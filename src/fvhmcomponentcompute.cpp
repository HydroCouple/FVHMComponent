#include "stdafx.h"
#include "fvhmcomponent.h"
#include "controlvolume.h"
#include "inletoutletflowbc.h"
#include "wsebc.h"
#include "initialwsebc.h"
#include "outletwseslope.h"
#include "precipbc.h"
#include "core/abstractoutput.h"
#include "criticaldepthoutflowbc.h"

#include <math.h>
#include <vector>

#ifdef USE_OPENMP
#include <omp.h>
#endif

#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef USE_SUITESPARSE
#include "SuiteSparseQR.hpp"
#endif


#include "netcdf"
#include "HYPRE.h"
#include "HYPRE_utilities.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_krylov.h"
#include "HYPRE_seq_mv.h"
#include "HYPRE_parcsr_mv.h"
#include "_hypre_parcsr_mv.h"

using namespace std;
using namespace HydroCouple::Temporal;

double FVHMComponent::defaultInletOutletWSESlope() const
{
  return m_inletOutletWSESlope;
}

double FVHMComponent::viscousSubLayerDepth() const
{
  return m_viscousSubLayerDepth;
}

double FVHMComponent::getNextTimeStep()
{
  double estimatedTimeStep = m_minTimeStep;

  if(m_initTimeStepCycle < m_numFixedTimeSteps)
  {
    estimatedTimeStep = m_minTimeStep;
    m_initTimeStepCycle ++;
  }
  else if(m_useAdaptiveTimeStep)
  {
    switch (m_adaptiveTSMode)
    {
      case AdaptiveTSMode::MaxCourantNumber:
        {
          double maxCFactor = 0;

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
          for(int i = 0 ; i < m_numCells  ; i++)
          {
            TriCV *cv = m_controlVolumes[i];

            if(cv->wetIndex)
            {
              double cFactor = cv->getCourantFactor();

              if(cFactor > maxCFactor)
              {

#ifdef USE_OPENMP
#pragma omp critical
#endif
                {
                  maxCFactor = cFactor;
                }
              }
            }
          }

          estimatedTimeStep = maxCFactor ? m_maxCourantNumber * m_timeStepRelaxFactor / maxCFactor : m_maxTimeStep;
        }
        break;
      case AdaptiveTSMode::RMSCourantNumber:
        {
          double count = 0;
          double tStep = 0;

          for(int i = 0 ; i < m_numCells ; i++)
          {
            TriCV *cv = m_controlVolumes[i];

            if(cv->wetIndex )
            {
              double cfactor = cv->getCourantFactor();
              tStep += cfactor * cfactor;
              count ++;
            }
          }

          if(count)
          {
            estimatedTimeStep = m_RMSCourantNumber / sqrt(tStep/count);
          }
        }
        break;
    }

    double dtStep = estimatedTimeStep - m_timeStep;
    double fract = dtStep / m_timeStep;

    if(dtStep > 0 && fabs(fract) > m_maxTimeStepCF)
    {
      estimatedTimeStep = m_timeStep + m_maxTimeStepCF * fract * m_timeStep / fabs(fract);
    }
  }

  if(estimatedTimeStep < m_minTimeStep)
  {
    estimatedTimeStep = m_minTimeStep;
  }
  else if(estimatedTimeStep > m_maxTimeStep)
  {
    estimatedTimeStep = m_maxTimeStep;
  }

  double nextDateTime = m_currentDateTime + estimatedTimeStep / 86400.0;

  if(nextDateTime > m_nextOutputTime  && m_nextOutputTime > m_currentDateTime)
  {
    estimatedTimeStep = (m_nextOutputTime - m_currentDateTime) * 86400.0;
    estimatedTimeStep = max(estimatedTimeStep, m_minTimeStep);
  }

  return estimatedTimeStep;
}

void FVHMComponent::applyInitialBoundaryConditions()
{

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0 ; i < m_numCells ; i++)
  {
    TriCV *cv = m_controlVolumes[i];
    cv->inflow = 0.0;
    cv->externalForce->v[0] = 0.0;
    cv->externalForce->v[1] = 0.0;

    for(int j = 0 ; j < cv->numEdges ;  j++)
    {
      int cvnIndex = cv->nTris[j];

      if(cvnIndex < 0 == false && cv->faceNormalVels[j].isBC)
      {
        FaceNormVelBC &faceVel = cv->faceNormalVels[j];
        faceVel.value = 0.0;
        faceVel.associatedValue = 0.0;
        faceVel.isBC = true;
        faceVel.velocityCalculateMode = FaceNormVelBC::PreCalculated;
        faceVel.initialiazeVelocityVariable();
      }
    }
  }

  readRestartFile();

  m_initialWSEBC->applyBoundaryConditions(m_currentDateTime, m_timeStep);

  int numWetCells;
  getWetCells(numWetCells);
  calculateCellWSEGradients();
  calculateCellVelocityGradients();
  calculateCellEddyViscosities();

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0 ; i < m_numCells ; i++)
  {
    TriCV *cv = m_controlVolumes[i];
    cv->copyVariablesToPrev();
  }
}

void FVHMComponent::applyBoundaryConditions(double time)
{

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0 ; i < m_numCells ; i++)
  {
    TriCV *cv = m_controlVolumes[i];
    cv->inflow = 0.0;
    cv->externalForce->v[0] = 0.0;
    cv->externalForce->v[1] = 0.0;
  }

#ifdef USE_OPENMP
#pragma omp parallel sections
#endif
  {
#ifdef USE_OPENMP
#pragma omp section
    {
#endif
      m_outletWSEBC->applyBoundaryConditions(time, m_timeStep);

      m_outletWSESlope->applyBoundaryConditions(time, m_timeStep);

      m_criticalDepthOutflowBC->applyBoundaryConditions(time, m_timeStep);

#ifdef USE_OPENMP
    }
#endif

#ifdef USE_OPENMP
#pragma omp section
#endif
    m_inletWSEBC->applyBoundaryConditions(time, m_timeStep);

#ifdef USE_OPENMP
#pragma omp section
#endif
    m_outletFlowBCArgument->applyBoundaryConditions(time, m_timeStep);

#ifdef USE_OPENMP
#pragma omp section
#endif
    m_inletFlowBCArgument->applyBoundaryConditions(time, m_timeStep);

#ifdef USE_OPENMP
#pragma omp section
#endif
    m_precipitationArgument->applyBoundaryConditions(time, m_timeStep);
  }

}

double FVHMComponent::getMinOutputTime(const QList<HydroCouple::IOutput *> &requiredOutputs)
{
  double minTime = std::numeric_limits<double>::max();

  for(const auto &output : requiredOutputs)
  {
    getMinIOutputTime(output, minTime);
  }

  return minTime;
}

void FVHMComponent::getMinIOutputTime(const HydroCouple::IOutput *output, double &minTime)
{
  for(HydroCouple::IInput *input : output->consumers())
  {
    ITimeComponentDataItem *tComponentDataItem =
        dynamic_cast<ITimeComponentDataItem*>(input);

    if(tComponentDataItem && tComponentDataItem->timeCount())
    {
      int lastTime = tComponentDataItem->timeCount() - 1;
      minTime = min(minTime,tComponentDataItem->time(lastTime)->modifiedJulianDay());
    }
  }

  for(HydroCouple::IAdaptedOutput *adaptedOutput : output->adaptedOutputs())
  {
    getMinIOutputTime(adaptedOutput, minTime);
  }
}

void FVHMComponent::prepareForNextTimeStep()
{

  double maxV = std::numeric_limits<double>::lowest();
  double minV = std::numeric_limits<double>::max();

  double maxU = std::numeric_limits<double>::lowest();
  double minU = std::numeric_limits<double>::max();

  double maxH = std::numeric_limits<double>::lowest();
  double minH = std::numeric_limits<double>::max();

  double maxZ = std::numeric_limits<double>::lowest();
  double minZ = std::numeric_limits<double>::max();

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0 ; i < m_numCells ; i++)
  {
    TriCV *cv = m_controlVolumes[i];
    cv->contResidualTotal += cv->contResidualIter;
    cv->velResidual[0] = cv->velResidualIter[0];
    cv->velResidual[1] = cv->velResidualIter[1];

    if(m_verbose && m_printFrequencyCounter <= 0)
    {
#ifdef USE_OPENMP
#pragma omp critical
#endif
      {
        maxU = max(maxU , cv->vel[0].value);
        maxV = max(maxV , cv->vel[1].value);

        minU = min(minU , cv->vel[0].value);
        minV = min(minV , cv->vel[1].value);

        maxH = max(maxH, cv->h->value);
        minH = min(minH, cv->h->value);

        maxZ = max(maxZ, cv->z->value);
        minZ = min(minZ, cv->z->value);
      }
    }

    cv->copyVariablesToPrev();
  }

  if(m_verbose && m_printFrequencyCounter <= 0)
    printf("Time: %f, MinZ: %f, MaxZ: %f, MinH: %f, MaxH: %f, MinU: %f, MaxU: %f, MinV: %f, MaxV: %f\n" , m_currentDateTime, minZ, maxZ, minH, maxH, minU, maxU, minV, maxV);

}

void FVHMComponent::getWetCells(int &numWetCells)
{
  numWetCells = 0;

  for(int i = 0; i < m_numCells; i++)
  {
    TriCV* cv = m_controlVolumes[i];

    if(cv->h->value > m_wetCellDepth)
    {
      cv->wetIndex = 1;
      cv->wetCellIndex = numWetCells;
      numWetCells++;
    }
    else
    {
      cv->wetIndex = 0;
      cv->wetCellIndex = -1;

      for(int c = 0 ; c < 2; c++)
      {
        cv->vel[c].value = 0.0;
        //        cv->grad_vel[c].v[0] = 0.0;
        //        cv->grad_vel[c].v[1] = 0.0;
        cv->prevIterVel[c] = 0.0;
        cv->velResidualIter[c] = 0.0;

      }

      //      cv->grad_z->v[0] = 0.0;
      //      cv->grad_z->v[1] = 0.0;

      //      for(int f = 0; f < cv->numEdges; f++)
      //      {
      //        FaceNormVelBC &faceVel = cv->faceNormalVels[f];

      //        if(faceVel.isBC == false && faceVel.value > 0.0)
      //        {
      //          faceVel.value = 0.0;
      //        }
      //      }
    }
  }
}

bool FVHMComponent::hasEdgeFlow(TriCV *cv)
{
  int cvnIndex = -1;
  for(int i = 0; i < cv->numEdges; i++)
  {
    if(cv->faceNormalVels[i].value ||  ((cvnIndex = cv->nTris[i]) > -1 && m_controlVolumes[cvnIndex]->wetIndex))
      return true;
  }

  return false;
}


FVHMComponent::ErrorCode FVHMComponent::performSimpleTimeStep(double tStep, int &numIterations, int &numWetCells,
                                                              double &uvelRelResidualNormInit, double &uvelRelResidualNormFin,
                                                              double &vvelRelResidualNormInit, double &vvelRelResidualNormFin,
                                                              double &pressureRelRisdualNormInit, double &pressureRelRisdualNormFin,
                                                              double &continuityRelResidualNormInit, double &continuityRelResidualNormFin,
                                                              bool &converged, QString &errorMessage)
{

  ErrorCode error = ErrorCode::NoError;

  numIterations = 0;

  uvelRelResidualNormFin = vvelRelResidualNormFin = pressureRelRisdualNormFin =
      continuityRelResidualNormFin = uvelRelResidualNormInit = vvelRelResidualNormInit =
      pressureRelRisdualNormInit = continuityRelResidualNormInit = 0.0;

  converged = false;

  numWetCells = m_numCells;

  getWetCells(numWetCells);

#ifdef USE_OPENMP
#pragma omp parallel sections
#endif
  {
#ifdef USE_OPENMP
#pragma omp section
#endif
    calculateFriction(0);

#ifdef USE_OPENMP
#pragma omp section
#endif
    calculateFriction(1);
  }


  while(numIterations < m_itersPerTimeStep)
  {

    m_currentCoeff = m_currentCoeff == 0 ? 0 : 0;

    ErrorCode uerror = ErrorCode::NoError;
    ErrorCode verror = ErrorCode::NoError;
    QString uErrorMessage;
    QString vErrorMessage;
    uvelRelResidualNormFin = vvelRelResidualNormFin =
        pressureRelRisdualNormFin = continuityRelResidualNormFin = 0.0;


#ifdef USE_OPENMP
#pragma omp parallel sections
#endif
    {
#ifdef USE_OPENMP
#pragma omp section
#endif
      uerror = solveMomentumEquations(tStep, uvelRelResidualNormFin, numWetCells, 0, uErrorMessage);

#ifdef USE_OPENMP
#pragma omp section
#endif
      verror = solveMomentumEquations(tStep, vvelRelResidualNormFin, numWetCells, 1, vErrorMessage);
    }

    if(uerror == ErrorCode::CriticalFailure)
    {
      error = uerror;
      errorMessage = uErrorMessage;
      break;
    }

    if(verror == ErrorCode::CriticalFailure)
    {
      error = verror;
      errorMessage = vErrorMessage;
      break;
    }

    calculateFaceVelocities();

    error = solvePressureCorrection(tStep, pressureRelRisdualNormFin, 0.0, errorMessage);

    if(error == ErrorCode::CriticalFailure)
    {
      break;
    }

#ifdef USE_OPENMP
#pragma omp parallel sections
#endif
    {
#ifdef USE_OPENMP
#pragma omp section
#endif
      applyPressureCorrections();

#ifdef USE_OPENMP
#pragma omp section
#endif
      applyVelocityCorrections();
    }

    calculateCellWSEGradients();

    calculateFaceVelocities();

    calculateContinuityResiduals(tStep, continuityRelResidualNormFin);

    if(numIterations == 0)
    {
      uvelRelResidualNormInit = uvelRelResidualNormFin;
      vvelRelResidualNormInit = vvelRelResidualNormFin;
      pressureRelRisdualNormInit = pressureRelRisdualNormFin;
      continuityRelResidualNormInit = continuityRelResidualNormFin;
    }

    //check for convergence
    if(uvelRelResidualNormFin < m_uvelConvergenceTol &&
       vvelRelResidualNormFin < m_vvelConvergenceTol &&
       fabs(pressureRelRisdualNormFin) < m_pressureConvergenceTol
       && continuityRelResidualNormFin < m_contConvergenceTol
       )
    {
      numIterations++;
      converged = true;
      break;
    }

    calculateCellVelocityGradients();

    calculateCellEddyViscosities();

    numIterations++;

  }

  if(error == ErrorCode::NoError)
  {
    if(!converged)
    {
      errorMessage = "FVHMComponent Id:" + id() + " failed to converged! at time " + m_currentDateTime;
      error = ErrorCode::FailedToConverge;
    }
    else
    {
      errorMessage = "";
      return ErrorCode::NoError;
    }
  }

  return error;
}


FVHMComponent::ErrorCode FVHMComponent::performSimpleCTimeStep(double tStep, int &numIterations, int &numWetCells,
                                                               double &uvelRelResidualNormInit, double &uvelRelResidualNormFin,
                                                               double &vvelRelResidualNormInit, double &vvelRelResidualNormFin,
                                                               double &pressureRelRisdualNormInit, double &pressureRelRisdualNormFin,
                                                               double &continuityRelResidualNormInit, double &continuityRelResidualNormFin,
                                                               bool &converged, QString &errorMessage)
{
  ErrorCode error = ErrorCode::NoError;

  numIterations = 0;

  uvelRelResidualNormFin = vvelRelResidualNormFin = pressureRelRisdualNormFin =
      continuityRelResidualNormFin = uvelRelResidualNormInit = vvelRelResidualNormInit =
      pressureRelRisdualNormInit = continuityRelResidualNormInit = 0.0;

  converged = false;

  numWetCells = m_numCells;

  getWetCells(numWetCells);


#ifdef USE_OPENMP
#pragma omp parallel sections
#endif
  {

#ifdef USE_OPENMP
#pragma omp section
#endif
    calculateFriction(0);

#ifdef USE_OPENMP
#pragma omp section
#endif
    calculateFriction(1);
  }


  while(numIterations < m_itersPerTimeStep)
  {
    uvelRelResidualNormFin = vvelRelResidualNormFin = pressureRelRisdualNormFin = continuityRelResidualNormFin = 0.0;

    ErrorCode uerror = ErrorCode::NoError;
    ErrorCode verror = ErrorCode::NoError;
    QString uErrorMessage;
    QString vErrorMessage;

#ifdef USE_OPENMP
#pragma omp parallel sections
#endif
    {
#ifdef USE_OPENMP
#pragma omp section
#endif
      uerror = solveMomentumEquations(tStep, uvelRelResidualNormFin, numWetCells, 0, uErrorMessage);

#ifdef USE_OPENMP
#pragma omp section
#endif
      verror = solveMomentumEquations(tStep, vvelRelResidualNormFin, numWetCells, 1, vErrorMessage);
    }

    if(uerror == ErrorCode::CriticalFailure)
    {
      error = uerror;
      errorMessage = uErrorMessage;
      break;
    }

    if(verror == ErrorCode::CriticalFailure)
    {
      error = verror;
      errorMessage = vErrorMessage;
      break;
    }

    calculateFaceVelocities();

    error = solvePressureCorrection(tStep, pressureRelRisdualNormFin, 1.0, errorMessage);

    if(error == ErrorCode::CriticalFailure)
    {
      break;
    }

#ifdef USE_OPENMP
#pragma omp parallel sections
#endif
    {

#ifdef USE_OPENMP
#pragma omp section
#endif
      applyPressureCorrections();

#ifdef USE_OPENMP
#pragma omp section
#endif
      applyVelocityCorrections();

    }

    calculateCellWSEGradients();

    calculateFaceVelocities();

    calculateContinuityResiduals(tStep, continuityRelResidualNormFin);

    if(numIterations == 0)
    {
      uvelRelResidualNormInit = uvelRelResidualNormFin;
      vvelRelResidualNormInit = vvelRelResidualNormFin;
      pressureRelRisdualNormInit = pressureRelRisdualNormFin;
      continuityRelResidualNormInit = continuityRelResidualNormFin;
    }

    //check for convergence
    if(uvelRelResidualNormFin < m_uvelConvergenceTol &&
       vvelRelResidualNormFin < m_vvelConvergenceTol &&
       fabs(pressureRelRisdualNormFin) < m_pressureConvergenceTol &&
       continuityRelResidualNormFin < m_contConvergenceTol)
    {
      numIterations++;
      converged = true;
      break;
    }

    calculateCellVelocityGradients();

    calculateCellEddyViscosities();

    numIterations++;

  }

  if(error == ErrorCode::NoError)
  {
    if(!converged)
    {
      errorMessage = "FVHMComponent Id:" + id() + " failed to converged! at time " + m_currentDateTime;
      error = ErrorCode::FailedToConverge;
    }
    else
    {
      errorMessage = "";
      return ErrorCode::NoError;
    }
  }

  return error;
}

FVHMComponent::ErrorCode FVHMComponent::performPISOTimeStep(double tStep, int &numIterations, int &numWetCells,
                                                            double &uvelRelResidualNormInit, double &uvelRelResidualNormFin,
                                                            double &vvelRelResidualNormInit, double &vvelRelResidualNormFin,
                                                            double &pressureRelRisdualNormInit, double &pressureRelRisdualNormFin,
                                                            double &continuityRelResidualNormInit, double &continuityRelResidualNormFin,
                                                            bool &converged, QString &errorMessage)
{
  ErrorCode error = ErrorCode::NoError;

  numIterations = 0;

  uvelRelResidualNormFin = vvelRelResidualNormFin = pressureRelRisdualNormFin =
      continuityRelResidualNormFin = uvelRelResidualNormInit = vvelRelResidualNormInit =
      pressureRelRisdualNormInit = continuityRelResidualNormInit = 0.0;

  converged = false;

  numWetCells = m_numCells;

  getWetCells(numWetCells);


#ifdef USE_OPENMP
#pragma omp parallel sections
#endif
  {

#ifdef USE_OPENMP
#pragma omp section
#endif
    calculateFriction(0);

#ifdef USE_OPENMP
#pragma omp section
#endif
    calculateFriction(1);
  }

  while(numIterations < m_itersPerTimeStep)
  {
    uvelRelResidualNormFin = vvelRelResidualNormFin =
        pressureRelRisdualNormFin = continuityRelResidualNormFin = 0.0;

    ErrorCode uerror = ErrorCode::NoError;
    ErrorCode verror = ErrorCode::NoError;
    QString uErrorMessage;
    QString vErrorMessage;

#ifdef USE_OPENMP
#pragma omp parallel sections
#endif
    {
#ifdef USE_OPENMP
#pragma omp section
#endif
      uerror = solveMomentumEquations(tStep, uvelRelResidualNormFin, numWetCells, 0, uErrorMessage);

#ifdef USE_OPENMP
#pragma omp section
#endif
      verror = solveMomentumEquations(tStep, vvelRelResidualNormFin, numWetCells, 1, vErrorMessage);
    }

    if(uerror == ErrorCode::CriticalFailure)
    {
      error = uerror;
      errorMessage = uErrorMessage;
      break;
    }

    if(verror == ErrorCode::CriticalFailure)
    {
      error = verror;
      errorMessage = vErrorMessage;
      break;
    }

    calculateFaceVelocities();

    for(int i = 0 ; i < 2; i++)
    {
      error = solvePressureCorrection(tStep, pressureRelRisdualNormFin, 0.0, errorMessage);

#ifdef USE_OPENMP
#pragma omp parallel sections
#endif
      {

#ifdef USE_OPENMP
#pragma omp section
#endif
        applyPressureCorrections();

#ifdef USE_OPENMP
#pragma omp section
#endif
        applyVelocityCorrections();

      }

      calculateCellWSEGradients();

      calculateFaceVelocities();
    }

    if(error == ErrorCode::CriticalFailure)
    {
      break;
    }

    calculateContinuityResiduals(tStep, continuityRelResidualNormFin);

    if(numIterations == 0)
    {
      uvelRelResidualNormInit = uvelRelResidualNormFin;
      vvelRelResidualNormInit = vvelRelResidualNormFin;
      pressureRelRisdualNormInit = pressureRelRisdualNormFin;
      continuityRelResidualNormInit = continuityRelResidualNormFin;
    }

    //check for convergence
    if(uvelRelResidualNormFin < m_uvelConvergenceTol &&
       vvelRelResidualNormFin < m_vvelConvergenceTol &&
       fabs(pressureRelRisdualNormFin) < m_pressureConvergenceTol &&
       continuityRelResidualNormFin < m_contConvergenceTol)
    {
      numIterations++;
      converged = true;
      break;
    }

    calculateCellVelocityGradients();

    //calculate turbulent viscousities
    calculateCellEddyViscosities();

    numIterations++;
  }

  if(error == ErrorCode::NoError)
  {
    if(!converged)
    {
      errorMessage = "FVHMComponent Id:" + id() + " failed to converged! at time " + m_currentDateTime;
      error = ErrorCode::FailedToConverge;
    }
    else
    {
      errorMessage = "";
      return ErrorCode::NoError;
    }
  }

  return error;
}

FVHMComponent::ErrorCode FVHMComponent::solveMomentumEquations(double tStep, double &uvelRelResidualNorm, int numWetCells, int uv, QString &errorMessage)
{
  //momentum
  double *u_b = new double[numWetCells];
  SparseMatrix coeffMatrix(numWetCells, 4);
  double *x = new double[numWetCells];

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int index = 0 ; index < m_numCells ; index++)
  {

    TriCV *cv = m_controlVolumes[index];

    for(int u = 0; u < cv->numEdges + 1 ; u++)
    {
      cv->velCoeffs[uv][u] = 0.0;
    }

    double a_p = cv->h->value * cv->area / tStep;

    if(cv->wetIndex)
    {
      double b = 0.0;

      for(int i = 0; i < cv->numEdges; i++)
      {
        FaceNormVelBC &faceVel = cv->faceNormalVels[i];
        double evel = faceVel.value;

        if(evel)
        {
          TriCV *cvn = nullptr;
          int cvIndex = cv->nTris[i];
          double edgeDepth = cv->faceDepths[i].value;
          double tempCoeff = m_useAdvection * evel * cv->r_eta[i] * edgeDepth;

          if(cvIndex > -1 && (cvn = m_controlVolumes[cvIndex])->wetIndex)
          {

            double a_n = 0.0;
            double wl_n = cv->w_l[i];
            double wl_p = 1.0 - wl_n;

            Vect &r_xi = cv->r_xi[i];

            //TVD
            if(evel > 0)
            {
              double ru = calculateFluxLimiter(cv, i, wl_n, cv->grad_vel[uv], cvn->grad_vel[uv], cv->vel[uv].value, cvn->vel[uv].value, r_xi.x() , r_xi.y());
              a_p += tempCoeff * (1.0 - ru * wl_n);
              a_n += tempCoeff * ru * wl_n;
            }
            else
            {
              double ru = calculateFluxLimiter(cv, i, wl_p, cvn->grad_vel[uv], cv->grad_vel[uv], cvn->vel[uv].value, cv->vel[uv].value, -r_xi.x() , -r_xi.y());
              a_p += tempCoeff *  ru * wl_p;
              a_n += tempCoeff *  (1.0 - ru * wl_p);
            }

            //            double gradVelX = wl_p * cv->grad_vel[uv].v[0] + wl_n * cvn->grad_vel[uv].v[0];
            //            double gradVelY = wl_p * cv->grad_vel[uv].v[1] + wl_n * cvn->grad_vel[uv].v[1];
            //            b -= tempCoeff * Vect::dotProduct(cv->df[i].v[0], cv->df[i].v[1], 0.0, gradVelX, gradVelY, 0.0);

            double faceTurbulentVisco = wl_p * cv->eddyViscosity + wl_n * cvn->eddyViscosity;

            //direct diffusion
            a_n -= m_useEddyViscosity * edgeDepth * (faceTurbulentVisco + m_viscosity) * cv->dir_diff_term[i];
            a_p += m_useEddyViscosity * edgeDepth * (faceTurbulentVisco + m_viscosity) * cv->dir_diff_term[i];

            //cross diffusion
            b += m_useEddyViscosity * edgeDepth * (faceTurbulentVisco + m_viscosity) *
                 cv->cross_diff_term[i] * (cv->nodeVels[uv][i + 1] - cv->nodeVels[uv][i]);

            cv->velCoeffs[uv][i + 1] = a_n;

          }
          else if(faceVel.isBC)
          {
            b -= tempCoeff * faceVel.vel->v[uv];
          }
        }
      }

      cv->velCoeffs[uv][0] = a_p;

      for(int i = 0; i < cv->numEdges + 1; i++)
      {
        int cvnIndex = cv->orderedNTris[i];
        double coeff = 0.0;

        if(cvnIndex > -1 && (coeff = cv->velCoeffs[uv][cv->orderedNTrisIndexes[i]]))
        {
          TriCV *cvn = m_controlVolumes[cvnIndex];
          coeffMatrix.appendValue(cv->wetCellIndex, cvn->wetCellIndex, coeff);
        }
      }

      //const

      //del_v/del_t
      double dv_dt = cv->prevH->value * cv->prevVel[uv].value * cv->area / tStep;

      //pressure gradient
      double gradp = -g * cv->grad_z->v[uv] * cv->area * cv->h->value;

      //external force
      double extforce = cv->externalForce->v[uv] * cv->area * cv->h->value;

      double friction = 0.0;

      //friction
      {
        friction = -cv->friction[uv];

        for(int f = 0 ; f < cv->numEdges; f++)
        {
          friction += -m_useWall * cv->wallShearFriction[uv][f];
        }
      }

      b += dv_dt + gradp + friction + extforce;


      //const coefficients
      u_b[cv->wetCellIndex] = b;

      //initial guess
      x[cv->wetCellIndex] = cv->vel[uv].value;

    }
    else
    {
      cv->velCoeffs[uv][0] = a_p;
    }
  }

  //set up momentum equations calculate u and solve
  double *residuals = new double[numWetCells]();
  ErrorCode error = solve(coeffMatrix, u_b, x, residuals, uvelRelResidualNorm,errorMessage);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int index = 0; index < m_numCells ; index++)
  {
    TriCV *cv = m_controlVolumes[index];

    if(cv->wetIndex)
    {
      cv->prevIterVel[uv] = cv->vel[uv].value;
      cv->vel[uv].value = x[cv->wetCellIndex];
      cv->velResidualIter[uv] = residuals[cv->wetCellIndex];
    }
  }

  delete[] u_b;
  delete[] x;
  delete[] residuals;

  return error;
}

void FVHMComponent::interpolateWettedCellVelocities()
{
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int c = 0 ; c < m_numCells ; c++)
  {
    TriCV *cv = m_controlVolumes[c];

    if((cv->wetIndex == 0 || cv->wetIndex == 2) )
    {
      if(cv->h->value > m_wetCellDepth)
      {
        cv->wetIndex = 2;
        interpolateWettedCellVelocity(cv);
      }
      else
      {
        cv->wetIndex = 0;
        cv->vel[0].value = 0.0;
        cv->vel[1].value = 0.0;
      }
    }
  }
}

void FVHMComponent::interpolateWettedCellVelocity(TriCV *cv)
{
  int index = -1;
  double maxH = std::numeric_limits<double>::lowest();

  for(int i = 0 ; i < cv->numEdges; i++)
  {
    int cvnIndex = cv->nTris[i];
    TriCV* cvn = nullptr;

    if(cvnIndex > -1 && (cvn = m_controlVolumes[cvnIndex])->wetIndex == 1 && cvn->h->value > maxH)
    {
      index  = i;
      maxH = cvn->h->value;
    }
  }

  if(index > -1)
  {
    TriCV* cvn = this->m_controlVolumes[cv->nTris[index]];
    double u = cvn->vel[0].value + Vect::dotProduct(cvn->r_xi[index] * 1.0 , cvn->grad_vel[0]);
    double v = cvn->vel[1].value + Vect::dotProduct(cvn->r_xi[index] * 1.0 , cvn->grad_vel[1]);
    cv->vel[0].value = u;
    cv->vel[1].value = v;
  }
  else
  {
    cv->vel[0].value = 0;
    cv->vel[1].value = 0;
  }
}

void FVHMComponent::calculateFaceVelocities()
{

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int c = 0 ; c < m_numCells ; c++)
  {
    m_controlVolumes[c]->calculateFaceVelocities();
  }
}

void FVHMComponent::calculateBoundaryFaceVelocities()
{
  //#ifdef USE_OPENMP
  //#pragma omp parallel for
  //#endif
  //  for(int c = 0 ; c < m_numCells ; c++)
  //  {
  //    TriCV *cv = m_controlVolumes[c];
  //    calculateBoundaryFaceVelocities(cv);
  //  }
}

void FVHMComponent::calculateBoundaryFaceVelocities(TriCV *cv)
{
  for(int i = 0 ; i < cv->numEdges ; i++)
  {
    FaceNormVelBC &faceNormVel =  cv->faceNormalVels[i];

    if(faceNormVel.isBC)
    {
      if(faceNormVel.calculateFromQ)
      {
        double area = cv->faceDepths[i].value * cv->r_eta[i];

        if(fabs(area - 0.0) > m_epsilon)
        {
          faceNormVel.value = faceNormVel.associatedValue / area;
        }
      }

      switch (faceNormVel.velocityCalculateMode)
      {
        case FaceNormVelBC::PreCalculated:
          {
          }
          break;
        case FaceNormVelBC::DotProductEdgeVel:
          {
            faceNormVel.vel->v[0] = faceNormVel.value * cv->e_n[i].v[0];
            faceNormVel.vel->v[1] = faceNormVel.value * cv->e_n[i].v[1];
          }
          break;
        case FaceNormVelBC::CellVel:
          {
            faceNormVel.vel->v[0]  = cv->vel[0].value;
            faceNormVel.vel->v[1]  = cv->vel[1].value;
          }
          break;
        case FaceNormVelBC::CellGradient:
          {
            faceNormVel.vel->v[0] = cv->vel[0].value + Vect::dotProduct(cv->grad_vel[0], cv->r_e[i]);
            faceNormVel.vel->v[1] = cv->vel[1].value + Vect::dotProduct(cv->grad_vel[1], cv->r_e[i]);
          }
          break;
        case FaceNormVelBC::ZeroGradient:
          {
            std::tuple<double, double, double> vel = calculateZeroGradientFaceVelocity(cv, i);
            faceNormVel.vel->v[0]  = std::get<0>(vel);
            faceNormVel.vel->v[1]  = std::get<1>(vel);
            faceNormVel.value = std::get<2>(vel);
            faceNormVel.associatedValue = faceNormVel.value * cv->r_eta[i] * cv->faceDepths[i].value;
          }
          break;
      }
    }
  }
}

/*!
 * \brief FVHMComponent::calculateFluxLimiter
 * \param cv
 * \param faceIndex
 * \param uvel
 * \brief Skewness correction added from Denner, F. and B.G.M. van Wachem, 2015. TVD Differencing on Three-Dimensional Unstructured
 * Meshes with Monotonicity-Preserving Correction of Mesh Skewness. Journal of Computational Physics 298:466–479.
 * \return
 */
double FVHMComponent::calculateFluxLimiter(const TriCV* cv, int faceIndex, double idwfactor, const Vect &grad_c, const Vect& grad_d, double phi_c, double phi_d, double r_cdX, double r_cdY)
{
  //upwind default
  double ru = 0.0;

  switch (m_advectionScheme)
  {
    //Central difference
    case 1:
      {
        ru = 1.0;
      }
      break;
      //Quick
    case 2:
      {
        double du = phi_d - phi_c;

        if(fabs(du - 0.0) >= m_epsilon)
        {
          double rf = 0.0;

          //Darwish and Moukalled (2002)
          {
            rf = Vect::dotProduct(2.0 * grad_c.v[0], 2.0 * grad_c.v[1], 0.0, r_cdX, r_cdY, 0.0) / du;
            rf = rf - 1.0;
          }

          //Ubbink and Issa (1999)
          //          {
          //            double gradFaceX = (1.0 - idwfactor)* grad_c.v[0] + idwfactor * grad_d.v[0];
          //            double gradFaceY = (1.0 - idwfactor)* grad_c.v[1] + idwfactor * grad_d.v[1];
          //            double phi_u = phi_d - 2.0 * Vect::dotProduct(gradFaceX, gradFaceY, 0.0, dir,dir, 0.0);
          //            rf = (phi_c - phi_u)/(phi_d - phi_c);
          //          }

          ru = max({0.0, min({2.0 * rf, (3.0 + rf) / 4.0 , 2.0})});

          //Skewness correction FabianDenner∗, Berend G.M.vanWachem
          {
            //            ru = ru + Vect::dotProduct(gradFaceX / idwfactor, gradFaceY / idwfactor, 0.0, cv->df[faceIndex]) / du;

            //            //1st order
            //            {
            //              double maxRu = max(0.0, min(rf/idwfactor, 1.0/idwfactor));
            //              double minRu = 0.0;
            //              ru = min(max(minRu, ru), maxRu);
            //            }

            //2nd order
            //          {
            //            double maxRu = max({0.0, min(rf/idwfactor, 1.0), min(rf, 1.0/ wl)});
            //            double minRu = max(0.0, min(rf,1.0));
            //            ru = min(max(minRu, ru), maxRu);
            //          }
          }

        }
      }
      break;
      //Min Mod
    case 3:
      {
        double du = phi_d - phi_c;

        if(fabs(du - 0.0) >= m_epsilon)
        {
          double rf = 0.0;

          //Darwish and Moukalled (2002)
          {
            rf = Vect::dotProduct(2.0 * grad_c.v[0], 2.0 * grad_c.v[1], 0.0, r_cdX, r_cdY, 0.0) / du;
            rf = rf - 1.0;
          }

          //Ubbink and Issa (1999)
          //          {
          //            double gradFaceX = (1.0 - idwfactor)* grad_c.v[0] + idwfactor * grad_d.v[0];
          //            double gradFaceY = (1.0 - idwfactor)* grad_c.v[1] + idwfactor * grad_d.v[1];
          //            double phi_u = phi_d - 2.0 * Vect::dotProduct(gradFaceX, gradFaceY, 0.0, dir,dir, 0.0);
          //            rf = (phi_c - phi_u)/(phi_d - phi_c);
          //          }

          ru = rf > 0 ? min(rf,1.0) : 0.0;

          //Skewness correction FabianDenner∗, Berend G.M.vanWachem
          {
            //            ru = ru + Vect::dotProduct(gradFaceX / idwfactor, gradFaceY / idwfactor, 0.0, cv->df[faceIndex]) / du;

            //            //1st order
            //            {
            //              double maxRu = max(0.0, min(rf/idwfactor, 1.0/idwfactor));
            //              double minRu = 0.0;
            //              ru = min(max(minRu, ru), maxRu);
            //            }

            //2nd order
            //          {
            //            double maxRu = max({0.0, min(rf/idwfactor, 1.0), min(rf, 1.0/ wl)});
            //            double minRu = max(0.0, min(rf,1.0));
            //            ru = min(max(minRu, ru), maxRu);
            //          }
          }

        }
      }
      break;
      //Van Leer
    case 4:
      {
        double du = phi_d - phi_c;

        if(fabs(du - 0.0) >= m_epsilon)
        {
          double rf = 0.0;

          //Darwish and Moukalled (2002)
          {
            rf = Vect::dotProduct(2.0 * grad_c.v[0], 2.0 * grad_c.v[1], 0.0, r_cdX, r_cdY, 0.0) / du;
            rf = rf - 1.0;
          }

          //Ubbink and Issa (1999)
          //          {
          //          double gradFaceX = (1.0 - idwfactor)* grad_c.v[0] + idwfactor * grad_d.v[0];
          //          double gradFaceY = (1.0 - idwfactor)* grad_c.v[1] + idwfactor * grad_d.v[1];
          //            double phi_u = phi_d - 2.0 * Vect::dotProduct(gradFaceX, gradFaceY, 0.0, dir,dir, 0.0);
          //            rf = (phi_c - phi_u)/(phi_d - phi_c);
          //          }

          ru = (rf + fabs(rf))/(1.0 + fabs(rf));

          //Skewness correction FabianDenner∗, Berend G.M.vanWachem
          {
            //            ru = ru + Vect::dotProduct(gradFaceX / idwfactor, gradFaceY / idwfactor, 0.0, cv->df[faceIndex]) / du;

            //1st order
            //            {
            //              double maxRu = max(0.0, min(rf/idwfactor, 1.0/idwfactor));
            //              double minRu = 0.0;
            //              ru = min(max(minRu, ru), maxRu);
            //            }

            //2nd order
            //          {
            //            double maxRu = max({0.0, min(rf/idwfactor, 1.0), min(rf, 1.0/ wl)});
            //            double minRu = max(0.0, min(rf,1.0));
            //            ru = min(max(minRu, ru), maxRu);
            //          }
          }

        }
      }
      break;
      //Van Albada
    case 5:
      {
        double du = phi_d - phi_c;

        if(fabs(du - 0.0) >= m_epsilon)
        {
          double rf = 0.0;

          //Darwish and Moukalled (2002)
          {
            rf = Vect::dotProduct(2.0 * grad_c.v[0], 2.0 * grad_c.v[1], 0.0, r_cdX, r_cdY, 0.0) / du;
            rf = rf - 1.0;
          }

          //Ubbink and Issa (1999)
          //          {
          //          double gradFaceX = (1.0 - idwfactor)* grad_c.v[0] + idwfactor * grad_d.v[0];
          //          double gradFaceY = (1.0 - idwfactor)* grad_c.v[1] + idwfactor * grad_d.v[1];
          //            double phi_u = phi_d - 2.0 * Vect::dotProduct(gradFaceX, gradFaceY, 0.0, dir,dir, 0.0);
          //            rf = (phi_c - phi_u)/(phi_d - phi_c);
          //          }

          ru = (rf + rf * rf)/(1 + rf * rf);

          //Skewness correction FabianDenner∗, Berend G.M.vanWachem
          {
            //            ru = ru + Vect::dotProduct(gradFaceX / idwfactor, gradFaceY / idwfactor, 0.0, cv->df[faceIndex]) / du;

            //            //1st order
            //            {
            //              double maxRu = max(0.0, min(rf/idwfactor, 1.0/idwfactor));
            //              double minRu = 0.0;
            //              ru = min(max(minRu, ru), maxRu);
            //            }

            //2nd order
            //          {
            //            double maxRu = max({0.0, min(rf/idwfactor, 1.0), min(rf, 1.0/ wl)});
            //            double minRu = max(0.0, min(rf,1.0));
            //            ru = min(max(minRu, ru), maxRu);
            //          }
          }
        }
      }
      break;
      //Superbee
    case 6:
      {
        double du = phi_d - phi_c;

        if(fabs(du - 0.0) >= m_epsilon)
        {
          double rf = 0.0;

          //Darwish and Moukalled (2002)
          {
            rf = Vect::dotProduct(2.0 * grad_c.v[0], 2.0 * grad_c.v[1], 0.0, r_cdX, r_cdY, 0.0) / du;
            rf = rf - 1.0;
          }


          //Ubbink and Issa (1999)
          //          {
          //          double gradFaceX = (1.0 - idwfactor)* grad_c.v[0] + idwfactor * grad_d.v[0];
          //          double gradFaceY = (1.0 - idwfactor)* grad_c.v[1] + idwfactor * grad_d.v[1];
          //            double phi_u = phi_d - 2.0 * Vect::dotProduct(gradFaceX, gradFaceY, 0.0, dir,dir, 0.0);
          //            rf = (phi_c - phi_u)/(phi_d - phi_c);
          //          }

          ru = max(0.0, max(min(2.0 * rf , 1.0), min(rf,2.0)));

          //Skewness correction FabianDenner∗, Berend G.M.vanWachem
          {
            //            ru = ru + Vect::dotProduct(gradFaceX / idwfactor, gradFaceY / idwfactor, 0.0, cv->df[faceIndex]) / du;

            //            //1st order
            //            {
            //              double maxRu = max(0.0, min(rf/idwfactor, 1.0/idwfactor));
            //              double minRu = 0.0;
            //              ru = min(max(minRu, ru), maxRu);
            //            }

            //2nd order
            //          {
            //            double maxRu = max({0.0, min(rf/idwfactor, 1.0), min(rf, 1.0/ wl)});
            //            double minRu = max(0.0, min(rf,1.0));
            //            ru = min(max(minRu, ru), maxRu);
            //          }
          }

        }
      }
      break;
      //Umist
    case 7:
      {
        double du = phi_d - phi_c;

        if(fabs(du - 0.0) >= m_epsilon)
        {
          double rf = 0.0;

          //Darwish and Moukalled (2002)
          {
            rf = Vect::dotProduct(2.0 * grad_c.v[0], 2.0 * grad_c.v[1], 0.0, r_cdX, r_cdY, 0.0) / du;
            rf = rf - 1.0;
          }


          //Ubbink and Issa (1999)
          //          {
          //          double gradFaceX = (1.0 - idwfactor)* grad_c.v[0] + idwfactor * grad_d.v[0];
          //          double gradFaceY = (1.0 - idwfactor)* grad_c.v[1] + idwfactor * grad_d.v[1];
          //            double phi_u = phi_d - 2.0 * Vect::dotProduct(gradFaceX, gradFaceY, 0.0, dir,dir, 0.0);
          //            rf = (phi_c - phi_u)/(phi_d - phi_c);
          //          }

          ru = max(0.0,min(min( min(2.0 * rf, (1.0 + 3.0 * rf)/4.0) ,(3.0 + rf)/4.0),2.0));

          //Skewness correction FabianDenner∗, Berend G.M.vanWachem
          {
            //            ru = ru + Vect::dotProduct(gradFaceX / idwfactor, gradFaceY / idwfactor, 0.0, cv->df[faceIndex]) / du;

            //            //1st order
            //            {
            //              double maxRu = max(0.0, min(rf/idwfactor, 1.0/idwfactor));
            //              double minRu = 0.0;
            //              ru = min(max(minRu, ru), maxRu);
            //            }

            //2nd order
            //          {
            //            double maxRu = max({0.0, min(rf/idwfactor, 1.0), min(rf, 1.0/ wl)});
            //            double minRu = max(0.0, min(rf,1.0));
            //            ru = min(max(minRu, ru), maxRu);
            //          }
          }
        }
      }
      break;
      //Quadratic
    case 8:
      {
        double du = phi_d - phi_c;

        if(fabs(du - 0.0) >= m_epsilon)
        {
          double rf = 0.0;

          //Darwish and Moukalled (2002)
          {
            rf = Vect::dotProduct(2.0 * grad_c.v[0], 2.0 * grad_c.v[1], 0.0, r_cdX, r_cdY, 0.0) / du;
            rf = rf - 1.0;
          }

          //Ubbink and Issa (1999)
          //          {
          //            double gradFaceX = (1.0 - idwfactor)* grad_c.v[0] + idwfactor * grad_d.v[0];
          //            double gradFaceY = (1.0 - idwfactor)* grad_c.v[1] + idwfactor * grad_d.v[1];
          //            double phi_u = phi_d - 2.0 * Vect::dotProduct(gradFaceX, gradFaceY, 0.0, dir,dir, 0.0);
          //            rf = (phi_c - phi_u)/(phi_d - phi_c);
          //          }

          ru = rf <= 2.0 ? (2.0 * rf + rf*rf)/(2.0 + rf + rf*rf) : 1.0 ;

          //Skewness correction FabianDenner∗, Berend G.M.vanWachem
          {
            //            ru = ru + Vect::dotProduct(gradFaceX / idwfactor, gradFaceY / idwfactor, 0.0, cv->df[faceIndex]) / du;

            //            //1st order
            //            {
            //              double maxRu = max(0.0, min(rf/idwfactor, 1.0/idwfactor));
            //              double minRu = 0.0;
            //              ru = min(max(minRu, ru), maxRu);
            //            }

            //2nd order
            //          {
            //            double maxRu = max({0.0, min(rf/idwfactor, 1.0), min(rf, 1.0/ wl)});
            //            double minRu = max(0.0, min(rf,1.0));
            //            ru = min(max(minRu, ru), maxRu);
            //          }
          }
        }
      }
      break;
      //Cubic
    case 9:
      {
        double du = phi_d - phi_c;

        if(fabs(du - 0.0) >= m_epsilon)
        {
          double rf = 0.0;

          //Darwish and Moukalled (2002)
          {
            rf = Vect::dotProduct(2.0 * grad_c.v[0], 2.0 * grad_c.v[1], 0.0, r_cdX, r_cdY, 0.0) / du;
            rf = rf - 1.0;
          }

          //Ubbink and Issa (1999)
          //          {
          //            double gradFaceX = (1.0 - idwfactor)* grad_c.v[0] + idwfactor * grad_d.v[0];
          //            double gradFaceY = (1.0 - idwfactor)* grad_c.v[1] + idwfactor * grad_d.v[1];
          //            double phi_u = phi_d - 2.0 * Vect::dotProduct(gradFaceX, gradFaceY, 0.0, dir,dir, 0.0);
          //            rf = (phi_c - phi_u)/(phi_d - phi_c);
          //          }

          ru = rf <= 2.0 ? (4.0*rf + rf*rf*rf)/(4.0 + rf*rf + rf*rf*rf) : 1.0 ;

          //Skewness correction FabianDenner∗, Berend G.M.vanWachem
          {
            //            ru = ru + Vect::dotProduct(gradFaceX / idwfactor, gradFaceY / idwfactor, 0.0, cv->df[faceIndex]) / du;

            //            //1st order
            //            {
            //              double maxRu = max(0.0, min(rf/idwfactor, 1.0/idwfactor));
            //              double minRu = 0.0;
            //              ru = min(max(minRu, ru), maxRu);
            //            }

            //2nd order
            //          {
            //            double maxRu = max({0.0, min(rf/idwfactor, 1.0), min(rf, 1.0/ wl)});
            //            double minRu = max(0.0, min(rf,1.0));
            //            ru = min(max(minRu, ru), maxRu);
            //          }
          }
        }
      }
      break;
      //MUSCL
    case 10:
      {
        double du = phi_d - phi_c;

        if(fabs(du - 0.0) >= m_epsilon)
        {
          double rf = 0.0;

          //Darwish and Moukalled (2002)
          {
            rf = Vect::dotProduct(2.0 * grad_c.v[0], 2.0 * grad_c.v[1], 0.0, r_cdX, r_cdY, 0.0) / du;
            rf = rf - 1.0;
          }


          //Ubbink and Issa (1999)
          //          {
          //            double gradFaceX = (1.0 - idwfactor)* grad_c.v[0] + idwfactor * grad_d.v[0];
          //            double gradFaceY = (1.0 - idwfactor)* grad_c.v[1] + idwfactor * grad_d.v[1];
          //            double phi_u = phi_d - 2.0 * Vect::dotProduct(gradFaceX, gradFaceY, 0.0, dir,dir, 0.0);
          //            rf = (phi_c - phi_u)/(phi_d - phi_c);
          //          }

          ru = (rf  + fabs(rf))/(1 + fabs(rf));

          //Skewness correction FabianDenner∗, Berend G.M.vanWachem
          {
            //            ru = ru + Vect::dotProduct(gradFaceX / idwfactor, gradFaceY / idwfactor, 0.0, cv->df[faceIndex]) / du;

            //            //1st order
            //            {
            //              double maxRu = max(0.0, min(rf/idwfactor, 1.0/idwfactor));
            //              double minRu = 0.0;
            //              ru = min(max(minRu, ru), maxRu);
            //            }

            //2nd order
            //          {
            //            double maxRu = max({0.0, min(rf/idwfactor, 1.0), min(rf, 1.0/ wl)});
            //            double minRu = max(0.0, min(rf,1.0));
            //            ru = min(max(minRu, ru), maxRu);
            //          }
          }
        }
      }
      break;
      //FROMM
    case 11:
      {
        double du = phi_d - phi_c;

        if(fabs(du - 0.0) >= m_epsilon)
        {
          double rf = 0.0;

          //Darwish and Moukalled (2002)
          {
            rf = Vect::dotProduct(2.0 * grad_c.v[0], 2.0 * grad_c.v[1], 0.0, r_cdX, r_cdY, 0.0) / du;
            rf = rf - 1.0;
          }

          //Ubbink and Issa (1999)
          //          {
          //            double gradFaceX = (1.0 - idwfactor)* grad_c.v[0] + idwfactor * grad_d.v[0];
          //            double gradFaceY = (1.0 - idwfactor)* grad_c.v[1] + idwfactor * grad_d.v[1];
          //            double phi_u = phi_d - 2.0 * Vect::dotProduct(gradFaceX, gradFaceY, 0.0, dir,dir, 0.0);
          //            rf = (phi_c - phi_u)/(phi_d - phi_c);
          //          }

          ru = (1.0 + rf)/(2.0);

          //Skewness correction FabianDenner∗, Berend G.M.vanWachem
          {
            //            ru = ru + Vect::dotProduct(gradFaceX / idwfactor, gradFaceY / idwfactor, 0.0, cv->df[faceIndex]) / du;

            //            //1st order
            //            {
            //              double maxRu = max(0.0, min(rf/idwfactor, 1.0/idwfactor));
            //              double minRu = 0.0;
            //              ru = min(max(minRu, ru), maxRu);
            //            }

            //2nd order
            //          {
            //            double maxRu = max({0.0, min(rf/idwfactor, 1.0), min(rf, 1.0/ wl)});
            //            double minRu = max(0.0, min(rf,1.0));
            //            ru = min(max(minRu, ru), maxRu);
            //          }
          }
        }
      }
      break;
      //Sweby
    case 12:
      {
        double du = phi_d - phi_c;

        if(fabs(du - 0.0) >= m_epsilon)
        {
          double rf = 0.0;

          //Darwish and Moukalled (2002)
          {
            rf = Vect::dotProduct(2.0 * grad_c.v[0], 2.0 * grad_c.v[1], 0.0, r_cdX, r_cdY, 0.0) / du;
            rf = rf - 1.0;
          }

          //Ubbink and Issa (1999)
          //          {
          //            double gradFaceX = (1.0 - idwfactor)* grad_c.v[0] + idwfactor * grad_d.v[0];
          //            double gradFaceY = (1.0 - idwfactor)* grad_c.v[1] + idwfactor * grad_d.v[1];
          //            double phi_u = phi_d - 2.0 * Vect::dotProduct(gradFaceX, gradFaceY, 0.0, dir,dir, 0.0);
          //            rf = (phi_c - phi_u)/(phi_d - phi_c);
          //          }

          ru = (1.0 + rf)/(2.0);

          //          Vect &df = cv->df[faceIndex];
          //          ru +=  Vect::dotProduct(gradFace / wl, df) / (du);


          //          //1st order
          //          {
          //            double maxRu = max(0.0, min(rf/idwfactor, 1.0/idwfactor));
          //            double minRu = 0.0;
          //            ru = min(max(minRu, ru), maxRu);
          //          }

          //2nd order
          //          {
          //            double maxRu = max({0.0, min(rf/idwfactor, 1.0), min(rf, 1.0/ wl)});
          //            double minRu = max(0.0, min(rf,1.0));
          //            ru = min(max(minRu, ru), maxRu);
          //          }
        }
      }
      break;
    case 13:
      {
        ru = 0.3;
      }
      break;
  }

  return ru;
}

FVHMComponent::ErrorCode FVHMComponent::solvePressureCorrection(double tStep, double &pressureRelResidualNorm, double minorTermsCoeff , QString &errorMessage)
{

  ErrorCode error;
  double *x = new double[m_numCells];

  for(int iter = 0 ; iter < m_numPressureCorrectionIterations ; iter++)
  {

    //#ifdef USE_OPENMP
    //#pragma omp parallel for
    //#endif
    //    for(int rc = 0 ; rc < m_numCells ; rc++)
    //    {
    //      TriCV * cv = m_controlVolumes[rc];
    //      cv->calculateZCorrectionGradients();
    //    }

    double *h_b = new double[m_numCells];
    SparseMatrix coeffMatrix(m_numCells, 4);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(int index = 0 ; index < m_numCells ; index++)
    {
      x[index] = 0.0;

      TriCV *cv = m_controlVolumes[index];

      double constTerm  = 0;

      double h_p = cv->area / tStep;

      double sumEdgeCoeffs = 0;

      for(int e = 0; e < cv->numEdges; e++)
      {
        sumEdgeCoeffs += cv->velCoeffs[m_currentCoeff][e+1];
      }

      double coeffucv = cv->velCoeffs[m_currentCoeff][0] + minorTermsCoeff * sumEdgeCoeffs;

      double *coefficients = new double[cv->numEdges + 1]();
      coefficients[0] = 0.0;

      double zcorrest = 0.0;

      for(int i = 0; i < cv->numEdges; i++)
      {
        coefficients[i+1] = 0.0;
        int cvnIndex = cv->nTris[i];
        double edgeDepth = cv->faceDepths[i].value;
        double eta = cv->r_eta[i];

        TriCV *cvn = nullptr;

        double wl_n = cv->w_l[i];
        double wl_p = 1.0 - wl_n;

        if(cvnIndex > -1 && (cvn = m_controlVolumes[cvnIndex]))
        {
          sumEdgeCoeffs = 0;

          for(int e = 0; e < cvn->numEdges; e++)
          {
            sumEdgeCoeffs += cvn->velCoeffs[m_currentCoeff][e+1];
          }

          double coeffucvn = cvn->velCoeffs[m_currentCoeff][0] + minorTermsCoeff * sumEdgeCoeffs;

          double A_e = (wl_p * cv->area / coeffucv) + (wl_n * cvn->area / coeffucvn);

          double velCorrection = -g * A_e *  cv->faceDepths[i].value / cv->r_xi_dot_e_n[i];

          h_p += -1.0 * edgeDepth * eta * velCorrection * m_velocityRelaxFactor;
          double h_n = 1.0 * edgeDepth * eta * velCorrection * m_velocityRelaxFactor;

          h_p += wl_p * eta * cv->faceNormalVels[i].value;
          h_n += wl_n * eta * cv->faceNormalVels[i].value;

          //          constTerm -= eta * cv->faceNormalVels[i].value * Vect::dotProduct(cv->df[i], cv->grad_z_corr[0]);

          //          double facCV = wl_p * g * cv->h->value;
          //          double facCVN = wl_n * g * cvn->h->value;

          //          double tv2x = facCV * cv->grad_z_corr->v[0] + facCVN * cvn->grad_z_corr->v[0];
          //          double tv2y = facCV * cv->grad_z_corr->v[1] + facCVN * cvn->grad_z_corr->v[1];

          //          constTerm -= A_e * Vect::dotProduct(tv2x,tv2y, 0, cv->r_xi[i]) / cv->r_xi_dot_e_n[i];


          //          double cvgradZX = cv->h->value * cv->grad_z_corr->v[0];
          //          double cvgradZY = cv->h->value * cv->grad_z_corr->v[1];

          //          double cvngradZX = cvn->h->value * cvn->grad_z_corr->v[0];
          //          double cvngradZY = cvn->h->value * cvn->grad_z_corr->v[1];

          //          constTerm += g * A_e * (Vect::dotProduct(cvngradZX, cvngradZY, 0.0, cv->r_n_corr[i]) -
          //                                  Vect::dotProduct(cvgradZX , cvgradZY , 0.0, cv->r_p_corr[i])) / cv->r_xi_dot_e_n[i];

          coefficients[i+1] = h_n;
        }

        double infest = eta * cv->faceDepths[i].value * cv->faceNormalVels[i].value;
        constTerm -= infest;
        zcorrest += infest * tStep / cv->area;

      }

      coefficients[0] = h_p;

      for(int i = 0; i < cv->numEdges + 1; i++)
      {
        int cvnIndex = cv->orderedNTris[i];
        double coeff = 0.0;

        if(cvnIndex > -1 && (coeff = coefficients[cv->orderedNTrisIndexes[i]]))
        {
          coeffMatrix.appendValue(cv->index, m_controlVolumes[cvnIndex]->index, coeff);
        }
      }

      delete[] coefficients;


      //constTerm
      constTerm += cv->inflow;

      zcorrest += cv->inflow * tStep / cv->area;

      x[cv->index] = zcorrest;

      //tstep terms
      constTerm += (cv->prevH->value - cv->h->value) * cv->area  / tStep;

      //set const value
      h_b[cv->index] = constTerm;
    }

    double *residuals = new double[m_numCells];
    error = solve(coeffMatrix, h_b, x, residuals, pressureRelResidualNorm, errorMessage);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(int rc = 0 ; rc < m_numCells ; rc++)
    {
      TriCV * cv = m_controlVolumes[rc];
      double zCorr = x[cv->index];
      cv->zCorrection = zCorr;
    }

    delete[] residuals;
    delete[] h_b;
  }

  pressureRelResidualNorm = l2Norm(x,m_numCells);

  delete[] x;

  return error;
}

void FVHMComponent::applyPressureCorrections()
{

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int c = 0 ; c < m_numCells ; c++)
  {

    TriCV *cv = m_controlVolumes[c];

    if(cv->z->isBC == false)
    {
      double zCorr = cv->zCorrection;

      //      double zNew = cv->z->value + m_pressureRelaxFactor * zCorr;
      //      cv->setVFRWSE(zNew);

      double hNew = max(m_viscousSubLayerDepth, cv->h->value + m_pressureRelaxFactor * zCorr);
      cv->setVFRDepth(hNew);
    }
  }
}

void FVHMComponent::applyVelocityCorrections()
{

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int c = 0 ; c < m_numCells ; c++)
  {
    TriCV *cv = m_controlVolumes[c];
    cv->vel[0].value = m_velocityRelaxFactor * cv->vel[0].value  + (1.0 - m_velocityRelaxFactor) * cv->prevIterVel[0];
    cv->vel[1].value = m_velocityRelaxFactor * cv->vel[1].value  + (1.0 - m_velocityRelaxFactor) * cv->prevIterVel[1];
  }
}

void FVHMComponent::calculateContinuityResiduals(double tStep, double &contRelResidualNorm)
{
  double *residuals = new double[m_numCells];

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int rc = 0 ; rc < m_numCells ; rc++)
  {
    TriCV *cv = m_controlVolumes[rc];

    double sumflow = -cv->inflow * tStep;

    for(int i = 0; i < cv->numEdges ; i++)
    {
      sumflow += cv->faceDepths[i].value * cv->r_eta[i] * cv->faceNormalVels[i].value * tStep;
    }

    sumflow += (cv->h->value - cv->prevH->value) * cv->area;
    residuals[rc] = sumflow;
    cv->contResidualIter = sumflow;
  }

  contRelResidualNorm = l2Norm(residuals,m_numCells);

  delete[] residuals;
}

void FVHMComponent::updateExternalInflowOutflowTotals(double tStep)
{
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int rc = 0 ; rc < m_numCells ; rc++)
  {
    TriCV *cv = m_controlVolumes[rc];

    if(cv->inflow > 0)
    {
      cv->totalExternalInflow += cv->inflow * tStep;
    }
    else if(cv->inflow)
    {
      cv->totalExternalOutflow += -cv->inflow * tStep;
    }

    for(int f = 0 ; f < cv->numEdges; f++)
    {
      FaceNormVelBC &faceVel = cv->faceNormalVels[f];

      if(faceVel.isBC)
      {
        if(faceVel.associatedValue > 0.0)
        {
          cv->totalExternalOutflow += faceVel.associatedValue * tStep;
        }
        else
        {
          cv->totalExternalInflow += -faceVel.associatedValue * tStep;
        }
      }
    }
  }
}

FVHMComponent::ErrorCode FVHMComponent::solveConstituentEquations(int cIndex, double tStep, double &vvelRelResidualNorm, QString &errorMessage)
{
  errorMessage ="";
  vvelRelResidualNorm = cIndex * tStep;

  return ErrorCode::NoError;
}

void FVHMComponent::calculateFriction(int uv)
{

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0; i < m_numCells; i++)
  {
    TriCV *cv = m_controlVolumes[i];

    if(cv->wetIndex)
    {
      calculateFriction(cv, uv);
    }
    else
    {
      cv->friction[uv] = 0.0;

      for(int f = 0; f < 2; f++)
      {
        for(int j = 0 ; j < cv->numEdges ; j++)
        {
          cv->wallShearFriction[f][j] = 0.0;
        }
      }
    }
  }
}

void FVHMComponent::calculateFriction(TriCV *cv, int uv)
{
  double cf = g * cv->h->value * cv->mannings * cv->mannings / pow(cv->h->value, 4.0/3.0);
  double umag = hypot(cv->vel[0].value, cv->vel[1].value);
  double frictionVel = sqrt(cf * umag * umag);

  if(fabs(cv->vel[uv].value - 0.0) > m_epsilon)
  {
    cv->friction[uv] = cf * umag * cv->vel[uv].value *  cv->area;

    for(int i = 0 ; i < cv->numEdges ; i++)
    {
      int xy = (uv + 1) % 2;
      double relp = cv->r_e_l_p[i].v[xy];

      if(cv->faceNormalVels[i].calculateWallShearStress && relp > 0.0)
      {
        double yp = frictionVel * relp  / m_viscosity;
        double wallShear = frictionVel * 0.41 * cv->vel[uv].value * cv->faceDepths[i].value * cv->r_eta[i] / log(9.758 * yp);
        cv->wallShearFriction[uv][i] = wallShear;
      }
      else
      {
        cv->wallShearFriction[uv][i] = 0.0;
      }
    }
  }
  else
  {
    cv->friction[uv] = 0.0;

    for(int f = 0; f < 2 ; f++)
    {
      for(int i = 0 ; i < cv->numEdges; i++)
      {
        cv->wallShearFriction[f][i] = 0.0;
      }
    }
  }
}

void FVHMComponent::interpolateFaceVelocityFromGradient(TriCV *cv, int faceIndex)
{
  Vect edgeVel;
  Vect &rn =  cv->r_e[faceIndex];
  edgeVel.v[0] =  cv->vel[0].value + Vect::dotProduct(rn,cv->grad_vel[0]);
  edgeVel.v[1] =  cv->vel[1].value + Vect::dotProduct(rn,cv->grad_vel[1]);
  double ev = Vect::dotProduct(cv->e_n[faceIndex],edgeVel);
  cv->faceNormalVels[faceIndex].value = ev ;
}

void FVHMComponent::calculateCellVelocityGradients()
{

  double* residualsUX = new double[m_numCells];
  double* residualsUY = new double[m_numCells];
  double* residualsVX = new double[m_numCells];
  double* residualsVY = new double[m_numCells];

  double rux = 0.0;
  double ruy = 0.0;

  double rvx = 0.0;
  double rvy = 0.0;

  int iters = 0;


  do
  {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(int i = 0 ; i < m_numCells ; i++)
    {
      TriCV *cv = m_controlVolumes[i];

      double drux = cv->grad_vel[0].v[0];
      double druy = cv->grad_vel[0].v[1];

      double drvx = cv->grad_vel[1].v[0];
      double drvy = cv->grad_vel[1].v[1];

      cv->calculateVelocityGradient();

      residualsUX[i] = drux - cv->grad_vel[0].v[0];
      residualsUY[i] = druy - cv->grad_vel[0].v[1];

      residualsVX[i] = drvx - cv->grad_vel[1].v[0];
      residualsVY[i] = drvy - cv->grad_vel[1].v[1];

    }

    rux = l2Norm(residualsUX, m_numCells);
    ruy = l2Norm(residualsUY, m_numCells);

    rvx = l2Norm(residualsVX, m_numCells);
    rvy = l2Norm(residualsVY, m_numCells);

    iters ++;

  } while (rux > 1e-6 && ruy > 1e-6 && rvx > 1e-6 && rvy > 1e-6 && iters < 500);


#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0 ; i < m_numCells ; i++)
  {
    TriCV *cv = m_controlVolumes[i];
//    cv->interpolateNodeVelocities();
    cv->calculateNodeVelocities();
  }


  delete[] residualsUX;
  delete[] residualsUY;
  delete[] residualsVX;
  delete[] residualsVY;

}

void FVHMComponent::calculateCellEddyViscosities()
{
  if(m_useEddyViscosity)
  {
    int scheme = m_turbulenceScheme;

    switch (scheme)
    {
      case 1:
        {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
          for(int i = 0 ; i < m_numCells ; i++)
          {
            //Smargorinsky. Need to check to make sure correct.
            TriCV *cv = m_controlVolumes[i];
            double grad_ux_2 = cv->grad_vel[0].v[0] * cv->grad_vel[0].v[0];
            double grad_vy_2 = cv->grad_vel[1].v[1] * cv->grad_vel[1].v[1];
            double grad_uv_2 = cv->grad_vel[1].v[0] + cv->grad_vel[0].v[1];
            grad_uv_2 = grad_uv_2 * grad_uv_2;

            double mag = grad_ux_2 + grad_vy_2 + 0.5 * grad_uv_2;

            cv->eddyViscosity = mag ?  m_smargorinskyCoefficient * m_smargorinskyCoefficient * cv->area * sqrt(mag) : 0.0;
          }
        }
        break;
      default:
        {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
          // parabolic eddy viscosity
          for(int i = 0 ; i < m_numCells ; i++)
          {
            TriCV *cv = m_controlVolumes[i];
            double usquared = cv->vel[0].value * cv->vel[0].value + cv->vel[1].value * cv->vel[1].value;
            double frictionVel = sqrt(g * cv->h->value * cv->mannings * cv->mannings * usquared / pow(cv->h->value, 4.0/3.0));
            cv->eddyViscosity = m_parabolicEddyViscosityConstant * frictionVel * cv->h->value ;
          }
        }
        break;
    }


#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    // parabolic eddy viscosity
    for(int i = 0 ; i < m_numCells ; i++)
    {
      TriCV *cv = m_controlVolumes[i];
      cv->interpolateNodeEddyViscosities();
    }
  }
}

void FVHMComponent::calculateCellWSEGradients()
{

  double* residualsX = new double[m_numCells];
  double* residualsY = new double[m_numCells];

  double rx = 0.0;
  double ry = 0.0;

  int iters = 0;

  do
  {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(int i = 0 ; i < m_numCells ; i++)
    {
      TriCV *cv = m_controlVolumes[i];

      if(!cv->grad_z->isBC)
      {
        double gradDx = cv->grad_z->v[0];
        double gradDy = cv->grad_z->v[1];

        cv->calculateWSEGradient();

        residualsX[i] = cv->grad_z->v[0] - gradDx;
        residualsY[i] = cv->grad_z->v[1] - gradDy;
      }
      else
      {
        residualsX[i] = 0.0;
        residualsY[i] = 0.0;
      }
    }

    rx = l2Norm(residualsX, m_numCells);
    ry = l2Norm(residualsY, m_numCells);

    iters ++;

  }while ( rx > 1e-6 && ry > 1e-6 && iters < 500) ;


#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0 ; i < m_numCells ; i++)
  {
    TriCV *cv = m_controlVolumes[i];
    //    cv->interpolateNodeElevFromNeighbors();
    cv->calculateNodeElevations();
   cv->calculateEdgeDepths();
  }


  delete[] residualsX;
  delete[] residualsY;

}

FVHMComponent::ErrorCode FVHMComponent::solve(const SparseMatrix &A, const double b[], double x[], double residuals[], double &relativeResidualNorm, QString &errorMessage)
{
  ErrorCode solved = ErrorCode::NoError;

  switch (m_solverType)
  {
    //PCG with AMG
    case 0:
      {

        HYPRE_IJMatrix mA;
        HYPRE_IJVector mb;
        HYPRE_IJVector mx;

        int numRows = A.numRows();
        int *colsPerRow = A.colsPerRow();
        int *rowIndexes = A.rows();

        HYPRE_IJMatrixCreate(MPI_COMM_WORLD, 0, numRows - 1, 0, numRows - 1, &mA);
        HYPRE_IJMatrixSetObjectType(mA, HYPRE_PARCSR);
        HYPRE_IJMatrixInitialize(mA);

#ifdef USE_OPENMP
        HYPRE_IJMatrixSetOMPFlag(mA,1);
#endif

        HYPRE_IJVectorCreate(MPI_COMM_WORLD,0, numRows - 1,&mb);
        HYPRE_IJVectorSetObjectType(mb, HYPRE_PARCSR);
        HYPRE_IJVectorInitialize(mb);


        HYPRE_IJVectorCreate(MPI_COMM_WORLD,0,numRows - 1,&mx);
        HYPRE_IJVectorSetObjectType(mx, HYPRE_PARCSR);
        HYPRE_IJVectorInitialize(mx);

        int dataSize = A.getDataSize();
        double* values = new double[dataSize];
        int* colIndexes = new int[dataSize];

        A.getDataByRow(values);
        A.getColumnIndexes(colIndexes);

        HYPRE_IJMatrixSetValues(mA,numRows,colsPerRow, rowIndexes ,colIndexes,values);
        HYPRE_IJMatrixAssemble(mA);

        HYPRE_IJVectorSetValues(mb,numRows, rowIndexes,b);
        HYPRE_IJVectorAssemble(mb);

        HYPRE_IJVectorSetValues(mx,numRows, rowIndexes,x);
        HYPRE_IJVectorAssemble(mx);

        HYPRE_ParCSRMatrix parcsr_A;
        HYPRE_IJMatrixGetObject(mA, (void**) &parcsr_A);

        HYPRE_ParVector par_b;
        HYPRE_IJVectorGetObject(mb,(void**) &par_b);

        HYPRE_ParVector par_x;
        HYPRE_IJVectorGetObject(mx,(void**) &par_x);

        HYPRE_Solver solver;
        HYPRE_Solver precond;

        int num_iterations;

        HYPRE_ParVector par_residualVector;


        /* Create solver */
        HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);

        /* Set some parameters (See Reference Manual for more parameters) */
        HYPRE_PCGSetMaxIter(solver, m_solverMaximumNumberOfIterations); /* max iterations */
        HYPRE_PCGSetTol(solver, m_solverConvergenceTol); /* conv. tolerance */
        HYPRE_PCGSetTwoNorm(solver, 1); /* use the two norm as the stopping criteria */
        HYPRE_PCGSetPrintLevel(solver, 0); /* print solve info */
        HYPRE_PCGSetLogging(solver, 1); /* needed to get run info later */

        /* Now set up the AMG preconditioner and specify any parameters */
        HYPRE_BoomerAMGCreate(&precond);
        HYPRE_BoomerAMGSetPrintLevel(precond, 0); /* print amg solution info */
        HYPRE_BoomerAMGSetCoarsenType(precond, 6);

#ifdef USE_OPENMP
        HYPRE_BoomerAMGSetRelaxType(precond, 5); /* Sym G.S./Jacobi hybrid */
#else
        HYPRE_BoomerAMGSetRelaxType(precond, 6);
#endif

        HYPRE_BoomerAMGSetNumSweeps(precond, 1);
        HYPRE_BoomerAMGSetMaxIter(precond, m_solverAMGPreconditionerNumberOfIterations); /* do only one iteration! */
        HYPRE_BoomerAMGSetTol(precond, 0.0); /* conv. tolerance zero */

        /* Set the FlexGMRES preconditioner */
        HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
                            (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);

        /* Now setup and solve! */
        HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_b, par_x);

        HYPRE_Int result = HYPRE_ParCSRPCGSolve(solver, parcsr_A, par_b, par_x);

        /*Now copy solution! */
        //    HYPRE_IJVectorPrint(mx,"/Users/calebbuahin/Documents/Projects/HydroCouple/FVHMComponent/examples/XVector.txt");
        HYPRE_IJVectorGetValues(mx,numRows,rowIndexes,x);

        /*Copy residuals! */
        HYPRE_FlexGMRESGetResidual(solver, (void**)&par_residualVector);
        double* residualValues = hypre_VectorData(hypre_ParVectorLocalVector(par_residualVector));

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
        for(int i = 0; i < numRows; i++)
        {
          residuals[i] = residualValues[i];
        }

        /* Run info - needed logging turned on */
        HYPRE_PCGGetNumIterations(solver, &num_iterations);
        HYPRE_PCGGetFinalRelativeResidualNorm(solver, &relativeResidualNorm);

        switch (result)
        {
          case 0:
            {
              errorMessage ="";
            }
            break;
          case 1:
          case 2:
          case 3:
            {
              if(m_outputDir.exists())
              {
                HYPRE_IJMatrixPrint(mA, qPrintable(m_outputDir.absolutePath() + "/AMatrix.txt"));
                HYPRE_IJVectorPrint(mb, qPrintable(m_outputDir.absolutePath() + "/BVector.txt"));
                HYPRE_IJVectorPrint(mx, qPrintable(m_outputDir.absolutePath() + "/XVector.txt"));
              }

              char message[100];
              HYPRE_DescribeError(result,message);
              errorMessage = QString(message);
              solved = ErrorCode::CriticalFailure;
            }
            break;
          case 256:
            {
              if(m_outputDir.exists())
              {
                HYPRE_IJMatrixPrint(mA, qPrintable(m_outputDir.absolutePath() + "/AMatrix.txt"));
                HYPRE_IJVectorPrint(mb, qPrintable(m_outputDir.absolutePath() + "/BVector.txt"));
                HYPRE_IJVectorPrint(mx, qPrintable(m_outputDir.absolutePath() + "/XVector.txt"));
              }

              char message[100];
              HYPRE_DescribeError(result,message);
              errorMessage = QString(message);
              solved = ErrorCode::SolverFailedToConverge;
            }
            break;
        }

        /* Destroy solver and preconditioner */
        HYPRE_ParCSRPCGDestroy(solver);
        HYPRE_BoomerAMGDestroy(precond);

        HYPRE_IJMatrixDestroy(mA);
        HYPRE_IJVectorDestroy(mb);
        HYPRE_IJVectorDestroy(mx);

        delete[] values;
        delete[] colIndexes;
      }
      break;
      //GMRES and AMG
    case 1:
      {
        HYPRE_IJMatrix mA;
        HYPRE_IJVector mb;
        HYPRE_IJVector mx;

        int numRows = A.numRows();
        int *colsPerRow = A.colsPerRow();
        int *rowIndexes = A.rows();

        HYPRE_IJMatrixCreate(MPI_COMM_WORLD, 0, numRows - 1, 0, numRows - 1, &mA);
        HYPRE_IJMatrixSetObjectType(mA, HYPRE_PARCSR);
        HYPRE_IJMatrixInitialize(mA);

#ifdef USE_OPENMP
        HYPRE_IJMatrixSetOMPFlag(mA,1);
#endif

        HYPRE_IJVectorCreate(MPI_COMM_WORLD,0, numRows - 1,&mb);
        HYPRE_IJVectorSetObjectType(mb, HYPRE_PARCSR);
        HYPRE_IJVectorInitialize(mb);


        HYPRE_IJVectorCreate(MPI_COMM_WORLD,0,numRows - 1,&mx);
        HYPRE_IJVectorSetObjectType(mx, HYPRE_PARCSR);
        HYPRE_IJVectorInitialize(mx);

        int dataSize = A.getDataSize();
        double* values = new double[dataSize];
        int* colIndexes = new int[dataSize];

        A.getDataByRow(values);
        A.getColumnIndexes(colIndexes);

        HYPRE_IJMatrixSetValues(mA,numRows,colsPerRow, rowIndexes ,colIndexes,values);
        HYPRE_IJMatrixAssemble(mA);

        HYPRE_IJVectorSetValues(mb,numRows, rowIndexes,b);
        HYPRE_IJVectorAssemble(mb);

        HYPRE_IJVectorSetValues(mx,numRows, rowIndexes,x);
        HYPRE_IJVectorAssemble(mx);

        HYPRE_ParCSRMatrix parcsr_A;
        HYPRE_IJMatrixGetObject(mA, (void**) &parcsr_A);

        HYPRE_ParVector par_b;
        HYPRE_IJVectorGetObject(mb,(void**) &par_b);

        HYPRE_ParVector par_x;
        HYPRE_IJVectorGetObject(mx,(void**) &par_x);

        HYPRE_Solver solver;
        HYPRE_Solver precond;

        int num_iterations;

        HYPRE_ParVector par_residualVector;

        HYPRE_ParCSRFlexGMRESCreate(MPI_COMM_WORLD, &solver);

        /* Set some parameters (See Reference Manual for more parameters) */
        //HYPRE_FlexGMRESSetKDim(solver, 30);
        HYPRE_FlexGMRESSetMaxIter(solver, m_solverMaximumNumberOfIterations); /* max iterations */
        HYPRE_FlexGMRESSetTol(solver, m_solverConvergenceTol); /* conv. tolerance */
        HYPRE_FlexGMRESSetPrintLevel(solver, 0); /* print solve info */
        HYPRE_FlexGMRESSetLogging(solver, 1); /* needed to get run info later */

        /* Now set up the AMG preconditioner and specify any parameters */
        HYPRE_BoomerAMGCreate(&precond);
        HYPRE_BoomerAMGSetPrintLevel(precond, 0); /* print amg solution info */
        HYPRE_BoomerAMGSetCoarsenType(precond, 6);

#ifdef USE_OPENMP
        HYPRE_BoomerAMGSetRelaxType(precond, 5);
#else
        HYPRE_BoomerAMGSetRelaxType(precond, 6);/* Sym G.S./Jacobi hybrid */
#endif

        HYPRE_BoomerAMGSetNumSweeps(precond, 1);
        HYPRE_BoomerAMGSetMaxIter(precond, m_solverAMGPreconditionerNumberOfIterations); /* do only one iteration! */
        HYPRE_BoomerAMGSetTol(precond, 0.0); /* conv. tolerance zero */

        /* Set the FlexGMRES preconditioner */
        HYPRE_FlexGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
                                  (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);

        /* Now setup and solve! */
        HYPRE_ParCSRFlexGMRESSetup(solver, parcsr_A, par_b, par_x);

        HYPRE_Int result = HYPRE_ParCSRFlexGMRESSolve(solver, parcsr_A, par_b, par_x);


        /*Now copy solution! */
        HYPRE_IJVectorGetValues(mx,numRows,rowIndexes,x);

        /*Copy residuals! */
        HYPRE_FlexGMRESGetResidual(solver, (void**)&par_residualVector);
        double* residualValues = hypre_VectorData(hypre_ParVectorLocalVector(par_residualVector));

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
        for(int i = 0; i < numRows; i++)
        {
          residuals[i] = residualValues[i];
        }

        /* Run info - needed logging turned on */
        HYPRE_FlexGMRESGetNumIterations(solver, &num_iterations);
        HYPRE_FlexGMRESGetFinalRelativeResidualNorm(solver, &relativeResidualNorm);

        switch (result)
        {
          case 0:
            {
              errorMessage ="";
            }
            break;
          case 1:
          case 2:
          case 3:
            {

              if(m_outputDir.exists())
              {
                HYPRE_IJMatrixPrint(mA, qPrintable(m_outputDir.absolutePath() + "/AMatrix.txt"));
                HYPRE_IJVectorPrint(mb, qPrintable(m_outputDir.absolutePath() + "/BVector.txt"));
                HYPRE_IJVectorPrint(mx, qPrintable(m_outputDir.absolutePath() + "/XVector.txt"));
              }

              char message[100];
              HYPRE_DescribeError(result,message);
              errorMessage = QString(message);
              solved = ErrorCode::CriticalFailure;
            }
            break;
          case 256:
            {
              if(m_outputDir.exists())
              {
                HYPRE_IJMatrixPrint(mA, qPrintable(m_outputDir.absolutePath() + "/AMatrix.txt"));
                HYPRE_IJVectorPrint(mb, qPrintable(m_outputDir.absolutePath() + "/BVector.txt"));
                HYPRE_IJVectorPrint(mx, qPrintable(m_outputDir.absolutePath() + "/XVector.txt"));
              }

              char message[100];
              HYPRE_DescribeError(result,message);
              errorMessage = QString(message);
              solved = ErrorCode::SolverFailedToConverge;
            }
            break;
        }

        /* Destroy solver and preconditioner */
        HYPRE_ParCSRFlexGMRESDestroy(solver);
        HYPRE_BoomerAMGDestroy(precond);

        HYPRE_IJMatrixDestroy(mA);
        HYPRE_IJVectorDestroy(mb);
        HYPRE_IJVectorDestroy(mx);

        delete[] values;
        delete[] colIndexes;
      }
      break;

#ifdef USE_SUITESPARSE
      //SUITE SPARSE
    case 2:
      {
        //        cholmod_common Common, *cc ;
        //        cholmod_sparse *AA ;
        //        cholmod_dense *XX, *BB, *Residual ;

        //        int dataSize = A.getDataSize();
        //        int numRows = A.numRows();

        //        double rnorm, one [2] = {1,0}, minusone [2] = {-1,0};

        //        // start CHOLMOD
        //        cc = &Common ;
        //        cholmod_l_start (cc) ;

        //        // load A
        //        AA = cholmod_allocate_sparse(numRows, numRows, dataSize, FALSE, FALSE, 0, CHOLMOD_REAL, cc);
        //        int *cpoints = AA->p;
        //        int *rowIndexes = AA->i;
        //        int *nzLoc = AA->nz;
        //        double* values = AA->x;

        //        A.getRowIndexes(rowIndexes,cpoints,nzLoc);
        //        A.getDataByColumn(values);


        //        BB = cholmod_allocate_dense(numRows, 1, numRows, CHOLMOD_REAL, cc);

        //#ifdef USE_OPENMP
        //#pragma omp parallel for
        //#endif
        //        for(int i = 0; i < numRows; i++)
        //        {
        //          ((double*)BB->x)[i] = b[i];
        //        }


        //        // X = A\B
        //        XX = SuiteSparseQR<double>(AA, BB, cc) ;


        //        // rnorm = norm (B-A*X)
        //        Residual = cholmod_l_copy_dense (BB, cc) ;
        //        cholmod_l_sdmult (AA, 0, minusone, one, XX, Residual, cc);

        //        double *castRes = (double*) Residual->x;
        //        double *castX = (double*) XX->x;

        //#ifdef USE_OPENMP
        //#pragma omp parallel for
        //#endif
        //        for(int i = 0; i < numRows; i++)
        //        {
        //          residuals[i] = castRes[i];
        //          x[i] = castX[i];
        //        }


        //        rnorm = cholmod_l_norm_dense (Residual, 2, cc);
        //        relativeResidualNorm = rnorm;

        //        printf ("2-norm of residual: %8.1e\n", rnorm) ;
        //        printf ("rank %ld\n", cc->SPQR_istat [4]);
        //        // free everything and finish CHOLMOD
        //        cholmod_l_free_dense (&Residual, cc) ;
        //        cholmod_l_free_sparse (&AA, cc) ;
        //        cholmod_l_free_dense (&XX, cc) ;
        //        cholmod_l_free_dense (&BB, cc) ;
      }
      break;
      //KLU
    case 3:
      {

      }
      break;
      //CHOLMOD
    case 4:
      {


      }
      break;
      //UMFPACK
    case 5:
      {


      }
      break;
#endif
      //GMRES and AMG
    default:
      {
        HYPRE_IJMatrix mA;
        HYPRE_IJVector mb;
        HYPRE_IJVector mx;

        int numRows = A.numRows();
        int *colsPerRow = A.colsPerRow();
        int *rowIndexes = A.rows();

        HYPRE_IJMatrixCreate(MPI_COMM_WORLD, 0, numRows - 1, 0, numRows - 1, &mA);
        HYPRE_IJMatrixSetObjectType(mA, HYPRE_PARCSR);
        HYPRE_IJMatrixInitialize(mA);

#ifdef USE_OPENMP
        HYPRE_IJMatrixSetOMPFlag(mA,1);
#endif

        HYPRE_IJVectorCreate(MPI_COMM_WORLD,0, numRows - 1,&mb);
        HYPRE_IJVectorSetObjectType(mb, HYPRE_PARCSR);
        HYPRE_IJVectorInitialize(mb);


        HYPRE_IJVectorCreate(MPI_COMM_WORLD,0,numRows - 1,&mx);
        HYPRE_IJVectorSetObjectType(mx, HYPRE_PARCSR);
        HYPRE_IJVectorInitialize(mx);

        int dataSize = A.getDataSize();
        double* values = new double[dataSize];
        int* colIndexes = new int[dataSize];

        A.getDataByRow(values);
        A.getColumnIndexes(colIndexes);

        HYPRE_IJMatrixSetValues(mA,numRows,colsPerRow, rowIndexes ,colIndexes,values);
        HYPRE_IJMatrixAssemble(mA);

        HYPRE_IJVectorSetValues(mb,numRows, rowIndexes,b);
        HYPRE_IJVectorAssemble(mb);

        HYPRE_IJVectorSetValues(mx,numRows, rowIndexes,x);
        HYPRE_IJVectorAssemble(mx);

        HYPRE_ParCSRMatrix parcsr_A;
        HYPRE_IJMatrixGetObject(mA, (void**) &parcsr_A);

        HYPRE_ParVector par_b;
        HYPRE_IJVectorGetObject(mb,(void**) &par_b);

        HYPRE_ParVector par_x;
        HYPRE_IJVectorGetObject(mx,(void**) &par_x);

        HYPRE_Solver solver;
        HYPRE_Solver precond;

        int num_iterations;

        HYPRE_ParVector par_residualVector;

        HYPRE_ParCSRFlexGMRESCreate(MPI_COMM_WORLD, &solver);

        /* Set some parameters (See Reference Manual for more parameters) */
        HYPRE_FlexGMRESSetMaxIter(solver, m_solverMaximumNumberOfIterations); /* max iterations */
        HYPRE_FlexGMRESSetTol(solver, m_solverConvergenceTol); /* conv. tolerance */
        HYPRE_FlexGMRESSetPrintLevel(solver, 0); /* print solve info */
        HYPRE_FlexGMRESSetLogging(solver, 1); /* needed to get run info later */

        /* Now set up the AMG preconditioner and specify any parameters */
        HYPRE_BoomerAMGCreate(&precond);
        HYPRE_BoomerAMGSetPrintLevel(precond, 0); /* print amg solution info */
        HYPRE_BoomerAMGSetCoarsenType(precond, 6);

#ifdef USE_OPENMP
        HYPRE_BoomerAMGSetRelaxType(precond, 5); /* Sym G.S./Jacobi hybrid */
#else
        HYPRE_BoomerAMGSetRelaxType(precond, 6);
#endif

        HYPRE_BoomerAMGSetNumSweeps(precond, 1);
        HYPRE_BoomerAMGSetMaxIter(precond, m_solverAMGPreconditionerNumberOfIterations); /* do only one iteration! */
        HYPRE_BoomerAMGSetTol(precond, 0.0); /* conv. tolerance zero */

        /* Set the FlexGMRES preconditioner */
        HYPRE_FlexGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
                                  (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);

        /* Now setup and solve! */
        HYPRE_ParCSRFlexGMRESSetup(solver, parcsr_A, par_b, par_x);

        HYPRE_Int result = HYPRE_ParCSRFlexGMRESSolve(solver, parcsr_A, par_b, par_x);


        /*Now copy solution! */
        //    HYPRE_IJVectorPrint(mx,"/Users/calebbuahin/Documents/Projects/HydroCouple/FVHMComponent/examples/basic_tests/sloping/XVector.txt");
        HYPRE_IJVectorGetValues(mx,numRows,rowIndexes,x);

        /*Copy residuals! */
        HYPRE_FlexGMRESGetResidual(solver, (void**)&par_residualVector);
        double* residualValues = hypre_VectorData(hypre_ParVectorLocalVector(par_residualVector));

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
        for(int i = 0; i < numRows; i++)
        {
          residuals[i] = residualValues[i];
        }

        /* Run info - needed logging turned on */
        HYPRE_FlexGMRESGetNumIterations(solver, &num_iterations);
        HYPRE_FlexGMRESGetFinalRelativeResidualNorm(solver, &relativeResidualNorm);

        switch (result)
        {
          case 0:
            {
              errorMessage ="";
            }
            break;
          case 1:
          case 2:
          case 3:
            {
              if(m_outputDir.exists())
              {
                HYPRE_IJMatrixPrint(mA, qPrintable(m_outputDir.absolutePath() + "/AMatrix.txt"));
                HYPRE_IJVectorPrint(mb, qPrintable(m_outputDir.absolutePath() + "/BVector.txt"));
                HYPRE_IJVectorPrint(mx, qPrintable(m_outputDir.absolutePath() + "/XVector.txt"));
              }

              char message[100];
              HYPRE_DescribeError(result,message);
              errorMessage = QString(message);
              solved = ErrorCode::CriticalFailure;
            }
            break;
          case 256:
            {
              if(m_outputDir.exists())
              {
                HYPRE_IJMatrixPrint(mA, qPrintable(m_outputDir.absolutePath() + "/AMatrix.txt"));
                HYPRE_IJVectorPrint(mb, qPrintable(m_outputDir.absolutePath() + "/BVector.txt"));
                HYPRE_IJVectorPrint(mx, qPrintable(m_outputDir.absolutePath() + "/XVector.txt"));
              }

              char message[100];
              HYPRE_DescribeError(result,message);
              errorMessage = QString(message);
              solved = ErrorCode::SolverFailedToConverge;
            }
            break;
        }

        /* Destroy solver and preconditioner */
        HYPRE_ParCSRFlexGMRESDestroy(solver);
        HYPRE_BoomerAMGDestroy(precond);

        HYPRE_IJMatrixDestroy(mA);
        HYPRE_IJVectorDestroy(mb);
        HYPRE_IJVectorDestroy(mx);

        delete[] values;
        delete[] colIndexes;
      }
      break;
  }



  return solved;
}

double FVHMComponent::sign(double value)
{
  return (value >= 0) ? 1 : ((value < 0) ? -1 : 0);
}

double FVHMComponent::l2Norm(const double residuals[], int length)
{
  double value = 0;

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0; i < length ; i++)
  {
    double v = residuals[i];
    v *=v;
#ifdef USE_OPENMP
#pragma omp critical
#endif
    value += v;
  }

  return value > 0.0 ? sqrt(value) : 0.0;
}

double FVHMComponent::relativeResidualNorm(const double residuals[], const double values[], int length)
{
  double denom = 0;
  double numer = 0;

  for(int i = 0 ; i < length; i++)
  {
    double r = residuals[i];
    numer += r * r;
    denom += values[i];
  }

  return sqrt(length * numer) / denom;
}

bool FVHMComponent::isWet(const TriCV *cv)
{
  return cv->h->value >= m_wetCellDepth;
}

std::tuple<double, double, double> FVHMComponent::calculateZeroGradientFaceVelocity(TriCV *cv, int face)
{
  Vect &unorm = cv->e_n[face];
  Vect &rdist = cv->r_e[face];

  Vect gradUEdge = (cv->grad_vel[0]) - Vect::dotProduct(cv->grad_vel[0], unorm) * unorm;
  Vect gradVEdge = (cv->grad_vel[1]) - Vect::dotProduct(cv->grad_vel[1], unorm) * unorm;

  Vect v1;
  v1.v[0]  = cv->vel[0].value + Vect::dotProduct(gradUEdge, rdist);
  v1.v[1]  = cv->vel[1].value + Vect::dotProduct(gradVEdge, rdist);

  double outvel = Vect::dotProduct(v1,unorm);

  return std::make_tuple(v1.v[0] , v1.v[1] , outvel);
}
