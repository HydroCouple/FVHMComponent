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
          for(int i = 0 ; i < m_numContCells  ; i++)
          {
            int cvIndex = m_contCells[i];
            TriCV *cv = m_controlVolumes[cvIndex];

            double cFactor = cv->getCourantFactor();

            if(cFactor > maxCFactor)
            {

#ifdef USE_OPENMP
#pragma omp atomic read
#endif
              maxCFactor = cFactor;
            }
          }

          estimatedTimeStep = maxCFactor ?  m_maxCourantNumber * m_timeStepRelaxFactor / maxCFactor : m_maxTimeStep;
        }
        break;
      case AdaptiveTSMode::RMSCourantNumber:
        {
          double count = 0;
          double tStep = 0;

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
          for(int i = 0 ; i < m_numWetCells ; i++)
          {
            int cvIndex = m_wetCells[i];
            TriCV *cv = m_controlVolumes[cvIndex];

            double cfactor = cv->getCourantFactor();

#ifdef USE_OPENMP
#pragma omp atomic
#endif
            tStep += cfactor * cfactor;

#ifdef USE_OPENMP
#pragma omp atomic
#endif
            count ++;

          }

          if(count && tStep)
          {
            estimatedTimeStep =  m_RMSCourantNumber * m_timeStepRelaxFactor / sqrt(tStep/count) ;
          }
        }
        break;
    }

    double dtStep = estimatedTimeStep - m_timeStep;
    double fract = dtStep / m_timeStep;

    if(dtStep > 0)
    {
      if(fabs(fract) > m_maxTimeStepIncreaseCF)
      {
        estimatedTimeStep = m_timeStep + m_maxTimeStepIncreaseCF * fract * m_timeStep / fabs(fract);
      }
    }
    else
    {
      estimatedTimeStep = estimatedTimeStep - estimatedTimeStep * m_maxTimeStepDecreaseCF;
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
    cv->inflowOutflow = 0.0;
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
    cv->inflowOutflow = 0.0;
    cv->externalForce->v[0] = 0.0;
    cv->externalForce->v[1] = 0.0;
  }

#ifdef USE_OPENMP
#pragma omp parallel sections
#endif
  {

#ifdef USE_OPENMP
#pragma omp section
#endif
    {
      m_outletWSEBC->applyBoundaryConditions(time, m_timeStep);

      m_outletWSESlope->applyBoundaryConditions(time, m_timeStep);

      m_criticalDepthOutflowBC->applyBoundaryConditions(time, m_timeStep);
    }

#ifdef USE_OPENMP
#pragma omp section
#endif
    {
      m_inletWSEBC->applyBoundaryConditions(time, m_timeStep);
    }

#ifdef USE_OPENMP
#pragma omp section
#endif
    {
      m_outletFlowBCArgument->applyBoundaryConditions(time, m_timeStep);
    }

#ifdef USE_OPENMP
#pragma omp section
#endif
    {
      m_inletFlowBCArgument->applyBoundaryConditions(time, m_timeStep);
    }

#ifdef USE_OPENMP
#pragma omp section
#endif
    {
      m_precipitationArgument->applyBoundaryConditions(time, m_timeStep);
    }
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

void FVHMComponent::prepareForNextTimeStep(double &minU, double &maxU, double &minV, double &maxV, double &minH, double &maxH, double &minZ, double &maxZ)
{

  maxV = std::numeric_limits<double>::lowest();
  minV = std::numeric_limits<double>::max();

  maxU = std::numeric_limits<double>::lowest();
  minU = std::numeric_limits<double>::max();

  maxH = std::numeric_limits<double>::lowest();
  minH = std::numeric_limits<double>::max();

  maxZ = std::numeric_limits<double>::lowest();
  minZ = std::numeric_limits<double>::max();

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0 ; i < m_numCells ; i++)
  {
    TriCV *cv = m_controlVolumes[i];
    cv->contResidualTotal += cv->contResidualIter;
    cv->velResidual[0] = cv->velResidualIter[0];
    cv->velResidual[1] = cv->velResidualIter[1];

    maxU = max(maxU , cv->vel[0].value);
    maxV = max(maxV , cv->vel[1].value);

    minU = min(minU , cv->vel[0].value);
    minV = min(minV , cv->vel[1].value);

    maxH = max(maxH, cv->h->value);
    minH = min(minH, cv->h->value);

    maxZ = max(maxZ, cv->z->value);
    minZ = min(minZ, cv->z->value);

    cv->copyVariablesToPrev();
  }
}

void FVHMComponent::setWetAndContCells()
{
  m_numWetCells = 0;
  m_numContCells = 0;

  for(int i = 0; i < m_numCells; i++)
  {
    TriCV* cv = m_controlVolumes[i];
    double coeff = cv->h->value * cv->area / m_timeStep;
    cv->velCoeffs[0][0] = cv->velCoeffs[1][0] = coeff;
    cv->wetCellIndex = -1;
    cv->contCellIndex = -1;
    cv->contResidualIter = cv->velResidualIter[0] = cv->velResidualIter[1] =
        cv->isWetOrHasWetNeigh = cv->wetIndex = cv->contIndex = 0.0;

    if(cv->h->value > m_wetCellDepth)
    {
      cv->isWetOrHasWetNeigh = 1;
      cv->wetIndex = 1;
      cv->contIndex = 1;
      cv->wetCellIndex = m_numWetCells; m_wetCells[m_numWetCells] = i; m_numWetCells++;
      cv->contCellIndex = m_numContCells; m_contCells[m_numContCells] = i; m_numContCells++;
    }
    else
    {
      for(int c = 0 ; c < 2; c++)
      {
        cv->prevIterVel[c] = 0.0;
        cv->velResidualIter[c] = 0.0;
        cv->vel[c].value = 0.0;
      }

      bool cont = false;

      if(cv->inflowOutflow)
      {
        cont = true;
      }

      for(int j = 0; j < cv->numEdges; j++)
      {

        int cvnIndex = cv->nTris[j];

        if(cvnIndex > -1)
        {
          TriCV *cvn = m_controlVolumes[cvnIndex];

          if(cvn->h->value > m_wetCellDepth)
          {
            cv->isWetOrHasWetNeigh = 2;
            cont = true;
          }
        }

        if(cv->faceNormalVels[j].isBC && cv->faceNormalVels[j].value)
        {
          cont = true;
        }
      }

      if(cont)
      {
        cv->contIndex = 1;
        cv->contCellIndex = m_numContCells; m_contCells[m_numContCells] = i; m_numContCells++;
      }
      else
      {
        cv->grad_h->zero();
        cv->grad_z->zero();
        cv->grad_vel[0].zero();
        cv->grad_vel[1].zero();
        cv->zCorrection = 0.0;
        cv->grad_zcorr->zero();
      }
    }
  }
}

bool FVHMComponent::neighbourIsWet(TriCV *cv, double wetCellDepth)
{
  for(int i = 0; i < cv->numEdges; i++)
  {
    int cvIndex = cv->nTris[i];

    if(cvIndex > -1)
    {
      TriCV *cvn = m_controlVolumes[cvIndex];

      if(cvn->h->value > wetCellDepth)
      {
        return true;
      }
    }
  }

  return false;
}

int FVHMComponent::performSimpleTimeStep(double tStep, int &numIterations,
                                         double &uvelRelResidualNormInit, double &uvelRelResidualNormFin, int &minUVelSolvIters, int &maxUVelSolvIters,
                                         double &vvelRelResidualNormInit, double &vvelRelResidualNormFin, int &minVVelSolvIters, int &maxVVelSolvIters,
                                         double &pressureRelRisdualNormInit, double &pressureRelRisdualNormFin, int &minPressSolvIters,int &maxPressSolvIters,
                                         double &continuityRelResidualNormInit, double &continuityRelResidualNormFin,
                                         bool &converged, QString &errorMessage)
{

  ErrorCode error = ErrorCode::NoError;

  numIterations = 0;

  uvelRelResidualNormFin = vvelRelResidualNormFin = pressureRelRisdualNormFin =
      continuityRelResidualNormFin = uvelRelResidualNormInit = vvelRelResidualNormInit =
      pressureRelRisdualNormInit = continuityRelResidualNormInit = 0.0;

  converged = false;

  double fractStep = tStep / m_numFractionalSteps;
  double currentStep = fractStep;

  m_currentCoeff = m_currentCoeff == 0 ? 0 : 0;

  //Friction
  {
#ifdef USE_OPENMP
#pragma omp parallel sections
#endif
    {

#ifdef USE_OPENMP
#pragma omp section
#endif
      {
        calculateFriction(0);
      }

#ifdef USE_OPENMP
#pragma omp section
#endif
      {
        calculateFriction(1);
      }
    }
  }

  while (numIterations < m_itersPerTimeStep)
  {
    ErrorCode uerror = ErrorCode::NoError;
    ErrorCode verror = ErrorCode::NoError;
    QString uErrorMessage;
    QString vErrorMessage;
    uvelRelResidualNormFin = vvelRelResidualNormFin =
        pressureRelRisdualNormFin = continuityRelResidualNormFin = 0.0;

    int numUVelSolvIters = 0;
    int numVVelSolvIters = 0;
    int numPSolvIters  = 0;

    currentStep += fractStep;
    currentStep = min(currentStep, tStep);

    if (m_numWetCells)
    {
      //Solve Momentum
      {
#ifdef USE_OPENMP
#pragma omp parallel sections
#endif
        {

#ifdef USE_OPENMP
#pragma omp section
#endif
          {
            uerror = solveMomentumEquations(currentStep, 0, uvelRelResidualNormFin, numUVelSolvIters, uErrorMessage);
          }

#ifdef USE_OPENMP
#pragma omp section
#endif
          {
            verror = solveMomentumEquations(currentStep, 1, vvelRelResidualNormFin, numVVelSolvIters, vErrorMessage);
          }
        }

        if (uerror == ErrorCode::CriticalFailure)
        {
          error = uerror;
          errorMessage = uErrorMessage;
          break;
        }

        if (verror == ErrorCode::CriticalFailure)
        {
          error = verror;
          errorMessage = vErrorMessage;
          break;
        }
      }
    }

    calculateFaceVelocities();

    error = solvePressureCorrection(currentStep, pressureRelRisdualNormFin, 0.0, numPSolvIters, errorMessage);

    minUVelSolvIters = min(minUVelSolvIters, numUVelSolvIters);
    minVVelSolvIters = min(minVVelSolvIters, numVVelSolvIters);
    minPressSolvIters = min(minPressSolvIters,numPSolvIters);

    maxUVelSolvIters = max(maxUVelSolvIters, numUVelSolvIters);
    maxVVelSolvIters = max(maxVVelSolvIters, numVVelSolvIters);
    maxPressSolvIters = max(maxPressSolvIters,numPSolvIters);

    if (error == ErrorCode::CriticalFailure)
    {
      break;
    }

    //Velocity & Pressure Correction
    {
#ifdef USE_OPENMP
#pragma omp parallel sections
#endif
      {

#ifdef USE_OPENMP
#pragma omp section
#endif
        {
          applyPressureCorrections();
          calculateCellWSEGradients();
        }

#ifdef USE_OPENMP
#pragma omp section
#endif
        {
          applyVelocityCorrections();
          calculateCellVelocityGradients();
          calculateCellEddyViscosities();
        }
      }
    }

    calculateContinuityResiduals(currentStep, continuityRelResidualNormFin);

    if (numIterations == 0)
    {
      uvelRelResidualNormInit = uvelRelResidualNormFin;
      vvelRelResidualNormInit = vvelRelResidualNormFin;
      pressureRelRisdualNormInit = pressureRelRisdualNormFin;
      continuityRelResidualNormInit = continuityRelResidualNormFin;
    }

    //check for convergence
    {
      if (uvelRelResidualNormFin < m_uvelConvergenceTol &&
          vvelRelResidualNormFin < m_vvelConvergenceTol &&
          fabs(pressureRelRisdualNormFin) < m_pressureConvergenceTol &&
          continuityRelResidualNormFin < m_contConvergenceTol &&
          currentStep == tStep)
      {
        numIterations++;
        converged = true;
        break;
      }
    }

    numIterations++;
  }

  if (error == ErrorCode::NoError)
  {
    if (!converged)
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


int FVHMComponent::performSimpleCTimeStep(double tStep, int &numIterations,
                                          double &uvelRelResidualNormInit, double &uvelRelResidualNormFin, int &minUVelSolvIters, int &maxUVelSolvIters,
                                          double &vvelRelResidualNormInit, double &vvelRelResidualNormFin, int &minVVelSolvIters, int &maxVVelSolvIters,
                                          double &pressureRelRisdualNormInit, double &pressureRelRisdualNormFin, int &minPressSolvIters, int &maxPressSolvIters,
                                          double &continuityRelResidualNormInit, double &continuityRelResidualNormFin,
                                          bool &converged, QString &errorMessage)
{

  ErrorCode error = ErrorCode::NoError;

  numIterations = 0;

  uvelRelResidualNormFin = vvelRelResidualNormFin = pressureRelRisdualNormFin =
      continuityRelResidualNormFin = uvelRelResidualNormInit = vvelRelResidualNormInit =
      pressureRelRisdualNormInit = continuityRelResidualNormInit =0;
  converged = false;

  double fractStep = tStep / m_numFractionalSteps;
  double currentStep = fractStep;

  m_currentCoeff = m_currentCoeff == 0 ? 0 : 0;

  //Friction
  {
#ifdef USE_OPENMP
#pragma omp parallel sections
#endif
    {

#ifdef USE_OPENMP
#pragma omp section
#endif
      {
        calculateFriction(0);
      }

#ifdef USE_OPENMP
#pragma omp section
#endif
      {
        calculateFriction(1);
      }
    }
  }


  while (numIterations < m_itersPerTimeStep)
  {

    ErrorCode uerror = ErrorCode::NoError;
    ErrorCode verror = ErrorCode::NoError;
    QString uErrorMessage;
    QString vErrorMessage;
    uvelRelResidualNormFin = vvelRelResidualNormFin =
        pressureRelRisdualNormFin = continuityRelResidualNormFin = 0.0;

    int numUVelSolvIters = 0;
    int numVVelSolvIters = 0;
    int numPSolvIters  = 0;

    currentStep += fractStep;
    currentStep = min(currentStep, tStep);

    if (m_numWetCells)
    {


      //Solve Momentum
      {
#ifdef USE_OPENMP
#pragma omp parallel sections
#endif
        {

#ifdef USE_OPENMP
#pragma omp section
#endif
          {
            uerror = solveMomentumEquations(currentStep, 0, uvelRelResidualNormFin, numUVelSolvIters, uErrorMessage);
          }

#ifdef USE_OPENMP
#pragma omp section
#endif
          {
            verror = solveMomentumEquations(currentStep, 1, vvelRelResidualNormFin, numVVelSolvIters, vErrorMessage);
          }
        }

        if (uerror == ErrorCode::CriticalFailure)
        {
          error = uerror;
          errorMessage = uErrorMessage;
          break;
        }

        if (verror == ErrorCode::CriticalFailure)
        {
          error = verror;
          errorMessage = vErrorMessage;
          break;
        }
      }
    }

    calculateFaceVelocities();

    error = solvePressureCorrection(currentStep, pressureRelRisdualNormFin, 1.0, numPSolvIters, errorMessage);

    minUVelSolvIters = min(minUVelSolvIters, numUVelSolvIters);
    minVVelSolvIters = min(minVVelSolvIters, numVVelSolvIters);
    minPressSolvIters = min(minPressSolvIters,numPSolvIters);

    maxUVelSolvIters = max(maxUVelSolvIters, numUVelSolvIters);
    maxVVelSolvIters = max(maxVVelSolvIters, numVVelSolvIters);
    maxPressSolvIters = max(maxPressSolvIters, numPSolvIters);

    if (error == ErrorCode::CriticalFailure)
    {
      break;
    }

    //Velocity & Pressure Correction
    {
#ifdef USE_OPENMP
#pragma omp parallel sections
#endif
      {

#ifdef USE_OPENMP
#pragma omp section
#endif
        {
          applyPressureCorrections();
          calculateCellWSEGradients();
        }

#ifdef USE_OPENMP
#pragma omp section
#endif
        {
          applyVelocityCorrections();
          calculateCellVelocityGradients();
          calculateCellEddyViscosities();
        }
      }
    }

    calculateContinuityResiduals(currentStep, continuityRelResidualNormFin);

    if (numIterations == 0)
    {
      uvelRelResidualNormInit = uvelRelResidualNormFin;
      vvelRelResidualNormInit = vvelRelResidualNormFin;
      pressureRelRisdualNormInit = pressureRelRisdualNormFin;
      continuityRelResidualNormInit = continuityRelResidualNormFin;
    }

    //check for convergence
    {
      if (uvelRelResidualNormFin < m_uvelConvergenceTol &&
          vvelRelResidualNormFin < m_vvelConvergenceTol &&
          fabs(pressureRelRisdualNormFin) < m_pressureConvergenceTol &&
          continuityRelResidualNormFin < m_contConvergenceTol &&
          currentStep == tStep)
      {
        numIterations++;
        converged = true;
        break;
      }
    }

    numIterations++;
  }

  if (error == ErrorCode::NoError)
  {
    if (!converged)
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

int FVHMComponent::performPISOTimeStep(double tStep, int &numIterations,
                                       double &uvelRelResidualNormInit, double &uvelRelResidualNormFin, int &minUVelSolvIters, int &maxUVelSolvIters,
                                       double &vvelRelResidualNormInit, double &vvelRelResidualNormFin, int &minVVelSolvIters, int &maxVVelSolvIters,
                                       double &pressureRelRisdualNormInit, double &pressureRelRisdualNormFin, int &minPressSolvIters, int &maxPressSolvIters,
                                       double &continuityRelResidualNormInit, double &continuityRelResidualNormFin,
                                       bool &converged, QString &errorMessage)
{
  ErrorCode error = ErrorCode::NoError;

  numIterations = 0;

  uvelRelResidualNormFin = vvelRelResidualNormFin = pressureRelRisdualNormFin =
      continuityRelResidualNormFin = uvelRelResidualNormInit = vvelRelResidualNormInit =
      pressureRelRisdualNormInit = continuityRelResidualNormInit = 0.0;

  converged = false;

  double fractStep = tStep / m_numFractionalSteps;
  double currentStep = fractStep;

  m_currentCoeff = m_currentCoeff == 0 ? 0 : 0;


  //Friction
  {
#ifdef USE_OPENMP
#pragma omp parallel sections
#endif
    {

#ifdef USE_OPENMP
#pragma omp section
#endif
      {
        calculateFriction(0);
      }

#ifdef USE_OPENMP
#pragma omp section
#endif
      {
        calculateFriction(1);
      }
    }
  }

  while (numIterations < m_itersPerTimeStep)
  {

    ErrorCode uerror = ErrorCode::NoError;
    ErrorCode verror = ErrorCode::NoError;
    QString uErrorMessage;
    QString vErrorMessage;
    uvelRelResidualNormFin = vvelRelResidualNormFin =
        pressureRelRisdualNormFin = continuityRelResidualNormFin = 0.0;

    int numUVelSolvIters = 0;
    int numVVelSolvIters = 0;
    int numPSolvIters  = 0;

    currentStep += fractStep;
    currentStep = min(currentStep, tStep);

    if (m_numWetCells)
    {

      //Solve Momentum
      {
#ifdef USE_OPENMP
#pragma omp parallel sections
#endif
        {

#ifdef USE_OPENMP
#pragma omp section
#endif
          {
            uerror = solveMomentumEquations(currentStep, 0, uvelRelResidualNormFin, numUVelSolvIters, uErrorMessage);
          }

#ifdef USE_OPENMP
#pragma omp section
#endif
          {
            verror = solveMomentumEquations(currentStep, 1, vvelRelResidualNormFin, numVVelSolvIters, vErrorMessage);
          }
        }

        if (uerror == ErrorCode::CriticalFailure)
        {
          error = uerror;
          errorMessage = uErrorMessage;
          break;
        }

        if (verror == ErrorCode::CriticalFailure)
        {
          error = verror;
          errorMessage = vErrorMessage;
          break;
        }
      }

    }

    for(int m = 0; m < 2; m++)
    {
      calculateFaceVelocities();

      error = solvePressureCorrection(currentStep, pressureRelRisdualNormFin, 0.0, numPSolvIters, errorMessage);

      minUVelSolvIters = min(minUVelSolvIters, numUVelSolvIters);
      minVVelSolvIters = min(minVVelSolvIters, numVVelSolvIters);
      minPressSolvIters = min(minPressSolvIters,numPSolvIters);

      maxUVelSolvIters = max(maxUVelSolvIters, numUVelSolvIters);
      maxVVelSolvIters = max(maxVVelSolvIters, numVVelSolvIters);
      maxPressSolvIters = max(maxPressSolvIters, numPSolvIters);

      if (error == ErrorCode::CriticalFailure)
      {
        break;
      }

      //Velocity & Pressure Correction
      {
#ifdef USE_OPENMP
#pragma omp parallel sections
#endif
        {

#ifdef USE_OPENMP
#pragma omp section
#endif
          {
            applyPressureCorrections();
            calculateCellWSEGradients();
          }

#ifdef USE_OPENMP
#pragma omp section
#endif
          {
            applyVelocityCorrections();
            calculateCellVelocityGradients();
            calculateCellEddyViscosities();
          }
        }
      }
    }


    calculateContinuityResiduals(currentStep, continuityRelResidualNormFin);

    if (numIterations == 0)
    {
      uvelRelResidualNormInit = uvelRelResidualNormFin;
      vvelRelResidualNormInit = vvelRelResidualNormFin;
      pressureRelRisdualNormInit = pressureRelRisdualNormFin;
      continuityRelResidualNormInit = continuityRelResidualNormFin;
    }

    //check for convergence
    {
      if (uvelRelResidualNormFin < m_uvelConvergenceTol &&
          vvelRelResidualNormFin < m_vvelConvergenceTol &&
          fabs(pressureRelRisdualNormFin) < m_pressureConvergenceTol &&
          continuityRelResidualNormFin < m_contConvergenceTol &&
          currentStep == tStep)
      {
        numIterations++;
        converged = true;
        break;
      }
    }

    numIterations++;
  }

  if (error == ErrorCode::NoError)
  {
    if (!converged)
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

FVHMComponent::ErrorCode FVHMComponent::solveMomentumEquations(double tStep, int uv, double &uvelRelResidualNorm, int &numSolvIters, QString &errorMessage)
{
  //momentum
  double *u_b = new double[m_numWetCells];
  SparseMatrix coeffMatrix(0,m_numWetCells-1,m_numWetCells, 4);
  double *x = new double[m_numWetCells];

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int index = 0 ; index < m_numWetCells ; index++)
  {
    int cvIndex = m_wetCells[index];
    TriCV *cv = m_controlVolumes[cvIndex];

    for(int u = 0; u < cv->numEdges + 1 ; u++)
    {
      cv->velCoeffs[uv][u] = 0.0;
    }

    double a_p = cv->h->value * cv->area / tStep;

    double b = 0.0;

    for(int i = 0; i < cv->numEdges; i++)
    {
      FaceNormVelBC &faceVel = cv->faceNormalVels[i];
      double evel = faceVel.value;

      if(evel)
      {

        int cvIndex = cv->nTris[i];
        double edgeDepth = cv->faceDepths[i].value;
        double tempCoeff = m_useAdvection * evel * cv->r_eta[i] * edgeDepth;

        if(cvIndex > -1)
        {
          TriCV *cvn = cvn = m_controlVolumes[cvIndex];

          if(cvn->wetIndex)
          {
            double a_n = 0.0;
            double wl_n = cv->w_l[i];
            double wl_p = 1.0 - wl_n;

            Vect &r_xi = cv->r_xi[i];

            //TVD
            if(evel > 0)
            {
              double ru = calculateFluxLimiter(cv, i, wl_n, cv->grad_vel[uv], cvn->grad_vel[uv], cv->vel[uv].value, cvn->vel[uv].value, r_xi.v[0] , r_xi.v[1]);
              a_p += tempCoeff * (1.0 - ru * wl_n);
              a_n += tempCoeff * ru * wl_n;
            }
            else
            {
              double ru = calculateFluxLimiter(cv, i, wl_p, cvn->grad_vel[uv], cv->grad_vel[uv], cvn->vel[uv].value, cv->vel[uv].value, -r_xi.v[0] , -r_xi.v[1]);
              a_p += tempCoeff *  ru * wl_p;
              a_n += tempCoeff *  (1.0 - ru * wl_p);
            }

            double faceTurbulentVisco = wl_p * cv->eddyViscosity + wl_n * cvn->eddyViscosity;

            //direct diffusion
            a_n -= m_useEddyViscosity * edgeDepth * (faceTurbulentVisco + m_viscosity) * cv->dir_diff_term[i];
            a_p += m_useEddyViscosity * edgeDepth * (faceTurbulentVisco + m_viscosity) * cv->dir_diff_term[i];

            //cross diffusion
            b += m_useEddyViscosity * edgeDepth * (faceTurbulentVisco + m_viscosity) *
                 cv->cross_diff_term[i] * (cv->nodeVels[uv][i + 1] - cv->nodeVels[uv][i]);

            cv->velCoeffs[uv][i + 1] = a_n;
          }
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

    //del_v/del_t
    double dv_dt = cv->prevH->value * cv->prevVel[uv].value * cv->area / tStep;

    //pressure gradient
    double gradp = -g * cv->grad_z->v[uv] * cv->h->value * cv->area;

    //external force
    double extforce = cv->externalForce->v[uv] * cv->h->value * cv->area;

    double friction = 0.0;

    //friction

    //if(!isInfOrNan(cv->friction[uv]))
    {
      friction = -cv->friction[uv];

      if(m_useWall)
      {
        for(int f = 0 ; f < cv->numEdges; f++)
        {
          friction -= cv->wallShearFriction[uv][f];
        }
      }
    }

    b += dv_dt + gradp + friction + extforce;

    //const coefficients
    u_b[cv->wetCellIndex] = b;

    //initial guess
    x[cv->wetCellIndex] = cv->vel[uv].value;

  }

  //set up momentum equations calculate u and solve
  double *residuals = new double[m_numWetCells]();

  ErrorCode error = ErrorCode::NoError;

#ifdef USE_OPENMP
#pragma omp critical
#endif
  {
    error = mpiSolve(coeffMatrix, u_b, x, residuals, uvelRelResidualNorm, numSolvIters, errorMessage);
  }

  if(error)
  {
    for(int index = 0; index < m_numWetCells ; index++)
    {
      int cvIndex = m_wetCells[index];
      TriCV *cv = m_controlVolumes[cvIndex];
      cv->printDetails();
    }
  }

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int index = 0; index < m_numWetCells ; index++)
  {
    int cvIndex = m_wetCells[index];
    TriCV *cv = m_controlVolumes[cvIndex];
    cv->prevIterVel[uv] = cv->vel[uv].value;
    cv->vel[uv].value = x[cv->wetCellIndex];
    cv->velResidualIter[uv] = residuals[cv->wetCellIndex];
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
  for(int c = 0; c < m_numContCells; c++)
  {
    int cvIndex = m_contCells[c];
    TriCV *cv =  m_controlVolumes[cvIndex];
    cv->calculateFaceVelocities();
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
          //          {
          //            rf = Vect::dotProduct(2.0 * grad_c.v[0], 2.0 * grad_c.v[1], 0.0, r_cdX, r_cdY, 0.0) / du;
          //            rf = rf - 1.0;
          //          }

          //Ubbink and Issa (1999)
          {
            double gradFaceX = (1.0 - idwfactor)* grad_c.v[0] + idwfactor * grad_d.v[0];
            double gradFaceY = (1.0 - idwfactor)* grad_c.v[1] + idwfactor * grad_d.v[1];
            double phi_u = phi_d - 2.0 * Vect::dotProduct(gradFaceX, gradFaceY, 0.0, r_cdX, r_cdY, 0.0);
            rf = (phi_c - phi_u)/(phi_d - phi_c);
          }

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
          //          {
          //            rf = Vect::dotProduct(2.0 * grad_c.v[0], 2.0 * grad_c.v[1], 0.0, r_cdX, r_cdY, 0.0) / du;
          //            rf = rf - 1.0;
          //          }

          //Ubbink and Issa (1999)
          {
            double gradFaceX = (1.0 - idwfactor)* grad_c.v[0] + idwfactor * grad_d.v[0];
            double gradFaceY = (1.0 - idwfactor)* grad_c.v[1] + idwfactor * grad_d.v[1];
            double phi_u = phi_d - 2.0 * Vect::dotProduct(gradFaceX, gradFaceY, 0.0, r_cdX, r_cdY, 0.0);
            rf = (phi_c - phi_u)/(phi_d - phi_c);
          }

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
          //          {
          //            rf = Vect::dotProduct(2.0 * grad_c.v[0], 2.0 * grad_c.v[1], 0.0, r_cdX, r_cdY, 0.0) / du;
          //            rf = rf - 1.0;
          //          }

          //Ubbink and Issa (1999)
          {
            double gradFaceX = (1.0 - idwfactor)* grad_c.v[0] + idwfactor * grad_d.v[0];
            double gradFaceY = (1.0 - idwfactor)* grad_c.v[1] + idwfactor * grad_d.v[1];
            double phi_u = phi_d - 2.0 * Vect::dotProduct(gradFaceX, gradFaceY, 0.0, r_cdX, r_cdY, 0.0);
            rf = (phi_c - phi_u)/(phi_d - phi_c);
          }

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
          //          double phi_u = phi_d - 2.0 * Vect::dotProduct(gradFaceX, gradFaceY, 0.0, r_cdX, r_cdY, 0.0);
          //          rf = (phi_c - phi_u)/(phi_d - phi_c);
          //          }

          ru = (rf + rf * rf)/(1 + rf * rf);

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
      //Superbee
    case 6:
      {
        double du = phi_d - phi_c;

        if(fabs(du - 0.0) >= m_epsilon)
        {
          double rf = 0.0;

          //Darwish and Moukalled (2002)
          //          {
          //            rf = Vect::dotProduct(2.0 * grad_c.v[0], 2.0 * grad_c.v[1], 0.0, r_cdX, r_cdY, 0.0) / du;
          //            rf = rf - 1.0;
          //          }


          //Ubbink and Issa (1999)
          {
            double gradFaceX = (1.0 - idwfactor)* grad_c.v[0] + idwfactor * grad_d.v[0];
            double gradFaceY = (1.0 - idwfactor)* grad_c.v[1] + idwfactor * grad_d.v[1];
            double phi_u = phi_d - 2.0 * Vect::dotProduct(gradFaceX, gradFaceY, 0.0, r_cdX, r_cdY, 0.0);
            rf = (phi_c - phi_u)/(phi_d - phi_c);
          }

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
          //          {
          //            rf = Vect::dotProduct(2.0 * grad_c.v[0], 2.0 * grad_c.v[1], 0.0, r_cdX, r_cdY, 0.0) / du;
          //            rf = rf - 1.0;
          //          }


          //Ubbink and Issa (1999)
          {
            double gradFaceX = (1.0 - idwfactor)* grad_c.v[0] + idwfactor * grad_d.v[0];
            double gradFaceY = (1.0 - idwfactor)* grad_c.v[1] + idwfactor * grad_d.v[1];
            double phi_u = phi_d - 2.0 * Vect::dotProduct(gradFaceX, gradFaceY, 0.0, r_cdX, r_cdY, 0.0);
            rf = (phi_c - phi_u)/(phi_d - phi_c);
          }

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
          //          {
          //            rf = Vect::dotProduct(2.0 * grad_c.v[0], 2.0 * grad_c.v[1], 0.0, r_cdX, r_cdY, 0.0) / du;
          //            rf = rf - 1.0;
          //          }

          //Ubbink and Issa (1999)
          {
            double gradFaceX = (1.0 - idwfactor)* grad_c.v[0] + idwfactor * grad_d.v[0];
            double gradFaceY = (1.0 - idwfactor)* grad_c.v[1] + idwfactor * grad_d.v[1];
            double phi_u = phi_d - 2.0 * Vect::dotProduct(gradFaceX, gradFaceY, 0.0, r_cdX, r_cdY, 0.0);
            rf = (phi_c - phi_u)/(phi_d - phi_c);
          }

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
          //          {
          //            rf = Vect::dotProduct(2.0 * grad_c.v[0], 2.0 * grad_c.v[1], 0.0, r_cdX, r_cdY, 0.0) / du;
          //            rf = rf - 1.0;
          //          }

          //Ubbink and Issa (1999)
          {
            double gradFaceX = (1.0 - idwfactor)* grad_c.v[0] + idwfactor * grad_d.v[0];
            double gradFaceY = (1.0 - idwfactor)* grad_c.v[1] + idwfactor * grad_d.v[1];
            double phi_u = phi_d - 2.0 * Vect::dotProduct(gradFaceX, gradFaceY, 0.0, r_cdX, r_cdY, 0.0);
            rf = (phi_c - phi_u)/(phi_d - phi_c);
          }

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
          //          {
          //            rf = Vect::dotProduct(2.0 * grad_c.v[0], 2.0 * grad_c.v[1], 0.0, r_cdX, r_cdY, 0.0) / du;
          //            rf = rf - 1.0;
          //          }


          //Ubbink and Issa (1999)
          {
            double gradFaceX = (1.0 - idwfactor)* grad_c.v[0] + idwfactor * grad_d.v[0];
            double gradFaceY = (1.0 - idwfactor)* grad_c.v[1] + idwfactor * grad_d.v[1];
            double phi_u = phi_d - 2.0 * Vect::dotProduct(gradFaceX, gradFaceY, 0.0, r_cdX, r_cdY, 0.0);
            rf = (phi_c - phi_u)/(phi_d - phi_c);
          }

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
          //          {
          //            rf = Vect::dotProduct(2.0 * grad_c.v[0], 2.0 * grad_c.v[1], 0.0, r_cdX, r_cdY, 0.0) / du;
          //            rf = rf - 1.0;
          //          }

          //Ubbink and Issa (1999)
          {
            double gradFaceX = (1.0 - idwfactor)* grad_c.v[0] + idwfactor * grad_d.v[0];
            double gradFaceY = (1.0 - idwfactor)* grad_c.v[1] + idwfactor * grad_d.v[1];
            double phi_u = phi_d - 2.0 * Vect::dotProduct(gradFaceX, gradFaceY, 0.0, r_cdX, r_cdY, 0.0);
            rf = (phi_c - phi_u)/(phi_d - phi_c);
          }

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
          //          {
          //            rf = Vect::dotProduct(2.0 * grad_c.v[0], 2.0 * grad_c.v[1], 0.0, r_cdX, r_cdY, 0.0) / du;
          //            rf = rf - 1.0;
          //          }

          //Ubbink and Issa (1999)
          {
            double gradFaceX = (1.0 - idwfactor)* grad_c.v[0] + idwfactor * grad_d.v[0];
            double gradFaceY = (1.0 - idwfactor)* grad_c.v[1] + idwfactor * grad_d.v[1];
            double phi_u = phi_d - 2.0 * Vect::dotProduct(gradFaceX, gradFaceY, 0.0, r_cdX, r_cdY, 0.0);
            rf = (phi_c - phi_u)/(phi_d - phi_c);
          }

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

FVHMComponent::ErrorCode FVHMComponent::solvePressureCorrection(double tStep, double &pressureRelResidualNorm, double minorTermsCoeff, int &numSolvIters, QString &errorMessage)
{
  numSolvIters = 0;

  ErrorCode error = ErrorCode::NoError;
  double *x = new double[m_numContCells];
  double *h_b = new double[m_numContCells];
  double *residuals = new double[m_numContCells];

  SparseMatrix coeffMatrix(0,m_numContCells - 1,m_numContCells, 4);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int index = 0 ; index < m_numContCells ; index++)
  {
    int cvIndex = m_contCells[index];
    TriCV *cv = m_controlVolumes[cvIndex];

    x[cv->contCellIndex] = 0.0;

    double zCorrEst = 0;

    double constTerm  = 0;

    double h_p = cv->area / tStep;

    double sumEdgeCoeffs = 0;

    for(int e = 0; e < cv->numEdges; e++)
    {
      sumEdgeCoeffs += cv->velCoeffs[0][e+1];
    }

    double coeffucv = cv->velCoeffs[0][0] + minorTermsCoeff * sumEdgeCoeffs;

    double *coefficients = new double[cv->numEdges + 1]();

    for(int i = 0; i < cv->numEdges; i++)
    {
      coefficients[i+1] = 0.0;
      double edgeDepth = cv->faceDepths[i].value;
      double eta = cv->r_eta[i];
      double evel = cv->faceNormalVels[i].value;

      int cvnIndex = cv->nTris[i];

      if(cvnIndex > -1)
      {
        TriCV *cvn = m_controlVolumes[cvnIndex];
        double wl_n = cv->w_l[i];
        double wl_p = 1.0 - wl_n;
        double h_n = 0.0;

        if(cv->wetIndex && cvn->wetIndex)
        {
          sumEdgeCoeffs = 0;

          for(int e = 0; e < cvn->numEdges; e++)
          {
            sumEdgeCoeffs += cvn->velCoeffs[0][e+1];
          }

          double coeffucvn = cvn->velCoeffs[0][0] + minorTermsCoeff * sumEdgeCoeffs;

          double A_e = (wl_p * cv->area / coeffucv + wl_n * cvn->area / coeffucvn);

          double velCorrection = -g * A_e *  cv->h->value / cv->r_xi_l[i];

          h_p -= edgeDepth * eta * velCorrection * m_velocityRelaxFactor;
          h_n += edgeDepth * eta * velCorrection * m_velocityRelaxFactor;

        }

        h_p += wl_p * evel * eta;
        h_n += wl_n * evel * eta;

        double zcorrX = wl_p * cv->grad_zcorr->v[0] + wl_n * cvn->grad_zcorr->v[0];
        double zcorrY = wl_p * cv->grad_zcorr->v[1] + wl_n * cvn->grad_zcorr->v[1];
        constTerm -= evel * eta * Vect::dotProduct(zcorrX, zcorrY, 0.0, cv->df[i]);

        coefficients[i+1] = h_n;
      }

      double edgeFlow = eta * edgeDepth * evel;
      constTerm -= edgeFlow;
      zCorrEst += edgeFlow;
    }

    coefficients[0] = h_p;

    for(int i = 0; i < cv->numEdges + 1; i++)
    {
      int cvnIndex = cv->orderedNTris[i];
      double coeff = 0.0;

      if(cvnIndex > -1 && (coeff = coefficients[cv->orderedNTrisIndexes[i]]))
      {
        coeffMatrix.appendValue(cv->contCellIndex, m_controlVolumes[cvnIndex]->contCellIndex, coeff);
      }
    }

    delete[] coefficients;

    //constTerm
    constTerm += cv->inflowOutflow;

    //tstep terms
    constTerm += (cv->prevH->value - cv->h->value) * cv->area  / tStep;

    //set const value
    h_b[cv->contCellIndex] = constTerm;

    zCorrEst += constTerm;
    zCorrEst *= tStep / cv->area;
    x[cv->contCellIndex] = zCorrEst;
  }

  int numIters = 0;
  error = mpiSolve(coeffMatrix, h_b, x, residuals, pressureRelResidualNorm, numIters, errorMessage);
  numSolvIters = max(numIters, numSolvIters);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int rc = 0 ; rc < m_numContCells ; rc++)
  {
    int cvIndex = m_contCells[rc];
    TriCV * cv = m_controlVolumes[cvIndex];
    double zCorr = x[cv->contCellIndex];
    cv->zCorrection = zCorr;
  }

  if(error)
  {
    for(int rc = 0 ; rc < m_numContCells ; rc++)
    {
      int cvIndex = m_contCells[rc];
      TriCV * cv = m_controlVolumes[cvIndex];
      cv->printDetails();
    }
  }

  //Cheap
  //  for(int it = 0; it < m_numPressureCorrectionIterations; it++)
  //  {
  //#ifdef USE_OPENMP
  //#pragma omp parallel for
  //#endif
  //    for(int index = 0 ; index < m_numContCells ; index++)
  //    {
  //      int cvnIndex  = m_contCells[index];
  //      TriCV *cv = m_controlVolumes[cvnIndex];

  //      double zCorrection = 0.0;

  //      zCorrection = cv->inflowOutflow + (cv->prevH->value - cv->h->value) * cv->area / tStep;

  //      for(int i = 0; i < cv->numEdges; i++)
  //      {
  //        zCorrection -= cv->faceNormalVels[i].value * cv->faceDepths[i].value * cv->r_eta[i];
  //      }

  //      zCorrection =  zCorrection * tStep / cv->area;
  //      x[cv->index] = zCorrection;
  //      cv->zCorrection = zCorrection;
  //    }
  //  }

  pressureRelResidualNorm = l2Norm(x,m_numContCells);

  delete[] x;
  delete[] residuals;
  delete[] h_b;

  return error;
}

void FVHMComponent::applyPressureCorrections()
{

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int c = 0 ; c < m_numContCells ; c++)
  {
    int cvIndex = m_contCells[c];
    TriCV *cv = m_controlVolumes[cvIndex];

    if(cv->z->isBC == false)
    {
      double hNew = cv->h->value + m_pressureRelaxFactor * cv->zCorrection;
      cv->setVFRDepth(hNew);
    }
  }
}

void FVHMComponent::applyVelocityCorrections()
{

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int c = 0 ; c < m_numWetCells ; c++)
  {
    int cvIndex = m_wetCells[c];
    TriCV *cv = m_controlVolumes[cvIndex];
    cv->vel[0].value = m_velocityRelaxFactor * cv->vel[0].value  + (1.0 - m_velocityRelaxFactor) * cv->prevIterVel[0];
    cv->vel[1].value = m_velocityRelaxFactor * cv->vel[1].value  + (1.0 - m_velocityRelaxFactor) * cv->prevIterVel[1];
  }
}

void FVHMComponent::calculateContinuityResiduals(double tStep, double &contRelResidualNorm)
{
  double *residuals = new double[m_numContCells];

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int rc = 0 ; rc < m_numContCells ; rc++)
  {
    int cvIndex = m_contCells[rc];
    TriCV *cv = m_controlVolumes[cvIndex];

    double sumflow = -cv->inflowOutflow * tStep;

    for(int i = 0; i < cv->numEdges ; i++)
    {
      sumflow += cv->faceDepths[i].value * cv->r_eta[i] * cv->faceNormalVels[i].value * tStep;
    }

    sumflow += (cv->h->value - cv->prevH->value) * cv->area;
    residuals[cv->contCellIndex] = sumflow;
    cv->contResidualIter = sumflow;
  }

  contRelResidualNorm = l2Norm(residuals,m_numContCells);

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

    if(cv->inflowOutflow > 0)
    {
      cv->totalExternalInflow += cv->inflowOutflow * tStep;
    }
    else if(cv->inflowOutflow)
    {
      cv->totalExternalOutflow += -cv->inflowOutflow * tStep;
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
  for(int i = 0; i < m_numWetCells; i++)
  {
    int cvIndex = m_wetCells[i];
    TriCV *cv = m_controlVolumes[cvIndex];
    calculateFriction(cv, uv);
  }
}

void FVHMComponent::calculateFriction(TriCV *cv, int uv)
{
  double cf = g * cv->mannings * cv->mannings / pow(cv->h->value, 1.0/3.0);
  double umag = hypot(cv->vel[0].value, cv->vel[1].value);
  double frictionVel = sqrt(cf * umag * umag);

  double friction =  cf * umag * cv->vel[uv].value * cv->area;
  cv->friction[uv] = friction;

  if(m_useWall)
  {
    for(int i = 0 ; i < cv->numEdges ; i++)
    {
      //int xy = (uv + 1) % 2;
      double relp = cv->r_e_l_p[i];//.v[xy];

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

  double rux = 0.0;
  double ruy = 0.0;

  double rvx = 0.0;
  double rvy = 0.0;

  int iters = 0;

  do
  {

    rux = 0.0;
    ruy = 0.0;

    rvx = 0.0;
    rvy = 0.0;

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < m_numContCells; i++)
    {
      int cvIndex = m_contCells[i];
      TriCV *cv = m_controlVolumes[cvIndex];

      double drux = cv->grad_vel[0].v[0];
      double druy = cv->grad_vel[0].v[1];

      double drvx = cv->grad_vel[1].v[0];
      double drvy = cv->grad_vel[1].v[1];

      cv->calculateVelocityGradient();

      drux = fabs(drux - cv->grad_vel[0].v[0]);
      druy = fabs(druy - cv->grad_vel[0].v[1]);

      drvx = fabs(drvx - cv->grad_vel[1].v[0]);
      drvy = fabs(drvy - cv->grad_vel[1].v[1]);

      if (drux > rux)
      {
#ifdef USE_OPENMP
#pragma omp atomic read
#endif
        rux = drux;
      }

      if (druy > ruy)
      {
#ifdef USE_OPENMP
#pragma omp atomic read
#endif
        ruy = druy;
      }

      if (drvx > rvx)
      {
#ifdef USE_OPENMP
#pragma omp atomic read
#endif
        rvx = drvx;
      }

      if (drvy > rvy)
      {
#ifdef USE_OPENMP
#pragma omp atomic read
#endif
        rvy = drvy;
      }
    }

    iters++;

  } while (rux > 1e-6 && ruy > 1e-6 && rvx > 1e-6 && rvy > 1e-6 && iters < 1000);


#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < m_numContCells; i++)
  {
    int cvIndex = m_contCells[i];
    TriCV *cv = m_controlVolumes[cvIndex];
    cv->calculateNodeVelocities();
  }

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
          for(int i = 0 ; i < m_numWetCells ; i++)
          {
            int cvIndex = m_wetCells[i];
            //Smargorinsky. Need to check to make sure correct.
            TriCV *cv = m_controlVolumes[cvIndex];
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
          for(int i = 0 ; i < m_numWetCells ; i++)
          {
            int cvIndex = m_wetCells[i];
            TriCV *cv = m_controlVolumes[cvIndex];
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
    for(int i = 0 ; i < m_numWetCells ; i++)
    {
      int cvIndex = m_wetCells[i];
      TriCV *cv = m_controlVolumes[cvIndex];
      cv->interpolateNodeEddyViscosities();
    }
  }
}

void FVHMComponent::calculateCellWSEGradients()
{

#ifdef USE_OPENMP
#pragma omp parallel sections
#endif
  {

#ifdef USE_OPENMP
#pragma omp section
#endif
    {
      double zx = 0.0;
      double zy = 0.0;

      int iters = 0;

      do
      {

        zx = 0.0;
        zy = 0.0;


#ifdef USE_OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < m_numContCells; i++)
        {
          int cvIndex = m_contCells[i];
          TriCV *cv = m_controlVolumes[cvIndex];

          if (!cv->grad_z->isBC)
          {
            double gradzx = cv->grad_z->v[0];
            double gradzy = cv->grad_z->v[1];

            cv->calculateWSEGradient();

            double tmpzx = fabs(cv->grad_z->v[0] - gradzx);
            double tmpzy = fabs(cv->grad_z->v[1] - gradzy);

            if (tmpzx > zx)
            {
#ifdef USE_OPENMP
#pragma omp atomic read
#endif
              zx = tmpzx;
            }

            if (tmpzy > zy)
            {
#ifdef USE_OPENMP
#pragma omp atomic read
#endif
              zy = tmpzy;
            }

          }
        }

        iters++;

      } while (zx > 1e-6 && zy > 1e-6 && iters < 1000);

      if (iters >= 1000)
      {
        printf("dh_x:%f\tdh_y:%f\tNum iters: %i\n", zx, zy, iters);
      }
    }

#ifdef USE_OPENMP
#pragma omp section
#endif
    {
      double hx = 0.0;
      double hy = 0.0;

      int iters = 0;

      do
      {

        hx = 0.0;
        hy = 0.0;


#ifdef USE_OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < m_numContCells; i++)
        {
          int cvIndex = m_contCells[i];
          TriCV *cv = m_controlVolumes[cvIndex];

          double gradhx = cv->grad_h->v[0];
          double gradhy = cv->grad_h->v[1];

          cv->calculateDepthGradient();

          double tmphx = fabs(cv->grad_h->v[0] - gradhx);
          double tmphy = fabs(cv->grad_h->v[1] - gradhy);

          if (tmphx > hx)
          {
#ifdef USE_OPENMP
#pragma omp atomic read
#endif
            hx = tmphx;
          }

          if (tmphy > hy)
          {
#ifdef USE_OPENMP
#pragma omp atomic read
#endif
            hy = tmphy;
          }
        }

        iters++;

      } while (hx > 1e-6 && hy > 1e-6 && iters < 1000);

      if (iters >= 1000)
      {
        printf("dz_x: %f\tdz_y: %f\tNum iters: %i\n", hx, hy, iters);
      }
    }

#ifdef USE_OPENMP
#pragma omp section
#endif
    {
      calculateCellZCorrGradients();
    }
  }


#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < m_numCells; i++)
  {
    TriCV *cv = m_controlVolumes[i];
    cv->calculateEdgeDepths();
    cv->calculateNodeElevations();
  }

}

void FVHMComponent::calculateCellZCorrGradients()
{
  double zx = 0.0;
  double zy = 0.0;

  int iters = 0;

  do
  {

    zx = 0.0;
    zy = 0.0;


#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < m_numContCells; i++)
    {
      int cvIndex = m_contCells[i];
      TriCV *cv = m_controlVolumes[cvIndex];

      double gradzx = cv->grad_zcorr->v[0];
      double gradzy = cv->grad_zcorr->v[1];

      cv->calculateZCorrGradient();

      double tmpzx = fabs(cv->grad_zcorr->v[0] - gradzx);
      double tmpzy = fabs(cv->grad_zcorr->v[1] - gradzy);


      {
        if (tmpzx > zx)
        {
#ifdef USE_OPENMP
#pragma omp atomic read
#endif
          zx = tmpzx;
        }

        if (tmpzy > zy)
        {
#ifdef USE_OPENMP
#pragma omp atomic read
#endif
          zy = tmpzy;
        }

      }
    }

    iters++;

  } while (zx > 1e-6 && zy > 1e-6 && iters < 1000);

  if (iters >= 1000)
  {
    printf("dzcorr_x: %f\tdzcorr_y: %f\tNum iters: %i\n", zx, zy, iters);
  }
}

FVHMComponent::ErrorCode FVHMComponent::mpiSolve(const SparseMatrix &A, const double b[], double x[],
                                                 double residuals[], double &relativeResidualNorm, int &numIterations, QString &errorMessage)
{
  ErrorCode solved = ErrorCode::NoError;
  int result = 0;

  if(mpiProcessRank() == 0)
  {
    int numMPIProcs = m_mpiAllocatedProcessesArray.size();
    int rowsPerProc = A.rowCount() / numMPIProcs;
    int remRows = A.rowCount() - rowsPerProc * numMPIProcs;

    if(numMPIProcs == 1 || rowsPerProc == 0 || A.rowCount() < m_mpiSolverSplitThreshold)
    {
      result = solve(MPI_COMM_SELF, A,A.ilower(), A.iupper(),b,x,residuals,relativeResidualNorm, numIterations);
    }
    //divide and conquer
    else
    {

      int start = rowsPerProc + remRows;

      for(int i = 1; i < m_mpiAllocatedProcessesArray.size() ; i++)
      {
        int procRank =  m_mpiAllocatedProcessesArray[i];

        int ilower = start + (i-1) * rowsPerProc;
        int iupper = ilower + rowsPerProc - 1;

        int size = 4 + (A.getDataSize(ilower,iupper) * 2) + ((iupper - ilower + 1) * 3);
        double *serializedData = new double[size];

        int counter = 0;
        A.serializeRows(ilower,iupper, serializedData, counter);

        for(int r = ilower; r <= iupper; r++)
        {
          serializedData[counter] = b[r]; counter++;
        }

        for(int r = ilower; r <= iupper; r++)
        {
          serializedData[counter] = x[r]; counter++;
        }

        MPI_Send(serializedData, size, MPI_DOUBLE, procRank, 1001, MPI_COMM_WORLD);

        delete[] serializedData;
      }

      result = solve(m_ComponentMPIComm, A, A.ilower(), start - 1, b, x, residuals, relativeResidualNorm, numIterations);

      //recieve solutions
      for(int i = 1; i < m_mpiAllocatedProcessesArray.size() ; i++)
      {
        int procRank =  m_mpiAllocatedProcessesArray[i];

        int ilower = start + (i-1) * rowsPerProc;
        int iupper = ilower + rowsPerProc - 1;

        int bufferSize = rowsPerProc * 2 + 2;
        double *values = new double[bufferSize];

        MPI_Status status;
        MPI_Recv(values, bufferSize, MPI_DOUBLE, procRank, 1001, MPI_COMM_WORLD, &status);

        result = max(result, (int)values[0]);
        relativeResidualNorm = max(relativeResidualNorm, values[1]);

        for(int j = ilower; j <= iupper; j++)
        {
          x[j] = values[j - ilower + 2];
          residuals[j] = values[rowsPerProc  + j - ilower + 2];
        }

        delete[] values;
      }
    }
  }
#ifdef USE_MPI
  else
  {
    result = solve(m_ComponentMPIComm, A, A.ilower(), A.iupper(), b, x, residuals, relativeResidualNorm, numIterations);
  }
#endif

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

        char message[100];
        HYPRE_DescribeError(result,message);
        errorMessage = QString(message);
        solved = ErrorCode::CriticalFailure;
      }
      break;
    case 256:
      {
        char message[100];
        HYPRE_DescribeError(result,message);
        errorMessage = QString(message);
        solved = ErrorCode::SolverFailedToConverge;
      }
      break;
  }


  return solved;
}

int FVHMComponent::solve(MPI_Comm communicator, const SparseMatrix &A, int ilower, int iupper,
                         const double b[], double x[], double residuals[], double &relativeResidualNorm, int &numIterations)
{
  HYPRE_Int result = 0;

  HYPRE_IJMatrix mA;
  HYPRE_IJVector mb;
  HYPRE_IJVector mx;

  int numRows = iupper - ilower + 1;
  int dataSize = A.getDataSize(ilower, iupper);

  int *colsPerRow = new int[numRows];
  int *rowIndexes = new int[numRows];
  int *colIndexes = new int[dataSize];
  double *values = new double[dataSize];

  A.getColsPerRow(colsPerRow, ilower, iupper);
  A.getRowIndexes(rowIndexes, ilower, iupper);
  A.getColumnIndexes(colIndexes, ilower, iupper);
  A.getValuesByRows(values, ilower, iupper);

  //A matrix
  HYPRE_IJMatrixCreate(communicator, ilower, iupper, ilower, iupper, &mA);

#ifdef USE_OPENMP
  HYPRE_IJMatrixSetOMPFlag(mA, 1);
#endif

  HYPRE_IJMatrixSetObjectType(mA, HYPRE_PARCSR);
  HYPRE_IJMatrixInitialize(mA);

  HYPRE_IJMatrixSetValues(mA, numRows, colsPerRow, rowIndexes, colIndexes, values);
  HYPRE_IJMatrixAssemble(mA);

  //B vector
  HYPRE_IJVectorCreate(communicator, ilower, iupper, &mb);
  HYPRE_IJVectorSetObjectType(mb, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(mb);

  HYPRE_IJVectorSetValues(mb, numRows, rowIndexes, b);
  HYPRE_IJVectorAssemble(mb);

  //X VECTOR
  HYPRE_IJVectorCreate(communicator, ilower, iupper, &mx);
  HYPRE_IJVectorSetObjectType(mx, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(mx);

  HYPRE_IJVectorSetValues(mx, numRows, rowIndexes, x);
  HYPRE_IJVectorAssemble(mx);

  switch (m_solverType)
  {
    //PCG with AMG
    case 0:
      {
        HYPRE_ParCSRMatrix parcsr_A;
        HYPRE_IJMatrixGetObject(mA, (void **)&parcsr_A);

        HYPRE_ParVector par_b;
        HYPRE_IJVectorGetObject(mb, (void **)&par_b);

        HYPRE_ParVector par_x;
        HYPRE_IJVectorGetObject(mx, (void **)&par_x);

        HYPRE_Solver solver;
        HYPRE_Solver precond;

        int num_iterations;

        HYPRE_ParVector par_residualVector;

        /* Create solver */
        HYPRE_ParCSRPCGCreate(communicator, &solver);

        /* Set some parameters (See Reference Manual for more parameters) */
        HYPRE_PCGSetMaxIter(solver, m_solverMaximumNumberOfIterations); /* max iterations */
        HYPRE_PCGSetTol(solver, m_solverConvergenceTol);                /* conv. tolerance */
        HYPRE_PCGSetTwoNorm(solver, 1);                                 /* use the two norm as the stopping criteria */
        HYPRE_PCGSetPrintLevel(solver, 0);                              /* print solve info */
        HYPRE_PCGSetLogging(solver, 1);                                 /* needed to get run info later */

        /* Now set up the AMG preconditioner and specify any parameters */
        HYPRE_BoomerAMGCreate(&precond);
        HYPRE_BoomerAMGSetPrintLevel(precond, 0); /* print amg solution info */
        HYPRE_BoomerAMGSetCoarsenType(precond, 8);

#ifdef USE_OPENMP
        HYPRE_BoomerAMGSetRelaxType(precond, 0); /* Sym G.S./Jacobi hybrid */
#else
        HYPRE_BoomerAMGSetRelaxType(precond, 6);
#endif

        HYPRE_BoomerAMGSetNumSweeps(precond, 1);
        HYPRE_BoomerAMGSetMaxIter(precond, m_solverAMGPreconditionerNumberOfIterations); /* do only one iteration! */
        HYPRE_BoomerAMGSetTol(precond, 0.0);                                             /* conv. tolerance zero */

        /* Set the FlexGMRES preconditioner */
        HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSolve,
                            (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSetup, precond);

        /* Now setup and solve! */
        HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_b, par_x);

        result = HYPRE_ParCSRPCGSolve(solver, parcsr_A, par_b, par_x);

        /*Get number of iterations used*/
        HYPRE_PCGGetNumIterations(solver, &numIterations);

        /*Now copy solution! */
        HYPRE_IJVectorGetValues(mx, numRows, rowIndexes, x);

        /*Copy residuals! */
        HYPRE_FlexGMRESGetResidual(solver, (void **)&par_residualVector);

        double *residualValues = hypre_VectorData(hypre_ParVectorLocalVector(par_residualVector));

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < numRows; i++)
        {
          residuals[i] = residualValues[i];
        }

        /* Run info - needed logging turned on */
        HYPRE_PCGGetNumIterations(solver, &num_iterations);
        HYPRE_PCGGetFinalRelativeResidualNorm(solver, &relativeResidualNorm);

        /* Destroy solver and preconditioner */
        HYPRE_ParCSRPCGDestroy(solver);
        HYPRE_BoomerAMGDestroy(precond);
      }
      break;
      //GMRES and AMG
    case 1:
      {
        HYPRE_ParCSRMatrix parcsr_A;
        HYPRE_IJMatrixGetObject(mA, (void **)&parcsr_A);

        HYPRE_ParVector par_b;
        HYPRE_IJVectorGetObject(mb, (void **)&par_b);

        HYPRE_ParVector par_x;
        HYPRE_IJVectorGetObject(mx, (void **)&par_x);

        HYPRE_Solver solver;
        HYPRE_Solver precond;

        int num_iterations;

        HYPRE_ParVector par_residualVector;

        HYPRE_ParCSRFlexGMRESCreate(communicator, &solver);

        /* Set some parameters (See Reference Manual for more parameters) */
        //HYPRE_FlexGMRESSetKDim(solver, 30);
        HYPRE_FlexGMRESSetMaxIter(solver, m_solverMaximumNumberOfIterations); /* max iterations */
        HYPRE_FlexGMRESSetTol(solver, m_solverConvergenceTol);                /* conv. tolerance */
        HYPRE_FlexGMRESSetPrintLevel(solver, 0);                              /* print solve info */
        HYPRE_FlexGMRESSetLogging(solver, 1);                                 /* needed to get run info later */

        /* Now set up the AMG preconditioner and specify any parameters */
        HYPRE_BoomerAMGCreate(&precond);
        HYPRE_BoomerAMGSetPrintLevel(precond, 0); /* print amg solution info */
        HYPRE_BoomerAMGSetCoarsenType(precond, 8);

#ifdef USE_OPENMP
        HYPRE_BoomerAMGSetRelaxType(precond, 5); /* Sym G.S./Jacobi hybrid */
#else
        HYPRE_BoomerAMGSetRelaxType(precond, 6);
#endif

        HYPRE_BoomerAMGSetNumSweeps(precond, 1);
        HYPRE_BoomerAMGSetMaxIter(precond, m_solverAMGPreconditionerNumberOfIterations); /* do only one iteration! */
        HYPRE_BoomerAMGSetTol(precond, 0.0);                                             /* conv. tolerance zero */

        /* Set the FlexGMRES preconditioner */
        HYPRE_FlexGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSolve,
                                  (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSetup, precond);

        /* Now setup and solve! */
        HYPRE_ParCSRFlexGMRESSetup(solver, parcsr_A, par_b, par_x);

        result = HYPRE_ParCSRFlexGMRESSolve(solver, parcsr_A, par_b, par_x);


        /*Get number of iterations used*/
        HYPRE_FlexGMRESGetNumIterations(solver, &numIterations);

        /*Now copy solution! */
        HYPRE_IJVectorGetValues(mx, numRows, rowIndexes, x);

        /*Copy residuals! */
        HYPRE_FlexGMRESGetResidual(solver, (void **)&par_residualVector);
        double *residualValues = hypre_VectorData(hypre_ParVectorLocalVector(par_residualVector));

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < numRows; i++)
        {
          residuals[i] = residualValues[i];
        }

        /* Run info - needed logging turned on */
        HYPRE_FlexGMRESGetNumIterations(solver, &num_iterations);
        HYPRE_FlexGMRESGetFinalRelativeResidualNorm(solver, &relativeResidualNorm);

        /* Destroy solver and preconditioner */
        HYPRE_ParCSRFlexGMRESDestroy(solver);
        HYPRE_BoomerAMGDestroy(precond);
      }
      break;
      //GMRES and AMG
    default:
      {

        HYPRE_ParCSRMatrix parcsr_A;
        HYPRE_IJMatrixGetObject(mA, (void **)&parcsr_A);

        HYPRE_ParVector par_b;
        HYPRE_IJVectorGetObject(mb, (void **)&par_b);

        HYPRE_ParVector par_x;
        HYPRE_IJVectorGetObject(mx, (void **)&par_x);

        HYPRE_Solver solver;
        HYPRE_Solver precond;

        int num_iterations;

        HYPRE_ParVector par_residualVector;

        HYPRE_ParCSRFlexGMRESCreate(communicator, &solver);

        /* Set some parameters (See Reference Manual for more parameters) */
        HYPRE_FlexGMRESSetMaxIter(solver, m_solverMaximumNumberOfIterations); /* max iterations */
        HYPRE_FlexGMRESSetTol(solver, m_solverConvergenceTol);                /* conv. tolerance */
        HYPRE_FlexGMRESSetPrintLevel(solver, 0);                              /* print solve info */
        HYPRE_FlexGMRESSetLogging(solver, 1);                                 /* needed to get run info later */

        /* Now set up the AMG preconditioner and specify any parameters */
        HYPRE_BoomerAMGCreate(&precond);
        HYPRE_BoomerAMGSetPrintLevel(precond, 0); /* print amg solution info */
        HYPRE_BoomerAMGSetCoarsenType(precond, 8);

#ifdef USE_OPENMP
        HYPRE_BoomerAMGSetRelaxType(precond, 5); /* Sym G.S./Jacobi hybrid */
#else
        HYPRE_BoomerAMGSetRelaxType(precond, 6);
#endif
        HYPRE_BoomerAMGSetNumSweeps(precond, 1);
        HYPRE_BoomerAMGSetMaxIter(precond, m_solverAMGPreconditionerNumberOfIterations); /* do only one iteration! */
        HYPRE_BoomerAMGSetTol(precond, 0.0);                                             /* conv. tolerance zero */

        /* Set the FlexGMRES preconditioner */
        HYPRE_FlexGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSolve,
                                  (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSetup, precond);

        /* Now setup and solve! */
        HYPRE_ParCSRFlexGMRESSetup(solver, parcsr_A, par_b, par_x);

        result = HYPRE_ParCSRFlexGMRESSolve(solver, parcsr_A, par_b, par_x);

        /*Get number of iterations used*/
        HYPRE_FlexGMRESGetNumIterations(solver, &numIterations);

        /*Now copy solution! */
        //    HYPRE_IJVectorPrint(mx,"/Users/calebbuahin/Documents/Projects/HydroCouple/FVHMComponent/examples/basic_tests/sloping/XVector.txt");
        HYPRE_IJVectorGetValues(mx, numRows, rowIndexes, x);

        /*Copy residuals! */
        HYPRE_FlexGMRESGetResidual(solver, (void **)&par_residualVector);
        double *residualValues = hypre_VectorData(hypre_ParVectorLocalVector(par_residualVector));

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < numRows; i++)
        {
          residuals[i] = residualValues[i];
        }

        /* Run info - needed logging turned on */
        HYPRE_FlexGMRESGetNumIterations(solver, &num_iterations);
        HYPRE_FlexGMRESGetFinalRelativeResidualNorm(solver, &relativeResidualNorm);

        /* Destroy solver and preconditioner */
        HYPRE_ParCSRFlexGMRESDestroy(solver);
        HYPRE_BoomerAMGDestroy(precond);
      }
      break;
  }

  switch (result)
  {
    case 1:
    case 2:
    case 3:
    case 256:
      {
        if (m_outputDir.exists())
        {
          HYPRE_IJMatrixPrint(mA, qPrintable(m_outputDir.absolutePath() + "/AMatrix_" + QString::number(mpiProcessRank()) + ".txt"));
          HYPRE_IJVectorPrint(mb, qPrintable(m_outputDir.absolutePath() + "/BVector_" + QString::number(mpiProcessRank()) + ".txt"));
          HYPRE_IJVectorPrint(mx, qPrintable(m_outputDir.absolutePath() + "/XVector_" + QString::number(mpiProcessRank()) + ".txt"));
        }
      }
      break;
  }

  HYPRE_IJMatrixDestroy(mA);
  HYPRE_IJVectorDestroy(mb);
  HYPRE_IJVectorDestroy(mx);

  delete[] colsPerRow;
  delete[] rowIndexes;
  delete[] colIndexes;
  delete[] values;

  //  printf("exiting solver \n");

  return result;
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

#ifdef USE_OPENMP
#pragma omp atomic
#endif
    value += v * v;
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
