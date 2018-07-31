#include "stdafx.h"
#include "fvhmcomponent.h"
#include "fvhmcomponentinfo.h"
#include "controlvolume.h"
#include "core/idbasedargument.h"
#include "core/valuedefinition.h"
#include "core/idbasedinputs.h"
#include "core/componentstatuschangeeventargs.h"
#include "core/unit.h"
#include "core/idbasedargument.h"
#include "core/dimension.h"
#include "spatial/octree.h"
#include "spatial/tinargument.h"
#include "spatial/polyhedralsurface.h"
#include "spatial/polygon.h"
#include "spatial/edge.h"
#include "spatial/geometryfactory.h"
#include "temporal/timeseriesargument.h"
#include "temporal/timedata.h"
#include "inletoutletflowbc.h"
#include "wsebc.h"
#include "initialwsebc.h"
#include "outletwseslope.h"
#include "precipbc.h"
#include "cvwseoutput.h"
#include "cvdepthoutput.h"
#include "inflowinput.h"
#include "criticaldepthoutflowbc.h"
#include "spatial/geometryexchangeitems.h"
#include "core/abstractadaptedoutput.h"
#include "progresschecker.h"
#include "edgefluxesio.h"

#include <QDebug>
#include <QDir>
#include <QDateTime>
#include <math.h>

#include "netcdf"


#include "HYPRE.h"
#include "HYPRE_utilities.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_krylov.h"
#include "HYPRE_seq_mv.h"
#include "HYPRE_parcsr_mv.h"
#include "_hypre_parcsr_mv.h"

using namespace netCDF;
using namespace netCDF::exceptions;
using namespace std;
using namespace HydroCouple;
using namespace HydroCouple::Spatial;
using namespace HydroCouple::Temporal;
using namespace HydroCouple::SpatioTemporal;


FVHMComponent::FVHMComponent(const QString &id, FVHMComponentInfo* componentInfo)
  :AbstractTimeModelComponent(id,componentInfo),
    m_TINMeshArgument(nullptr)
{

  performStep = &FVHMComponent::performSimpleTimeStep;
  m_epsilon = std::numeric_limits<double>::epsilon();
  m_octree = new Octree(Octree::Octree2D, Octree::AlongEnvelopes,10,1000);
  createDimensions();
  createArguments();
}

FVHMComponent::FVHMComponent(const QString &id, const QString &caption, FVHMComponentInfo* componentInfo)
  :AbstractTimeModelComponent(id, caption, componentInfo),
    m_TINMeshArgument(nullptr)
{
  performStep = &FVHMComponent::performSimpleTimeStep;
  m_epsilon = std::numeric_limits<double>::epsilon();
  m_octree = new Octree(Octree::Octree2D, Octree::AlongEnvelopes,10,1000);
  createDimensions();
  createArguments();
}

FVHMComponent::~FVHMComponent()
{
  delete m_octree;
  deleteControlVolumes();

  if(m_outputNetCDF)
  {
    m_outputNetCDF->sync();
    m_outputNetCDF->close();
    delete m_outputNetCDF;
    m_outputNetCDF = nullptr;
  }
}

QList<QString> FVHMComponent::validate()
{
  //check if mesh available

  //check time step options

  //check

  return QList<QString>();
}

void FVHMComponent::prepare()
{
  if(!isPrepared())
  {

    m_timeStepCount = 0;
    m_convergedCount = 0;

    if(mpiProcessRank() == 0)
    {
      initializeAdaptedOutputs();

      applyInitialBoundaryConditions();

      applyBoundaryConditions(m_startDateTime);

      updateOutputValues(QList<HydroCouple::IOutput*>());

      setWetAndContCells();

#ifdef USE_OPENMP
#pragma omp parallel sections
#endif
      {

#ifdef USE_OPENMP
#pragma omp section
#endif
        {
          calculateCellWSEGradients();
        }

#ifdef USE_OPENMP
#pragma omp section
#endif
        {
          calculateCellVelocityGradients();
          calculateCellEddyViscosities();
        }
      }

      writeOutputs();

      m_nextOutputTime = m_nextOutputTime + (m_outputTimeStep / 1440.00);

      m_initTimeStepCycle = 0;
      m_printFrequencyCounter = 0;
    }

    setPrepared(true);
    setStatus(IModelComponent::Updated,"");
  }
}

void FVHMComponent::update(const QList<HydroCouple::IOutput*> &requiredOutputs)
{
  if(mpiProcessRank() == 0)
  {
    if(status() == IModelComponent::Updated)
    {
      setStatus(IModelComponent::Updating);

      m_printFrequencyCounter++;

      //get minimum time being requested from outputs
      double outputDataItemsMinTime = getMinOutputTime(requiredOutputs);

      //check if iteration or move to next time step
      bool moveToNextTimeStep = outputDataItemsMinTime <= m_currentDateTime &&
                                m_currentDateTime != m_startDateTime ? false : true;

      if(moveToNextTimeStep)
      {
        //calculate next time step;
        m_timeStep = getNextTimeStep();

        //update current date time
        m_currentDateTime = m_currentDateTime + m_timeStep / 86400.0 ;

        m_timeStepCount++;

        //apply new boundary conditions for  next next timestep/current if iteration
        applyBoundaryConditions(m_currentDateTime);

        currentDateTimeInternal()->setJulianDay(m_currentDateTime);
      }
      //or perform iteration
      else
      {
        moveToNextTimeStep = false;
      }


      m_qtDateTime = SDKTemporal::DateTime(m_currentDateTime).dateTime();

      //retrive input boundary conditions from exchangedItems overrides internal boundaries
      applyInputValues();

      setWetAndContCells();

#ifdef USE_OPENMP
#pragma omp parallel sections
#endif
      {

#ifdef USE_OPENMP
#pragma omp section
#endif
        {
          calculateCellWSEGradients();
        }

#ifdef USE_OPENMP
#pragma omp section
#endif
        {
          calculateCellVelocityGradients();
          calculateCellEddyViscosities();
        }
      }

      //performtimeStep
      int iterations = 0;

      QString message = "";

      //residuals
      double uvelRelResidualNormFin = 0.0, vvelRelResidualNormFin = 0.0,
          pressureRelRisdualNormFin = 0.0, continuityRelResidualNormFin = 0.0,
          uvelRelResidualNormInit = 0.0, vvelRelResidualNormInit = 0.0,
          pressureRelRisdualNormInit = 0.0, continuityRelResidualNormInit = 0.0;

      int minUVelSolvIters = 100000000, maxUVelSolvIters = 0, minVVelSolvIters = 100000000,
          maxVVelSolvIters = 0, minPressSolvIters = 100000000, maxPressSolvIters = 0;

      bool converged = false;

      int error = ErrorCode::NoError;


      error = (this->*performStep) (m_timeStep, iterations,
                                    uvelRelResidualNormInit, uvelRelResidualNormFin, minUVelSolvIters, maxUVelSolvIters,
                                    vvelRelResidualNormInit, vvelRelResidualNormFin, minVVelSolvIters, maxVVelSolvIters,
                                    pressureRelRisdualNormInit, pressureRelRisdualNormFin, minPressSolvIters, maxPressSolvIters,
                                    continuityRelResidualNormInit, continuityRelResidualNormFin,
                                    converged, message);



      double minU, maxU, minV, maxV, minH, maxH, minZ, maxZ;

      //copy new values to previous
      if(moveToNextTimeStep)
      {
#ifdef USE_OPENMP
#pragma omp parallel sections
#endif
        {

#ifdef USE_OPENMP
#pragma omp section
#endif
          {
            prepareForNextTimeStep(minU, maxU, minV, maxV, minH, maxH, minZ, maxZ);
          }

#ifdef USE_OPENMP
#pragma omp section
#endif
          {
            updateExternalInflowOutflowTotals(m_timeStep);
          }
        }
      }

      if(m_printFrequencyCounter >= m_printFrequency || error)
      {
        if(m_verbose || error)
        {
          printf("\nResiduals\n"
                 "dt: %s \t dt: %f \t ts: %f\n"
                 "u_min: %f \t u_max: %f \t u_i: %e \t u_f: %e \t u_min_iters: %i \t u_max_iters: %i\n"
                 "v_min: %f \t v_max: %f \t v_i: %e \t v_f: %e \t v_min_iters: %i \t v_max_iters: %i\n"
                 "h_min: %f \t h_max: %f \t p_i: %e \t p_f: %e \t p_min_iters: %i \t p_max_iters: %i\n"
                 "z_min: %f \t z_max: %f \t c_i: %e \t c_f: %e \n"
                 "c_cells: %i/%i \t m_cells: %i/%i \t iters: %i/%i \t conv?: %s\n\n",
                 qPrintable(m_qtDateTime.toString(Qt::ISODate)), m_currentDateTime, m_timeStep,
                 minU, maxU, uvelRelResidualNormInit, uvelRelResidualNormFin, minUVelSolvIters, maxUVelSolvIters,
                 minV, maxV, vvelRelResidualNormInit, vvelRelResidualNormFin, minVVelSolvIters, maxVVelSolvIters,
                 minH, maxH, pressureRelRisdualNormInit, pressureRelRisdualNormFin, minPressSolvIters, maxPressSolvIters,
                 minZ, maxZ, continuityRelResidualNormInit, continuityRelResidualNormFin,
                 m_numContCells, m_numCells , m_numWetCells, m_numCells,  iterations, m_itersPerTimeStep, converged ? "true" : "false");
        }

        m_printFrequencyCounter = 0;
      }

      switch (error)
      {
        case ErrorCode::SolverFailedToConverge:
          {
            setStatus(IModelComponent::Updated , "FVHM simulation with component id " + id() + " solver failed to converge!",
                      (int)((m_currentDateTime - m_startDateTime) * 100.0 /(m_endDateTime - m_startDateTime)));
          }
          break;
        case ErrorCode::FailedToConverge:
          {
            setStatus(IModelComponent::Updated , "Did not converge for current time step.",
                      (int)((m_currentDateTime - m_startDateTime) * 100.0 /(m_endDateTime - m_startDateTime)));
          }
          break;
        case ErrorCode::CriticalFailure:
          {
            setStatus(IModelComponent::Failed , "Critical failure. Try smaller timestep or relaxation factors.",
                      (int)((m_currentDateTime - m_startDateTime) * 100.0 /(m_endDateTime - m_startDateTime)));
            return;
          }
        case ErrorCode::NoError:
          {
            if(moveToNextTimeStep)
            {
              m_convergedCount++;
            }
          }
          break;
      }

      //logging
      {
        writeToLogFile(iterations,
                       uvelRelResidualNormInit,uvelRelResidualNormFin,
                       vvelRelResidualNormInit,vvelRelResidualNormFin,
                       continuityRelResidualNormInit,continuityRelResidualNormFin,
                       pressureRelRisdualNormInit,pressureRelRisdualNormFin);
      }

      updateOutputValues(QList<IOutput*>({}));

      //write outputs
      if(m_currentDateTime != m_previouslyWrittenDateTime && m_currentDateTime >= m_nextOutputTime)
      {

        m_qtDateTime = SDKTemporal::DateTime(m_nextOutputTime).dateTime();
        m_previouslyWrittenDateTime = m_currentDateTime;
        setStatus(IModelComponent::Updating, "Writing output for " + m_qtDateTime.toString(Qt::ISODate) + " ...", progressChecker()->progress());

        writeOutputs();

        m_nextOutputTime = m_nextOutputTime + m_outputTimeStep / 1440.00;
      }


      //check if simulation has been completed
      if(m_currentDateTime >=  m_endDateTime)
      {
        setStatus(IModelComponent::Done , "Simulation Finished successfully!",100);
      }
      else if(progressChecker()->performStep(m_currentDateTime))
      {
        setStatus(IModelComponent::Updated, "Current DateTime: " + m_qtDateTime.toString(Qt::ISODate), progressChecker()->progress());
      }
      else
      {
        setStatus(IModelComponent::Updated);
      }
    }
  }
#ifdef USE_MPI
  else
  {
    if(status() == IModelComponent::Updated)
    {
      setStatus(IModelComponent::Updating);

      MPI_Status status;
      int result = 0;

      if((result = MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status)) == MPI_SUCCESS)
      {
        switch (status.MPI_TAG)
        {
          case 1001:
            {

              int dataSize = 0;
              MPI_Get_count(&status, MPI_DOUBLE, &dataSize);

              if(dataSize)
              {
                double *values = new double[dataSize];
                MPI_Recv(values, dataSize, MPI_DOUBLE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);

                double *b = nullptr;
                SparseMatrix *A = nullptr;
                double *x = nullptr;
                int counter  = 0;
                SparseMatrix::deserialize(values, A, b, x, counter);

                int rowCount = A->rowCount();
                double *residuals = new double[rowCount];
                double residualNorm = 0;
                int numIterations = 0;
                //                printf("Solving on slave\n");
                int result = solve(mpiCommunicator(), *A, A->ilower(), A->iupper(), b, x, residuals, residualNorm, numIterations);

                int size = rowCount * 2 + 2;
                double* solution = new double[size];
                solution[0] = result;
                solution[1] = residualNorm;

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
                for(int i = 0; i < rowCount; i++)
                {
                  solution[i + 2] = x[i];
                }

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
                for(int i = 0; i < rowCount; i++)
                {
                  solution[rowCount + 2 + i] = residuals[i];
                }

                MPI_Send(&solution[0], size, MPI_DOUBLE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD);

                delete[] values;
                delete[] x;
                delete[] residuals;
                delete[] b;
                delete[] solution;
                delete A;
              }
            }
            break;
        }
      }

      setStatus(IModelComponent::Updated);
    }
  }
#endif
}

void FVHMComponent::finish()
{
  if(isPrepared())
  {
    setStatus(IModelComponent::Finishing , "");

    if(mpiProcessRank() == 0)
    {
      if(m_outputNetCDF)
      {
        m_outputNetCDF->sync();
        m_outputNetCDF->close();
        delete m_outputNetCDF;
        m_outputNetCDF = nullptr;
      }

      closeOutputCSVFile();
      closeLogFile();
    }

    setPrepared(false);
    setInitialized(false);

    setStatus(IModelComponent::Finished , "");
    setStatus(IModelComponent::Created , "");
  }
}

double FVHMComponent::startDateTime() const
{
  return m_startDateTime;
}

double FVHMComponent::endDateTime() const
{
  return m_endDateTime;
}

void FVHMComponent::initializeFailureCleanUp()
{
  closeLogFile();
}

void FVHMComponent::createDimensions()
{
  m_patchDimension = new Dimension("PatchDimension", "Dimension for control volume patches",this);
  m_edgeDimension = new Dimension("EdgeDimension","Dimension for control volume edges",this);
  m_nodeDimension = new Dimension("NodeDimension","Dimension for control volume nodes",this);
  m_idDimension = new Dimension("IdDimension","Dimension for identifiers",this);
  m_timeDimension = new Dimension("TimeSeriesDimension","Dimension for times", this);
  m_geometryDimension = new Dimension("GeometryDimension","Dimension for geometries", this);
}

void FVHMComponent::createArguments()
{
  createMeshArguments();
  createTimeArguments();
  createInputFilesArguments();
  createOutputFilesArguments();
  createEdgeFluxesIOArguments();
  createConvergenceArguments();
  createSolverOptionsArguments();
  createAdvectionDiffusionSchemeArguments();
  createPressureVelCouplingArguments();
  createWettingAndDryingArguments();
  createInitialUniformBCArguments();
  createMiscOptionsArguments();
  createInletOutletFlowBCArguments();
  createWSEBCArguments();
  createCriticalDepthOutflowBCArguments();
  createInitialWSEBCArguments();
  createWSESlopeArguments();
  createSourceSinkFlowsBCArguments();
  createPrecipitationArguments();
}

void FVHMComponent::createMeshArguments()
{
  Quantity *TINQuantity = Quantity::unitLessValues("TINValuesQuantity",QVariant::Double,this);

  m_TINMeshArgument = new TINArgumentDouble("TINMesh", HydroCouple::Spatial::Centroid,
                                            m_patchDimension,m_edgeDimension,m_nodeDimension,TINQuantity,this);

  m_TINMeshArgument->setCaption("Computational TIN Mesh");
  m_TINMeshArgument->setDescription("<h1>Computational TIN Mesh</h1>"
                                    "</hr>"
                                    "<p>"
                                    "The Triangular Irregular Network (TIN) mesh used for the computations"
                                    "</p>");

  addArgument(m_TINMeshArgument);
}

void FVHMComponent::createTimeArguments()
{
  QStringList identifiers;
  identifiers.append("StartDateTime");
  identifiers.append("EndDateTime");

  Quantity* quantity = Quantity::unitLessValues("DateTimeQuantity", QVariant::DateTime , this);

  m_simulationDurationArgument = new IdBasedArgumentDateTime("SimulationDuration",identifiers,m_idDimension,quantity,this);
  m_simulationDurationArgument->setCaption("Simulation Duration");
  m_simulationDurationArgument->setMatchIdentifiersWhenReading(true);

  QDate date = QDate::currentDate();

  (*m_simulationDurationArgument)["StartDateTime"] = QDateTime(date);
  (*m_simulationDurationArgument)["EndDateTime"] = QDateTime(date.addDays(1));

  m_simulationDurationArgument->setDescription("<h1>Simulation Duration</h1>"
                                               "<hr>"
                                               "<p>"
                                               "The start and end times for the current simulation."
                                               "</p>");

  addArgument(m_simulationDurationArgument);

  QStringList timeStepIdentifiers;
  timeStepIdentifiers.append("Time Step Mode");
  timeStepIdentifiers.append("Min Time Step");
  timeStepIdentifiers.append("Max Time Step");
  timeStepIdentifiers.append("Max Time Step Change Factor");
  timeStepIdentifiers.append("Max Time Decrease Step Change Factor");
  timeStepIdentifiers.append("Max Courant Number");
  timeStepIdentifiers.append("RMS Courant Number");
  timeStepIdentifiers.append("Time Step Relaxation Factor");
  timeStepIdentifiers.append("Output Time Step");
  timeStepIdentifiers.append("Number of Fixed Time Steps");

  Quantity* timeStepQuantity = Quantity::unitLessValues("TimeStepQuantity", QVariant::Double , this);
  m_timeStepArguments = new IdBasedArgumentDouble("TimeStepOptions",timeStepIdentifiers, m_idDimension, timeStepQuantity , this);
  m_timeStepArguments->setCaption("Time Step Options");
  m_timeStepArguments->setMatchIdentifiersWhenReading(true);

  (*m_timeStepArguments)["Time Step Mode"] = m_useAdaptiveTimeStep ? (m_adaptiveTSMode == MaxCourantNumber ? 1 : 2) : 0;
  (*m_timeStepArguments)["Min Time Step"] = m_minTimeStep;
  (*m_timeStepArguments)["Max Time Step"] = m_maxTimeStep;
  (*m_timeStepArguments)["Max Time Step Change Factor"] = m_maxTimeStepIncreaseCF;
  (*m_timeStepArguments)["Max Time Step Decrease Change Factor"] = m_maxTimeStepDecreaseCF;
  (*m_timeStepArguments)["Max Courant Number"] = m_maxCourantNumber;
  (*m_timeStepArguments)["RMS Courant Number"] = m_RMSCourantNumber;
  (*m_timeStepArguments)["Time Step Relaxation Factor"] = m_timeStepRelaxFactor;
  (*m_timeStepArguments)["Output Time Step"] = m_outputTimeStep;
  (*m_timeStepArguments)["Number of Fixed Time Steps"] = m_numFixedTimeSteps;


  QString description = "<div>"
                        "<h1>Time Step Options</h1>"
                        "<hr>"
                        "<h2>Time Step Mode</h2>"
                        "<p>"
                        "<dl>"
                        "<dt>0</dt>"
                        "<dd><p>Constant time step</p></dd>"
                        "<dt>1</dt>"
                        "<dd><p>Adaptive time step using maximum courant number for all cells. If time step "
                        "calculated is less than specified time step, the specified time step is used</p></dd>"
                        "<dt>2</dt>"
                        "<dd><p>Adaptive time step using root mean square of the courant number for all cells. If time step "
                        "calculated is less than specified time step, the specified time step is used</p></dd>"
                        "</dl>"
                        "</p>"
                        "<h1>Time Step</h1>"
                        "<hr>"
                        "<p>Time step specified in seconds</p>"
                        "<h1>Time Step Relaxation Factor</h1>"
                        "<hr>"
                        "<p>Time step relaxation factor for adaptive time-step</p>"
                        "<h1>Output Time Step</h1>"
                        "<hr>"
                        "<p>Output Time step specified in minutes</p>"
                        "</div>";

  m_timeStepArguments->setDescription(description);
  addArgument(m_timeStepArguments);
}

void FVHMComponent::createInputFilesArguments()
{
  Quantity* quantity = Quantity::unitLessValues("OutputFileQuality", QVariant::String , this);
  quantity->setDefaultValue(QString(""));
  quantity->setMissingValue(QString(""));

  QStringList options;
  options << "Restart File";

  m_inputFiles = new IdBasedArgumentString("Input Files", options, m_idDimension, quantity, this);
  m_inputFiles->setCaption("Input Files");
  m_inputFiles->setMatchIdentifiersWhenReading(true);
  addArgument(m_inputFiles);
}

void FVHMComponent::createOutputFilesArguments()
{
  Quantity* quantity = Quantity::unitLessValues("OutputFileQuality", QVariant::String , this);
  quantity->setDefaultValue(QString(""));
  quantity->setMissingValue(QString(""));

  QStringList options;
  options << "Output NetCDF File";
  options << "Output Shapefile";
  options << "Output Edge Fluxes File";
  options << "WriteActiveCellsOnly";
  options << "Log File";

  m_outputFilesArgument = new IdBasedArgumentString("Output Options", options, m_idDimension,quantity,this);
  m_outputFilesArgument->setCaption("Output Options");
  m_outputFilesArgument->setMatchIdentifiersWhenReading(true);

  addArgument(m_outputFilesArgument);
}

void FVHMComponent::createEdgeFluxesIOArguments()
{
  Quantity* quantity = Quantity::flowInCMS(this);
  quantity->setCaption("Outflow (cms)");
  quantity->setMissingValue(QString(""));
  quantity->setDefaultValue(QString(""));

  m_edgeFluxesIO = new EdgeFluxesIO("Output Edge Fluxes",
                                    m_geometryDimension,
                                    quantity,this);

  m_edgeFluxesIO->setCaption("Output Edge Fluxes");
  m_edgeFluxesIO->setIsOptional(true);
  addArgument(m_edgeFluxesIO);
}

void FVHMComponent::createConvergenceArguments()
{
  QStringList identifiers;
  identifiers.append("Pressure Rel. Norm Criteria");
  identifiers.append("Continuity Rel. Norm Criteria");
  identifiers.append("X-Velocity Rel. Norm Criteria");
  identifiers.append("Y-Velocity Rel. Norm Criteria");

  Quantity* quantity = Quantity::unitLessValues("ConvergenceQuantity", QVariant::Double , this);
  quantity->setMissingValue(-99999999);

  m_convergenceArguments = new IdBasedArgumentDouble("Convergence Options",identifiers,m_idDimension,quantity,this);
  m_convergenceArguments->setCaption("Convergence Options");
  m_convergenceArguments->setMatchIdentifiersWhenReading(true);



  (*m_convergenceArguments)["Pressure Rel. Norm Criteria"] = m_pressureConvergenceTol;
  (*m_convergenceArguments)["Continuity Rel. Norm Criteria"] = m_contConvergenceTol;
  (*m_convergenceArguments)["X-Velocity Rel. Norm Criteria"] = m_uvelConvergenceTol;
  (*m_convergenceArguments)["Y-Velocity Rel. Norm Criteria"] = m_vvelConvergenceTol;

  m_convergenceArguments->setDescription("<h1>Convergence Options</h1>"
                                         "<hr>"
                                         "<p>"
                                         "Convergence relative residual norm criteria  for convergence."
                                         "</p>");

  addArgument(m_convergenceArguments);
}

void FVHMComponent::createSolverOptionsArguments()
{
  Dimension *solverOptionsDimension = new Dimension("SolverOptionDimension","Dimension for solver option identifiers", this);
  QStringList identifiers;
  identifiers.append("Solver Type");
  identifiers.append("Solver Max No. Iterations");
  identifiers.append("AMG Max No.Iterations");
  identifiers.append("Convergence Criteria");

  Quantity* quantity = Quantity::unitLessValues("ConvergenceQuantity", QVariant::Double , this);

  m_solverOptionsArguments = new IdBasedArgumentDouble("Solver Options",identifiers,solverOptionsDimension,quantity,this);
  m_solverOptionsArguments->setCaption("Solver Options");
  m_solverOptionsArguments->setMatchIdentifiersWhenReading(true);

  (*m_solverOptionsArguments)["Solver Type"] = m_solverType;
  (*m_solverOptionsArguments)["Solver Max No. Iterations"] = m_solverMaximumNumberOfIterations;
  (*m_solverOptionsArguments)["AMG Max No.Iterations"] = m_solverAMGPreconditionerNumberOfIterations;
  (*m_solverOptionsArguments)["Convergence Criteria"] = m_solverConvergenceTol;

  addArgument(m_solverOptionsArguments);
}

void FVHMComponent::createPressureVelCouplingArguments()
{
  QStringList identifiers;
  identifiers.append("Coupling Type");
  identifiers.append("Iterations per TimeStep");
  identifiers.append("Pressure Relaxation Factor");
  identifiers.append("Velocity Relaxation Factor");

  Quantity* quantity = Quantity::unitLessValues("PressureVelocityCoupleQuantity", QVariant::Double , this);
  quantity->setMissingValue(-99999999);

  m_pressureVelCouplingArguments = new IdBasedArgumentDouble("Pressure-Velocity Coupling Options",
                                                             identifiers,m_idDimension,quantity,this);
  m_pressureVelCouplingArguments->setCaption("Pressure-Velocity Coupling Options");
  m_pressureVelCouplingArguments->setMatchIdentifiersWhenReading(true);

  (*m_pressureVelCouplingArguments)["Coupling Type"] = m_pressureVelocityCouplingType;
  (*m_pressureVelCouplingArguments)["Iterations per TimeStep"] = m_itersPerTimeStep;
  (*m_pressureVelCouplingArguments)["Pressure Relaxation Factor"] = m_pressureRelaxFactor;
  (*m_pressureVelCouplingArguments)["Velocity Relaxation Factor"] = m_velocityRelaxFactor;

  addArgument(m_pressureVelCouplingArguments);
}

void FVHMComponent::createWettingAndDryingArguments()
{
  QStringList identifiers;
  identifiers.append("Wet Cell Depth");
  identifiers.append("Wet Cell Node Depth");
  identifiers.append("Minimum Viscous Sublayer Depth");
  identifiers.append("Wet Dry Factor");

  Quantity* quantity = Quantity::lengthInMeters(this);
  quantity->setCaption("Wetting and Drying Depths (m)");
  quantity->setMissingValue(-99999999);

  m_wettingAndDryingArgument = new IdBasedArgumentDouble("Wetting and Drying Parameters",identifiers,m_idDimension,quantity,this);
  m_wettingAndDryingArgument->setCaption("Wetting and Drying Parameters");
  m_wettingAndDryingArgument->setMatchIdentifiersWhenReading(true);


  (*m_wettingAndDryingArgument)["Wet Cell Depth"] = m_wetCellDepth;
  (*m_wettingAndDryingArgument)["Wet Cell Node Depth"] = m_wetCellNodeDepth;
  (*m_wettingAndDryingArgument)["Minimum Viscous Sublayer Depth"] = m_viscousSubLayerDepth;
  (*m_wettingAndDryingArgument)["Wet Dry Factor"] = m_wetDryFactor;

  addArgument(m_wettingAndDryingArgument);
}

void FVHMComponent::createInitialUniformBCArguments()
{
  QStringList identifiers;

  identifiers.append("Uniform WSE");
  identifiers.append("Uniform Depth");
  identifiers.append("Uniform U-Velocity");
  identifiers.append("Uniform V-Velocity");
  identifiers.append("Uniform Mannings");
  identifiers.append("Viscosity");
  identifiers.append("Min Outlet WSE Slope");
  identifiers.append("Inlet WSEBC Vel Relaxation Factor");
  identifiers.append("Outlet WSEBC Vel Relaxation Factor");

  Quantity* quantity = Quantity::unitLessValues("Initial Model Variables",QVariant::Double, this);
  quantity->setDefaultValue(0.0);
  quantity->setMissingValue(-99999999);

  m_initialParamArguments = new IdBasedArgumentDouble("Initial Model Variables",identifiers,m_idDimension,quantity,this);
  m_initialParamArguments->setCaption("Initial Model Variables");
  m_initialParamArguments->setMatchIdentifiersWhenReading(true);

  (*m_initialParamArguments)["Uniform WSE"] = m_uniWSE;
  (*m_initialParamArguments)["Uniform Depth"] = m_uniDepth;
  (*m_initialParamArguments)["Uniform U-Velocity"] = m_uniUVel;
  (*m_initialParamArguments)["Uniform V-Velocity"] = m_uniVVel;
  (*m_initialParamArguments)["Uniform Mannings"] = m_uniMannings;
  (*m_initialParamArguments)["Viscosity"] = m_viscosity;
  (*m_initialParamArguments)["Min Outlet WSE Slope"] = m_inletOutletWSESlope;
  (*m_initialParamArguments)["Inlet WSEBC Vel Relaxation Factor"] = m_inletWSEBCVRelFactor;
  (*m_initialParamArguments)["Outlet WSEBC Vel Relaxation Factor"] = m_outletWSEBCVRelFactor;

  addArgument(m_initialParamArguments);
}

void FVHMComponent::createInletOutletFlowBCArguments()
{

  Quantity* quantity = Quantity::flowInCMS(this);
  quantity->setMissingValue(-99999999);

  m_inletFlowBCArgument = new InletFlowBCArgument("Inlet Flow BC", m_timeDimension,
                                                  m_geometryDimension,quantity,this);
  m_inletFlowBCArgument->setCaption("Inlet Flow BC");
  m_inletFlowBCArgument->setIsOptional(true);

  addArgument(m_inletFlowBCArgument);


  m_outletFlowBCArgument = new OutletFlowBCArgument("Outlet Flow BC", m_timeDimension,
                                                    m_geometryDimension,quantity,this);
  m_outletFlowBCArgument->setCaption("Outlet Flow BC");
  m_outletFlowBCArgument->setIsOptional(true);

  addArgument(m_outletFlowBCArgument);
}

void FVHMComponent::createSourceSinkFlowsBCArguments()
{

  Quantity* quantity = Quantity::flowInCMS(this);
  quantity->setMissingValue(-99999999);

  m_sourceSinkFlowBCArgument = new TimeGeometryArgumentDouble("Sources & Sinks Flow BC",
                                                              IGeometry::Polygon,
                                                              m_timeDimension,
                                                              m_geometryDimension, quantity,this);
  m_sourceSinkFlowBCArgument->setCaption("Sources and Sinks Flow BC");
  m_sourceSinkFlowBCArgument->setIsOptional(true);

  addArgument(m_sourceSinkFlowBCArgument);
}

void FVHMComponent::createWSEBCArguments()
{
  Quantity* quantity = Quantity::lengthInMeters(this);
  quantity->setCaption("Water Surface Elevation (m)");
  quantity->setMissingValue(-99999999);

  m_outletWSEBC = new WSEBC("Outlet WSE BC", WSEBC::Outlet,
                            m_outletWSEBCVRelFactor,
                            m_timeDimension,
                            m_geometryDimension,quantity,this);

  m_outletWSEBC->setCaption("Outlet Water Surface Elevation BC");
  m_outletWSEBC->setIsOptional(true);


  m_inletWSEBC = new WSEBC("Inlet WSE BC",
                           WSEBC::Inlet,
                           m_inletWSEBCVRelFactor,
                           m_timeDimension,
                           m_geometryDimension,quantity,this);

  m_inletWSEBC->setCaption("Inlet Water Surface Elevation BC");
  m_inletWSEBC->setIsOptional(true);

  addArgument(m_outletWSEBC);
  addArgument(m_inletWSEBC);
}

void FVHMComponent::createCriticalDepthOutflowBCArguments()
{
  Quantity* quantity = Quantity::lengthInMeters(this);
  quantity->setCaption("Critical Depth (m)");
  quantity->setMissingValue(-99999999);

  m_criticalDepthOutflowBC = new CriticalDepthOutflowBC("Critical Depth Outflow BC",
                                                        m_geometryDimension,
                                                        quantity,this);

  m_criticalDepthOutflowBC->setCaption("Critical Depth Outflow BC");
  m_criticalDepthOutflowBC->setIsOptional(true);
  addArgument(m_criticalDepthOutflowBC);
}

void FVHMComponent::createInitialWSEBCArguments()
{
  //  Quantity* quantity = Quantity::lengthInMeters(this);
  //  quantity->setCaption("Water Surface Elevation (m)");
  //  quantity->setMissingValue(-99999999);

  m_initialWSEBC = new InitialWSEBC("Initial WSE BC",this);
  m_initialWSEBC->setCaption("Initial Water Surface Elevation BC");
  ValueDefinition *quantity = dynamic_cast<ValueDefinition*>(m_initialWSEBC->valueDefinition());
  quantity->setCaption("Water Surface Elevation (m)");
  quantity->setMissingValue(-99999999);

  addArgument(m_initialWSEBC);
}

void FVHMComponent::createWSESlopeArguments()
{
  Quantity* quantity = Quantity::lengthInMeters(this);
  quantity->setCaption("Water Surface Elevation (m)");
  quantity->setMissingValue(-99999999);

  m_outletWSESlope = new OutletWSESlope("Outlet WSE Slope BC",
                                        m_geometryDimension,
                                        quantity,this);

  m_outletWSESlope->setCaption("Outlet Water Surface Elevation Slope BC");
  m_outletWSESlope->setIsOptional(true);
  addArgument(m_outletWSESlope);
}

void FVHMComponent::createPrecipitationArguments()
{

  //fix
  Quantity* quantity = Quantity::lengthInMeters(this);
  quantity->setCaption("Rainfall Rate (mm/hr)");
  quantity->setMissingValue(-99999999);

  m_precipitationArgument = new PrecipBC("Precipitation Rate BC", m_timeDimension,m_geometryDimension, quantity,this);
  m_precipitationArgument->setCaption("Precipitation Rate BC");
  m_precipitationArgument->setDescription("Rainfall Rate (mm/hr)");
  m_precipitationArgument->setIsOptional(true);

  addArgument(m_precipitationArgument);
}

void FVHMComponent::createAdvectionDiffusionSchemeArguments()
{
  QStringList identifiers;
  identifiers.append("Use Turbulence");
  identifiers.append("Turbulence Scheme");
  identifiers.append("Smagorinsky-Lilly Coefficient");
  identifiers.append("Parabolic Eddy Viscosity Constant");
  identifiers.append("Advection Scheme");

  Quantity* quantity = Quantity::unitLessValues("Indexes for Advec/Diffusion Modes", QVariant::Double, this);


  m_advectionDiffusionSchemeArgument = new IdBasedArgumentDouble("Advection Diffusion Scheme Options",identifiers,m_idDimension,quantity,this);
  m_advectionDiffusionSchemeArgument->setCaption("Advection Diffusion Scheme Options");
  m_advectionDiffusionSchemeArgument->setMatchIdentifiersWhenReading(true);


  (*m_advectionDiffusionSchemeArgument)["Use Turbulence"] = m_useEddyViscosity ? 1.0 : 0.0 ;
  (*m_advectionDiffusionSchemeArgument)["Turbulence Scheme"] = m_turbulenceScheme;
  (*m_advectionDiffusionSchemeArgument)["Smagorinsky-Lilly Coefficient"] = m_smargorinskyCoefficient;
  (*m_advectionDiffusionSchemeArgument)["Parabolic Eddy Viscosity Constant"] = m_parabolicEddyViscosityConstant;
  (*m_advectionDiffusionSchemeArgument)["Advection Scheme"] = m_advectionScheme;

  addArgument(m_advectionDiffusionSchemeArgument);
}

void FVHMComponent::createMiscOptionsArguments()
{
  QStringList identifiers;

  identifiers.append("Verbose");
  identifiers.append("Print Frequency");
  identifiers.append("Output Write Flush Frequency");
  identifiers.append("Gradient Calculation Mode");
  identifiers.append("Use Wall Friction");
  identifiers.append("MPI Solver Split Threshold");

  Quantity* quantity = Quantity::unitLessValues("Misc Options Values", QVariant::Double, this);
  quantity->setDefaultValue(0.0);
  quantity->setMissingValue(-99999999);

  m_miscOptionsArgument = new IdBasedArgumentDouble("Misc Options Argument",identifiers,m_idDimension,quantity,this);
  m_miscOptionsArgument->setCaption("Misc Options");
  m_miscOptionsArgument->setMatchIdentifiersWhenReading(true);

  (*m_miscOptionsArgument)["Verbose"] = m_verbose ? 1.0 : 0.0;
  (*m_miscOptionsArgument)["Print Frequency"] = m_printFrequency;
  (*m_miscOptionsArgument)["Output Write Flush Frequency"] = m_writeFrequency;
  (*m_miscOptionsArgument)["Gradient Calculation Mode"] = TriCV::gradientCalcMode;
  (*m_miscOptionsArgument)["Use Wall Friction"] = m_useWall;
  (*m_miscOptionsArgument)["MPI Solver Split Threshold"] = m_mpiSolverSplitThreshold;

  addArgument(m_miscOptionsArgument);
}

bool FVHMComponent::initializeArguments(QString &message)
{
  if(mpiProcessRank() == 0)
  {
    bool initialized = initializeWettingAndDryingArguments(message) &&
                       initializeMeshArguments(message) &&
                       initializeInitialUniformBCArguments(message) &&
                       initializeTimeArguments(message) &&
                       initializeOutputFilesArguments(message) &&
                       initializeEdgeFluxesIOArguments(message) &&
                       initializeConvergenceArguments(message) &&
                       initializeSolverOptionsArguments(message) &&
                       initializePressureVelCouplingArguments(message) &&
                       initializeInletOutletFlowBCArguments(message) &&
                       initializeSourceSinkFlowsBCArguments(message) &&
                       initializeWSEBCArguments(message) &&
                       initializeCriticalDepthOutflowBCArguments(message) &&
                       initializeInitialWSEBCArguments(message) &&
                       initializeWSESlopeBCArguments(message) &&
                       initializePrecipitationArguments(message) &&
                       initializeAdvectionDiffusionSchemeArguments(message) &&
                       initializeMiscOptionsArguments(message);

    return initialized;
  }
  else
  {
    bool initialized = initializeWettingAndDryingArguments(message) &&
                       //                       initializeMeshArguments(message) &&
                       //                       initializeInitialUniformBCArguments(message) &&
                       //                       initializeTimeArguments(message) &&
                       initializeOutputFilesArguments(message) &&
                       //                       initializeEdgeFluxesIOArguments(message) &&
                       initializeConvergenceArguments(message) &&
                       initializeSolverOptionsArguments(message) &&
                       initializePressureVelCouplingArguments(message) &&
                       //                       initializeInletOutletFlowBCArguments(message) &&
                       //                       initializeSourceSinkFlowsBCArguments(message) &&
                       //                       initializeWSEBCArguments(message) &&
                       //                       initializeCriticalDepthOutflowBCArguments(message) &&
                       //                       initializeInitialWSEBCArguments(message) &&
                       //                       initializeWSESlopeBCArguments(message) &&
                       //                       initializePrecipitationArguments(message) &&
                       //                       initializeAdvectionDiffusionSchemeArguments(message) &&
                       initializeMiscOptionsArguments(message) ;

    return initialized;
  }
}

bool FVHMComponent::initializeMeshArguments(QString &message)
{
  deleteControlVolumes();

  m_sharedTriangles.clear();

  if(m_TINMeshArgument->TINInternal() != nullptr)
  {
    if((m_numCells = m_TINMeshArgument->TINInternal()->patchCount()))
    {
      m_wetCells.clear();
      m_wetCells.reserve(m_numCells);

      m_contCells.clear();
      m_contCells.reserve(m_numCells);

      m_controlVolumes.clear();
      m_controlVolumes.reserve(m_numCells);

      for(int i = 0; i < m_numCells ; i++)
      {
        HCTriangle *triangle = m_TINMeshArgument->TINInternal()->triangleInternal(i);
        triangle->setIndex(i);

        TriCV *cv = new TriCV(triangle, this);

        cv->index = i;
        cv->wetIndex = 1;
        cv->wetCellIndex = i;

        m_controlVolumes.push_back(cv);
        m_wetCells.push_back(i);
        m_contCells.push_back(i);

        m_sharedTriangles.append(QSharedPointer<HCGeometry>(triangle, deleteLater));
      }

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
      for(int i = 0 ; i < m_numCells ; i++)
      {
        m_controlVolumes[i]->calculateAdjacentCellParams();
      }

      resetOctree();
      return true;
    }
  }

  m_numCells = 0;
  message = "Mesh has not been set or does not have any elements";
  resetOctree();

  return false;

}

bool FVHMComponent::initializeTimeArguments(QString &message)
{
  QDateTime startDateTime = (*m_simulationDurationArgument)["StartDateTime"];
  QDateTime endDateTime = (*m_simulationDurationArgument)["EndDateTime"];

  if(startDateTime >= endDateTime)
  {
    message =  "Start datetime must be earlier than end time";
    return false;
  }
  else
  {
    m_startDateTime =  SDKTemporal::DateTime::toJulianDays(startDateTime);
    m_endDateTime = SDKTemporal::DateTime::toJulianDays(endDateTime);
    m_currentDateTime = m_startDateTime;
    progressChecker()->reset(m_startDateTime, m_endDateTime);
  }

  m_minTimeStep = (*m_timeStepArguments)["Min Time Step"];
  m_maxTimeStep = (*m_timeStepArguments)["Max Time Step"];

  if(m_startDateTime + (m_minTimeStep/(60*60*24)) > m_endDateTime)
  {
    message = "Time step must be less than specified simulation duration interval";
    return false;
  }


  m_maxTimeStepIncreaseCF = (*m_timeStepArguments)["Max Time Step Change Factor"];

  if(m_maxTimeStepIncreaseCF < 0)
    m_maxTimeStepIncreaseCF = 0.05;

  m_maxTimeStepDecreaseCF = (*m_timeStepArguments)["Max Time Step Decrease Change Factor"];

  if(m_maxTimeStepDecreaseCF < 0)
    m_maxTimeStepDecreaseCF = 0.1;

  m_outputTimeStep = (*m_timeStepArguments)["Output Time Step"];
  m_previouslyWrittenDateTime = m_currentDateTime;
  m_nextOutputTime = m_currentDateTime;
  m_qtDateTime = SDKTemporal::DateTime(m_currentDateTime).dateTime();

  timeHorizonInternal()->setJulianDay(m_startDateTime);
  timeHorizonInternal()->setDuration(m_endDateTime - m_startDateTime);
  currentDateTimeInternal()->setJulianDay(m_currentDateTime);

  if(m_currentDateTime + (m_outputTimeStep / 1440.00) > m_endDateTime)
  {
    message = "Output time step must be less than specified simulation duration interval";
    return false;
  }

  int timeStepMode = (*m_timeStepArguments)["Time Step Mode"];

  m_numFixedTimeSteps = (*m_timeStepArguments)["Number of Fixed Time Steps"];

  if(m_numFixedTimeSteps < 2)
    m_numFixedTimeSteps = 2;

  switch (timeStepMode)
  {

    case 1:
      {
        m_useAdaptiveTimeStep = true;
        m_adaptiveTSMode = AdaptiveTSMode::MaxCourantNumber;
        m_maxCourantNumber = (*m_timeStepArguments)["Max Courant Number"];

        if(m_maxCourantNumber < 0.0)
          m_maxCourantNumber = 1.0;

      }
      break;
    case 2:
      {
        m_useAdaptiveTimeStep = true;
        m_adaptiveTSMode = AdaptiveTSMode::RMSCourantNumber;
        m_RMSCourantNumber = (*m_timeStepArguments)["RMS Courant Number"];

        if(m_RMSCourantNumber < 0.0)
          m_RMSCourantNumber = 1.0;
      }
      break;
    default:
      {
        m_useAdaptiveTimeStep = false;
      }
      break;
  }

  return true;
}

bool FVHMComponent::initializeConvergenceArguments(QString &message)
{
  message += "";

  m_pressureConvergenceTol = (*m_convergenceArguments)["Pressure Rel. Norm Criteria"];
  m_contConvergenceTol = (*m_convergenceArguments)["Continuity Rel. Norm Criteria"];
  m_uvelConvergenceTol = (*m_convergenceArguments)["X-Velocity Rel. Norm Criteria"];
  m_vvelConvergenceTol = (*m_convergenceArguments)["Y-Velocity Rel. Norm Criteria"];

  return true;
}

bool FVHMComponent::initializeSolverOptionsArguments(QString &message)
{
  message += "";

  m_solverType = (*m_solverOptionsArguments)["Solver Type"];
  m_solverMaximumNumberOfIterations = (*m_solverOptionsArguments)["Solver Max No. Iterations"];
  m_solverAMGPreconditionerNumberOfIterations = (*m_solverOptionsArguments)["AMG Max No.Iterations"];
  m_solverConvergenceTol = (*m_solverOptionsArguments)["Convergence Criteria"];

  return true;
}

bool FVHMComponent::initializePressureVelCouplingArguments(QString &message)
{
  message += "";

  // m_solverType = (*m_pressureVelCouplingArguments)["Coupling Type"];
  m_itersPerTimeStep = (*m_pressureVelCouplingArguments)["Iterations per TimeStep"];

  if(m_itersPerTimeStep < 1)
  {
    message = "Number of iterations for the pressure velocity coupling must be greater than 0";
    return false;
  }

  m_pressureVelocityCouplingType = (*m_pressureVelCouplingArguments)["Coupling Type"];

  switch (m_pressureVelocityCouplingType)
  {
    case 1:
      {
        performStep = &FVHMComponent::performSimpleCTimeStep;
      }
      break;
    case 2:
      {
        performStep = &FVHMComponent::performPISOTimeStep;
      }
      break;
    default:
      {
        performStep = &FVHMComponent::performSimpleTimeStep;
      }
      break;
  }


  m_pressureRelaxFactor = (*m_pressureVelCouplingArguments)["Pressure Relaxation Factor"];

  if(m_pressureRelaxFactor <= 0 || m_pressureRelaxFactor > 1)
  {
    message = "Pressure relaxation factor must be greater than 0 and less than or equal to 1";
    return false;
  }

  m_velocityRelaxFactor = (*m_pressureVelCouplingArguments)["Velocity Relaxation Factor"];

  if(m_velocityRelaxFactor <= 0 || m_velocityRelaxFactor > 1)
  {
    message = "Velocity relaxation factor must be greater than 0 and less than or equal to 1";
    return false;
  }

  return true;
}

bool FVHMComponent::initializeInletOutletFlowBCArguments(QString &message)
{
  message = "";

  m_inletFlowBCArgument->clear();
  m_inletFlowBCArgument->findAssociatedCVGeometries();
  m_inletFlowBCArgument->prepare();

  m_outletFlowBCArgument->clear();
  m_outletFlowBCArgument->findAssociatedCVGeometries();
  m_outletFlowBCArgument->prepare();

  return true;
}

bool FVHMComponent::initializeSourceSinkFlowsBCArguments(QString &message)
{
  message = "";
  return true;
}

bool FVHMComponent::initializeWSEBCArguments(QString &message)
{
  message = "";

  m_outletWSEBC->setVelocityRelaxationFactor(m_outletWSEBCVRelFactor);
  m_outletWSEBC->clear();
  m_outletWSEBC->findAssociatedCVGeometries();
  m_outletWSEBC->prepare();

  m_inletWSEBC->setVelocityRelaxationFactor(m_inletWSEBCVRelFactor);
  m_inletWSEBC->clear();
  m_inletWSEBC->findAssociatedCVGeometries();
  m_inletWSEBC->prepare();

  return true;
}

bool FVHMComponent::initializeCriticalDepthOutflowBCArguments(QString &message)
{
  message = "";
  m_criticalDepthOutflowBC->clear();
  m_criticalDepthOutflowBC->findAssociatedCVGeometries();
  m_criticalDepthOutflowBC->prepare();

  return true;
}

bool FVHMComponent::initializeInitialWSEBCArguments(QString &message)
{
  message = "";

  m_initialWSEBC->clear();
  m_initialWSEBC->findAssociatedCVGeometries();
  m_initialWSEBC->prepare();

  return true;
}

bool FVHMComponent::initializeWSESlopeBCArguments(QString &message)
{
  message = "";
  //  m_outletWSESlope->setVelocityRelaxationFactor(m_outletWSEBCVRelFactor);
  m_outletWSESlope->clear();
  m_outletWSESlope->findAssociatedCVGeometries();
  m_outletWSESlope->prepare();

  return true;
}

bool FVHMComponent::initializePrecipitationArguments(QString &message)
{
  message = "";
  m_precipitationArgument->clear();
  m_precipitationArgument->findAssociatedCVGeometries();
  m_precipitationArgument->prepare();
  return true;
}

bool FVHMComponent::initializeInitialUniformBCArguments(QString &message)
{
  message = "";

  m_uniDepth=(*m_initialParamArguments)["Uniform Depth"];
  m_uniUVel = (*m_initialParamArguments)["Uniform U-Velocity"];
  m_uniVVel = (*m_initialParamArguments)["Uniform V-Velocity"];
  m_uniMannings  = (*m_initialParamArguments)["Uniform Mannings"];
  m_viscosity = (*m_initialParamArguments)["Viscosity"];
  m_inletOutletWSESlope = (*m_initialParamArguments)["Min Outlet WSE Slope"];
  m_uniWSE = (*m_initialParamArguments)["Uniform WSE"];
  m_inletWSEBCVRelFactor = (*m_initialParamArguments)["Inlet WSEBC Vel Relaxation Factor"];
  m_outletWSEBCVRelFactor = (*m_initialParamArguments)["Outlet WSEBC Vel Relaxation Factor"];

  double missingData =  dynamic_cast<HydroCouple::IQuantity*>(m_initialParamArguments->valueDefinition())->missingValue().toDouble();

  if(m_uniWSE != missingData)
  {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(int i = 0; i < m_numCells ; i++)
    {
      TriCV *cv = m_controlVolumes[i];
      cv->maxH = std::numeric_limits<double>::lowest();
      cv->maxZ = cv->maxH;
      cv->maxVel[0] = cv->maxH;
      cv->maxVel[1] = cv->maxH;
      cv->setVFRWSE(m_uniWSE,true);
      cv->setVFRWSE(m_uniWSE);
    }
  }

  if(m_uniDepth != missingData)
  {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(int i = 0; i < m_numCells ; i++)
    {
      TriCV *cv = m_controlVolumes[i];
      cv->setVFRWSE(m_uniDepth,true);
      cv->setVFRWSE(m_uniDepth);
    }
  }

  if(m_uniUVel != missingData)
  {

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(int i = 0; i < m_numCells ; i++)
    {
      TriCV *cv = m_controlVolumes[i];
      cv->prevVel[0].value = m_uniUVel;
      cv->vel[0].value = m_uniUVel;
      cv->prevIterVel[0] = m_uniUVel;
    }
  }

  if(m_uniVVel != missingData)
  {

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(size_t i = 0; i < m_controlVolumes.size() ; i++)
    {
      TriCV *cv = m_controlVolumes[i];
      cv->prevVel[1].value = m_uniVVel;
      cv->vel[1].value = m_uniVVel;
      cv->prevIterVel[1] = m_uniVVel;
    }
  }

  if(m_uniMannings != missingData)
  {

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(size_t i = 0; i < m_controlVolumes.size() ; i++)
    {
      TriCV *cv = m_controlVolumes[i];
      cv->mannings = m_uniMannings;
    }
  }

  return true;
}

bool FVHMComponent::initializeWettingAndDryingArguments(QString &message)
{
  message = "";

  m_wetCellDepth = (*m_wettingAndDryingArgument)["Wet Cell Depth"];
  m_wetCellNodeDepth = (*m_wettingAndDryingArgument)["Wet Cell Node Depth"];
  m_viscousSubLayerDepth = (*m_wettingAndDryingArgument)["Minimum Viscous Sublayer Depth"];
  m_wetDryFactor = (*m_wettingAndDryingArgument)["Wet Dry Factor"];

  if(m_wetDryFactor < 0 || m_wetDryFactor >= 1.0)
    m_wetDryFactor = 0.025;

  return true;
}

bool FVHMComponent::initializeAdvectionDiffusionSchemeArguments(QString &message)
{
  message = "";

  //Use switch later for differen turbulence options
  m_useEddyViscosity = (*m_advectionDiffusionSchemeArgument)["Use Turbulence"] ? 1.0 : 0.0;
  m_turbulenceScheme = (*m_advectionDiffusionSchemeArgument)["Turbulence Scheme"];
  m_smargorinskyCoefficient = (*m_advectionDiffusionSchemeArgument)["Smagorinsky-Lilly Coefficient"];
  m_parabolicEddyViscosityConstant = (*m_advectionDiffusionSchemeArgument)["Parabolic Eddy Viscosity Constant"];
  m_advectionScheme =  (*m_advectionDiffusionSchemeArgument)["Advection Scheme"];

  return true;
}

bool FVHMComponent::initializeMiscOptionsArguments(QString &message)
{
  message = "";

  m_verbose = (*m_miscOptionsArgument)["Verbose"] ? true : false;
  m_printFrequency = (*m_miscOptionsArgument)["Print Frequency"];
  m_useWall = (*m_miscOptionsArgument)["Use Wall Friction"];

  if(m_useWall)
    m_useWall = 1.0;

  if(m_printFrequency < 1)
  {
    m_printFrequency = 1;
  }

  m_writeFrequency = (*m_miscOptionsArgument)["Output Write Flush Frequency"];

  if(m_writeFrequency < 1)
  {
    m_writeFrequency = 1;
  }

  if(m_miscOptionsArgument->containsIdentifier("Gradient Calculation Mode"))
  {
    TriCV::gradientCalcMode = (*m_miscOptionsArgument)["Gradient Calculation Mode"];

    if(TriCV::gradientCalcMode < 0 )
    {
      TriCV::gradientCalcMode = 0;
    }
  }

  m_mpiSolverSplitThreshold = (*m_miscOptionsArgument)["MPI Solver Split Threshold"];

  if(m_mpiSolverSplitThreshold < 100)
    m_mpiSolverSplitThreshold = 100;

  return true;
}

void FVHMComponent::createInputs()
{
  createInflowInput();
}

void FVHMComponent::createInflowInput()
{

  Quantity *flowQuantity = Quantity::flowInCMS(this);

  m_inflowInput = new InflowInput("CVInflows", m_TINMeshArgument->sharedTIN(),
                                  m_timeDimension, m_patchDimension, m_edgeDimension,
                                  m_nodeDimension, flowQuantity,this);
  m_inflowInput->setCaption("Inflows (m^3/s)");

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_startDateTime- 1.0/10000000000.0, m_inflowInput);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_startDateTime, m_inflowInput);
  m_inflowInput->addTime(dt1);
  m_inflowInput->addTime(dt2);


  addInput(m_inflowInput);
}

void FVHMComponent::createOutputs()
{
  createWSEOutput();
  createDepthOutput();
  createCVAreaOutput();
}

void FVHMComponent::createWSEOutput()
{

  Quantity *wseQuantity = Quantity::lengthInMeters(this);
  wseQuantity->setCaption("Water Surface Elevation (m)");

  m_WSEOutput = new CVWSEOutput("CVNodeWSE", m_TINMeshArgument->sharedTIN(),
                                m_timeDimension, m_patchDimension, m_edgeDimension,
                                m_nodeDimension, wseQuantity,this);
  m_WSEOutput->setCaption("Water Surface Elevation (m)");

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_startDateTime - 1.0/10000000000.0, m_WSEOutput);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_startDateTime, m_WSEOutput);
  m_WSEOutput->addTime(dt1);
  m_WSEOutput->addTime(dt2);

  addOutput(m_WSEOutput);
}

void FVHMComponent::createDepthOutput()
{

  Quantity *depthQuantity = Quantity::lengthInMeters(this);
  depthQuantity->setCaption("Water Depth (m)");

  m_depthOutput = new CVDepthOutput("CVWaterDepth", m_TINMeshArgument->sharedTIN(),
                                    m_timeDimension, m_patchDimension, m_edgeDimension,
                                    m_nodeDimension, depthQuantity,this);
  m_depthOutput->setCaption("Water Depth (m)");

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_startDateTime - 1.0/10000000000.0, m_depthOutput);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_startDateTime, m_depthOutput);
  m_depthOutput->addTime(dt1);
  m_depthOutput->addTime(dt2);

  addOutput(m_depthOutput);
}

void FVHMComponent::createCVAreaOutput()
{

  Quantity *areaQuantity = Quantity::areaInSquareMeters(this);
  m_CVAreasOutput = new GeometryOutputDouble("CVAreas", IGeometry::TriangleZ, m_geometryDimension,
                                             areaQuantity,this);
  m_CVAreasOutput->setCaption("Control Volume Areas (m^2)");


  if(m_sharedTriangles.length())
  {
    m_CVAreasOutput->addGeometries(m_sharedTriangles);
  }

  for(int i = 0 ; i < m_sharedTriangles.length(); i++)
  {
    m_CVAreasOutput->setValue(i, &m_controlVolumes[i]->area);
  }

  addOutput(m_CVAreasOutput);
}

bool FVHMComponent::isInfOrNan(double number)
{
#ifdef UTAH_CHPC
  return std::isinf(number) || std::isnan(number);
#else
  return isinf(number) || isnan(number);
#endif
}

bool FVHMComponent::getBoolFromString(const QString &boolString)
{
  bool convertedToInt = false;
  int value = boolString.toInt(&convertedToInt);

  if(!boolString.trimmed().compare("True", Qt::CaseInsensitive) ||
     !boolString.trimmed().compare("Yes", Qt::CaseInsensitive) ||
     !boolString.trimmed().compare("1", Qt::CaseInsensitive) ||
     (convertedToInt && value)
     )
  {
    return true;
  }

  return false;
}

void FVHMComponent::deleteLater(HCGeometry *geometry)
{
  geometry->index();
}

