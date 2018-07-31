/*!
 *  \file    fvhmcompopnent.h
 *  \author  Caleb Amoa Buahin <caleb.buahin@gmail.com>
 *  \version 1.0.0
 *  \section Description
 *  fvhmcompopnent.h, associated files and libraries are free software;
 *  you can redistribute it and/or modify it under the terms of the
 *  Lesser GNU Lesser General Public License as published by the Free Software Foundation;
 *  either version 3 of the License, or (at your option) any later version.
 *  fvhmcompopnent.h its associated files is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.(see <http://www.gnu.org/licenses/> for details)
 *  \date 2014-2016
 *  \pre
 *  \bug
 *  \todo
 *  \warning
 */

#ifndef FVHMCOMPONENT_H
#define FVHMCOMPONENT_H

#include "fvhmcomponent_global.h"
#include "sparsemat.h"
#include "hydrocouplespatial.h"
#include "temporal/abstracttimemodelcomponent.h"

#include <QDateTime>
#include <math.h>
#include <QTextStream>
#include <QVector>
#include <tuple>

#ifdef USE_OPENMP
#include <omp.h>
#endif

class FVHMComponentInfo;
class FVHMInputObjectItem;
class FVHMOutputObjectItem;
class IdBasedArgumentString;
class IdBasedArgumentDouble;
class IdBasedArgumentDateTime;
class TINArgumentDouble;
class HCTriangle;
class IArgument;
class InletFlowBCArgument;
class OutletFlowBCArgument;
class WSEBC;
class InitialWSEBC;
class OutletWSESlope;
class TimeSeriesArgumentDouble;
class TimeGeometryArgumentDouble;
struct TriCV;
class Octree;
class HCGeometry;
class Edge;
class HCVertex;
class PrecipBC;
class InflowInput;
class CVWSEOutput;
class CVDepthOutput;
class GeometryOutputDouble;
class Dimension;
struct Vect;
class CriticalDepthOutflowBC;
class EdgeFluxesIO;


namespace netCDF
{
  class NcFile;
}

class FVHMComponent;

typedef int (FVHMComponent::*PerformStepFunction)(double tStep, int &numIterations,
                                                  double &uvelRelResidualNormInit, double &uvelRelResidualNormFin, int &minUVelSolvIters, int &maxUVelSolvIters,
                                                  double &vvelRelResidualNormInit, double &vvelRelResidualNormFin, int &minVVelSolvIters, int &maxVVelSolvIters,
                                                  double &pressureRelRisdualNormInit, double &pressureRelRisdualNormFin, int &minPressSolvIters, int &maxPressSolvIters,
                                                  double &continuityRelResidualNormInit, double &continuityRelResidualNormFin,
                                                  bool &converged, QString &errorMessage);

class FVHMCOMPONENT_EXPORT FVHMComponent : public AbstractTimeModelComponent
{
    friend class FVHMComponentInfo;
    friend struct TriCV;
    friend class WSEBC;
    friend class InletFlowBCArgument;
    friend class OutletWSESlope;
    friend class PrecipBC;
    friend class CriticalDepthOutflowBC;
    friend class EdgeFluxesIO;
    friend class InletOutletFlowBCArgument;
    friend class InitialWSEBC;
    friend class OutletFlowBCArgument;

    Q_OBJECT

  public:

    enum AdaptiveTSMode
    {
      MaxCourantNumber,
      RMSCourantNumber
    };

    enum ErrorCode
    {
      NoError,
      SolverFailedToConverge,
      FailedToConverge,
      CriticalFailure,
    };

    FVHMComponent(const QString &id, FVHMComponentInfo* componentInfo = nullptr);

    FVHMComponent(const QString &id, const QString &caption, FVHMComponentInfo* componentInfo = nullptr);

    virtual ~FVHMComponent();

    QList<QString> validate() override;

    void prepare() override;

    void update(const QList<HydroCouple::IOutput*> &requiredOutputs = QList<HydroCouple::IOutput*>()) override;

    void finish() override;

    TINArgumentDouble *meshArgument() const;

    const std::vector<TriCV*> &controlVolumes() const;

    HCVertex *findVertex(HydroCouple::Spatial::IPoint *startPoint);

    Edge *findNextEdge(HCVertex* origin, HydroCouple::Spatial::IPoint* destPoint);

    TriCV *findCoincidentCV(HydroCouple::Spatial::IPoint *point);

    Octree *triangleOctree() const;

    double defaultInletOutletWSESlope() const;

    double viscousSubLayerDepth() const;

    double startDateTime() const;

    double endDateTime() const;

  protected:

    void initializeFailureCleanUp() override;

  private:

    void createDimensions();

    void createArguments() override;

    void createMeshArguments();

    void createTimeArguments();

    void createInputFilesArguments();

    void createOutputFilesArguments();

    void createEdgeFluxesIOArguments();

    void createConvergenceArguments();

    void createSolverOptionsArguments();

    void createPressureVelCouplingArguments();

    void createInitialUniformBCArguments();

    void createWettingAndDryingArguments();

    void createInletOutletFlowBCArguments();

    void createSourceSinkFlowsBCArguments();

    void createWSEBCArguments();

    void createCriticalDepthOutflowBCArguments();

    void createInitialWSEBCArguments();

    void createWSESlopeArguments();

    void createPrecipitationArguments();

    void createAdvectionDiffusionSchemeArguments();

    void createMiscOptionsArguments();

    bool initializeArguments(QString &message) override;

    bool initializeMeshArguments(QString &message);

    bool initializeTimeArguments(QString &message);

    bool initializeOutputFilesArguments(QString &message);

    bool initializeNetCDFOutputFile(QString &message);

    bool initializeShapeOutputFile(QString &message);

    bool initializeLogFile(QString &message);

    bool initializeEdgeFluxesIOArguments(QString &message);

    bool initializeConvergenceArguments(QString &message);

    bool initializeSolverOptionsArguments(QString &message);

    bool initializePressureVelCouplingArguments(QString &message);

    bool initializeInletOutletFlowBCArguments(QString &message);

    bool initializeSourceSinkFlowsBCArguments(QString &message);

    bool initializeWSEBCArguments(QString &message);

    bool initializeCriticalDepthOutflowBCArguments(QString &message);

    bool initializeInitialWSEBCArguments(QString &message);

    bool initializeWSESlopeBCArguments(QString &message);

    bool initializePrecipitationArguments(QString &message);

    bool initializeInitialUniformBCArguments(QString &message);

    bool initializeWettingAndDryingArguments(QString &message);

    bool initializeAdvectionDiffusionSchemeArguments(QString &message);

    bool initializeMiscOptionsArguments(QString &message);

    void createInputs() override;

    void createInflowInput();

    void createOutputs() override;

    void createWSEOutput();

    void createDepthOutput();

    void createCVAreaOutput();

    double getNextTimeStep();

    void applyInitialBoundaryConditions();

    void readRestartFile();

    void applyBoundaryConditions(double time);

    double getMinOutputTime(const QList<HydroCouple::IOutput*> &requiredOutputs);

    void getMinIOutputTime(const HydroCouple::IOutput *output, double &minTime);

    void prepareForNextTimeStep(double &minU, double &maxU, double &minV, double &maxV, double &minH, double &maxH, double &minZ, double &maxZ);

    void setWetAndContCells();


    int performSimpleTimeStep(double tStep, int &numIterations,
                              double &uvelRelResidualNormInit, double &uvelRelResidualNormFin, int &minUVelSolvIters, int &maxUVelSolvIters,
                              double &vvelRelResidualNormInit, double &vvelRelResidualNormFin, int &minVVelSolvIters, int &maxVVelSolvIters,
                              double &pressureRelRisdualNormInit, double &pressureRelRisdualNormFin, int &minPressSolvIters, int &maxPressSolvIters,
                              double &continuityRelResidualNormInit, double &continuityRelResidualNormFin,
                              bool &converged, QString &errorMessage);

    int performSimpleCTimeStep(double tStep, int &numIterations,
                               double &uvelRelResidualNormInit, double &uvelRelResidualNormFin, int &minUVelSolvIters, int &maxUVelSolvIters,
                               double &vvelRelResidualNormInit, double &vvelRelResidualNormFin, int &minVVelSolvIters, int &maxVVelSolvIters,
                               double &pressureRelRisdualNormInit, double &pressureRelRisdualNormFin, int &minPressSolvIters, int &maxPressSolvIters,
                               double &continuityRelResidualNormInit, double &continuityRelResidualNormFin,
                               bool &converged, QString &errorMessage);

    int performPISOTimeStep(double tStep, int &numIterations,
                            double &uvelRelResidualNormInit, double &uvelRelResidualNormFin, int &minUVelSolvIters, int &maxUVelSolvIters,
                            double &vvelRelResidualNormInit, double &vvelRelResidualNormFin, int &minVVelSolvIters, int &maxVVelSolvIters,
                            double &pressureRelRisdualNormInit, double &pressureRelRisdualNormFin, int &minPressSolvIters, int &maxPressSolvIters,
                            double &continuityRelResidualNormInit, double &continuityRelResidualNormFin,
                            bool &converged, QString &errorMessage);


    bool neighbourIsWet(TriCV *cv, double wetCellDepth);

    ErrorCode solveMomentumEquations(double tStep, int uv, double &uvelRelResidualNorm, int &numSolvIters, QString &errorMessage);

    void interpolateWettedCellVelocities();

    void interpolateWettedCellVelocity(TriCV *cv);

    double calculateFluxLimiter(const TriCV* cv, int faceIndex, double idwfactor,
                                const Vect &grad_c, const Vect &grad_d,
                                double phi_c, double phi_d, double r_cdX, double r_cdY);

    ErrorCode solvePressureCorrection(double tStep, double &pressureRelResidualNorm, double minorTermsCoeff, int &numSolvIters, QString &errorMessage);

    void applyPressureCorrections();

    void applyVelocityCorrections();

    void calculateContinuityResiduals(double tStep, double &contRelResidualNorm);

    void updateExternalInflowOutflowTotals(double tStep);

    ErrorCode solveConstituentEquations(int cIndex, double tStep, double &vvelRelResidualNorm, QString &errorMessage);

    void calculateFriction(int uv = 0);

    void calculateFriction(TriCV *cv, int uv);

    void calculateFaceVelocities();

    void calculateBoundaryFaceVelocities();

    void calculateBoundaryFaceVelocities(TriCV *cv);

    void interpolateFaceVelocityFromGradient(TriCV *cv, int faceIndex);

    void calculateCellVelocityGradients();

    void calculateCellEddyViscosities();

    void calculateCellWSEGradients();

    void calculateCellZCorrGradients();

    ErrorCode mpiSolve(const SparseMatrix &A, const double b[], double x[],
                       double residuals[], double &relativeResidualNorm, int &numIterations, QString &errorMessage);

    int solve(MPI_Comm communicator, const SparseMatrix &A, int ilower, int iupper,
              const double b[], double x[], double residuals[], double &relativeResidualNorm, int &numIterations);

    static double sign(double value);

    static double l2Norm(const double residuals[], int length);

    static double relativeResidualNorm(const double residuals[], const double values[], int length);

    bool isWet(const TriCV* cv);

    static bool isInfOrNan(double number);

    static bool getBoolFromString(const QString &boolString);

    static std::tuple<double, double, double> calculateZeroGradientFaceVelocity(TriCV* cv, int face);

    void writeOutputs();

    void writeToNetCDF();

    void writeToCSV();

    void writeToLogFile(int iters,
                        double uInitRes, double uFinRes,
                        double vInitRes, double vFinRes,
                        double contInitRes, double contFinRes,
                        double pressInitRes, double pressFinRes);

    void deleteControlVolumes();

    void closeOutputCSVFile();

    void closeLogFile();

    void resetOctree();

    static void deleteLater(HCGeometry *geometry);

  private:
    IdBasedArgumentString *m_inputFiles = nullptr;
    IdBasedArgumentString *m_outputFilesArgument = nullptr;
    TINArgumentDouble *m_TINMeshArgument = nullptr;
    IdBasedArgumentDouble *m_transportConstituentArguments = nullptr;
    IdBasedArgumentDouble *m_timeStepArguments = nullptr;
    IdBasedArgumentDouble *m_convergenceArguments = nullptr;
    IdBasedArgumentDouble *m_solverOptionsArguments = nullptr;
    IdBasedArgumentDouble *m_pressureVelCouplingArguments = nullptr;
    IdBasedArgumentDouble *m_initialParamArguments = nullptr;
    IdBasedArgumentDouble *m_octreeArguments = nullptr;
    IdBasedArgumentDateTime *m_simulationDurationArgument = nullptr;
    IdBasedArgumentDouble *m_wettingAndDryingArgument = nullptr;
    IdBasedArgumentDouble *m_advectionDiffusionSchemeArgument = nullptr;
    IdBasedArgumentDouble *m_miscOptionsArgument = nullptr;
    InletFlowBCArgument *m_inletFlowBCArgument = nullptr;
    OutletFlowBCArgument *m_outletFlowBCArgument = nullptr;

    WSEBC *m_outletWSEBC = nullptr;
    WSEBC *m_inletWSEBC = nullptr;
    CriticalDepthOutflowBC *m_criticalDepthOutflowBC = nullptr;
    InitialWSEBC *m_initialWSEBC = nullptr;
    EdgeFluxesIO *m_edgeFluxesIO = nullptr;


    OutletWSESlope *m_outletWSESlope = nullptr;
    TimeGeometryArgumentDouble *m_sourceSinkFlowBCArgument= nullptr;
    PrecipBC *m_precipitationArgument = nullptr;

    CVWSEOutput *m_WSEOutput = nullptr;
    CVDepthOutput *m_depthOutput = nullptr;

    InflowInput *m_inflowInput = nullptr;
    GeometryOutputDouble *m_CVAreasOutput = nullptr;

    Dimension *m_patchDimension = nullptr,
    *m_edgeDimension = nullptr,
    *m_nodeDimension = nullptr,
    *m_idDimension = nullptr,
    *m_geometryDimension = nullptr,
    *m_timeDimension = nullptr;

    bool m_useAdaptiveTimeStep = true,
    m_verbose = true,
    m_writeShapefile = false,
    m_writeActiveCellsOnly = true;

    double m_startDateTime = 0,
    m_endDateTime = 1.0,
    m_currentDateTime = 0.0,
    m_previouslyWrittenDateTime = std::numeric_limits<double>::lowest(),
    m_minTimeStep = 0.1, m_maxTimeStep = 10.0, m_maxTimeStepIncreaseCF = 0.1, m_maxTimeStepDecreaseCF = 0.1,
    m_outputTimeStep = 30.0,
    m_maxCourantNumber = 1.0,
    m_nextOutputTime = 0.0,
    m_RMSCourantNumber = 1.0,
    m_timeStepRelaxFactor = 0.9,
    m_wetCellDepth = 10e-3,
    m_wetCellNodeDepth = 10e-7,
    m_viscousSubLayerDepth = 10e-10,
    m_viscosity = 1.004 *10e-6,
    m_uvelConvergenceTol = 1e-7,
    m_vvelConvergenceTol = 1e-7,
    m_pressureConvergenceTol = 1e-5,
    m_contConvergenceTol = 1e-5,
    m_pressureRelaxFactor = 0.7,
    m_velocityRelaxFactor = 0.7,
    m_solverConvergenceTol =1e-10,
    m_timeStep = 0.01,
    m_uniDepth = -99999999,
    m_uniWSE =-99999999,
    m_uniUVel = -99999999,
    m_uniVVel = -99999999,
    m_uniMannings = 0.03,
    m_useEddyViscosity = 1,
    m_inletOutletWSESlope = -0.0001,
    m_inletWSEBCVRelFactor = 1.00,
    m_outletWSEBCVRelFactor = 1.00,
    m_turbulenceScheme = 1.0,
    m_smargorinskyCoefficient = 0.2,
    m_parabolicEddyViscosityConstant = 0.7,
    m_advectionDampeningFactor = 0.3,
    m_timeStepCount = 0.0,
    m_convergedCount = 0.0,
    m_useAdvection = 1.0,
    m_useWall = 1.0,
    m_epsilon = 0.0,
    m_maxWSEGradient = 3.0,
    m_numFractionalSteps = 3.0,
    m_wetDryFactor = 0.05;

    QDateTime m_qtDateTime;

    const double g = 9.80665;

    int m_itersPerTimeStep  = 150,
    m_innerItersPerTimeStep = 2,
    m_solverMaximumNumberOfIterations = 150,
    m_solverAMGPreconditionerNumberOfIterations = 1,
    m_solverType = 0,
    m_initTimeStepCycle = 0,
    m_numFixedTimeSteps = 100,
    m_pressureVelocityCouplingType = 0,
    m_numCells = 0,
    m_advectionScheme = 0,
    m_printFrequency = 5,
    m_printFrequencyCounter = 0,
    m_writeFrequency = 20,
    m_writeFrequencyCounter = 0,
    m_numPressureCorrectionIterations = 5,
    m_currentCoeff = 0,
    m_mpiSolverSplitThreshold = 100,
    m_numWetCells, m_numContCells;

    AdaptiveTSMode m_adaptiveTSMode = AdaptiveTSMode::MaxCourantNumber;

    std::vector<TriCV*> m_controlVolumes;
    std::vector<int> m_wetCells, m_contCells;
    QFileInfo m_outputNetCDFFile, m_outputShapefile, m_outputShapefileJoin, m_logFile, m_edgeFluxesFile;
    QDir m_outputDir;
    QList<QSharedPointer<HCGeometry>> m_sharedTriangles;
    QFile m_logFileIO, m_CSVOutputIO;
    netCDF::NcFile *m_outputNetCDF = nullptr;
    QTextStream m_logFileTextStream, m_CSVOutputTextStream;
    Octree *m_octree;
    PerformStepFunction performStep = nullptr;

};

Q_DECLARE_METATYPE(FVHMComponent*)

#endif // FVHMCOMPONENT_H
