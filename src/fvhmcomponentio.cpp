#include "stdafx.h"
#include "fvhmcomponent.h"
#include "controlvolume.h"
#include "temporal/timedata.h"
#include "core/idbasedargument.h"

#include "spatial/polygon.h"
#include "spatial/geometryfactory.h"
#include "spatial/tinargument.h"
#include "cvwseoutput.h"
#include "inflowinput.h"
#include "edgefluxesio.h"

#include "netcdf"

using namespace netCDF;
using namespace netCDF::exceptions;

bool FVHMComponent::initializeOutputFilesArguments(QString &message)
{
  if(mpiProcessRank() == 0)
  {
    if(initializeNetCDFOutputFile(message) &&
       initializeShapeOutputFile(message) &&
       initializeLogFile(message))
    {
      return true;
    }
  }
  else
  {
    if(initializeNetCDFOutputFile(message))
    {
      return true;
    }
  }

  return false;
}

bool FVHMComponent::initializeNetCDFOutputFile(QString &message)
{

  QString file = (*m_outputFilesArgument)["Output NetCDF File"];

  if(file.isEmpty())
  {
    return true;
  }
  else if(!file.isEmpty() && !file.isNull() && !(m_outputNetCDFFile = getAbsoluteFilePath(file)).absoluteDir().exists())
  {
    message = "NetCDF output file directory does not exist: " + file;
    return false;
  }

  if(mpiProcessRank() == 0)
  {
    if(m_outputNetCDF)
    {
      m_outputNetCDF->sync();
      m_outputNetCDF->close();
      delete m_outputNetCDF;
      m_outputNetCDF = nullptr;
    }



    GeometryFactory::writeTINToNetCDF(m_TINMeshArgument->TINInternal(), m_outputNetCDFFile.absoluteFilePath(), message);

    try
    {

      m_outputNetCDF = new NcFile(m_outputNetCDFFile.absoluteFilePath().toStdString() , NcFile::write);

      NcDim edge = m_outputNetCDF->getDim("edges");

      if(edge.isNull())
        m_outputNetCDF->addDim("edges",3);

      //time variable
      NcDim timeDim = m_outputNetCDF->addDim("time");
      // NcDim vertDim = m_outputNetCDF->getDim("vertices");

      NcVar timeVar = m_outputNetCDF->addVar("time", NcType::nc_DOUBLE , timeDim);
      timeVar.putAtt("time:long_name","time");
      timeVar.putAtt("time:units","days since 1858-11-17 0:0:0");
      timeVar.putAtt("time:calendar","modified_julian");

      //add variables

      NcVar inflowTotalVar = m_outputNetCDF->addVar("externalInflowTotalVolume","double", std::vector<std::string>({"time"}));
      inflowTotalVar.putAtt("externalInflowTotalVolume:long_name","inflow");
      inflowTotalVar.putAtt("externalInflowTotalVolume:units","m^3");

      NcVar outflowTotalVar = m_outputNetCDF->addVar("externalOutflowTotalVolume","double", std::vector<std::string>({"time"}));
      outflowTotalVar.putAtt("externalOutflowTotalVolume:long_name","inflow");
      outflowTotalVar.putAtt("externalOutflowTotalVolume:units","m^3");

      NcVar inflowVar = m_outputNetCDF->addVar("externalInflowTotalCellVolume","double", std::vector<std::string>({"time","triangles"}));
      inflowVar.putAtt("externalInflowTotalCellVolume:long_name","inflow");
      inflowVar.putAtt("externalInflowTotalCellVolume:units","m^3");

      NcVar outflowVar = m_outputNetCDF->addVar("externalOutflowTotalCellVolume","double", std::vector<std::string>({"time","triangles"}));
      outflowVar.putAtt("externalOutflowTotalCellVolume:long_name","inflow");
      outflowVar.putAtt("externalOutflowTotalCellVolume:units","m^3");

      NcVar depthVar = m_outputNetCDF->addVar("depth","double", std::vector<std::string>({"time","triangles"}));
      depthVar.putAtt("depth:long_name","depth");
      depthVar.putAtt("depth:units","m");

      NcVar maxdepthVar = m_outputNetCDF->addVar("maxdepth","double", std::vector<std::string>({"triangles"}));
      maxdepthVar.putAtt("maxdepth:long_name","depth");
      maxdepthVar.putAtt("maxdepth:units","m");

      NcVar elevVar = m_outputNetCDF->addVar("elev","double", std::vector<std::string>({"time","triangles"}));
      elevVar.putAtt("elev:long_name","elevation");
      elevVar.putAtt("elev:units","m");

      NcVar maxElevVar = m_outputNetCDF->addVar("maxelev","double", std::vector<std::string>({"triangles"}));
      maxElevVar.putAtt("maxelev:long_name","elevation");
      maxElevVar.putAtt("maxelev:units","m");

      NcVar uVar = m_outputNetCDF->addVar("uvel","double", std::vector<std::string>({"time","triangles"}));
      uVar.putAtt("uvel:long_name","u_velocity");
      uVar.putAtt("uvel:units","m/s");

      NcVar maxUVar = m_outputNetCDF->addVar("maxuvel","double", std::vector<std::string>({"triangles"}));
      maxUVar.putAtt("maxuvel:long_name","u_velocity");
      maxUVar.putAtt("maxuvel:units","m/s");

      NcVar vVar = m_outputNetCDF->addVar("vvel","double", std::vector<std::string>({"time","triangles"}));
      vVar.putAtt("vvel:long_name","v_velocity");
      vVar.putAtt("vvel:units","m/s");

      NcVar maxVVar = m_outputNetCDF->addVar("maxvvel","double", std::vector<std::string>({"triangles"}));
      maxVVar.putAtt("maxvvel:long_name","v_velocity");
      maxVVar.putAtt("maxvvel:units","m/s");

      NcVar uVarRes = m_outputNetCDF->addVar("uvelresidual","double", std::vector<std::string>({"time","triangles"}));
      uVarRes.putAtt("uvelresidual:long_name","uvel_residual_velocity");
      uVarRes.putAtt("uvelresidual:units","m/s");

      NcVar vVarRes = m_outputNetCDF->addVar("vvelresidual","double", std::vector<std::string>({"time","triangles"}));
      vVarRes.putAtt("vvelresidual:long_name","vvel_residual_velocity");
      vVarRes.putAtt("vvelresidual:units","m/s");

      NcVar continuityRes = m_outputNetCDF->addVar("continuityresidual","double", std::vector<std::string>({"time","triangles"}));
      continuityRes.putAtt("continuityresidual:long_name","continuity_residual");
      continuityRes.putAtt("continuityresidual:units","m3/s");

      NcVar continuityResTot = m_outputNetCDF->addVar("totalcontinuityresidual","double", std::vector<std::string>({"time"}));
      continuityResTot.putAtt("totalcontinuityresidual:long_name","totalcontinuity_residual");
      continuityResTot.putAtt("totalcontinuityresidual:units","m3/s");

      NcVar pressurRes = m_outputNetCDF->addVar("pressureresidual","double", std::vector<std::string>({"time","triangles"}));
      pressurRes.putAtt("pressureresidual:long_name","pressure_residual");
      pressurRes.putAtt("pressureresidual:units","m");


      NcVar edgeFlow = m_outputNetCDF->addVar("edgeflow","double",std::vector<std::string>({"time","triangles","edges"}));
      edgeFlow.putAtt("edgeflow:long_name","control_volume_edge_flow");
      edgeFlow.putAtt("edgeflow:units","m^3/s");

      //add other variables

      m_outputDir = m_outputNetCDFFile.absoluteDir();
    }
    catch(NcException &e)
    {
      message = QString(e.what());
      return false;
    }
  }
  else
  {
    m_outputDir = m_outputNetCDFFile.absoluteDir();
  }

  return true;
}

bool FVHMComponent::initializeShapeOutputFile(QString &message)
{
  QString file = (*m_outputFilesArgument)["Output Shapefile"];

  if(!file.isEmpty() && !file.isNull() && !(m_outputShapefile = getAbsoluteFilePath(file)).absoluteDir().exists())
  {
    message = "Output shapefile directory does not exist: " + file;
    m_writeShapefile = false;
    return false;
  }

  if(!file.isEmpty() && !file.isNull())
  {
    QString activeOnly = (*m_outputFilesArgument)["WriteActiveCellsOnly"];
    m_writeActiveCellsOnly = getBoolFromString(activeOnly);

    closeOutputCSVFile();

    m_outputShapefileJoin = QFileInfo(m_outputShapefile.absoluteFilePath().replace(m_outputShapefile.suffix(),"csv"));
    m_CSVOutputIO.setFileName(m_outputShapefileJoin.absoluteFilePath());

    if(GeometryFactory::writeTINPolygons(m_TINMeshArgument->TINInternal(), m_outputShapefile.absoluteFilePath(),"ESRI Shapefile"))
    {
      if(m_CSVOutputIO.open(QIODevice::WriteOnly | QIODevice::Truncate))
      {
        QFile vrtFile(m_outputShapefile.absoluteFilePath().replace(m_outputShapefile.suffix(),"vrt"));
        vrtFile.open(QIODevice::WriteOnly);

        QTextStream vrtStream(&vrtFile);

        vrtStream << "<OGRVRTDataSource>\n"
                  << "<OGRVRTLayer name=\""
                  << m_outputShapefile.baseName()
                  << "\">\n"
                  << "<SrcDataSource>"
                  << m_outputShapefile.fileName().replace(m_outputShapefile.suffix(),"csv")
                  << "</SrcDataSource>\n"
                  << "<GeometryType>wkbPolygon25D</GeometryType>\n"
                  << "<GeometryField encoding=\"WKT\" field=\"PolygonWKT\"/>\n"
                  << "<Field name=\"TriIndex\" src=\"TriIndex\" type=\"Integer\"/>\n"
                  << "<Field name=\"Cx\" src=\"Cx\" type=\"Real\"/>\n"
                  << "<Field name=\"Cy\" src=\"Cy\" type=\"Real\"/>\n"
                  << "<Field name=\"Cz\" src=\"Cz\" type=\"Real\"/>\n"
                  << "<Field name=\"DateTime\" src=\"DateTime\" type=\"String\"/>\n"
                  << "<Field name=\"DTIndex\" src=\"DTIndex\" type=\"Integer\"/>\n"
                  << "<Field name=\"Depth\" src=\"Depth\" type=\"Real\"/>\n"
                  << "<Field name=\"Elevation\" src=\"Elevation\" type=\"Real\"/>\n"
                  << "<Field name=\"UVel\" src=\"UVel\" type=\"Real\"/>\n"
                  << "<Field name=\"VVel\" src=\"VVel\" type=\"Real\"/>\n"
                  << "<Field name=\"EdgeVel1\" src=\"EdgeVel1\" type=\"Real\"/>\n"
                  << "<Field name=\"EdgeVel2\" src=\"EdgeVel2\" type=\"Real\"/>\n"
                  << "<Field name=\"EdgeVel3\" src=\"EdgeVel3\" type=\"Real\"/>\n"
                  << "<Field name=\"Inflow\" src=\"Inflow\" type=\"Real\"/>\n"
                  << "<Field name=\"URes\" src=\"URes\" type=\"Real\"/>\n"
                  << "<Field name=\"VRes\" src=\"VRes\" type=\"Real\"/>\n"
                  << "<Field name=\"ContRes\" src=\"ContRes\" type=\"Real\"/>\n"
                  << "</OGRVRTLayer>\n"
                  << "</OGRVRTDataSource>";

        vrtFile.close();

        m_CSVOutputTextStream.setDevice(&m_CSVOutputIO);
        m_CSVOutputTextStream.setRealNumberPrecision(10);
        m_CSVOutputTextStream.setRealNumberNotation(QTextStream::SmartNotation);

        m_CSVOutputTextStream << "PolygonWKT," << "TriIndex," << "Cx," << "Cy,"<<"Cz," << "DateTime," << "DTIndex,"
                              << "Depth," << "Elevation," << "UVel," << "VVel,"
                              << "EdgeVel1," << "EdgeVel2," << "EdgeVel3," << "Inflow,"
                              << "URes," << "VRes," << "ContRes" << endl;

        m_CSVOutputIO.flush();
      }

      m_writeShapefile = true;
      return true;
    }
  }
  else
  {
    m_writeShapefile = false;
    return true;
  }

  return false;
}

bool FVHMComponent::initializeLogFile(QString &message)
{

  QString file = (*m_outputFilesArgument)["Log File"];

  if(!file.isEmpty() && !file.isNull() && !(m_logFile = getAbsoluteFilePath(file)).absoluteDir().exists())
  {
    message = "Output log file directory does not exist: " + file;;
    return false;
  }

  if(!file.isEmpty() && !file.isNull())
  {
    closeLogFile();

    m_logFileIO.setFileName(m_logFile.absoluteFilePath());

    if(m_logFileIO.open(QIODevice::WriteOnly | QIODevice::Truncate))
    {
      m_logFileTextStream.setDevice(&m_logFileIO);
      m_logFileTextStream.setRealNumberPrecision(10);
      m_logFileTextStream.setRealNumberNotation(QTextStream::SmartNotation);
      m_logFileTextStream << "Time,"<< "TimeStep," << "NumIters,"
                          << "U-Res-Initial," << "U-Res-Final,"
                          << "V-Res-Initial," << "V-Res-Final,"
                          << "Cont-Res-Initial," << "Cont-Res-Final,"
                          <<  "P-Res-Initial," <<  "P-Res-Final" << endl;

      return true;
    }
  }
  else
  {
    return true;
  }

  return false;

}

bool FVHMComponent::initializeEdgeFluxesIOArguments(QString &message)
{
  QString file = (*m_outputFilesArgument)["Output Edge Fluxes File"];

  if(!file.isEmpty() && !file.isNull() && !(m_edgeFluxesFile = getAbsoluteFilePath(file)).absoluteDir().exists())
  {
    message = "Output log file directory does not exist: " + file;;
    return false;
  }

  m_edgeFluxesIO->clear();

  if(!file.isEmpty() && !file.isNull())
  {
    m_edgeFluxesIO->setOutputFile(m_edgeFluxesFile);
    m_edgeFluxesIO->findAssociatedCVGeometries();
    m_edgeFluxesIO->prepare();

    return true;
  }
  else
  {
    return true;
  }

  return false;
}

void FVHMComponent::readRestartFile()
{
  QString file = (*m_inputFiles)["Restart File"];
  QFileInfo restartFile;

  if(!file.isEmpty() && !file.isNull() && (restartFile = getAbsoluteFilePath(file)).absoluteDir().exists())
  {
    try
    {
      NcFile restartNetCDF(restartFile.absoluteFilePath().toStdString() , NcFile::read);

      NcDim timeDim = restartNetCDF.getDim("time");
      NcVar timeVar = restartNetCDF.getVar("time");

      if(!timeDim.isNull() && !timeVar.isNull())
      {
        size_t timeSize = timeDim.getSize();

        double *dateTimes = new double[timeSize];
        timeVar.getVar(dateTimes);

        for(size_t i = 0 ; i < timeSize ; i++)
        {
          double dateTime = dateTimes[i];

          if(m_currentDateTime == dateTime)
          {
            NcVar uvelsVar = restartNetCDF.getVar("uvel");
            NcVar vvelsVar = restartNetCDF.getVar("vvel");
            NcVar elevsVar = restartNetCDF.getVar("elev");

            double *uvels = new double[m_numCells];
            double *vvels = new double[m_numCells];
            double *elevs = new double[m_numCells];

            uvelsVar.getVar(std::vector<size_t>({i,0}), std::vector<size_t>({1,(size_t)m_numCells}), uvels);
            vvelsVar.getVar(std::vector<size_t>({i,0}), std::vector<size_t>({1,(size_t)m_numCells}), vvels);
            elevsVar.getVar(std::vector<size_t>({i,0}), std::vector<size_t>({1,(size_t)m_numCells}), elevs);


            for(size_t j = 0; j < m_controlVolumes.size() ; j++)
            {
              TriCV *cv = m_controlVolumes[j];
              cv->vel[0].value = uvels[j];
              cv->vel[1].value = vvels[j];
              cv->setVFRWSE(elevs[i]);
            }

            delete[] uvels;
            delete[] vvels;
            delete[] elevs;

          }
          else if(m_currentDateTime > dateTime && i < timeSize - 1)
          {
            double nextTime = dateTimes[i+1];
            double factor =  (m_currentDateTime - dateTime) / (nextTime - dateTime);


            NcVar uvelsVar = restartNetCDF.getVar("uvel");
            NcVar vvelsVar = restartNetCDF.getVar("vvel");
            NcVar elevsVar = restartNetCDF.getVar("elev");

            double *uvels = new double[m_numCells * 2];
            double *vvels = new double[m_numCells * 2];
            double *elevs = new double[m_numCells * 2];

            uvelsVar.getVar(std::vector<size_t>({i,0}), std::vector<size_t>({2,(size_t) m_numCells}), uvels);
            vvelsVar.getVar(std::vector<size_t>({i,0}), std::vector<size_t>({2,(size_t) m_numCells}), vvels);
            elevsVar.getVar(std::vector<size_t>({i,0}), std::vector<size_t>({2,(size_t) m_numCells}), elevs);


            for(size_t j = 0; j < m_controlVolumes.size() ; j++)
            {
              TriCV *cv = m_controlVolumes[j];

              cv->vel[0].value = uvels[j] + factor * (uvels[m_controlVolumes.size() + j] - uvels[j]);
              cv->vel[1].value = uvels[j] + factor * (vvels[m_controlVolumes.size() + j] - vvels[j]);

              double elev = elevs[j] + factor * (elevs[m_controlVolumes.size() + j] - elevs[j]);

              cv->setVFRWSE(elev);
            }

            delete[] uvels;
            delete[] vvels;
            delete[] elevs;
          }
        }

        delete[] dateTimes;
      }

    }
    catch(NcException &e)
    {
      printf("%s\n",e.what());
      return;
    }
  }
}

void FVHMComponent::writeOutputs()
{
  writeToNetCDF();
  writeToCSV();
  m_edgeFluxesIO->applyBoundaryConditions(m_currentDateTime,m_timeStep);
}

void FVHMComponent::writeToNetCDF()
{
  try
  {




    if(m_outputNetCDF && !m_outputNetCDF->isNull())
    {
      //time variable
      NcDim timeDim = m_outputNetCDF->getDim("time");
      size_t currentTime = timeDim.getSize();

      NcVar timeVar = m_outputNetCDF->getVar("time");
      timeVar.putVar(std::vector<size_t>({currentTime}), m_nextOutputTime);


      double *externalInflowValues = new double[m_numCells];
      double *externalOutflowValues = new double[m_numCells];
      double *depthValues = new double[m_numCells];
      double *maxDepthValues = new double[m_numCells];
      double *elevValues = new double[m_numCells];
      double *maxElevValues = new double[m_numCells];
      double *uVelValues = new double[m_numCells];
      double *maxUVelValues = new double[m_numCells];
      double *uVelResValues = new double[m_numCells];
      double *vVelValues = new double[m_numCells];
      double *maxVVelValues = new double[m_numCells];
      double *vVelResValues = new double[m_numCells];
      double *contResValues = new double[m_numCells];
      double *pressResValues = new double[m_numCells];
      double *edgeFlowValues = new double[m_numCells * 3];

      double contResTot = 0, inflowTotal = 0.0, outflowTotal = 0.0;

      //Depth
      NcVar inflowTotalVar = m_outputNetCDF->getVar("externalInflowTotalVolume");
      NcVar outflowTotalVar = m_outputNetCDF->getVar("externalOutflowTotalVolume");
      NcVar externalInflowVar = m_outputNetCDF->getVar("externalInflowTotalCellVolume");
      NcVar externalOutflowVar = m_outputNetCDF->getVar("externalOutflowTotalCellVolume");
      NcVar depthVar = m_outputNetCDF->getVar("depth");
      NcVar maxdepthVar = m_outputNetCDF->getVar("maxdepth");
      NcVar elevVar = m_outputNetCDF->getVar("elev");
      NcVar maxelevVar = m_outputNetCDF->getVar("maxelev");
      NcVar uVar = m_outputNetCDF->getVar("uvel");
      NcVar maxUVar = m_outputNetCDF->getVar("maxuvel");
      NcVar uResVar = m_outputNetCDF->getVar("uvelresidual");
      NcVar vVar = m_outputNetCDF->getVar("vvel");
      NcVar maxVVar = m_outputNetCDF->getVar("maxvvel");
      NcVar vResVar = m_outputNetCDF->getVar("vvelresidual");
      NcVar continuityResVar = m_outputNetCDF->getVar("continuityresidual");
      NcVar continuityResTotalVar = m_outputNetCDF->getVar("totalcontinuityresidual");
      NcVar pressureResVar = m_outputNetCDF->getVar("pressureresidual");
      NcVar edgeFlowVar = m_outputNetCDF->getVar("edgeflow");

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
      for(int i = 0; i < m_numCells ; i++)
      {
        TriCV *cv = m_controlVolumes[i];

        depthValues[i] = cv->h->value;
        maxDepthValues[i] = cv->maxH;

        elevValues[i] = cv->z->value;
        maxElevValues[i] = cv->maxZ;

        uVelValues[i] = cv->vel[0].value;
        maxUVelValues[i] = cv->maxVel[0];
        uVelResValues[i] = cv->velResidual[0];

        vVelValues[i] = cv->vel[1].value;
        maxVVelValues[i] = cv->maxVel[1];
        vVelResValues[i] = cv->velResidual[1];

        contResValues[i] = cv->contResidualTotal;
        pressResValues[i] = cv->zCorrection;

        externalInflowValues[i] = cv->totalExternalInflow;
        externalOutflowValues[i] = cv->totalExternalOutflow;

#ifdef USE_OPENMP
#pragma omp atomic
#endif
        contResTot += cv->contResidualTotal;

#ifdef USE_OPENMP
#pragma omp atomic
#endif
        inflowTotal += cv->totalExternalInflow;

#ifdef USE_OPENMP
#pragma omp atomic
#endif
        outflowTotal += cv->totalExternalOutflow;
      }

      externalInflowVar.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, (size_t)m_numCells}), externalInflowValues);
      externalOutflowVar.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, (size_t)m_numCells}), externalOutflowValues);

      depthVar.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, (size_t)m_numCells}), depthValues);
      maxdepthVar.putVar(std::vector<size_t>({0}), std::vector<size_t>({(size_t)m_numCells}), maxDepthValues);

      elevVar.putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1,(size_t)m_numCells}), elevValues);
      maxelevVar.putVar(std::vector<size_t>({0}), std::vector<size_t>({(size_t)m_numCells}), maxElevValues);

      uVar.putVar(std::vector<size_t>({currentTime, 0}) , std::vector<size_t>({1,(size_t)m_numCells}) , uVelValues);
      uResVar.putVar(std::vector<size_t>({currentTime, 0}) , std::vector<size_t>({1,(size_t)m_numCells}) , uVelResValues);
      maxUVar.putVar(std::vector<size_t>({0}), std::vector<size_t>({(size_t)m_numCells}), maxUVelValues);

      vVar.putVar(std::vector<size_t>({currentTime, 0}) , std::vector<size_t>({1,(size_t)m_numCells}) , vVelValues);
      vResVar.putVar(std::vector<size_t>({currentTime, 0}) , std::vector<size_t>({1,(size_t)m_numCells}) , vVelResValues);
      maxVVar.putVar(std::vector<size_t>({0}) , std::vector<size_t>({(size_t)m_numCells}),maxVVelValues);

      continuityResVar.putVar(std::vector<size_t>({currentTime, 0}) , std::vector<size_t>({1,(size_t)m_numCells}) , contResValues);
      pressureResVar.putVar(std::vector<size_t>({currentTime, 0}) , std::vector<size_t>({1,(size_t)m_numCells}) , pressResValues);

      continuityResTotalVar.putVar(std::vector<size_t>({currentTime}) , std::vector<size_t>({1}) , &contResTot);
      inflowTotalVar.putVar(std::vector<size_t>({currentTime}) , std::vector<size_t>({1}) , &inflowTotal);
      outflowTotalVar.putVar(std::vector<size_t>({currentTime}) , std::vector<size_t>({1}) , &outflowTotal);

      edgeFlowVar.putVar(std::vector<size_t>({currentTime,0,0}) , std::vector<size_t>({1, (size_t)m_numCells, 3}), edgeFlowValues);

      delete[] externalInflowValues;
      delete[] edgeFlowValues;
      delete[] depthValues;
      delete[] maxDepthValues;
      delete[] elevValues;
      delete[] maxElevValues;
      delete[] uVelValues;
      delete[] maxUVelValues;
      delete[] uVelResValues;
      delete[] vVelValues;
      delete[] maxVVelValues;
      delete[] vVelResValues;
      delete[] contResValues;
      delete[] pressResValues;

      m_writeFrequencyCounter++;

      //    if(m_writeFrequencyCounter >= m_writeFrequency)
      {
        m_outputNetCDF->sync();
      }
    }

  }
  catch(NcException &e)
  {
    printf("NetCDF FileIO Error: %s", e.what());
    return ;
  }
}

void FVHMComponent::writeToCSV()
{
  if(m_CSVOutputIO.isOpen())
  {

    for(int i = 0; i < m_numCells ; i++)
    {
      TriCV *cv = m_controlVolumes[i];

      if((m_writeActiveCellsOnly && cv->wetIndex) || !m_writeActiveCellsOnly)
      {
        double inflow = cv->inflowOutflow;

        for(int j = 0; j < cv->numEdges; j++)
        {
          double flow = cv->faceNormalVels[j].value * cv->faceDepths[j].value * cv->r_eta[j] ;

          if(cv->faceNormalVels[j].isBC)
          {
            inflow -= flow;
          }
        }

        HCVertex *v1 = cv->vertices[0];
        HCVertex *v2 = cv->vertices[1];
        HCVertex *v3 = cv->vertices[2];

        double z1 = cv->nWSE[0];
        double z2 = cv->nWSE[1];
        double z3 = cv->nWSE[2];

        //        if(cv->wetIndex == 0)
        {
          z1 = std::max(cv->nz[0], z1);
          z2 = std::max(cv->nz[1], z2);
          z3 = std::max(cv->nz[2], z3);
        }

        m_CSVOutputTextStream << "\"POLYGON Z ((" << v1->x() << " " << v1->y() << " " << z1 << ","
                              << v2->x() << " " << v2->y() << " " << z2 << ","
                              << v3->x() << " " << v3->y() << " " << z3 << ","
                              << v1->x() << " " << v1->y() << " " << z1 << "))\","
                              << cv->cell->index() << ","
                              << cv->center->v[0] << ","
                              << cv->center->v[1] << ","
                              << cv->cz << ","
                              << "\"" << m_qtDateTime.toString(Qt::ISODate).replace("-","").replace(":","").replace("T","") << "\","
                              << m_timeStepCount << ","
                              << cv->h->value << ","
                              << cv->z->value << ","
                              << cv->vel[0].value << ","
                              << cv->vel[1].value << ","
                              << cv->faceNormalVels[0].value << ","
                              << cv->faceNormalVels[1].value << ","
                              << cv->faceNormalVels[2].value << ","
                              << inflow << ","
                              << cv->velResidual[0] <<  ","
                              << cv->velResidual[1] << ","
                              << cv->contResidualTotal <<  endl;

      }
    }

    {
      m_CSVOutputTextStream.flush();
    }
  }
}

void FVHMComponent::writeToLogFile(int iters,
                                   double uInitRes, double uFinRes,
                                   double vInitRes, double vFinRes,
                                   double contInitRes, double contFinRes,
                                   double pressInitRes, double pressFinRes)
{
  if(m_logFileIO.isOpen())
  {

    m_logFileTextStream << QString::number(m_currentDateTime,'e',15) << ","
                        << m_timeStep << ","
                        << iters << ","
                        << uInitRes << ","
                        << uFinRes << ","
                        << vInitRes << ","
                        << vFinRes << ","
                        << contInitRes << ","
                        << contFinRes << ","
                        << pressInitRes << ","
                        << pressFinRes << ","
                        << endl;

    m_logFileIO.flush();
  }
}

void FVHMComponent::closeOutputCSVFile()
{
  if(m_CSVOutputIO.isOpen())
  {
    m_CSVOutputIO.flush();
    m_CSVOutputIO.close();
  }
}

void FVHMComponent::closeLogFile()
{
  if(m_logFileIO.isOpen())
  {
    m_logFileIO.flush();
    m_logFileIO.close();
  }
}




