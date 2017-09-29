#include "stdafx.h"
#include "edgefluxesio.h"
#include "spatial/geometry.h"
#include "spatial/edge.h"
#include "spatial/polygon.h"
#include "core/valuedefinition.h"
#include "fvhmcomponent.h"
#include "controlvolume.h"
#include "spatial/point.h"
#include <QFile>
#include <QTextStream>

using namespace HydroCouple;
using namespace HydroCouple::Spatial;

using namespace std;

EdgeFluxesIO::EdgeFluxesIO(const QString &id,
                           Dimension *geometryDimension,
                           ValueDefinition *valueDefinition,
                           FVHMComponent *modelComponent)
  :EdgeBC(id, geometryDimension,valueDefinition,modelComponent)
{

}

EdgeFluxesIO::~EdgeFluxesIO()
{
  if(m_edgeFluxesFile.isOpen())
  {
    m_edgeFluxesFile.flush();
    m_edgeFluxesFile.close();
  }
}


void EdgeFluxesIO::prepare()
{
  m_missingData = valueDefinition()->missingValue().toDouble();
  m_defaultValue = valueDefinition()->defaultValue().toDouble();

  const std::vector<TriCV*>& controlVolumes = m_modelComponent->m_controlVolumes;

  for(unordered_map<HCGeometry*,set<Edge*>>::iterator it = m_edges.begin() ; it != m_edges.end() ; it++)
  {
    for(Edge *edge : (*it).second)
    {
      HCPolygon *tri = edge->leftInternal();
      TriCV *cv = controlVolumes[tri->index()];
      int index = edge->marker();
      cv->faceDepths[index].isBC = true;
      cv->hasEdgeDepthBC = true;
      FaceNormVelBC &faceVel = cv->faceNormalVels[index];
      faceVel.isBC = true;
      faceVel.calculateWallShearStress = false;
      faceVel.initialiazeVelocityVariable();
      faceVel.velocityCalculateMode = FaceNormVelBC::PreCalculated;
    }
  }

  if(!m_outputFile.fileName().isEmpty() &&
     !m_outputFile.fileName().isNull() && m_outputFile.absoluteDir().exists())
  {
    m_edgeFluxesFile.setFileName(m_outputFile.absoluteFilePath());


    if(m_edgeFluxesFile.open(QIODevice::WriteOnly | QIODevice::Truncate))
    {
      m_edgeFluxesTextStream.setDevice(&m_edgeFluxesFile);
      m_edgeFluxesTextStream.setRealNumberPrecision(10);
      m_edgeFluxesTextStream.setRealNumberNotation(QTextStream::SmartNotation);

      m_edgeFluxesTextStream << "DateTime";

      int count = 0;

      for(unordered_map<HCGeometry*,set<Edge*>>::iterator it = m_edges.begin() ; it != m_edges.end() ; it++)
      {
        QString name = "Geom" + QString::number(count);
        m_edgeFluxesTextStream << "," << name << "Flow" ;
        count++;
      }

      m_edgeFluxesTextStream << endl;
    }
  }
}

void EdgeFluxesIO::applyBoundaryConditions(double dateTime, double prevTimeStep)
{
  if(m_edgeFluxesFile.isOpen())
  {
    const std::vector<TriCV*>& controlVolumes = m_modelComponent->m_controlVolumes;

    m_edgeFluxesTextStream << m_modelComponent->m_qtDateTime.toString(Qt::ISODate);

    for(int i = 0; i < m_geometries.length() ; i++)
    {
      double totalQ = 0.0;

      HCGeometry* geometry = m_geometries[i].data();
      set<Edge*> & edges = m_edges[geometry];

      for(Edge *edge : edges)
      {
        HCPolygon *tri = edge->leftInternal();
        TriCV *cv = controlVolumes[tri->index()];
        int face = edge->marker();
        totalQ += cv->faceNormalVels[face].value * cv->r_eta[face] * cv->faceDepths[face].value;
      }

      m_edgeFluxesTextStream << "," << totalQ;
    }

    m_edgeFluxesTextStream << endl;
  }
}

void EdgeFluxesIO::setOutputFile(const QFileInfo &file)
{
  m_outputFile = file;
}

QFileInfo EdgeFluxesIO::outputFile() const
{
  return m_outputFile;
}

