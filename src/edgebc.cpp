#include "stdafx.h"
#include "edgebc.h"
#include "fvhmcomponent.h"

#include "core/valuedefinition.h"

#include "spatial/geometry.h"
#include "spatial/edge.h"
#include "spatial/polygon.h"
#include "temporal/timedata.h"


using namespace HydroCouple;
using namespace HydroCouple::Spatial;
using namespace SDKTemporal;


EdgeBC::EdgeBC(const QString &id,
               Dimension *geometryDimension,
               ValueDefinition *valueDefinition,
               FVHMComponent *modelComponent)
  : GeometryArgumentDouble(id, IGeometry::LineString,
                           geometryDimension, valueDefinition, modelComponent),
    m_modelComponent(modelComponent)
{

}

EdgeBC::~EdgeBC()
{
}

void EdgeBC::findAssociatedCVGeometries()
{
  m_edges.clear();

  int count = 0;
  for(QSharedPointer<HCGeometry> geometry : m_geometries)
  {
    count++;

    ILineString *lineString = dynamic_cast<ILineString*>(geometry.data());

    if(lineString->pointCount() > 1)
    {
      IPoint* startPoint = lineString->point(0);

      HCVertex *startVertex = m_modelComponent->findVertex(startPoint);

      if(startVertex)
      {
        for(int i = 1 ; i < lineString->pointCount() ; i++)
        {
          IPoint* endPoint = lineString->point(i);
          Edge *edge = m_modelComponent->findNextEdge(startVertex, endPoint);

          if(edge)
          {
            m_edges[geometry.data()].insert(edge);
            startVertex = edge->destInternal();
          }
          else
          {

          }
        }
      }
    }
    else
    {

    }
  }
}

void EdgeBC::clear()
{
  m_edges.clear();
}



TimeEdgeBC::TimeEdgeBC(const QString &id,
               Dimension *timeDimension,
               Dimension *geometryDimension,
               ValueDefinition *valueDefinition,
               FVHMComponent *modelComponent)
  : TimeGeometryArgumentDouble(id, IGeometry::LineString,
                               timeDimension, geometryDimension,
                               valueDefinition,modelComponent),
    m_modelComponent(modelComponent)
{

}

TimeEdgeBC::~TimeEdgeBC()
{
}

void TimeEdgeBC::findAssociatedCVGeometries()
{
  m_edges.clear();

  for(QSharedPointer<HCGeometry> geometry : m_geometries)
  {
    ILineString *lineString = dynamic_cast<ILineString*>(geometry.data());

    if(lineString->pointCount() > 1)
    {
      IPoint* startPoint = lineString->point(0);

      HCVertex *startVertex = m_modelComponent->findVertex(startPoint);

      if(startVertex)
      {
        for(int i = 1 ; i < lineString->pointCount() ; i++)
        {
          IPoint* endPoint = lineString->point(i);
          Edge *edge = m_modelComponent->findNextEdge(startVertex, endPoint);

          if(edge)
          {
            m_edges[geometry.data()].insert(edge);
            startVertex = edge->destInternal();
          }
        }
      }
    }
  }

  timeCursor()->resetCursor();
}

int TimeEdgeBC::findDateTimeIndex(double dateTime)
{
  int index = -1;
  int beginCursor = timeCursor()->index();

  do
  {
    if(m_times[timeCursor()->index()]->julianDay() >= dateTime)
    {
      index = timeCursor()->index();
      break;
    }

  }while( timeCursor()->moveNext() != beginCursor);

  return index;
}

void TimeEdgeBC::clear()
{
  m_edges.clear();
}
