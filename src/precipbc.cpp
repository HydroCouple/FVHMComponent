#include "stdafx.h"
#include "precipbc.h"
#include "fvhmcomponent.h"
#include "spatial/octree.h"
#include "temporal/timedata.h"
#include "spatial/polygon.h"
#include "spatial/point.h"
#include "controlvolume.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif

using namespace HydroCouple;
using namespace HydroCouple::Spatial;
using namespace SDKTemporal;
using namespace std;

PrecipBC::PrecipBC(const QString &id,
                   Dimension *timeDimension,
                   Dimension *geometryDimension,
                   ValueDefinition *valueDefinition,
                   FVHMComponent *modelComponent)
  : TimeGeometryArgumentDouble(id, IGeometry::Polygon,
                               timeDimension, geometryDimension,
                               valueDefinition,modelComponent),
    m_modelComponent(modelComponent)

{

}

PrecipBC::~PrecipBC()
{

}

void PrecipBC::findAssociatedCVGeometries()
{
  Octree *octree = new Octree(Octree::Octree2D, Octree::AlongEnvelopes, 7, 1000);

  for(QSharedPointer<HCGeometry> geometry : m_geometries)
  {
    HCPolygon *polygon = dynamic_cast<HCPolygon*>(geometry.data());
    m_controlVolumes[polygon] = vector<int>();
    octree->addGeometry(polygon);
  }

  std::vector<TriCV*> & controlVolumes = m_modelComponent->m_controlVolumes;

  for(TriCV *cv : controlVolumes)
  {
    HCPoint center(cv->center->v[0], cv->center->v[1]);

    std::vector<IGeometry*> geometries =  octree->findCollidingGeometries(&center);

    for(IGeometry* geometry : geometries)
    {
      HCPolygon* polygon = dynamic_cast<HCPolygon*>(geometry);

      if(polygon->contains(&center))
      {
        vector<int> & cvols = m_controlVolumes[polygon];
        cvols.push_back(cv->index);
        break;
      }
    }
  }

  delete octree;
}

void PrecipBC::prepare()
{

}

void PrecipBC::clear()
{
  m_controlVolumes.clear();
  timeCursor()->resetCursor();
}

void PrecipBC::applyBoundaryConditions(double dateTime, double prevTimeStep)
{

  if(m_controlVolumes.size())
  {
    std::vector<TriCV*> & controlVolumes = m_modelComponent->m_controlVolumes;

    if(m_times.size())
    {
      int index = findDateTimeIndex(dateTime);

      if(index > -1)
      {
        double dtu = m_times[index]->modifiedJulianDay();

        if(dtu == dateTime)
        {
          double precipTotalInflow = 0;


          for(size_t i = 0; i < m_geometries.size() ; i++)
          {
            HCPolygon* polygon = dynamic_cast<HCPolygon*>(m_geometries[i].data());

            double value = (*this)(index ,i);

            vector<int> & cvols = m_controlVolumes[polygon];


#ifdef USE_OPENMP
#pragma omp parallel for
#endif
            for(size_t g = 0; g < cvols.size(); g++)
            {
              TriCV *cv = controlVolumes[cvols[g]];
              double inflowRate = value * cv->area * 2.7777777778e-7;
              cv->inflowOutflow += inflowRate;

#ifdef USE_OPENMP
#pragma omp atomic
#endif
              precipTotalInflow += inflowRate;
            }
          }

          if(m_modelComponent->m_printFrequencyCounter >= m_modelComponent->m_printFrequency)
          {
            printf("FVHM Precip Total Q: %f\n",precipTotalInflow);
          }

        }
        else if(index - 1 > -1)
        {
          double dtl = m_times[index - 1]->modifiedJulianDay();
          double factor  = (dateTime - dtl) / (dtu - dtl);
          double precipTotalInflow = 0.0;

          for(size_t i = 0; i < m_geometries.size() ; i++)
          {
            HCPolygon *polygon = dynamic_cast<HCPolygon*>(m_geometries[i].data());
            double valueu = (*this)(index,i);
            double valuel = (*this)(index - 1,i);

            double value = valuel + factor * (valueu - valuel);

            vector<int> & cvols = m_controlVolumes[polygon];

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
            for(size_t g = 0; g < cvols.size(); g++)
            {
              TriCV *cv = controlVolumes[cvols[g]];
              double inflowRate = value * cv->area * 2.7777777778e-7;
              cv->inflowOutflow += inflowRate;

#ifdef USE_OPENMP
#pragma omp atomic
#endif
              precipTotalInflow += inflowRate;
            }
          }

          if(m_modelComponent->m_printFrequencyCounter >= m_modelComponent->m_printFrequency)
          {
            printf("FVHM Precip Total Q: %f\n",precipTotalInflow);
          }

        }
      }
    }
  }
}

int PrecipBC::findDateTimeIndex(double dateTime)
{
  int index = -1;
  int beginCursor = timeCursor()->index();

  do
  {
    if(m_times[timeCursor()->index()]->modifiedJulianDay() >= dateTime)
    {
      index = timeCursor()->index();
      break;
    }

  }while( timeCursor()->moveNext() != beginCursor);

  return index;
}
