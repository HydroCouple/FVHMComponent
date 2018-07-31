#include "stdafx.h"
#include "inflowinput.h"
#include "hydrocouplespatiotemporal.h"
#include "controlvolume.h"
#include "spatial/polyhedralsurface.h"
#include "temporal/timedata.h"
#include "spatial/octree.h"
#include "spatial/polygon.h"

using namespace HydroCouple;
using namespace HydroCouple::Spatial;
using namespace HydroCouple::SpatioTemporal;
using namespace std;

InflowInput::InflowInput(const QString &id,
                         const QSharedPointer<HCTIN> &tinSurface,
                         Dimension *timeDimension,
                         Dimension *cellDimension,
                         Dimension *edgeDimension,
                         Dimension *nodeDimension,
                         ValueDefinition *valueDefinition,
                         FVHMComponent *modelComponent)
  :TimeTINMultiInputDouble(id,
                           tinSurface,
                           HydroCouple::Spatial::Centroid,
                           timeDimension,
                           cellDimension,
                           edgeDimension,
                           nodeDimension,
                           valueDefinition,
                           modelComponent),
    m_FVHMComponent(modelComponent)
{
}

InflowInput::~InflowInput()
{
}

bool InflowInput::addProvider(HydroCouple::IOutput *provider)
{
  if(AbstractMultiInput::addProvider(provider))
  {
    ITimeGeometryComponentDataItem* geometryItem = dynamic_cast<ITimeGeometryComponentDataItem*>(provider);

    int geomCount = geometryItem->geometryCount();

    Octree *octree = new Octree(Octree::Octree2D, Octree::AlongEnvelopes,6,500);
    std::map<IGeometry*,int> indexes ;

    for(int i = 0; i < geomCount; i++)
    {
      IPoint *point = dynamic_cast<IPoint*>(geometryItem->geometry(i));
      octree->addGeometry(point);
      indexes[point] = i;

//      TriCV *cv = m_FVHMComponent->findCoincidentCV(point);

//      if(cv)
//      {
//        m_controlVolumeMapping[cv->index][geometryItem] = i;
//      }
    }

    const std::vector<TriCV*>& controlVolumes = m_FVHMComponent->controlVolumes();

    for(size_t i = 0; i < controlVolumes.size(); i++)
    {
      TriCV *cv = controlVolumes[i];
      std::vector<IGeometry*> collidingGeometries = octree->findCollidingGeometries(cv->cell);

      for(size_t j = 0; j < collidingGeometries.size(); j++)
      {
        IGeometry *geometry = collidingGeometries[j];

        if(cv->cell->contains(geometry))
        {
          m_controlVolumeMapping[cv->index][geometryItem] = indexes[geometry];
          break;
        }
      }
    }


    delete octree;

    return true;
  }

  return false;
}

bool InflowInput::removeProvider(HydroCouple::IOutput *provider)
{
  if(AbstractMultiInput::removeProvider(provider))
  {
    const std::vector<TriCV*> &controlVolumes = m_FVHMComponent->controlVolumes();
    ITimeGeometryComponentDataItem* geometryItem = dynamic_cast<ITimeGeometryComponentDataItem*>(provider);

    for(size_t i = 0; controlVolumes.size(); i++)
    {
      TriCV* cv = controlVolumes[i];

      if(m_controlVolumeMapping.find(cv->index) != m_controlVolumeMapping.end() )
      {
        unordered_map<ITimeGeometryComponentDataItem*, int> & geomIndexMap =  m_controlVolumeMapping[cv->index];

        if(geomIndexMap.find(geometryItem) != geomIndexMap.end())
          geomIndexMap.erase(geometryItem);
      }
    }

    return true;
  }

  return false;
}

bool InflowInput::canConsume(HydroCouple::IOutput *provider, QString &message) const
{
  message = "";

  ITimeGeometryComponentDataItem* geometryItem = dynamic_cast<ITimeGeometryComponentDataItem*>(provider);

  if(geometryItem &&
     (geometryItem->geometryType() == IGeometry::Point ||
      geometryItem->geometryType() == IGeometry::PointM ||
      geometryItem->geometryType() == IGeometry::PointZ ||
      geometryItem->geometryType() == IGeometry::PointZM) &&
     dynamic_cast<HydroCouple::IQuantity*>(geometryItem->valueDefinition()) &&
     geometryItem->valueDefinition()->type() == QVariant::Double)
  {
    return true;
  }

  return false;
}

void InflowInput::retrieveValuesFromProvider()
{
  moveDataToPrevTime();

  size_t currentTimeIndex = m_times.size() - 1;
  m_times[currentTimeIndex]->setJulianDay(m_FVHMComponent->currentDateTime()->julianDay());

  for(IOutput* output : m_providers)
  {
    output->updateValues(this);
  }

  double missingValue = valueDefinition()->missingValue().toDouble();

  const std::vector<TriCV*> &controlVolumes = m_FVHMComponent->controlVolumes();

  unordered_map<int, unordered_map<HydroCouple::SpatioTemporal::ITimeGeometryComponentDataItem*,int> >::iterator it;

  for(it = m_controlVolumeMapping.begin(); it != m_controlVolumeMapping.end(); it++)
  {
    TriCV* cv = controlVolumes[(*it).first];
    unordered_map<HydroCouple::SpatioTemporal::ITimeGeometryComponentDataItem*,int>  &mapping = m_controlVolumeMapping[cv->index];
    std::vector<int> flags(mapping.size(), 0);

    for(int i = 0 ; i < m_providers.length(); i++)
    {
      ITimeGeometryComponentDataItem* geometryItem = dynamic_cast<ITimeGeometryComponentDataItem*>(m_providers[i]);

      if(mapping.find(geometryItem) != mapping.end())
      {
        int lastTimeIndex = geometryItem->timeCount() - 1;
        int srcId = mapping[geometryItem];

        double value = 0;
        geometryItem->getValue(lastTimeIndex,srcId, &value);

        if(flags[i] == 0)
        {
          setValue(currentTimeIndex, cv->index, 0, 0, &value);
          flags[i] = 1;
        }
        else
        {
          double currValue = 0;
          getValue(currentTimeIndex, cv->index, 0, 0, &currValue);

          currValue += value;
          setValue(currentTimeIndex, cv->index, 0, 0, &currValue);
        }
      }
      else
      {
        setValue(currentTimeIndex, i, 0, 0, &missingValue);
      }
    }
  }
}

void InflowInput::applyData()
{
  size_t currentTimeIndex = m_times.size() - 1;
  double missingValue = valueDefinition()->missingValue().toDouble();

  const std::vector<TriCV*> &controlVolumes = m_FVHMComponent->controlVolumes();

  unordered_map<int, unordered_map<HydroCouple::SpatioTemporal::ITimeGeometryComponentDataItem*,int> >::iterator it;

  for(it = m_controlVolumeMapping.begin(); it != m_controlVolumeMapping.end(); it++)
  {
    TriCV* cv = controlVolumes[(*it).first];

    double value = 0.0;
    getValue(currentTimeIndex, cv->index,0,0,&value);

    if(value != missingValue && (value > 0.0 || (value < 0.0 && cv->wetIndex)))
    {
      cv->inflowOutflow += value;
    }
  }
}
