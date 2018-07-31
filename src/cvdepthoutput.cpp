#include "stdafx.h"
#include "cvdepthoutput.h"
#include "hydrocoupletemporal.h"
#include "temporal/timedata.h"
#include "controlvolume.h"
#include "fvhmcomponent.h"

using namespace HydroCouple;
using namespace HydroCouple::Temporal;


CVDepthOutput::CVDepthOutput(const QString &id,
                             const QSharedPointer<HCTIN> &tinSurface,
                             Dimension *timeDimension,
                             Dimension *cellDimension,
                             Dimension *edgeDimension,
                             Dimension *nodeDimension,
                             ValueDefinition *valueDefinition,
                             FVHMComponent *modelComponent)
  : TimeTINOutputDouble(id,
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

CVDepthOutput::~CVDepthOutput()
{

}

void CVDepthOutput::updateValues(HydroCouple::IInput *querySpecifier)
{
  ITimeComponentDataItem* timeExchangeItem = dynamic_cast<ITimeComponentDataItem*>(querySpecifier);

  if(timeExchangeItem)
  {
   double queryTime = timeExchangeItem->time(timeExchangeItem->timeCount() - 1)->julianDay();

    while (m_FVHMComponent->currentDateTime()->julianDay() < queryTime &&
           m_FVHMComponent->status() == IModelComponent::Updated)
    {
      m_FVHMComponent->update(QList<IOutput*>({this}));
    }
  }
  else
  {
    if(m_FVHMComponent->status() == IModelComponent::Updated)
    {
      m_FVHMComponent->update(QList<IOutput*>({this}));
    }
  }

  refreshAdaptedOutputs();
}

void CVDepthOutput::updateValues()
{
  moveDataToPrevTime();

  int currentTimeIndex = timeCount() - 1;
  m_times[currentTimeIndex]->setJulianDay(m_FVHMComponent->currentDateTime()->julianDay());
  resetTimeSpan();

  const std::vector<TriCV*>& controlVolumes = m_FVHMComponent->controlVolumes();

  for(size_t i = 0; i < controlVolumes.size(); i++)
  {
    TriCV* cv = controlVolumes[i];
    setValue(currentTimeIndex,i,0,0,&cv->h->value);
  }

  refreshAdaptedOutputs();
}
