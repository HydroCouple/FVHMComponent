#include "stdafx.h"
#include "cvwseoutput.h"
#include "hydrocoupletemporal.h"
#include "temporal/timedata.h"
#include "controlvolume.h"
#include "fvhmcomponent.h"

using namespace HydroCouple;
using namespace HydroCouple::Temporal;


CVWSEOutput::CVWSEOutput(const QString &id,
                     const QSharedPointer<HCTIN> &tinSurface,
                     Dimension *timeDimension,
                     Dimension *cellDimension,
                     Dimension *edgeDimension,
                     Dimension *nodeDimension,
                     ValueDefinition *valueDefinition,
                     FVHMComponent *modelComponent)
  : TimeTINOutputDouble(id,
                        tinSurface,
                        HydroCouple::Spatial::Node,
                        timeDimension,
                        cellDimension,
                        edgeDimension,
                        nodeDimension,
                        valueDefinition,
                        modelComponent),
    m_FVHMComponent(modelComponent)
{

}

CVWSEOutput::~CVWSEOutput()
{

}

void CVWSEOutput::updateValues(HydroCouple::IInput *querySpecifier)
{
  ITimeComponentDataItem* timeExchangeItem = dynamic_cast<ITimeComponentDataItem*>(querySpecifier);

  if(timeExchangeItem)
  {
    double queryTime = timeExchangeItem->time(timeExchangeItem->timeCount() - 1)->modifiedJulianDay();

    while (m_FVHMComponent->currentDateTime() < queryTime &&
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

void CVWSEOutput::updateValues()
{
  moveDataToPrevTime();

  int currentTimeIndex = timeCount() - 1;
  m_times[currentTimeIndex]->setModifiedJulianDay(m_FVHMComponent->currentDateTime());
  resetTimeSpan();

  const std::vector<TriCV*>& controlVolumes = m_FVHMComponent->controlVolumes();

  for(size_t i = 0; i < controlVolumes.size(); i++)
  {
    TriCV* cv = controlVolumes[i];

    for(int j = 0; j < 3; j++)
    {
      for(int k = 0; k < 2; k++)
      {
        setValue(currentTimeIndex,cv->index,j,k,&cv->nWSE[j + k]);
      }
    }
  }

  refreshAdaptedOutputs();
}
