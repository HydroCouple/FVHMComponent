#include "stdafx.h"
#include "fvhmtimeseriesexchangeitems.h"
#include "fvhmcomponent.h"
#include "core/dimension.h"
#include "temporal/timedata.h"
#include <QDebug>

using namespace SDKTemporal;

FVHMNodeWSETimeSeriesOutput::FVHMNodeWSETimeSeriesOutput(TNode* node,
                                                         Dimension *dimension,
                                                         const QList<SDKTemporal::Time*>& times,
                                                         ValueDefinition *valueDefinition,
                                                         FVHMComponent *component)
  : TimeSeriesOutputDouble(QString(node->ID) + "-N-WSE-Out", times, dimension , valueDefinition, component)
{
  m_component = component;
  m_node = node;
}

FVHMNodeWSETimeSeriesOutput::~FVHMNodeWSETimeSeriesOutput()
{

}

void FVHMNodeWSETimeSeriesOutput::update(HydroCouple::IInput *querySpecifier)
{
  ITimeExchangeItem* timeExchangeItem = dynamic_cast<ITimeExchangeItem*>(querySpecifier);

  if(timeExchangeItem)
  {
    QList<HydroCouple::Temporal::ITime*> inpTimes = timeExchangeItem->times();
    HydroCouple::Temporal::ITime *lastTime = inpTimes[inpTimes.length() -1];

    while (m_component->currentDateTime()->dateTime() < lastTime->dateTime() &&
           m_component->status() == HydroCouple::Updated)
    {
      m_component->update();
    }
  }
  else
  {
    while (m_component->status() != HydroCouple::Done &&
           m_component->status() != HydroCouple::Failed &&
           m_component->status() != HydroCouple::Finished)
    {
      m_component->update();
    }
  }
}

void FVHMNodeWSETimeSeriesOutput::retrieveDataFromModel()
{
  int timeDimLength = timeDimension()->length();
  QList<Time*> ctimes = timesInternal();
  Time *lastTime = ctimes[timeDimLength -1];
  double elev = m_node->invertElev + m_node->newDepth;

  if(m_component->currentDateTime()->dateTime() > lastTime->dateTime())
  {
    if(timeDimLength > 1)
    {
      double values[timeDimLength -1];

      getValues(1,timeDimLength - 1,values);
      setValues(0,timeDimLength - 1,values);

      for(int i = 0 ; i < timeDimLength -1; i++)
      {
        ctimes[i]->setDateTime(ctimes[i+1]->qDateTime());
      }
    }
    setValues(timeDimLength -1, 1, &elev);
    lastTime->setDateTime(m_component->currentDateTime()->qDateTime());
  }
  else if(m_component->currentDateTime()->dateTime() == lastTime->dateTime())
  {
    setValues(timeDimLength -1, 1, &elev);
  }

  QList<HydroCouple::IAdaptedOutput*> tadaptedOutputs = adaptedOutputs();

  for(HydroCouple::IAdaptedOutput* adaptedOutput :tadaptedOutputs)
  {
    adaptedOutput->refresh();
  }
}

//======================================================================================================================================================================

FVHMNodeWSETimeSeriesInput::FVHMNodeWSETimeSeriesInput(TNode* node,
                                                       Dimension *dimension,
                                                       const QList<SDKTemporal::Time*>& times,
                                                       ValueDefinition *valueDefinition,
                                                       FVHMComponent *component)
  : TimeSeriesInputDouble(QString(node->ID) + "-N-WSE-Inp", times, dimension, valueDefinition, component)
{
  m_component = component;
  m_node = node;
}

FVHMNodeWSETimeSeriesInput::~FVHMNodeWSETimeSeriesInput()
{

}


bool FVHMNodeWSETimeSeriesInput::canConsume(HydroCouple::IOutput* provider, QString& message) const
{
  return true;
}


void FVHMNodeWSETimeSeriesInput::retrieveOuputItemData()
{
  HydroCouple::IOutput *provider = this->provider();

  int timeDimLength = timeDimension()->length();
  QList<Time*> ctimes = timesInternal();
  Time *lastTime = ctimes[timeDimLength -1];

  if(m_component->currentDateTime()->dateTime() > lastTime->dateTime())
  {
    if(timeDimLength > 1)
    {
      double values[timeDimLength -1];
      getValues(1,timeDimLength - 1,values);
      setValues(0,timeDimLength - 1,values);

      for(int i = 0 ; i < timeDimLength -1; i++)
      {
        ctimes[i]->setDateTime(ctimes[i+1]->qDateTime());
      }
    }

    lastTime->setDateTime(m_component->currentDateTime()->qDateTime());
  }

  provider->update(this);


  double value = 0;

  ITimeSeriesExchangeItem *tsoutput = dynamic_cast<ITimeSeriesExchangeItem*>(provider);

  if(tsoutput)
  {
    tsoutput->getValues(tsoutput->timeDimension()->length() -1, 1 , &value);
  }

  int index = project_findObject(m_component->SWMMProject(), NODE, m_node->ID);

  if(value && value > m_node->invertElev)
  {
    setValues(timeDimLength-1,1,&value);
    qDebug() << "depth" << value - m_node->invertElev;
    addNodeDepth(m_component->SWMMProject(),index,value - m_node->invertElev);
  }
  else
  {
    setValues(timeDimLength-1,1,&m_node->invertElev);
    removeNodeLateralInflow(m_component->SWMMProject(),index);
  }
}

//======================================================================================================================================================================

FVHMNodeLatInflowTimeSeriesInput::FVHMNodeLatInflowTimeSeriesInput(TNode* node,
                                                                   Dimension *dimension,
                                                                   const QList<SDKTemporal::Time*>& times,
                                                                   ValueDefinition *valueDefinition,
                                                                   FVHMComponent *component)
  : TimeSeriesMultiInputDouble(QString(node->ID)+"-N-Lat-Inf", times, dimension, valueDefinition, component)
{
  m_component = component;
  m_node = node;
}

FVHMNodeLatInflowTimeSeriesInput::~FVHMNodeLatInflowTimeSeriesInput()
{

}

bool FVHMNodeLatInflowTimeSeriesInput::canConsume(HydroCouple::IOutput* provider, QString& message) const
{
  return true;
}

void FVHMNodeLatInflowTimeSeriesInput::retrieveOuputItemData()
{
  QList<HydroCouple::IOutput*> inproviders = providers();
  int timeDimLength = timeDimension()->length();
  QList<Time*> ctimes = timesInternal();
  Time *lastTime = ctimes[timeDimLength -1];

  if(m_component->currentDateTime()->dateTime() > lastTime->dateTime())
  {
    if(timeDimLength > 1)
    {
      double values[timeDimLength -1];
      getValues(1,timeDimLength - 1,values);
      setValues(0,timeDimLength - 1,values);

      for(int i = 0 ; i < timeDimLength -1; i++)
      {
        ctimes[i]->setDateTime(ctimes[i+1]->dateTime());
      }
    }

    lastTime->setDateTime(m_component->currentDateTime()->dateTime());
  }

  for(HydroCouple::IOutput *output : inproviders)
  {
    output->update(this);
  }

  double value = 0;

  for(HydroCouple::IOutput *inpProvider : inproviders)
  {
    ITimeExchangeItem* inpTProvider = dynamic_cast<ITimeExchangeItem*>(inpProvider);

    if(inpTProvider)
    {
      double ovalue = 0;
      int ind[1]={inpProvider->dimensions()[0]->length() - 1};
      int str[1] = {1};
      inpProvider->getValues(ind, str, &ovalue);
      value = value + ovalue;
    }
  }

  int index = project_findObject(m_component->SWMMProject(), NODE, m_node->ID);

  if(value)
  {
    setValues(timeDimLength-1,1,&value);
    addNodeLateralInflow(m_component->SWMMProject(),index,value);
  }
  else
  {
    removeNodeLateralInflow(m_component->SWMMProject(),index);
  }
}

//======================================================================================================================================================================

FVHMLinkDischargeTimeSeriesOutput::FVHMLinkDischargeTimeSeriesOutput(TLink* link,
                                                                     Dimension *dimension,
                                                                     const QList<SDKTemporal::Time*>& times,
                                                                     ValueDefinition *valueDefinition,
                                                                     FVHMComponent *component)
  : TimeSeriesOutputDouble(QString(link->ID) + "-L-F", times, dimension, valueDefinition, component)
{
  m_component = component;
  m_link = link;
}

FVHMLinkDischargeTimeSeriesOutput::~FVHMLinkDischargeTimeSeriesOutput()
{

}

void FVHMLinkDischargeTimeSeriesOutput::update(HydroCouple::IInput *querySpecifier)
{
  ITimeExchangeItem* timeExchangeItem = dynamic_cast<ITimeExchangeItem*>(querySpecifier);

  if(timeExchangeItem)
  {
    QList<HydroCouple::Temporal::ITime*> inptimes = timeExchangeItem->times();
    HydroCouple::Temporal::ITime *lastTime = inptimes[inptimes.length() -1];

    while (m_component->currentDateTime()->dateTime() < lastTime->dateTime() &&
           m_component->status() == HydroCouple::Updated)
    {
      m_component->update();
    }
  }
  else //run to completion
  {
    while (m_component->status() != HydroCouple::Done &&
           m_component->status() != HydroCouple::Failed &&
           m_component->status() != HydroCouple::Finished)
    {
      m_component->update();
    }
  }

  QList<HydroCouple::IAdaptedOutput*> tadaptedOutputs = adaptedOutputs();

  for(HydroCouple::IAdaptedOutput* adaptedOutput :tadaptedOutputs)
  {
    adaptedOutput->refresh();
  }
  //otherwise run to end;
}

void FVHMLinkDischargeTimeSeriesOutput::retrieveDataFromModel()
{
  int timeDimLength = timeDimension()->length();
  QList<HydroCouple::Temporal::ITime*> times = this->times();
  HydroCouple::Temporal::ITime *lastTime = times[timeDimLength -1];

  if(m_component->currentDateTime()->dateTime() > lastTime->dateTime())
  {
    if(timeDimLength > 1)
    {
      double values[timeDimLength -1];

      getValues(1,timeDimLength - 1,values);
      setValues(0,timeDimLength - 1,values);

      for(int i = 0 ; i < timeDimLength -1; i++)
      {
        times[i]->setDateTime(times[i+1]->dateTime());
      }
    }

    double flow = m_link->newFlow;
    setValues(timeDimLength -1, 1, &flow);
    lastTime->setDateTime(m_component->currentDateTime()->dateTime());
  }
  else if(m_component->currentDateTime()->dateTime() == lastTime->dateTime())
  {
    setValues(timeDimLength -1, 1, &m_link->newFlow);
  }
}
