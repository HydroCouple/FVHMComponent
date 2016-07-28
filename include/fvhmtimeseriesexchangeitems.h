#ifndef FVHMNODETIMESERIESEXCHANGEITEM
#define FVHMNODETIMESERIESEXCHANGEITEM

#include "fvhmcomponent_global.h"
#include "core/headers.h"
#include "temporal/timeseriesexchangeitem.h"
#include "fvhmobjectitems.h"

class FVHMComponent ;


class FVHMCOMPONENT_EXPORT FVHMNodeWSETimeSeriesOutput : public TimeSeriesOutputDouble,
    public virtual FVHMOutputObjectItem
{
    Q_OBJECT

  public:

    FVHMNodeWSETimeSeriesOutput(TNode* node,
                                Dimension *dimension,
                                const QList<SDKTemporal::Time*>& times,
                                ValueDefinition *valueDefinition,
                                FVHMComponent *component);

    virtual ~FVHMNodeWSETimeSeriesOutput();

    void update(HydroCouple::IInput *querySpecifier) override;

    void retrieveDataFromModel() override;

  private:
    TNode *m_node;
    FVHMComponent *m_component;

};

class FVHMCOMPONENT_EXPORT FVHMNodeWSETimeSeriesInput : public TimeSeriesInputDouble,
    public virtual FVHMInputObjectItem
{
    Q_OBJECT

  public:

    FVHMNodeWSETimeSeriesInput(TNode* node,
                               Dimension *dimension,
                               const QList<SDKTemporal::Time*>& times,
                               ValueDefinition *valueDefinition,
                               FVHMComponent *component);

    virtual ~FVHMNodeWSETimeSeriesInput();

    bool canConsume(HydroCouple::IOutput* provider, QString &message) const override;

    void retrieveOuputItemData() override;

  private:
    TNode *m_node;
    FVHMComponent *m_component;

};

class FVHMCOMPONENT_EXPORT FVHMNodeLatInflowTimeSeriesInput : public TimeSeriesMultiInputDouble,
    public virtual FVHMInputObjectItem
{
    Q_OBJECT

  public:

    FVHMNodeLatInflowTimeSeriesInput(TNode* node,
                                     Dimension *dimension,
                                     const QList<SDKTemporal::Time*>& times,
                                     ValueDefinition *valueDefinition,
                                     FVHMComponent *component);

    virtual ~FVHMNodeLatInflowTimeSeriesInput();

    bool canConsume(HydroCouple::IOutput* provider, QString &message) const override;

    void retrieveOuputItemData() override;

  private:
    TNode *m_node;
    FVHMComponent *m_component;
};


class FVHMCOMPONENT_EXPORT FVHMLinkDischargeTimeSeriesOutput : public TimeSeriesOutputDouble,
    public virtual FVHMOutputObjectItem
{
    Q_OBJECT

  public:

    FVHMLinkDischargeTimeSeriesOutput(TLink* link,
                                      Dimension *dimension,
                                      const QList<SDKTemporal::Time*>& times,
                                      ValueDefinition *valueDefinition,
                                      FVHMComponent *component);

    virtual ~FVHMLinkDischargeTimeSeriesOutput();

    void update(HydroCouple::IInput *querySpecifier) override;

    void retrieveDataFromModel() override;

  private:
    TLink *m_link;
    FVHMComponent *m_component;
};



#endif // FVHMNODETIMESERIESEXCHANGEITEM

