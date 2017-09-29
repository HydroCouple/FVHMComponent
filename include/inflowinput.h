#ifndef INFLOWINPUT_H
#define INFLOWINPUT_H

#include "fvhmcomponent.h"
#include "fvhmcomponent_global.h"
#include "spatiotemporal/timetinexchangeitem.h"

#include <unordered_map>

class HYDROCOUPLESDK_EXPORT InflowInput: public TimeTINMultiInputDouble
{
    Q_OBJECT


  public:

    InflowInput(const QString &id,
                const QSharedPointer<HCTIN> &tinSurface,
                Dimension *timeDimension,
                Dimension *cellDimension,
                Dimension *edgeDimension,
                Dimension *nodeDimension,
                ValueDefinition *valueDefinition,
                FVHMComponent *modelComponent);

    virtual ~InflowInput();

    bool addProvider(HydroCouple::IOutput *provider) override;

    bool removeProvider(HydroCouple::IOutput *provider) override;

    bool canConsume(HydroCouple::IOutput *provider, QString &message) const override;

    void retrieveValuesFromProvider() override;

    void applyData() override;

  private:
    std::unordered_map<int, std::unordered_map<HydroCouple::SpatioTemporal::ITimeGeometryComponentDataItem*, int> > m_controlVolumeMapping;
    FVHMComponent *m_FVHMComponent;
};

#endif // INFLOWINPUT_H
