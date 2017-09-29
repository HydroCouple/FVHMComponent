#ifndef PRECIPBC_H
#define PRECIPBC_H

#include "fvhmcomponent_global.h"
#include "boundarycondition.h"
#include "spatiotemporal/timegeometryargument.h"

#include <unordered_map>
#include <set>

class FVHMComponent;
class HCPolygon;

class FVHMCOMPONENT_EXPORT PrecipBC : public TimeGeometryArgumentDouble,
    public virtual IBoundaryCondition
{
    Q_OBJECT
    Q_INTERFACES(IBoundaryCondition)

  public:

    PrecipBC(const QString &id,
             Dimension *timeDimension,
             Dimension *geometryDimension,
             ValueDefinition *valueDefinition,
             FVHMComponent *modelComponent);

    virtual ~PrecipBC();

    void findAssociatedCVGeometries() override;

    void prepare() override;

    void clear() override;

    void applyBoundaryConditions(double dateTime, double prevTimeStep) override;

    int findDateTimeIndex(double dateTime);

  protected:
    std::unordered_map<HCPolygon*, std::vector<int>> m_controlVolumes;
    FVHMComponent *m_modelComponent;

};

#endif // PRECIPBC_H
