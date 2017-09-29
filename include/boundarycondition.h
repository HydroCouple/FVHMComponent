#ifndef BOUNDARYCONDITION_H
#define BOUNDARYCONDITION_H

#include <QObject>
#include "fvhmcomponent_global.h"


class IBoundaryCondition
{
  public:

    virtual ~IBoundaryCondition(){}

    virtual void findAssociatedCVGeometries() = 0;

    virtual void prepare() = 0;

    virtual void applyBoundaryConditions(double dateTime, double prevTimeStep) = 0;

    virtual void clear() = 0;

};

Q_DECLARE_INTERFACE(IBoundaryCondition, "FVHMComponent::IBoundaryCondition")
#endif // BOUNDARYCONDITION_H
