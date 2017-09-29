#ifndef OUTLETWSESLOPE_H
#define OUTLETWSESLOPE_H

#include "edgebc.h"

struct  TriCV;

class FVHMCOMPONENT_EXPORT OutletWSESlope : public EdgeBC
{
    Q_OBJECT

  public:

    OutletWSESlope(const QString &id,
                   Dimension *geometryDimension,
                   ValueDefinition *valueDefinition,
                   FVHMComponent *modelComponent);

    virtual ~OutletWSESlope();

    void prepare() override;

    void applyBoundaryConditions(double dateTime, double prevTimeStep) override;

    void calculateFaceVelocity(TriCV *cv, int face);

  private:
    double m_missingData;
    double m_defaultValue;
};

#endif // OUTLETWSESLOPE_H
