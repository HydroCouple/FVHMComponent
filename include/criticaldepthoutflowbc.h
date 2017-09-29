#ifndef CRITICALDEPTHOUTFLOWBC_H
#define CRITICALDEPTHOUTFLOWBC_H


#include "fvhmcomponent_global.h"
#include "edgebc.h"
#include "controlvolume.h"

#include <QHash>
#include <QSet>
#include <tuple>

class FVHMComponent;
class HCGeometry;
class Edge;
struct TriCV;
struct Vect;

class FVHMCOMPONENT_EXPORT CriticalDepthOutflowBC  : public EdgeBC
{
    Q_OBJECT

  public:

    CriticalDepthOutflowBC(const QString &id,
                           Dimension *geometryDimension,
                           ValueDefinition *valueDefinition,
                           FVHMComponent *modelComponent);

    virtual ~CriticalDepthOutflowBC();

    void prepare() override;

    void applyBoundaryConditions(double dateTime, double prevTimeStep) override;

  private:
    double m_missingData;
    double m_defaultValue;
};





#endif // CRITICALDEPTHOUTFLOWBC_H
