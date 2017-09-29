#ifndef INITIALWSEBC_H
#define INITIALWSEBC_H

#include "fvhmcomponent_global.h"
#include "boundarycondition.h"
#include "spatial/geometryargument.h"

#include <set>
#include <unordered_map>

class FVHMComponent;
class HCGeometry;
struct TriCV;
class HCTIN;
class HCTriangle;
class Octree;
class HCPoint;

class FVHMCOMPONENT_EXPORT InitialWSEBC : public GeometryArgumentDouble,
    public virtual IBoundaryCondition
{
    Q_OBJECT
    Q_INTERFACES(IBoundaryCondition)

  public:

    InitialWSEBC(const QString &id,
                 FVHMComponent *modelComponent);

    virtual ~InitialWSEBC();

    void findAssociatedCVGeometries() override;

    void prepare() override;

    void applyBoundaryConditions(double dateTime, double prevTimeStep) override;

    void clear() override;

  private:

    static HCTriangle *findInOctree(HCPoint *point, Octree *octree);

  private:

    HCTIN *m_waterSurfaceTIN;
    std::unordered_map<HCTriangle*,std::set<TriCV*>> m_controlVolumeMapping;
    FVHMComponent *m_modelComponent;

};

#endif // INITIALWSEBC_H
