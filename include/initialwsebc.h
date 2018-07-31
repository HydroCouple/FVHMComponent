/*! \file   initialwsebc.h
 *  \author Caleb Amoa Buahin <caleb.buahin@gmail.com>
 *  \version   1.0.0
 *  \section   Description
 *  \section License
 *  cvwseoutput.h, associated files and libraries are free software;
 *  you can redistribute it and/or modify it under the terms of the
 *  Lesser GNU Lesser General Public License as published by the Free Software Foundation;
 *  either version 3 of the License, or (at your option) any later version.
 *  The FVHMComponent library and its associated files is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.(see <http://www.gnu.org/licenses/> for details)
 *  \date 2014-2016
 *  \pre
 *  \bug
 *  \warning
 *  \todo
 */

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
