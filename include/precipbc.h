/*! \file   precipbc.h
 *  \author Caleb Amoa Buahin <caleb.buahin@gmail.com>
 *  \version   1.0.0.0
 *  \section   Description
 *  \section License
 *  cvwseoutput.h, associated files and libraries are free software;
 *  you can redistribute it and/or modify it under the terms of the
 *  Lesser GNU General Public License as published by the Free Software Foundation;
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
