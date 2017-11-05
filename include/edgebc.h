/*! \file   edgebc.h
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


#ifndef EDGEBC_H
#define EDGEBC_H

#include "fvhmcomponent_global.h"
#include "spatial/geometryargument.h"
#include "spatiotemporal/timegeometryargument.h"
#include "boundarycondition.h"

#include <unordered_map>
#include <set>

class FVHMComponent;
class HCGeometry;
class Edge;

class FVHMCOMPONENT_EXPORT EdgeBC : public GeometryArgumentDouble,
    public virtual IBoundaryCondition
{
    Q_OBJECT
    Q_INTERFACES(IBoundaryCondition)

  public:

    EdgeBC(const QString &id,
           Dimension *geometryDimension,
           ValueDefinition *valueDefinition,
           FVHMComponent *modelComponent);

    virtual ~EdgeBC();

    void findAssociatedCVGeometries() override;

    void clear() override;

  protected:
    std::unordered_map<HCGeometry*, std::set<Edge*>> m_edges;
    FVHMComponent *m_modelComponent;
};


class FVHMCOMPONENT_EXPORT TimeEdgeBC : public TimeGeometryArgumentDouble,
    public virtual IBoundaryCondition
{
    Q_OBJECT
    Q_INTERFACES(IBoundaryCondition)

  public:

    TimeEdgeBC(const QString &id,
               Dimension *timeDimension,
               Dimension *geometryDimension,
               ValueDefinition *valueDefinition,
               FVHMComponent *modelComponent);

    virtual ~TimeEdgeBC();

    void findAssociatedCVGeometries() override;

    void clear() override;

    int findDateTimeIndex(double dateTime);

  protected:
    std::unordered_map<HCGeometry*, std::set<Edge*>> m_edges;
    FVHMComponent *m_modelComponent;
};




#endif // EDGEBC_H
