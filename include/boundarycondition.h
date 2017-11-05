/*! \file   boundarycondition.h
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
