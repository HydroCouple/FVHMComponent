/*!
 *  \file    inletoutletflowbc.h
 *  \author  Caleb Amoa Buahin <caleb.buahin@gmail.com>
 *  \version 1.0.0
 *  \section Description
 * Inflow and outflow boundary condition.
 *  \section License
 *  inletoutletflowbc.h, associated files and libraries are free software;
 *  you can redistribute it and/or modify it under the terms of the
 *  Lesser GNU Lesser General Public License as published by the Free Software Foundation;
 *  either version 3 of the License, or (at your option) any later version.
 *  inletoutletflowbc.h its associated files is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.(see <http://www.gnu.org/licenses/> for details)
 *  \date 2014-2016
 *  \pre
 *  \bug
 *  \todo
 *  \warning
 */

#ifndef INLETOUTLETFLOWBC_H
#define INLETOUTLETFLOWBC_H

#include "fvhmcomponent_global.h"
#include "edgebc.h"

class FVHMComponent;
class HCGeometry;
class HCVertex;

class FVHMCOMPONENT_EXPORT InletOutletFlowBCArgument : public TimeEdgeBC
{
    Q_OBJECT

  public:

    InletOutletFlowBCArgument(const QString &id,
                              Dimension *timeDimension,
                              Dimension *geometryDimension,
                              ValueDefinition *valueDefinition,
                              FVHMComponent *modelComponent);

    virtual ~InletOutletFlowBCArgument();

    void prepare() override;

  protected:

    double m_missingData;
    double m_defaultValue;
};


class FVHMCOMPONENT_EXPORT InletFlowBCArgument : public InletOutletFlowBCArgument
{
    Q_OBJECT

  public:

    InletFlowBCArgument(const QString &id,
                        Dimension *timeDimension,
                        Dimension *geometryDimension,
                        ValueDefinition *valueDefinition,
                        FVHMComponent *modelComponent);

    virtual ~InletFlowBCArgument();

    void prepare() override;

    void applyBoundaryConditions(double dateTime, double prevTimeStep) override;
};

class FVHMCOMPONENT_EXPORT OutletFlowBCArgument : public InletOutletFlowBCArgument
{
    Q_OBJECT

  public:

    OutletFlowBCArgument(const QString &id,
                        Dimension *timeDimension,
                        Dimension *geometryDimension,
                        ValueDefinition *valueDefinition,
                        FVHMComponent *modelComponent);

    virtual ~OutletFlowBCArgument();

    void prepare() override;

    void applyBoundaryConditions(double dateTime, double prevTimeStep) override;
};


#endif // INLETOUTLETFLOWBC_H
