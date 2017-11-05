/*!
 *  \file    wsebc.h
 *  \author  Caleb Amoa Buahin <caleb.buahin@gmail.com>
 *  \version 1.0.0.0
 *  \section Description
 * Water surface elevation boundary condition.
 *  \section License
 *  wsebc.h, associated files and libraries are free software;
 *  you can redistribute it and/or modify it under the terms of the
 *  Lesser GNU General Public License as published by the Free Software Foundation;
 *  either version 3 of the License, or (at your option) any later version.
 *  wsebc.h its associated files is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.(see <http://www.gnu.org/licenses/> for details)
 *  \date 2014-2016
 *  \pre
 *  \bug
 *  \todo
 *  \warning
 */

#ifndef WSEBC_H
#define WSEBC_H

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

class FVHMCOMPONENT_EXPORT WSEBC  : public TimeEdgeBC
{
    Q_OBJECT

  public:

    enum InletOutletType
    {
      None,
      Outlet,
      Inlet,
      OutletInlet
    };

    WSEBC(const QString &id,
          InletOutletType inletOutletType,
          double velRelaxationFactor,
          Dimension *timeDimension,
          Dimension *geometryDimension,
          ValueDefinition *valueDefinition,
          FVHMComponent *modelComponent);

    virtual ~WSEBC();

    void prepare() override;

    void applyBoundaryConditions(double dateTime, double prevTimeStep) override;

    double velocityRelaxationFactor() const;

    void setVelocityRelaxationFactor(double velocityRalaxationFactor);

  private:

    double m_missingData;
    double m_defaultValue;
    double m_velRelaxationFactor = 1.00;
    InletOutletType m_inletOutletType;
};

#endif // WSEBC_H
