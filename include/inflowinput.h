/*! \file   inflowinput.h
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

#ifndef INFLOWINPUT_H
#define INFLOWINPUT_H

#include "fvhmcomponent.h"
#include "fvhmcomponent_global.h"
#include "spatiotemporal/timetinexchangeitem.h"

#include <unordered_map>

class HYDROCOUPLESDK_EXPORT InflowInput: public TimeTINMultiInputDouble
{
    Q_OBJECT


  public:

    InflowInput(const QString &id,
                const QSharedPointer<HCTIN> &tinSurface,
                Dimension *timeDimension,
                Dimension *cellDimension,
                Dimension *edgeDimension,
                Dimension *nodeDimension,
                ValueDefinition *valueDefinition,
                FVHMComponent *modelComponent);

    virtual ~InflowInput();

    bool addProvider(HydroCouple::IOutput *provider) override;

    bool removeProvider(HydroCouple::IOutput *provider) override;

    bool canConsume(HydroCouple::IOutput *provider, QString &message) const override;

    void retrieveValuesFromProvider() override;

    void applyData() override;

  private:
    std::unordered_map<int, std::unordered_map<HydroCouple::SpatioTemporal::ITimeGeometryComponentDataItem*, int> > m_controlVolumeMapping;
    FVHMComponent *m_FVHMComponent;
};

#endif // INFLOWINPUT_H
