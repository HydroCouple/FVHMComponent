/*! \file   cvdepthoutput.h
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

#ifndef CVDEPTHOUTPUT_H
#define CVDEPTHOUTPUT_H

#include "fvhmcomponent_global.h"
#include "spatiotemporal/timetinexchangeitem.h"

class FVHMComponent;

class HYDROCOUPLESDK_EXPORT CVDepthOutput: public TimeTINOutputDouble
{
    Q_OBJECT


  public:

    CVDepthOutput(const QString &id,
                  const QSharedPointer<HCTIN> &tinSurface,
                  Dimension *timeDimension,
                  Dimension *cellDimension,
                  Dimension *edgeDimension,
                  Dimension *nodeDimension,
                  ValueDefinition *valueDefinition,
                  FVHMComponent *modelComponent);

    virtual ~CVDepthOutput();

    void updateValues(HydroCouple::IInput *querySpecifier) override;

    void updateValues() override;

  private:

    FVHMComponent *m_FVHMComponent;
};


#endif // CVDEPTHOUTPUT_H
