/*!
 *  \file    fvhmcompopnentinfo.h
 *  \author  Caleb Amoa Buahin <caleb.buahin@gmail.com>
 *  \version 1.0.0
 *  \section Description
 * Main code for the FVHMComponent model.
 *  \section License
 *  fvhmcompopnentinfo.h, associated files and libraries are free software;
 *  you can redistribute it and/or modify it under the terms of the
 *  Lesser GNU Lesser General Public License as published by the Free Software Foundation;
 *  either version 3 of the License, or (at your option) any later version.
 *  fvhmcompopnentinfo.h its associated files is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.(see <http://www.gnu.org/licenses/> for details)
 *  \date 2014-2016
 *  \pre
 *  \bug
 *  \todo
 *  \warning
 */

#ifndef FVHMCOMPONENTINFO
#define FVHMCOMPONENTINFO

#include "fvhmcomponent_global.h"
#include "core/abstractmodelcomponentinfo.h"

class FVHMCOMPONENT_EXPORT FVHMComponentInfo : public AbstractModelComponentInfo
{
      Q_OBJECT
      Q_PLUGIN_METADATA(IID "FVHMComponentInfo")

   public:

      FVHMComponentInfo(QObject *parent = nullptr);

      virtual ~FVHMComponentInfo();

      HydroCouple::IModelComponent* createComponentInstance() override;

      void createAdaptedOutputFactories();

};

Q_DECLARE_METATYPE(FVHMComponentInfo*)

#endif // FVHMCOMPONENTINFO

