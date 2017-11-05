 /*!
 *  \file    fvhmcomponent_global.h
 *  \author  Caleb Amoa Buahin <caleb.buahin@gmail.com>
 *  \version 1.0.0.0
 *  \section Description
 * Main code for the FVHMComponent model.
 *  \section License
 *  fvhmcompopnentinfo.h, associated files and libraries are free software;
 *  you can redistribute it and/or modify it under the terms of the
 *  Lesser GNU General Public License as published by the Free Software Foundation;
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


#ifndef FVHMCOMPONENT_GLOBAL_H
#define FVHMCOMPONENT_GLOBAL_H

#include <QtCore/qglobal.h>

#if defined(FVHMCOMPONENT_LIBRARY)
#  define FVHMCOMPONENT_EXPORT Q_DECL_EXPORT
#else
#  define FVHMCOMPONENT_EXPORT Q_DECL_IMPORT
#endif

#endif // FVHMCOMPONENT_GLOBAL_H
