/*! \file   edgefluxesio.h
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

#ifndef EDGEFLUXESIO_H
#define EDGEFLUXESIO_H

#include "fvhmcomponent_global.h"
#include "edgebc.h"
#include "controlvolume.h"

#include <QHash>
#include <QSet>
#include <tuple>
#include <QTextStream>

class FVHMComponent;
class HCGeometry;
class Edge;
struct TriCV;
struct Vect;


class FVHMCOMPONENT_EXPORT EdgeFluxesIO  : public EdgeBC
{
    Q_OBJECT

  public:

    EdgeFluxesIO(const QString &id,
                           Dimension *geometryDimension,
                           ValueDefinition *valueDefinition,
                           FVHMComponent *modelComponent);

    virtual ~EdgeFluxesIO();

    void prepare() override;

    void applyBoundaryConditions(double dateTime, double prevTimeStep) override;

    void setOutputFile(const QFileInfo& file);

    QFileInfo outputFile() const;

  private:
    double m_missingData;
    double m_defaultValue;
    QFileInfo m_outputFile;
    QFile m_edgeFluxesFile;
    QTextStream m_edgeFluxesTextStream;
};

#endif // EDGEFLUXESIO_H
