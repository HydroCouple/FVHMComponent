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
