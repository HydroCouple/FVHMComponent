#ifndef WSEOUTPUT_H
#define WSEOUTPUT_H

#include "fvhmcomponent_global.h"
#include "spatiotemporal/timetinexchangeitem.h"

class FVHMComponent;

class HYDROCOUPLESDK_EXPORT CVWSEOutput: public TimeTINOutputDouble
{
    Q_OBJECT


  public:

    CVWSEOutput(const QString &id,
              const QSharedPointer<HCTIN> &tinSurface,
              Dimension *timeDimension,
              Dimension *cellDimension,
              Dimension *edgeDimension,
              Dimension *nodeDimension,
              ValueDefinition *valueDefinition,
              FVHMComponent *modelComponent);

    virtual ~CVWSEOutput();

    void updateValues(HydroCouple::IInput *querySpecifier) override;

    void updateValues() override;

  private:

    FVHMComponent *m_FVHMComponent;
};

#endif // WSEOUTPUT_H
