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
