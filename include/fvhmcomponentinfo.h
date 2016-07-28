#ifndef FVHMCOMPONENTINFO
#define FVHMCOMPONENTINFO

#include "fvhmcomponent_global.h"
#include "core/modelcomponentinfo.h"

class FVHMCOMPONENT_EXPORT FVHMComponentInfo : public ModelComponentInfo, public virtual HydroCouple::IModelComponentInfo
{
      Q_OBJECT
      Q_PLUGIN_METADATA(IID "FVHMComponentInfo")
      Q_INTERFACES(HydroCouple::IModelComponentInfo)

   public:
      FVHMComponentInfo(QObject *parent = nullptr);

      virtual ~FVHMComponentInfo();

      HydroCouple::IModelComponent* createComponentInstance() override;

};

Q_DECLARE_METATYPE(FVHMComponentInfo*)

#endif // FVHMCOMPONENTINFO

