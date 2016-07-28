#ifndef FVHMOBJECTITEMS
#define FVHMOBJECTITEMS

#include "fvhmcomponent_global.h"


class FVHMCOMPONENT_EXPORT FVHMOutputObjectItem
{

   public:
      virtual ~FVHMOutputObjectItem(){}

      virtual void retrieveDataFromModel() = 0;
};


class FVHMCOMPONENT_EXPORT FVHMInputObjectItem
{

   public:
      virtual ~FVHMInputObjectItem(){}

      virtual void retrieveOuputItemData() = 0;
};

#endif // FVHMOUTPUTITEM

