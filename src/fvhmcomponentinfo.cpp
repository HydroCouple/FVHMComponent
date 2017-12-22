#include "stdafx.h"
#include "fvhmcomponentinfo.h"
#include "fvhmcomponent.h"
#include "spatial/geometryfactory.h"
#include "temporal/temporalinterpolationfactory.h"
#include <QUuid>

FVHMComponentInfo::FVHMComponentInfo(QObject *parent)
   : AbstractModelComponentInfo(parent)
{
   GeometryFactory::registerGDAL();
   setId("FVHMComponent 1.0.0");
   setCaption("FVHM Model");
   setIconFilePath("./../../resources/images/fvhmicon.png");
   setDescription("A finite volume distributed hydrologic model");
   setCategory("Hydrology\\Distributed");
   setCopyright("Caleb Buahin 2016");
   setVendor("www.hydrocouple.org");
   setUrl("www.hydrocouple.org");
   setEmail("caleb.buahin@aggiemail.usu.edu");
   setVersion("1.0.0.0 ");
   createAdaptedOutputFactories();
}

FVHMComponentInfo::~FVHMComponentInfo()
{

}

HydroCouple::IModelComponent* FVHMComponentInfo::createComponentInstance()
{
   QString id =  QUuid::createUuid().toString();
   FVHMComponent* component = new FVHMComponent(id, "FVHM Model Instance",this);
   component->setDescription("FVHM Model Instance");
   return component;

  return nullptr;
}

void FVHMComponentInfo::createAdaptedOutputFactories()
{
  TemporalInterpolationFactory* tempInterpFactory = new TemporalInterpolationFactory("TimeInterpolationFactory",this);
  addAdaptedOutputFactory(tempInterpFactory);
}
