#include "stdafx.h"
#include "fvhmcomponentinfo.h"
#include "fvhmcomponent.h"
#include <QUuid>

FVHMComponentInfo::FVHMComponentInfo(QObject *parent)
   :ModelComponentInfo(parent)
{
   setId("FVHMComponent_1.0.0");
   setName("FVHM Component 1.0.0");
   setCaption("Finite Volume Hydrologic Model");
   setIconFilePath("./fvhmicon.png");
   setDescription("A finite volume distributed hydrologic model based on");
   setCategory("Hydrology\\Distributed");
   setCopyright("");
   setVendor("United States Environmental Protection Agency (USEPA)");
   setUrl("https://www.epa.gov/water-research/storm-water-management-model-swmm");
   setEmail("caleb.buahin@aggiemail.usu.edu");
   setVersion("1.0.0.0 ");

   QStringList publications;

   publications << "Gironás, J., L.A. Roesner, J. Davis, L.A. Rossman, and W. Supply, 2009. Storm Water Management Model Applications Manual."
                   "National Risk Management Research Laboratory, Office of Research and Development, US Environmental Protection Agency. "
                   "http://www.epa.gov/NRMRL/pubs/600r09077/600r09077.pdf. Accessed 20 Oct 2014."
                << "Roesner, L.A.; A., Camp, Dresser and McKee, Inc., Annandale, VA, Florida Univ., Gainesville. Dept. of Environmental Engineering Sciences, "
                   "Environmental Research Lab., Athens, GA, J.A.; Aldrich, R.E. Dickinson, and T.O. Barnwell, 1988. "
                   "Storm Water Management Model User’s Manual, Version 4, EXTRAN Addendum."
                   "Environmental Research Laboratory, Office of Research and Development, U.S. Environmental Protection Agency."
                << "Rossman, L.A., 2006. Storm Water Management Model, Quality Assurance Report: Dynamic Wave Flow Routing. "
                   "US Environmental Protection Agency, Office of Research and Development, National Research Management Research Laboratory."
                   " http://www.hydrolatinamerica.org/jahia/webdav/site/hydrolatinamerica/shared/Manuals/epaswmm5_qa/SWMM5_QA.pdf. Accessed 7 Sep 2014."
                << "Rossman, L.A., 2010. Storm Water Management Model User’s Manual, Version 5.0. National Risk Management Research Laboratory, "
                   "Office of Research and Development, US Environmental Protection Agency. "
                   "ftp://152.66.121.2/Oktatas/Epito2000/KozmuhalozatokTervezese-SP2/swmm/epaswmm5_manual.pdf. Accessed 13 Nov 2013.";

   setPublications(publications);

}

FVHMComponentInfo::~FVHMComponentInfo()
{

}

HydroCouple::IModelComponent* FVHMComponentInfo::createComponentInstance()
{
   QString id =  QUuid::createUuid().toString();
   FVHMComponent* component = new FVHMComponent(id, "FVHM Model Instance");
   component->setDescription("FVHM Model Instance");
   component->setComponentInfo(this);
   return component;

  return nullptr;
}

