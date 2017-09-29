#include "stdafx.h"
#include "initialwsebc.h"
#include "fvhmcomponent.h"
#include "controlvolume.h"

#include "core/dimension.h"
#include "core/unit.h"
#include "core/valuedefinition.h"
#include "spatial/polyhedralsurface.h"
#include "spatial/octree.h"
#include "spatial/point.h"
#include "spatial/polygon.h"

using namespace HydroCouple::Spatial;

InitialWSEBC::InitialWSEBC(const QString &id, FVHMComponent *modelComponent):
  GeometryArgumentDouble(id,
                         IGeometry::PointZ,
                         new Dimension("GeometryDimension",modelComponent),
                         new Quantity(QVariant::Double, Unit::lengthInMeters(modelComponent),modelComponent),
                         modelComponent),
  m_waterSurfaceTIN(nullptr),
  m_modelComponent(modelComponent)
{

}

InitialWSEBC::~InitialWSEBC()
{
  clear();
}

void InitialWSEBC::findAssociatedCVGeometries()
{
  if(m_waterSurfaceTIN != nullptr)
    clear();

  //void triangulate and find
  std::vector<HCPoint*> points(m_geometries.size());
  std::vector<HCLineString*> plsg(0);
  std::vector<HCPoint*> holes(0);

  if(m_geometries.length())
  {
    for(int i = 0; i <  m_geometries.length(); i++)
    {
      HCPoint *point = dynamic_cast<HCPoint*>(m_geometries[i].data());
      points[i] = point;
    }

    m_waterSurfaceTIN = HCTIN::triangulateSFS("zcDYQ",points, plsg, holes);

    Octree octree(Octree::Octree2D, Octree::AlongEnvelopes,8,1000);

    for(int i = 0; i < m_waterSurfaceTIN->patchCount(); i++)
    {
      octree.addGeometry(m_waterSurfaceTIN->triangle(i));
    }

    std::vector<TriCV*>& controlVolumes = m_modelComponent->m_controlVolumes;

    for(size_t i = 0; i < controlVolumes.size(); i++)
    {
      TriCV *cv = controlVolumes[i];
      HCPoint *p = new HCPoint(cv->center->v[0],cv->center->v[1]);
      HCTriangle *triangle = nullptr;

      if((triangle = findInOctree(p, &octree)) != nullptr)
      {
        m_controlVolumeMapping[triangle].insert(cv);
      }

      delete p;
    }
  }
}

void InitialWSEBC::prepare()
{

}

void InitialWSEBC::applyBoundaryConditions(double dateTime, double prevTimeStep)
{
  if(dateTime == m_modelComponent->currentDateTime())
  {
    std::unordered_map<HCTriangle*,std::set<TriCV*>>::iterator it;

    for(it = m_controlVolumeMapping.begin(); it != m_controlVolumeMapping.end(); it++)
    {
      HCTriangle *triangle = it->first;
      std::set<TriCV*> & cvs = it->second;

      Vect v1(*triangle->edge()->orig());
      Vect v2(*triangle->edge()->dest());
      Vect v3(*triangle->edge()->leftNext()->dest());

      Vect v11 = v3 - v2;
      Vect v22 = v1 - v2;

      Vect n =  Vect::crossProduct(v11,v22);
      n.normalize();
      double d = -n.v[0] * v1.v[0] - n.v[1] * v1.v[1] - n.v[2] * v1.v[2];

      for(TriCV *cv : cvs)
      {
        double z = (n.v[0] * cv->center->v[0] + n.v[1] * cv->center->v[1] + d) / -n.v[2];
        cv->setVFRWSE(z);
      }
    }
  }
}

void InitialWSEBC::clear()
{
  if(m_waterSurfaceTIN != nullptr)
  {
    delete m_waterSurfaceTIN;
    m_waterSurfaceTIN = nullptr;
  }

  m_controlVolumeMapping.clear();
}

HCTriangle *InitialWSEBC::findInOctree(HCPoint *point, Octree *octree)
{
  std::vector<IGeometry*> eligibleGeoms = octree->findCollidingGeometries(point);

  for(IGeometry *geom : eligibleGeoms)
  {
    if(geom->contains(point))
    {
      return dynamic_cast<HCTriangle*>(geom);
    }
  }

  return nullptr;
}

