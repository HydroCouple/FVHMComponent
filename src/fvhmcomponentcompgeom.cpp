#include "stdafx.h"
#include "fvhmcomponent.h"
#include "spatial/octree.h"
#include "spatial/tinargument.h"
#include "spatial/polygon.h"
#include "spatial/polyhedralsurface.h"
#include "spatial/edge.h"
#include "controlvolume.h"

#include <assert.h>

using namespace HydroCouple;
using namespace HydroCouple::Spatial;

TINArgumentDouble *FVHMComponent::meshArgument() const
{
  return m_TINMeshArgument;
}

const std::vector<TriCV *> &FVHMComponent::controlVolumes() const
{
  return m_controlVolumes;
}

HCVertex *FVHMComponent::findVertex(HydroCouple::Spatial::IPoint *startPoint)
{
  std::vector<IGeometry*> potential = m_octree->findCollidingGeometries(startPoint);

  for(IGeometry *geometry : potential)
  {
    HCTriangle *triangle = dynamic_cast<HCTriangle*>(geometry);

    Edge *start = triangle->edgeInternal();
    Edge *search = start;

    do
    {
      if(search->origInternal()->equals(startPoint,0.001))
      {
        return search->origInternal();
      }
      else if(search->destInternal()->equals(startPoint,0.001))
      {
        return search->destInternal();
      }

      search = search->leftNextInternal();

    }while(search != start);

  }

  return nullptr;
}

Edge *FVHMComponent::findNextEdge(HCVertex *origin, HydroCouple::Spatial::IPoint *destPoint)
{
  Edge* start = origin->edgeInternal();
  Edge* search = start;

  do
  {
    if(search->destInternal()->equals(destPoint,0.001))
    {
      return search;
    }

    search = search->origNextInternal();

  }while(search != start);

  return nullptr;
}

TriCV *FVHMComponent::findCoincidentCV(HydroCouple::Spatial::IPoint *point)
{
  std::vector<IGeometry*> potential = m_octree->findCollidingGeometries(point);

  for(IGeometry *geometry : potential)
  {
    HCTriangle *triangle = dynamic_cast<HCTriangle*>(geometry);

    if(triangle->contains(point))
    {
      TriCV *cv = m_controlVolumes[triangle->index()];
      assert(cv->cell == triangle);
      return cv;
    }
  }

  return nullptr;
}

Octree *FVHMComponent::triangleOctree() const
{
  return m_octree;
}

void FVHMComponent::resetOctree()
{
  m_octree->clearObjects();

  if(m_TINMeshArgument->TINInternal())
  {
    for(int i = 0 ; i < m_TINMeshArgument->TINInternal()->patchCount() ; i++)
    {
      m_octree->addGeometry(m_TINMeshArgument->TINInternal()->triangleInternal(i));
    }
  }
}

void FVHMComponent::deleteControlVolumes()
{
  m_octree->clearObjects();

  for(TriCV *cv : m_controlVolumes)
    delete cv;

  m_controlVolumes.clear();
}
