#include "stdafx.h"
#include "criticaldepthoutflowbc.h"
#include "spatial/geometry.h"
#include "spatial/edge.h"
#include "spatial/polygon.h"
#include "core/valuedefinition.h"
#include "fvhmcomponent.h"
#include "controlvolume.h"
#include "spatial/point.h"


using namespace HydroCouple;
using namespace HydroCouple::Spatial;

using namespace std;

CriticalDepthOutflowBC::CriticalDepthOutflowBC(const QString &id,
                                               Dimension *geometryDimension,
                                               ValueDefinition *valueDefinition,
                                               FVHMComponent *modelComponent)
  :EdgeBC(id, geometryDimension,valueDefinition,modelComponent)
{

}

CriticalDepthOutflowBC::~CriticalDepthOutflowBC()
{

}


void CriticalDepthOutflowBC::prepare()
{
  m_missingData = valueDefinition()->missingValue().toDouble();
  m_defaultValue = valueDefinition()->defaultValue().toDouble();

  const std::vector<TriCV*>& controlVolumes = m_modelComponent->controlVolumes();

  for(unordered_map<HCGeometry*,set<Edge*>>::iterator it = m_edges.begin() ; it != m_edges.end() ; it++)
  {
    for(Edge *edge : (*it).second)
    {
      HCPolygon *tri = edge->leftInternal();
      TriCV *cv = controlVolumes[tri->index()];
      int index = edge->marker();
      cv->faceDepths[index].isBC = true;
      cv->hasEdgeDepthBC = true;
      FaceNormVelBC &faceVel = cv->faceNormalVels[index];
      faceVel.isBC = true;
      faceVel.calculateWallShearStress = false;
      faceVel.initialiazeVelocityVariable();
      faceVel.velocityCalculateMode = FaceNormVelBC::PreCalculated;
    }
  }
}

void CriticalDepthOutflowBC::applyBoundaryConditions(double dateTime, double prevTimeStep)
{
  const std::vector<TriCV*>& controlVolumes = m_modelComponent->controlVolumes();

  for(int i = 0; i < m_geometries.length() ; i++)
  {
    HCGeometry* geometry = m_geometries[i].data();
    set<Edge*> & edges = m_edges[geometry];

    //double slopeValue = (*this)[i];

    for(Edge *edge : edges)
    {
      HCPolygon *tri = edge->leftInternal();
      TriCV *cv = controlVolumes[tri->index()];
      int face = edge->marker();

      std::tuple<double, double, double> fv = FVHMComponent::calculateZeroGradientFaceVelocity(cv,face);

      VarBC &faceDepth = cv->faceDepths[face];
      FaceNormVelBC &faceVel = cv->faceNormalVels[face];

      double outflowVel = std::max(0.0, std::get<2>(fv));
      double depth = outflowVel * outflowVel / m_modelComponent->g;

      if(outflowVel > 0 && depth > 1e-3 && depth  < 1e4)
      {
        cv->hasEdgeDepthBC = true;
        faceDepth.isBC = true;
        faceDepth.associatedValue = cv->ecz[face] + depth;
        faceDepth.value = depth;

        faceVel.value = std::get<2>(fv);
        faceVel.vel->v[0] = std::get<0>(fv);
        faceVel.vel->v[1] = std::get<1>(fv);
        faceVel.associatedValue = faceVel.value * depth * cv->r_eta[face];

        if(m_modelComponent->m_printFrequencyCounter >= m_modelComponent->m_printFrequency)
        {
          printf("FVHM QOut: %f, VOut: %f, EdgeDepth: %f\n",faceVel.associatedValue,faceVel.value,faceDepth.value);
        }
      }
      else
      {
        cv->hasEdgeDepthBC = false;
        faceDepth.isBC = false;
        cv->calculateWSEGradient();

        faceVel.vel->v[0] = 0.0;
        faceVel.vel->v[1] = 0.0;
        faceVel.value = 0.0;
        faceVel.associatedValue = 0.0;

        if(m_modelComponent->m_printFrequencyCounter >= m_modelComponent->m_printFrequency)
        {
          printf("FVHM QOut: %f, VOut: %f, EdgeDepth: %f\n",faceVel.associatedValue,faceVel.value,faceDepth.value);
        }
      }
    }
  }
}


