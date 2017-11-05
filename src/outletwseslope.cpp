#include "stdafx.h"
#include "outletwseslope.h"
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

OutletWSESlope::OutletWSESlope(const QString &id,
                               Dimension *geometryDimension,
                               ValueDefinition *valueDefinition,
                               FVHMComponent *modelComponent)
  :EdgeBC(id, geometryDimension,valueDefinition,modelComponent)
{

}

OutletWSESlope::~OutletWSESlope()
{

}

void OutletWSESlope::prepare()
{
  m_missingData = valueDefinition()->missingValue().toDouble();
  m_defaultValue = valueDefinition()->defaultValue().toDouble();

  std::vector<TriCV*>& controlVolumes = m_modelComponent->m_controlVolumes;

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
      faceVel.calculateFromQ = true;
      faceVel.initialiazeVelocityVariable();
      faceVel.velocityCalculateMode = FaceNormVelBC::PreCalculated;
    }
  }
}

void OutletWSESlope::applyBoundaryConditions(double dateTime, double prevTimeStep)
{

  std::vector<TriCV*>& controlVolumes = m_modelComponent->m_controlVolumes;

  for(int i = 0; i < m_geometries.length() ; i++)
  {
    HCGeometry* geometry = m_geometries[i].data();
    set<Edge*> & edges = m_edges[geometry];

    double slopeValue = (*this)[i];

    for(Edge *edge : edges)
    {
      HCPolygon *tri = edge->leftInternal();
      TriCV *cv = controlVolumes[tri->index()];
      int face = edge->marker();

      if(cv->wetIndex)
      {
        double edgeZ = (slopeValue * cv->r_e_l[face]) + cv->z->value;
        cv->setFaceElevation(face, edgeZ);
        std::tuple<double, double, double> fv = FVHMComponent::calculateZeroGradientFaceVelocity(cv,face);

        FaceNormVelBC &faceVel = cv->faceNormalVels[face];
        VarBC &faceDepth = cv->faceDepths[face];

        faceVel.vel->v[0] = 0.0;
        faceVel.vel->v[1] = 0.0;
        faceVel.value = 0.0;
        faceVel.associatedValue = 0.0;

        if(faceDepth.value > 5e-2)
        {
          if(std::get<2>(fv) > 0.0)
          {
            double factor = faceDepth.value * cv->r_eta[face];

            faceVel.value = std::get<2>(fv);
            faceVel.vel->v[0] = std::get<0>(fv);
            faceVel.vel->v[1] = std::get<1>(fv);
            faceVel.associatedValue = faceVel.value * factor;

            if(m_modelComponent->m_printFrequencyCounter >= m_modelComponent->m_printFrequency)
            {
              printf("FVHM Q: %f\tV: %f\tVx: %f\tVy: %f\tEdgeDepth: %f\n",faceVel.associatedValue,faceVel.value,faceVel.vel->v[0], faceVel.vel->v[1], faceDepth.value);
            }
          }
          else
          {
            faceVel.vel->v[0] = 0.0;
            faceVel.vel->v[1] = 0.0;
            faceVel.value = 0.0;
            faceVel.associatedValue = 0.0;

            if(m_modelComponent->m_printFrequencyCounter >= m_modelComponent->m_printFrequency)
            {
              printf("FVHM Q: %f\tV: %f\tVx: %f\tVy: %f\tEdgeDepth: %f\n",faceVel.associatedValue,faceVel.value,faceVel.vel->v[0], faceVel.vel->v[1], faceDepth.value);
            }
          }
        }
      }
      else
      {
        double edgeZ =  (slopeValue * cv->r_e_l[face]) + cv->z->value;
        cv->setFaceElevation(face, edgeZ);
      }
    }
  }
}

