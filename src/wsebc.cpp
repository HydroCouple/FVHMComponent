#include "stdafx.h"
#include "wsebc.h"
#include "fvhmcomponent.h"
#include "core/valuedefinition.h"
#include "spatial/geometry.h"
#include "spatial/edge.h"
#include "spatial/polygon.h"
#include "temporal/timedata.h"
#include "spatial/point.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include <math.h>

using namespace HydroCouple;
using namespace HydroCouple::Spatial;
using namespace SDKTemporal;

using namespace std;

WSEBC::WSEBC(const QString &id,
             InletOutletType inletOutletType,
             double velRelaxationFactor,
             Dimension *timeDimension,
             Dimension *geometryDimension,
             ValueDefinition *valueDefinition,
             FVHMComponent *modelComponent)
  : TimeEdgeBC(id,timeDimension, geometryDimension,valueDefinition,modelComponent),
    m_velRelaxationFactor(velRelaxationFactor),
    m_inletOutletType(inletOutletType)
{
}

WSEBC::~WSEBC()
{

}

void WSEBC::prepare()
{
  m_missingData = valueDefinition()->missingValue().toDouble();
  m_defaultValue = valueDefinition()->defaultValue().toDouble();

  std::vector<TriCV*>& controlVolumes = m_modelComponent->m_controlVolumes;

  for(unordered_map<HCGeometry*,set<Edge*>>::iterator it = m_edges.begin() ; it != m_edges.end() ; it++)
  {
    //    double totalLength = 0;

    for(Edge *edge : it->second)
    {
      HCPolygon *tri = edge->leftInternal();
      TriCV *cv = controlVolumes[tri->index()];
      int index = edge->marker();

      cv->faceDepths[index].isBC = true;
      cv->hasEdgeDepthBC = true;

      FaceNormVelBC &faceVel = cv->faceNormalVels[index];
      faceVel.calculateWallShearStress = false;

      switch (m_inletOutletType)
      {
        case Inlet:
          {
            faceVel.isBC = true;
            faceVel.calculateFromQ = false;
            faceVel.velocityCalculateMode = FaceNormVelBC::ZeroGradient;
            faceVel.initialiazeVelocityVariable();
          }
          break;
        case Outlet:
          {
            faceVel.isBC = true;
            faceVel.calculateFromQ = false;
            faceVel.velocityCalculateMode = FaceNormVelBC::PreCalculated;
            faceVel.initialiazeVelocityVariable();
          }
          break;
        case OutletInlet:
          {
            faceVel.isBC = true;
            faceVel.calculateFromQ = false;
            faceVel.velocityCalculateMode = FaceNormVelBC::ZeroGradient;
            faceVel.initialiazeVelocityVariable();
          }
          break;
        case None:
          break;
      }
    }
  }

  timeCursor()->resetCursor();
}

void WSEBC::applyBoundaryConditions(double dateTime, double prevTimeStep)
{
  bool applied = false;
  std::vector<TriCV*> & controlVolumes = m_modelComponent->m_controlVolumes;

  if(m_times.size())
  {
    int index = findDateTimeIndex(dateTime);

    if(index > -1)
    {
      double dtu = m_times[index]->modifiedJulianDay();

      if(dtu == dateTime)
      {
        for(size_t i = 0; i < m_geometries.size() ; i++)
        {
          HCGeometry* geometry = m_geometries[i].data();
          set<Edge*> & edges = m_edges[geometry];

          double value = (*this)(index ,i);

          for(Edge *edge : edges)
          {
            HCPolygon *tri = edge->leftInternal();
            TriCV *cv = controlVolumes[tri->index()];
            int faceIndex = edge->marker();

            double depth = value - cv->ecz[faceIndex] ;
            depth = max(0.0, depth);

            VarBC &edgeDepth = cv->faceDepths[faceIndex];
            edgeDepth.value = depth;
            edgeDepth.associatedValue = max(value , cv->snz[0]);

            FaceNormVelBC &faceVel = cv->faceNormalVels[faceIndex];

            faceVel.vel->v[0] = 0.0;
            faceVel.vel->v[1] = 0.0;
            faceVel.value = 0.0;
            faceVel.associatedValue = 0.0;

            cv->calculateWSEGradient();

            if(depth > 1e-4)
            {
              double r_eta = cv->r_eta[faceIndex];
              double factor = edgeDepth.value * r_eta;

              faceVel.value = 0;

              switch (m_inletOutletType)
              {
                case Outlet:
                  {
                    std::tuple<double,double, double> fv = FVHMComponent::calculateZeroGradientFaceVelocity(cv,faceIndex);

                    faceVel.value = max(std::get<2>(fv),0.0);
                    faceVel.vel->v[0] = std::get<0>(fv);
                    faceVel.vel->v[1] = std::get<1>(fv);
                    faceVel.associatedValue = faceVel.value * factor;

                    if(m_modelComponent->m_printFrequencyCounter >= m_modelComponent->m_printFrequency)
                      printf("FVHM QOut: %f, VOut: %f, EdgeDepth: %f\n",faceVel.associatedValue,faceVel.value,edgeDepth.value);
                  }
                  break;
                case Inlet:
                  {
                    std::tuple<double,double, double> fv = FVHMComponent::calculateZeroGradientFaceVelocity(cv,faceIndex);

                    faceVel.value = min(std::get<2>(fv),0.0);
                    faceVel.vel->v[0] = std::get<0>(fv);
                    faceVel.vel->v[1] = std::get<1>(fv);
                    faceVel.associatedValue = faceVel.value * factor;

                    if(m_modelComponent->m_printFrequencyCounter >= m_modelComponent->m_printFrequency)
                      printf("FVHM Qin: %f, Vin: %f, EdgeDepth: %f\n",faceVel.associatedValue,faceVel.value,edgeDepth.value);
                  }
                  break;
                case OutletInlet:
                  {
                    std::tuple<double,double, double> fv = FVHMComponent::calculateZeroGradientFaceVelocity(cv,faceIndex);

                    faceVel.value = std::get<2>(fv);
                    faceVel.vel->v[0] = std::get<0>(fv);
                    faceVel.vel->v[1] = std::get<1>(fv);
                    faceVel.associatedValue = faceVel.value * factor;

                    if(m_modelComponent->m_printFrequencyCounter >= m_modelComponent->m_printFrequency)
                      printf("FVHM Qin: %f, Vin: %f, EdgeDepth: %f\n",faceVel.associatedValue,faceVel.value,edgeDepth.value);
                  }
                  break;
                case None:
                  {
                    faceVel.vel->v[0] = 0.0;
                    faceVel.vel->v[1] = 0.0;
                    faceVel.value = 0.0;
                    faceVel.associatedValue = 0.0;
                  }
                  break;
              }
            }
          }
        }
      }
      else if(index - 1 > -1)
      {
        double dtl = m_times[index - 1]->modifiedJulianDay();
        double factor  = (dateTime - dtl) / (dtu - dtl);

        for(size_t i = 0; i < m_geometries.size() ; i++)
        {
          HCGeometry *geometry = m_geometries[i].data();
          set<Edge*> & edges = m_edges[geometry];

          double valueu = (*this)(index,i);
          double valuel = (*this)(index - 1,i);

          double value = valuel + factor * (valueu - valuel);

          for(Edge *edge : edges)
          {
            HCPolygon *tri = edge->leftInternal();
            TriCV *cv = controlVolumes[tri->index()];
            int faceIndex = edge->marker();

            double depth = value - cv->ecz[faceIndex] ;
            depth = max(0.0, depth);

            VarBC &edgeDepth = cv->faceDepths[faceIndex];
            edgeDepth.value = depth;
            edgeDepth.associatedValue = max(value , cv->snz[0]);

            cv->calculateWSEGradient();

            double factor = edgeDepth.value * cv->r_eta[faceIndex];

            FaceNormVelBC &faceVel = cv->faceNormalVels[faceIndex];
            faceVel.vel->v[0] = 0.0;
            faceVel.vel->v[1] = 0.0;
            faceVel.value = 0.0;
            faceVel.associatedValue = 0.0;


            if(edgeDepth.value > 1e-4)
            {
              switch (m_inletOutletType)
              {
                case Outlet:
                  {
                    std::tuple<double,double, double> fv = FVHMComponent::calculateZeroGradientFaceVelocity(cv,faceIndex);

                    faceVel.value = max(std::get<2>(fv),0.0);
                    faceVel.vel->v[0] = std::get<0>(fv);
                    faceVel.vel->v[1] = std::get<1>(fv);
                    faceVel.associatedValue = faceVel.value * factor;

                    if(m_modelComponent->m_printFrequencyCounter >= m_modelComponent->m_printFrequency)
                      printf("FVHM QOut: %f, VOut: %f, EdgeDepth: %f\n",faceVel.associatedValue,faceVel.value,edgeDepth.value);
                  }
                  break;
                case Inlet:
                  {
                    std::tuple<double,double, double> fv = FVHMComponent::calculateZeroGradientFaceVelocity(cv,faceIndex);

                    faceVel.value = min(std::get<2>(fv),0.0);
                    faceVel.vel->v[0] = std::get<0>(fv);
                    faceVel.vel->v[1] = std::get<1>(fv);
                    faceVel.associatedValue = faceVel.value * factor;

                    if(m_modelComponent->m_printFrequencyCounter >= m_modelComponent->m_printFrequency)
                      printf("FVHM Qin: %f, Vin: %f, EdgeDepth: %f\n",faceVel.associatedValue,faceVel.value,edgeDepth.value);
                  }
                  break;
                case OutletInlet:
                  {
                    std::tuple<double,double, double> fv = FVHMComponent::calculateZeroGradientFaceVelocity(cv,faceIndex);

                    faceVel.value = std::get<2>(fv);
                    faceVel.vel->v[0] = std::get<0>(fv);
                    faceVel.vel->v[1] = std::get<1>(fv);
                    faceVel.associatedValue = faceVel.value * factor;

                    if(m_modelComponent->m_printFrequencyCounter >= m_modelComponent->m_printFrequency)
                      printf("FVHM Qin: %f, Vin: %f, EdgeDepth: %f\n",faceVel.associatedValue,faceVel.value,edgeDepth.value);
                  }
                  break;
                case None:
                  {
                    faceVel.vel->v[0] = 0.0;
                    faceVel.vel->v[1] = 0.0;
                    faceVel.value = 0.0;
                    faceVel.associatedValue = 0.0;
                  }
                  break;
              }
            }
          }
        }
      }

      applied = true;

    }
  }
  else
  {
    for(size_t i = 0; i < m_geometries.size() ; i++)
    {
      HCGeometry *geometry = m_geometries[i].data();
      set<Edge*> & edges = m_edges[geometry];

      for(Edge *edge : edges)
      {
        HCPolygon *tri = edge->leftInternal();
        TriCV *cv = controlVolumes[tri->index()];
        int faceIndex = edge->marker();

        FaceNormVelBC &faceVel = cv->faceNormalVels[faceIndex];
        faceVel.value = 0.0;
      }
    }
  }
}

double WSEBC::velocityRelaxationFactor() const
{
  return m_velRelaxationFactor;
}

void WSEBC::setVelocityRelaxationFactor(double velocityRalaxationFactor)
{
  m_velRelaxationFactor = velocityRalaxationFactor;
}

