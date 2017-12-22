#include "stdafx.h"
#include "inletoutletflowbc.h"
#include "fvhmcomponent.h"

#include "core/valuedefinition.h"

#include "spatial/geometry.h"
#include "spatial/edge.h"
#include "spatial/polygon.h"

#include "controlvolume.h"
#include "temporal/timedata.h"

#include <math.h>
#include <unordered_map>

using namespace HydroCouple;
using namespace HydroCouple::Spatial;
using namespace SDKTemporal;
using namespace std;

InletOutletFlowBCArgument::InletOutletFlowBCArgument(const QString &id,
                                                     Dimension *timeDimension,
                                                     Dimension *geometryDimension,
                                                     ValueDefinition *valueDefinition,
                                                     FVHMComponent *modelComponent)
  :TimeEdgeBC(id,timeDimension, geometryDimension, valueDefinition,modelComponent)
{
}

InletOutletFlowBCArgument::~InletOutletFlowBCArgument()
{
}

void InletOutletFlowBCArgument::prepare()
{
  std::vector<TriCV*>& controlVolumes = m_modelComponent->m_controlVolumes;

  m_missingData = valueDefinition()->missingValue().toDouble();
  m_defaultValue = valueDefinition()->defaultValue().toDouble();

  for(unordered_map<HCGeometry*, set<Edge*>>::iterator it = m_edges.begin() ; it != m_edges.end() ; it++)
  {
    set<Edge*> edgesForG = (*it).second;

    for(Edge *edge :edgesForG)
    {
      HCPolygon *tri = edge->leftInternal();
      TriCV *cv = controlVolumes[tri->index()];
      int faceIndex = edge->marker();

      FaceNormVelBC &faceVel = cv->faceNormalVels[faceIndex];
      faceVel.calculateWallShearStress = false;
      faceVel.isBC = true;
      faceVel.calculateFromQ = true;
      faceVel.velocityCalculateMode = FaceNormVelBC::PreCalculated;
      faceVel.initialiazeVelocityVariable();
    }
  }

  timeCursor()->resetCursor();

}


InletFlowBCArgument::InletFlowBCArgument(const QString &id,
                                         Dimension *timeDimension,
                                         Dimension *geometryDimension,
                                         ValueDefinition *valueDefinition,
                                         FVHMComponent *modelComponent)
  :InletOutletFlowBCArgument(id,timeDimension, geometryDimension,valueDefinition,modelComponent)
{
}

InletFlowBCArgument::~InletFlowBCArgument()
{

}

void InletFlowBCArgument::prepare()
{
  const std::vector<TriCV*>& controlVolumes = m_modelComponent->m_controlVolumes;

  m_missingData = valueDefinition()->missingValue().toDouble();
  m_defaultValue = valueDefinition()->defaultValue().toDouble();

  for(unordered_map<HCGeometry*, set<Edge*>>::iterator it = m_edges.begin() ; it != m_edges.end() ; it++)
  {
    set<Edge*> edgesForG = (*it).second;

    for(Edge *edge :edgesForG)
    {
      HCPolygon *tri = edge->leftInternal();
      TriCV *cv = controlVolumes[tri->index()];
      int faceIndex = edge->marker();

      FaceNormVelBC &faceVel = cv->faceNormalVels[faceIndex];
      faceVel.isBC = true;
      faceVel.calculateFromQ = true;
      faceVel.velocityCalculateMode = FaceNormVelBC::DotProductEdgeVel;
      faceVel.calculateWallShearStress = false;
      faceVel.initialiazeVelocityVariable();
    }
  }

  timeCursor()->resetCursor();
}

void InletFlowBCArgument::applyBoundaryConditions(double dateTime, double prevTimeStep)
{
  if(m_geometries.size() && m_edges.size())
  {
    const std::vector<TriCV*> & controlVolumes = m_modelComponent->controlVolumes();

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
            double count = (double)edges.size();

            if(count > 0.0)
            {
              value = min(- value / count , 0.0);

              for(Edge *edge : edges)
              {
                HCPolygon *tri = edge->leftInternal();
                TriCV *cv = controlVolumes[tri->index()];

                int faceIndex = edge->marker();

                VarBC &edgeDepth = cv->faceDepths[faceIndex];

                double denom = edgeDepth.value * cv->r_eta[faceIndex];
                FaceNormVelBC &faceVel = cv->faceNormalVels[faceIndex];

                double vel = value/denom;

                if(cv->wetIndex && fabs(denom - 0.0) > std::numeric_limits<double>::epsilon() && fabs(vel) <= 4.0)
                {
                  faceVel.value =  vel;
                  faceVel.vel->v[0] = vel * cv->e_n[faceIndex].v[0];
                  faceVel.vel->v[1] = vel * cv->e_n[faceIndex].v[1];
                  faceVel.associatedValue = value;

                  if(m_modelComponent->m_printFrequencyCounter >= m_modelComponent->m_printFrequency)
                    printf("FVHM Q: %f\tV: %f\tVx: %f\tVy: %f\tEdgeDepth: %f\n",value,vel,faceVel.vel->v[0], faceVel.vel->v[1], cv->faceDepths[faceIndex].value);
                }
                else
                {
                  faceVel.vel->v[0] = 0.0;
                  faceVel.vel->v[1] = 0.0;
                  faceVel.value = 0.0;
                  faceVel.associatedValue = 0.0;
                  cv->inflowOutflow -= value;
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

            double count = (double)edges.size();

            if(count > 0)
            {
              value = min(-1.0 * value / count , 0.0);

              for(Edge *edge : edges)
              {
                HCPolygon *tri = edge->leftInternal();
                TriCV *cv = controlVolumes[tri->index()];
                int faceIndex = edge->marker();

                VarBC &edgeDepth = cv->faceDepths[faceIndex];
                double denom = edgeDepth.value * cv->r_eta[faceIndex];
                FaceNormVelBC &faceVel = cv->faceNormalVels[faceIndex];

                double vel = value / denom;

                if(cv->wetIndex && fabs(denom - 0.0) > std::numeric_limits<double>::epsilon() && fabs(vel) <= 4.0)
                {
                  faceVel.value =  vel * 1.00000;
                  faceVel.vel->v[0] = vel * cv->e_n[faceIndex].v[0] * 1.00000;
                  faceVel.vel->v[1] = vel * cv->e_n[faceIndex].v[1] * 1.00000;
                  faceVel.associatedValue = value;

                  if(m_modelComponent->m_printFrequencyCounter >= m_modelComponent->m_printFrequency)
                    printf("FVHM Q: %f\tV: %f\tVx: %f\tVy: %f\tEdgeDepth: %f\n",value,vel,faceVel.vel->v[0], faceVel.vel->v[1], cv->faceDepths[faceIndex].value);
                }
                else
                {
                  vel = 0.0;
                  faceVel.vel->v[0] = 0.0;
                  faceVel.vel->v[1] = 0.0;
                  faceVel.value = 0.0;
                  faceVel.associatedValue = 0.0;
                  cv->inflowOutflow -= value;

                  if(m_modelComponent->m_printFrequencyCounter >= m_modelComponent->m_printFrequency)
                    printf("FVHM Q: %f\tV: %f\tVx: %f\tVy: %f\tEdgeDepth: %f\n",value,vel,faceVel.vel->v[0], faceVel.vel->v[1], cv->faceDepths[faceIndex].value);
                }
              }
            }
          }
        }
      }
      else
      {
        for(size_t i = 0; i < m_geometries.size() ; i++)
        {
          HCGeometry* geometry = m_geometries[i].data();
          set<Edge*> & edges = m_edges[geometry];

          for(Edge *edge : edges)
          {
            HCPolygon *tri = edge->leftInternal();
            TriCV *cv = controlVolumes[tri->index()];
            int faceIndex = edge->marker();

            VarBC &edgeDepth = cv->faceDepths[faceIndex];
            double denom = edgeDepth.value * cv->r_eta[faceIndex];

            FaceNormVelBC &faceVel = cv->faceNormalVels[faceIndex];

            if(cv->wetIndex && edgeDepth.value > 1e-3)
            {
              double vel = min(m_defaultValue / denom , 0.0);

              faceVel.value = vel * 1.00000;
              faceVel.vel->v[0] = vel * cv->e_n[faceIndex].v[0] * 1.00000;
              faceVel.vel->v[1] = vel * cv->e_n[faceIndex].v[1] * 1.00000;
              faceVel.associatedValue = min(m_defaultValue,0.0);
            }
            else
            {
              faceVel.vel->v[0] = 0.0;
              faceVel.vel->v[1] = 0.0;
              faceVel.value = 0.0;
              faceVel.associatedValue = 0.0;
              cv->inflowOutflow += m_defaultValue;
            }
          }
        }
      }
    }
    else
    {
      for(size_t i = 0; i < m_geometries.size() ; i++)
      {
        HCGeometry* geometry = m_geometries[i].data();
        set<Edge*> & edges = m_edges[geometry];

        for(Edge *edge : edges)
        {
          HCPolygon *tri = edge->leftInternal();
          TriCV *cv = controlVolumes[tri->index()];
          int faceIndex = edge->marker();

          VarBC &edgeDepth = cv->faceDepths[faceIndex];
          double denom = edgeDepth.value * cv->r_eta[faceIndex];

          FaceNormVelBC &faceVel = cv->faceNormalVels[faceIndex];

          if(cv->wetIndex && edgeDepth.value > 1e-3)
          {
            double vel = min(m_defaultValue / denom , 0.0);

            faceVel.value = vel * 1.00000;
            faceVel.vel->v[0] = vel * cv->e_n[faceIndex].v[0] * 1.00000;
            faceVel.vel->v[1] = vel * cv->e_n[faceIndex].v[1] * 1.00000;
            faceVel.associatedValue = min(m_defaultValue,0.0);
          }
          else
          {
            faceVel.vel->v[0] = 0.0;
            faceVel.vel->v[1] = 0.0;
            faceVel.value = 0.0;
            faceVel.associatedValue = 0.0;
            cv->inflowOutflow += m_defaultValue;
          }
        }
      }
    }
  }
}


OutletFlowBCArgument::OutletFlowBCArgument(const QString &id,
                                           Dimension *timeDimension,
                                           Dimension *geometryDimension,
                                           ValueDefinition *valueDefinition,
                                           FVHMComponent *modelComponent)
  :InletOutletFlowBCArgument(id,timeDimension, geometryDimension,valueDefinition,modelComponent)
{

}

OutletFlowBCArgument::~OutletFlowBCArgument()
{

}

void OutletFlowBCArgument::prepare()
{
  InletOutletFlowBCArgument::prepare();
}

void OutletFlowBCArgument::applyBoundaryConditions(double dateTime, double prevTimeStep)
{

  if(m_geometries.size() && m_edges.size())
  {
    const std::vector<TriCV*> & controlVolumes = m_modelComponent->controlVolumes();

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
            double count = (double)edges.size();

            if(count > 0.0)
            {
              value = max(value / count , 0.0);

              for(Edge *edge : edges)
              {
                HCPolygon *tri = edge->leftInternal();
                TriCV *cv = controlVolumes[tri->index()];

                int faceIndex = edge->marker();
                double denom = cv->faceDepths[faceIndex].value * cv->r_eta[faceIndex];

                if(cv->faceDepths[faceIndex].value >= m_modelComponent->m_wetCellDepth)
                {
                  double vel = value / denom;
                  cv->faceNormalVels[faceIndex].value =  vel;
                  cv->faceNormalVels[faceIndex].associatedValue = value;
                }
                else
                {
                  cv->faceNormalVels[faceIndex].value = 0.0;
                  cv->faceNormalVels[faceIndex].associatedValue = 0.0;
                  cv->inflowOutflow += value;
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

            double count = (double)edges.size();


            if(count > 0)
            {
              value = max(value / count , 0.0);


              for(Edge *edge : edges)
              {
                HCPolygon *tri = edge->leftInternal();
                TriCV *cv = controlVolumes[tri->index()];
                int faceIndex = edge->marker();

                double denom = cv->faceDepths[faceIndex].value * cv->r_eta[faceIndex];

                if(cv->faceDepths[faceIndex].value >= m_modelComponent->m_wetCellDepth)
                {
                  double vel = value / denom;
                  cv->faceNormalVels[faceIndex].value = vel;
                  cv->faceNormalVels[faceIndex].associatedValue = value;
                }
                else
                {
                  cv->faceNormalVels[faceIndex].value = 0.0;
                  cv->faceNormalVels[faceIndex].associatedValue = 0.0;
                  cv->inflowOutflow += value;
                }
              }
            }
          }
        }
      }
      else
      {
        for(size_t i = 0; i < m_geometries.size() ; i++)
        {
          HCGeometry* geometry = m_geometries[i].data();
          set<Edge*> & edges = m_edges[geometry];

          for(Edge *edge : edges)
          {
            HCPolygon *tri = edge->leftInternal();
            TriCV *cv = controlVolumes[tri->index()];
            int faceIndex = edge->marker();

            double denom = cv->faceDepths[faceIndex].value * cv->r_eta[faceIndex];

            if(cv->faceDepths[faceIndex].value >= m_modelComponent->m_wetCellDepth)
            {
              double vel = max(m_defaultValue / denom , 0.0);
              cv->faceNormalVels[faceIndex].value = vel;
              cv->faceNormalVels[faceIndex].associatedValue = max(m_defaultValue,0.0);
            }
            else
            {
              cv->faceNormalVels[faceIndex].value = 0.0;
              cv->faceNormalVels[faceIndex].associatedValue = 0.0;
              cv->inflowOutflow += m_defaultValue;
            }
          }
        }
      }
    }
    else
    {
      for(size_t i = 0; i < m_geometries.size() ; i++)
      {
        HCGeometry* geometry = m_geometries[i].data();
        set<Edge*> & edges = m_edges[geometry];

        for(Edge *edge : edges)
        {
          HCPolygon *tri = edge->leftInternal();
          TriCV *cv = controlVolumes[tri->index()];
          int faceIndex = edge->marker();

          double denom = cv->faceDepths[faceIndex].value * cv->r_eta[faceIndex];

          if(cv->faceDepths[faceIndex].value >= m_modelComponent->m_wetCellDepth)
          {
            double vel = max(m_defaultValue / denom , 0.0);
            cv->faceNormalVels[faceIndex].value = vel;
            cv->faceNormalVels[faceIndex].associatedValue = max(m_defaultValue, 0.0);
          }
          else
          {
            cv->faceNormalVels[faceIndex].value = 0.0;
            cv->faceNormalVels[faceIndex].associatedValue = 0.0;
            cv->inflowOutflow += m_defaultValue;
          }
        }
      }
    }
  }
}

