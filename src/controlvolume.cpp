#include "stdafx.h"
#include "controlvolume.h"
#include "fvhmcomponent.h"
#include "spatial/polygon.h"
#include "spatial/edge.h"
#include "spatial/point.h"

#include "qrsolve.h"

#include <math.h>

#ifdef USE_OPENMP
#include <omp.h>
#endif

using namespace std;

int TriCV::numEdges = 3;
int TriCV::gradientCalcMode = 1;

TriCV::TriCV(HCTriangle *cell, FVHMComponent *modelComponent)
  : modelComponent(modelComponent)
{
  this->cell = cell;
  index = cell->index();
  area = cell->area3D();

  int numNodes = numEdges + 1;

  minNodeDepth = 0.0;

  numNodeTriangles = new int[numEdges];
  nodeTriangleIndexes = new int*[numEdges];
  nodeTriangleDistances = new double*[numEdges];
  //  nodeTriangleDistancesV = new Vect*[numEdges];

  vertices = new HCVertex*[numNodes];

  snz = new double[numEdges];
  snzMinEdge = new double[numEdges];
  snzMaxEdge = new double[numEdges];

  ecz = new double[numEdges];
  nz = new double[numNodes];
  nWSE = new double[numNodes];

  orderedNTris = new int[numNodes];
  orderedNTrisIndexes = new int[numNodes];

  r_e = new Vect[numEdges];
  r_e_cvn = new Vect[numEdges];
  r_e_l = new double[numEdges];
  r_e_l_p = new double[numEdges];
  r_n = new Vect[numNodes];
  e_n = new Vect[numNodes];
  df = new Vect[numEdges];
  r_eta = new double[numEdges];
  nTris = new int[numEdges];

  cz = 0.0;

  prevZ = new VarBC();
  z = new VarBC();
  prevH = new VarBC();
  h = new VarBC();

  prevH->value = modelComponent->m_viscousSubLayerDepth;
  h->value = modelComponent->m_viscousSubLayerDepth;

  vel = new VarBC[2];
  prevVel = new VarBC[2];

  prevIterVel = new double[2]();

  r_xi = new Vect[numEdges];
  r_xi_l = new double[numEdges]();
  r_xi_dot_e_n = new double[numEdges]();
  friction = new double[2]();

  wallShearFriction = new double*[2];
  velCoeffs = new double*[2];
  nodeVels = new double*[2];

  for(int i = 0; i < 2; i++)
  {
    wallShearFriction[i] = new double[numEdges]();
    velCoeffs[i] = new double[numNodes]();
    nodeVels[i] = new double[numNodes]();
  }

  nodalEddyViscosity = new double[numEdges]();

  faceDepths = new VarBC[numEdges];

  grad_z = new VectBC();
  grad_h = new Vect();
  grad_zcorr = new Vect();

  grad_vel = new VectBC[2];
  faceNormalVels = new FaceNormVelBC[numEdges];

  dir_diff_term = new double[numEdges]();
  cross_diff_term = new double[numEdges]();

  externalForce = new Vect();

  w_l = new double[numEdges]();
  e_xi = new Vect[numEdges];

  maxVel = new double[2]{std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest()};

  velResidual = new double[2]();
  velResidualIter = new double[2]();
  prevIterVel = new double[2]();

  Edge *edge = cell->edgeInternal();

  std::vector<double> tz(3);

  double cx = 0.0, cy = 0.0;
  cz = 0.0;

  for(int i = 0; i < numEdges ; i++)
  {
    HCVertex* p1 = edge->origInternal();
    HCVertex* p2 = edge->destInternal();

    edge->setMarker(i);
    edge = edge->leftNextInternal();

    HCVertex* p3 = edge->destInternal();

    vertices[i] =  p1;
    vertices[i+1] = p2;

    nz[i] =  nWSE[i] = p1->z();
    nz[i+1] = nWSE[i+1] = p2->z();

    cz += p1->z();
    cx += p1->x();
    cy += p1->y();

    tz[i] = p1->z();

    ecz[i] = (p1->z() + p2->z()) / 2.0;
    snzMinEdge[i] = min(p1->z(), p2->z());
    snzMaxEdge[i] = max(p1->z(), p2->z());

    Vect currEdge = Vect(*p2) - Vect(*p1);
    Vect nextEdge = Vect(*p3) - Vect(*p2);

    currEdge.v[2] = 0.0;
    r_eta[i] = currEdge.length();

    e_n[i] = currEdge.unitNormal2dToVector();

    nextEdge.v[2] = 0.0;
    e_n[i+1] = nextEdge.unitNormal2dToVector();

    nodeTriangleIndexes[i] = nullptr;
    nodeTriangleDistances[i] = nullptr;
    nTris[i] = -1;
  }

  std::sort(tz.begin(),tz.end());

  cx = cx / numEdges;
  cy = cy / numEdges;

  center = new Vect(cx,cy,0.0);
  cz = cz / numEdges;

  for(int i = 0 ; i < numEdges + 1; i++)
  {
    HCPoint *p1 = vertices[i];
    HCPoint *p2 = vertices[i+1];

    r_n[i] = Vect(p1->x(), p1->y()) -  Vect(*center);

    if(i < numEdges)
    {
      r_e[i] = Vect((p1->x() + p2->x())/2.0, (p1->y() + p2->y())/2.0, 0.0) - (*center);
      r_e_l[i] =  r_e[i].length();

      //      double fac = Vect::dotProduct(e_n[i], r_e[i]);
      //      r_e_l_p[i].v[0] = fabs(fac * e_n[i].v[0]);
      //      r_e_l_p[i].v[1] = fabs(fac * e_n[i].v[1]);
      r_e_l_p[i] = fabs(Vect::dotProduct(e_n[i], r_e[i]));
      snz[i] = tz[i];
    }

    nWSE[i] = tz[0];
    orderedNTris[i] = -1;
    orderedNTrisIndexes[i] = -1;
  }

  z->value = prevZ->value = snz[0];
  h->value = prevH->value = modelComponent->m_viscousSubLayerDepth;

  calculateInitialWSEGradient();
}

TriCV::~TriCV()
{

  delete[] nTris;
  delete[] orderedNTris;
  delete[] orderedNTrisIndexes;
  delete[] numNodeTriangles;

  for(int i = 0 ; i < numEdges; i++)
  {
    if(nodeTriangleIndexes[i])
      delete[]  nodeTriangleIndexes[i];

    if(nodeTriangleDistances[i])
    {
      delete[] nodeTriangleDistances[i];
      //      delete[] nodeTriangleDistancesV[i];
    }
  }

  delete[] nodeTriangleIndexes;
  delete[] nodeTriangleDistances;
  delete[] vertices;

  delete z;
  delete grad_z;

  if(grad_z_init)
    delete grad_z_init;

  delete grad_h;
  delete grad_zcorr;
  delete prevZ;

  delete h;
  delete prevH;

  delete[] vel;
  delete[] grad_vel;
  delete[] prevVel;
  delete[] maxVel;

  delete[] e_n;
  delete[] e_xi;
  delete[] r_e;
  delete[] r_e_cvn;
  delete[] r_e_l;
  delete[] r_e_l_p;
  delete[] r_n;
  delete[] df;
  delete[] r_eta;
  delete[] r_xi;
  delete[] r_xi_dot_e_n;
  delete[] r_xi_l;
  delete[] w_l;

  delete[] faceNormalVels;

  for(int i = 0; i < 2; i++)
  {
    delete[] velCoeffs[i];
    delete[] wallShearFriction[i];
    delete[] nodeVels[i];
  }

  delete[] nodeVels;

  delete center;

  delete[] snz;
  delete[] snzMinEdge;
  delete[] snzMaxEdge;

  delete[] nz;
  delete[] ecz;
  delete[] nWSE;
  //  delete[] nZCorr;

  delete[] faceDepths;

  delete[] dir_diff_term;
  delete[] cross_diff_term;
  //  delete[] r_p_corr;
  //  delete[] r_n_corr;

  delete externalForce;

  delete[] nodalEddyViscosity;
  delete[] velResidual;
  delete[] velResidualIter;

  delete[] velCoeffs;
  delete[] prevIterVel;
  delete[] friction;
  delete[] wallShearFriction;

}

void TriCV::calculateAdjacentCellParams()
{
  std::vector<Vect*> points(numEdges+1);
  Edge *edge = cell->edgeInternal();

  std::map<int,int> orderedTris;
  orderedTris[cell->index()] = 0;
  totalnTris = 0;

  for(int i = 0; i < numEdges + 1 ; i++)
  {
    HCVertex *vertex = edge->origInternal();
    Vect *p = new Vect(*vertex);
    p->v[2] = 0.0;

    if(i < numEdges)
    {
      if(edge->rightInternal())
      {
        int nIndex = edge->rightInternal()->index();
        nTris[i] = nIndex;
        orderedTris[nIndex] = i + 1;
        faceNormalVels[i].calculateWallShearStress = false;
      }
      else
      {
        faceNormalVels[i].calculateWallShearStress = true;
      }

      std::vector<int> indexes;
      std::vector<double> distances;
      Edge *vedge = edge;

      do
      {
        if(vedge->leftInternal())
        {
          int index = vedge->leftInternal()->index();
          indexes.push_back(index);
          TriCV *cvn = modelComponent->m_controlVolumes[index];
          Vect dist = Vect(cvn->center->x() , cvn->center->y()) - (*p);
          distances.push_back(dist.length());
        }

        vedge = vedge->origNextInternal();

      } while(vedge != edge);

      numNodeTriangles[i] = indexes.size();
      nodeTriangleIndexes[i] = new int[indexes.size()];
      nodeTriangleDistances[i] = new double[indexes.size()];

      totalnTris += indexes.size();

      for(size_t d = 0 ; d < indexes.size(); d++)
      {
        nodeTriangleIndexes[i][d] = indexes[d];
        nodeTriangleDistances[i][d] = distances[d];
      }
    }

    points[i] = p;
    edge = edge->leftNextInternal();
  }

  int i = 0;

  for(const auto & orderedTri : orderedTris)
  {
    orderedNTris[i] = orderedTri.first;
    orderedNTrisIndexes[i] = orderedTri.second;
    i++;
  }

  Vect *e_eta = new Vect[numEdges];
  Vect *eta = new Vect[numEdges];
  Vect *edgesCenter = new Vect[numEdges];

  for(int i = 0; i < numEdges; i++)
  {
    eta[i] = Vect(*points[i+1]) - Vect(*points[i]);
    e_eta[i] = eta[i].normalized();
    edgesCenter[i] = ((*points[i]) + (*points[i+1])) / 2.0;
  }

  for(int i = 0 ; i < numEdges; i++)
  {
    int index = nTris[i];

    if(index > -1)
    {
      TriCV * cvn = modelComponent->m_controlVolumes[index];
      Vect v = (*cvn->center) - (*center);
      r_xi[i] = v;

      double lpn = v.length();

      r_xi_l[i] = lpn;

      r_xi_dot_e_n[i] = Vect::dotProduct(r_xi[i], e_n[i]);
      e_xi[i] = v.normalized();

      Vect &edgeC = edgesCenter[i];

      Vect f;
      intersectPoint((*points[i]), (*points[i+1]), *center, *cvn->center, f);

      df[i] = edgeC - f;
      Vect r2 =  edgeC - (*cvn->center);

      r_e_cvn[i] = r2;

      //      r_p_corr[i] = (edgeC - Vect::dotProduct(edgeC - (*center) , e_n[i]) * e_n[i]) - (*center);
      //      r_n_corr[i] = (edgeC - Vect::dotProduct(edgeC - (*cvn->center) , e_n[i]) * e_n[i]) - (*cvn->center);

      //      double sig1 = r_e_l[i];
      //      double sg1 = Vect::dotProduct(r_e[i], e_n[i]);

      double sig2 = r2.length();
      //      double sg2 = Vect::dotProduct(-1.0 * r2, e_n[i]);

      //w_l[i] = (1.0 / (sig2 * sig2)) / ((1.0 / (sig2 * sig2)) + (1.0 / (sig1 * sig1)));
      //      w_l[i] = sig1 / (sig1 + r2.length());
      w_l[i] = 1.0 - sig2 / lpn;
      //cx[i] = Vect::dotProduct(-sg1 * r2 - sg2 * r_e[i], eta[i]) / ((sg1 + sg2) * r_eta[i] * r_eta[i]);

      //      r_e_n[i] = r2;

      //turbulence terms
      dir_diff_term[i] = Vect::dotProduct(e_n[i], e_n[i] * r_eta[i]) / (Vect::dotProduct(e_n[i],e_xi[i]) * r_xi_l[i]);
      cross_diff_term[i] = Vect::dotProduct(e_xi[i], e_eta[i] * r_eta[i])  / (Vect::dotProduct(e_n[i],e_xi[i]) * r_eta[i]);
    }
    else
    {
      w_l[i] = 0.0;
    }
  }

  delete[] e_eta;
  delete[] eta;
  delete[] edgesCenter;

  for(Vect *v : points)
    delete v;

  points.clear();
}

double TriCV::getCourantFactor() const
{
  double factor = 0;

  for(int i = 0; i < 3 ; i++)
  {
    if(faceNormalVels[i].value > 0.0)
    {
      factor +=  faceNormalVels[i].value * r_eta[i] * faceDepths[i].value;
    }
  }

  if(inflowOutflow < 0)
  {
    factor -= inflowOutflow / area;
  }

  factor = max(0.0, factor);

  double actualH = max(0.0, h->value - modelComponent->m_wetCellDepth);

  if(actualH)
  {
    factor = factor / (area * actualH);
  }
  else
  {
    factor = 0.0;
  }

  return factor;
}

void TriCV::setVFRWSE(double wse, bool prev)
{
  double zz1 = snz[0];
  double zz2 = snz[1];
  double zz3 = snz[2];
  wse = max(wse, zz1);

  if(prev)
  {
    if(wse > zz3)
    {
      prevH->value = wse - cz;
      prevZ->value = wse;
    }
    else if(wse > zz2 && wse <= zz3)
    {
      double dp = (wse * wse + wse * zz3 - 3.0 * wse * zz1 - zz3 * zz2 + zz1 * zz2 + zz1 * zz1) /(3.0 * (zz3 - zz1));
      prevH->value = dp;
      prevZ->value = wse;
    }
    else  if(wse > zz1 && wse <= zz2)
    {
      double num = wse - zz1;
      double dp = (num * num * num) / (3.0 * (zz2 - zz1) * (zz3 - zz1));
      prevH->value = dp;
      prevZ->value = wse;
    }
    else
    {
      prevZ->value = zz1;
      prevH->value = modelComponent->m_viscousSubLayerDepth;
    }
  }
  else
  {
    if(wse > zz3)
    {
      h->value = wse - cz;
      z->value = wse;
    }
    else if(wse > zz2 && wse <= zz3)
    {
      double dp = (wse * wse + wse * zz3 - 3.0 * wse * zz1 - zz3 * zz2 + zz1 * zz2 + zz1 * zz1)/(3.0 * (zz3 - zz1));
      h->value = dp;
      z->value = wse;
    }
    else if(wse > zz1 && wse <= zz2)
    {
      double num = wse - zz1;
      double dp = (num * num * num) / (3.0 * (zz2 - zz1) * (zz3 - zz1));
      h->value = dp;
      z->value = wse;
    }
    else
    {
      z->value =  zz1;
      h->value = modelComponent->m_viscousSubLayerDepth;
    }
  }
}

void TriCV::setVFRDepth(double depth, bool prev)
{

  double zz1 = snz[0];
  double zz2 = snz[1];
  double zz3 = snz[2];

  depth = max(modelComponent->m_viscousSubLayerDepth , depth);

  if(prev)
  {
    prevH->value = depth;
    double wse = prevZ->value;

    if(wse > zz3)
    {
      prevZ->value = cz + depth;
    }
    else if(wse > zz2 && wse <= zz3)
    {
      double gamma1 = zz3 - 3.0 * zz1;
      double gamma2 = 3.0 * depth * zz1 - 3.0 * depth * zz3 - zz3 * zz2 +  zz1 * zz2 + zz1 * zz1;
      prevZ->value = 0.5 * (-gamma1 + sqrt(gamma1 * gamma1 - 4.0 * gamma2));
    }
    else if(wse > zz1 && wse <= zz2)
    {
      prevZ->value = zz1 + cbrt(3.0 * depth * (zz2 - zz1) * (zz3 - zz1)) +  modelComponent->m_viscousSubLayerDepth ;
    }
    else
    {
      prevZ->value = zz1 + cbrt(3.0 * depth * (zz2 - zz1) * (zz3 - zz1)) +  modelComponent->m_viscousSubLayerDepth  ;
    }
  }
  else
  {
    double wse = z->value;
    h->value = depth;

    if(wse > zz3)
    {
      z->value = cz + depth;
    }
    else if(wse > zz2 && wse <= zz3)
    {
      double gamma1 = zz3 - 3.0 * zz1;
      double gamma2 = 3.0 * depth * zz1 - 3.0 * depth * zz3 - zz3*zz2 +  zz1*zz2 + zz1*zz1;
      z->value = 0.5 *(-gamma1 + sqrt(gamma1*gamma1 - 4.0 * gamma2));
    }
    else if(wse > zz1 && wse <= zz2)
    {
      z->value = zz1 + cbrt(3.0 * depth * (zz2 - zz1) * (zz3 - zz1)) + modelComponent->m_viscousSubLayerDepth;
    }
    else
    {
      z->value = zz1 + cbrt(3.0 * depth * (zz2 - zz1) * (zz3 - zz1)) + modelComponent->m_viscousSubLayerDepth;
    }
  }
}

void TriCV::copyVariablesToPrev()
{
  prevZ->copy(*z);
  prevH->copy(*h);

  for(int i = 0; i < 2; i++)
  {
    prevVel[i].copy(vel[i]);
    prevIterVel[i] = vel[i].value;
    maxVel[i] = max(maxVel[i],vel[i].value);
  }

  maxZ = max(maxZ,z->value);
  maxH = max(maxH,h->value);

}

void TriCV::calculateEdgeDepths()
{
  minNodeDepth = 1e40;

  for(int i = 0; i < numEdges; i++)
  {
    VarBC &faceDepth = faceDepths[i];

    if(!faceDepth.isBC)
    {
      int cvnIndex = nTris[i];
      double edgeDepth = 0.0;
      double edgeZ = 0.0;


      if(cvnIndex > -1)
      {
        TriCV *cvn = modelComponent->m_controlVolumes[cvnIndex];
        double wl_n = w_l[i];
        double wl_p = 1.0 - wl_n;

        //        double gradh_x = wl_p * grad_h->v[0] + wl_n * cvn->grad_h->v[0];
        //        double gradh_y = wl_p * grad_h->v[1] + wl_n * cvn->grad_h->v[1];
        //        edgeDepth = wl_p * h->value + wl_n * cvn->h->value + Vect::dotProduct(gradh_x, gradh_y, 0.0, df[i]);

        double gradz_x = wl_p * grad_z->v[0] + wl_n * cvn->grad_z->v[0];
        double gradz_y = wl_p * grad_z->v[1] + wl_n * cvn->grad_z->v[1];
        edgeZ = wl_p * z->value + wl_n * cvn->z->value + Vect::dotProduct(gradz_x, gradz_y, 0.0, df[i]);
      }
      else
      {
        //        edgeDepth =  h->value + Vect::dotProduct(*grad_h, r_e[i]);
        edgeZ =  z->value + Vect::dotProduct(*grad_z, r_e[i]);
      }

      //      edgeDepth = max(0.0, edgeDepth);
      //      faceDepth.value = edgeDepth;
      //      faceDepth.associatedValue = edgeDepth + ecz[i];

      double sn1 = snzMinEdge[i];
      double sn2 = snzMaxEdge[i];

      if(edgeZ <= sn1)
      {
        faceDepth.value = modelComponent->m_viscousSubLayerDepth;
        faceDepth.associatedValue = sn1 + modelComponent->m_viscousSubLayerDepth;
      }
      else if(edgeZ > sn1 && edgeZ <= sn2)
      {
        double num = edgeZ - sn1;
        faceDepth.value = (num * num) / (2.0 *(sn2 - sn1));
        faceDepth.associatedValue = edgeZ;
      }
      else
      {
        faceDepth.value = edgeZ - ecz[i];
        faceDepth.associatedValue = edgeZ;
      }
    }

    minNodeDepth = min(faceDepth.value, minNodeDepth);
  }
}

void TriCV::setFaceElevation(int face, double elevation)
{
  VarBC &faceDepth = faceDepths[face];

  double sn1 = snzMinEdge[face];
  double sn2 = snzMaxEdge[face];

  if(elevation <= sn1)
  {
    faceDepth.value = modelComponent->m_viscousSubLayerDepth;
    faceDepth.associatedValue = sn1 + modelComponent->m_viscousSubLayerDepth;
  }
  else if(elevation > sn1 && elevation <= sn2)
  {
    double num = elevation - sn1;
    faceDepth.value = (num * num) / (2.0 *(sn2 - sn1));
    faceDepth.associatedValue = elevation;
  }
  else
  {
    faceDepth.value = elevation - ecz[face];
    faceDepth.associatedValue = elevation;
  }
}

void TriCV::calculateWSEGradient()
{
  if(grad_z->isBC == false)
  {

    double dz_dx = 0.0;
    double dz_dy = 0.0;

    for(int i = 0 ; i < numEdges; i++)
    {
      int cvnIndex = nTris[i];

      if(faceDepths[i].isBC)
      {
        double z_e = faceDepths[i].associatedValue;

        dz_dx += z_e * r_eta[i] * e_n[i].v[0];
        dz_dy += z_e * r_eta[i] * e_n[i].v[1];

      }
      else if(cvnIndex  > -1)
      {
        TriCV *cvn = modelComponent->m_controlVolumes[cvnIndex];

        double wl_n = w_l[i];
        double wl_p = 1.0 - wl_n;

        double z_e = wl_p * z->value + wl_n * cvn->z->value;

        double gradz_x = wl_p * grad_z->v[0] + wl_n * cvn->grad_z->v[0];
        double gradz_y = wl_p * grad_z->v[1] + wl_n * cvn->grad_z->v[1];

        z_e += Vect::dotProduct(gradz_x, gradz_y, 0, df[i]);

        dz_dx += z_e * r_eta[i] * e_n[i].v[0];
        dz_dy += z_e * r_eta[i] * e_n[i].v[1];

      }
      else
      {
        double z_e =  z->value + Vect::dotProduct(*grad_z, r_e[i]);
        dz_dx += z_e * r_eta[i] * e_n[i].v[0];
        dz_dy += z_e * r_eta[i] * e_n[i].v[1];
      }
    }

    dz_dx /= area;
    dz_dy /= area;

    grad_z->v[0] = dz_dx;
    grad_z->v[1] = dz_dy;
  }
}

void TriCV::calculateDepthGradient()
{

  double dh_dx = 0.0;
  double dh_dy = 0.0;

  for(int i = 0 ; i < numEdges; i++)
  {
    int cvnIndex = nTris[i];

    if(faceDepths[i].isBC)
    {
      double h_e = max(0.0, faceDepths[i].value);

      dh_dx += h_e * r_eta[i] * e_n[i].v[0];
      dh_dy += h_e * r_eta[i] * e_n[i].v[1];
    }
    else if(cvnIndex  > -1)
    {
      TriCV *cvn = modelComponent->m_controlVolumes[cvnIndex];

      double wl_n = w_l[i];
      double wl_p = 1.0 - wl_n;

      double gradh_x = wl_p * grad_h->v[0] + wl_n * cvn->grad_h->v[0];
      double gradh_y = wl_p * grad_h->v[1] + wl_n * cvn->grad_h->v[1];

      double h_e = max(0.0, wl_p * h->value + wl_n * cvn->h->value + Vect::dotProduct(gradh_x, gradh_y, 0, df[i]));

      dh_dx += h_e * r_eta[i] * e_n[i].v[0];
      dh_dy += h_e * r_eta[i] * e_n[i].v[1];

    }
    else
    {
      double h_e = max(0.0, h->value + Vect::dotProduct(*grad_h, r_e[i]));
      dh_dx += h_e * r_eta[i] * e_n[i].v[0];
      dh_dy += h_e * r_eta[i] * e_n[i].v[1];
    }
  }

  dh_dx /= area;
  dh_dy /= area;

  grad_h->v[0] = dh_dx;
  grad_h->v[1] = dh_dy;
}

void TriCV::calculateZCorrGradient()
{

  double dz_dx = 0.0;
  double dz_dy = 0.0;

  for(int i = 0 ; i < numEdges; i++)
  {
    int cvnIndex = nTris[i];

    if(cvnIndex  > -1)
    {
      TriCV *cvn = modelComponent->m_controlVolumes[cvnIndex];

      double wl_n = w_l[i];
      double wl_p = 1.0 - wl_n;

      double z_e = wl_p * zCorrection + wl_n * cvn->zCorrection;

      double gradz_x = wl_p * grad_zcorr->v[0] + wl_n * cvn->grad_zcorr->v[0];
      double gradz_y = wl_p * grad_zcorr->v[1] + wl_n * cvn->grad_zcorr->v[1];

      z_e += Vect::dotProduct(gradz_x, gradz_y, 0, df[i]);

      dz_dx += z_e * r_eta[i] * e_n[i].v[0];
      dz_dy += z_e * r_eta[i] * e_n[i].v[1];

    }
    else
    {
      double z_e =  zCorrection + Vect::dotProduct(*grad_zcorr, r_e[i]);
      dz_dx += z_e * r_eta[i] * e_n[i].v[0];
      dz_dy += z_e * r_eta[i] * e_n[i].v[1];
    }
  }

  dz_dx /= area;
  dz_dy /= area;

  grad_zcorr->v[0] = dz_dx;
  grad_zcorr->v[1] = dz_dy;

}

void TriCV::calculateNodeElevations()
{

  for(int i = 0; i < numEdges ; i++)
  {
    double totalValue = 0;
    double totalDistance = 0;

    int numTri = numNodeTriangles[i];

    for(int j = 0; j < numTri ; j++ )
    {
      TriCV *cvn = modelComponent->m_controlVolumes[nodeTriangleIndexes[i][j]];
      double dist = nodeTriangleDistances[i][j];
      totalValue += cvn->z->value / (dist * dist);
      totalDistance += 1.0 / (dist * dist);
    }

    if(totalDistance)
    {
      nWSE[i] = max(snz[0], totalValue / totalDistance);
    }
    else
    {
      nWSE[i] =  max(snz[0], z->value + Vect::dotProduct(*grad_z, r_n[i]));
    }
  }


  //  minNodeDepth = max(0.0, nWSE[0]- nz[0]);

  //  for(int i = 1; i < numEdges; i++)
  //  {
  //    minNodeDepth = min(minNodeDepth, nWSE[i] - nz[i]);
  //  }

  //  nWSE[numEdges] = nWSE[0];
}

void TriCV::calculateVelocityGradient()
{
  //  //Green-Gauss
  for(int f = 0 ; f < 2 ; f++)
  {
    VectBC &gradVel = grad_vel[f];

    if(gradVel.isBC == false)
    {
      double du_dx = 0.0;
      double du_dy = 0.0;

      for(int i = 0;  i < numEdges; i++)
      {
        int cvnIndex = nTris[i];

        if(faceNormalVels[i].isBC && faceNormalVels[i].value)
        {
          double u_e = faceNormalVels[i].vel->v[f];
          du_dx += u_e * r_eta[i] * e_n[i].v[0];
          du_dy += u_e * r_eta[i] * e_n[i].v[1];
        }
        else if(cvnIndex > -1)
        {
          TriCV *cvn = modelComponent->m_controlVolumes[cvnIndex];
          double u_e = 0;

          double wl_n = w_l[i];
          double wl_p = 1.0 - wl_n;

          double gradVX = wl_p * gradVel.v[0] + wl_n * cvn->grad_vel[f].v[0];
          double gradVY = wl_p * gradVel.v[1] + wl_n * cvn->grad_vel[f].v[1];

          u_e = wl_p * vel[f].value + wl_n * cvn->vel[f].value + Vect::dotProduct(gradVX, gradVY, 0.0, this->df[i]);

          du_dx += u_e * r_eta[i] * e_n[i].v[0];
          du_dy += u_e * r_eta[i] * e_n[i].v[1];
        }
        else
        {
          double u_e = vel[f].value + Vect::dotProduct(gradVel, r_e[i]);
          du_dx += u_e * r_eta[i] * e_n[i].v[0];
          du_dy += u_e * r_eta[i] * e_n[i].v[1];
        }
      }

      du_dx /= area;
      du_dy /= area;

      gradVel.v[0] = du_dx;
      gradVel.v[1] = du_dy;
    }
  }
}

void TriCV::calculateNodeVelocities()
{
  for(int f = 0; f < 2 ; f++)
  {

    for(int i = 0; i < numEdges ; i++)
    {
      //      double totalValue = 0;
      //      double totalDistance = 0;

      //      int numTri = numNodeTriangles[i];

      //      for(int j = 0; j < numTri ; j++ )
      //      {
      //        TriCV *cvn = modelComponent->m_controlVolumes[nodeTriangleIndexes[i][j]];
      //        double dist = nodeTriangleDistances[i][j];
      //        totalValue = cvn->vel[f].value / dist;
      //        totalDistance = 1.0 / dist;
      //      }

      //      if(totalDistance)
      //      {
      //        nodeVels[f][i] = totalValue / totalDistance;
      //      }
      //      else
      //      {
      //        nodeVels[f][i] = vel[f].value;
      //      }

      for(int i = 0; i < numEdges; i++)
      {
        nodeVels[f][i] = vel[f].value + Vect::dotProduct(r_n[i],this->grad_vel[f]);
      }
    }

    nodeVels[f][numEdges] = nodeVels[f][0];

  }
}

/*!
 * \brief Rhie, C.M. and W.L. Chow, 1983. Numerical Study of the Turbulent
 * Flow Past an Airfoil with Trailing Edge Separation. AIAA Journal 21:1525–1532.
 * \param cv
 */
void TriCV::calculateFaceVelocities()
{
  double u_p_x = vel[0].value;
  double u_p_y = vel[1].value;

  for(int i = 0 ; i < numEdges ; i++)
  {
    FaceNormVelBC &faceNormVel =  faceNormalVels[i];

    if(faceNormVel.isBC == false)
    {
      faceNormVel.value = 0.0;
      faceNormVel.associatedValue = 0.0;

      int cvnIndex = nTris[i];

      if(cvnIndex > -1)
      {
        double v1 = 0.0, v2 = 0.0, v3 = 0.0, evel = 0.0 ,
            A_e_cv = 0.0, A_e_cvn = 0.0,
            A_e_cvy = 0.0, A_e_cvny = 0.0;

        TriCV* cvn = modelComponent->m_controlVolumes[cvnIndex];

        if(wetIndex && cvn->wetIndex)
        {

          double wl_n = w_l[i];
          double wl_p = 1.0 - wl_n;
          A_e_cv  = wl_p * area / velCoeffs[0][0];
          A_e_cvn  = wl_n * cvn->area / cvn->velCoeffs[0][0];

          A_e_cvy  = wl_p * area / velCoeffs[1][0];
          A_e_cvny  = wl_n * cvn->area / cvn->velCoeffs[1][0];


          //v1
          {
            double v1x = wl_p * u_p_x + wl_n * cvn->vel[0].value;
            double v1y = wl_p * u_p_y + wl_n * cvn->vel[1].value;

            double gradUx = wl_p * grad_vel[0].v[0] + wl_n * cvn->grad_vel[0].v[0];
            double gradUy = wl_p * grad_vel[0].v[1] + wl_n * cvn->grad_vel[0].v[1];

            v1x += Vect::dotProduct(gradUx, gradUy, 0.0, df[i]);

            double gradVx = wl_p * grad_vel[1].v[0] + wl_n * cvn->grad_vel[1].v[0];
            double gradVy = wl_p * grad_vel[1].v[1] + wl_n * cvn->grad_vel[1].v[1];

            v1y += Vect::dotProduct(gradVx, gradVy, 0.0, df[i]);

            v1 = Vect::dotProduct(v1x, v1y, 0, e_n[i]);
          }

          //v2
          {
            double tv2x = A_e_cv * h->value * grad_z->v[0] + A_e_cvn * cvn->h->value * cvn->grad_z->v[0];
            double tv2y = A_e_cvy * h->value * grad_z->v[1] + A_e_cvny * cvn->h->value * cvn->grad_z->v[1];
            v2 =  -modelComponent->g * Vect::dotProduct(tv2x, tv2y, 0, e_xi[i]);
          }

          //v3
          {
            v3  = -modelComponent->g * (A_e_cv + A_e_cvn) * h->value * (cvn->z->value - z->value) / r_xi_l[i];
          }

          evel = v1 - v2 + v3;
        }
        else if(cvn->wetIndex)
        {

          //          double vv1 = vel[0].value + Vect::dotProduct(grad_vel[0], r_e[i]);
          //          double vv2 = vel[1].value + Vect::dotProduct(grad_vel[1], r_e[i]);
          //          v1 = Vect::dotProduct(vv1, vv2, 0.0, e_n[i]);
          //          evel = min(0.0, v1);

          double wl_n = w_l[i];
          double wl_p = 1.0 - wl_n;

          double gradUx = wl_p * grad_vel[0].v[0] + wl_n * cvn->grad_vel[0].v[0];
          double gradUy = wl_p * grad_vel[0].v[1] + wl_n * cvn->grad_vel[0].v[1];

          double gradVx = wl_p * grad_vel[1].v[0] + wl_n * cvn->grad_vel[1].v[0];
          double gradVy = wl_p * grad_vel[1].v[1] + wl_n * cvn->grad_vel[1].v[1];

          double ux = wl_p * vel[0].value + wl_n * cvn->vel[0].value + Vect::dotProduct(gradUx, gradUy, 0.0, df[i]);
          double uy = wl_p * vel[1].value + wl_n * cvn->vel[1].value + Vect::dotProduct(gradVx, gradVy, 0.0, df[i]);

          v1 = Vect::dotProduct(ux,uy,0.0, e_n[i]);
          evel = min(0.0, v1);
        }
        else if(wetIndex)
        {

          //          double vv1 = vel[0].value + Vect::dotProduct(grad_vel[0], r_e[i]);
          //          double vv2 = vel[1].value + Vect::dotProduct(grad_vel[1], r_e[i]);
          //          v1 = Vect::dotProduct(vv1, vv2, 0.0, e_n[i]);
          //          evel = max(0.0, v1);


          double wl_n = w_l[i];
          double wl_p = 1.0 - wl_n;

          double gradUx = wl_p * grad_vel[0].v[0] + wl_n * cvn->grad_vel[0].v[0];
          double gradUy = wl_p * grad_vel[0].v[1] + wl_n * cvn->grad_vel[0].v[1];

          double gradVx = wl_p * grad_vel[1].v[0] + wl_n * cvn->grad_vel[1].v[0];
          double gradVy = wl_p * grad_vel[1].v[1] + wl_n * cvn->grad_vel[1].v[1];

          double ux = wl_p * vel[0].value + wl_n * cvn->vel[0].value + Vect::dotProduct(gradUx, gradUy, 0.0, df[i]);
          double uy = wl_p * vel[1].value + wl_n * cvn->vel[1].value + Vect::dotProduct(gradVx, gradVy, 0.0, df[i]);

          v1 = Vect::dotProduct(ux,uy,0.0, e_n[i]);
          evel = max(0.0, v1);
        }

        faceNormVel.value = evel;
        faceNormVel.associatedValue = evel * faceDepths[i].value * r_eta[i];
      }
    }
  }
}

double TriCV::verifyFaceVelocity(TriCV *cv, TriCV *cvn, int faceIndex, double faceVelocity)
{
  if (cv->wetIndex && cvn->wetIndex)
    return faceVelocity;
  else if (cv->wetIndex && cvn->wetIndex == 0)
    return max(faceVelocity,0.0);
  else if (cv->wetIndex == 0 && cvn->wetIndex)
    return min(faceVelocity,0.0);

  return 0.0;
}

void TriCV::interpolateNodeEddyViscosities()
{
  for(int i = 0; i < numEdges ; i++)
  {
    double valTotal = 0.0;
    double distTotal = 0.0;

    int numTri = numNodeTriangles[i];

    for(int j = 0; j < numTri ; j++)
    {
      TriCV *cvn = modelComponent->m_controlVolumes[nodeTriangleIndexes[i][j]];
      double dist = nodeTriangleDistances[i][j];
      valTotal += cvn->eddyViscosity / dist;
      distTotal += 1.0 / dist;
    }

    nodalEddyViscosity[i] = valTotal / distTotal;

  }
}

double TriCV::getVFRDepth(double wse)
{
  double zz1 = snz[0];
  double zz2 = snz[1];
  double zz3 = snz[2];

  if(wse <= zz1)
  {
    return  modelComponent->m_viscousSubLayerDepth;
  }
  else
  {
    if(wse > zz1 && wse <= zz2)
    {
      double num = wse - zz1;
      double dp = num * num * num / (3.0 *(zz2-zz1)*(zz3-zz1));
      return dp;
    }
    else if(wse > zz2 && wse <= zz3)
    {
      double dp = (wse*wse + wse*zz3 - 3.0*wse*zz1 - zz3*zz2 + zz1*zz2 + zz1*zz1) / (3.0*(zz3-zz1));
      return dp;
    }
    else
    {
      return (wse - cz);
    }
  }
}

double TriCV::getBedZ(double x, double y)
{
  double v1[3] = { vertices[1]->x() - vertices[0]->x(),
                   vertices[1]->y() - vertices[0]->y(),
                   nz[1] - nz[0]};

  double v2[3] = { vertices[2]->x() - vertices[0]->x(),
                   vertices[2]->y() - vertices[0]->y(),
                   nz[2] - nz[0]};

  double norm[3];

  TriCV::crossProduct(v1,v2,norm);
  TriCV::normalizeVector(norm);

  double d = -norm[0] * vertices[0]->x()
             -norm[1] * vertices[0]->y()
             -norm[2] * vertices[0]->z();


  double z = (norm[0] * x + norm[1] * y + d) /(-norm[2]);

  return z;
}

double TriCV::getWSEZ(double x, double y)
{
  //  double v1[3] = { vertices[1]->x() - vertices[0]->x(),
  //                   vertices[1]->y() - vertices[0]->y(),
  //                   nWSE[1] - nWSE[0]};

  //  double v2[3] = { vertices[2]->x() - vertices[0]->x(),
  //                   vertices[2]->y() - vertices[0]->y(),
  //                   nWSE[2] - nWSE[0]};

  //  double norm[3];

  //  TriCV::crossProduct(v1,v2,norm);
  //  TriCV::normalizeVector(norm);

  //  double d = -norm[0] * vertices[0]->x()
  //             -norm[1] * vertices[0]->y()
  //             -norm[2] * nWSE[0];


  //  double z = (norm[0] * x + norm[1] * y + d) /(-norm[2]);

  //  return z;
  return 0.0;
}

bool TriCV::intersectPoint(double x1, double y1, double x2, double y2,
                           double x3, double y3, double x4, double y4,
                           double &x, double &y)
{
  double denom = ((y4-y3)*(x2-x1) - (x4-x3)*(y2-y1));

  if(fabs(denom - 0.0) > std::numeric_limits<double>::epsilon())
  {
    double u = ((x4-x3)*(y1-y3) - (y4-y3)*(x1-x3)) / denom;

    x = x1 + u * (x2 - x1);
    y = y1 + u * (y2 - y1);

    if(u >= 0 && u <= 1.0)
    {
      return true;
    }
  }

  return false;
}

bool TriCV::intersectPoint(const Vect& v1, const Vect &v2,
                           const Vect &v3, const Vect &v4, Vect &p)
{
  double denom = ((v4.y() - v3.y()) * (v2.x() - v1.x()) - (v4.x() - v3.x()) * (v2.y() - v1.y()));

  if(fabs(denom - 0.0) > std::numeric_limits<double>::epsilon())
  {
    double u = ((v4.x() - v3.x()) * (v1.y() - v3.y()) - (v4.y() - v3.y()) * (v1.x() - v3.x())) / denom;

    p.v[0] = v1.x() + u * (v2.x() - v1.x());
    p.v[1] = v1.y() + u * (v2.y() - v1.y());
    p.v[2] = v1.z() + u * (v2.z() - v1.z());

    if(u >= 0 && u <= 1.0)
    {
      return true;
    }
  }

  return false;
}

void TriCV::lsGradReconstruction(double value, Vect *distances[],
                                 const double values[], int rowCount, double output[])
{

  output[0] = 0;
  output[1] = 0;

  double *a = new double[2 * rowCount];
  double *b = new double[rowCount] ;

  for(int i = 0 ; i  < rowCount ; i++)
  {
    const Vect *dist = distances[i];
    b[i] = (values[i] - value);
    a[i + rowCount * 0] = dist->v[0];
    a[i + rowCount * 1] = dist->v[1];
  }

  qr_solve(output,rowCount,2,a,b);

  delete[] a;
  delete[] b;
}

void TriCV::transpose(int m, int n, double *a, double *a_t)
{
  for(int i = 0; i < m ; i++)
  {
    for(int j = 0; j < n ; j++)
    {
      a_t[j + i * n] = a[i + j * m];
    }
  }
}

void TriCV::multiply(int am, int an, double *a, int bm, int bn, double *b, double *x)
{
  assert(an == bm);

  for(int i = 0 ; i < am; i++)
  {
    for(int j = 0; j < bn; j++)
    {
      x[i + am*j]  = 0.0;

      for(int k = 0; k < an; k++)
      {
        x[i + am*j] += a[i + k * am] * b[k + an*j];
      }
    }
  }

}

void TriCV::normalizeVector(double v[])
{
  double x= v[0];
  double y= v[1];
  double z= v[2];

  double denom = sqrt(x*x + y*y + z*z);

  v[0] = x /denom;
  v[1] = y /denom;
  v[2] = z /denom;
}

double TriCV::dotProduct(double u[], double v[])
{
  return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
}

void TriCV::crossProduct(double u[], double v[], double out[])
{
  out[0] = u[1]*v[2] - u[2]*v[1];
  out[1] = u[2]*v[0] - u[0]*v[2];
  out[2] = u[0]*v[1] - u[1]*v[0];
}

void TriCV::printDetails()
{
  printf("\ncv index: %i\n"
         "wetIndex: %i\tcontIndex: %i\n"
         "inflow: %f\tfriction_x: %f\tfriction_y: %f\n"
         "h: %f\tgradh_x: %f\tgradh_y: %f\n"
         "z: %f\tgradz_x: %f\tgradz_y: %f\n"
         "u: %f\tgradu_x: %f\tgradu_y :%f\n"
         "v :%f\tgradv_x: %f\tgradv_y: %f\n"
         "fvel_1: %f\tfvel_2: %f\tfvel_3: %f\n",
         index,
         wetIndex, contIndex,
         inflowOutflow, friction[0], friction[1],
      h->value, grad_h->v[0], grad_h->v[1],
      z->value, grad_z->v[0], grad_z->v[1],
      vel[0].value, grad_vel[0].v[0], grad_vel[0].v[1],
      vel[1].value, grad_vel[1].v[0], grad_vel[1].v[1],
      faceNormalVels[0].value, faceNormalVels[1].value, faceNormalVels[2].value);
}

void TriCV::calculateInitialWSEGradient()
{
  if(grad_z_init == nullptr)
  {
    grad_z_init = new Vect();
    double *values = new double[numEdges];
    Vect  **distances = new Vect*[numEdges];

    for(int i = 0; i < numEdges; i++)
    {
      values[i] = ecz[i];
      distances[i] = &r_e[i];
    }

    lsGradReconstruction(cz, distances, values, numEdges, grad_z_init->v);

    delete[] values;
    delete[] distances;
  }

  grad_z->v[0] = grad_z_init->v[0];
  grad_z->v[1] = grad_z_init->v[1];
  grad_z->v[2] = 0;

}
