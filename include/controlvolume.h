/*! \file   controlvolume.h
 *  \author Caleb Amoa Buahin <caleb.buahin@gmail.com>
 *  \version   1.0.0.0
 *  \section   Description
 *  \section License
 *  cvwseoutput.h, associated files and libraries are free software;
 *  you can redistribute it and/or modify it under the terms of the
 *  Lesser GNU General Public License as published by the Free Software Foundation;
 *  either version 3 of the License, or (at your option) any later version.
 *  The FVHMComponent library and its associated files is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.(see <http://www.gnu.org/licenses/> for details)
 *  \date 2014-2016
 *  \pre
 *  \bug
 *  \warning
 *  \todo
 */

#ifndef GEOMDATAOBJECTS_H
#define GEOMDATAOBJECTS_H

#include "matrix.h"
#include "fvhmcomponent_global.h"
#include <spatial/point.h>

class HCTriangle;
class FVHMComponent;

#define btoa(x) ((x)?"true":"false")

struct FVHMCOMPONENT_EXPORT VectBC : public Vect
{

    VectBC(const Vect &vect):
      Vect(vect)
    {

    }

    VectBC(const HCPoint &point)
      :Vect(point)
    {

    }

    VectBC(double x = 0, double y = 0, double z = 0)
      :Vect(x,y,z)
    {

    }

    VectBC &operator =(const Vect &vect)
    {
      Vect::operator =(vect);
      return *this;
    }



    bool isBC = false;
    bool hasValue = true;
};

struct FVHMCOMPONENT_EXPORT VarBC
{
    double value = 0;

    double associatedValue = 0;

    bool isBC = false;

    void copy(const VarBC & var)
    {
      value = var.value;
      isBC = var.isBC;
    }

    void print()
    {
      printf("Value:%f, isBC=%i\n",value,isBC);
    }
};

struct FVHMCOMPONENT_EXPORT FaceNormVelBC : public VarBC
{
  public:

    enum VelCalculationMode
    {
      PreCalculated,
      DotProductEdgeVel,
      CellVel,
      CellGradient,
      ZeroGradient
    };

    bool calculateFromQ = false;

    bool calculateWallShearStress = false;

    VelCalculationMode velocityCalculateMode = VelCalculationMode::PreCalculated;

    Vect *vel = nullptr;

    FaceNormVelBC():
      VarBC()
    {

    }

    ~FaceNormVelBC()
    {
      if(vel != nullptr)
      {
        delete vel;
        vel = nullptr;
      }
    }

    void initialiazeVelocityVariable()
    {
      if(vel == nullptr)
      {
        vel = new Vect(0,0,0);
      }
    }

};

struct FVHMCOMPONENT_EXPORT TriCV
{
    int index = 0;

    int constCount = 0;

    static int numEdges;

    int totalnTris = 0;

    //Neighbouring triangles
    int *nTris = nullptr;

    int *orderedNTris = nullptr;

    int *orderedNTrisIndexes = nullptr;

    int *numNodeTriangles = nullptr;

    int **nodeTriangleIndexes = nullptr;

    double **nodeTriangleDistances = nullptr;

    //    Vect **nodeTriangleDistancesV = nullptr;

    HCTriangle *cell = nullptr;

    HCVertex **vertices = nullptr;

    double area = 0.0;

    VarBC *z = nullptr;

    VarBC *prevZ = nullptr;

    double maxZ =  std::numeric_limits<double>::lowest();

    VarBC *h = nullptr;

    VarBC *prevH = nullptr;

    double maxH =  std::numeric_limits<double>::lowest();

    VarBC *vel = nullptr;

    VectBC *grad_z = nullptr;

    Vect *grad_z_init = nullptr;

    Vect *grad_h = nullptr;

    Vect *grad_zcorr = nullptr;

    VectBC *grad_vel = nullptr;

    double *maxVel =  nullptr;

    VarBC *prevVel = nullptr;

    //1 for wet, 0 for dry
    int wetIndex = 0;

    int wetCellIndex = 0;

    int contIndex = 0;

    int contCellIndex = 0;

    //Unit normal vectors for edge
    Vect *e_n = nullptr;

    //Unit normal vectors from this cells centroid to adjacent cell centroid.
    Vect *e_xi = nullptr;

    //Vector distances from edge centers to this cells centroid.
    Vect *r_e =  nullptr;

    Vect *r_e_cvn = nullptr;

    //distance to edge centroid
    double *r_e_l = nullptr;

    //normal distance from center of cell to  wall
    double *r_e_l_p = nullptr;

    //Distances from cell nodes to this cells centroid
    Vect *r_n = nullptr;

    //Edge Skewness
    Vect *df = nullptr;

    //Lengths of cell edges
    double *r_eta = nullptr;

    //Lengths of adjacent cell centroids to this cells centroid
    Vect *r_xi = nullptr;

    //
    double *r_xi_dot_e_n = nullptr;

    //Lengths of adjacent cell centroids to this cells centroid
    double *r_xi_l = nullptr;

    //inverse weighted length to adjancent cell with distances calculated from shared edge centroids.
    double *w_l = nullptr;

    //Edge normal velocities
    FaceNormVelBC *faceNormalVels = nullptr;

    //Node velocities
    double **nodeVels = nullptr;

    //center z
    double cz = 0;

    //center
    Vect *center = nullptr;

    //sorted node zs
    double *snz = nullptr;

    double *snzMinEdge = nullptr;

    double *snzMaxEdge = nullptr;

    //unsorted node elevations
    double *nz = nullptr;

    //edge center elevations
    double *ecz = nullptr;

    //current node water surface elevation
    double *nWSE = nullptr;

    //Edge depths. Associated data is the elevation estimated for the edge
    VarBC *faceDepths = nullptr;

    //Direct diffusion const term
    double *dir_diff_term = nullptr;

    //Cross diffusion const term
    double *cross_diff_term = nullptr;

    //Inflow boundary conditions
    double inflowOutflow = 0;

    //External force
    Vect *externalForce = nullptr;

    //Min node depth
    double minNodeDepth = 0;

    double mannings = 0.03;

    //Turbulent viscosity
    double eddyViscosity = 0.0;

    double *nodalEddyViscosity = nullptr;

    double *elevationFactors = nullptr;

    double *velResidual = nullptr;

    double *velResidualIter = nullptr;

    double contResidualIter = 0;

    double contResidualTotal = 0;

    double zCorrection = 0.0;

    double **velCoeffs = nullptr;

    //Previous iteration velocities
    double *prevIterVel = nullptr;

    //Storage depths
    double storageDepth = 0.001;

    double *friction = nullptr;

    double **wallShearFriction = nullptr;

    FVHMComponent *modelComponent = nullptr;

    bool hasEdgeDepthBC = false;

    double totalExternalInflow = 0.0;

    double totalExternalOutflow = 0.0;

    int fillMode;

    static int gradientCalcMode;

    int isWetOrHasWetNeigh = 0;

  public:

    TriCV(HCTriangle *cell, FVHMComponent* modelComponent);

    ~TriCV();

    void calculateAdjacentCellParams();

    double getCourantFactor() const;

    void setVFRWSE(double wse, bool prev = false);

    void setVFRDepth(double depth, bool prev = false);

    void copyVariablesToPrev();

    void calculateEdgeDepths();

    void setFaceElevation(int face, double elevation);

    void calculateWSEGradient();

    void calculateDepthGradient();

    void calculateZCorrGradient();

    void calculateNodeElevations();

    void calculateVelocityGradient();

    void calculateNodeVelocities();

    void calculateFaceVelocities();

    static double verifyFaceVelocity(TriCV *cv, TriCV *cvn, int faceIndex, double faceVelocity);

    void interpolateNodeFaceVelocityFactors();

    void interpolateNodeEddyViscosities();

    double getVFRDepth(double wse);

    double getBedZ(double x, double y);

    double getWSEZ(double x, double y);

    static bool intersectPoint(double x1, double y1, double x2,  double y2, double x3, double y3, double x4, double y4, double &x, double &y);

    static bool intersectPoint(const Vect& v1, const Vect &v2, const Vect &v3, const Vect &v4, Vect &p);

    static void lsGradReconstruction(double value, Vect *distances[], const double values[], int rowCount, double output[]);

    static void transpose(int m, int n, double* a, double *a_t);

    static void multiply(int am, int an, double* a, int bm, int bn, double *b, double *x);

    static void normalizeVector(double v[]);

    static double dotProduct(double u[] , double v[]);

    static void crossProduct(double u[], double v[], double out[]);

    void printDetails();

    void calculateInitialWSEGradient();
};

#endif // GEOMDATAOBJECTS_H
