// Responsible for creating a Mesh (1D for now)
#ifndef _MESHGEN_H_
#define _MESHGEN_H_
#include "ExactNozzle.h"
#include <fstream>
#include <iostream>

class MeshGen1D { //creates a uniform mesh
  double xmin,xmax;
  int cellnumber;

  public:
  MeshGen1D(double &a, double &b, int &c);

  static double GetCellVolume(int &loc,double &dx,vector<double> &xcoords);
 
  void GenerateMesh(vector<double> &xcoords);

  void OutputNozzleAreas(vector<double> &xcoords,const char *filename);

  ~MeshGen1D();


};

class MeshGen2D { //reads in a non-uniform 2D mesh
  double xmin,xmax;
  double ymin,ymax;
  int cellnumber;
  const char* filename;

  public:
  MeshGen2D(const char* name);

  vector<double> xcoords,ycoords,zcoords;
  int imax,jmax,kmax;

  void ReadMeshFile();

  double GetCellVolume(int &i,double &dx,int &j,double &dy,vector<double> &xcoords,vector<double> &ycoords);


};

#endif
