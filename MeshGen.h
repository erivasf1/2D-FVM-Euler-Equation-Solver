// Responsible for creating a Mesh (1D for now)
#ifndef _MESHGEN_H_
#define _MESHGEN_H_
#include "ExactNozzle.h"
#include <fstream>
#include <iostream>

class MeshGen1D {
  double xmin,xmax;
  int cellnumber;

  public:
  MeshGen1D(double &a, double &b, int &c);

  static double GetCellVolume(int &loc,double &dx,vector<double> &xcoords);
 
  void GenerateMesh(vector<double> &xcoords);

  void OutputNozzleAreas(vector<double> &xcoords,const char *filename);

  ~MeshGen1D();


};




#endif
