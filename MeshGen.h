// Responsible for creating a Mesh (1D for now)
#ifndef _MESHGEN_H_
#define _MESHGEN_H_
#include "ExactNozzle.h"

class MeshGen1D {
  double xmin,xmax;
  int cellnumber;

  public:
  MeshGen1D(double &a, double &b, int &c);
 
  void GenerateMesh(vector<double> &xcoords);

  ~MeshGen1D();


};




#endif
