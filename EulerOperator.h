// Responsible for creating a Mesh (1D for now)
#ifndef _EULEROPERATOR_H_
#define _EULEROPERATOR_H_
#include "ExactNozzle.h"
#include "MeshGen.h"

class Euler1D {
  vector<double> &xcoords;

  public:
  Euler1D(vector<double> &coords);
 

  ~Euler1D();


};




#endif
