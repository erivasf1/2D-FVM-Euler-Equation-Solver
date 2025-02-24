// Responsible for creating a Mesh (1D for now)
#ifndef _EULEROPERATOR_H_
#define _EULEROPERATOR_H_
#include "ExactNozzle.h"
#include "MeshGen.h"

class Euler1D {
  vector<double> &xcoords;

  public:
  Euler1D(vector<double> &coords);

  void SetBoundaryConditions(); //TODO
  void SetInitialConditions(array<double,3> &init_val,array<double,3>* &field); //Complete (tested)
  void FluxComputation();//TODO
  void JamesonDamping();//TODO
 

  ~Euler1D();


};




#endif
