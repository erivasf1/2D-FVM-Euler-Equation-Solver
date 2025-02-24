//User-defined functions
#include "EulerOperator.h" 

// EULER1D DEFINITIONS

Euler1D::Euler1D(vector<double> &coords)
  : xcoords(coords) {}

//-----------------------------------------------------------
void Euler1D::SetInitialConditions(array<double,3> &init_val,array<double,3>* &field){

  for (int i=0;i<(int)xcoords.size();i++){
    for (int n=0;n<(int)init_val.size();n++) {
      field[i][n] = init_val[n]; //setting the flow quantities (primitive) to the initial conditions

    }    
  }

  return;

}
//-----------------------------------------------------------

Euler1D::~Euler1D(){}
