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
void Euler1D::SetBoundaryConditions(vector<array<double,3>> &Field,array<double,3> &init){

  //Inserting ghost cells for inflow and outflow locations
  //iterator it = Field.begin();
  Field.insert(Field.begin(),init); //!< setting ghost cells to initial conditions 
  Field.push_back(init); 
  
}




//-----------------------------------------------------------

Euler1D::~Euler1D(){}
