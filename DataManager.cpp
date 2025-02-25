//User-defined functions
#include "DataManager.h" 

// SPACEVARIABLE1D DEFINITIONS

SpaceVariables1D::SpaceVariables1D(int &c,vector<array<double,3>> &Field,array<double,3>* &field)
  : cell_number(c) {

  //TODO: Account for ghost cells!
  //Field(cell_number); //creating vector that stores all values of conservative variables
  field = Field.data(); //assinging array pointer (field) to point to 1st element of vector which is an array
}
//---------------------------------------------------------

void SpaceVariables1D::ConvertToConservative() {



}

//---------------------------------------------------------
void SpaceVariables1D::ConvertToPrimitive() {


}

//---------------------------------------------------------
SpaceVariables1D::~SpaceVariables1D(){}
