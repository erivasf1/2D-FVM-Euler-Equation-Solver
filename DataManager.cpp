//User-defined functions
#include "DataManager.h" 

// SPACEVARIABLE1D DEFINITIONS

SpaceVariables1D::SpaceVariables1D(int &c,array<double,3>* &field)
  : cell_number(c) {

  //TODO: Account for ghost cells!
  vector<array<double,3>> Field(c); //creating vector that stores all values of conservative variables
  field = Field.data(); //assinging array pointer (field) to point to 1st array element of the 1st point 
}
//---------------------------------------------------------

void SpaceVariables1D::ConvertToConservative() {



}

//---------------------------------------------------------
void SpaceVariables1D::ConvertToPrimitive() {


}

//---------------------------------------------------------
SpaceVariables1D::~SpaceVariables1D(){}
