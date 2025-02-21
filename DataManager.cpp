//User-defined functions
#include "DataManager.h" 

// SPACEVARIABLE1D DEFINITIONS

SpaceVariable1D::SpaceVariable1D(int &c)
  : cell_number(c) {

  //TODO:
  vector<array<double,3>> Field(c); //creating vector that stores all values of conservative variables
  //Setup();
}
//---------------------------------------------------------

void SpaceVariable1D::Setup() {

  //vector<array<double,3>> Field(cell_number); //creating vector that stores all values of conservative variables
  //field = Field.data(); //assigning pointer to point to each coord

}


//---------------------------------------------------------
SpaceVariable1D::~SpaceVariable1D(){}
