//User-defined functions
#include "DataManager.h" 

// SPACEVARIABLE1D DEFINITIONS

SpaceVariables1D::SpaceVariables1D(int &c,vector<array<double,3>> &Field,array<double,3>* &field)
  : cell_number(c) {

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
array<double,3> SpaceVariables1D::ComputeSolutionNorms(array<double,3>* &field){

  array<double,3> norm{0.0,0.0,0.0};

  //using L2 norm  
  for (int i=0;i<cell_number;i++){
    Tools::print("Cell number: %d\n",i);
    //continuity
    Tools::print("continuity res.: %e\n",field[i][0]);
    norm[0]+= pow(field[i][0],2);
    //x-mom
    Tools::print("x-mom. res.: %e\n",field[i][1]);
    norm[1]+= pow(field[i][1],2);
    //energy
    Tools::print("energy res.: %e\n",field[i][2]);
    norm[2]+= pow(field[i][2],2);

  }
  norm[0] = sqrt(norm[0]);
  norm[1] = sqrt(norm[1]);
  norm[2] = sqrt(norm[2]);

  return norm;
}

//---------------------------------------------------------
SpaceVariables1D::~SpaceVariables1D(){}
