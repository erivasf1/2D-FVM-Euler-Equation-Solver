//User-defined functions
#include "Test.h" 

// TEST DEFINITIONS

//-----------------------------------------------------------
void Test::Test2ndOrderDamping(array<double,3>* field,int& loc){

  Euler1D Euler;
  Tools::print("Testing 2nd Order Damping fcn.\n");
  array<double,3> D = Euler.Compute2ndOrderDamping(field,loc);

  Tools::print("2ndOrderDamping Term: [%f,%f,%f]^T",D[0],D[1],D[2]);

  return;


}


//-----------------------------------------------------------
void Test::TestComputeResidual(vector<array<double,3>> &Test){




}


//-----------------------------------------------------------

Test::~Test(){}
