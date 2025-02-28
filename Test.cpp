//User-defined functions
#include "Test.h" 

// TEST DEFINITIONS

//-----------------------------------------------------------
void Test::Test2ndOrderDamping(array<double,3>* field,int& loc){

  Tools::print("Testing 2nd Order Damping fcn.\n");
  array<double,3> D = Euler1D::Compute2ndOrderDamping(field,loc);

  print("2ndOrderDamping Term: [%f,%f,%f]^T",D[0],D[1],D[2]);

  return;


}


//-----------------------------------------------------------

Test::~Test(){}
