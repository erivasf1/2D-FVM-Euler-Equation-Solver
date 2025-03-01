// Responsible for creating unit tests for functions used in the solver
#ifndef _TEST_H_
#define _TEST_H_
#include "ExactNozzle.h"
#include "MeshGen.h"
#include "EulerOperator.h"

class Test {
  array<double,3>* &field;

  public:
  Test();
  
  void Test2ndOrderDamping(array<double,3>* field,int& loc);
  void Test4thOrderDamping();
  void TestComputeResidual();
 

  ~Test();


};




#endif
