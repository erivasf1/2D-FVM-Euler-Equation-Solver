// Responsible for creating a Mesh (1D for now)
#ifndef _DATAMANAGER_TEST_H_
#define _DATAMANAGER_TEST_H_
#include "MeshGen.h"
#include "EulerOperator_TEST.h"
#include <iostream>
#include <fstream>

using namespace std;

//SpaceVariable Class stores the values of ALL flow quantities at all points
class SpaceVariables1D {
  int cell_number;

  public:

  SpaceVariables1D();
 
  array<double,3> ComputeSolutionNorms(vector<array<double,3>> &Field);

  void OutputPrimitiveVariables(vector<array<double,3>> &Field,Euler1D &Euler,const char *filename);

  void OutputResidualTerms(array<double,3> F_right,array<double,3> F_left,double S,array<double,3> D2_right,array<double,3> D2_left,array<double,3> D4_right,array<double,3> D4_left,const char* filename);

  ~SpaceVariables1D();


};




#endif
