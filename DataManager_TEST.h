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

  void AllOutputPrimitiveVariables(vector<array<double,3>> &Field,Euler1D &Euler,string filename,bool cond,int iter,vector<double> &xcoords);

  void OutputLocalResiduals(vector<array<double,3>> &Resid,const char *filename); //TODO

  void OutputResidualTerms(array<double,3> F_right,array<double,3> F_left,double S,array<double,3> D2_right,array<double,3> D2_left,array<double,3> D4_right,array<double,3> D4_left,const char* filename);

  void ComputeCellAveragedSol(vector<array<double,3>> &SolFace,vector<array<double,3>> &SolCell,vector<double> &xcoords,double dx);

  ~SpaceVariables1D();


};




#endif
