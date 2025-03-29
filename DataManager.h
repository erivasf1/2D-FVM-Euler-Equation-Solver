// Responsible for creating a Mesh (1D for now)
#ifndef _DATAMANAGER_H_
#define _DATAMANAGER_H_
#include "MeshGen.h"
#include "TimeIntegrator.h"


//SpaceVariable Class stores the values of ALL flow quantities at all points
class SpaceVariables1D {
  int cell_number;

  public:

  SpaceVariables1D();

  array<double,3> ComputeSolutionNorms(vector<array<double,3>>* &resid);

  void OutputPrimitiveVariables(vector<array<double,3>>* &field,Euler1D* &euler,const char *filename);

  void AllOutputPrimitiveVariables(vector<array<double,3>>* &field,Euler1D* &euler,string filename,bool cond,int iter,vector<double> &xcoords);

  void OutputLocalResiduals(vector<array<double,3>> &Resid,const char *filename); 

  void OutputResidualTerms(array<double,3> F_right,array<double,3> F_left,double S,array<double,3> D2_right,array<double,3> D2_left,array<double,3> D4_right,array<double,3> D4_left,const char* filename); //Not needed for now

  void ComputeCellAveragedSol(vector<array<double,3>>* &cell_faces,vector<array<double,3>>* &cell_sols,vector<double> &xcoords,double dx);

  ~SpaceVariables1D();

};




#endif
