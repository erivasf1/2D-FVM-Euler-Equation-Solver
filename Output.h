// Responsible for outputting results in a clean format
#ifndef _OUTPUT_H_
#define _OUTPUT_H_
#include <fstream> 
#include <cmath> 

#include "ExactNozzle.h"
#include "MeshGen.h"
#include "DataManager.h"
#include "TimeIntegrator.h"

using namespace std;

class Output {
//  array<double,3>* &field;

  public:
  //Output(array<double,3>* &field);
  Output();
  
  void PrintResidualNorm(int &cellnum,int &n);

  void DiscretizationErrorNorms(vector<array<double,3>>* &field,vector<array<double,3>>* &exact_field,vector<array<double,3>>* &errors,SpaceVariables1D* &sols);

  void CalculateOrderofAccuracy(const char *filename_read,const char *filename_write); //creates a new file containing the order of accuracy value given the discretization error file.txt
 
  void OutputPrimitiveVariables(vector<array<double,4>>* &field,string filename,bool cond,int iter,vector<double> &xcoords,vector<double> &ycoords,int cell_number,int imax,int jmax); //outputs primitive variables in tecplot format

  void OutputManufacturedSourceTerms(vector<array<double,4>>* &field,string filename,bool cond,int iter,vector<double> &xcoords,vector<double> &ycoords,int cell_number,int imax,int jmax); //outputs primitive variables in tecplot format
  //void ConvertToDatFile(const char*filename_read,const char *filename_write); //TODO: creates a .dat file of a given .txt file

  ~Output();


};




#endif
