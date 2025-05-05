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

  void DiscretizationErrorNorms(vector<array<double,3>>* &field,vector<array<double,3>>* &exact_field,vector<array<double,3>>* &errors,SpaceVariablesBASE* &sols);

  void CalculateOrderofAccuracy(const char *filename_read,const char *filename_write); //creates a new file containing the order of accuracy value given the discretization error file.txt
 
  //void ConvertToDatFile(const char*filename_read,const char *filename_write); //TODO: creates a .dat file of a given .txt file

  ~Output();


};




#endif
