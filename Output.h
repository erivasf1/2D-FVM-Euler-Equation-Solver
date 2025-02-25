// Responsible for creating a Mesh (1D for now)
#ifndef _OUTPUT_H_
#define _OUTPUT_H_
#include "ExactNozzle.h"
#include "MeshGen.h"

class Output {
  array<double,3>* &field;

  public:
  Output(array<double,3>* &field);
  
  void PrintResidualNorm(int &cellnum,int &n);
 

  ~Output();


};




#endif
