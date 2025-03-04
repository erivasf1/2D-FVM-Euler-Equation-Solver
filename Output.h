// Responsible for outputting results in a clean format
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
