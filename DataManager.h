// Responsible for creating a Mesh (1D for now)
#ifndef _DATAMANAGER_H_
#define _DATAMANAGER_H_
#include "MeshGen.h"

//SpaceVariable Class stores the values of ALL flow quantities at all points
class SpaceVariables1D {
  int cell_number;

  public:

  SpaceVariables1D(int &c,vector<array<double,3>> &Field,array<double,3>* &field);
 
  //void Pointer(vector<array<double,3>> &Field); //sets up domain by assigning each point a 3D array to store conservative variables
  void ConvertToConservative(); //converts conservative variable values to primitive variable values
  void ConvertToPrimitive(); //converts primitive variable values to primitive variable values
  void OutputSolution();

  ~SpaceVariables1D();


};




#endif
