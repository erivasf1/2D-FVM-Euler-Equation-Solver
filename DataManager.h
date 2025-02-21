// Responsible for creating a Mesh (1D for now)
#ifndef _DATAMANAGER_H_
#define _DATAMANAGER_H_
#include "ExactNozzle.h"

//SpaceVariable Class stores the values of ALL flow quantities at all points
class SpaceVariable1D {
  int cell_number;

  public:

  //array<double,3>* field;

  SpaceVariable1D(int &c);
 
  void Setup(); //sets up domain by assigning each point a 3D array to store conservative variables
  void ConvertToPrimitive(); //converts conservative variable values to primitive variable values
  void OutputSolution();

  ~SpaceVariable1D();


};




#endif
