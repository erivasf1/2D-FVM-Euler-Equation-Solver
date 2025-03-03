// Responsible for creating a Mesh (1D for now)
#ifndef _TIMEINTEGRATOR_H_
#define _TIMEINTEGRATOR_H_
#include "ExactNozzle.h"
#include "EulerOperator.h"

class EulerExplicit {
  //double xmin,xmax;
  int cellnumber;

  public:
  EulerExplicit(int &c);
 
  vector<double> ComputeLocalTimeStep(array<double,3>* &field,Euler1D &Euler,const double &CFL,double &dx); // Outputting the local time step for every cell in the domain
  double ComputeGlobalTimeStep(const double &CFL,double &dx,double &lambda_max); // Computing a single time step for all cells in the domain(i.e. min(all local timesteps) 

  

  ~EulerExplicit();


};




#endif
