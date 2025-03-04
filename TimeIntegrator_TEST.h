// Responsible for creating a Mesh (1D for now)
#ifndef _TIMEINTEGRATOR_TEST_H_
#define _TIMEINTEGRATOR_TEST_H_
#include "ExactNozzle.h"
#include "EulerOperator_TEST.h"

class EulerExplicit {
  //double xmin,xmax;
  int cellnumber;

  public:
  EulerExplicit(int &c);
 
  vector<double> ComputeLocalTimeStep(vector<array<double,3>> &Field,Euler1D &Euler,const double &CFL,double &dx); // Outputting the local time step for every cell in the domain
  double ComputeGlobalTimeStep(const double &CFL,double &dx,double &lambda_max); // Computing a single time step for all cells in the domain(i.e. min(all local timesteps) 

  void FWDEulerAdvance(vector<array<double,3>> &Field,vector<array<double,3>> &Resid,vector<double> &time_steps,vector<double> &xcoords,double &dx); //Computing the new solution at the next time step

  

  ~EulerExplicit();


};




#endif
