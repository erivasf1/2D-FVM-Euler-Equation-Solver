// Responsible for creating a Mesh (1D for now)
#ifndef _TIMEINTEGRATOR_TEST_H_
#define _TIMEINTEGRATOR_TEST_H_
#include "ExactNozzle.h"
#include "EulerOperator_TEST.h"

class EulerExplicit {
  //double xmin,xmax;
  int cellnumber; // cellnumber in domain
  double Density_min = 1.0e-8;
  double Velocity_min = 1.0e-8;
  double Pressure_min = 2.0;

  double Density_max = 1.0e8;
  double Velocity_max = 1.0e8;
  double Pressure_max = 1.0e8;

  public:
  EulerExplicit(int &c);
 
  vector<double> ComputeLocalTimeStep(vector<array<double,3>> &Field,Euler1D &Euler,const double &CFL,double &dx); // Outputting the local time step for every cell in the domain
  double ComputeGlobalTimeStep(const double &CFL,double &dx,double &lambda_max); // Computing a single time step for all cells in the domain(i.e. min(all local timesteps) 


  void FWDEulerAdvance(vector<array<double,3>> &Field,vector<array<double,3>> &Resid,Euler1D &Euler,vector<double> &time_steps,vector<double> &xcoords,double &dx); //Computing the new solution at the next time step

  void SolutionLimiter(array<double,3> &Sol); //Re-assigning primitive variables to specified max and min limits
  

  ~EulerExplicit();


};




#endif
