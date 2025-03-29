// Responsible for creating a Mesh (1D for now)
#ifndef _TIMEINTEGRATOR_TEST_H_
#define _TIMEINTEGRATOR_TEST_H_
#include "ExactNozzle.h"
#include "EulerOperator_TEST.h"

//1st order accurate time step
class EulerExplicit {
  //double xmin,xmax;
  int cellnumber; // cellnumber in domain
  double Density_min = 1.0e-7;
  double Velocity_min = 1.0e-7;
  double Pressure_min = 1.0e-7;
  //double Pressure_min = 1.0e-7;

  double Density_max = 1.0e6;
  double Velocity_max = 1.0e6;
  double Pressure_max = 1.0e9;

  public:
  EulerExplicit(int &c);
 
  vector<double> ComputeLocalTimeStep(vector<array<double,3>> &Field,Euler1D &Euler,const double &CFL,double &dx); // Outputting the local time step for every cell in the domain
  vector<double> ComputeGlobalTimeStep(vector<array<double,3>> &Field,Euler1D &Euler,const double &CFL,double &dx); // Computing a single time step for all cells in the domain(i.e. min(all local timesteps) 


  void FWDEulerAdvance(vector<array<double,3>> &Field,vector<array<double,3>> &Resid,Euler1D &Euler,vector<double> &time_steps,vector<double> &xcoords,double &dx,array<double,3> &Omega); //Computing the new solution at the next time step

  void SolutionLimiter(vector<array<double,3>> &Sol); //Re-assigning primitive variables to specified max and min limits

  void UnderRelaxationCheck(array<double,3> ResidPrevNorm,array<double,3> ResidNorm,double C,array<bool,3> &check); //looks for residual norm increase to mark needed under-relaxation for each equation
  

  ~EulerExplicit();


};




#endif
