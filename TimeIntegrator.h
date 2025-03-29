// Responsible for creating a Mesh (1D for now)
#ifndef _TIMEINTEGRATOR_H_
#define _TIMEINTEGRATOR_H_
#include "ExactNozzle.h"
#include "EulerOperator.h"

class EulerExplicit {
  //double xmin,xmax;
  int cellnumber;
  double Density_min = 1.0e-7;
  double Velocity_min = 1.0e-7;
  double Pressure_min = 1.0e-7;

  double Density_max = 1.0e6;
  double Velocity_max = 1.0e6;
  double Pressure_max = 1.0e9;

  public:
  EulerExplicit(int &c);
 
  vector<double> ComputeLocalTimeStep(vector<array<double,3>>* &field,Euler1D* &euler,const double &CFL,double &dx); // Outputting the local time step for every cell in the domain
  vector<double> ComputeGlobalTimeStep(vector<array<double,3>>* &field,Euler1D* &euler,const double &CFL,double &dx); // Outputting the local time step for every cell in the domain

  void FWDEulerAdvance(vector<array<double,3>>* &field,vector<array<double,3>>* &resid,vector<double> &time_steps,vector<double> &xcoords,double &dx); //Computing the new solution at the next time step

  void SolutionLimiter(vector<array<double,3>>* &field); //Re-assigning primitive variables to specified max and min limits

  void UnderRelaxationCheck(array<double,3> ResidPrevNorm,array<double,3> ResidNorm,double C,array<bool,3> &check); //looks for residual norm increase to mark needed under-relaxation for each equation 

  ~EulerExplicit();


};




#endif
