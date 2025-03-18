// Responsible for Euler Eqs. Computations
#ifndef _EULEROPERATOR_TEST_H_
#define _EULEROPERATOR_TEST_H_
#include "ExactNozzle.h"
#include "MeshGen.h"

using namespace std;

class Euler1D {
  //vector<double> &xcoords;
  double stag_pressure; //stagnation pressure
  double stag_temperature; //stagnation temperature
  double gamma;
  const double Ru = 8314.0; // J/(kmol*K) -- universal gas constant   
  const double MolMass = 28.96; // kg/kmol

  double Density_min = 1.0e-4; //primitive variable sol. limits
  double Velocity_min = 1.0e-7;
  double Pressure_min = 1.0e-4;

  double Density_max = 1.0e6;
  double Velocity_max = 1.0e6;
  double Pressure_max = 1.0e9;
 

  public:
  double dx; //cell thickness -- gets assigned value in Source Term fcn. (P.S. not sure if this is the best way in doing this)
  int interior_cellnum; //holds the cell num (interior only)
  int total_cellnum; //holds the total cell num (including boundary conditions) -- gets assigned value in SetBoundaryConditions fcn.
  double R = Ru / MolMass; //specific gas constant
  

  Euler1D(); //empty constructor for unit testing

  Euler1D(int &cellnum,double &P0,double &T0,double &g); //constructor for Main file

  // Primitive & Conserved variables fcns.
  array<double,3> ComputeConserved(vector<array<double,3>> &Field,int loc);

  void ComputePrimitive(vector<array<double,3>> &Field,array<double,3> &Conserved,int loc);

  // Boundary + Initial Conditions Fcns.
  void SetInitialConditions(vector<array<double,3>> &Field,vector<double> &xcoords); //Complete (tested)
  void SetBoundaryConditions(vector<array<double,3>> &Field,bool &cond); //ADDS ghost cells nodes and computes their values
  void ComputeTotalBoundaryConditions(vector<array<double,3>> &Field,bool &cond); //only computes their values
  void ComputeInflowBoundaryConditions(vector<array<double,3>> &Field); //only computes their values
  void ComputeOutflowBoundaryConditions(vector<array<double,3>> &Field,bool cond); //only computes their values

  // Residual
  void ComputeResidual(vector<array<double,3>> &Resid,vector<array<double,3>> &Field,vector<double> &xcoords,double &dx); //TODO: Computes the residual vector (uses artificial viscosity and dampening

  // Spatial Fluxes Fcns. (including source term)
  array<double,3> ComputeSpatialFlux(vector<array<double,3>> &Field,int loc,int nbor);//TODO
  double ComputeSourceTerm(vector<array<double,3>> &Field,int loc,vector<double> &xcoords);//TODO: May have to look into Pi

  // Artificial Dissipaton Fcns. (using JST Dampening)
  array<double,3> Compute2ndOrderDamping(vector<array<double,3>> &Field,int loc); // viscous term for shocks (c(2))
  array<double,3> Compute4thOrderDamping(vector<array<double,3>> &Field,int loc); // prevents odd-even decoupling (c(4))
  
  double GetEpsilon2(vector<array<double,3>> &Field,int loc); //sensor that detects "sharp" gradients
  double GetEpsilon4(vector<array<double,3>> &Field,int loc); 
 
  double GetLambda(vector<array<double,3>> &Field,int loc);
  double GetNu(vector<array<double,3>> &Field,int loc); //switching fcn.
  double GetMachNumber(vector<array<double,3>> &Field,int loc); //computing Mach Number given known primitve varialbes
  //double ComputeIsentropicMachNumber(vector<array<double,3>> &Field,int loc); //computing Mach Number using isentropic conditions -- only used for
  
  // Supplemental Fcns. (may be used for other fcns. of other classes)
  double GetLambdaMax(vector<array<double,3>> &Field,int loc); //extracts largest eigenvalue for a given cell
  static double GetCellAverageSol(double &A_left,double &A_right,double &dx,array<double,3> &sol_left,array<double,3>&sol_right); //testing x-velocity for now

  ~Euler1D();


};




#endif
