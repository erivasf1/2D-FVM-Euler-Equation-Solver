// Responsible for Euler Eqs. Computations
#ifndef _EULEROPERATOR_H_
#define _EULEROPERATOR_H_
#include "ExactNozzle.h"
#include "MeshGen.h"

class Euler1D {
  //vector<double> &xcoords;
  double dx; //cell thickness
  double stag_pressure; //stagnation pressure
  double stag_temperature; //stagnation temperature
  double gamma;
  const double Ru = 8314.0; // J/(kmol*K) -- universal gas constant   
  const double MolMass = 28.96; // kg/kmol
  double R = Ru / MolMass; //specific gas constant
 

  public:
  int interior_cellnum; //holds the cell num (interior only)
  int total_cellnum; //holds the total cell num (including boundary conditions)

  Euler1D(); //empty constructor for unit testing

  Euler1D(int &cellnum,double &P0,double &T0,double &g); //constructor for Main file

  // Boundary + Initial Conditions Fcns.
  void SetInitialConditions(array<double,3>* &field,vector<double> &xcoords); //Complete (tested)
  void SetBoundaryConditions(vector<array<double,3>> &Field,array<double,3>* &field,bool &cond); //ADDS ghost cells nodes and computes their values
  void ComputeTotalBoundaryConditions(array<double,3>* &field,bool &cond); //only computes their values
  void ComputeInflowBoundaryConditions(array<double,3>* &field); //only computes their values
  void ComputeOutflowBoundaryConditions(array<double,3>* &field,bool& cond); //only computes their values

  // Residual
  void ComputeResidual(array<double,3>* &resid,array<double,3>* &field,vector<double> &xcoords,double &dx); //TODO: Computes the residual vector (uses artificial viscosity and dampening

  // Spatial Fluxes Fcns. (including source term)
  array<double,3> ComputeSpatialFlux(array<double,3>* &field,int loc,int nbor);//TODO
  double ComputeSourceTerm(array<double,3>* &field,int &loc,vector<double> &xcoords);//TODO: May have to look into Pi

  // Artificial Dissipaton Fcns. (using JST Dampening)
  array<double,3> Compute2ndOrderDamping(array<double,3>* &field,int loc); // viscous term for shocks (c(2))
  array<double,3> Compute4thOrderDamping(array<double,3>* &field,int loc); // prevents odd-even decoupling (c(4))
  
  double GetEpsilon2(array<double,3>* &field,int &loc); //sensor that detects "sharp" gradients
  double GetEpsilon4(array<double,3>* &field,int &loc); 
 
  double GetLambda(array<double,3>* &field,int &loc);
  double GetNu(array<double,3>* &field,int loc); //switching fcn.
  double GetMachNumber(array<double,3>* field,int loc); //used for GetLambda fcn.
  
  // Supplemental Fcns. (may be used for other fcns. of other classes)
  double GetLambdaMax(array<double,3>* &field,int &loc); //extracts largest eigenvalue for a given cell

  ~Euler1D();


};




#endif
