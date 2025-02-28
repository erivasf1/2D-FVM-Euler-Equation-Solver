// Responsible for Euler Eqs. Computations
#ifndef _EULEROPERATOR_H_
#define _EULEROPERATOR_H_
#include "ExactNozzle.h"
#include "MeshGen.h"

class Euler1D {
  vector<double> &xcoords;
  double dx; //cell thickness
  double stag_pressure; //stagnation pressure
  double stag_temperature; //stagnation temperature
  double gamma;
  const double Ru = 8314.0; // J/(kmol*K) -- universal gas constant   
  const double MolMass = 28.96; // kg/kmol
  double R = Ru / MolMass; //specific gas constant

  static Euler1D* instance; //static pointer to hold instance so static function can use non-static functions

  public:
  Euler1D(vector<double> &coords,double &P0,double &T0,double &g);

  // Boundary + Initial Conditions Fcns.
  void SetInitialConditions(array<double,3>* &field); //Complete (tested)
  void SetBoundaryConditions(vector<array<double,3>> &Field,array<double,3> &init); 
  void ComputeBoundaryConditions(vector<array<double,3>> &Field,array<double,3> &init); 

  // Residual
  void ComputeResidual(); //TODO: Computes the residual vector (uses artificial viscosity and dampening

  // Spatial Fluxes Fcns. (including source term)
  array<double,3> ComputeSpatialFlux(array<double,3>* &field,int &loc,int &nbor);//TODO
  double ComputeSourceTerm(array<double,3>* &field,int &loc);//TODO: May have to look into Pi

  // Artificial Dissipaton Fcns. (using JST Dampening)
  static array<double,3> Compute2ndOrderDamping(array<double,3>* &field,int loc); // viscous term for shocks (c(2))
  array<double,3> Compute4thOrderDamping(array<double,3>* &field,int loc); // prevents odd-even decoupling (c(4))
  
  double GetEpsilon2(array<double,3>* &field,int &loc); //sensor that detects "sharp" gradients
  double GetEpsilon4(array<double,3>* &field,int &loc); 
 
  double GetLambda(array<double,3>* &field,int &loc);
  double GetNu(array<double,3>* &field,int loc); //switching fcn.
  double GetMachNumber(array<double,3>* field,int &loc);
  

  ~Euler1D();


};




#endif
