// Responsible for creating a Mesh (1D for now)
#ifndef _EULEROPERATOR_H_
#define _EULEROPERATOR_H_
#include "ExactNozzle.h"
#include "MeshGen.h"

class Euler1D {
  vector<double> &xcoords;
  double dx; //cell thickness

  public:
  Euler1D(vector<double> &coords);

  // Boundary Conditions Fcns.
  void SetInitialConditions(array<double,3> &init_val,array<double,3>* &field); //Complete (tested)
  void SetBoundaryConditions(vector<array<double,3>> &Field,array<double,3> &init); 
  void ComputeBoundaryConditions(vector<array<double,3>> &Field,array<double,3> &init); 

  // Residual
  void ComputeResidual(); //TODO: Computes the residual vector (uses artificial viscosity and dampening

  // Spatial Fluxes Fcns. (including source term)
  array<double,3> ComputeSpatialFlux(array<double,3>* &field,int &loc,int &nbor);//TODO
  double ComputeSourceTerm(array<double,3>* &field,int &loc);//TODO: May have to look into Pi

  // Artificial Dissipaton Fcns. (using JST Dampening)
  array<double,3> Compute2ndOrderDamping(array<double,3>* &field,int loc); // viscous term for shocks (c(2))
  array<double,3> Compute4thOrderDamping(array<double,3>* &field,int loc); // prevents odd-even decoupling (c(4))
  
  double GetEpsilon2(array<double,3>* &field,int &loc); //sensor that detects "sharp" gradients
  double GetEpsilon4(array<double,3>* &field,int &loc); 
 
  double GetLambda(array<double,3>* &field,int &loc);
  double GetNu(array<double,3>* &field,int loc); //switching fcn.
  

  ~Euler1D();


};




#endif
