// Responsible for Euler Eqs. Computations
#ifndef _EULEROPERATOR_H_
#define _EULEROPERATOR_H_
#include "ExactNozzle.h"
#include "MeshGen.h"

class Euler1D {
  //vector<double> &xcoords;
  double stag_pressure; //stagnation pressure
  double back_pressure; //stagnation pressure
  double stag_temperature; //stagnation temperature
  double gamma;
  const double Ru = 8314.0; // J/(kmol*K) -- universal gas constant   
  const double MolMass = 28.96; // kg/kmol
 
  double Density_max = 1.0e9;
  double Velocity_max = 1.0e9;
  double Pressure_max = 1.0e9;

  public:
  double dx; //cell thickness
  int interior_cellnum; //holds the cell num (interior only)
  int total_cellnum; //holds the total cell num (including boundary conditions)
  double R = Ru / MolMass; //specific gas constant

  Euler1D(); //empty constructor for unit testing

  Euler1D(int &cellnum,double &P0,double &BP,double &T0,double &g); //constructor for Main file

  // Primitive & Conserved variables fcns.
  array<double,3> ComputeConserved(vector<array<double,3>>* &field,int loc);

  void ComputePrimitive(vector<array<double,3>>* &field,array<double,3> &Conserved,int loc);

  // Boundary + Initial Conditions Fcns.
  void SetInitialConditions(vector<array<double,3>>* &field,vector<double> &xcoords); //Complete (tested)
  void SetBoundaryConditions(vector<array<double,3>>* &field,bool &cond); //ADDS ghost cells nodes and computes their values
  void ComputeTotalBoundaryConditions(vector<array<double,3>>* &field,bool &cond); //only computes their values
  void ComputeInflowBoundaryConditions(vector<array<double,3>>* &field); //only computes their values
  void ComputeOutflowBoundaryConditions(vector<array<double,3>>* &field,bool& cond); //only computes their values

  // Residual Fcns.
  void ComputeResidual(vector<array<double,3>>* &resid,vector<array<double,3>>* &field,vector<double> &xcoords,double &dx);

  // Spatial Fluxes Fcns. (including source term)
  // Central Difference using Central quadrature fcn.
  array<double,3> ComputeSpatialFlux_BASE(vector<array<double,3>>* &field,int loc,int nbor);
  //Upwind Schemes
  array<double,3> ComputeSpatialFlux_UPWIND1stOrder(vector<array<double,3>>* &field,bool &method,int loc,int rnbor); //1st order upwind schemes
  array<double,3> ComputeSpatialFlux_UPWIND2ndOrder(bool &method); //2nd order upwind schemes
  //VanLeer Fcns.
  array<double,3> VanLeerCompute(vector<array<double,3>>* &field,int loc,bool sign); //returns convective+pressure flux of specified state (either right or left)
  double GetC(double M,bool sign); //c value
  double GetAlpha(double M,bool sign); //alpha value
  double GetBeta(double M); //beta value
  double GetVanLeerM(double M,bool sign); //Van Leer MachNumber
  double GetD(double M,bool sign); //D value
  double GetP2Bar(double M,bool sign); //Pressure double bar
  // TODO: Add Upwind Schemes here: Van Leer and Roe's Method
  //Upwind fcn. -- this will sum right and left states
  //VanLeerCompute -- evaluates either specified left or right state via Van Leer's Method
  //RoeCompute -- evaluates either specified left or right state via Roe's Method

  double ComputeSourceTerm(vector<array<double,3>>* &field,int loc,vector<double> &xcoords);

  // Artificial Dissipaton Fcns. (using JST Dampening)
  array<double,3> Compute2ndOrderDamping(vector<array<double,3>>* &field,int loc); // viscous term for shocks (c(2))
  array<double,3> Compute4thOrderDamping(vector<array<double,3>>* &field,int loc); // prevents odd-even decoupling (c(4))
  
  double GetEpsilon2(vector<array<double,3>>* &field,int loc); //sensor that detects "sharp" gradients
  double GetEpsilon4(vector<array<double,3>>* &field,int loc); 
 
  double GetLambda(vector<array<double,3>>* &field,int loc);
  double GetNu(vector<array<double,3>>* &field,int loc); //switching fcn.
  double GetMachNumber(vector<array<double,3>>* &field,int loc); //used for GetLambda fcn.
  
  // Supplemental Fcns. (may be used for other fcns. of other classes)
  double GetLambdaMax(vector<array<double,3>>* &field,int loc); //extracts largest eigenvalue for a given cell
  static double GetCellAverageSol(double &A_left,double &A_right,double &dx,array<double,3> &sol_left,array<double,3> &sol_right); //testing x-velocity for now

  ~Euler1D();


};




#endif
