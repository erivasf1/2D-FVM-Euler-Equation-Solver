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

  // BOUNDARY + INITIAL CONDITIONS FCNS.
  void SetInitialConditions(vector<array<double,3>>* &field,vector<double> &xcoords); //Complete (tested)
  void SetBoundaryConditions(vector<array<double,3>>* &field,bool &cond); //ADDS ghost cells nodes and computes their values
  void ComputeTotalBoundaryConditions(vector<array<double,3>>* &field,bool &cond); //only computes their values
  void ComputeInflowBoundaryConditions(vector<array<double,3>>* &field); //only computes their values
  void ComputeOutflowBoundaryConditions(vector<array<double,3>>* &field,bool& cond); //only computes their values

  // RESIDUAL FCNS.
  void ComputeResidual(vector<array<double,3>>* &resid,vector<array<double,3>>* &field,vector<double> &xcoords,double &dx,bool flux_scheme,bool flux_accuracy,bool upwind_scheme,double epsilon);

  // SPATIAL FLUXES FCNS. (INCLUDING SOURCE TERM)
  // Central Difference using Central quadrature fcn.
  array<double,3> ComputeSpatialFlux_CELL(array<double,3> &field_state); //at cell
  array<double,3> ComputeSpatialFlux_BASE(vector<array<double,3>>* &field,int loc,int nbor); //at cell interface
  // Upwind Schemes
  array<double,3> ComputeSpatialFlux_UPWIND1stOrder(vector<array<double,3>>* &field,bool method,int loc,int nbor); //1st order upwind schemes
  array<double,3> ComputeSpatialFlux_UPWIND2ndOrder(vector<array<double,3>>* &field,bool method,int loc, int nbor,double epsilon); //2nd order upwind schemes
  //VanLeer Fcns.
  array<double,3> VanLeerCompute(array<double,3> &field_state,bool sign); //returns convective+pressure flux of specified state (either right or left)
  double GetC(double M,bool sign); //c value
  double GetAlpha(double M,bool sign); //alpha value
  double GetBeta(double M); //beta value
  double GetVanLeerM(double M,bool sign); //Van Leer MachNumber
  double GetD(double M,bool sign); //D value
  double GetP2Bar(double M,bool sign); //Pressure double bar
  //Roe Fcns.
  array<double,3> ComputeRoeFlux(array<double,3> &field_ltstate,array<double,3> &field_rtstate); //rho-avg eigenvalues
  array<double,3> ComputeRoeWaveAmps(array<double,3> &roe_vars,array<double,3> &field_ltstate,array<double,3> &field_rtstate,double abar); //rho-avg eigenvalues
  array<double,3> ComputeRoeEigenVals(array<double,3> &rho_vars,double abar); //rho-avg eigenvalues
  array<array<double,3>,3> ComputeRoeEigenVecs(array<double,3> &roe_vars,double abar); //rho-avg eigenvectors
  array<double,3> ComputeRoeAvgVars(array<double,3> &field_ltstate,array<double,3> &field_rtstate,double &abar); //rho-avg. vars
  double ComputeSourceTerm(vector<array<double,3>>* &field,int loc,vector<double> &xcoords);

  // MUSCL extrapolation + Flux Limiters
  array<array<double,3>,2> MUSCLApprox(vector<array<double,3>>* &field,int loc,int nbor,double epsilon); //outputs the left and right state primitive variables
  array<array<double,3>,2> ComputeBetaLimiter(vector<array<double,3>>* &field,int loc,int nbor,int r_nbor,int l_nbor,double beta); //computes limiter using the beta limiter method
  array<double,3> ComputeVanLeerLimiter(vector<array<double,3>>* &field,array<double,3> &r_vec); //computes limiter using the Van Leer method
  array<array<double,3>,2> ComputeRVariation(vector<array<double,3>>* &field,int loc,int nbor,int r_nbor,int l_nbor); //consecutive variation

  array<double,3> ComputeRPlusVariation(vector<array<double,3>>* &field,int loc,int r_nbor,int nbor); //plus part of consecutive variation
  array<double,3> ComputeRMinusVariation(vector<array<double,3>>* &field,int loc,int l_nbor,int nbor); //minus part of consecutive variation

  void FreezeLimiter(bool freeze,double resid_current,double resid_prev,int freeze_count);

  // ARTIFICIAL DISSIPATON FCNS. (USING JST DAMPENING)
  array<double,3> Compute2ndOrderDamping(vector<array<double,3>>* &field,int loc); // viscous term for shocks (c(2))
  array<double,3> Compute4thOrderDamping(vector<array<double,3>>* &field,int loc); // prevents odd-even decoupling (c(4))
  
  double GetEpsilon2(vector<array<double,3>>* &field,int loc); //sensor that detects "sharp" gradients
  double GetEpsilon4(vector<array<double,3>>* &field,int loc); 
 
  double GetLambda(vector<array<double,3>>* &field,int loc);
  double GetNu(vector<array<double,3>>* &field,int loc); //switching fcn.
  double GetMachNumber(vector<array<double,3>>* &field,int loc); //used for GetLambda fcn.
  
  // SUPPLEMENTAL FCNS. (MAY BE USED FOR OTHER FCNS. OF OTHER CLASSES)
  double GetLambdaMax(vector<array<double,3>>* &field,int loc); //extracts largest eigenvalue for a given cell
  static double GetCellAverageSol(double &A_left,double &A_right,double &dx,array<double,3> &sol_left,array<double,3> &sol_right); //testing x-velocity for now
  double ComputeMachNumber(array<double,3> &sols); //computes Mach number for solution vars. that are not part of the field -- same formula as GetMachNumber

  ~Euler1D();


};




#endif
