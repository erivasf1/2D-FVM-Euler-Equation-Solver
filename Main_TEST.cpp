//Quasi-1D Nozzle, Euler Eqs. of Converging-to-Diverging Nozzle - Erick Rivas
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdarg.h>

#include "ExactNozzle.h"
#include "MeshGen.h"
#include "EulerOperator_TEST.h"
#include "DataManager_TEST.h"
#include "Output.h"
#include "TimeIntegrator_TEST.h" 

using namespace std;


int main() {

  // Domain and Fluid Parameters
  double xmin = -1.0;
  double xmax = 1.0;
  double stag_pressure = 300.0; //kPa
  double stag_temp = 600.0; //K
  double gamma = 1.4; //specific heat ratio
  double area;
  double area_star = Tools::AreaVal(0.5*(xmin+xmax)); //area at throat (midpoint of xmin and xmax)
  bool cond{false}; //true for subsonic & false for supersonic

  //Mesh Specifications
  int cellnum = 80; //recommending an even number for cell face at the throat of nozzle
  vector<double> xcoords; //!< stores the coords of the cell FACES!!! (i.e. size of xcoords is cellnum+1)!


  //Temporal Specifications
  const int iter_max = 1e2; //max number of iterations
  const int iterout = 1; //number of iterations per solution output
  const double CFL = 0.5; //CFL number (must <= 1 for Euler Explicit integration)
  bool timestep{true}; //true = local time stepping; false = global time stepping

  //Governing Eq. Residuals
  double cont_tol = 1e-3;
  double xmom_tol = 1e-3;
  double energy_tol = 1e-3;

  // ALGORITHM:
  // Create Mesh (verified) -- may have to add ghost cells
  // Set Initial Conditions (verified)
  // Calculate Exact Sol. for comparison (verified)
  // Set Boundary Conditions (verified)
  // Output Initial Residual Norm (verified)
  // Output Initial Solution (verified)
  // Begin Main Loop to iteratively solve Euler Eqs.
  //   Set Time Step
  //   Solve Flow quantity at the new time step
  //   Calc. primitive variables from conserved variables //   Reset boundary conditions (ghost node approach)
  //   Output sol. at every "iterout" iterations
  //   calculate residual norms (normalized by initial value)
  //   check for convergence (if converged, exit main loop)
  // End Main Loop
  // Output solution
  // TODO: Evaluate discreization for error norms for grid convergence (isentropic)

  //double M; //used for debugging M


  //Object Initializations
  vector<array<double,3>> Field(cellnum); //stores the primitive variable sols. at all cells
  vector<array<double,3>> ExactField(cellnum); //stores the exact sol. at all cells
  vector<array<double,3>> Residual(cellnum); //stores the residuals for all interior cells
  vector<array<double,3>> InitResidual(cellnum); //stores the initial residuals for all interior cells


  MeshGen1D Mesh(xmin,xmax,cellnum); //mesh

  //TODO: Only need 1 of these
  // for outputting primitive variables and sol. norms
  SpaceVariables1D Sols; 
  SpaceVariables1D ExactSols; 
  SpaceVariables1D ResidSols; 
  SpaceVariables1D InitResidSols; 

  array<double,3> ResidualNorms; //stores the norms of the residuals

  Tools tool; //used as utilities object

  EulerExplicit Time(cellnum); //for computing time steps

  // Generating Mesh
  Mesh.GenerateMesh(xcoords); //stores all coords in xcoords list
  double dx = abs(xcoords[0]-xcoords[1]); //delta x distance
  Tools::print("dx:%f\n",dx);
  Euler1D Euler(cellnum,stag_pressure,stag_temp,gamma); //for solving Euler eqs.

  // Printing out Temporal Stats
  Tools::print("-Temporal Statistics:\n");
  Tools::print("--CFL: %f\n",CFL);
  Tools::print("Time-stepping method: ");
  (timestep == true) ? Tools::print("Local time-stepping\n") : Tools::print("Global time-stepping\n");  


  //debugging:
  /*array<double,3> init{10,50,100};
  Euler.SetInitialConditions(init,field);
  Tools::print("First,middle,& end pt. of rho: %f,%f,%f\n",field[0][0],field[cellnum/2-1][0],field[cellnum-1][0]);
  Tools::print("First,middle,& end pt. of velocity: %f,%f,%f\n",field[0][1],field[cellnum/2-1][1],field[cellnum-1][1]);
  Tools::print("First,middle,& end pt. of pressure: %f,%f,%f\n",field[0][2],field[cellnum/2-1][2],field[cellnum-1][2]);

  */

  //!!! Solution format: [rho,velocity,pressure]^T

  // SETTING INITIAL CONDITIONS
  //Tools::print("At initial conditions\n");
  //Euler.SetInitialConditions(Field,xcoords);

  //Debug: printing initial conditions
  //const char* filename = "InitialSolutions.txt"; 
  //Sols.OutputPrimitiveVariables(Field,Euler,filename);

  // COMPUTING EXACT SOLUTION -- (should be outputted to a file)
  vector<array<double,3>> ExactSol_Faces(cellnum+1);
  //Tools::print("Exact Solution Output\n");
  //Computing and storing exact sol. at face
  for (int n=0;n<(int)ExactSol_Faces.size();n++) {
    area = tool.AreaVal(xcoords[n]);
    cond = (xcoords[n] < 0) ? true:false; 
    SuperSonicNozzle Nozzle(area,area_star,stag_pressure,stag_temp,cond);
    Nozzle.ComputeExactSol(ExactSol_Faces[n]); 
 
    //Tools::print("Point %f\n",xcoords[i]);
    //Tools::print("Density,Velocity,& Pressure: %f,%f,%f\n",field[i][0],field[i][1],field[i][2]);
  }

  //Computing cell-average sol. for all cells
  ExactSols.ComputeCellAveragedSol(ExactSol_Faces,ExactField,xcoords,dx);

  //TODO: Temporarily set initial conditions to exact solutions
  Field = ExactField;
  //Debug: printing initial conditions w/ no BCs
  const char* filename = "InitialSolutions.txt"; 
  Sols.OutputPrimitiveVariables(Field,Euler,filename);
  

  // SETTING BOUNDARY CONDITIONS
  // Note: Was found that extrapolating values at boundary gave a negative pressure
  Euler.SetBoundaryConditions(Field,cond);

  //for (int i=0;i<(int)Field.size();i++) //!< Applying sol. limiter for every cell
  Time.SolutionLimiter(Field);

  

  const char* filename2 = "InitSolutionswBCs.txt"; 
  Sols.OutputPrimitiveVariables(Field,Euler,filename2);

  //debug
  /*
  array<double,3> init{10,50,100};
  Tools::print("Size of Field before set BC: %d\n",Field.size());
  Euler.SetBoundaryConditions(Field,init);
  Tools::print("Size of Field after set BC: %d\n",Field.size());
  */

  // COMPUTING INITIAL RESIDUAL NORMS
  // using ResidSols spacevariable
  array<double,3> InitNorms;
  Euler.ComputeResidual(InitResidual,Field,xcoords,dx); //computing residuals per cell
  InitNorms = InitResidSols.ComputeSolutionNorms(InitResidual); //computing L2 norm of residuals
  Tools::print("-Initial Residual Norms\n");
  Tools::print("--Continuity:%e\n",InitNorms[0]);
  Tools::print("--X-Momentum:%e\n",InitNorms[1]);
  Tools::print("--Energy:%e\n",InitNorms[2]);
  
  const char* filename3 = "InitialLocalResiduals.txt";
  InitResidSols.OutputLocalResiduals(InitResidual,filename3);


  // BEGIN OF MAIN LOOP
  vector<double> time_steps;
  Residual = InitResidual;


  string it,name,resid; //used for outputting file name
  int iter; //iteration number

  //Opening File that stores residuals
  ofstream myresids;
  myresids.open("SolResids.txt");
  myresids<<"Iteration"<<"  "<<"Contintuity"<<"  "<<"X-Momentum"<<"  "<<"Energy"<<endl;

  for (iter=1;iter<iter_max;iter++){

    //debugging only
    Tools::print("------Iteration #: %d----------\n",iter);
  
    //COMPUTE TIME STEP
    // if global time step, chosen then create a vector<double> of the smallest time step
    time_steps = Time.ComputeLocalTimeStep(Field,Euler,CFL,dx);

    //COMPUTE NEW SOL. VALUES 
    Time.FWDEulerAdvance(Field,Residual,Euler,time_steps,xcoords,dx);
    Time.SolutionLimiter(Field); //applies solution limiter to all cells (including ghost cells)

    //COMPUTE RESIDUAL NORMS
    //ResidSols.ComputeSolutionNorms(Residual);

    //COMPUTE BOUNDARY CONDITIONS
    Euler.ComputeTotalBoundaryConditions(Field,cond);
    Time.SolutionLimiter(Field); //temporarily reapplying the limiter

    //OUTPUT SOL. IN TEXT FILE EVERY "ITEROUT" STEPS
    if (iter % iterout == 0) {
      it = to_string(iter);
      name = "SolResults/Iteration";
      name += it;
      name += ".txt";
      const char* filename_iter = name.c_str(); 
      Sols.OutputPrimitiveVariables(Field,Euler,filename_iter);
      
    }

    //COMPUTE NEW RESIDUALS 
    Euler.ComputeResidual(Residual,Field,xcoords,dx); //computing residuals per cell


    //COMPUTE RESIDUAL NORMS & CHECK FOR CONVERGENCE
    ResidualNorms = ResidSols.ComputeSolutionNorms(Residual);
    Tools::print("Norms\n");
    Tools::print("Continuity:%e\nX-Momentum:%e\nEnergy:%e\n",ResidualNorms[0],ResidualNorms[1],ResidualNorms[2]);
    myresids<<iter<<"  "<<ResidualNorms[0]<<"  "<<ResidualNorms[1]<<"  "<<ResidualNorms[2]<<endl;

    if (ResidualNorms[0]/InitNorms[0] <= cont_tol || ResidualNorms[1]/InitNorms[1] <= xmom_tol || ResidualNorms[2]/InitNorms[2] <= energy_tol)
      break;
    
    //debug:
    //resid = Residual.data();
    //Tools::print("(Before)1st resid of continuity:%e\n",resid[0][0]);
    //Tools::print("(Before)1st Resid of continuity:%e\n",Residual[0][0]);
    //Tools::print("(After)1st resid of continuity:%e\n",resid[0][0]);

    //Compute new Boundary Conditions
    //Output sol. every "iterout" steps

    //Cal. residual norms and see if it is within iterative convergence tolerance

  }

  //Final Output of Solution

  if (iter==iter_max)
    Tools::print("Failed to converge!\n");

  else {
    const char* filename_final = "ConvergedSolution.txt" ; 
    Sols.OutputPrimitiveVariables(Field,Euler,filename_final);
  }

  //Closing Residuals file
  myresids.close();

  //Evaluate discretization norms for grid convergence

  return 0;
}
