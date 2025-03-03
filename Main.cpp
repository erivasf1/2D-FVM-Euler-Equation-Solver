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
#include "EulerOperator.h"
#include "DataManager.h"
#include "Output.h"
#include "TimeIntegrator.h"

using namespace std;


int main() {

  // Initializing Parameters
  double xmin = -1.0;
  double xmax = 1.0;
  double stag_pressure = 300.0; //kPa
  double stag_temp = 600.0; //K
  double gamma = 1.4; //specific heat ratio
  int pt_num = 5; //# of evenly-spaced requested points (including xmin and xmax)
  double area;
  double area_star; //area at throat
  bool cond{false}; //true for subsonic & false for supersonic

  //Mesh Specifications
  int cellnum = 8; //recommending an even number for cell face at the throat of nozzle
  vector<double> xcoords; //!< stores the coords of the cell FACES!!! (i.e. size of xcoords is cellnum+1)!


  //Temporal Specifications
  const int iter_max = 1000;
  const int iterout = 5; //number of iterations per solution output
  const double CFL = 0.8; //CFL number (must <= 1 for Euler Explicit integration)
  bool timestep{true}; //true = local time stepping; false = global time stepping

  // ALGORITHM:
  // Create Mesh (verified) -- may have to add ghost cells
  // Set Initial Conditions (verified)
  // Calculate Exact Sol. for comparison (verified)
  // Set Boundary Conditions (verified)
  // TODO: Output Initial Residual Norm
  // TODO: Output Initial Solution
  // TODO: Begin Main Loop to iteratively solve Euler Eqs.
  //   Set Time Step
  //   Solve Flow quantity at the new time step
  //   Calc. primitive variables from conserved variables //   Reset boundary conditions (ghost node approach)
  //   Output sol. at every "iterout" iterations
  //   calculate residual norms (normalized by initial value)
  //   check for convergence (if converged, exit main loop)
  // End Main Loop
  // TODO: Output solution
  // TODO: Evaluate discreization for error norms for grid convergence (isentropic)

  //double M; //used for debugging M


  //Object Initializations
  vector<array<double,3>> Field(cellnum);
  vector<array<double,3>> ExactField(cellnum);
  vector<array<double,3>> Residual(cellnum);
  vector<array<double,3>> InitResidual(cellnum);

  array<double,3>* field; //pointer to Field solutions
  array<double,3>* exact_sols; //pointer to exact solution field values
  array<double,3>* resid; //pointer to residual field values per cell
  array<double,3>* init_resid; //pointer to residual field values per cell

  MeshGen1D Mesh(xmin,xmax,cellnum); //mesh

  SpaceVariables1D Sols(cellnum,Field,field); //for storing solutions
  SpaceVariables1D ExactSols(cellnum,ExactField,exact_sols); //for storing exact solutions
  SpaceVariables1D ResidSols(cellnum,Residual,resid); //for storing residuals for every cell
  SpaceVariables1D InitResidSols(cellnum,InitResidual,init_resid); //for storing initial residuals for every cell

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
  Euler.SetInitialConditions(field,xcoords);

  // COMPUTING EXACT SOLUTION -- (should be outputted to a file)
  array<double,3> sol;
  area_star = tool.AreaVal(0.0); //area at throat
  //Tools::print("Exact Solution Output\n");
  /*for (int i=0;i<cellnum;i++) {
    area = tool.AreaVal(xcoords[i]);
    cond = (xcoords[i] < 0) ? true:false; 
    SuperSonicNozzle Nozzle(area,area_star,stag_pressure,stag_temp,cond);
    Nozzle.ComputeExactSol(sol);
 
    exact_sols[i] = sol; //storing solution values to exact sol.
    
    //Tools::print("Point %f\n",xcoords[i]);
    //Tools::print("Density,Velocity,& Pressure: %f,%f,%f\n",field[i][0],field[i][1],field[i][2]);

  }*/
  //area = tool.AreaVal(xcoord[i]);

  // SETTING BOUNDARY CONDITIONS
  Euler.SetBoundaryConditions(Field,field,cond);

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
  Euler.ComputeResidual(init_resid,field,xcoords,dx); //computing residuals per cell
  InitNorms = InitResidSols.ComputeSolutionNorms(init_resid); //computing L2 norm of residuals
  Tools::print("-Initial Residual Norms\n");
  Tools::print("--Continuity:%e\n",InitNorms[0]);
  Tools::print("--X-Momentum:%e\n",InitNorms[1]);
  Tools::print("--Energy:%e\n",InitNorms[2]);


  // BEGIN OF MAIN LOOP
  vector<double> time_steps;
  resid = init_resid; //assining residuals to initial residuals

  /*time_steps = Time.ComputeLocalTimeStep(field,Euler,CFL,dx);
  Tools::print("Local time step list:\n");
  for(int i=0;i<cellnum;i++)
    Tools::print("cell index:%d &time step: %f\n",i,time_steps[i]);
  */
  for (int iter=1;iter<iter_max;iter++){

    //debugging only
    Tools::print("Iteration #: %d\n",iter);
  
    //COMPUTE TIME STEP
    // if global time step, chosen then create a vector<double> of the smallest time step
    time_steps = Time.ComputeLocalTimeStep(field,Euler,CFL,dx);

    //COMPUTE NEW SOL. VALUES 
    Time.FWDEulerAdvance(field,resid,time_steps,xcoords,dx);
    Euler.ComputeResidual(resid,field,xcoords,dx); //computing residuals per cell
    ResidSols.ComputeSolutionNorms(resid);

    //Compute new Boundary Conditions
    //Output sol. every "iterout" steps

    //Cal. residual norms and see if it is within iterative convergence tolerance

  }

  //Output solution
  //Evaluate discretization norms for grid convergence

  return 0;
}
