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

using namespace std;

Euler1D* Euler1D::instance = nullptr; //used to assign the static pointer member variable, so that static functions can access non static functions

int main() {

  // Initializing Parameters
  Tools tool;
  double xmin = -1.0;
  double xmax = 1.0;
  double stag_pressure = 300.0; //kPa
  double stag_temp = 600.0; //K
  double gamma = 1.4; //specific heat ratio
  int pt_num = 5; //# of evenly-spaced requested points (including xmin and xmax)
  double area;
  double area_star = tool.AreaVal(0.0); //area at throat
  bool cond; //true for subsonic & false for supersonic

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

  //Mesh Specifications
  int cellnum = 8; //recommending an even number for cell face at the throat of nozzle
  vector<double> xcoords;

  //Object Initializations
  vector<array<double,3>> Field(cellnum);
  vector<array<double,3>> ExactField(cellnum);
  array<double,3>* field; //pointer to Field solutions
  array<double,3>* exact_sols; //pointer to exact solution field values
  MeshGen1D Mesh(xmin,xmax,cellnum); //mesh
  Euler1D Euler(xcoords,stag_pressure,stag_temp,gamma); //for solving Euler eqs.
  SpaceVariables1D Sols(cellnum,Field,field); //for storing solutions
  SpaceVariables1D ExactSols(cellnum,ExactField,exact_sols); //for storing exact solutions

  
  // Generating Mesh
  Mesh.GenerateMesh(xcoords); //stores all coords in xcoords list

  //debugging:
  /*array<double,3> init{10,50,100};
  Euler.SetInitialConditions(init,field);
  Tools::print("First,middle,& end pt. of rho: %f,%f,%f\n",field[0][0],field[cellnum/2-1][0],field[cellnum-1][0]);
  Tools::print("First,middle,& end pt. of velocity: %f,%f,%f\n",field[0][1],field[cellnum/2-1][1],field[cellnum-1][1]);
  Tools::print("First,middle,& end pt. of pressure: %f,%f,%f\n",field[0][2],field[cellnum/2-1][2],field[cellnum-1][2]);

  */

  //!!! Solution format: [rho,velocity,pressure]^T
  // Setting Initial Conditions
  Euler.SetInitialConditions(field);

  // Computing Exact Solution
  array<double,3> sol;
  for (int i=0;i<cellnum;i++) {
    area = tool.AreaVal(xcoords[i]);
    cond = (xcoords[i] < 0) ? true:false; 
    SuperSonicNozzle Nozzle(area,area_star,stag_pressure,stag_temp,cond);
    Nozzle.ComputeExactSol(sol);
 
    exact_sols[i] = sol; //storing solution values to exact sol.
    
    Tools::print("Point %f\n",xcoords[i]);
    Tools::print("Density,Velocity,& Pressure: %f,%f,%f\n",field[i][0],field[i][1],field[i][2]);

  }
  //area = tool.AreaVal(xcoord[i]);

  //TODO: Setting Boundary Conditions
  //Euler.SetBoundaryConditions(Field,init);

  //debug
  /*
  array<double,3> init{10,50,100};
  Tools::print("Size of Field before set BC: %d\n",Field.size());
  Euler.SetBoundaryConditions(Field,init);
  Tools::print("Size of Field after set BC: %d\n",Field.size());
  */

  // Printing Initial Residual norms











  return 0;
}
