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

using namespace std;

int main() {

  // Initializing Parameters
  Tools tool;
  double xmin = -1.0;
  double xmax = 1.0;
  double stag_pressure = 300.0; //kPa
  double stag_temp = 600.0; //K
  int pt_num = 5; //# of evenly-spaced requested points (including xmin and xmax)
  double area;
  double area_star = tool.AreaVal(0.0); //area at throat
  bool cond; //true for subsonic & false for supersonic

  // ALGORITHM:
  // Create Mesh (verified) -- may have to add ghost cells
  // TODO: Set Initial Conditions (In progress)
  // TODO: Calculate Exact Sol. for comparison  
  // TODO: Set Boundary Conditions (ghost node approach) -- just setting values
  // TODO: Output Initial Residual Norm
  // TODO: Output Initial Solution
  // TODO: Begin Main Loop to iteratively solve Euler Eqs.
  //   Set Time Step
  //   Solve Flow quantity at the new time step
  //   Calc. primitive variables from conserved variables
  //   Reset boundary conditions (ghost node approach)
  //   Output sol. at every "iterout" iterations
  //   calculate residual norms (normalized by initial value)
  //   check for convergence (if converged, exit main loop)
  // End Main Loop
  // TODO: Output solution
  // TODO: Evaluate discreization for error norms for grid convergence (isentropic)

  //double M; //used for debugging M

  //Mesh Specifications
  int cellnum = 10;
  vector<double> xcoords;

  //Object Initializations
  array<double,3>* field; //pointer to solution field solutions
  MeshGen1D Mesh(xmin,xmax,cellnum); //mesh
  Euler1D Euler(xcoords); //for solving Euler eqs.
  SpaceVariables1D Sols(cellnum,field); //for storing solutions
  
  Mesh.GenerateMesh(xcoords); //stores all coords in xcoords list
  //debugging:
  /*array<double,3> init{10,50,100};
  Euler.SetInitialConditions(init,field);
  Tools::print("First,middle,& end pt. of rho: %f,%f,%f\n",field[0][0],field[cellnum/2-1][0],field[cellnum-1][0]);
  Tools::print("First,middle,& end pt. of velocity: %f,%f,%f\n",field[0][1],field[cellnum/2-1][1],field[cellnum-1][1]);
  Tools::print("First,middle,& end pt. of pressure: %f,%f,%f\n",field[0][2],field[cellnum/2-1][2],field[cellnum-1][2]);

  */





  return 0;
}
