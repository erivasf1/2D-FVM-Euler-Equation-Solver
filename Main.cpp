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
  // Create Mesh (verified)
  // TODO: Set Initial Conditions
  // TODO: Calculate Exact Sol. for comparison  
  // TODO: Set Boundary Conditions (ghost node approach)
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
  MeshGen1D Mesh(xmin,xmax,cellnum);
  Mesh.GenerateMesh(xcoords);
  
  //debugging:
  /*tool.print("Xcoords:\n");
  tool.print("Xcoords size:%d\n",(int)xcoords.size());
  for (int n=0;n<(int)xcoords.size();n++){
    tool.print("%e\n",xcoords[n]);
  }
  */

  tool.print("---------\n");
  tool.print("RESULTS\n");
  tool.print("---------\n");
  //SuperSonicNozzle nozzle(xmin,xmax,pt_num,stag_pressure,stag_temp);  

  //nozzle.ComputeExactSol();
  
  //nozzle.RetrievePoints(); 
  //debug compute mach number
  /*nozzle.print("X. loc.: -1.0\n");
  double a = nozzle.AreaVal(-1.0);
  double m = nozzle.ComputeMachNumber(a);
  nozzle.print("Area:%f & Mach #:%f\n",a,m);*/




  return 0;
}
