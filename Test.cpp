// Unit testing file that uses the Catch2 header library
#define CATCH_CONFIG_MAIN
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdarg.h>
#include <catch2/catch.hpp> //for unit testing

#include "ExactNozzle.h"
#include "MeshGen.h"
#include "EulerOperator.h"
#include "DataManager.h"
#include "Output.h"
#include "TimeIntegrator.h"

using namespace std;


//Start using Testcases as cases to test whole classes or segments of the Main file; not individual functions


//EulerOperator Fcns.

TEST_CASE(" EulerOperator " ){

  //Mesh Specifications
  int cellnum = 10.0;
  double xmin = -1.0; double xmax = 1.0;
  vector<double> xcoords;

  //Fluid properties
  double P0 = 10.0 * 1000.0;
  double Pb = 10.0 * 1000.0;
  double T0 = 10.0;
  double gamma = 1.4; //specific heat ratio
  bool cond{false}; //flow condition (false = supersonic & true = subsonic)

  // Flux Specifications
  const bool flux_scheme{false}; //true for JST Damping & false for Upwind
  const bool upwind_scheme{false}; //true for Van Leer & false for Rhoe
  const bool flux_accuracy{true}; //true for 1st order & false for 2nd order

  //Field variables
  vector<array<double,3>> Field(cellnum); //stores primitive variable sols.
  vector<array<double,3>> ExpectedField(cellnum); //stores primitive variable sols.

  //Pointers to Field variables
  vector<array<double,3>>* field = &Field; //pointer to Field solutions
  vector<array<double,3>>* expected_field = &ExpectedField; //pointer to Field solutions

  //Object Initializations
  MeshGen1D Mesh(xmin,xmax,cellnum); //mesh

  SpaceVariables1D Sols; //for operating on Field variables

  Tools tool; //utilities object

  Euler1D Euler(cellnum,P0,Pb,T0,gamma); //for performing Euler Eq. operations 

  //Pointers to Objects
  MeshGen1D* mesh = &Mesh;
  Euler1D* euler = &Euler;

  //Generating 1D Mesh
  Mesh.GenerateMesh(xcoords);


  double expected_density,expected_velocity,expected_pressure;

  //Initialize Field with trivial solutions
  for (int i=0;i<(int)Field.size();i++){
    for (int j=0;j<3;j++)
      Field[i][j] = 10.0*j + 1; 
  }  

  //Testing location
  int loc = cellnum/2; int nbor = (cellnum/2) + 1;

  SECTION( "Rho-Avg. Computation" ){
    
    //Testing fcn.
    double abar_test;
    array<double,3> tested_roeavg = euler->ComputeRoeAvgVars(field,loc,nbor,abar_test);
      
    //Expected fcn.
    array<double,3> expected_roeavg;
    double rho_loc = (*field)[loc][0];double rho_nbor = (*field)[nbor][0];
    double u_loc = (*field)[loc][1];double u_nbor = (*field)[nbor][1];
    double p_loc = (*field)[loc][2];double p_nbor = (*field)[nbor][2];
    double ht_loc = (gamma/(gamma-1.0))*(p_loc/rho_loc); //pressure work
    ht_loc += pow(u_loc,2.0) / 2.0; //kinetic energy
    double ht_nbor = (gamma/(gamma-1.0))*(p_nbor/rho_nbor); //pressure work
    ht_nbor += pow(u_nbor,2.0) / 2.0; //kinetic energy

    double R_ihalf = sqrt(rho_nbor/rho_loc);

    expected_roeavg[0] = R_ihalf * rho_loc;
    expected_roeavg[1] = ((R_ihalf*u_nbor) + u_loc) / (R_ihalf + 1.0);
    expected_roeavg[2] = ((R_ihalf*ht_nbor) + ht_loc) / (R_ihalf + 1.0);
    
    for (int i=0;i<3;i++) {
      REQUIRE(tested_roeavg[i] == Approx(expected_roeavg[i]));
    }
    

  }

  double abar;
  array<double,3> roe_avg = euler->ComputeRoeAvgVars(field,loc,nbor,abar);

  SECTION( "Roe Eigen-values" ){

    //Tested fcn.
    array<double,3> test_eigenvals = euler->ComputeRoeEigenVals(roe_avg,abar);

    //Expected fcn.
    array<double,3> expected_eigenvals;
    double ubar = roe_avg[1];
    expected_eigenvals[0] = ubar;
    expected_eigenvals[1] = ubar + abar;
    expected_eigenvals[2] = ubar - abar;

    for (int i=0;i<3;i++) {
      REQUIRE(test_eigenvals[i] == Approx(expected_eigenvals[i]));
    }
    


  }

}

