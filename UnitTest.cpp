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
#include "EulerOperator_TEST.h"
#include "DataManager_TEST.h"
#include "Output.h"
#include "TimeIntegrator_TEST.h"

using namespace std;


//Start using Testcases as cases to test whole classes or segments of the Main file; not individual functions

//MeshGen Fcns.
//double dx=0.0;

/*TEST_CASE( "CellVolumeComputation" ) {
// Initializing Parameters
double xmin = -1.0;
double xmax = 1.0;
double stag_pressure = 300.0; //kpa
double stag_temp = 600.0; //k
double gamma = 1.4; //specific heat ratio
int pt_num = 5; //# of evenly-spaced requested points (including xmin and xmax)
double area;
double area_star; //area at throat
bool cond{false}; //true for subsonic & false for supersonic

//Mesh Specifications
int cellnum = 8; //recommending an even number for cell face at the throat of nozzle
vector<double> xcoords; //!< stores the coords of the cell FACES!!! (i.e. size of xcoords is cellnum+1)!
MeshGen1D Mesh(xmin,xmax,cellnum);
Mesh.GenerateMesh(xcoords);
int loc = 2;
double dx = 0.0;

  REQUIRE( MeshGen1D::GetCellVolume(loc,dx,xcoords) == 0.0);


}
*/


//EulerOperator Fcns.

TEST_CASE(" EulerOperator " ){

  //Mesh Specifications
  double xmin = -1.0; double xmax = 1.0;
  int cellnum = 6;
  vector<double> xcoords;
  MeshGen1D Mesh(xmin,xmax,cellnum);
  Mesh.GenerateMesh(xcoords);

  //Fluid properties
  double P0 = 300.0;
  double T0 = 600.0;
  double gamma = 1.4; //specific heat ratio
  bool cond{false}; //flow condition (false = supersonic & true = subsonic)

  Euler1D Euler(cellnum,P0,T0,gamma);
  array<double,3> empty{0.0,0.0,0.0};
  vector<array<double,3>> Field(cellnum,empty); //initilizing Field vecor
  vector<array<double,3>> expected = Field; //stores the expected values of the whole field
  double expected_density,expected_velocity,expected_pressure;

  SECTION( "Setting Initial Conditions" ) { //Verified
    for (int i=0;i<cellnum;i++) {
      double M = (9.0/10.0)*xcoords[i] + 1.0; //testing for 1st grid cell.
      double psi = 1.0+ (gamma-1.0)*0.5 * pow(M,2.0);


      //expected pressure calc.
      expected_pressure = pow(psi,gamma/(gamma-1.0));
      expected_pressure = P0/expected_pressure;
      expected[i][2] = expected_pressure;
    
      //expected density calc.
      double T = T0 / psi; // local temperature
      expected_density = expected_pressure / (Euler.R*T); 
      expected[i][0] = expected_density;

      //expected velocity calc.
      double a = sqrt(gamma*Euler.R*T); // local temperature
      expected_velocity = abs(M*a);
      expected[i][1] = expected_velocity;


      /*Euler.SetInitialConditions(Field,xcoords);
      REQUIRE(Field[i][0] == Approx(expected_density)); 
      REQUIRE(Field[i][1] == Approx(expected_velocity)); 
      REQUIRE(Field[i][2] == Approx(expected_pressure)); 
      */

  }

    Euler.SetInitialConditions(Field,xcoords);

    REQUIRE(Field.size() == expected.size()); 
    for (int n=0;n<(int)Field.size();n++){
      for (int i=0;i<3;i++) {
        REQUIRE(Field[n][i] == Approx(expected[n][i]));
      }
    }


  }


  SECTION( " Get Mach Number" ){ //Verified
 
    double M_expected; //expected values (initial condition)
    double M; //from fcn.
    Euler.SetInitialConditions(Field,xcoords); //primitive variables

    for(int i=0;i<cellnum;i++){
      M_expected = (9.0/10.0)*xcoords[i] + 1.0;
      M = Euler.GetMachNumber(Field,i);
      //Tools::print("M_expected is %f\n",M_expected);
      //Tools::print("M is %f\n",M);
      CAPTURE(M,M_expected); //used for printing out values if Test fails
      REQUIRE(M_expected == Approx(M));

    }

  }

    Euler.SetInitialConditions(Field,xcoords); //!<initializing the domain

  //Adding the ghost cells to Field
  //Inflow
  Field.insert(Field.begin(),empty); //!< temporarily setting ghost cells to empty arrays 
  Field.insert(Field.begin(),empty); 

  //Outflow
  Field.push_back(empty); 
  Field.push_back(empty); 



  SECTION( "Computing Inflow Boundary Conditions" ){ //Verified
    //Note: Look into solely extrapolating the primitive variables instead of the Mach Number
    
    //Checking if Field now includes ghost cells
    REQUIRE(Field.size() == cellnum+4);

    //Extrapolated Mach Number for Ghost Cells
    double M1,M2;
    M1 = Euler.GetMachNumber(Field,2);
    M2 = Euler.GetMachNumber(Field,3);
    CAPTURE(M1,M2);
    double MG1 = 2.0*M1 - M2;//closest cell from interior; index:1

    M1 = abs(MG1);
    M2 = Euler.GetMachNumber(Field,2);
    CAPTURE(M1,M2);
    double MG2 = 2.0*M1 - M2; //furthest cell from interior; index:0
  
    //Tools::print("MG1:%f & MG2:%f\n",MG1,MG2);

    double expected_density,expected_velocity,expected_pressure;
    double psi;
    double T,a;
    //Extrapolating ghost cell closest to interior cells
    psi = 1.0+ (gamma-1.0)*0.5 * pow(MG1,2.0);
    CAPTURE(psi); //used for printing out values if Test fails

    //Pressure calc.
    expected_pressure = pow(psi,gamma/(gamma-1.0));
    expected_pressure = P0 / expected_pressure; 

    //Density calc.
    T = T0 / psi; // local temperature
    expected_density = expected_pressure / (Euler.R*T);

    //Velocity calc.
    a = sqrt(gamma*Euler.R*T);
    expected_velocity = abs(MG1 * a);

    Euler.ComputeInflowBoundaryConditions(Field);

    REQUIRE(expected_density == Approx(Field[1][0]));
    REQUIRE(expected_velocity == Approx(Field[1][1]));
    REQUIRE(expected_pressure == Approx(Field[1][2]));

    //Extrapolating ghost cell furthest to interior cells
    psi = 1.0+ (gamma-1.0)*0.5 * pow(MG2,2.0);

    //Pressure calc.
    expected_pressure = pow(psi,gamma/(gamma-1.0));
    expected_pressure = P0 / expected_pressure; 

    //Density calc.
    T = T0 / psi; // local temperature
    expected_density = expected_pressure / (Euler.R*T);

    //Velocity calc.
    a = sqrt(gamma*Euler.R*T);
    expected_velocity = abs(MG2 * a);

    CAPTURE(Field[0][1],expected_velocity); //used for printing out values if Test fails
    CAPTURE(MG1,MG2); //used for printing out values if Test fails
    CAPTURE(psi); //used for printing out values if Test fails

    REQUIRE(expected_density == Approx(Field[0][0]));
    REQUIRE(expected_velocity == Approx(Field[0][1]));
    REQUIRE(expected_pressure == Approx(Field[0][2]));

  }


  SECTION( "Computing Outflow Boundary Conditions" ){ //Verified only for supersonic outflow

    //Checking if Field now includes ghost cells
    REQUIRE(Field.size() == cellnum+4);

    //Only extrapolating primitive solution variables
    vector<array<double,3>> V = Field; //primitive sol. field
    double V1,V2,V3; //primitive variables vector 
    int index = V.size()-2; //index of first interior cell
    for (int i=0;i<2;i++) {
      V[index+i][0] = 2.0*V[index+i-1][0] - V[index+i-2][0];
      V[index+i][1] = 2.0*V[index+i-1][1] - V[index+i-2][1];
      V[index+i][2] = 2.0*V[index+i-1][2] - V[index+i-2][2];
    }

    Euler.ComputeOutflowBoundaryConditions(Field,false);
    for (int n=0;n<Field.size();n++){
      for (int i=0;i<3;i++){
      CAPTURE(n,i);
      CAPTURE(cellnum);
      REQUIRE(V[n][i] == Approx(Field[n][i]));
      }

    }

  }

  
  Euler.ComputeTotalBoundaryConditions(Field,cond); //Calls both Inflow and Outflow Boundary Conditions

  SECTION( "Compute Spatial Flux" ){ //Verified

    int cell_index = (int)Field.size()/2; //test cell
    array<double,3> Flux;
    array<double,3> Expected_Flux = Flux;

    array<double,3> Uright; //conserved variable vector to right of test cell
    array<double,3> U; //conserved variable vector at test cell
    //Computing right face flux
    Flux = Euler.ComputeSpatialFlux(Field,cell_index,cell_index+1);
    
    Uright = Euler.ComputeConserved(Field,cell_index+1); 
    U = Euler.ComputeConserved(Field,cell_index); 
    for (int n=0;n<3;n++){ 
      Expected_Flux[n] = (Uright[n]+U[n]) * 0.5;

    }

    for (int n=0;n<3;n++){ 
      CAPTURE(Expected_Flux[n],Flux[n]);
      REQUIRE(Expected_Flux[n] == Approx(Flux[n])); 
    }

  }


}


/*
//EulerOperator Fcns.
TEST_CASE(" GetCellAverageSol" ){
  
  //Fluid Properties
  double stag_pressure = 300.0; //kpa
  double stag_temp = 600.0; //k
  double gamma = 1.4; //specific heat ratio

  //Acquiring xcoords list
  double xmin = -1.0;
  double xmax = 1.0;
  int cellnum = 6;
  vector<double> xcoords; //!< stores the coords of the cell FACES!!! (i.e. size of xcoords is cellnum+1)!
  MeshGen1D Mesh(xmin,xmax,cellnum);
  Mesh.GenerateMesh(xcoords);
  double dx = abs(xcoords[1]-xcoords[0]);

  //Computing exact solution at faces of cell to right of throat 
  //Area
  double area_star = Tools::AreaVal(0.0); //area at throat
  double area_left = Tools::AreaVal(xcoords[cellnum/2]);
  double area_right = Tools::AreaVal(xcoords[(cellnum/2)+1]);
  bool cond{false};
  SuperSonicNozzle NozzleLeft(area_left,area_star,stag_pressure,stag_temp,cond);
  SuperSonicNozzle NozzleRight(area_right,area_star,stag_pressure,stag_temp,cond);
  
  array<double,3> sol_left;
  NozzleLeft.ComputeExactSol(sol_left);

  array<double,3> sol_right;
  NozzleRight.ComputeExactSol(sol_right);

  //Computing Cell Average Sol
  //double sol = Euler1D::GetCellAverageSol(area_left,area_right,dx,sol_left,sol_right);

  //Testing function
  //REQUIRE( MeshGen1D::GetCellVolume(loc,dx,xcoords) == 0.0);
  REQUIRE(Euler1D::GetCellAverageSol(area_left,area_right,dx,sol_left,sol_right) == 3896.9);

}
*/

/*TEST_CASE( " ExactSolutionResiduals" ) {


  //Fluid Property constants
  const double stag_pressure = 300.0; //kpa
  const double stag_temp = 600.0; //k
  const double gamma = 1.4; //specific heat ratio

  //MeshSpecs definition
  double xmin = -1.0;
  double xmax = 1.0;
  int cellnum = 8; //recommending an even number for cell face at the throat of nozzle
  vector<double> xcoords; //!< stores the coords of the cell FACES!!! (i.e. size of xcoords is cellnum+1)!
  double dx = abs(xcoords[1]-xcoords[0]);
  MeshGen1D Mesh(xmin,xmax,cellnum);
  Mesh.GenerateMesh(xcoords);


  //Exact Solution Parameters & Computation
  vector<array<double,3>> ExactField(cellnum);
  array<double,3>* exact_sols; //pointer to exact solution field values
  SpaceVariables1D ExactSols(cellnum,ExactField,exact_sols); //for storing exact solutions

  double area;
  array<double,3> sol;
  area_star = tool.AreaVal(0.0); //area at throat
  //Tools::print("Exact Solution Output\n");
  for (int i=0;i<cellnum;i++) {
    area = tool.AreaVal(xcoords[i]);
    cond = (xcoords[i] < 0) ? true:false; 
    SuperSonicNozzle Nozzle(area,area_star,stag_pressure,stag_temp,cond);
    Nozzle.ComputeExactSol(sol);
 
    exact_sols[i] = sol; //storing solution values to exact sol.
    
    //Tools::print("Point %f\n",xcoords[i]);
    //Tools::print("Density,Velocity,& Pressure: %f,%f,%f\n",field[i][0],field[i][1],field[i][2]);

  }
  //Compute cell-averaged exact solutions





}
*/
//EulerOperator Fcns.

/*TEST_CASE( "ComputeResidualFunction" ) {

  //Fluid Property constants
  const double stag_pressure = 300.0; //kpa
  const double stag_temp = 600.0; //k
  const double gamma = 1.4; //specific heat ratio

  //MeshSpecs definition
  double xmin = -1.0;
  double xmax = 1.0;
  int cellnum = 8; //recommending an even number for cell face at the throat of nozzle
  vector<double> xcoords; //!< stores the coords of the cell FACES!!! (i.e. size of xcoords is cellnum+1)!
  double dx = abs(xcoords[1]-xcoords[0]);
  MeshGen1D Mesh(xmin,xmax,cellnum);
  Mesh.GenerateMesh(xcoords);


  //Exact Solution Parameters & Computation
  vector<array<double,3>> ExactField(cellnum);
  array<double,3>* exact_sols; //pointer to exact solution field values
  SpaceVariables1D ExactSols(cellnum,ExactField,exact_sols); //for storing exact solutions

  double area;
  array<double,3> sol;
  area_star = tool.AreaVal(0.0); //area at throat
  //Tools::print("Exact Solution Output\n");
  for (int i=0;i<cellnum;i++) {
    area = tool.AreaVal(xcoords[i]);
    cond = (xcoords[i] < 0) ? true:false; 
    SuperSonicNozzle Nozzle(area,area_star,stag_pressure,stag_temp,cond);
    Nozzle.ComputeExactSol(sol);
 
    exact_sols[i] = sol; //storing solution values to exact sol.
    
    //Tools::print("Point %f\n",xcoords[i]);
    //Tools::print("Density,Velocity,& Pressure: %f,%f,%f\n",field[i][0],field[i][1],field[i][2]);

  }

  //Residual Computation

  
  Euler1D Euler(cellnum,stag_pressure,stag_temp,gamma); //for solving Euler eqs.

  Euler.ComputeResidual(resid,field,scoords,dx);




}
*/
