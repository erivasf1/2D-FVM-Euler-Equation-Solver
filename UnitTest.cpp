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

//MeshGen Fcns.
//double dx=0.0;

TEST_CASE( "CellVolumeComputation" ) {
// Initializing Parameters
double xmin = -1.0;
double xmax = 1.0;
/*double stag_pressure = 300.0; //kpa
double stag_temp = 600.0; //k
double gamma = 1.4; //specific heat ratio
int pt_num = 5; //# of evenly-spaced requested points (including xmin and xmax)
double area;
double area_star; //area at throat
bool cond{false}; //true for subsonic & false for supersonic
*/
//Mesh Specifications
int cellnum = 8; //recommending an even number for cell face at the throat of nozzle
vector<double> xcoords; //!< stores the coords of the cell FACES!!! (i.e. size of xcoords is cellnum+1)!
MeshGen1D Mesh(xmin,xmax,cellnum);
Mesh.GenerateMesh(xcoords);
int loc = 2;
double dx = 0.0;

  REQUIRE( MeshGen1D::GetCellVolume(loc,dx,xcoords) == 0.0);


}


//DataManager Functions
/*
TEST_CASE(" ConvertToConservative" ){

  vector<array<double,3>> Field(cellnum);
  array<double,3>* field; //pointer to Field solutions
  SpaceVariables1D Sols(cellnum,Field,field); //for storing solutions







}
*/


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
