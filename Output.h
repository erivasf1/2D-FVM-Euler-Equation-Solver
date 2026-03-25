// Responsible for outputting results in a clean format
#ifndef _OUTPUT_H_
#define _OUTPUT_H_
#include <fstream> 
#include <cmath> 

#include "ExactNozzle.h"
#include "MeshGen.h"
#include "DataManager.h"
#include "EulerOperator.h"
//#include "TimeIntegrator.h"

using namespace std;

//Forward Declarations
class SpaceVariables1D;
class SpaceVariables2D;
class EulerBASE;

class Output {

  string results_prefix, resids_prefix, mms_error_prefix;
  int iterout;
  MeshGenBASE* mesh;

  public:
  Output();
  Output(string &fresults, string &fresids, string &fmmserror,int iout, MeshGenBASE* m);
  
  void PrintResidualNorm(int &cellnum,int &n);

  void OutputResidualNorms(const char* &filename,int iter,array<double,4> ResidualNorms);

  void DiscretizationErrorNorms(vector<array<double,4>>* &field,vector<array<double,4>>* &ms_field,vector<array<double,4>>* &errors,SpaceVariables2D* &sols,const char* &filename);

  void CalculateOrderofAccuracy(const char *filename_read,const char *filename_write); //creates a new file containing the order of accuracy value given the discretization error file.txt
 
  void OutputPrimitiveVariables(vector<array<double,4>>* &field,string filename,bool cond,int iter,vector<double> &xcoords,vector<double> &ycoords,int cell_number,int imax,int jmax); //outputs primitive variables in tecplot format

  void OutputPrimitiveVariables_VTS(string filename,vector<array<double,4>>* &field); //outputs primitive variables in .vts format for ParaView visualization

  void OutputManufacturedSourceTerms(string filename,vector<array<double,4>>* &field); //outputs primitive variables in .vts format for ParaView visualization

  void OutputGhostCoords(string filename,vector<double> &xcoords,vector<double> &ycoords,int Nx,int Ny); //for visualizing the ghost cells
  void OutputGhostCells(vector<array<double,4>>* &ghost_cell,string filename,vector<double> &xcoords,vector<double> &ycoords,vector<double> &ghost_xcoords,vector<double> &ghost_ycoords,int Nx,int Ny,int ghost_Nx,int ghost_Ny,int side); //for visualizing the ghost cells

  void WriteSolutions(int iter,vector<array<double,4>>* &field,vector<array<double,4>>* &resid,[[maybe_unused]] vector<array<double,4>>* &field_ms,[[maybe_unused]] vector<array<double,4>>* &mms_error,array<double,4> &ResidualNorms,EulerBASE* &euler,int scenario,bool &resid_stall,vector<string> &iter_visuals_primitive,vector<string> &iter_visuals_resid,[[maybe_unused]] vector<string> &iter_visuals_MMSerror); //outputs: 1) prim. vars. 2)residuals 3) MMS error

  void WritePVDFile(const char* &filename,vector<string> &iter_visuals);
  //void ConvertToDatFile(const char*filename_read,const char *filename_write); //TODO: creates a .dat file of a given .txt file

  string zeroPad(int number, int padWidth);

  ~Output();


};




#endif
