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

  //Timing Purposes
  double start_time,stop_time;
  start_time = MPI_Wtime();

  //! INITIALIZATION
  // Constants
  double xmin = -1.0;
  double xmax = 1.0;
  double stag_pressure = 300.0 * 1000.0; //kPa -> Pa
  double back_pressure = 120.0 * 1000.0; //kPa -> Pa (for subsonic outflow cond.)
  double stag_temp = 600.0; //K
  double gamma = 1.4; //specific heat ratio
  double area_star = Tools::AreaVal(0.5*(xmin+xmax)); //area at throat

  // Case Specification
  bool cond_loc{false}; //true for subsonic & false for supersonic (FOR EXACT SOL.)
  bool cond_bc{false}; //true for subsonic & false for supersonic (FOR OUTFLOW BC)

  // Mesh Specifications
  int cellnum = 50; //recommending an even number for cell face at the throat of nozzle
  vector<double> xcoords; //stores the coords of the cell FACES!!! (i.e. size of xcoords is cellnum+1)!

  // Temporal Specifications
  const int iter_max = 1e5;
  const int iterout = 100; //number of iterations per solution output
  const double CFL = 0.3; //CFL number (must <= 1 for Euler Explicit integration)
  //const double CFL = 2.9e-4; //CFL number (must <= 1 for Euler Explicit integration)
  bool timestep{false}; //true = local time stepping; false = global time stepping

  // Flux Specifications
  const bool flux_scheme{false}; //true for JST Damping & false for Upwind
  const bool upwind_scheme{true}; //true for Van Leer & false for Rhoe
  const bool flux_accuracy{false}; //true for 1st order & false for 2nd order

  // Under-Relaxation Parameters
  double C = 1.2; //residual norm check
  array<double,3> Omega{1.0,1.0,1.0}; //FWD Advance Limiter
  //int subiter_max = 0; //max number of relaxation sub-iterations
  int subiter_max = 1e2; //max number of relaxation sub-iterations

  // Governing Eq. Residuals
  double cont_tol = 1.0e-10;
  double xmom_tol = 1.0e-10;
  double energy_tol = 1.0e-10;

  //! DATA ALLOCATION
  //Field variables
  vector<array<double,3>> Field(cellnum); //stores primitive variable sols.
  vector<array<double,3>> FieldStar(cellnum); //stores intermediate primitive variable sols.

  vector<array<double,3>> ExactField(cellnum); //stores exact cell-averaged primitve variable sols.
  vector<array<double,3>> ExactFaces(cellnum+1); //stores exact primitve variable sols. at cell faces

  vector<array<double,3>> Residual(cellnum); //stores the local residuals per eq.
  vector<array<double,3>> ResidualStar(cellnum); //stores the intermediate stage of primtive variables
  vector<array<double,3>> InitResidual(cellnum); //stores the initial residual

  vector<double> TimeSteps; //for storing the time step (delta_t) for each cell
  array<double,3> ResidualNorms; //for storing the global residual norms

  //Pointers to Field variables
  vector<array<double,3>>* field = &Field; //pointer to Field solutions
  vector<array<double,3>>* field_star = &FieldStar; //pointer to intermediate Field solutions
  vector<array<double,3>>* exact_sols = &ExactField; //pointer to exact solution field values
  vector<array<double,3>>* exact_faces = &ExactFaces; //pointer to exact solution field values
  vector<array<double,3>>* resid = &Residual; //pointer to residual field values per cell
  vector<array<double,3>>* resid_star = &ResidualStar; //pointer to intermediate residual field values per cell
  vector<array<double,3>>* init_resid = &InitResidual; //pointer to residual field values per cell
  vector<double>* time_steps = &TimeSteps;

  //Object Initializations
  MeshGen1D Mesh(xmin,xmax,cellnum); //mesh

  SpaceVariables1D Sols; //for operating on Field variables

  Tools tool; //utilities object

  Euler1D Euler(cellnum,stag_pressure,back_pressure,stag_temp,gamma); //for performing Euler Eq. operations 

  EulerExplicit Time(cellnum); //for computing time steps

  Output Error; //for discretization error operations

  //Pointers to Objects
  MeshGen1D* mesh = &Mesh;
  SpaceVariables1D* sols = &Sols;
  Euler1D* euler = &Euler;
  EulerExplicit* time = &Time;
  Output* error = &Error;

  //Intermediate variables for under-relaxation
  array<double,3> ResidualStarNorms; //stores the intermediate global residual norms
  array<bool,3> check{false,false,false}; //false by default to check if under-relaxation is needed


  //! GENERATING UNIFORM MESH
  mesh->GenerateMesh(xcoords); //stores all coords in xcoords list
  double dx = abs(xcoords[0]-xcoords[1]); //delta x distance


  //! PRINTING OUT SIMULATION INFO
  // Title
  Tools::print("1D EULER EQ. SOLVER\n");
  // Case Spec
  Tools::print("-Case Selected: ");
  (cond_bc == true) ? Tools::print("Shock Wave Case\n") : Tools::print("Isentropic Case\n");
  // Spatial Stats
  Tools::print("-Spatial Statistics:\n");
  Tools::print("--Cell Number:%d\n",cellnum);
  Tools::print("--Delta x:%f\n",dx);
  // Temporal Stats
  Tools::print("-Temporal Statistics:\n");
  Tools::print("--CFL: %f\n",CFL);
  Tools::print("--Time-stepping method: ");
  (timestep == true) ? Tools::print("Local time-stepping\n") : Tools::print("Global time-stepping\n");  
  // Temporal Stats
  Tools::print("-Flux Statistics:");
  if (flux_scheme == false){ //Upwind Schemes
    if (flux_accuracy == true)//1st order
      (upwind_scheme == true) ? Tools::print(" 1st Order Van Leer Scheme\n") : Tools::print(" 1st Order Roe's Scheme\n");
    else if (flux_accuracy == false)//2nd order
      (upwind_scheme == true) ? Tools::print(" 2nd Order Van Leer Scheme\n") : Tools::print(" 2nd Order Roe's Scheme\n");
  }
  else
    Tools::print(" JST Damping\n");


  //! SETTING INITIAL CONDITIONS
  //Tools::print("At initial conditions\n");
  euler->SetInitialConditions(field,xcoords);

  //! COMPUTING EXACT SOLUTION -- (should be outputted to a file)
  if (cond_bc == false){ //Compute Exact Solution if isentropic case is selected
    array<double,3> sol;
    double area;
    for (int i=0;i<(int)exact_faces->size();i++) {
      area = tool.AreaVal(xcoords[i]);
      cond_loc = (xcoords[i] < 0) ? true:false; 
      SuperSonicNozzle Nozzle(area,area_star,stag_pressure,stag_temp,cond_loc);
      Nozzle.ComputeExactSol(sol);
 
      (*exact_faces)[i] = sol; //assigning to exact faces vector
    
    }
  }

  // Computing cell-average sol. for all cells
  sols->ComputeCellAveragedSol(exact_faces,exact_sols,xcoords,dx);

  //Debug: Temporarily set initial conditions to exact solutions
  //Field = ExactField;
  //Debug: printing initial conditions w/ no BCs
  const char* filename = "InitialSolutions.txt";
  sols->OutputPrimitiveVariables(field,euler,filename);


  // SETTING BOUNDARY CONDITIONS
  euler->SetBoundaryConditions(field,cond_bc);

  time->SolutionLimiter(field);

  //!< Outputting initial solutions with BC's
  const char* filename2 = "InitSolutionswBCs.txt";
  sols->OutputPrimitiveVariables(field,euler,filename2);

  // COMPUTING INITIAL RESIDUAL NORMS
  // using ResidSols spacevariable
  array<double,3> InitNorms;
  euler->ComputeResidual(init_resid,field,xcoords,dx,flux_scheme,flux_accuracy,upwind_scheme); //computing residuals per cell
  InitNorms = sols->ComputeSolutionNorms(init_resid); //computing L2 norm of residuals
  Tools::print("-Initial Residual Norms\n");
  Tools::print("--Continuity:%e\n",InitNorms[0]);
  Tools::print("--X-Momentum:%e\n",InitNorms[1]);
  Tools::print("--Energy:%e\n\n",InitNorms[2]);


  (*resid) = (*init_resid);//!< setting initial residual to intermediate
  ResidualNorms = InitNorms;

  string it,name; //used for outputting file name
  int iter; //iteration number

  //Opening File that stores residuals
  ofstream myresids;
  myresids.open("SolResids.dat");
  myresids<<"variables= \"Iteration num.\" \"Continuity\" \"Momentum\"  \"Energy\""<<endl;
  myresids<<"zone T= "<<"\""<<0<<"\""<<endl;
  myresids<<"DATAPACKING=POINT"<<endl;
  myresids<<"DT=(DOUBLE DOUBLE DOUBLE DOUBLE )"<<endl;
  //myresids<<"Iteration"<<"  "<<"Contintuity"<<"  "<<"X-Momentum"<<"  "<<"Energy"<<endl;

  myresids<<0<<"  "<<InitNorms[0]<<"  "<<InitNorms[1]<<"  "<<InitNorms[2]<<endl; //printing out the initial residuals first

  //Printing to TECPLOT
  std::string filename_totalsols = "AllSolutions.dat";
  sols->AllOutputPrimitiveVariables(field,euler,filename_totalsols,false,0,xcoords);

  //Assigning Intermediate Field to Initial Field (including residuals)
  (*field_star) = (*field);
  (*resid_star) = (*resid);

  //! BEGIN OF MAIN LOOP
  for (iter=1;iter<iter_max;iter++){

    //debugging only
    //Tools::print("Iteration #: %d\n",iter);
  
    //! COMPUTE TIME STEP
    // if global time step, chosen then create a vector<double> of the smallest time step
    //time_steps = Time.ComputeLocalTimeStep(field,Euler,CFL,dx);//TESTING
    (*time_steps) = (timestep == true) ? time->ComputeLocalTimeStep(field,euler,CFL,dx) : time->ComputeGlobalTimeStep(field,euler,CFL,dx);
    //time_steps = Time.ComputeLocalTimeStep(field,Euler,CFL,dx);

    //! COMPUTE NEW SOL. VALUES 
    time->FWDEulerAdvance(field_star,resid_star,euler,time_steps,xcoords,dx,Omega);//TESTING
    time->SolutionLimiter(field_star); //applies solution limiter to all cells (including ghost cells)

    //! COMPUTE BOUNDARY CONDITIONS
    euler->ComputeTotalBoundaryConditions(field_star,cond_bc);
    time->SolutionLimiter(field_star); //temporarily reapplying the limiter


    euler->ComputeResidual(resid_star,field_star,xcoords,dx,flux_scheme,flux_accuracy,upwind_scheme);
    ResidualStarNorms = sols->ComputeSolutionNorms(resid_star);
    time->UnderRelaxationCheck(ResidualNorms,ResidualStarNorms,C,check);


    //! UNDER-RELAXATION CHECK
    if (check[0]==true || check[1] == true || check[2] == true){ //perform under-relaxation if any of these are true
      for (int j=0;j<subiter_max;j++){
        for (int i=0;i<3;i++) //!< reassigns omega to half of current value if under-relaxation detected
          Omega[i] = (check[i] == true) ?  Omega[i] /= 2.0 : Omega[i] = 1.0;

        (*field_star) = (*field); //resetting primitive variables to previous time step values
        time->FWDEulerAdvance(field_star,resid_star,euler,time_steps,xcoords,dx,Omega); //advancing intermediate solution w/ under-relaxation factor 
        euler->ComputeResidual(resid_star,field_star,xcoords,dx,flux_scheme,flux_accuracy,upwind_scheme); //compute under-relaxed residual
        ResidualStarNorms = sols->ComputeSolutionNorms(resid_star);
        time->UnderRelaxationCheck(ResidualNorms,ResidualStarNorms,C,check);

        if (check[0]==false && check[1] == false && check[2] == false){ //checking if new residuals now do not need under-relaxation
        for (int i=0;i<3;i++) //!< resetting omega to 1
          Omega[i] = 1.0;
        break;
        }

        if (j == subiter_max)
          Tools::print("Under-relaxation Failed!\n");
      }
    }

    //Assinging New Time Step Values to Intermediate Values
    (*field) = (*field_star);
    (*resid) = (*resid_star); ResidualNorms = ResidualStarNorms;


    //! OUTPUT SOL. IN TEXT FILE EVERY "ITEROUT" STEPS
    if (iter % iterout == 0) {
      it = to_string(iter);
      name = "SolResults/Iteration";
      name += it;
      name += ".txt";
      const char* filename_iter = name.c_str();
      sols->OutputPrimitiveVariables(field,euler,filename_iter);

      //Printing to TECPLOT
      Sols.AllOutputPrimitiveVariables(field,euler,filename_totalsols,true,iter,xcoords);

      //Printing Residual Norms to Screen
      Tools::print("------Iteration #: %d----------\n",iter);
      Tools::print("Continuity:%e\nX-Momentum:%e\nEnergy:%e\n",ResidualNorms[0],ResidualNorms[1],ResidualNorms[2]);

      // Writing Residuals history to "SolResids.txt" file
      myresids<<iter<<"  "<<ResidualNorms[0]<<"  "<<ResidualNorms[1]<<"  "<<ResidualNorms[2]<<endl;

    }


    //! CHECK FOR CONVERGENCE (w/ respect to the intial residual norms)
    if (ResidualNorms[0]/InitNorms[0] <= cont_tol && ResidualNorms[1]/InitNorms[1] <= xmom_tol && ResidualNorms[2]/InitNorms[2] <= energy_tol)
      break;

  }

  //! FINAL OUTPUT OF SOLUTION

  if (iter==iter_max)
    Tools::print("Failed to converge!\n");

  else {
    Tools::print("\n");
    Tools::print("------------------------------------------------------------\n");
    Tools::print("CONGRATS you converged!\n");
    Tools::print("Continuity: %e\nX-Momentum: %e\nEnergy: %e\n",ResidualNorms[0],ResidualNorms[1],ResidualNorms[2]);
    const char* filename_final = "ConvergedSolution.txt" ;
    sols->OutputPrimitiveVariables(field,euler,filename_final);
  }

  //Closing Residuals file
  myresids.close();


  //! EVALUATE DISCRETIZATION NORMS FOR GRID CONVERGENCE AND PRINT OUT TO FILE
  if (cond_bc == false){
    field->erase(field->begin()); field->erase(field->begin()); //!< erasing ghost cells
    field->erase(field->end()); field->erase(field->end());

    vector<array<double,3>> Errors(Field);
    vector<array<double,3>>* errors = &Errors;

    error->DiscretizationErrorNorms(field,exact_sols,errors,sols);
  }

  stop_time = MPI_Wtime();
  Tools::print("Elapsed time: %fs\n",stop_time-start_time);

  return 0;
}
