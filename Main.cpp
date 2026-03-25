// 2D FVM Euler Eq. Solver - Erick Rivas
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
#include "MeshAccess.hpp"
#include "MeshGen.h"
#include "EulerOperator.h"
#include "DataManager.h"
#include "Output.h"
#include "TimeIntegrator.h" 

using namespace std;

//For ease of defining parameters
enum CASE_2D { 
  INLET,AIRFOIL1,AIRFOIL2};

enum CASE_MMS {
  SUPERSONIC,SUBSONIC};

enum BOUNDARY_COND {
  INFLOW,OUTFLOW,SLIP_WALL,PERIODIC};


int main() {

  //Timing Purposes
  double start_time,stop_time;
  start_time = MPI_Wtime();

  //! INITIALIZATION
  // Scenario
  int scenario = 4; //1 = 1D, 2 = 2D, 3 = 2D MMS, 4 = TRUE CARTESIAN of MMS
  CASE_2D case_2d = INLET;
  CASE_MMS case_mms = SUPERSONIC;

  // Constants for 1D case or True Cartesian 2D MMS case
  [[maybe_unused]]double xmin = 0.0; [[maybe_unused]]double xmax = 1.0;
  [[maybe_unused]]double ymin = 0.0; [[maybe_unused]]double ymax = 1.0;
  [[maybe_unused]]int Nx = 17; [[maybe_unused]]int Ny = 17;

  // Boundary Conditions Specification
  //FOR NOW: MMS case is initialized with MS & boundaries are set to outflow for extrapolating to ghost cells
  BOUNDARY_COND top_cond,btm_cond,left_cond,right_cond;
  if (scenario == 3 || scenario == 4){ //MMS case
    top_cond = INFLOW; //OUTFLOW; 
    btm_cond = INFLOW; //OUTFLOW; 
    left_cond = INFLOW; //OUTFLOW; 
    right_cond = INFLOW; //OUTFLOW;  
  }
  else if (scenario == 2 && case_2d == INLET) { //Inlet case
    top_cond = INFLOW; 
    btm_cond = SLIP_WALL;
    left_cond = SLIP_WALL;
    right_cond = OUTFLOW;
  
  }
  else if (scenario == 2 && (case_2d == AIRFOIL1 || case_2d == AIRFOIL2) ){ //Airfoil case
    top_cond = INFLOW;
    btm_cond = SLIP_WALL;
    left_cond = OUTFLOW;
    right_cond = OUTFLOW;

  }

  else{
    cerr<<"Scenario not supported yet!"<<endl;
    return 0;
  }

  [[maybe_unused]]bool cond_loc{false}; //true for subsonic & false for supersonic (FOR EXACT SOL.)
  [[maybe_unused]]bool cond_bc{true}; //true for subsonic & false for supersonic (FOR OUTFLOW BC)

  // Mesh Specifications
  //AIRFOIL GRIDS
  //[[maybe_unused]]const char* meshfile = "Grids/AirfoilGrids/NACA64A006.extra-coarse.27x14.grd"; //name of 2D file -- Note: set to NULL if 1D case is to be ran
  //[[maybe_unused]]const char* meshfile = "Grids/AirfoilGrids/NACA64A006.medium.193x53.grd"; //name of 2D file -- Note: set to NULL if 1D case is to be ran
  //[[maybe_unused]]const char* meshfile = "Grids/AirfoilGrids/NACA64A006.fine.385x105.grd"; //name of 2D file -- Note: set to NULL if 1D case is to be ran

  // INLET GRIDS
  //[[maybe_unused]]const char* meshfile = "Grids/InletGrids/Inlet.53x17.grd"; //name of 2D file -- Note: set to NULL if 1D case is to be ran
  //[[maybe_unused]]const char* meshfile = "Grids/InletGrids/Inlet.105x33.grd"; //name of 2D file -- Note: set to NULL if 1D case is to be ran
  //[[maybe_unused]]const char* meshfile = "Grids/InletGrids/Inlet.209x65.grd"; //name of 2D file -- Note: set to NULL if 1D case is to be ran
  //[[maybe_unused]]const char* meshfile = "Grids/InletGrids/Inlet.417x129.grd"; //name of 2D file -- Note: set to NULL if 1D case is to be ran

  //MMS GRIDS
  //[[maybe_unused]]const char* meshfile = "Grids/CurvilinearGrids/curv2d9.grd"; //name of 2D file -- Note: set to NULL if 1D case is to be ran
  //[[maybe_unused]]const char* meshfile = "Grids/CurvilinearGrids/curv2d17.grd"; //name of 2D file -- Note: set to NULL if 1D case is to be ran
  //[[maybe_unused]]const char* meshfile = "Grids/CurvilinearGrids/curv2d33.grd"; //name of 2D file -- Note: set to NULL if 1D case is to be ran
  //[[maybe_unused]]const char* meshfile = "Grids/CurvilinearGrids/curv2d65.grd"; //name of 2D file -- Note: set to NULL if 1D case is to be ran
  //[[maybe_unused]]const char* meshfile = "Grids/CurvilinearGrids/curv2d129.grd"; //name of 2D file -- Note: set to NULL if 1D case is to be ran
  //[[maybe_unused]]const char* meshfile = "Grids/CurvilinearGrids/curv2d257.grd"; //name of 2D file -- Note: set to NULL if 1D case is to be ran

  [[maybe_unused]]const char* meshfile = NULL; //! UNCOMMENT IF RUNNING 2D MMS EXACT CARTESIAN CASE
  [[maybe_unused]]int cellnum = 100; //recommending an even number for cell face at the throat of nozzle (NOTE: will get reassigned val. if mesh is provided)

  // Temporal Specifications
  const int iter_max = 10;
  int iterout = 1; //number of iterations per solution output
  const double CFL = 0.98; //CFL number (must <= 1 for Euler Explicit integration)
  //const double CFL = 1e-2; //CFL number (must <= 1 for Euler Explicit integration)
  bool timestep{false}; //true = local time stepping; false = global time stepping
  int time_scheme = 1; //0 for Euler Explicit, 1 for RungeKutta2, 2 for RungeKutta4

  // Flux Specifications
  int flux_scheme{1}; //0=JST, 1=Van Leer, 2 = Roe 
  double epsilon = 1.0; //0 for 1st order and 1 for 2nd order
  bool epsilon_ramp = false; //true to enable ramp from 2nd order to 1st
  [[maybe_unused]] int ramp_start = 1.2e4; 
  [[maybe_unused]] int ramp_stop = 1.2e4;

  int flux_limiter = 0; //0 for disable & 1 for Van Leer
  bool freeze_limiter = false; //true/false for enabling/disabling limiter freeze
  [[maybe_unused]]bool resid_stall{false};//for detecting if residuals have stalled, NOTE: leave as false!
  [[maybe_unused]]int stall_count = 0;

  // Under-Relaxation Parameters
  double C = 1.2; //residual norm check
  array<double,4> Omega{1.0,1.0,1.0,1.0}; //FWD Advance Limiter
  int subiter_max = 0; //max number of relaxation sub-iterations
  array<bool,4> check{false,false,false,false}; //false by default to check if under-relaxation is needed

  // Outputting results parameters
  string results_prefix = "Results/";
  string resids_prefix = "Residuals/";
  string mmserror_prefix = "MMSError/";

  // Governing Eq. Residuals
  double cont_tol = 1.0e-15;
  double xmom_tol = 1.0e-15;
  double ymom_tol = 1.0e-15;
  double energy_tol = 1.0e-15;

  //! GENERATING MESH 
  MeshGenBASE* mesh; 
  if (scenario == 1){ //1D Nozzle Mesh Case -- TODO:After 2D is done
    mesh = new MeshGenNozzle(xmin,xmax,mesh->cellnumber);
    //mesh->GenerateMesh();
  }
  else if ((scenario == 2) || (scenario == 3)) {
    mesh = new MeshGen2D(meshfile);
    mesh->OutputMesh(); //!< outputs mesh file for visualization
  }
  else if ((scenario == 4) && (!meshfile)) {
    mesh = new MeshGen2D(xmin,xmax,ymin,ymax,Nx,Ny);
  }
  else{
    cerr<<"Unknown scenario number!"<<endl;
    return 0;
  }
  
  //! DATA ALLOCATION
  //Field variables
  vector<array<double,4>> Field(mesh->cellnumber); //stores primitive variable sols.
  vector<array<double,4>> FieldStar(mesh->cellnumber); //stores intermediate primitive variable sols.
  vector<array<double,4>> FieldStall(mesh->cellnumber); //stores primitive variable sols. before stall (if detected)
  vector<array<double,4>> FieldMS(mesh->cellnumber); //stores manufactured sol.
  [[maybe_unused]] vector<array<double,4>> FieldMS_Source(mesh->cellnumber); //stores manufactured source term for all cells
  vector<array<double,4>> FieldMS_Error(mesh->cellnumber); //stores manufactured sol. error

  vector<array<double,4>> ExactField(mesh->cellnumber); //stores exact cell-averaged primitve variable sols.
  vector<array<double,4>> ExactFaces(mesh->cellnumber+1); //stores exact primitve variable sols. at cell faces

  vector<array<double,4>> Residual(mesh->cellnumber); //stores the local residuals per eq.
  vector<array<double,4>> ResidualStar(mesh->cellnumber); //stores the intermediate stage of primtive variables
  vector<array<double,4>> InitResidual(mesh->cellnumber); //stores the initial residual

  vector<double> TimeSteps(mesh->cellnumber); //for storing the time step (delta_t) for each cell
  array<double,4> ResidualNorms; //for storing the global residual norms
  array<double,4> ResidualStarNorms; //stores the intermediate global residual norms
  [[maybe_unused]] array<double,4> Prev_ResidualNorms; //for storing the previous global residual norms

  //Pointers to Field variables
  vector<array<double,4>>* field = &Field; //pointer to Field solutions
  vector<array<double,4>>* field_star = &FieldStar; //pointer to intermediate Field solutions
  vector<array<double,4>>* field_stall = &FieldStall; //pointer to stalled field sols
  vector<array<double,4>>* field_ms = &FieldMS; //pointer to manufactured sol. field
  [[maybe_unused]] vector<array<double,4>>* field_ms_source = &FieldMS_Source; //pointer to intermediate Field solutions
  vector<array<double,4>>* field_ms_error = &FieldMS_Error; //pointer to intermediate Field solutions
  [[maybe_unused]] vector<array<double,4>>* exact_sols = &ExactField; //pointer to exact solution field values
  [[maybe_unused]] vector<array<double,4>>* exact_faces = &ExactFaces; //pointer to exact solution field values
  vector<array<double,4>>* resid = &Residual; //pointer to residual field values per cell
  vector<array<double,4>>* resid_star = &ResidualStar; //pointer to intermediate residual field values per cell
  vector<array<double,4>>* init_resid = &InitResidual; //pointer to residual field values per cell
  vector<double>* time_steps = &TimeSteps;

  vector<string> iter_visuals_primitive, iter_visuals_resid; 
  [[maybe_unused]] vector<string> iter_visuals_MMSerror;

  //!OBJECT INITIALIZATIONS

  //Euler Operator spec.
  // Specifying EulerOperator
  EulerBASE* euler;
  //Temp -- will add scenario == 1 once 1D section is fixed!
  if (scenario == 2) 
    euler = new Euler2D(case_2d,mesh->cell_imax,mesh->cell_jmax,flux_scheme,flux_limiter,epsilon,top_cond,btm_cond,left_cond,right_cond,mesh);
  else if ((scenario == 3) || (scenario == 4))
    euler = new Euler2DMMS(mesh->cell_imax,mesh->cell_jmax,flux_scheme,flux_limiter,epsilon,top_cond,btm_cond,left_cond,right_cond,mesh,field_ms_source,field_ms,case_mms);
  else{
    cerr<<"Error: scenario # not recognized!"<<endl;
    return 0;
  }
  
  //Time Integrator spec.
  EulerExplicitBASE* time; //for computing time steps
  if (time_scheme == 0)
    time = new EulerExplicit(mesh,euler,CFL);
  else if (time_scheme == 1)
    time = new RungeKutta2(mesh,euler,CFL);
  else if (time_scheme == 2)
    time = new RungeKutta4(mesh,euler,CFL);
  else 
    cerr<<"Error: unknown time scheme spec.!"<<endl;

  

  SpaceVariables2D Sols; //for operating on Field variables

  Output Error(results_prefix,resids_prefix,mmserror_prefix,iterout,mesh); 

  //Pointers to Objects
  SpaceVariables2D* sols = &Sols;

  Output* output = &Error;

  //! PRINTING OUT SIMULATION INFO
  // TITLE
  if (meshfile)
    Tools::print("2D EULER EQ. SOLVER\n");
  else
    Tools::print("1D EULER EQ. SOLVER\n");
  // CASE SPEC
  if (meshfile){
    Tools::print("-Mesh Selected: ");
    Tools::print("%s\n",meshfile);
  }
  else if (scenario == 4){
    Tools::print("-Case Selected: ");
    Tools::print("True Cartesian MMS Case\n");
  }
  else{
    Tools::print("-Case Selected: ");
    (cond_bc == true) ? Tools::print("Shock Wave Case\n") : Tools::print("Isentropic Case\n");
  }
  // SPATIAL STATS
  Tools::print("-Spatial Statistics:\n");
  Tools::print("--Cell Number: %d\n",mesh->cellnumber);
  //Tools::print("--Delta x: %f\n",dx);

  // TEMPORAL STATS
  Tools::print("-Temporal Statistics:\n");
  Tools::print("--CFL: %f\n",CFL);
  Tools::print("--Time Scheme: ");
  if (time_scheme == 0) //Euler Explicit Time Scheme
    Tools::print("Euler Explicit\n");
  else if (time_scheme == 1) //Runge-Kutta2 Time Scheme
    Tools::print("Runge-Kutta2\n");
  else  //Runge-Kutta4 Time Scheme
    Tools::print("Runge-Kutta4\n");

  // FLUX STATS
  Tools::print("-Flux Statistics:");
  if (flux_scheme == 1) //Van Leer Flux 
    (epsilon == 0.0) ? Tools::print(" 1st Order Van Leer Scheme\n") : Tools::print(" 2nd Order Van Leer Scheme\n");
  else if (flux_scheme == 2) //Roe Flux
    (epsilon == 0.0) ? Tools::print(" 1st Order Roe's Scheme\n") : Tools::print(" 2nd Order Roe's Scheme\n");
  
  else
    Tools::print(" JST Damping\n");

  Tools::print("-Flux Limiter: ");
  if ((flux_limiter == 1) && (epsilon > 0.0) ) 
    Tools::print("Van Leer\n");
  else
    Tools::print("N/A\n");
  
  // AIRFOIL STATS
  if (scenario == 2 && (case_2d == AIRFOIL1 || case_2d == AIRFOIL2)){
    Tools::print("-Airfoil Stats:\n");
    if (case_2d == AIRFOIL1)
      Tools::print("--AOA = 0 degrees\n");
    else if (case_2d == AIRFOIL2)
      Tools::print("--AOA = 8 degrees\n");
    else {}
  }

  // MMS STATS
  if (scenario == 3){
    Tools::print("-MMS Stats:\n");
    if (case_mms == SUPERSONIC)
      Tools::print("--Supersonic Conds. specified\n");
    else if (case_mms == SUBSONIC)
      Tools::print("--Subsonic Conds. specified\n");
    else {}

  }

  //! COMPUTING MANUFACTURED SOLUTION AND SOURCE TERMS (MMS ONLY)
  if ((scenario == 3) || (scenario == 4)){
    euler->SetCoefficients(); //setting supersonic or subsonic coefficients here
    string mms_sol_filename = "ManufacturedSols.vts"; string mms_source_filename = "SourceTerms.vts";
    euler->ManufacturedPrimitiveSols(field_ms,sols); //!< computing manufactured sol.
    euler->EvalSourceTerms(sols); //!< computing manufactured source terms
    output->OutputPrimitiveVariables_VTS(mms_sol_filename,field_ms);
    output->OutputManufacturedSourceTerms(mms_source_filename,field_ms_source);

  }

  //! SETTING INITIAL CONDITIONS
  euler->SetInitialConditions(field);

  string val = output->zeroPad(0,4);
  string init_name = "Results/Iteration_";
  init_name += val;
  init_name += ".vts";
  const char* filename_init = init_name.c_str();
  output->OutputPrimitiveVariables_VTS(filename_init,field);
  iter_visuals_primitive.push_back(val); //inserting initial sol. to iter_visuals for visualization

  //initial error w/ Manufacturd Sol.
  if ( (scenario == 3) || (scenario == 4) ){ 
    string mms_val = output->zeroPad(0,4);
    string text = "MMSError/Iteration_";
    text += mms_val;
    text += ".vts";
    const char* filename_mms_error = text.c_str();
    euler->ComputeMSError(field_ms_error,field,field_ms);
    output->OutputPrimitiveVariables_VTS(filename_mms_error,field_ms_error);
    iter_visuals_MMSerror.push_back(mms_val);
  }

  //!TODO: COMPUTING EXACT SOLUTION -- ONLY FOR 1D QUASI-STEADY NOZZLE
  /*
  if ((cond_bc == false) && (!meshfile)){ //Compute Exact Solution if isentropic case is selected
    array<double,3> sol;
    double area;
    for (int i=0;i<(int)exact_faces->size();i++) {
      area = tool.AreaVal(xcoords[i]);
      cond_loc = (xcoords[i] < 0) ? true:false; 
      SuperSonicNozzle Nozzle(area,area_star,stag_pressure,stag_temp,cond_loc);
      Nozzle.ComputeExactSol(sol);
 
      (*exact_faces)[i] = sol; //assigning to exact faces vector
    
    }
    // Computing cell-average sol. for all cells
    sols->ComputeCellAveragedSol(exact_faces,exact_sols,xcoords);
  }
  */

  // SETTING EPSILON RAMP (IF SPECIFIED)
  if (epsilon_ramp == true && epsilon == 1.0)
    time->RampEpsilon(ramp_start,ramp_stop,0);

  // SETTING BOUNDARY CONDITIONS
  euler->Setup2DBoundaryConditions(field,output); 

  //time->SolutionLimiter(field_test);


  // COMPUTING INITIAL RESIDUAL NORMS
  // using ResidSols spacevariable
  array<double,4> InitNorms;


  euler->ComputeResidual(init_resid,field,field_stall,resid_stall);

  //debug: Residual
  const char* resid_file = "Residuals/InitialResiduals.vts"; 
  output->OutputPrimitiveVariables(init_resid,resid_file,false,0,mesh->xcoords,mesh->ycoords,mesh->cellnumber,mesh->Nx,mesh->Ny);

  InitNorms = sols->ComputeL2SolutionNorms(init_resid); //computing L2 norm of residuals

  Tools::print("-Initial Residual Norms\n");
  Tools::print("--Continuity:%e\n",InitNorms[0]);
  Tools::print("--X-Momentum:%e\n",InitNorms[1]);
  Tools::print("--Y-Momentum:%e\n",InitNorms[2]);
  Tools::print("--Energy:%e\n\n",InitNorms[3]);

  output->OutputResidualNorms(resid_file,0,InitNorms);


  string it,name; //used for outputting file name
  int iter; //iteration number

  //Opening file that stores residuals -- TODO: make this into a fcn.
  ofstream myresids;
  myresids.open("SolResids.dat");
  myresids<<"variables= \"Iteration num.\" \"Continuity\" \"X-Momentum\" \"Y-Momentum\" \"Energy\""<<endl;
  myresids<<"zone T= "<<"\""<<0<<"\""<<endl;
  myresids<<"DATAPACKING=POINT"<<endl;
  myresids<<"DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE )"<<endl;
  //myresids<<"Iteration"<<"  "<<"Contintuity"<<"  "<<"X-Momentum"<<"  "<<"Energy"<<endl;

  myresids<<0<<"  "<<InitNorms[0]<<"  "<<InitNorms[1]<<"  "<<InitNorms[2]<<"  "<<InitNorms[3]<<endl; //printing out the initial residuals first

  //Printing primitive vars. for TECPLOT visualization
  std::string filename_totalsols = "AllSolutions.dat";
  output->OutputPrimitiveVariables(field,filename_totalsols,false,0,mesh->xcoords,mesh->ycoords,mesh->cellnumber,mesh->Nx,mesh->Ny);

  //Assigning Intermediate Field to Initial Field (including residuals)
  (*resid) = (*init_resid); //!< setting residual to initial
  ResidualNorms = InitNorms;

  (*field_star) = (*field); //!< setting intermediate to current
  (*resid_star) = (*resid);

  //! BEGIN OF MAIN LOOP
  for (iter=1;iter<iter_max;iter++){

    // SET EPSILON RAMP (IF SPECIFIED)
    if (epsilon_ramp == true && epsilon == 1.0)
      time->RampEpsilon(ramp_start,ramp_stop,iter);

    // CHECK FOR STALL RESIDUALS (IF LIMITER FREEZE IS SPECIFIED)
    if (euler->epsilon == 1.0 && freeze_limiter == true && resid_stall == false){
      resid_stall = time->CheckStallResids(field,field_stall,ResidualNorms,Prev_ResidualNorms,sols);
    }
  
    //! COMPUTE TIME STEP
    // if global time step, chosen then create a vector<double> of the smallest time step
    (*time_steps) = (timestep == true) ? time->ComputeLocalTimeStep(field) : time->ComputeGlobalTimeStep(field);

    //! COMPUTE NEW SOL. VALUES 
    time->ComputeNewSolution(field_star,resid_star,time_steps,Omega,field_stall,resid_stall);//TESTING
    time->SolutionLimiter(field_star); //applies solution limiter to all cells (including ghost cells)

    //! ENFORCE BOUNDARY CONDITIONS
    euler->Enforce2DBoundaryConditions(field_star,false);
    //time->SolutionLimiter(field_star); //temporarily reapplying the limiter


    euler->ComputeResidual(resid_star,field_star,field_stall,resid_stall);
    ResidualStarNorms = sols->ComputeL2SolutionNorms(resid_star);
    time->UnderRelaxationCheck(ResidualNorms,ResidualStarNorms,C,check);


    //! UNDER-RELAXATION CHECK
    if (check[0]==true || check[1] == true || check[2] == true || check[3] == true){ //perform under-relaxation if any of these are true
      for (int j=0;j<subiter_max;j++){
        for (int i=0;i<4;i++) //!< reassigns omega to half of current value if under-relaxation detected
          Omega[i] = (check[i] == true) ?  Omega[i] /= 2.0 : Omega[i] = 1.0;

        (*field_star) = (*field); //resetting primitive variables to previous time step values

        time->ComputeNewSolution(field_star,resid_star,time_steps,Omega,field_stall,resid_stall); //advancing intermediate solution w/ under-relaxation factor 
        time->SolutionLimiter(field_star); //temporarily reapplying the limiter

        euler->ComputeResidual(resid_star,field_star,field_stall,resid_stall);
        ResidualStarNorms = sols->ComputeL2SolutionNorms(resid_star);

        time->UnderRelaxationCheck(ResidualNorms,ResidualStarNorms,C,check);

        //checking if new residuals do not need under-relaxation
        if (check[0]==false && check[1] == false && check[2] == false && check[3] == false){ 
          for (int i=0;i<4;i++) //!< resetting omega to 1
            Omega[i] = 1.0;
          break;
        }

        if (j == subiter_max)
          Tools::print("Under-relaxation Failed!\n");
      }
    }

    //! SAVING PREV. RESIDUALS (ONLY IF LIMITER FREEZE SPECIFIED)
    if (freeze_limiter == true && resid_stall == false)
      Prev_ResidualNorms = ResidualNorms;
    
    if (resid_stall == true)
      Tools::print("Limiters are now frozen!\n");


    //! ASSINGING NEW TIME STEP VALUES TO INTERMEDIATE VALUES
    (*field) = (*field_star);
    (*resid) = (*resid_star); ResidualNorms = ResidualStarNorms;


    //! OUTPUTING SOL.
    output->WriteSolutions(iter,field,resid,field_ms,field_ms_error,ResidualNorms,euler,scenario,resid_stall,iter_visuals_primitive,iter_visuals_resid,iter_visuals_MMSerror);


    //! CHECK FOR CONVERGENCE 
    //based on intial residual norms
    if (ResidualNorms[0]/InitNorms[0] <= cont_tol && ResidualNorms[1]/InitNorms[1] <= xmom_tol && 
        ResidualNorms[2]/InitNorms[2] <= ymom_tol && ResidualNorms[3]/InitNorms[3] <= energy_tol)
      break;

  }

  //! FINAL OUTPUT OF SOLUTION
  if (iter==iter_max)
    Tools::print("Failed to converge!\n");

  else {
    Tools::print("\n");
    Tools::print("------------------------------------------------------------\n");
    Tools::print("CONGRATS you converged!\n");
    Tools::print("Continuity: %e\nX-Momentum: %e\nY-Momentum: %e\nEnergy: %e\n",ResidualNorms[0],ResidualNorms[1],ResidualNorms[2],ResidualNorms[3]);
    const char* filename_final = "ConvergedSolution.dat" ;
    output->OutputPrimitiveVariables(field,filename_final,false,0,mesh->xcoords,mesh->ycoords,mesh->cellnumber,mesh->Nx,mesh->Ny);

  }

  //Writing PVD file
  const char* primitive_pvd = "Results/solution.pvd";
  const char* resid_pvd = "Residuals/resids.pvd";
  const char* mms_error_pvd = "MMSError/mms_error.pvd";

  output->WritePVDFile(primitive_pvd,iter_visuals_primitive);
  output->WritePVDFile(resid_pvd,iter_visuals_resid);
  output->WritePVDFile(mms_error_pvd,iter_visuals_MMSerror);

  //! INLET ONLY: COMPUTE TOTAL PRESSURE LOSS (SRQ)
  if (case_2d == 0 && scenario == 2){
    const char* file_pressureloss = "InletPressureLoss.txt"; //writing to file to save pressure loss
    double P_loss = euler->ComputePressureLoss(file_pressureloss,field);
    P_loss /= 1000.0;
    Tools::print("Inlet Pressure loss = %f kPa\n",P_loss);

  }

  //! MMS ONLY: EVALUATING ERROR NORMS FOR OBSERVED ORDER OF ACCURACY (SRQ)
  if (scenario == 3 || scenario == 4){

    vector<array<double,4>> Errors(mesh->cellnumber);
    vector<array<double,4>>* errors = &Errors;
    const char* file_errornorms = "DiscretizationErrorNorms.txt";

    output->DiscretizationErrorNorms(field,field_ms,errors,sols,file_errornorms); //writing to file saving error norms
  }

  //AIRFOIL ONLY: COMPUTING LIFT AND DRAG COEFFICIENT (SRQ)
  if ((case_2d == 1 || case_2d == 2) && scenario == 2){
    array<double,2> Airfoil_coeffs = euler->ComputeLiftAndDragCoefficient(field);
    double C_D = Airfoil_coeffs[0];
    double C_L = Airfoil_coeffs[1];
    Tools::print("Airfoil Drag Coefficient = %f \n",C_D);
    Tools::print("Airfoil Lift Coefficient = %f \n",C_L);

    const char* file_surfacep = "AirfoilPressureDistribution.txt";
    euler->ComputeSurfacePresssureCoefficient(file_surfacep,field);
  }

  //Closing Residuals file
  myresids.close();



  stop_time = MPI_Wtime();
  double elapsed_min = (stop_time-start_time) / 60.0;
  double elapsed_hr = elapsed_min / 60.0;
  Tools::print("Elapsed time: %f s | %f min. | %f hrs.\n",stop_time-start_time,elapsed_min,elapsed_hr);

  //! CLEANUP
  delete euler;
  delete time;
  delete mesh;

  return 0;
}
