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
#include "EulerOperator_TEST.h"
#include "DataManager_TEST.h"
#include "Output_TEST.h"
#include "TimeIntegrator_TEST.h" 

using namespace std;


int main() {

  // Domain and Fluid Parameters
  double xmin = -1.0;
  double xmax = 1.0;
  double stag_pressure = 300.0 * 1000.0; //kPa -> Pa
  double back_pressure = 120.0 * 1000.0; //kPa -> Pa (for subsonic outflow cond.)
  double stag_temp = 600.0; //K
  double gamma = 1.4; //specific heat ratio
  double area;
  double area_star = Tools::AreaVal(0.5*(xmin+xmax)); //area at throat (midpoint of xmin and xmax)
  bool cond{false}; //true for subsonic & false for supersonic (FOR EXACT SOL.)
  bool cond_bc{false}; //true for subsonic & false for supersonic (FOR OUTFLOW BC)

  //Mesh Specifications
  int cellnum = 100; //recommending an even number for cell face at the throat of nozzle
  vector<double> xcoords; //!< stores the coords of the cell FACES!!! (i.e. size of xcoords is cellnum+1)!

  //Tools::print("DNE xcoords val:%f\n",xcoords[10]);


  //Temporal Specifications
  const int iter_max = 1e6; //max number of iterations
  const int iterout = 50; //number of iterations per solution output
  const double CFL = 0.1; //CFL number (must <= 1 for Euler Explicit integration)
  bool timestep_cond{false}; //true = local time stepping; false = global time stepping

  //Under-Relaxation Parameters
  double C = 1.2; //residual norm check
  array<double,3> Omega{1.0,1.0,1.0}; //FWD Advance Limiter
  int subiter_max = 1e2; //max number of relaxation sub-iterations

  //Governing Eq. Residuals
  double cont_tol = 1.0e-10;
  double xmom_tol = 1.0e-10;
  double energy_tol = 1.0e-10;

  // ALGORITHM:
  // Create Mesh (verified) -- may have to add ghost cells
  // Set Initial Conditions (verified)
  // Calculate Exact Sol. for comparison (verified)
  // Set Boundary Conditions (verified)
  // Output Initial Residual Norm (verified)
  // Output Initial Solution (verified)
  // Begin Main Loop to iteratively solve Euler Eqs.
  //   Set Time Step
  //   Solve Flow quantity at the new time step
  //   Calc. primitive variables from conserved variables //   Reset boundary conditions (ghost node approach)
  //   Output sol. at every "iterout" iterations
  //   calculate residual norms (normalized by initial value)
  //   check for convergence (if converged, exit main loop)
  // End Main Loop
  // Output solution
  // TODO: Evaluate discreization for error norms for grid convergence (isentropic)

  //double M; //used for debugging M


  //Object Initializations
  vector<array<double,3>> Field(cellnum); //stores the primitive variable sols. at all cells
  vector<array<double,3>> ExactField(cellnum); //stores the exact sol. at all cells
  vector<array<double,3>> Residual(cellnum); //stores the residuals for all interior cells
  vector<array<double,3>> InitResidual(cellnum); //stores the initial residuals for all interior cells

  //intermediate variables for under-relaxation 
  vector<array<double,3>> ResidualStar(cellnum); 
  array<double,3> ResidualStarNorms;
  array<bool,3> check{false,false,false}; //false by default


  MeshGen1D Mesh(xmin,xmax,cellnum); //mesh

  //TODO: Only need 1 of these
  // for outputting primitive variables and sol. norms
  SpaceVariables1D Sols; 
  SpaceVariables1D ExactSols; 
  SpaceVariables1D ResidSols; 
  SpaceVariables1D InitResidSols; 

  //for calculating the discretization error
  Output Error;

  array<double,3> ResidualNorms; //stores the norms of the residuals

  Tools tool; //used as utilities object

  EulerExplicit Time(cellnum); //for computing time steps

  // Generating Mesh
  Mesh.GenerateMesh(xcoords); //stores all coords in xcoords list
  double dx = abs(xcoords[0]-xcoords[1]); //delta x distance
  Tools::print("dx:%f\n",dx);
  Euler1D Euler(cellnum,stag_pressure,back_pressure,stag_temp,gamma); //for solving Euler eqs.

  //DEBUG: Displaying Areas
  const char* filename_area = "Areas.txt";
  Mesh.OutputNozzleAreas(xcoords,filename_area);

  // Printing out Temporal Stats
  Tools::print("-Temporal Statistics:\n");
  Tools::print("--CFL: %f\n",CFL);
  Tools::print("Time-stepping method: ");
  (timestep_cond == true) ? Tools::print("Local time-stepping\n") : Tools::print("Global time-stepping\n");  


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
  Euler.SetInitialConditions(Field,xcoords);

  //Debug: printing initial conditions
  //const char* filename = "InitialSolutions.txt"; 
  //Sols.OutputPrimitiveVariables(Field,Euler,filename);

  // COMPUTING EXACT SOLUTION -- (should be outputted to a file)
  vector<array<double,3>> ExactSol_Faces(cellnum+1);
  //Tools::print("Exact Solution Output\n");
  //Computing and storing exact sol. at face
  for (int n=0;n<(int)ExactSol_Faces.size();n++) {
    area = tool.AreaVal(xcoords[n]);
    cond = (xcoords[n] < 0.0) ? true:false; 
    SuperSonicNozzle Nozzle(area,area_star,stag_pressure,stag_temp,cond);
    Nozzle.ComputeExactSol(ExactSol_Faces[n]); 

    //if (n==17){ //debugging supersonic mach number at subsonic section
      //Tools::print("here\n");
      //const char* Mach_filename = "MachNumberVariation.txt";
      //Nozzle.OutputAllMachNumbers(Mach_filename,500);
    //}
 
    //Tools::print("Point %f\n",xcoords[i]);
    //Tools::print("Density,Velocity,& Pressure: %f,%f,%f\n",field[i][0],field[i][1],field[i][2]);
  }

  //Computing cell-average sol. for all cells
  ExactSols.ComputeCellAveragedSol(ExactSol_Faces,ExactField,xcoords,dx);

  //Debug: Temporarily set initial conditions to exact solutions
  //Field = ExactField;
  //Debug: printing initial conditions w/ no BCs
  const char* filename = "InitialSolutions.txt"; 
  Sols.OutputPrimitiveVariables(Field,Euler,filename);
  

  // SETTING BOUNDARY CONDITIONS
  // Note: Was found that extrapolating values at boundary gave a negative pressure
  Euler.SetBoundaryConditions(Field,cond_bc);

  Time.SolutionLimiter(Field);

  
  //!< outputting initial solutions with BC's
  const char* filename2 = "InitSolutionswBCs.txt"; 
  Sols.OutputPrimitiveVariables(Field,Euler,filename2);


  // COMPUTING INITIAL RESIDUAL NORMS
  // using ResidSols spacevariable
  array<double,3> InitNorms;
  Euler.ComputeResidual(InitResidual,Field,xcoords,dx); //computing residuals per cell
  InitNorms = InitResidSols.ComputeSolutionNorms(InitResidual); //computing L2 norm of residuals
  Tools::print("-Initial Residual Norms\n");
  Tools::print("--Continuity:%e\n",InitNorms[0]);
  Tools::print("--X-Momentum:%e\n",InitNorms[1]);
  Tools::print("--Energy:%e\n",InitNorms[2]);
  
  const char* filename3 = "InitialLocalResiduals.txt";
  InitResidSols.OutputLocalResiduals(InitResidual,filename3);


  vector<double> time_steps;
  Residual = InitResidual;
  ResidualNorms = InitNorms;

  string it,name,resid; //used for outputting file name
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
  Sols.AllOutputPrimitiveVariables(Field,Euler,filename_totalsols,false,0,xcoords);

  //Setting up intermediate values
  vector<array<double,3>> FieldStar(Field); 
  ResidualStar = Residual;

  // BEGIN OF MAIN LOOP
  for (iter=1;iter<iter_max;iter++){

    //debugging only
    //Tools::print("------Iteration #: %d----------\n",iter);
  
    //COMPUTE TIME STEP
    // if global time step, chosen then create a vector<double> of the smallest time step
    time_steps = (timestep_cond == true) ? Time.ComputeLocalTimeStep(Field,Euler,CFL,dx) : Time.ComputeGlobalTimeStep(Field,Euler,CFL,dx);

    //COMPUTE NEW SOL. VALUES 
    Time.FWDEulerAdvance(FieldStar,ResidualStar,Euler,time_steps,xcoords,dx,Omega);
    Time.SolutionLimiter(FieldStar); //applies solution limiter to all cells (including ghost cells)
    //Time.FWDEulerAdvance(Field,Residual,Euler,time_steps,xcoords,dx);
    //Time.SolutionLimiter(Field); //applies solution limiter to all cells (including ghost cells)

    //COMPUTE BOUNDARY CONDITIONS
    Euler.ComputeTotalBoundaryConditions(FieldStar,cond_bc);
    Time.SolutionLimiter(FieldStar); //temporarily reapplying the limiter


    //UNDER-RELAXATION CHECK
    // compute residuals and norms
    // check if under-relaxation is needed
    // if "good" then set residual to residul new and same for Field
    // if "bad", then redo Euler FWD Advance with under-relaxation factor & set new values to Field
    // Euler.ComputeResidual(ResidStar,FieldStar,xcoords,dx);
    // ResidStarNorms = ResidSols.ComputeSolutionNorms(ResidStar);
    // Time.UnderRelaxationCheck();
    // if (any in check array is true) -- only checking for bad residual
    //   FieldStar = Field; --> resetting to previous timestep sol.
    //   for max number of relaxation sub-iterations 
    //     for (i=0;i<2;i++) --> setting omega individually
    //       if (check[i] == true) Omega[i] /= 2
    //     Euler.FWDEulerAdvance(FieldStar,ResidStar,Omega);
    //     Euler.ComputeResidual(ResidStar,FieldStar,xcoords,dx);
    //     ResidStarNorms = ResidSols.ComputeSolutionNorms(ResidStar);
    //     Time.UnderRelaxationCheck(ResidPrevNorm,ResidStarNorm,C,check);
    //     if (check==false) {
    //       for (i=0;i<2;i++) --> resetting omega back to 1
    //         Omega[i] = 1; 
    //       break;
    //     }
    //     if (subiter = maxnumsubiter) print("Under-relaxation failed");
    //   end
    //  Residual = ResidStar; Field = FieldStar; ResidualNorms = ResidStarNorms --> assigning the star values to next time step values
    Euler.ComputeResidual(ResidualStar,FieldStar,xcoords,dx);
    ResidualStarNorms = ResidSols.ComputeSolutionNorms(ResidualStar);
    Time.UnderRelaxationCheck(ResidualNorms,ResidualStarNorms,C,check); //C refers to factor of resid increase

    if (check[0]==true || check[1] == true || check[2] == true){ //perform under-relaxation if any of these are true
      for (int j=0;j<subiter_max;j++){
        for (int i=0;i<3;i++) //!< reassigns omega to half of current value if under-relaxation detected
          Omega[i] = (check[i] == true) ?  Omega[i] /= 2.0 : Omega[i] = 1.0;

        FieldStar = Field; //resetting primitive variables to previous time step values
        Time.FWDEulerAdvance(FieldStar,ResidualStar,Euler,time_steps,xcoords,dx,Omega); //advancing intermediate solution w/ under-relaxation factor
        Euler.ComputeResidual(ResidualStar,FieldStar,xcoords,dx);
        ResidualStarNorms = ResidSols.ComputeSolutionNorms(ResidualStar);
        Time.UnderRelaxationCheck(ResidualNorms,ResidualStarNorms,C,check);

        if (check[0]==false && check[1] == false && check[2] == false){ //checking if new residuals now do not need under-relaxation
        for (int i=0;i<3;i++) //!< resetting omega to 1
          Omega[i] = 1.0; 
        break;
        }

        if (j == subiter_max)
          Tools::print("Under-relaxation Failed!\n");
      }
    }
    
    //assigning new time step values to intermediate values
    Field = FieldStar; 
    Residual = ResidualStar; ResidualNorms = ResidualStarNorms;

    //OUTPUT SOL. IN TEXT FILE EVERY "ITEROUT" STEPS
    if (iter % iterout == 0) {
      it = to_string(iter);
      name = "SolResults/Iteration";
      name += it;
      name += ".txt";
      const char* filename_iter = name.c_str(); 
      Sols.OutputPrimitiveVariables(Field,Euler,filename_iter);
 
      //Printing to TECPLOT
      Sols.AllOutputPrimitiveVariables(Field,Euler,filename_totalsols,true,iter,xcoords);
      
      //Printing Residual Norms to Screen
      Tools::print("------Iteration #: %d----------\n",iter);
      Tools::print("Continuity:%e\nX-Momentum:%e\nEnergy:%e\n",ResidualNorms[0],ResidualNorms[1],ResidualNorms[2]);

      // Writing Residuals history to "SolResids.txt" file
      myresids<<iter<<"  "<<ResidualNorms[0]<<"  "<<ResidualNorms[1]<<"  "<<ResidualNorms[2]<<endl;

    }


    //CHECK FOR CONVERGENCE (w/ respect to the intial residual norms)
    if (ResidualNorms[0]/InitNorms[0] <= cont_tol && ResidualNorms[1]/InitNorms[1] <= xmom_tol && ResidualNorms[2]/InitNorms[2] <= energy_tol)
      break;
    
  }

  //Final Output of Solution

  if (iter==iter_max)
    Tools::print("Failed to converge!\n");

  else {
    Tools::print("\n");
    Tools::print("------------------------------------------------------------\n");
    Tools::print("CONGRATS you converged!\n");
    Tools::print("Continuity: %e\nX-Momentum: %e\nEnergy: %e\n",ResidualNorms[0],ResidualNorms[1],ResidualNorms[2]);
    const char* filename_final = "ConvergedSolution.txt" ; 
    Sols.OutputPrimitiveVariables(Field,Euler,filename_final);
  }

  //Closing Residuals file
  myresids.close();


  //Evaluate discretization norms for grid convergence and print out to file
  if (cond_bc == false){
    Field.erase(Field.begin()); Field.erase(Field.begin()); //!< erasing ghost cells
    Field.erase(Field.end()); Field.erase(Field.end());
    
    vector<array<double,3>> Errors(Field);

    Error.DiscretizationErrorNorms(Field,ExactField,Errors,Sols);   
  }


  return 0;
}
