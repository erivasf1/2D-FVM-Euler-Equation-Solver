//User-defined functions
#include "TimeIntegrator.h" 

// EULEREXPLICIT DEFINITIONS

//-----------------------------------------------------------
EulerExplicit::EulerExplicit(int &c)
: cellnumber(c)
{}

//-----------------------------------------------------------
vector<double> EulerExplicit::ComputeLocalTimeStep(vector<array<double,3>>* &field,Euler1D* &euler,const double &CFL,double &dx){

  double lambda_max;
  vector<double> time_steps(cellnumber);

  //array<double,cellnumber> time_steps; //stores the computed time steps per interior cell
  for (int i=0;i<euler->total_cellnum;i++){ //looping through all interior no|int nbor)
    if (i==0 | i==1 | i==euler->total_cellnum-2 | i==euler->total_cellnum-1) //skiping the ghost cells
      continue;

    lambda_max = euler->GetLambdaMax(field,i); //obtaining largest eigenvalue per cell
    if (std::isnan(lambda_max)){
      Tools::print("Infinitiy detected!\n");
      Tools::print("Velocity at loc(%d): %f\n",i,(*field)[i][0]);
      Tools::print("Mach # at loc(%d): %f\n",i,euler->GetMachNumber(field,i));
    }
    time_steps[i-2] = CFL * (dx/lambda_max);
    //Tools::print("CFL:%f,dx: %f,lambda_max:%f\n",CFL,dx,lambda_max);
     
  }

  return time_steps;
}

//-----------------------------------------------------------
vector<double> EulerExplicit::ComputeGlobalTimeStep(vector<array<double,3>>* &field,Euler1D* &euler,const double &CFL,double &dx){

  //extracting smallest local time step of all cells
  vector<double> time_steps = ComputeLocalTimeStep(field,euler,CFL,dx); //local time steps list for all cells
  double min_time_step = 1.0e5; //temp. value for min time_step

  for (int n=0;n<cellnumber;n++){
    if (time_steps[n] < min_time_step)
      min_time_step = time_steps[n];
  }

  vector<double> global_time_steps(cellnumber,min_time_step);

  return global_time_steps;

}
//-----------------------------------------------------------
void EulerExplicit::FWDEulerAdvance(vector<array<double,3>>* &field,vector<array<double,3>>* &resid,Euler1D* &euler,vector<double>* &time_steps,vector<double> &xcoords,double &dx,array<double,3> &Omega){

  double vol;
  array<double,3> conserve;
  //use indexing of interior cells for Resid!
  // Field still has ghost cells
  //First, convert to conservative to compute conservative values at new time step
  //Second, extract primitive variables from newly calculated conservative variables
  for (int n=0;n<cellnumber;n++){ //i+2 to skip inflow ghost cells

    vol = MeshGen1D::GetCellVolume(n,dx,xcoords); //acquiring cell vol
    //Tools::print("Volume of cell %d:%f\n",i,vol);
    conserve = euler->ComputeConserved(field,n+2); //!< computing conservative values

    for (int i=0;i<3;i++) // advancing to new timestep of conservative variable
      conserve[i] -= Omega[i]*((*time_steps)[i] / vol) * (*resid)[n][i];

    euler->ComputePrimitive(field,conserve,n+2); //!< extracting new primitive variables

  }
  
  return;

}

//-----------------------------------------------------------
void EulerExplicit::SolutionLimiter(vector<array<double,3>>* &field){

  for (int n=0;n<(int)field->size();n++){
    //Density
    (*field)[n][0] = std::min(Density_max,std::max(Density_min,(*field)[n][0]));

    //Velocity
    (*field)[n][1] = std::min(Velocity_max,std::max(Velocity_min,(*field)[n][1]));

    //Pressure
    (*field)[n][2] = std::min(Pressure_max,std::max(Pressure_min,(*field)[n][2]));

    //Printing out message if limiter kicks in
    if ((*field)[n][0] == Density_max || (*field)[n][0] == Density_min)
      Tools::print("Limiter was hit for density at cell %d | val is now:%e\n",n,(*field)[n][0]);

    if ((*field)[n][1] == Velocity_max || (*field)[n][1] == Velocity_min)
      Tools::print("Limiter was hit for velocity at cell %d | val is now:%e\n",n,(*field)[n][1]);

    if ((*field)[n][2] == Pressure_max || (*field)[n][2] == Pressure_min)
      Tools::print("Limiter was hit for pressure at cell %d | val is now:%e\n",n,(*field)[n][2]);

  }

  return;

}

//-----------------------------------------------------------
void EulerExplicit::UnderRelaxationCheck(array<double,3> ResidPrevNorm,array<double,3> ResidNorm,double C,array<bool,3> &check){

  for (int i=0;i<3;i++){
    if (ResidNorm[i] > C*ResidPrevNorm[i])
      check[i] = true; //assigning true to corresponding equation
  }


  return;
}

//-----------------------------------------------------------
bool EulerExplicit::CheckStallResids(int &count,array<double,3> &ResidNorms,array<double,3> &Prev_ResidualNorms,SpaceVariables1D* &sol){

  int count_tol = 1e5; double diff_tol = 1.0e-2; 
  double resid_avg = sol->ComputeNormAvg(ResidNorms); 
  double prev_resid_avg = sol->ComputeNormAvg(Prev_ResidualNorms); 

  count = (abs(resid_avg-prev_resid_avg) <= diff_tol) ? count + 1 : count;

  bool stall = (count > count_tol) ? true : false;

  return stall;

}

//-----------------------------------------------------------
EulerExplicit::~EulerExplicit(){}
