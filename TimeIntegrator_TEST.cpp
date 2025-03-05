//User-defined functions
#include "TimeIntegrator_TEST.h" 

// EULEREXPLICIT DEFINITIONS

//-----------------------------------------------------------
EulerExplicit::EulerExplicit(int &c)
: cellnumber(c)
{}

//-----------------------------------------------------------
vector<double> EulerExplicit::ComputeLocalTimeStep(vector<array<double,3>> &Field,Euler1D &Euler,const double &CFL,double &dx){

  double lambda_max;
  vector<double> time_steps(cellnumber);

  //array<double,cellnumber> time_steps; //stores the computed time steps per interior cell
  for (int i=0;i<Euler.total_cellnum;i++){ //looping through all interior no|int nbor)
    if (i==0 | i==1 | i==Euler.total_cellnum-2 | i==Euler.total_cellnum-1) //skiping the ghost cells
      continue;

    lambda_max = Euler.GetLambdaMax(Field,i); //obtaining largest eigenvalue per cell
    if (std::isnan(lambda_max)){
      Tools::print("Infinitiy detected!\n");
      Tools::print("Velocity at loc(%d): %f\n",i,Field[i][0]);
      Tools::print("Mach # at loc(%d): %f\n",i,Euler.GetMachNumber(Field,i));
    }
    time_steps[i-2] = CFL * (dx/lambda_max);
    //Tools::print("CFL:%f,dx: %f,lambda_max:%f\n",CFL,dx,lambda_max);
     
  }

  return time_steps;
}

//-----------------------------------------------------------
double EulerExplicit::ComputeGlobalTimeStep(const double &CFL,double &dx,double &lambda_max){

  //TODO

  double dt;
  return dt;

}
//-----------------------------------------------------------
void EulerExplicit::FWDEulerAdvance(vector<array<double,3>> &Field,vector<array<double,3>> &Resid,Euler1D &Euler,vector<double> &time_steps,vector<double> &xcoords,double &dx){

  double vol;
  array<double,3> conserve;
  //use indexing of interior cells for Resid!
  // Field still has ghost cells
  //First, convert to conservative to compute conservative values at new time step
  //Second, extract primitive variables from newly calculated conservative variables
  for (int i=0;i<cellnumber;i++){ //i+2 to skip inflow ghost cells
    
    vol = MeshGen1D::GetCellVolume(i,dx,xcoords); //acquiring cell vol
    //Tools::print("Volume of cell %d:%f\n",i,vol);
    conserve = Euler.ComputeConserved(Field,i); //!< computing conservative values

    //new continuity
    //Tools::print("previous density :%f\n",Field[i+2][0]);
    conserve[0] -= (time_steps[i] / vol) * Resid[i][0];
    //Tools::print("time step :%f\n",time_steps[i]);
    //Tools::print("continuity Resid :%f\n",Resid[i][0]);
    //Tools::print("new density :%f\n",Field[i+2][0]);

    //new velocity
    //Tools::print("previous velocity :%f\n",Field[i+2][1]);
    conserve[1] -= (time_steps[i] / vol) * Resid[i][1];
    //Tools::print("X-mom. Resid :%f\n",Resid[i][1]);

    //new pressure
    //Tools::print("previous pressure :%f\n",Field[i+2][2]);
    conserve[2] -= (time_steps[i] / vol) * Resid[i][2];
    //Tools::print("Energy Resid :%f\n",Resid[i][2]);

    Euler.ComputePrimitive(Field,conserve,i); //!< extracting new primitive variables

    //Applying Solution Limiter
    SolutionLimiter(Field[i+2]);


    /*
    Tools::print("new density :%f\n",Field[i+2][0]);
    Tools::print("new velocity :%f\n",Field[i+2][1]);
    Tools::print("new pressure :%f\n",Field[i+2][2]);
    Tools::print("------\n");
    */ 

  }

}

//-----------------------------------------------------------
void EulerExplicit::SolutionLimiter(array<double,3> &Sol){

  for(int i=0;i<cellnumber;i++){
    //Density
    Sol[0] = std::min(Density_max,std::max(Density_min,Sol[0]));

    //Velocity
    Sol[1] = std::min(Velocity_max,std::max(Velocity_min,Sol[1]));

    //Pressure
    Sol[2] = std::min(Pressure_max,std::max(Pressure_min,Sol[2]));

    //Printing out message if limiter kicks in
    if (Sol[0] == Density_max || Sol[0] == Density_min)
      Tools::print("Limiter was hit for density at cell %d\n",i);

    if (Sol[1] == Velocity_max || Sol[1] == Velocity_min)
      Tools::print("Limiter was hit for velocity at cell %d\n",i);

    if (Sol[2] == Pressure_max || Sol[2] == Pressure_min)
      Tools::print("Limiter was hit for pressure at cell %d\n",i);
  }
  return;

}

//-----------------------------------------------------------
EulerExplicit::~EulerExplicit(){}
