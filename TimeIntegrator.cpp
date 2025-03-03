//User-defined functions
#include "TimeIntegrator.h" 

// EULEREXPLICIT DEFINITIONS

//-----------------------------------------------------------
EulerExplicit::EulerExplicit(int &c)
: cellnumber(c)
{}

//-----------------------------------------------------------
vector<double> EulerExplicit::ComputeLocalTimeStep(array<double,3>* &field,Euler1D &Euler,const double &CFL,double &dx){

  double lambda_max;
  vector<double> time_steps(cellnumber);

  //array<double,cellnumber> time_steps; //stores the computed time steps per interior cell
  for (int i=0;i<Euler.total_cellnum;i++){ //looping through all interior no|int nbor)
    if (i==0 | i==1 | i==Euler.total_cellnum-2 | i==Euler.total_cellnum-1) //skiping the ghost cells
      continue;

    lambda_max = Euler.GetLambdaMax(field,i); //obtaining largest eigenvalue per cell
    time_steps[i-2] = CFL * (dx/lambda_max);
    //Tools::print("CFL:%f,dx: %f,lambda_max:%f\n",CFL,dx,lambda_max);
     
  }

  return time_steps;
}

//-----------------------------------------------------------
double EulerExplicit::ComputeGlobalTimeStep(const double &CFL,double &dx,double &lambda_max){


  double dt;
  return dt;

}
//-----------------------------------------------------------
void EulerExplicit::FWDEulerAdvance(array<double,3>* &field,array<double,3>* &resid,vector<double> &time_steps,vector<double> &xcoords,double &dx){

  double vol;
  //use indexing of interior cells!
  for (int i=0;i<cellnumber;i++){ //i+2 to skip inflow ghost cells
    
    vol = MeshGen1D::GetCellVolume(i,dx,xcoords); //acquiring cell vol
    Tools::print("Volume of cell %d:%f\n",i,vol);
    //new density
    Tools::print("previous density :%f\n",field[i+2][0]);
    field[i+2][0] -= (time_steps[i] / vol) * resid[i][0];
    Tools::print("time step :%f\n",time_steps[i]);
    Tools::print("Original resid :%f\n",resid[i][0]);
    Tools::print("new density :%f\n",field[i+2][0]);

    //new velocity
    Tools::print("previous velocity :%f\n",field[i+2][1]);
    field[i+2][1] -= (time_steps[i] / vol) * resid[i][1];
    Tools::print("Original resid :%f\n",resid[i][1]);
    Tools::print("new velocity :%f\n",field[i+2][1]);

    //new pressure
    Tools::print("previous pressure :%f\n",field[i+2][1]);
    field[i+2][2] -= (time_steps[i] / vol) * resid[i][2];
    Tools::print("Original resid :%f\n",resid[i][2]);
    Tools::print("new pressure :%f\n",field[i+2][2]);

  }

}

//-----------------------------------------------------------
EulerExplicit::~EulerExplicit(){}
