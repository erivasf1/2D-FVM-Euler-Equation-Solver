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
EulerExplicit::~EulerExplicit(){}
