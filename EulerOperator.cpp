//User-defined functions
#include "EulerOperator.h" 

// EULER1D DEFINITIONS

Euler1D::Euler1D(vector<double> &coords)
  : xcoords(coords) {

  dx = abs(xcoords[1] - xcoords[0]); //computing dx (assuming uniform mesh)
}

//-----------------------------------------------------------
void Euler1D::SetInitialConditions(array<double,3> &init_val,array<double,3>* &field){

  for (int i=0;i<(int)xcoords.size();i++){
    for (int n=0;n<(int)init_val.size();n++) {
      field[i][n] = init_val[n]; //setting the flow quantities (primitive) to the initial conditions

    }    
  }

  return;

}
//-----------------------------------------------------------
void Euler1D::SetBoundaryConditions(vector<array<double,3>> &Field,array<double,3> &init){

  //Inserting ghost cells for inflow and outflow locations
  //iterator it = Field.begin();
  Field.insert(Field.begin(),init); //!< setting ghost cells to initial conditions 
  Field.push_back(init); 
  
}

//-----------------------------------------------------------
array<double,3> Euler1D::ComputeSpatialFlux(array<double,3>* &field,int &loc,int &nbor){

  // flux computed using central quadrature
  double res_rho = (field[loc][0]+field[nbor][0]) / 2.0;
  double res_velocity = (field[loc][1]+field[nbor][1]) / 2.0;
  double res_pressure = (field[loc][2]+field[nbor][2]) / 2.0;

  array<double,3> res = {res_rho,res_velocity,res_pressure};
  return res; 


}

//-----------------------------------------------------------
double Euler1D::ComputeSourceTerm(array<double,3>* &field,int &loc) {

  //Get area for i+1/2 and i-1/2 locations (really just area at i=loc and i=loc-1)
  double A_rface = Tools::AreaVal(xcoords[loc]);
  double A_lface = Tools::AreaVal(xcoords[loc-1]);
  // pressure at i = loc
  double p = field[loc][2];
  // source = P_i * (A_i+1/2 - A_i-1/2) / dx * dx
  double res = p * A_rface - A_lface;

  return res;

}

//-----------------------------------------------------------
double Euler1D::GetEpsilon2(array<double,3>* &field,int &loc) {

  double Nu = GetNu(field,loc);
  double Nuleft = GetNu(field,loc-1);
  double Nuright = GetNu(field,loc+1);
  double Nuright2 = GetNu(field,loc+2);

  double kappa2 = 0.35; //typically from 1/4<kappa2<1/2
  
  double max = std::max({Nu,Nuleft,Nuright,Nuright2}); 
  double res = kappa2 * max;
  return res;

}


//-----------------------------------------------------------
double Euler1D::GetNu(array<double,3>* &field,int loc){

  double res = field[loc][2] - (2.0*field[loc][2]) +  field[loc][2];
  res /= field[loc][2] + (2.0*field[loc][2]) +  field[loc][2];
  res = abs(res);

  return res;

}


//-----------------------------------------------------------
double Euler1D::GetLambda(array<double,3>* &field,int &loc){

  //\bar{lambda_i}
  //double M = sqrt(gamma*R*T); //TODO:speed of sound
  double M;
  double lambda_i = abs(field[loc][1]) + M;

  //\bar{lambda_i+1}
  double lambda_iright = abs(field[loc+1][1]) + M;

  double res = (lambda_i + lambda_iright) / 2.0;
  return res;

}


//-----------------------------------------------------------
array<double,3> Euler1D::Compute2ndOrderDamping(array<double,3>* &field,int loc){

  double lambda = GetLambda(field,loc);
  double epsilon = GetEpsilon2(field,loc);

  double res_rho = lambda*epsilon*(field[loc+1][0] + field[loc+1][0]);
  double res_vel = lambda*epsilon*(field[loc+1][1] + field[loc+1][1]);
  double res_pressure = lambda*epsilon*(field[loc+1][2] + field[loc+1][2]);
  array<double,3> res = {res_rho,res_vel,res_pressure};

  return res;

}



//-----------------------------------------------------------
double Euler1D::GetEpsilon4(array<double,3>* &field,int &loc){

  double epsilon2 = GetEpsilon2(field,loc);
  double kappa4 = 1.0/50.0; //typically ranges from: 1/64<kappa4<1/32
  double res = std::max(0.0,(kappa4 - epsilon2));

  return res;

}


//-----------------------------------------------------------
array<double,3> Euler1D::Compute4thOrderDamping(array<double,3>* &field,int loc){

  double lambda = GetLambda(field,loc);
  double epsilon = GetEpsilon4(field,loc);

  double res_rho = lambda*epsilon*(field[loc+2][0] - 3.0*field[loc+1][0] + 3.0*field[loc][0] - field[loc-1][0]);
  double res_vel = lambda*epsilon*(field[loc+2][1] - 3.0*field[loc+1][1] + 3.0*field[loc][1] - field[loc-1][1]);
  double res_pressure = lambda*epsilon*(field[loc+2][2] - 3.0*field[loc+1][2] + 3.0*field[loc][2] - field[loc-1][2]);

  array<double,3> res = {res_rho,res_vel,res_pressure};

  return res;
}


//-----------------------------------------------------------
void Euler1D::ComputeResidual(){ //TODO



}


//-----------------------------------------------------------


Euler1D::~Euler1D(){}
