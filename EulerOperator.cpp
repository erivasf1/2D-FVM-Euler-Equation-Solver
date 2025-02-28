//User-defined functions
#include "EulerOperator.h" 

// EULER1D DEFINITIONS

Euler1D::Euler1D(vector<double> &coords,double &P0,double &T0,double &g)
  : xcoords(coords) , stag_pressure(P0) , stag_temperature(T0) 
  , gamma(g) {

  dx = abs(xcoords[1] - xcoords[0]); //computing dx (assuming uniform mesh)
  instance = this; //assings the instance to the object
}

//-----------------------------------------------------------
void Euler1D::SetInitialConditions(array<double,3>* &field){

  // Mach number is set to vary linearly via M(x) = 9/10(x) + 1
  // And flow quantities are calculated using isentropic conditions
  double M,psi,T,a;

  for (int i=0;i<(int)xcoords.size();i++){
    M = (9.0/10.0)*xcoords[i] + 1.0; //local Mach number
    psi = 1.0+(gamma-1.0)/2.0 * pow(M,2.0);
    
    //pressure calc.
    field[i][2] = pow(psi,gamma/(gamma-1.0));
    field[i][2] = stag_pressure / field[i][2]; 
    
    //density calc.
    T = stag_temperature / psi; // local temperature
    field[i][0] = field[i][2] / (R*T); 

    //velocity calc.
    a = sqrt(gamma*R*T); //local speed of sound
    field[i][1] = abs(M*a);
     
  }    
  
  return;

}
//-----------------------------------------------------------
void Euler1D::SetBoundaryConditions(vector<array<double,3>> &Field,array<double,3>* &field,bool &cond){

  //Inserting 4 ghost cells for inflow and outflow locations (2 for each)
  //iterator it = Field.begin();
  array<double,3> init{0.0,0.0,0.0};
  //Inserting empty solution vectors into inflow and outflow ghost cells
  //Inflow
  Field.insert(Field.begin(),init); //!< setting ghost cells to initial conditions 
  Field.insert(Field.begin(),init); 

  //Outflow
  Field.push_back(init); 
  Field.push_back(init); 

  //Calculating Boundary Condition values
  ComputeTotalBoundaryConditions(field,cond);

  return;
  
}

//-----------------------------------------------------------
void Euler1D::ComputeTotalBoundaryConditions(array<double,3>* &field,bool &cond){

  //Combines setting the inflow and outflow boundary conditions
  ComputeInflowBoundaryConditions(field);

  ComputeOutflowBoundaryConditions(field,cond);

  return;
}
//-----------------------------------------------------------
void Euler1D::ComputeInflowBoundaryConditions(array<double,3>* &field){

  // Domain:
  // [G1,G2,I1,...,Imax,G3,G4] --> READ THIS!!!
  //Inflow -- linear extrapolation of Mach Number (refer to class notes section 3 slide 34)
  // note: use the absolute velocity when computing Mach number to prevent negative velocities at the inflow
  double M0,M1,M2; //temp. values
  double psi;
  double T,a;
  for (int i=1;i>-1;i--){ //reverse for loop for ease of indexing!
    // i=0 (G1) & i=1 (G2)
    //acquiring Mach number from neighboring nodes
    M1 = GetMachNumber(field,i+1);
    M2 = GetMachNumber(field,i+2);
    M0 = 2.0*M1 - M2;
    psi = 1.0+(gamma-1.0)/2.0 * pow(M0,2.0);

    //pressure calc.
    field[i][2] = pow(psi,gamma/(gamma-1.0));
    field[i][2] = stag_pressure / field[i][2]; 
    
    //density calc.
    T = stag_temperature / psi; // local temperature
    field[i][0] = field[i][2] / (R*T); 

    //velocity calc.
    a = sqrt(gamma*R*T); //local speed of sound
    field[i][1] = abs(M0*a);
    }
   

  return;
}

//-----------------------------------------------------------
void Euler1D::ComputeOutflowBoundaryConditions(array<double,3>* &field,bool& cond){

  int n = (int)xcoords.size();
  if (cond == false){ //supersonic case 
    //using simple extrapolation from class notes section 3 slide 36
    for (int i=n;i<n+2;i++){ //for both outflow ghost cells
      field[i][0] = 2.0*field[i-1][0] - field[i-2][0]; //density
      field[i][1] = 2.0*field[i-1][1] - field[i-2][1]; //velocity
      field[i][2] = 2.0*field[i-1][2] - field[i-2][2]; //pressure
    }

  }

  else { //subsonic
    //fixing back pressure at boundary, not at ghost cell
    double Pb = 120.0; //kPa
    field[n][2] = 2.0*Pb - field[n-1][2]; //pressure at G3
    field[n+1][2] = 2.0*field[n][2] - field[n-1][2]; //extrapolated pressure for G4
    for (int i=n;i<n+2;i++){ //veloctiy and pressure extrapolations for both outflow ghost cells
      field[i][0] = 2.0*field[i-1][0] - field[i-2][0]; //density
      field[i][1] = 2.0*field[i-1][1] - field[i-2][1]; //velocity
    }



  }


  return;
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
double Euler1D::GetMachNumber(array<double,3>* field,int loc){

  // Using isentropic conditions
  double psi = (stag_pressure/field[loc][2]);
  psi = pow(psi,(gamma-1.0)/gamma);
  
  double M = (2.0*(psi-1.0)) / (gamma-1.0);
  M = sqrt(M);

  return M;

}


//-----------------------------------------------------------
double Euler1D::GetLambda(array<double,3>* &field,int &loc){

  //\bar{lambda_i}
  //double a = sqrt(gamma*R*T); //TODO:speed of sound (define a fcn. for this)
  // T
  double M = GetMachNumber(field,loc);
  double a = field[loc][1];
  double lambda_i = abs(field[loc][1]) + a;

  //\bar{lambda_i+1}
  double lambda_iright = abs(field[loc+1][1]) + a;

  double res = (lambda_i + lambda_iright) / 2.0;
  return res;

}


//-----------------------------------------------------------
array<double,3> Euler1D::Compute2ndOrderDamping(array<double,3>* &field,int loc){

  double lambda = instance->GetLambda(field,loc);
  double epsilon = instance->GetEpsilon2(field,loc);

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
