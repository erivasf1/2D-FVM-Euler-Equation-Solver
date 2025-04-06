//User-defined functions
#include "EulerOperator.h" 

// EULER1D DEFINITIONS

Euler1D::Euler1D(int& cellnum,double &P0,double &BP,double &T0,double &g)
  : interior_cellnum(cellnum), stag_pressure(P0) , back_pressure(BP), stag_temperature(T0) 
  , gamma(g) {}


//-----------------------------------------------------------
Euler1D::Euler1D() {}


//-----------------------------------------------------------
array<double,3> Euler1D::ComputeConserved(vector<array<double,3>>* &field,int loc){

  array<double,3> res; //to store conserved variables

  //density
  res[0] = (*field)[loc][0]; 

  //x-mom.
  res[1] = (*field)[loc][0] * (*field)[loc][1]; 

  //energy
  res[2] = (*field)[loc][2]/(gamma-1.0) + (0.5)* (*field)[loc][0]*pow((*field)[loc][1],2.0);

  return res;

}

//-----------------------------------------------------------
void Euler1D::ComputePrimitive(vector<array<double,3>>* &field,array<double,3> &Conserved,int loc) {

  //Computing Primitive variables given the conserved variables

  //density
  (*field)[loc][0] = Conserved[0];

  //x-velocity
  (*field)[loc][1] = Conserved[1] / Conserved[0];

  //pressure
  (*field)[loc][2] = (gamma-1.0) * (Conserved[2]-0.5*(pow(Conserved[1],2.0)/Conserved[0]));

  return;
}



//-----------------------------------------------------------
void Euler1D::SetInitialConditions(vector<array<double,3>>* &field,vector<double> &xcoords){

  // Mach number is set to vary linearly via M(x) = 9/10(x) + 1
  // And flow quantities are calculated using isentropic conditions
  // ASSUMPTION: Mach number at left face is equal to the cell averaged Mach number of a given cell (may be fine as an initial condition)
  double M,psi,T,a;

  //Tools::print("SetInitialConditions\n");
  for (int i=0;i<interior_cellnum;i++){
   // Tools::print("cell index: %d\n",i);

    M = (9.0/10.0)*xcoords[i] + 1.0; //local Mach number
    psi = 1.0+(gamma-1.0)/2.0 * pow(M,2.0);
    
    //pressure calc.
    (*field)[i][2] = pow(psi,gamma/(gamma-1.0));
    (*field)[i][2] = stag_pressure / (*field)[i][2]; 
    //Tools::print("pressure: %f\n",field[i][2]);
    
    //density calc.
    T = stag_temperature / psi; // local temperature
    (*field)[i][0] = (*field)[i][2] / (R*T); 
    //Tools::print("density: %f\n",field[i][0]);

    //velocity calc.
    a = sqrt(gamma*R*T); //local speed of sound
    (*field)[i][1] = abs(M*a);
    //Tools::print("velocity: %f\n",field[i][1]);
     
  }    
  
  return;

}
//-----------------------------------------------------------
void Euler1D::SetBoundaryConditions(vector<array<double,3>>* &field,bool &cond){

  //Inserting 4 ghost cells for inflow and outflow locations (2 for each)
  //iterator it = Field.begin();
  array<double,3> empty{0.0,0.0,0.0};
  //Inserting empty solution vectors into inflow and outflow ghost cells
  //TODO: Check if 2 empty lists were added to the start of Field domain instead of overwriting the 1st 2 elements
  //Inflow
  field->insert(field->begin(),empty); //!< temporarily setting ghost cells to empty arrays 
  field->insert(field->begin(),empty); 

  //Outflow
  field->push_back(empty); 
  field->push_back(empty); 

  total_cellnum = (int)field->size(); //!< saving the new total num. of cells (w/ ghost cells)

  //Calculating Boundary Condition values
  ComputeTotalBoundaryConditions(field,cond);

  return;
  
}

//-----------------------------------------------------------
void Euler1D::ComputeTotalBoundaryConditions(vector<array<double,3>>* &field,bool &cond){

  //Combines setting the inflow and outflow boundary conditions
  ComputeInflowBoundaryConditions(field);

  ComputeOutflowBoundaryConditions(field,cond);


  return;
}
//-----------------------------------------------------------
void Euler1D::ComputeInflowBoundaryConditions(vector<array<double,3>>* &field){

  // Domain:
  // [G1,G2,I1,...,Imax,G3,G4] --> READ THIS!!!
  //Inflow -- linear extrapolation of Mach Number (refer to class notes section 3 slide 34)
  // note: use the absolute velocity when computing Mach number to prevent negative velocities at the inflow
  double M0,M1,M2; //temp. values
  double psi;
  double T,a;
  double M_limit = 1.0e-3; //Mach num. limiter
  for (int i=1;i>-1;i--){ //reverse for loop for ease of indexing!
    // i=0 (G1) & i=1 (G2)
    //acquiring Mach number from neighboring nodes
    M1 = GetMachNumber(field,i+1);
    M2 = GetMachNumber(field,i+2);
    M0 = 2.0*M1 - M2;

    //applying limiter
    M0 = (M0<M_limit) ? M_limit : M0;

    //USING ISENTROPIC CONDITIONS for thermodynamic variables

    psi = 1.0+(gamma-1.0)/2.0 * pow(M0,2.0);
    //pressure calc.
    (*field)[i][2] = pow(psi,gamma/(gamma-1.0));
    (*field)[i][2] = stag_pressure / (*field)[i][2]; 
    
    //density calc.
    T = stag_temperature / psi; // local temperature
    (*field)[i][0] = (*field)[i][2] / (R*T); 

    //velocity calc.
    a = sqrt(gamma*R*T); //local speed of sound
    (*field)[i][1] = abs(M0*a);
    

    //Tools::print("B.C.\n");
    //Tools::print("Cell index:%d\n",i);
    //Tools::print("Density: %f, Velocity:%f,Pressure:%f\n",field[i][0],field[i][1],field[i][2]);
   }

  return;
}

//-----------------------------------------------------------
void Euler1D::ComputeOutflowBoundaryConditions(vector<array<double,3>>* &field,bool& cond){

  //TODO: Figure out why [i-1] node is not retrieving the node of interior node

  if (cond == false){ //supersonic case 
    //using simple extrapolation from class notes section 3 slide 36
    int c = (int)field->size() - 2; //index of the 1st outflow ghost cell
    for (int i=0;i<2;i++){ //for both outflow ghost cells
      (*field)[i+c][0] = 2.0* (*field)[(i+c)-1][0] - (*field)[(i+c)-2][0]; //density
      (*field)[i+c][1] = 2.0* (*field)[(i+c)-1][1] - (*field)[(i+c)-2][1]; //velocity
      (*field)[i+c][2] = 2.0* (*field)[(i+c)-1][2] - (*field)[(i+c)-2][2]; //pressure

    }

  }

  else { //subsonic
    //fixing back pressure at boundary, not at ghost cell
    int c = (int)field->size() - 2; //index of the 1st outflow ghost cell
    (*field)[c][2] = 2.0*back_pressure - (*field)[c-1][2]; //extrapolated pressure with back pressure enforced
    (*field)[c+1][2] = 2.0* (*field)[c][2] - (*field)[c-1][2]; //regular extrapolated pressure for last ghost cell

    for (int n=c;n<c+2;n++){ //reg. extrapolation
      for (int i=0;i<2;i++) //only iterating for density & velocity
        (*field)[n][i] = 2.0* (*field)[n-1][i] - (*field)[n-2][i]; //density
    }

  }


  return;
}

//-----------------------------------------------------------
array<double,3> Euler1D::ComputeSpatialFlux_BASE(vector<array<double,3>>* &field,int loc,int nbor){

  array<double,3> V_face; //primitive variable vector

  //RECONSTRUCTION: approximating primitive variables at face using central quadrature
  for (int i=0;i<3;i++)
    V_face[i] = 0.5*((*field)[loc][i]+(*field)[nbor][i]);

  array<double,3> flux;

  //Continuity Flux
  flux[0] = V_face[0]*V_face[1]; //rho*u

  //X-Momentum Flux
  flux[1] = V_face[0]*pow(V_face[1],2.0) + V_face[2]; //rho*u^2 + P

  //Energy Flux
  flux[2] =(gamma/(gamma-1.0))*V_face[2]*V_face[1] + (V_face[0]*pow(V_face[1],3.0))*0.5; //(gamma/gamma-1)*P*u + (rho*u^3)/2


  return flux; 

}

//-----------------------------------------------------------
array<double,3> Euler1D::ComputeSpatialFlux_UPWIND1stOrder(vector<array<double,3>>* &field,bool method,int loc,int nbor){

  array<double,3> flux; //total flux
  array<double,3> flux_rtstate; //right state flux (for upwinding in +c wave speed)
  array<double,3> flux_ltstate; //left state flux (for upwinding in -c wave speed)

  if (method == true){ //Van Leer Method
    flux_rtstate = VanLeerCompute(field,nbor,false); //false for negative c case
    flux_ltstate = VanLeerCompute(field,loc,true); //true for positive c case
   
  }
//  if (method == false) //Roe's Method
 //   flux = RoeCompute();

  for (int n=0;n<3;n++) //summing up left and right state fluxes
    flux[n] = flux_rtstate[n] + flux_ltstate[n];


  return flux;
}
//-----------------------------------------------------------
array<double,3> Euler1D::VanLeerCompute(vector<array<double,3>>* &field,int loc,bool sign){

  array<double,3> flux;
  //scalars
  double M = GetMachNumber(field,loc); //local Mach Number
  double a = (*field)[loc][1] / M; //local speed of sound
 
  //total energy term(h_t)
  double ht = (gamma/(gamma-1.0))*((*field)[loc][2]/(*field)[loc][0]); //pressure work
  ht += pow((*field)[loc][1],2.0) / 2.0; //kinetic energy

  //vectors
  array<double,3> convect_vec{1.0,(*field)[loc][1],ht}; //vector multiple for convective flux
  array<double,3> pressure_vec{0.0,(*field)[loc][2],0.0}; //vector multiple for pressure flux

  //Convective Flux
  double C = GetC(M,sign);
  for (int n=0;n<3;n++)
    flux[n] = (*field)[loc][0]*a*C*convect_vec[n]; 

  //Pressure Flux
  double D = GetD(M,sign);
  for (int n=0;n<3;n++)
    flux[n] += D*pressure_vec[n];
  

  return flux;
}

//-----------------------------------------------------------
double Euler1D::GetC(double M,bool sign){

  double alpha = GetAlpha(M,sign);
  double beta = GetBeta(M);
  double M_vl = GetVanLeerM(M,sign);

  double C = alpha*(1.0+beta)*M - beta*M_vl;
 
  return C;

}
//-----------------------------------------------------------
double Euler1D::GetAlpha(double M,bool sign){

  //sign = true for positive and false for negative
  double alpha = (sign == true) ? 0.5*(1.0+std::copysign(1.0,M)) : 0.5*(1.0-std::copysign(1.0,M));

  return alpha;

}
//-----------------------------------------------------------
double Euler1D::GetBeta(double M){

  double beta = -std::max(0.0,1.0-(int)abs(M)); //(int) -- explicit type conversion

  return beta;
}
//-----------------------------------------------------------
double Euler1D::GetVanLeerM(double M,bool sign){
  
  double M_VL = (sign == true) ? 0.25*pow((M+1.0),2.0) : -0.25*pow((M-1.0),2.0);
  return M_VL;
}
//-----------------------------------------------------------
double Euler1D::GetD(double M,bool sign){

  double alpha = GetAlpha(M,sign);
  double beta = GetBeta(M);
  double p2bar = GetP2Bar(M,sign);

  double D = alpha*(1.0+beta) - beta*p2bar;

  return D;
}
//-----------------------------------------------------------
double Euler1D::GetP2Bar(double M,bool sign){

  double M_vl = GetVanLeerM(M,sign);

  double p2bar = (sign == true) ? M_vl*(-M+2.0) : M_vl*(-M-2.0);

  return p2bar;
}
//-----------------------------------------------------------
double Euler1D::ComputeSourceTerm(vector<array<double,3>>* &field,int loc,vector<double> &xcoords) {

  //Note: loc is index of cell in Field -> xcoords index is i-2 of left face
  //Get area for i+1/2 and i-1/2 locations (refer to read me for data indexing)
  dx = abs(xcoords[0]-xcoords[1]); //assigning dx val. here
  double A_rface = Tools::AreaVal(xcoords[loc-1]);
  double A_lface = Tools::AreaVal(xcoords[loc-2]);
  // pressure at i = loc
  double p = (*field)[loc][2];
  double dAdx = (A_rface-A_lface) / dx;
  // source = P_i * (A_i+1/2 - A_i-1/2) / dx * dx
  double res = p * dAdx;

  return res;

}

//-----------------------------------------------------------
double Euler1D::GetEpsilon2(vector<array<double,3>>* &field,int loc) {

  double Nu = GetNu(field,loc);
  double Nuleft = GetNu(field,loc-1);
  double Nuright = GetNu(field,loc+1);
  double Nuright2 = GetNu(field,loc+2);

  double kappa2 = 1.0/2.0; //typically from 1/4<kappa2<1/2
  
  double max = std::max({Nu,Nuleft,Nuright,Nuright2}); //acquring max nu
  double res = kappa2 * max;
  return res;

}


//-----------------------------------------------------------
double Euler1D::GetNu(vector<array<double,3>>* &field,int loc){

  int total_cells = (int)field->size();
  double res;
  //Tools::print("loc:%d\n",loc);

  if (loc == total_cells-1){ //check for evaluating nu at last outflow ghost cell
    res = (*field)[loc][2] - (2.0* (*field)[loc][2]) +  (*field)[loc-1][2]; //assuming Pressure is same at cell loc
    res /= (*field)[loc][2] + (2.0* (*field)[loc][2]) +  (*field)[loc-1][2]; //assuming Pressure is same at cell loc
    res = 0.0;
    return res;
  }

  if (loc == 0){ //check for evaluating new at 1st inflow ghost cell
    res = (*field)[loc+1][2] - (2.0* (*field)[loc][2]) +  (*field)[loc][2]; //assuming Pressure is same at cell loc
    res /= (*field)[loc+1][2] + (2.0* (*field)[loc][2]) +  (*field)[loc][2]; //assuming Pressure is same at cell loc
    return res;
  }

  res = (*field)[loc+1][2] - (2.0* (*field)[loc][2]) +  (*field)[loc-1][2];
  res /= (*field)[loc+1][2] + (2.0* (*field)[loc][2]) +  (*field)[loc-1][2];
  res = abs(res);

  return res;

}


//-----------------------------------------------------------
double Euler1D::GetMachNumber(vector<array<double,3>>* &field,int loc){

  //Using insentropic conditions only for thermodynamic variables
  double T = (*field)[loc][2] / ((*field)[loc][0]*R); //if T is negative than M will be -nan!!!

  //using actual x-vel. value now
  double M = (*field)[loc][1] / sqrt(gamma * R * T); //M = u/a

  return M;

}


//-----------------------------------------------------------
double Euler1D::GetLambda(vector<array<double,3>>* &field,int loc){

  double lambda_i = GetLambdaMax(field,loc); //Lambda bar

  double lambda_iright = GetLambdaMax(field,loc+1); //Lambda_right bar

  double res = (lambda_i + lambda_iright) / 2.0;
  return res;

}


//-----------------------------------------------------------
array<double,3> Euler1D::Compute2ndOrderDamping(vector<array<double,3>>* &field,int loc){

  //following Roy's class notes nomenclature (section 3 slide 31)
  //returns solely \arrow{d^2} vector!
  //look into applying damping terms to boundary?
  array<double,3> conserved = ComputeConserved(field,loc);
  array<double,3> conserved_nbor = ComputeConserved(field,loc+1);
  
  double lambda = GetLambda(field,loc); //at cell face (i+1/2)
  double epsilon = GetEpsilon2(field,loc); //sensor for detecting shocks (will have to tweak the constant later)

  array<double,3> res;
   for (int i=0;i<3;i++) //computing 2nd order damping here
    res[i] = lambda*epsilon*(conserved_nbor[i]-conserved[i]);
  
  return res;

}



//-----------------------------------------------------------
double Euler1D::GetEpsilon4(vector<array<double,3>>* &field,int loc){

  double epsilon2 = GetEpsilon2(field,loc);
  double kappa4 = 1.0/30.0; //typically ranges from: 1/64<kappa4<1/32
  double res = std::max(0.0,(kappa4 - epsilon2));

  return res;

}


//-----------------------------------------------------------
array<double,3> Euler1D::Compute4thOrderDamping(vector<array<double,3>>* &field,int loc){

  double lambda = GetLambda(field,loc);
  double epsilon = GetEpsilon4(field,loc); //looks for gradients that cause odd-even decoupling (will have to tweak the constant later)
  //TODO: Convert to conservative variables HERE ONLY!
  array<double,3> conserved = ComputeConserved(field,loc);
  array<double,3> conserved_leftnbor = ComputeConserved(field,loc-1);
  array<double,3> conserved_rightnbor = ComputeConserved(field,loc+1);
  array<double,3> conserved_right2nbor = ComputeConserved(field,loc+2);

  array<double,3> res;
  for (int i=0;i<3;i++) //computing 4th order damping term here
    res[i] = lambda*epsilon*(conserved_right2nbor[i] - 3.0*conserved_rightnbor[i] + 3.0*conserved[i] - conserved_leftnbor[i]);


  return res;
}


//-----------------------------------------------------------
void Euler1D::ComputeResidual(vector<array<double,3>>* &resid,vector<array<double,3>>* &field,vector<double> &xcoords,double &dx,bool flux_scheme,bool flux_accuracy,bool upwind_scheme){ 

  //following nomenclature from class notes
  [[maybe_unused]] array<double,3> F_right,F_left; //left and right face spatial fluxes 
  [[maybe_unused]] array<double,3> D2_right,D2_left,D4_right,D4_left; //left and right face damping terms
  array<double,3> TotalF_right,TotalF_left; //left and right total fluxes (spatial + damping)
  array<double,3> Source{0.0,0.0,0.0}; //source term vector
  double S; //source term (only for x-mom. eq.)
  double A_left,A_right; //Area of corresponding faces

  for (int n=0;n<total_cellnum;n++){ //looping through all interior nodes
    if ((n==0) | (n==1) | (n==total_cellnum-2) | (n==total_cellnum-1)) //skipping the ghost cell nodes
      continue;


    //Fluxes Evaluation
    //note: \arrow{F}_(i-1/2) is same as \arrow{F}_(i+1/2) of cell to the left!
    if (flux_scheme == true){ //JST Damping Case -- central quadrature
      // Spatial Fluxes
       F_right = ComputeSpatialFlux_BASE(field,n,n+1);
       F_left = ComputeSpatialFlux_BASE(field,n-1,n);

      //JST Damping Terms (need a D2_left flux and D2_right flux vector; similar for D4)
      //note: \arrow{D}_(i-1/2) is same as \arrow{D}_(i+1/2) of cell to the left!
      // right face
      D2_right = Compute2ndOrderDamping(field,n);
      D4_right = Compute4thOrderDamping(field,n);
      // left face
      D2_left = Compute2ndOrderDamping(field,n-1);
      D4_left = Compute4thOrderDamping(field,n-1);

      //Total Flux Terms
      for (int i=0;i<3;i++){
        TotalF_right[i] = F_right[i] - (D2_right[i]-D4_right[i]);
        TotalF_left[i] = F_left[i] - (D2_left[i]-D4_left[i]);
      }
    }

    //TODO: Upwind Cases
    else{
      if (flux_accuracy == true){ //1st order accurate case
        TotalF_right = ComputeSpatialFlux_UPWIND1stOrder(field,upwind_scheme,n,n+1);
        TotalF_left = ComputeSpatialFlux_UPWIND1stOrder(field,upwind_scheme,n-1,n);

      }
    }

    //Source Term (external pressure) ONLY for x-mom. eq.
    // also, area already evaluated, but may need to be multiplied dx?
    S = ComputeSourceTerm(field,n,xcoords);
    Source[1] = S; //adding scalar to vector


    //Area Evaluations
    A_left = Tools::AreaVal(xcoords[n-2]);
    A_right = Tools::AreaVal(xcoords[n-1]);

    //Residual cal.
    for (int i=0;i<3;i++)
      (*resid)[n-2][i] = (TotalF_right[i]*A_right - TotalF_left[i]*A_left) - Source[i]*dx;
    

  }
  return;

}

//-----------------------------------------------------------
double Euler1D::GetLambdaMax(vector<array<double,3>>* &field,int loc){

  double M = GetMachNumber(field,loc); //cell-averaged Mach number
  double a = abs((*field)[loc][1]) / M; //cell-averaged speed of sound
  double lambda_max = abs((*field)[loc][1]) + a;

  return lambda_max;

}

//-----------------------------------------------------------
double Euler1D::GetCellAverageSol(double &A_left,double &A_right,double &dx,array<double,3> &sol_left,array<double,3> &sol_right){

  //using weighted geometric cell-average formula
  double num = 0.5*(A_right + A_left) * (sol_right[1]+sol_left[1]);
  double vol = 0.5*(A_right + A_left) * dx;

  double sol = num / vol;

  return sol;

}

//-----------------------------------------------------------


Euler1D::~Euler1D(){}
