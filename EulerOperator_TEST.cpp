//User-defined functions
#include "EulerOperator_TEST.h" 

// EULER1D DEFINITIONS

Euler1D::Euler1D(int& cellnum,double &P0,double &T0,double &g)
  : interior_cellnum(cellnum), stag_pressure(P0) , stag_temperature(T0) 
  , gamma(g) {}


//-----------------------------------------------------------
Euler1D::Euler1D() {}


//-----------------------------------------------------------
array<double,3> Euler1D::ComputeConserved(vector<array<double,3>> &Field,int loc){

  //TODO: Check for any negative values here
  array<double,3> res; //to store conserved variables

  //continuity
  res[0] = Field[loc][0]; 

  //x-mom.
  res[1] = Field[loc][0]*Field[loc][1]; 

  //total energy density (now changed to rho*u, used to be rho*p)
  res[2] = Field[loc][2]/(gamma-1.0) + 0.5*Field[loc][0]*pow(Field[loc][1],2);
  //res[2] = 1.0/(gamma-1.0) + (0.5)*Field[loc][0]*pow(Field[loc][1],2); -old

  return res;

}

//-----------------------------------------------------------
void Euler1D::ComputePrimitive(vector<array<double,3>> &Field,array<double,3> &Conserved,int loc) {

  // TODO: Check for negative values here
  //Computing Primitive variables given the conserved variables
  //density
  Field[loc][0] = Conserved[0];

  //x-velocity
  Field[loc][1] = Conserved[1] / Conserved[0];

  //pressure
  Field[loc][2] = (gamma-1.0) * (Conserved[2]-0.5*(pow(Conserved[1],2)/Conserved[0]));

  return;
}

//-----------------------------------------------------------
void Euler1D::SetInitialConditions(vector<array<double,3>> &Field,vector<double> &xcoords){

  // Mach number is set to vary linearly via M(x) = 9/10(x) + 1
  // And flow quantities are calculated using isentropic conditions
  // ASSUMPTION: Mach number at left face is equal to the cell averaged Mach number of a given cell (may be fine as an initial condition)
  // Also, initializing with isentropic conditions
  double M,psi,T,a;

  //Tools::print("SetInitialConditions\n");
  for (int i=0;i<interior_cellnum;i++){
   // Tools::print("cell index: %d\n",i);

    M = (9.0/10.0)*xcoords[i] + 1.0; //local Mach number
    psi = 1.0+ (gamma-1.0)*0.5 * pow(M,2.0);
    
    //pressure calc.
    Field[i][2] = pow(psi,gamma/(gamma-1.0));
    Field[i][2] = stag_pressure / Field[i][2]; 
    //Tools::print("pressure: %f\n",field[i][2]);
    
    //density calc.
    T = stag_temperature / psi; // local temperature
    Field[i][0] = Field[i][2] / (R*T); 
    //Tools::print("density: %f\n",field[i][0]);

    //velocity calc.
    a = sqrt(gamma*R*T); //local speed of sound
    Field[i][1] = abs(M*a);
    //Tools::print("velocity: %f\n",field[i][1]);
     
  }    
  
  return;

}
//-----------------------------------------------------------
void Euler1D::SetBoundaryConditions(vector<array<double,3>> &Field,bool &cond){

  //Inserting 4 ghost cells for inflow and outflow locations (2 for each)
  //iterator it = Field.begin();
  array<double,3> empty{0.0,0.0,0.0};
  //Inserting empty solution vectors into inflow and outflow ghost cells
  //TODO: Check if 2 empty lists were added to the start of Field domain instead of overwriting the 1st 2 elements
  //Inflow
  Field.insert(Field.begin(),empty); //!< temporarily setting ghost cells to empty arrays 
  Field.insert(Field.begin(),empty); 

  //Outflow
  Field.push_back(empty); 
  Field.push_back(empty); 

  total_cellnum = Field.size(); //!< saving the new total num. of cells (w/ ghost cells)
  //Tools::print("total number of cells after set BC:%d\n",total_cellnum);

  //Calculating Boundary Condition values
  ComputeTotalBoundaryConditions(Field,cond);

  return;
  
}

//-----------------------------------------------------------
void Euler1D::ComputeTotalBoundaryConditions(vector<array<double,3>> &Field,bool &cond){

  //Combines setting the inflow and outflow boundary conditions
  ComputeInflowBoundaryConditions(Field);

  ComputeOutflowBoundaryConditions(Field,cond);


  return;
}
//-----------------------------------------------------------
void Euler1D::ComputeInflowBoundaryConditions(vector<array<double,3>> &Field){

  // Domain:
  // [G1,G2,I1,...,Imax,G3,G4] --> READ THIS!!!
  //Inflow -- linear extrapolation of Mach Number (refer to class notes section 3 slide 34)
  // note: use the absolute velocity when computing Mach number to prevent negative velocities at the inflow

  // note[correction to previous note]: use the actual x-vel. value to compute M, but restrict it to a small positive value in the case it does want to go negative
  double M0,M1,M2; //temp. values
  double psi;
  double T,a;
  double M_limit = 1.0e-3; //Mach num. limiter
  //Tools::print("ComputeInflowBC\n");
  for (int i=1;i>-1;i--){ //reverse for loop for ease of indexing!
    // i=0 (G1) & i=1 (G2)
    //extrapolating Mach Number using slide 34, section 3 nomenclature
    M1 = GetMachNumber(Field,i+1);
    M2 = GetMachNumber(Field,i+2);
    M0 = 2.0*M1 - M2; //extrapolated Mach Num.

    //applying limiter
    M0 = (M0<M_limit) ? M_limit : M0;

    //using isentropic conditions for thermodynamic variables
    psi = 1.0 + (gamma-1.0)*0.5 * pow(M0,2.0);
    //pressure calc.
    Field[i][2] = pow(psi,gamma/(gamma-1.0));
    Field[i][2] = stag_pressure / Field[i][2]; 
    
    //density calc.
    T = stag_temperature / psi; // local temperature
    Field[i][0] = Field[i][2] / (R*T); 

    //velocity calc.
    a = sqrt(gamma*R*T); //local speed of sound
    Field[i][1] = M0 * a;
    

    //Tools::print("B.C.\n");
    //Tools::print("Cell index:%d\n",i);
    //Tools::print("Density: %f, Velocity:%f,Pressure:%f\n",field[i][0],field[i][1],field[i][2]);
   }

  return;
}

//-----------------------------------------------------------
void Euler1D::ComputeOutflowBoundaryConditions(vector<array<double,3>> &Field,bool cond){

  //TODO: TEST THIS!!!
  // From Class Notes Slide 36, Section 3

  if (cond == false){ //supersonic case 
    //using simple extrapolation from class notes section 3 slide 36
    int c = (int)Field.size() - 2; //index of the 1st outflow cell
    //int c = interior_cellnum; //index of the 1st outflow cell
    //Tools::print("Field.size()=%d\n",(int)Field.size());
    //NOTE: index i is use to advance index c to the next outflow cell
    for (int i=0;i<2;i++){ //for both outflow ghost cells
      Field[(i+c)][0] = 2.0*Field[(i+c)-1][0] - Field[(i+c)-2][0]; //density
      Field[(i+c)][1] = 2.0*Field[(i+c)-1][1] - Field[(i+c)-2][1]; //velocity
      Field[(i+c)][2] = 2.0*Field[(i+c)-1][2] - Field[(i+c)-2][2]; //pressure

      //debug:
      /*
      Tools::print("OutflowB.C.\n");
      Tools::print("Cell index:%d\n",i);
      Tools::print("[LeftNeighbor] Density: %f, Velocity:%f,Pressure:%f\n",field[i-1][0],field[i-1][1],field[i-1][2]);
      Tools::print("Density: %f, Velocity:%f,Pressure:%f\n",field[i][0],field[i][1],field[i][2]);
      */
    }

  }

  else { //subsonic
    //fixing back pressure at boundary, not at ghost cell
    int G3 = total_cellnum-2; //index for 1st outflow ghost cell
    double Pb = 120.0; //kPa
    Field[G3][2] = 2.0*Pb - Field[G3-1][2]; //pressure at G3
    Field[G3+1][2] = 2.0*Field[G3][2] - Field[G3-2][2]; //extrapolated pressure for G4

    for (int i=total_cellnum-2;i<total_cellnum;i++){ //reg. extrapolation
      Field[i][0] = 2.0*Field[i-1][0] - Field[i-2][0]; //density
      Field[i][1] = 2.0*Field[i-1][1] - Field[i-2][1]; //velocity
    }


  }


  return;
}

//-----------------------------------------------------------------------------------------
array<double,3> Euler1D::ComputeSpatialFlux(vector<array<double,3>> &Field,int loc,int nbor){

  //Now using Compute Conserved variable functions
  //array<double,3> U; //conserved variable vector
  //array<double,3> U_nbor; //conserved variable vector

  //array<double,3> U_face; //conserved variable at face
  array<double,3> V_face; //primitive variable at face
  //U = ComputeConserved(Field,loc);
  //U_nbor = ComputeConserved(Field,nbor);

  //RECONSTRUCTION: approximating primitive variables at face using central quadrature
  for (int i=0;i<3;i++)
    V_face[i] = 0.5*(Field[loc][i]+Field[nbor][i]);  

  array<double,3> flux;
  //double Fi,F_nbor; //fluxes at cells

  // Value at interface interpolated using central quadrature
  // flux notation from Section3, Slide 28
  /*double flux_continuity = (cv1+cv1_nbor) / 2.0;
  double flux_xmom = (cv2+cv2_nbor) / 2.0;
  double flux_energy = (cv3+cv3_nbor) / 2.0;
  */

  //Continuity Flux
  //Fi = U[1]; F_nbor = U_nbor[1]; //rho*u
  //flux[0] = 0.5*(F_nbor + Fi); 
  flux[0] = V_face[0]*V_face[1]; //rho*u

  //X-Momentum Flux
  //Fi = Field[loc][0]*pow(Field[loc][1],2) + Field[loc][2]; //rho*u^2 + P
  //F_nbor = Field[nbor][0]*pow(Field[nbor][1],2) + Field[nbor][2]; 
  //flux[1] = 0.5*(F_nbor + Fi); 
  flux[1] = V_face[0]*pow(V_face[1],2) + V_face[2]; //rho*u^2 + P

  //Energy Flux
  //Fi = (gamma/(gamma-1.0))*Field[loc][2]*Field[loc][1] + (Field[loc][0]*pow(Field[loc][1],3))*0.5; //(gamma/gamma-1)*P*u + (rho*u^3)/2
  //F_nbor = (gamma/(gamma-1.0))*Field[nbor][2]*Field[nbor][1] + (Field[nbor][0]*pow(Field[nbor][1],3))*0.5;
  //flux[2] = 0.5*(F_nbor + Fi); 
  flux[2] =(gamma/(gamma-1.0))*V_face[2]*V_face[1] + (V_face[0]*pow(V_face[1],3))*0.5; //(gamma/gamma-1)*P*u + (rho*u^3)/2

  return flux; 

}

//-----------------------------------------------------------------------------------------
double Euler1D::ComputeSourceTerm(vector<array<double,3>> &Field,int loc,vector<double> &xcoords) {


  //Note: loc is index of cell in Field -> xcoords index is i-2 of left face
  //Get area for i+1/2 and i-1/2 locations (refer to read me for data indexing)
  
  //double A_rface = Tools::AreaVal(xcoords[loc+1]);
  //double A_lface = Tools::AreaVal(xcoords[loc]);
  dx = abs(xcoords[0]-xcoords[1]); //assigning dx val. here
  double A_rface = Tools::AreaVal(xcoords[loc-1]);
  double A_lface = Tools::AreaVal(xcoords[loc-2]);
  // pressure at i = loc
  double p = Field[loc][2];
  // source = P_i * (A_i+1/2 - A_i-1/2) / dx * dx
  double res = p * ((A_rface - A_lface)/dx);

  return res;

}

//-----------------------------------------------------------
double Euler1D::GetEpsilon2(vector<array<double,3>> &Field,int loc) {

  double Nu = GetNu(Field,loc);
  double Nuleft = GetNu(Field,loc-1);
  double Nuright = GetNu(Field,loc+1);
  double Nuright2 = GetNu(Field,loc+2);

  double kappa2 = 1.0/2.0; //typically from 1/4<kappa2<1/2
  //kappa2 = 0.0;
  
  double max = std::max({Nu,Nuleft,Nuright,Nuright2}); //acquring max nu
  double res = kappa2 * max;
  return res;

}


//-----------------------------------------------------------
double Euler1D::GetNu(vector<array<double,3>> &Field,int loc){

  double res = Field[loc+1][2] - (2.0*Field[loc][2]) +  Field[loc-1][2];
  res /= Field[loc+1][2] + (2.0*Field[loc][2]) +  Field[loc-1][2];
  res = abs(res);

  return res;

}


//-----------------------------------------------------------
double Euler1D::GetMachNumber(vector<array<double,3>> &Field,int loc){

  /*// Using isentropic conditions
  Tools::print("stag_pressure: %f & gamma:%f\n",stag_pressure,gamma);
  double psi = (stag_pressure/Field[loc][2]);
  Tools::print("psi:%f\n",psi);
  Tools::print("pressure:%f\n",Field[loc][2]);
  //psi = psi(gamma-1.0)/gamma);
  psi = std::pow(psi,(gamma-1.0)/gamma);
  Tools::print("psi after pow:%f\n",psi);
  
  double M = (2.0*(psi-1.0)) / (gamma-1.0);
  Tools::print("M before sqrt:%f\n",M);
  M = sqrt(M);
  */

  //Using insentropic conditions only for thermodynamic variables
  double T = Field[loc][2] / (Field[loc][0]*R); //if T is negative than M will be -nan!!!

  //using actual x-vel. value now
  double M = Field[loc][1] / sqrt(gamma * R * T); //M = u/a

  /*if (M<0) { //uncomment for now

    Tools::print("Negative Mach Number Detected!\n");
    Tools::print("pressure at [%d]:%f\n",loc,Field[loc][2]);
    Tools::print("specific gas constant:%f\n",R);
    Tools::print("density:%f\n",Field[loc][0]);
    Tools::print("temp:%f\n",T);
    Tools::print("M before sqrt:%f\n",M);

  }*/

  return M;

}


//-----------------------------------------------------------
double Euler1D::GetLambda(vector<array<double,3>> &Field,int loc){

  //double M,a; //cell-averaged Mach number and speed of sound, respectively
  //\bar{lambda_i} at current cell
  //double a = sqrt(gamma*R*T); //TODO:speed of sound (define a fcn. for this)
  // T
  //Tools::print("In GetLambda fcn.\n");
  //Tools::print("Location: %d\n",loc);
  //M = GetMachNumber(Field,loc); 
  //a = Field[loc][1] * M;
  /*
  if (std::isnan(M) || std::isnan(a)){
    Tools::print("Mach Number:%f\n",M);
    Tools::print("Velocity:%f\n",Field[loc][1]);
    Tools::print("Pressure:%f\n",Field[loc][2]);
    Tools::print("Speed of Sound:%f\n",a);
    Tools::print("Cell Index:%d\n",loc);
  }
  */
  double lambda_i = GetLambdaMax(Field,loc); //Lambda bar
  //double lambda_i = abs(Field[loc][1]) + a;

  //\bar{lambda_i+1} at neighboring cell to the right
  //M = GetMachNumber(Field,loc+1); //cell-averaged Mach Number
  //a = Field[loc][1] / M; //speed of sound
  double lambda_iright = GetLambdaMax(Field,loc+1); //Lambda_right bar
  //double lambda_iright = abs(Field[loc+1][1]) + a;

  double res = (lambda_i + lambda_iright) / 2.0;
  return res;

}


//-----------------------------------------------------------
array<double,3> Euler1D::Compute2ndOrderDamping(vector<array<double,3>> &Field,int loc){

  //following Roy's class notes nomenclature (section 3 slide 31)
  //returns solely \arrow{d^2} vector!
  //look into applying damping terms to boundary?
  //TODO: Compute conservative variables HERE ONLY!
  array<double,3> conserved = ComputeConserved(Field,loc); 
  array<double,3> conserved_nbor = ComputeConserved(Field,loc+1);
  
  //Tools::print("Pressure before lambda:%f\n",Field[loc][2]);
  //Tools::print("Value of loc:%d\n",loc);
  double lambda = GetLambda(Field,loc); //at cell face (i+1/2)
  double epsilon = GetEpsilon2(Field,loc); //sensor for detecting shocks (will have to tweak the constant later)

  array<double,3> res;
  //Tools::print("Lambda: %f & Epsilon: %f\n",lambda,epsilon);
  for (int i=0;i<3;i++) //computing 2nd order damping here
    res[i] = lambda*epsilon*(conserved_nbor[i]-conserved[i]);
  /*
  double res_continuity = lambda*epsilon*(conserved_nbor[0]-conserved[0]);
  double res_xmom = lambda*epsilon*(conserved_nbor[1]-conserved[1]);
  double res_energy = lambda*epsilon*(conserved_nbor[2]-conserved[2]);
  */
  
  /*double res_rho = lambda*epsilon*(field[loc+1][0] - field[loc][0]);
  double res_vel = lambda*epsilon*(field[loc+1][1] + field[loc+1][1]);
  double res_pressure = lambda*epsilon*(field[loc+1][2] + field[loc+1][2]);
  */
  //array<double,3> res = {res_continuity,res_xmom,res_energy};
  //res[0]=0.0;res[1]=0.0;res[2]=0.0; //TODO: temporarily turning off

  return res;

}



//-----------------------------------------------------------
double Euler1D::GetEpsilon4(vector<array<double,3>> &Field,int loc){

  double epsilon2 = GetEpsilon2(Field,loc);
  double kappa4 = 1.0/32.0; //typically ranges from: 1/64<kappa4<1/32
  double res = std::max(0.0,(kappa4 - epsilon2));

  return res;

}


//-----------------------------------------------------------
array<double,3> Euler1D::Compute4thOrderDamping(vector<array<double,3>> &Field,int loc){

  double lambda = GetLambda(Field,loc);
  double epsilon = GetEpsilon4(Field,loc); //looks for gradients that cause odd-even decoupling (will have to tweak the constant later)
  //TODO: Convert to conservative variables HERE ONLY!
  array<double,3> conserved = ComputeConserved(Field,loc); 
  array<double,3> conserved_leftnbor = ComputeConserved(Field,loc-1);
  array<double,3> conserved_rightnbor = ComputeConserved(Field,loc+1);
  array<double,3> conserved_right2nbor = ComputeConserved(Field,loc+2);

  array<double,3> res;
  for (int i=0;i<3;i++) //computing 4th order damping term here
    res[i] = lambda*epsilon*(conserved_right2nbor[i] - 3.0*conserved_rightnbor[i] + 3.0*conserved[i] - conserved_leftnbor[i]);
    /*double res_continuity = lambda*epsilon*(conserved_right2nbor[0] - 3.0*conserved_rightnbor[0] + 3.0*conserved[0] - conserved_leftnbor[0]);
    double res_xmom = lambda*epsilon*(conserved_right2nbor[1] - 3.0*conserved_rightnbor[1] + 3.0*conserved[1] - conserved_leftnbor[1]);
    double res_energy = lambda*epsilon*(conserved_right2nbor[2] - 3.0*conserved_rightnbor[2] + 3.0*conserved[2] - conserved_leftnbor[2]);*/
  

  /*double res_rho = lambda*epsilon*(Field[loc+2][0] - 3.0*Field[loc+1][0] + 3.0*Field[loc][0] - Field[loc-1][0]);
  double res_vel = lambda*epsilon*(field[loc+2][1] - 3.0*field[loc+1][1] + 3.0*field[loc][1] - field[loc-1][1]);
  double res_pressure = lambda*epsilon*(field[loc+2][2] - 3.0*field[loc+1][2] + 3.0*field[loc][2] - field[loc-1][2]);
  */
  //array<double,3> res = {res_continuity,res_xmom,res_energy};

  return res;
}


//-----------------------------------------------------------
void Euler1D::ComputeResidual(vector<array<double,3>> &Resid,vector<array<double,3>> &Field,vector<double> &xcoords,double &dx){ 

  //Tools::print("I am here\n");
  //following nomenclature from class notes
  array<double,3> F_right,F_left; //left and right face spatial fluxes 
  array<double,3> D2_right,D2_left,D4_right,D4_left; //left and right face damping terms
  array<double,3> TotalF_right,TotalF_left; //left and right total fluxes (spatial + damping)
  array<double,3> Source{0.0,0.0,0.0}; //source term vector
  double S; //scalar source term (only for x-mom. eq.)
  double A_left,A_right; //Area of corresponding faces

  //Tools::print("------In Euler.ComputeResidual---------\n");
  total_cellnum = (int)Field.size();
  Tools::print("Total Size: %d\n",total_cellnum);
  for (int n=0;n<total_cellnum;n++){ //looping through all interior nodes
    if (n==0 | n==1 | n==total_cellnum-2 | n==total_cellnum-1) //skipping the ghost cell nodes
      continue;
    
    Tools::print("cell index:%d\n",n);


    //Spatial Flux Term
    //note: \arrow{F}_(i-1/2) is same as \arrow{F}_(i+1/2) of cell to the left!
    //Tools::print("Spatial Flux Energy\n");
    F_right = ComputeSpatialFlux(Field,n,n+1);
    F_left = ComputeSpatialFlux(Field,n-1,n);

    //Source Term (external pressure) ONLY for x-mom. eq.
    // also, area already evaluated, but may need to be multiplied dx?
    //Tools::print("Source Term\n");
    S = ComputeSourceTerm(Field,n,xcoords);//computing source term scalar
    Source[1] = S; //adding scalar to vector
    //Tools::print("S: %f\n",S);

    //JST Damping Terms (need a D2_left flux and D2_right flux vector; similar for D4)
    //note: \arrow{D}_(i-1/2) is same as \arrow{D}_(i+1/2) of cell to the left!
    // right face 
    D2_right = Compute2ndOrderDamping(Field,n);
    D4_right = Compute4thOrderDamping(Field,n);
    // left face
    D2_left = Compute2ndOrderDamping(Field,n-1);
    D4_left = Compute4thOrderDamping(Field,n-1);


    //Total Flux Terms
    for (int i=0;i<3;i++){
      TotalF_right[i] = F_right[i] - (D2_right[i]+D4_right[i]);
      TotalF_left[i] = F_left[i] - (D2_left[i]+D4_left[i]);
    }

    //Area Evaluations
    A_left = Tools::AreaVal(xcoords[n-2]);
    A_right = Tools::AreaVal(xcoords[n-1]);
    //Tools::print("A_left:%f & A_right:%f\n",A_left,A_right);
    //cout<<"A_left: "<<A_left<<"& "<<"A_right: "<<A_right<<endl;

    //Residual cal.
    for (int i=0;i<3;i++)
      Resid[n-2][i] = (TotalF_right[i]*A_right - TotalF_left[i]*A_left) - Source[i]*dx;

    
    //Tools::print("Cell Index: %d\n",i-2);
    
    /*
    //continuity residual (i-2 so that indexing is correct for resid spacevariable pointer) 
    Resid[i-2][0] = (TotalF_right[0]*A_right - TotalF_left[0]*A_left);
    
    //x-mom. Residual (w/ source term)
    Resid[i-2][1] =  (TotalF_right[1]*A_right - TotalF_left[1]*A_left) - S*dx;

    //energy Residual
    Resid[i-2][2] =  (TotalF_right[2]*A_right - TotalF_left[2]*A_left);
    */

  }
  return;

}

//-----------------------------------------------------------
double Euler1D::GetLambdaMax(vector<array<double,3>> &Field,int loc){

  double M = GetMachNumber(Field,loc); //cell-averaged Mach number
  double a = abs(Field[loc][1]) / M; //cell-averaged speed of sound
  //double a = absfield[loc][1] * M; //cell-averaged speed of sound
  double lambda_max = abs(Field[loc][1]) + a;

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
