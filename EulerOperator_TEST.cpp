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

  array<double,3> res; //to store conserved variables

  //continuity
  res[0] = Field[loc][0]; 

  //x-mom.
  res[1] = Field[loc][0]*Field[loc][1]; 

  //energy
  res[2] = 1.0/(gamma-1.0) + (0.5)*Field[loc][0]*pow(Field[loc][2],2);

  return res;

}

//-----------------------------------------------------------
void Euler1D::ComputePrimitive(vector<array<double,3>> &Field,array<double,3> &Conserved,int loc) {

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

  //TODO: TEST THIS!!
  // Domain:
  // [G1,G2,I1,...,Imax,G3,G4] --> READ THIS!!!
  //Inflow -- linear extrapolation of Mach Number (refer to class notes section 3 slide 34)
  // note: use the absolute velocity when computing Mach number to prevent negative velocities at the inflow
  double M0,M1,M2; //temp. values
  double psi;
  double T,a;
  //Tools::print("ComputeInflowBC\n");
  //for (int i=1;i>-1;i--){ //reverse for loop for ease of indexing!
  for (int i=1;i>-1;i--){ //reverse for loop for ease of indexing!
    // i=0 (G1) & i=1 (G2)
    //extrapolating Mach Number using slide 34, section 3 nomenclature
    M1 = GetMachNumber(Field,i+1);
    M2 = GetMachNumber(Field,i+2);
    M0 = 2.0*M1 - M2;
    psi = 1.0+ (gamma-1.0)*0.5 * pow(M0,2.0);
    Tools::print("at %d:M1=%f & M2=%f\n",i,M1,M2);

    //using isentropic conditions for thermodynamic variables
    //pressure calc.
    Field[i][2] = pow(psi,gamma/(gamma-1.0));
    Field[i][2] = stag_pressure / Field[i][2]; 
    
    //density calc.
    T = stag_temperature / psi; // local temperature
    Field[i][0] = Field[i][2] / (R*T); 

    //velocity calc.
    a = sqrt(gamma*R*T); //local speed of sound
    Field[i][1] = abs(M0*a);
    

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

//-----------------------------------------------------------
array<double,3> Euler1D::ComputeSpatialFlux(vector<array<double,3>> &Field,int loc,int nbor){

  //Now using Compute Conserved variable functions
  array<double,3> U; //conserved variable vector
  array<double,3> U_nbor; //conserved variable vector
  U = ComputeConserved(Field,loc);
  U_nbor = ComputeConserved(Field,nbor);
  array<double,3> flux;

  // Value at interface interpolated using central quadrature
  /*double flux_continuity = (cv1+cv1_nbor) / 2.0;
  double flux_xmom = (cv2+cv2_nbor) / 2.0;
  double flux_energy = (cv3+cv3_nbor) / 2.0;
  */
  for (int n=0;n<3;n++)
    flux[n] = (U[n]+U_nbor[n]) / 2.0; //central quadrature

  //array<double,3> flux = {flux_continuity,flux_xmom,flux_energy};
  return flux; 

}

//-----------------------------------------------------------
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
  double res = p * (A_rface - A_lface)/dx;

  return res;

}

//-----------------------------------------------------------
double Euler1D::GetEpsilon2(vector<array<double,3>> &Field,int loc) {

  double Nu = GetNu(Field,loc);
  double Nuleft = GetNu(Field,loc-1);
  double Nuright = GetNu(Field,loc+1);
  double Nuright2 = GetNu(Field,loc+2);

  double kappa2 = 1.0/3.0; //typically from 1/4<kappa2<1/2
  
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
  double T = Field[loc][2] / (Field[loc][0]*R);

  double M = abs(Field[loc][1]) / sqrt(gamma * R * T); //M = u/a

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

  double M,a; //cell-averaged Mach number and speed of sound, respectively
  //\bar{lambda_i} at current cell
  //double a = sqrt(gamma*R*T); //TODO:speed of sound (define a fcn. for this)
  // T
  //Tools::print("In GetLambda fcn.\n");
  //Tools::print("Location: %d\n",loc);
  M = GetMachNumber(Field,loc); 
  a = Field[loc][1] * M;
  if (std::isnan(M) || std::isnan(a)){
    Tools::print("Mach Number:%f\n",M);
    Tools::print("Velocity:%f\n",Field[loc][1]);
    Tools::print("Pressure:%f\n",Field[loc][2]);
    Tools::print("Speed of Sound:%f\n",a);
    Tools::print("Cell Index:%d\n",loc);
  }
  double lambda_i = abs(Field[loc][1]) + a;

  //\bar{lambda_i+1} at neighboring cell to the right
  M = GetMachNumber(Field,loc+1); //cell-averaged Mach Number
  a = Field[loc+1][1] * M;
  double lambda_iright = abs(Field[loc+1][1]) + a;

  double res = (lambda_i + lambda_iright) / 2.0;
  return res;

}


//-----------------------------------------------------------
array<double,3> Euler1D::Compute2ndOrderDamping(vector<array<double,3>> &Field,int loc){

  //following Roy's class notes nomenclature (section 3 slide 31)
  //returns solely \arrow{d^2} vector!
  //look into applying damping terms to boundary?
  //TODO: Compute conservative variables HERE ONLY!
  //Tools::print("In 2nd order damping fcn.\n");
  array<double,3> conserved = ComputeConserved(Field,loc);
  array<double,3> conserved_nbor = ComputeConserved(Field,loc+1);
  
  //Tools::print("Pressure before lambda:%f\n",Field[loc][2]);
  //Tools::print("Value of loc:%d\n",loc);
  double lambda = GetLambda(Field,loc); //at cell face (i+1/2)
  double epsilon = GetEpsilon2(Field,loc); //sensor for detecting shocks (will have to tweak the constant later)

  //Tools::print("Lambda: %f & Epsilon: %f\n",lambda,epsilon);
  double res_continuity = lambda*epsilon*(conserved_nbor[0]-conserved[0]);
  double res_xmom = lambda*epsilon*(conserved_nbor[1]-conserved[1]);
  double res_energy = lambda*epsilon*(conserved_nbor[2]-conserved[2]);
  
  /*double res_rho = lambda*epsilon*(field[loc+1][0] - field[loc][0]);
  double res_vel = lambda*epsilon*(field[loc+1][1] + field[loc+1][1]);
  double res_pressure = lambda*epsilon*(field[loc+1][2] + field[loc+1][2]);
  */
  array<double,3> res = {res_continuity,res_xmom,res_energy};

  return res;

}



//-----------------------------------------------------------
double Euler1D::GetEpsilon4(vector<array<double,3>> &Field,int loc){

  double epsilon2 = GetEpsilon2(Field,loc);
  double kappa4 = 1.0/50.0; //typically ranges from: 1/64<kappa4<1/32
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

  double res_continuity = lambda*epsilon*(conserved_right2nbor[0] - 3.0*conserved_rightnbor[0] + 3.0*conserved[0] - conserved_leftnbor[0]);
  double res_xmom = lambda*epsilon*(conserved_right2nbor[1] - 3.0*conserved_rightnbor[1] + 3.0*conserved[1] - conserved_leftnbor[1]);
  double res_energy = lambda*epsilon*(conserved_right2nbor[2] - 3.0*conserved_rightnbor[2] + 3.0*conserved[2] - conserved_leftnbor[2]);

  /*double res_rho = lambda*epsilon*(Field[loc+2][0] - 3.0*Field[loc+1][0] + 3.0*Field[loc][0] - Field[loc-1][0]);
  double res_vel = lambda*epsilon*(field[loc+2][1] - 3.0*field[loc+1][1] + 3.0*field[loc][1] - field[loc-1][1]);
  double res_pressure = lambda*epsilon*(field[loc+2][2] - 3.0*field[loc+1][2] + 3.0*field[loc][2] - field[loc-1][2]);
  */
  array<double,3> res = {res_continuity,res_xmom,res_energy};

  return res;
}


//-----------------------------------------------------------
void Euler1D::ComputeResidual(vector<array<double,3>> &Resid,vector<array<double,3>> &Field,vector<double> &xcoords,double &dx){ 

  //following nomenclature from class notes
  array<double,3> F_right,F_left; //left and right face spatial fluxes 
  array<double,3> D2_right,D2_left,D4_right,D4_left; //left and right face damping terms
  array<double,3> TotalF_right,TotalF_left; //left and right total fluxes (spatial + damping)
  double S; //source term (only for x-mom. eq.)
  double A_left,A_right; //Area of corresponding faces

  //Tools::print("------In Euler.ComputeResidual---------\n");
  //Tools::print("Total Size: %d\n",total_cellnum);
  for (int i=0;i<total_cellnum;i++){ //looping through all interior nodes
    if (i==0 | i==1 | i==total_cellnum-2 | i==total_cellnum-1) //skipping the ghost cell nodes
      continue;

    /*if (total_cellnum == (int)Field.size()){ //skipping the ghost cell nodes
      Tools::print("Field size w/ ghost cells size does not match with total_cellnum!\n");  
      Tools::print("Field size: %d & total_cellnum: %d\n",(int)Field.size(),total_cellnum);
      break;
    }*/

    //Tools::print("---Cell: %d---\n",i);
    //Spatial Flux Term
    //note: \arrow{F}_(i-1/2) is same as \arrow{F}_(i+1/2) of cell to the left!
    //Tools::print("Spatial Flux Energy\n");
    F_right = ComputeSpatialFlux(Field,i,i+1);
    F_left = ComputeSpatialFlux(Field,i-1,i);
    //F_left = ComputeSpatialFlux(Field,i,i-1);
    //Tools::print("F_right: %f\n",F_right[2]);
    //Tools::print("F_left: %f\n",F_left[2]);

    //Source Term (external pressure) ONLY for x-mom. eq.
    // also, area already evaluated, but may need to be multiplied dx?
    //Tools::print("Source Term\n");
    S = ComputeSourceTerm(Field,i,xcoords);
    //Tools::print("S: %f\n",S);

    //JST Damping Terms (need a D2_left flux and D2_right flux vector; similar for D4)
    //note: \arrow{D}_(i-1/2) is same as \arrow{D}_(i+1/2) of cell to the left!
    //Tools::print("Damping Flux Energy\n");
    // right face
    D2_right = Compute2ndOrderDamping(Field,i);
    D4_right = Compute4thOrderDamping(Field,i);
    // left face
    D2_left = Compute2ndOrderDamping(Field,i-1);
    D4_left = Compute4thOrderDamping(Field,i-1);

    /*
    Tools::print("D2_right: %f\n",D2_right[2]);
    Tools::print("D2_left: %f\n",D2_right[2]);
    Tools::print("D4_right: %f\n",D4_right[2]);
    Tools::print("D4_left: %f\n",D4_right[2]);
    */

    //Total Flux Terms
    //continuity
    TotalF_right[0] = F_right[0] - (D2_right[0]+D4_right[0]);
    TotalF_left[0] = F_left[0] - (D2_left[0]+D4_left[0]);
    //x-mom.
    TotalF_right[1] = F_right[1] - (D2_right[1]+D4_right[1]);
    TotalF_left[1] = F_left[1] - (D2_left[1]+D4_left[1]);
    //energy
    //Tools::print("Total Flux Energy\n");
    TotalF_right[2] = F_right[2] - (D2_right[2]+D4_right[2]);
    TotalF_left[2] = F_left[2] - (D2_left[2]+D4_left[2]);

    //Tools::print("TotalF_right: %f\n",TotalF_right[2]);
    //Tools::print("TotalF_left: %f\n",TotalF_left[2]);

    //Area Evaluations
    A_left = Tools::AreaVal(xcoords[i-2]);
    A_right = Tools::AreaVal(xcoords[i-1]);

    //Residual cal.
    //Tools::print("Cell Index: %d\n",i-2);
    
    //continuity residual (i-2 so that indexing is correct for resid spacevariable pointer) 
    //Tools::print("Continuity Residual\n");
    //Tools::print("TotalF_right:%f\n",TotalF_right[0]);
    //Tools::print("TotalF_left:%f\n",TotalF_left[0]);
    Resid[i-2][0] = (TotalF_right[0]*A_right - TotalF_left[0]*A_left);
    //Tools::print("Resid: %f\n",Resid[i-2][0]);
    
    //x-mom. Residual (w/ source term)
    //Tools::print("X-Mom. Residual\n");
    //Tools::print("TotalF_right:%f\n",TotalF_right[1]);
    //Tools::print("TotalF_left:%f\n",TotalF_left[1]);
    Resid[i-2][1] =  (TotalF_right[1]*A_right - TotalF_left[1]*A_left) - S*dx;
    //Tools::print("Resid: %f\n",Resid[i-2][1]);

    //energy Residual
    //Tools::print("Energy Residual\n");
    //Tools::print("TotalF_right:%f\n",TotalF_right[2]);
    //Tools::print("TotalF_left:%f\n",TotalF_left[2]);
    Resid[i-2][2] =  (TotalF_right[2]*A_right - TotalF_left[2]*A_left);
    //Tools::print("Resid: %f\n",Resid[i-2][2]);
    //Tools::print("---------------\n");

  }
  return;

}

//-----------------------------------------------------------
double Euler1D::GetLambdaMax(vector<array<double,3>> &Field,int loc){

  double M = GetMachNumber(Field,loc); //cell-averaged Mach number
  double a = abs(Field[loc][1]) * M; //cell-averaged speed of sound
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
