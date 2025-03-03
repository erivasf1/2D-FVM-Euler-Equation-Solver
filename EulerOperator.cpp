//User-defined functions
#include "EulerOperator.h" 

// EULER1D DEFINITIONS

Euler1D::Euler1D(int& cellnum,double &P0,double &T0,double &g)
  : interior_cellnum(cellnum), stag_pressure(P0) , stag_temperature(T0) 
  , gamma(g) {}


//-----------------------------------------------------------
Euler1D::Euler1D() {}


//-----------------------------------------------------------
void Euler1D::SetInitialConditions(array<double,3>* &field,vector<double> &xcoords){

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
    field[i][2] = pow(psi,gamma/(gamma-1.0));
    field[i][2] = stag_pressure / field[i][2]; 
    //Tools::print("pressure: %f\n",field[i][2]);
    
    //density calc.
    T = stag_temperature / psi; // local temperature
    field[i][0] = field[i][2] / (R*T); 
    //Tools::print("density: %f\n",field[i][0]);

    //velocity calc.
    a = sqrt(gamma*R*T); //local speed of sound
    field[i][1] = abs(M*a);
    //Tools::print("velocity: %f\n",field[i][1]);
     
  }    
  
  return;

}
//-----------------------------------------------------------
void Euler1D::SetBoundaryConditions(vector<array<double,3>> &Field,array<double,3>* &field,bool &cond){

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
  field = Field.data(); //reassign pointer to new Field (w/ ghost cells)
  //Tools::print("total number of cells after set BC:%d\n",total_cellnum);

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
    

    //Tools::print("B.C.\n");
    //Tools::print("Cell index:%d\n",i);
    //Tools::print("Density: %f, Velocity:%f,Pressure:%f\n",field[i][0],field[i][1],field[i][2]);
   }

  return;
}

//-----------------------------------------------------------
void Euler1D::ComputeOutflowBoundaryConditions(array<double,3>* &field,bool& cond){

  //TODO: Figure out why [i-1] node is not retrieving the node of interior node

  if (cond == false){ //supersonic case 
    //using simple extrapolation from class notes section 3 slide 36
    for (int i=total_cellnum-2;i<total_cellnum;i++){ //for both outflow ghost cells
      field[i][0] = 2.0*field[i-1][0] - field[i-2][0]; //density
      field[i][1] = 2.0*field[i-1][1] - field[i-2][1]; //velocity
      field[i][2] = 2.0*field[i-1][2] - field[i-2][2]; //pressure

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
    field[G3][2] = 2.0*Pb - field[G3-1][2]; //pressure at G3
    field[G3+1][2] = 2.0*field[G3][2] - field[G3-2][2]; //extrapolated pressure for G4

    for (int i=total_cellnum-2;i<total_cellnum;i++){ //reg. extrapolation
      field[i][0] = 2.0*field[i-1][0] - field[i-2][0]; //density
      field[i][1] = 2.0*field[i-1][1] - field[i-2][1]; //velocity
    }


  }


  return;
}

//-----------------------------------------------------------
array<double,3> Euler1D::ComputeSpatialFlux(array<double,3>* &field,int loc,int nbor){

  // Conversion into conservative variables(cv)
  //Rho * U conserved variable
  double cv1 = field[loc][0]*field[loc][1]; 
  double cv1_nbor = field[nbor][0]*field[nbor][1]; 
  //Rho * U^2 + P conserved variable
  double cv2 = field[loc][0]*pow(field[loc][1],2) + field[loc][2];
  double cv2_nbor = field[nbor][0]*pow(field[nbor][1],2) + field[nbor][2];
  //Rho*U*h_t conserved variable
  //Tools::print("rho: %f\n",field[loc][0]);
  double h_t = gamma/(gamma-1.0) * (field[loc][2]/field[loc][0]) + (pow(field[loc][1],2)/2.0);
  //Tools::print("h_t[i]:%f\n",h_t);
  double cv3 = field[loc][0]*field[loc][1]*h_t;
  //Tools::print("cv3[i]:%f\n",cv3);

  h_t = gamma/(gamma-1.0) * (field[nbor][2]/field[nbor][0]) + (pow(field[nbor][1],2)/2.0);
  //Tools::print("h_t[nbor]:%f\n",h_t);
  double cv3_nbor = field[nbor][0]*field[nbor][1]*h_t;
  //Tools::print("cv3[nbor]:%f\n",cv3_nbor);

  // Value at interface interpolated using central quadrature
  double flux_continuity = (cv1+cv1_nbor) / 2.0;
  double flux_xmom = (cv2+cv2_nbor) / 2.0;
  double flux_energy = (cv3+cv3_nbor) / 2.0;

  array<double,3> flux = {flux_continuity,flux_xmom,flux_energy};
  return flux; 

}

//-----------------------------------------------------------
double Euler1D::ComputeSourceTerm(array<double,3>* &field,int &loc,vector<double> &xcoords) {

  //Get area for i+1/2 and i-1/2 locations (refer to read me for data indexing)
  double A_rface = Tools::AreaVal(xcoords[loc-1]);
  double A_lface = Tools::AreaVal(xcoords[loc-2]);
  // pressure at i = loc
  double p = field[loc][2];
  // source = P_i * (A_i+1/2 - A_i-1/2) / dx * dx
  double res = p * (A_rface - A_lface);

  return res;

}

//-----------------------------------------------------------
double Euler1D::GetEpsilon2(array<double,3>* &field,int &loc) {

  double Nu = GetNu(field,loc);
  double Nuleft = GetNu(field,loc-1);
  double Nuright = GetNu(field,loc+1);
  double Nuright2 = GetNu(field,loc+2);

  double kappa2 = 0.35; //typically from 1/4<kappa2<1/2
  
  double max = std::max({Nu,Nuleft,Nuright,Nuright2}); //acquring max nu
  double res = kappa2 * max;
  return res;

}


//-----------------------------------------------------------
double Euler1D::GetNu(array<double,3>* &field,int loc){

  //Now changed to look at neighbors, originally everything was "loc"
  double res = field[loc+1][2] - (2.0*field[loc][2]) +  field[loc-1][2];
  res /= field[loc+1][2] + (2.0*field[loc][2]) +  field[loc-1][2];
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

  double M,a; //cell-averaged Mach number and speed of sound, respectively
  //\bar{lambda_i} at current cell
  //double a = sqrt(gamma*R*T); //TODO:speed of sound (define a fcn. for this)
  // T
  M = GetMachNumber(field,loc); 
  a = field[loc][1] * M;
  double lambda_i = abs(field[loc][1]) + a;

  //\bar{lambda_i+1} at neighboring cell to the right
  M = GetMachNumber(field,loc+1); //cell-averaged Mach Number
  a = field[loc+1][1] * M;
  double lambda_iright = abs(field[loc+1][1]) + a;

  double res = (lambda_i + lambda_iright) / 2.0;
  return res;

}


//-----------------------------------------------------------
array<double,3> Euler1D::Compute2ndOrderDamping(array<double,3>* &field,int loc){

  //following Roy's class notes nomenclature (section 3 slide 31)
  //returns solely \arrow{d^2} vector!
  //look into applying damping terms to boundary?
  // seems correct for now
  double lambda = GetLambda(field,loc); //at cell face (i+1/2)
  double epsilon = GetEpsilon2(field,loc); //sensor for detecting shocks (will have to tweak the constant later)

  double res_rho = lambda*epsilon*(field[loc+1][0] - field[loc][0]);
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
  double epsilon = GetEpsilon4(field,loc); //looks for gradients that cause odd-even decoupling (will have to tweak the constant later)

  double res_rho = lambda*epsilon*(field[loc+2][0] - 3.0*field[loc+1][0] + 3.0*field[loc][0] - field[loc-1][0]);
  double res_vel = lambda*epsilon*(field[loc+2][1] - 3.0*field[loc+1][1] + 3.0*field[loc][1] - field[loc-1][1]);
  double res_pressure = lambda*epsilon*(field[loc+2][2] - 3.0*field[loc+1][2] + 3.0*field[loc][2] - field[loc-1][2]);

  array<double,3> res = {res_rho,res_vel,res_pressure};

  return res;
}


//-----------------------------------------------------------
void Euler1D::ComputeResidual(array<double,3>* &resid,array<double,3>* &field,vector<double> &xcoords,double &dx){ //TODO

  //following nomenclature from class notes
  array<double,3> F_right,F_left; //left and right face spatial fluxes 
  array<double,3> D2_right,D2_left,D4_right,D4_left; //left and right face damping terms
  array<double,3> TotalF_right,TotalF_left; //left and right total fluxes (spatial + damping)
  double S; //source term (only for x-mom. eq.)
  double A_left,A_right; //Area of corresponding faces

  for (int i=0;i<total_cellnum;i++){ //looping through all interior nodes
    if (i==0 | i==1 | i==total_cellnum-2 | i==total_cellnum-1) //skipping the ghost cell nodes
      continue;

    //Tools::print("---Cell: %d\n---",i);
    //Spatial Flux Term
    //note: \arrow{F}_(i-1/2) is same as \arrow{F}_(i+1/2) of cell to the left!
    //Tools::print("Spatial Flux Energy\n");
    F_right = ComputeSpatialFlux(field,i,i+1);
    F_left = ComputeSpatialFlux(field,i-1,i);
    //F_left = ComputeSpatialFlux(field,i,i-1);
    //Tools::print("F_right: %f\n",F_right[2]);
    //Tools::print("F_left: %f\n",F_left[2]);

    //Source Term (external pressure) ONLY for x-mom. eq.
    // also, area already evaluated, but may need to be multiplied dx?
    //Tools::print("Source Term\n";
    S = ComputeSourceTerm(field,i,xcoords);
    //Tools::print("S: %f\n",S);

    //JST Damping Terms (need a D2_left flux and D2_right flux vector; similar for D4)
    //note: \arrow{D}_(i-1/2) is same as \arrow{D}_(i+1/2) of cell to the left!
    //Tools::print("Damping Flux Energy\n");
    // right face
    D2_right = Compute2ndOrderDamping(field,i);
    D4_right = Compute4thOrderDamping(field,i);
    // left face
    D2_left = Compute2ndOrderDamping(field,i-1);
    D4_left = Compute4thOrderDamping(field,i-1);

    //Tools::print("D2_right: %f\n",D2_right[2]);
    //Tools::print("D2_left: %f\n",D2_right[2]);
    //Tools::print("D4_right: %f\n",D4_right[2]);
    //Tools::print("D4_left: %f\n",D4_right[2]);

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
    
    //continuity residual (i-2 so that indexing is correct for resid spacevariable pointer) 
    resid[i-2][0] = (TotalF_right[0]*A_right - TotalF_left[0]*A_left);
    
    //x-mom. residual (w/ source term)
    resid[i-2][1] =  (TotalF_right[1]*A_right - TotalF_left[1]*A_left) - S*dx;

    //energy residual
    //Tools::print("Energy Residual\n");
    resid[i-2][2] =  (TotalF_right[2]*A_right - TotalF_left[2]*A_left);
    //Tools::print("resid: %f\n",resid[i-2][2]);

  }
  return;

}

//-----------------------------------------------------------
double Euler1D::GetLambdaMax(array<double,3>* &field,int &loc){

  double M = GetMachNumber(field,loc); //cell-averaged Mach number
  double a = abs(field[loc][1]) * M; //cell-averaged speed of sound
  //double a = absfield[loc][1] * M; //cell-averaged speed of sound
  double lambda_max = abs(field[loc][1]) + a;

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
