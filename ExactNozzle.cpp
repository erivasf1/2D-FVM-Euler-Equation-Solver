//User-defined functions
#include "ExactNozzle.h" 

// TOOLS DEFINITIONS
void Tools::print(const char format[],...) {

  va_list Argp;
  va_start(Argp, format);
  vprintf(format, Argp);
  va_end(Argp);
  return;

}

void Tools::printWithPrecision(int precision,const char* format,...) {
  //Code generated from ChatGPT
  std::ostringstream output;
    output << std::fixed << std::setprecision(precision);

    // Initialize the variadic argument list
    va_list args;
    va_start(args, format);

    const char* ptr = format;
    while (*ptr != '\0') {
        if (*ptr == '%' && *(ptr + 1) == 'f') { // Handle floating-point formatting
            double value = va_arg(args, double);
            output << value;
            ptr++; // Skip 'f'
        } else {
            output << *ptr;
        }
        ptr++;
    }

    va_end(args);

    // Print the final formatted string
    std::cout << output.str() << std::endl;

}





vector<double> Tools::RetrievePoints(double xmin,double xmax,int pt_num){

  vector<double> pts;
  //pts.push_back(xmin); -- used in HW1

  double dx = abs(xmax-xmin) / (pt_num-1.0); 
  for (int k=0;k<pt_num;k++)
    pts.push_back(xmin + k * dx); 
    //pts.push_back(xmin+(dx*(k+1))); -- used in HW1
   
  //debug:
  //print("Points:%f,%f,%f,%f,%f\n",pts[0],pts[1],pts[2],pts[3],pts[4]);

  return pts;
}

double Tools::AreaVal(double x){

  //A(x) = 0.2 + 0.4[1+sin(pi(x-0.5)]
  double PI = M_PI;
  double res = 0.2 + 0.4*(1+sin(PI*(x-0.5)));
  return res;
}

// SUPERSONICNOZZLE DEFINITIONS
SuperSonicNozzle::SuperSonicNozzle(double &a,double &b,double &c,double &d,bool &e) //constructor
    : area(a), area_star(b), stag_pressure(c), stag_temp(d), cond(e){}


SuperSonicNozzle::~SuperSonicNozzle(){} //destructor


double SuperSonicNozzle::GetPhi(double M) {

  double phi = 2.0/(gamma+1.0);
  phi *= 1.0 + ((gamma-1.0)/2.0) * pow(M,2.0);
  //phi *= (1.0 + ((gamma-1.0)/2.0)) * pow(M,2.0);
  return phi;
}

double SuperSonicNozzle::GetF(double Phi,double ABar,double M) {

  double f = pow(Phi,(gamma+1.0)/(gamma-1.0));
  f -= pow(ABar*M,2.0);
  return f;

}

double SuperSonicNozzle::GetFPrime(double Phi,double ABar,double M) {

  double fprime = 2.0 * M;
  fprime *= pow(Phi,2.0/(gamma-1.0)) - pow(ABar,2.0);
  return fprime;

}


double SuperSonicNozzle::ComputeMachNumber(){
//F[m] = [2/gamma+1(1+(gamma-1)/2)M^2]^[(gamma+1)/(gamma-1)]
// F'[M] = 2M[Phi^(2/gamma-1)-\Bar{A}^2] & Phi = [2/gamma+1(1+(gamma-1)/2)M^2]

  double M0,M1; //old (initial guess) and converged sol., respectively
  double ABar = area/area_star; //ratio of present area to throat area
  //print("Abar:%f\n",Abar);
  double F,Fprime; //fcn. and fcn. derivative
  double phi;
  double resid; // residual

  //debug:
  print("area: %f and tol.: %f\n",area,tol);
  
  M0 = (cond == true) ? 0.2:3.0; // initial guess of 0.2 for subsonic & 3.0 for supersonic region
  //print("M0 initial guess:%f\n",M0); //for debugging
  M1 = M0; //initial guess & resid.
  resid = 1.0e5;  

  //Computing initial residual
  double phi_init = GetPhi(M0);
  double F_init = GetF(phi_init,ABar,M0);
  
  //Employing Newton's Method to solve nonlinear equation
  for (int i=0;i<maxiter;i++) {
    if (resid <= tol) break;

    phi = GetPhi(M1);
    F = GetF(phi,ABar,M1);
    Fprime = GetFPrime(phi,ABar,M1);
    M1 -= F/Fprime;
    

    //checking residual
    phi = GetPhi(M1);
    F = GetF(phi,ABar,M1);
    resid = abs(F/F_init);
    //resid = abs(F);
    //F_init = F;

    //debug mode
    //print("M1: %f & residual:%f\n",M1,resid);
  }
  return M1;

}

void SuperSonicNozzle::ComputeExactSol(array<double,3> &sol){

  //vector<double> sol;
  double T,P,Rho,V;
  double M,Psi;
  double R = Ru/MolMass; //specific gas constant

  //print("Point: %f and Area: %f\n",pts[n],area);
  M = ComputeMachNumber();
  //print("Mach Number: %f\n",M);
  Psi = 1.0+(gamma-1.0)/2.0 * pow(M,2.0);

  //Temperature
  // Psi = 1 + (gamma-1/2)M^2
   T = stag_temp/Psi;
   //print("Temperature ratio: %f\n",T/stag_temp);
  
  //Pressure
  // P = P0/[Psi^(gamma/gamma-1)] 
   P = stag_pressure / pow(Psi,gamma/(gamma-1.0));
   //print("Pressure ratio: %f\n",P/stag_pressure);

  //Density
  // Rho = P / (RT)
   Rho = P / (R*T);
   double Rho_t = stag_pressure / (R*stag_temp);
   //print("Density ratio: %f\n",Rho/Rho_t);
  
  //Velocity
  // V = sqrt(gamma*R*T)*M
   V = abs(M * sqrt(gamma*R*T));

  //Appending flow quantities sol. vector of point
  /*sol.push_back(Rho);
  sol.push_back(V);
  sol.push_back(P);
  sol.push_back(M);
  */
  sol[0] = Rho;
  sol[1] = V;
  sol[2] = P;
    //debug:
    //print("Point: %f\n",pts[n]);
    //print("Rho[kg/m^3] = %f,Velocity[m/s] = %f,Pressure[kPa] = %f,Mach Number = %f\n",Rho,V,P,M);

  
  
  //PrintResults(sol);
}

void SuperSonicNozzle::PrintResults(vector<double> &sol) {

  /*print("---------\n");
  print("RESULTS\n");
  print("---------\n");*/


  //print("Rho[kg/m^3] = %f,Velocity[m/s] = %f,Pressure[kPa] = %f,Mach Number = %f\n",sol[0],sol[1],sol[2],sol[3]);
  printWithPrecision(14,"Rho[kg/m^3] = %f,Velocity[m/s] = %f,Pressure[kPa] = %f,Mach Number = %f\n",sol[0],sol[1],sol[2],sol[3]);
  print("-------\n");

}



