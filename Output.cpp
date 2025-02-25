//User-defined functions
#include "Output.h" 

// OUTPUT DEFINITIONS

Output::Output(array<double,3>* &f)
  : field(f) {}

//-----------------------------------------------------------

void Output::PrintResidualNorm(int &cellnum,int &n){

  if (n==1 || n==2 || n==3){ 
    array<double,3> norm = {0.0,0.0,0.0};
    for (int i=0;i<cellnum;i++){

      if (n==1){ //L1 norm case
        norm[0] += abs(field[i][0]); //density
        norm[1] += abs(field[i][1]); //velocity
        norm[2] += abs(field[i][2]); //pressure
      }
      else if (n==2){ //L2 norm case
        norm[0] += pow(field[i][0],2); //density
        norm[1] += pow(field[i][1],2); //velocity
        norm[2] += pow(field[i][2],2); //pressure
      }
      else{ //L inf. case
        if(abs(field[i][0]) > norm[0]) norm[0] = abs(field[i][0]); //density
        if(abs(field[i][1]) > norm[1]) norm[1] = abs(field[i][1]); //velocity
        if(abs(field[i][2]) > norm[2]) norm[2] = abs(field[i][2]); //pressure
      }

    }

    if (n==1){ //L2 norm case continued
      norm[0] = sqrt(norm[0]);
      norm[1] = sqrt(norm[1]);
      norm[2] = sqrt(norm[2]);
    }

    if (n==1 | n==2) //!< Printing out norms
      Tools::print("--L %d Norm Selected\n",n);
    else
      Tools::print("--L infinity Norm Selected\n");

    Tools::print("Density:%f,Velocity:%f,Pressure:%f\n",norm[0],norm[1],norm[2]);

    

  }
  

  else {
    Tools::print("Residual type unknown!\n");
    exit(0); 
  }
}

//-----------------------------------------------------------

Output::~Output(){}
