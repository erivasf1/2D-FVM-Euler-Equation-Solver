//User-defined functions
#include "DataManager_TEST.h" 

// SPACEVARIABLE1D DEFINITIONS

SpaceVariables1D::SpaceVariables1D()
{}

//---------------------------------------------------------
array<double,3> SpaceVariables1D::ComputeSolutionNorms(vector<array<double,3>> &Field){

  array<double,3> norm{0.0,0.0,0.0};

  //using L2 norm  
  Tools::print("Calculating Global Norm\n");
  for (int i=0;i<(int)Field.size();i++){
    //Tools::print("Cell number: %d\n",i);
    //continuity
    //Tools::print("continuity res.: %e\n",Field[i][0]);
    norm[0]+= pow(Field[i][0],2);
    //x-mom
    //Tools::print("x-mom. res.: %e\n",Field[i][1]);
    norm[1]+= pow(Field[i][1],2);
    //energy
    //Tools::print("energy res.: %e\n",Field[i][2]);
    norm[2]+= pow(Field[i][2],2);

  }
  norm[0] = sqrt(norm[0]);
  norm[1] = sqrt(norm[1]);
  norm[2] = sqrt(norm[2]);

  //Tools::print("Continuity res.: %e\n",norm[0]);
  //Tools::print("X-Momentum res.: %e\n",norm[1]);
  //Tools::print("Energy res.: %e\n",norm[2]);

  return norm;
}

//---------------------------------------------------------
void SpaceVariables1D::OutputPrimitiveVariables(vector<array<double,3>> &Field,Euler1D &Euler,const char *filename){


  ofstream myfile;
  myfile.open(filename);
  //ofstream myfile(filename);

  if (!myfile.is_open()){ //checking if file opened successfully
    cerr<<"Error: Could Not Open File "<<filename<<endl;
    return;
  }
    double M; //temp. variable for Mach number
    myfile<<"Primitive Variable Solutions (Including Ghost Cells)"<<endl;
    myfile<<"Cell#"<<"  "<<"Density"<<"  "<<"Velocity"<<"  "<<"Pressure"<<"  "<<"Mach Number"<<endl;
    for (int i=0;i<(int)Field.size();i++){
  
      //Computing Mach Number for checking the initial conditions
      M = Euler.GetMachNumber(Field,i);

      myfile<<i<<"  "<<Field[i][0]<<"  "<<Field[i][1]<<"  "<<Field[i][2]<<"  "<<M<<endl;

    }

  
  myfile.close(); //closing file writing to it
  //myfile.flush();

return;
}

//---------------------------------------------------------
void SpaceVariables1D::OutputResidualTerms(array<double,3> F_right,array<double,3> F_left,double S,array<double,3> D2_right,array<double,3> D2_left,array<double,3> D4_right,array<double,3> D4_left,const char* filename){

  ofstream myfile;
  myfile.open(filename);
  vector<string> names{"continuity","x-momentum","energy"};
  //ofstream myfile(filename);

  if (!myfile.is_open()){ //checking if file opened successfully
    cerr<<"Error: Could Not Open File "<<filename<<endl;
    return;
  }

   array<double,3> Source{0.0,S,0.0}; //array vector of sorce terms

    for (int n=0;n<3;n++){
      myfile<<names[0]<<" Residual Fluxes"<<endl;
      myfile<<"F_right: "<<F_right[n]<<"  "<<"F_left: "<<F_left[n]<<"  "<<"S: "<<Source[n]<<"  "<<"D2_right"<<D2_right[n]<<"  "<<"D2_left"<<D2_left[n]<<"  "<<"D4_right"<<D4_right[n]<<"  "<<"D4_left"<<D4_left[n]<<endl;
  
    }
  
  myfile.close(); //closing file writing to it
  //myfile.flush();
   



  return;
}

//---------------------------------------------------------
SpaceVariables1D::~SpaceVariables1D(){}
