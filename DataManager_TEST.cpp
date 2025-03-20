//User-defined functions
#include "DataManager_TEST.h" 

// SPACEVARIABLE1D DEFINITIONS

SpaceVariables1D::SpaceVariables1D()
{}

//---------------------------------------------------------
array<double,3> SpaceVariables1D::ComputeSolutionNorms(vector<array<double,3>> &Resid){

  array<double,3> norm{0.0,0.0,0.0};

  double imax = (double)Resid.size(); //number of interior nodes 

  //using L2 norm  
  //Tools::print("Calculating Global Norm\n");
  for (int i=0;i<(int)Resid.size();i++){
    //Tools::print("Cell number: %d\n",i);
    //continuity
    //Tools::print("continuity res.: %e\n",Resid[i][0]);
    norm[0]+= pow(Resid[i][0],2);
    //x-mom
    //Tools::print("x-mom. res.: %e\n",Resid[i][1]);
    norm[1]+= pow(Resid[i][1],2);
    //energy
    //Tools::print("energy res.: %e\n",Resid[i][2]);
    norm[2]+= pow(Resid[i][2],2);

  }
  norm[0] = sqrt(norm[0]/imax); //normalizing by imax
  norm[1] = sqrt(norm[1]/imax);
  norm[2] = sqrt(norm[2]/imax);

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
    myfile<<"Cell#"<<"  "<<"Density(kg/m^3)"<<"  "<<"Velocity(m/s)"<<"  "<<"Pressure(Pa)"<<"  "<<"Mach Number"<<endl;
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
void SpaceVariables1D::AllOutputPrimitiveVariables(vector<array<double,3>> &Field,Euler1D &Euler,string filename,bool cond,int iter,vector<double> &xcoords){

  cell_number = Field.size()-4; //number of interior cells
  std::ofstream myfile(filename,(cond==true) ? ios::app : ios::out); //true for append
  //myfile.open(filename);

  if (!myfile){ //checking if file opened successfully
    cerr<<"Error: Could Not Open File "<<filename<<endl;
    return;
  }

  if (cond==false)
    myfile<<"variables= \"cell index\" \"rho(kg/m^3)\" \"u(m/s)\"  \"Press(N/m^2)\" \"Mach\" \"Xcoords\""<<endl;

//Repeat the following each time you want to write out the solution
/*
write(40,*) 'zone T="',num_iter,'" '
write(40,*) 'I=',imax
write(40,*) 'DATAPACKING=POINT'
write(40,*) 'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE
& DOUBLE DOUBLE )'
  */

  myfile<<"zone T= "<<"\""<<iter<<"\""<<endl;
  myfile<<"I="<<cell_number<<endl;
  myfile<<"DATAPACKING=POINT"<<endl;
  myfile<<"DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE )"<<endl;

  double M;

  for (int n=0;n<cell_number;n++){
    //Computing Mach Number for checking the initial conditions
    M = Euler.GetMachNumber(Field,n+2);
    myfile<<n+1<<"  ";
    myfile<<Field[n][0]<<"  "<<Field[n][1]<<"  "<<Field[n][2]<<"  "<<M<<"  "<<xcoords[n]<<endl;

  }

  
  
  myfile.close(); //closing file writing to it
  //myfile.flush();

  return;
}

//---------------------------------------------------------
void SpaceVariables1D::OutputLocalResiduals(vector<array<double,3>> &Resid,const char *filename){


  ofstream myfile;
  myfile.open(filename);
  //ofstream myfile(filename);

  if (!myfile.is_open()){ //checking if file opened successfully
    cerr<<"Error: Could Not Open File "<<filename<<endl;
    return;
  }
  myfile<<"Local Residuals"<<endl;
  myfile<<"Cell"<<"  "<<"Continuity"<<"  "<<"Momentum"<<"  "<<"Energy"<<endl;
  for (int n=0;n<(int)Resid.size();n++){
    myfile<<n<<"  "<<Resid[n][0]<<"  "<<Resid[n][1]<<"  "<<Resid[n][2]<<endl;
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
void SpaceVariables1D::ComputeCellAveragedSol(vector<array<double,3>> &SolFace,vector<array<double,3>> &SolCell,vector<double> &xcoords,double dx){

  //SolFace has one more element than SolCell, due to more faces than cells
  // using trapezoidal to approximate integral
  int rght_face,lft_face; //indexes
  double A_right,A_left;
  for (int n=0;n<(int)SolCell.size();n++){ //computing sol. at each cell
    lft_face = n;
    rght_face = n + 1; 
    A_right = Tools::AreaVal(xcoords[rght_face]);
    A_left = Tools::AreaVal(xcoords[lft_face]);

    for (int i=0;i<3;i++){ //computing cell average for each primitive variable
      SolCell[n][i] = A_right*SolFace[rght_face][i] + A_left*SolFace[lft_face][i];
      SolCell[n][i] /= A_right + A_left;
      //SolCell[n][i] = 0.5*(A_right*SolFace[rght_face][i] + A_left*SolFace[lft_face][i]);
      //SolCell[n][i] /= 0.5*(A_right + A_left)*dx;
    }

  }

  return;
}

//---------------------------------------------------------
SpaceVariables1D::~SpaceVariables1D(){}
