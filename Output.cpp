//User-defined functions
#include "Output.h" 

// OUTPUT DEFINITIONS

Output::Output(){}

//Output::Output(array<double,3>* &f)
  //: field(f) {}

//-----------------------------------------------------------

/*
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
*/
//-----------------------------------------------------------
void Output::DiscretizationErrorNorms(vector<array<double,4>>* &field,vector<array<double,4>>* &exact_field,vector<array<double,4>>* &errors,SpaceVariables2D* &sols,MeshGenBASE* &mesh,const char* &filename){

  for (int n=0;n<(int)field->size();n++){ //calculating errors
    for (int i=0;i<4;i++)
      (*errors)[n][i] = (*field)[n][i] - (*exact_field)[n][i];
  }
  
  //L2 Norms of Error
  array<double,4> ErrorNorms = sols->ComputeSolutionNorms(errors);

  //PRINTING TO SCREEN
  Tools::print("-------------------------\n");
  Tools::print("Discretization Error Norms\n");
  Tools::print("Density: %e\n",ErrorNorms[0]);
  Tools::print("X-Velocity: %e\n",ErrorNorms[1]);
  Tools::print("Y-Velocity: %e\n",ErrorNorms[2]);
  Tools::print("Pressure: %e\n",ErrorNorms[3]);

  //WRITING TO FILE
  std::ofstream myfile(filename,ios::app); //true for append

  if (!myfile){ //checking if file opened successfully
    cerr<<"Error: Could Not Open File "<<filename<<endl;
    return;
  }

  myfile<<"-------------------------\n";
  myfile<<"Discretization Error Norms\n";
  myfile<<"Cell Number: "<<mesh->cellnumber<<endl;
  myfile<<"I: "<<mesh->cell_imax<<"  "<<"J: "<<mesh->cell_jmax<<endl;
  myfile<<"Density: "<<ErrorNorms[0]<<endl;
  myfile<<"X-Velocity: "<<ErrorNorms[1]<<endl;
  myfile<<"Y-Velocity: "<<ErrorNorms[2]<<endl;
  myfile<<"Pressure: "<<ErrorNorms[3]<<endl;
  
  myfile.close();

  return;
}

//-----------------------------------------------------------
void Output::CalculateOrderofAccuracy(const char *filename_read,const char *filename_write){

  //NOTE: order of error norms top(coarsest mesh) - btm(finest mesh)

  ifstream myfileread(filename_read);
  ofstream myfilewrite(filename_write);

  if (!myfileread){ //Error Handling
    cerr<<"Error Opening Error Norms File to Eval. Order of Accuracy!"<<endl;
    return; 
  }

  std::string line;

  vector<double> CellNumber;
  vector<double> I,J;
  vector<double> Density;
  vector<double> XVelocity;
  vector<double> YVelocity;
  vector<double> Pressure;

  // READING DISCREIZATION ERROR FILE(.TXT)
  while (std::getline(myfileread,line)){ //reading the line as a string

    std::stringstream ss(line); //reading the line as a string
    std::string label;
    double value,i,j;

    if (line.find("Cell Number:") != std::string::npos) { //found Cell Number

      ss >> label >> label >> value; 
      CellNumber.push_back(value);
    }

    else if (line.find("I:") != std::string::npos) { //found I and J number
      ss >> label >> i >> label >> j; 
      I.push_back(i);
      J.push_back(j);
    }

    else if (line.find("Density:") != std::string::npos) { //found Density
      ss >> label >> value; 
      Density.push_back(value);
    }

    else if (line.find("X-Velocity:") != std::string::npos) { //found X-Vel.
      ss >> label >> value; 
      XVelocity.push_back(value);
    }

    else if (line.find("Y-Velocity:") != std::string::npos) { //found Y-Vel.
      ss >> label >> value; 
      YVelocity.push_back(value);
    }
    else if (line.find("Pressure:") != std::string::npos) { //found Pressure
      ss >> label >> value; 
      Pressure.push_back(value);
    }

  }

  // CALCULATING OBSERVED ORDER OF ACCURACY
  //Reference: Section 4 Slide 31 Notes to calc. order of accuracy (p)
  // NOTE: arrangement of PHat lists start from coarsest and go to finest grids!

  vector<double> PHat_density((int)CellNumber.size()-1,0); //order of accuracy value
  vector<double> PHat_xvel((int)CellNumber.size()-1,0); //order of accuracy value
  vector<double> PHat_yvel((int)CellNumber.size()-1,0); //order of accuracy value
  vector<double> PHat_pressure((int)CellNumber.size()-1,0); //order of accuracy value

  double r = I[1] / I[0]; //mesh refinement factor

  for (int n=0;n<(int)CellNumber.size()-1;n++){ 
    PHat_density[n] = (log(Density[n]/Density[n+1])) / log(r);
    PHat_xvel[n] = (log(XVelocity[n]/XVelocity[n+1])) / log(r);
    PHat_yvel[n] = (log(YVelocity[n]/YVelocity[n+1])) / log(r);
    PHat_pressure[n] = (log(Pressure[n]/Pressure[n+1])) / log(r);
 }

  vector<double> h; //grid spacing 
  h.push_back(1.0); //1st element is the finest grid
  for (int i=1;i<=(int)CellNumber.size()-2;i++) //-2 b/c not evaluating coarsest mesh
    h.push_back(h[i-1]*r); //r times the previous mesh spacing


  reverse(h.begin(),h.end()); //reversing order to match with phat calc.


  // OUTPUTTING OBSERVED ORDER OF ACCURACY IN MATLAB-FRIENDLY FORMAT
  if (!myfilewrite){ //Error Handling
    cerr<<"Error Opening Output for Observed Order of Accuracy File!"<<endl;
    return; 
  }


  for (int n=0;n<(int)h.size();n++)
    myfilewrite<<h[n]<<" "<<PHat_density[n]<<" "<<PHat_xvel[n]<<" "<<PHat_yvel[n]<<" "<<PHat_pressure[n]<<endl;


  //.DAT FORMAT
  /*
  myfilewrite<<"variables= \"grid spacing(h)\" \"Phat(density)\" \"Phat(velocity)\"  \"Phat(Pressure)\""<<endl;

  myfilewrite<<"zone T= "<<"\""<<0<<"\""<<endl;
  myfilewrite<<"I="<<(int)h.size()<<endl;
  myfilewrite<<"DATAPACKING=POINT"<<endl;
  myfilewrite<<"DT=(DOUBLE DOUBLE DOUBLE DOUBLE )"<<endl;

  for (int n=0;n<(int)h.size();n++)
    myfilewrite<<h[n]<<"  "<<PHat_density[n]<<"  "<<PHat_velocity[n]<<"  "<<PHat_pressure[n]<<"  "<<endl;

  */
  
  myfilewrite.close();

  return;

}

//-----------------------------------------------------------
void Output::OutputPrimitiveVariables(vector<array<double,4>>* &field,string filename,bool cond,int iter,vector<double> &xcoords,vector<double> &ycoords,int cell_number,int imax,int jmax){

  std::ofstream myfile(filename,(cond==true) ? ios::app : ios::out); //true for append
  //myfile.open(filename);

  if (!myfile){ //checking if file opened successfully
    cerr<<"Error: Could Not Open File "<<filename<<endl;
    return;
  }

  if (cond==false){ //start of .dat file -- printing initial parameters
    myfile<<"TITLE = \" 2D Field Solutions \""<<endl;
    myfile<<"VARIABLES = \"X\",\"Y\",\"Rho\",\"U\",\"V\",\"P\""<<endl;
    //myfile<<"variables= \"cell index\" \"rho(kg/m^3)\" \"u(m/s)\"  \"Press(N/m^2)\" \"Mach\" \"Xcoords\""<<endl;
  }

//Repeat the following each time you want to write out the solution
/*
write(40,*) 'zone T="',num_iter,'" '
write(40,*) 'I=',imax
write(40,*) 'DATAPACKING=POINT'
write(40,*) 'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE
& DOUBLE DOUBLE )'
  */

  myfile<<"ZONE T="<<"\""<<iter<<"\""<<endl; //Now adding zone specific info.
  myfile<<"I="<<imax<<", "<<"J="<<jmax<<endl;
  myfile<<"DATAPACKING=BLOCK"<<endl;
  //myfile<<"DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE )"<<endl;
  myfile<<"VARLOCATION=([3-6]=CELLCENTERED)"<<endl; //-> tells Tecplot this is cell-centered val (must be size (imax-1)*(jmax-1) size


  // Saving all primitive variables in their own corresponding vector
  vector<double> all_rho,all_u,all_v,all_p;

  for (int i=0;i<cell_number;i++){
    all_rho.push_back((*field)[i][0]);
    all_u.push_back((*field)[i][1]);
    all_v.push_back((*field)[i][2]);
    all_p.push_back((*field)[i][3]);
  }

  //TODO: Reset count after each write of data type
  int count = 0;
  // Writing Xcoords
  for (int n=0;n<(int)xcoords.size();n++){
    count++;
    myfile<<std::setw(15)<<xcoords[n];
    if (count % 4 == 0)
      myfile<<endl;
  }

  // Writing Ycoords
  for (int n=0;n<(int)ycoords.size();n++){
    count++;
    myfile<<std::setw(15)<<ycoords[n];
    if (count % 4 == 0)
      myfile<<endl;
  }


  // Writing Rho
  for (int n=0;n<(int)all_rho.size();n++){
    count++;
    myfile<<std::setw(15)<<all_rho[n];
    if (count % 4 == 0)
      myfile<<endl;
  }
  
  // Writing U 
  for (int n=0;n<(int)all_u.size();n++){
    count++;
    myfile<<std::setw(15)<<all_u[n];
    if (count % 4 == 0)
      myfile<<endl;
  }

  // Writing V 
  for (int n=0;n<(int)all_v.size();n++){
    count++;
    myfile<<std::setw(15)<<all_v[n];
    if (count % 4 == 0)
      myfile<<endl;
  }

  // Writing P
  for (int n=0;n<(int)all_p.size();n++){
    count++;
    myfile<<std::setw(15)<<all_p[n];
    if (count % 4 == 0)
      myfile<<endl;
  }
  myfile<<endl;
  myfile.close(); //closing file writing to it
  //myfile.flush();

  return;
}

//-----------------------------------------------------------
void Output::OutputManufacturedSourceTerms(vector<array<double,4>>* &field,string filename,bool cond,int iter,vector<double> &xcoords,vector<double> &ycoords,int cell_number,int imax,int jmax){

  std::ofstream myfile(filename,(cond==true) ? ios::app : ios::out); //true for append
  //myfile.open(filename);

  if (!myfile){ //checking if file opened successfully
    cerr<<"Error: Could Not Open File "<<filename<<endl;
    return;
  }

  if (cond==false){ //start of .dat file -- printing initial parameters
    myfile<<"TITLE = \" 2D Field Solutions \""<<endl;
    myfile<<"VARIABLES = \"X\",\"Y\",\"Continuity\",\"X-Momentum\",\"Y-Momentum\",\"Energy\""<<endl;
  }

  myfile<<"ZONE T="<<"\""<<iter<<"\""<<endl; //Now adding zone specific info.
  myfile<<"I="<<imax<<", "<<"J="<<jmax<<endl;
  myfile<<"DATAPACKING=BLOCK"<<endl;
  //myfile<<"DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE )"<<endl;
  myfile<<"VARLOCATION=([3-6]=CELLCENTERED)"<<endl; //-> tells Tecplot this is cell-centered val (must be size (imax-1)*(jmax-1) size


  // Saving all primitive variables in their own corresponding vector
  vector<double> cont,xmom,ymom,energy;

  for (int i=0;i<cell_number;i++){
    cont.push_back((*field)[i][0]);
    xmom.push_back((*field)[i][1]);
    ymom.push_back((*field)[i][2]);
    energy.push_back((*field)[i][3]);
  }

  int count = 0;
  // Writing Xcoords
  for (int n=0;n<(int)xcoords.size();n++){
    count++;
    myfile<<std::setw(15)<<xcoords[n];
    if (count % 4 == 0)
      myfile<<endl;
  }

  // Writing Ycoords
  for (int n=0;n<(int)ycoords.size();n++){
    count++;
    myfile<<std::setw(15)<<ycoords[n];
    if (count % 4 == 0)
      myfile<<endl;
  }


  // Writing Continuity
  for (int n=0;n<(int)cont.size();n++){
    count++;
    myfile<<std::setw(15)<<cont[n];
    if (count % 4 == 0)
      myfile<<endl;
  }
  
  // Writing X-Momentum
  for (int n=0;n<(int)xmom.size();n++){
    count++;
    myfile<<std::setw(15)<<xmom[n];
    if (count % 4 == 0)
      myfile<<endl;
  }

  // Writing Y-Momentum
  for (int n=0;n<(int)ymom.size();n++){
    count++;
    myfile<<std::setw(15)<<ymom[n];
    if (count % 4 == 0)
      myfile<<endl;
  }

  // Writing Energy
  for (int n=0;n<(int)energy.size();n++){
    count++;
    myfile<<std::setw(15)<<energy[n];
    if (count % 4 == 0)
      myfile<<endl;
  }

  myfile.close(); //closing file writing to it
  //myfile.flush();

  return;
}
//-----------------------------------------------------------
void Output::OutputGhostCoords(string filename,vector<double> &xcoords,vector<double> &ycoords,int Nx,int Ny){

  std::ofstream myfile(filename); //true for append
  //myfile.open(filename);


  if (!myfile){ //checking if file opened successfully
    cerr<<"Error: Could Not Open File "<<filename<<endl;
    return;
  }

  //Title
  myfile<<"TITLE = \" 2D Ghost Cells \""<<endl;
  myfile<<"VARIABLES = \"X\",\"Y\""<<endl;
  myfile<<"ZONE I="<<Nx<<", "<<"J="<<Ny<<endl;
  myfile<<"DATAPACKING=BLOCK"<<endl;

  //Coords
  int count = 0;
  for (int n=0;n<(int)xcoords.size();n++){ //xcoords
    count++;
    myfile<<std::setw(15)<<xcoords[n];
    if (count % 4 == 0)
      myfile<<endl;
  }

  for (int n=0;n<(int)ycoords.size();n++){ //ycoords
    count++;
    myfile<<std::setw(15)<<ycoords[n];
    if (count % 4 == 0)
      myfile<<endl;
  }

  return;

}
//-----------------------------------------------------------
void Output::OutputGhostCells(vector<array<double,4>>* &ghost_cell,string filename,vector<double> &xcoords,vector<double> &ycoords,vector<double> &ghost_xcoords,vector<double> &ghost_ycoords,int Nx,int Ny,int ghost_Nx,int ghost_Ny,int side){

  std::ofstream myfile(filename); //true for append
  //myfile.open(filename);


  if (!myfile){ //checking if file opened successfully
    cerr<<"Error: Could Not Open File "<<filename<<endl;
    return;
  }

  //Title
  myfile<<"TITLE = \" 2D Ghost Cells \""<<endl;
  myfile<<"VARIABLES = \"X\",\"Y\",\"Rho\",\"U\",\"V\",\"P\""<<endl;
  myfile<<"ZONE I="<<ghost_Nx<<", "<<"J="<<ghost_Ny<<endl;
  myfile<<"DATAPACKING=BLOCK"<<endl;
  myfile<<"VARLOCATION=([3-6]=CELLCENTERED)"<<endl; //-> tells Tecplot this is cell-centered val (must be size (imax-1)*(jmax-1) size

  //Grouping interior coords w/ ghost coords
  vector<double> total_xcoords,total_ycoords;
  //interior
  if (side==0){ //for top ghost cells
    for (int i=0;i<Nx;i++){
      total_xcoords.push_back(xcoords[i+((Ny-1)*Nx)]);
      total_ycoords.push_back(ycoords[i+((Ny-1)*Nx)]);
    }

  }
  else if (side==1){ //for btm ghost cells
    for (int i=0;i<Nx;i++){
      total_xcoords.push_back(xcoords[i+(0*Nx)]);
      total_ycoords.push_back(ycoords[i+(0*Nx)]);
    }


  }
  else if (side==2){ //for left ghost cells
    for (int j=0;j<Ny;j++){
      total_xcoords.push_back(xcoords[0+(j*Nx)]);
      total_ycoords.push_back(ycoords[0+(j*Nx)]);
    }


  }
  else if (side==3){ //for right ghost cells
    for (int j=0;j<Ny;j++){
      total_xcoords.push_back(xcoords[(Nx-1)+(j*Nx)]);
      total_ycoords.push_back(ycoords[(Nx-1)+(j*Nx)]);
    }

  }

  else //error handling
    cerr<<"ERROR:Invalid side selection for outputting ghost cells"<<endl;
  //ghost coords
  for (int n=0;n<(int)ghost_xcoords.size();n++){
    total_xcoords.push_back(ghost_xcoords[n]);
    total_ycoords.push_back(ghost_ycoords[n]);
  }


  //Writing coords to file
  int count = 0;
  for (int n=0;n<(int)total_xcoords.size();n++){ //xcoords
    count++;
    myfile<<std::setw(15)<<total_xcoords[n];
    if (count % 4 == 0)
      myfile<<endl;
  }

  for (int n=0;n<(int)total_ycoords.size();n++){ //ycoords
    count++;
    myfile<<std::setw(15)<<total_ycoords[n];
    if (count % 4 == 0)
      myfile<<endl;
  }

  //Writing primitive variables
  vector<double> all_rho,all_u,all_v,all_p;

  for (int i=0;i<(int)(*ghost_cell).size();i++){
    all_rho.push_back((*ghost_cell)[i][0]);
    all_u.push_back((*ghost_cell)[i][1]);
    all_v.push_back((*ghost_cell)[i][2]);
    all_p.push_back((*ghost_cell)[i][3]);
  }

  // Writing Rho
  for (int n=0;n<(int)all_rho.size();n++){
    count++;
    myfile<<std::setw(15)<<all_rho[n];
    if (count % 4 == 0)
      myfile<<endl;
  }
  
  // Writing U 
  for (int n=0;n<(int)all_u.size();n++){
    count++;
    myfile<<std::setw(15)<<all_u[n];
    if (count % 4 == 0)
      myfile<<endl;
  }

  // Writing V 
  for (int n=0;n<(int)all_v.size();n++){
    count++;
    myfile<<std::setw(15)<<all_v[n];
    if (count % 4 == 0)
      myfile<<endl;
  }

  // Writing P
  for (int n=0;n<(int)all_p.size();n++){
    count++;
    myfile<<std::setw(15)<<all_p[n];
    if (count % 4 == 0)
      myfile<<endl;
  }

  myfile.close(); //closing file writing to it


  return;

}
//-----------------------------------------------------------
void Output::WritePVDFile(const char* &filename,vector<string> &iter_visuals){
ofstream myfile(filename);

  if (!myfile){ //checking if file opened successfully
    cerr<<"Error: Could Not Open File "<<filename<<endl;
    return;
  }

  //TITLE
  myfile<<"<?xml version=\"1.0\"?>"<<endl;
  myfile<<"<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
  myfile<<"  <Collection>"<<endl;
  //Writing file names
  for (int n=0;n<(int)iter_visuals.size();n++)
    myfile<<"    <DataSet timestep=\""<<iter_visuals[n]<<"\""<<" part=\"0\""<<" file=\"Iteration"<<iter_visuals[n]<<".dat\"/>"<<endl;

  myfile<<"  </Collection>"<<endl;
  myfile<<"</VTKFile>"<<endl;
/*
<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">
  <Collection>
    <DataSet timestep="0" part="0" file="solution_0.vtu"/>
    <DataSet timestep="1" part="0" file="solution_1.vtu"/>
    <DataSet timestep="2" part="0" file="solution_2.vtu"/>
    <DataSet timestep="3" part="0" file="solution_3.vtu"/>
    <!-- Add more timesteps here -->
  </Collection>
</VTKFile>
*/

  return;

}
//-----------------------------------------------------------
string Output::zeroPad(int number, int padWidth) {
  ostringstream ss;
  ss << std::setw(padWidth) << std::setfill('0') << number;
  return ss.str();
}
//-----------------------------------------------------------

Output::~Output(){}
