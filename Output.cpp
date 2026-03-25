//User-defined functions
#include "Output.h" 

// OUTPUT DEFINITIONS

//-----------------------------------------------------------
Output::Output()
{}
//-----------------------------------------------------------
Output::Output(string &fresults,string &fghostcells ,string &fresids, string &fmmserror,int iout, MeshGenBASE* m ) : results_prefix(fresults), ghostcells_prefix(fghostcells), resids_prefix(fresids), mms_error_prefix(fmmserror), iterout(iout), mesh(m)
{}

//-----------------------------------------------------------
void Output::OutputResidualNorms(const char* &filename,int iter,array<double,4> ResidualNorms){

  std::ofstream myfile(filename,ios::app);
  //Error handling
  if (!myfile){ //checking if file opened successfully
    cerr<<"Error: Could Not Open File "<<filename<<endl;
    return;
  }

  //WRITING RESIDUAL NORMS FOR PLOTTING
  //col1: iter, col2: continuity, col3: x-mom, col4: y-mom, col5: energy
  myfile<<iter<<"  "<<ResidualNorms[0]<<"  "<<ResidualNorms[1]<<"  "
                                            <<ResidualNorms[2]<<"  "
                                            <<ResidualNorms[3]<<endl;

  myfile.close();

  return;

}
//-----------------------------------------------------------
void Output::DiscretizationErrorNorms(vector<array<double,4>>* &field,vector<array<double,4>>* &exact_field,vector<array<double,4>>* &errors,SpaceVariables2D* &sols,const char* &filename){

  for (int n=0;n<(int)field->size();n++){ //calculating errors
    for (int i=0;i<4;i++)
      (*errors)[n][i] = (*field)[n][i] - (*exact_field)[n][i];
  }
  
  //L1 and L2 Norms of Error
  array<double,4> L1ErrorNorms = sols->ComputeL1SolutionNorms(errors);
  array<double,4> L2ErrorNorms = sols->ComputeL2SolutionNorms(errors);

  //PRINTING TO SCREEN
  Tools::print("-------------------------\n");
  Tools::print("Discretization Error Norms\n");
  Tools::print("Density: L1 = %e | L2 = %e\n",L1ErrorNorms[0],L2ErrorNorms[0]);
  Tools::print("X-Velocity: L1 = %e | L2 = %e\n",L1ErrorNorms[1],L2ErrorNorms[1]);
  Tools::print("Y-Velocity: L1 = %e | L2 = %e\n",L1ErrorNorms[2],L2ErrorNorms[2]);
  Tools::print("Pressure: L1 = %e | L2 = %e\n",L1ErrorNorms[3],L2ErrorNorms[3]);

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
  myfile<<"Density: "<<L1ErrorNorms[0]<<" "<<L2ErrorNorms[0]<<endl;
  myfile<<"X-Velocity: "<<L1ErrorNorms[1]<<" "<<L2ErrorNorms[1]<<endl;
  myfile<<"Y-Velocity: "<<L1ErrorNorms[2]<<" "<<L2ErrorNorms[2]<<endl;
  myfile<<"Pressure: "<<L1ErrorNorms[3]<<" "<<L2ErrorNorms[3]<<endl;
  
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
  vector<double> Density_L1, Density_L2;
  vector<double> XVelocity_L1, XVelocity_L2;
  vector<double> YVelocity_L1, YVelocity_L2;
  vector<double> Pressure_L1, Pressure_L2;

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
      ss >> label >> i >> j; 
      Density_L1.push_back(i);
      Density_L2.push_back(j);
    }

    else if (line.find("X-Velocity:") != std::string::npos) { //found X-Vel.
      ss >> label >> i >> j; 
      XVelocity_L1.push_back(i);
      XVelocity_L2.push_back(j);
    }

    else if (line.find("Y-Velocity:") != std::string::npos) { //found Y-Vel.
      ss >> label >> i >> j; 
      YVelocity_L1.push_back(i);
      YVelocity_L2.push_back(j);
    }
    else if (line.find("Pressure:") != std::string::npos) { //found Pressure
      ss >> label >> i >> j; 
      Pressure_L1.push_back(i);
      Pressure_L2.push_back(j);
    }

  }

  // CALCULATING OBSERVED ORDER OF ACCURACY
  //Reference: Section 4 Slide 31 Notes to calc. order of accuracy (p)
  // NOTE: arrangement of PHat lists start from coarsest and go to finest grids!

  vector<double> PHat_density_L1((int)CellNumber.size()-1,0); //order of accuracy value
  vector<double> PHat_density_L2 = PHat_density_L1;

  vector<double> PHat_xvel_L1((int)CellNumber.size()-1,0); //order of accuracy value
  vector<double> PHat_xvel_L2 = PHat_density_L1;

  vector<double> PHat_yvel_L1((int)CellNumber.size()-1,0); //order of accuracy value
  vector<double> PHat_yvel_L2 = PHat_density_L1;
  vector<double> PHat_pressure_L1((int)CellNumber.size()-1,0); //order of accuracy value
  vector<double> PHat_pressure_L2 = PHat_density_L1;

  double r = I[1] / I[0]; //mesh refinement factor

  //Observed order of accuracy calc. for L1
  for (int n=0;n<(int)CellNumber.size()-1;n++){ 
    PHat_density_L1[n] = (log(Density_L1[n]/Density_L1[n+1])) / log(r);
    PHat_xvel_L1[n] = (log(XVelocity_L1[n]/XVelocity_L1[n+1])) / log(r);
    PHat_yvel_L1[n] = (log(YVelocity_L1[n]/YVelocity_L1[n+1])) / log(r);
    PHat_pressure_L1[n] = (log(Pressure_L1[n]/Pressure_L1[n+1])) / log(r);
 }

  //Observed order of accuracy calc. for L2
  for (int n=0;n<(int)CellNumber.size()-1;n++){ 
    PHat_density_L2[n] = (log(Density_L2[n]/Density_L2[n+1])) / log(r);
    PHat_xvel_L2[n] = (log(XVelocity_L2[n]/XVelocity_L2[n+1])) / log(r);
    PHat_yvel_L2[n] = (log(YVelocity_L2[n]/YVelocity_L2[n+1])) / log(r);
    PHat_pressure_L2[n] = (log(Pressure_L2[n]/Pressure_L2[n+1])) / log(r);
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
    myfilewrite<<h[n]<<" "<<PHat_density_L1[n]<<" "<<PHat_density_L2[n]<<" "
                          <<PHat_xvel_L1[n]<<" "<<PHat_xvel_L2[n]<<" "
                          <<PHat_yvel_L1[n]<<" "<<PHat_yvel_L2[n]<<" "
                        <<PHat_pressure_L1[n]<<" "<<PHat_pressure_L2[n]<<endl;


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
void Output::OutputPrimitiveVariables_VTS(string filename,vector<array<double,4>>* &field){


  std::ofstream myfile(filename); 

  if (!myfile){ //checking if file opened successfully
    cerr<<"Error: Could Not Open File "<<filename<<endl;
    return;
  }

  int NI = mesh->cell_imax; int NJ = mesh->cell_jmax; //# of cells in each dir.
  int NPI = mesh->cell_imax + 1; int NPJ = mesh->cell_jmax + 1; //# of pts. in each dir.

  //-----------------------------HEADER-------------------------------
  myfile<<"<?xml version=\"1.0\"?>"<<endl;
  myfile<<"<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
  myfile<<"  <StructuredGrid WholeExtent=\"0 "<<NI<<" 0 "<<NJ<<" 0 1\">"<<endl;
  myfile<<"    <Piece Extent=\"0 "<<NI<<" 0 "<<NJ<<" 0 1\">"<<endl;


  //-----------------------------POINT DATA-------------------------------
  myfile<<"      <Points>"<<endl;
  myfile<<"        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"<<endl;;
  int id; //for indexing pts.
  for (int k=0;k<2;k++){
    for (int j=0;j<NPJ;j++){
      for (int i=0;i<NPI;i++){
       id = j * NPI + i;

        myfile<<"          "<<mesh->xcoords[id]<<" "<<mesh->ycoords[id]<<" "<<k<<endl;

      }
    }
  }
  myfile<<"        </DataArray>"<<endl;
  myfile<<"      </Points>"<<endl;

  //-----------------------------CELL DATA-------------------------------
  myfile<<"      <CellData Scalars=\"Density Pressure\" Vectors=\"Velocity\">"<<endl;

  //DENSITY
  myfile<<"        <DataArray type=\"Float64\" Name=\"Density\" NumberOfComponents=\"1\" format=\"ascii\">"<<endl;;
  
  for (int n=0;n<mesh->cellnumber;n++){
      myfile << (*field)[n][0] << " ";
      if ((n + 1) % 3 == 0) myfile << "\n";
  }
  if (mesh->cellnumber % 3 != 0) myfile << "\n";

  myfile << "        </DataArray>"<<endl;
  
  //PRESSURE
  myfile<<"        <DataArray type=\"Float64\" Name=\"Pressure\" NumberOfComponents=\"1\" format=\"ascii\">"<<endl;;
  
  for (int n=0;n<mesh->cellnumber;n++){
      myfile << (*field)[n][3] << " ";
      if ((n + 1) % 3 == 0) myfile << "\n";
  }
  if (mesh->cellnumber % 3 != 0) myfile << "\n";

  myfile << "        </DataArray>"<<endl;

  //VELOCITY
  myfile<<"        <DataArray type=\"Float64\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">"<<endl;;
  
  for (int n=0;n<mesh->cellnumber;n++){
      myfile << (*field)[n][1] << " " << (*field)[n][2]<<" "<<0.0<<endl;
      if ((n + 1) % 3 == 0) myfile << "\n";
  }
  if (mesh->cellnumber % 3 != 0) myfile << "\n";

  myfile << "        </DataArray>"<<endl;
  myfile << "      </CellData>"<<endl;

  //-----------------------------FOOTER-------------------------------
  myfile << "    </Piece>\n";
  myfile << "  </StructuredGrid>\n";
  myfile << "</VTKFile>\n";


}
//-----------------------------------------------------------
void Output::OutputManufacturedSourceTerms(string filename,vector<array<double,4>>* &field){


  std::ofstream myfile(filename); 

  if (!myfile){ //checking if file opened successfully
    cerr<<"Error: Could Not Open File "<<filename<<endl;
    return;
  }

  int NI = mesh->cell_imax; int NJ = mesh->cell_jmax; //# of cells in each dir.
  int NPI = mesh->cell_imax + 1; int NPJ = mesh->cell_jmax + 1; //# of pts. in each dir.

  //-----------------------------HEADER-------------------------------
  myfile<<"<?xml version=\"1.0\"?>"<<endl;
  myfile<<"<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
  myfile<<"  <StructuredGrid WholeExtent=\"0 "<<NI<<" 0 "<<NJ<<" 0 1\">"<<endl;
  myfile<<"    <Piece Extent=\"0 "<<NI<<" 0 "<<NJ<<" 0 1\">"<<endl;


  //-----------------------------POINT DATA-------------------------------
  myfile<<"      <Points>"<<endl;
  myfile<<"        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"<<endl;;
  int id; //for indexing pts.
  for (int k=0;k<2;k++){
    for (int j=0;j<NPJ;j++){
      for (int i=0;i<NPI;i++){
       id = j * NPI + i;

        myfile<<"          "<<mesh->xcoords[id]<<" "<<mesh->ycoords[id]<<" "<<k<<endl;

      }
    }
  }
  myfile<<"        </DataArray>"<<endl;
  myfile<<"      </Points>"<<endl;

  //-----------------------------CELL DATA-------------------------------
  myfile<<"      <CellData Scalars=\"Continuity Energy\" Vectors=\"Momentum\">"<<endl;

  //DENSITY
  myfile<<"        <DataArray type=\"Float64\" Name=\"Continuity\" NumberOfComponents=\"1\" format=\"ascii\">"<<endl;;
  
  for (int n=0;n<mesh->cellnumber;n++){
      myfile << (*field)[n][0] << " ";
      if ((n + 1) % 3 == 0) myfile << "\n";
  }
  if (mesh->cellnumber % 3 != 0) myfile << "\n";

  myfile << "        </DataArray>"<<endl;
  
  //PRESSURE
  myfile<<"        <DataArray type=\"Float64\" Name=\"Energy\" NumberOfComponents=\"1\" format=\"ascii\">"<<endl;;
  
  for (int n=0;n<mesh->cellnumber;n++){
      myfile << (*field)[n][3] << " ";
      if ((n + 1) % 3 == 0) myfile << "\n";
  }
  if (mesh->cellnumber % 3 != 0) myfile << "\n";

  myfile << "        </DataArray>"<<endl;

  //VELOCITY
  myfile<<"        <DataArray type=\"Float64\" Name=\"Momentum\" NumberOfComponents=\"3\" format=\"ascii\">"<<endl;;
  
  for (int n=0;n<mesh->cellnumber;n++){
      myfile << (*field)[n][1] << " " << (*field)[n][2]<<" "<<0.0<<endl;
      if ((n + 1) % 3 == 0) myfile << "\n";
  }
  if (mesh->cellnumber % 3 != 0) myfile << "\n";

  myfile << "        </DataArray>"<<endl;
  myfile << "      </CellData>"<<endl;

  //-----------------------------FOOTER-------------------------------
  myfile << "    </Piece>\n";
  myfile << "  </StructuredGrid>\n";
  myfile << "</VTKFile>\n";



}
//-----------------------------------------------------------
void Output::WriteAllGhostCellSolutions(const char* &filename_btm,const char* &filename_top,const char* &filename_left,const char* &filename_right){

  WriteGhostCellSolution_PVD(filename_btm,0); //bottom cells
  WriteGhostCellSolution_PVD(filename_top,1); //top cells
  WriteGhostCellSolution_PVD(filename_left,2); //left cells
  WriteGhostCellSolution_PVD(filename_right,3); //right cells

  return;
}
//-----------------------------------------------------------
void Output::WriteGhostCellSolution_PVD(const char* &filename,int tag){

  std::ofstream myfile(filename); 

  if (!myfile){ //checking if file opened successfully
    cerr<<"Error: Could Not Open File "<<filename<<endl;
    return;
  }

  int NI,NJ,NPI,NPJ;

  if (tag == 0 || tag == 1){ //top & btm cells
    NI = mesh->cell_imax; NJ = 2; //# of cells in each dir.
  }

  else if (tag == 2 || tag == 3) { //right & left cells
    NI = 2; NJ = mesh->cell_jmax; //# of cells in each dir.
  }

  else{
    cerr<<"Error: Unkown ID specified in WriteGhostCellSolution fcn."<<endl;
    return;
  }

  NPI = NI + 1; NPJ = NJ + 1; //# of pts. in each dir.
  int cellnum = NI * NJ; //# of cells

  //-----------------------------HEADER-------------------------------
  myfile<<"<?xml version=\"1.0\"?>"<<endl;
  myfile<<"<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
  myfile<<"  <StructuredGrid WholeExtent=\"0 "<<NI<<" 0 "<<NJ<<" 0 1\">"<<endl;
  myfile<<"    <Piece Extent=\"0 "<<NI<<" 0 "<<NJ<<" 0 1\">"<<endl;


  //-----------------------------POINT DATA-------------------------------
  myfile<<"      <Points>"<<endl;
  myfile<<"        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"<<endl;;
  int id{0}; //for indexing pts.
  for (int k=0;k<2;k++){
    for (int j=0;j<NPJ;j++){
      for (int i=0;i<NPI;i++){
       id = j * NPI + i;

        if (tag == 0) //btm cells 
          myfile<<"          "<<mesh->btm_xcoords[id]<<" "<<mesh->btm_ycoords[id]<<" "<<k<<endl;
        else if (tag == 1) //top cells 
          myfile<<"          "<<mesh->top_xcoords[id]<<" "<<mesh->top_ycoords[id]<<" "<<k<<endl;
        else if (tag == 2) //left cells
          myfile<<"          "<<mesh->left_xcoords[id]<<" "<<mesh->left_ycoords[id]<<" "<<k<<endl;
        else //right cells
          myfile<<"          "<<mesh->right_xcoords[id]<<" "<<mesh->right_ycoords[id]<<" "<<k<<endl;

      }
    }
  }
  myfile<<"        </DataArray>"<<endl;
  myfile<<"      </Points>"<<endl;

  //-----------------------------CELL DATA-------------------------------
  myfile<<"      <CellData Scalars=\"Density Pressure\" Vectors=\"Velocity\">"<<endl;

  //DENSITY
  myfile<<"        <DataArray type=\"Float64\" Name=\"Density\" NumberOfComponents=\"1\" format=\"ascii\">"<<endl;;
  
  for (int n=0;n<cellnum;n++){
    if (tag == 0)
      myfile << mesh->btm_cells[n][0] << " ";
    else if (tag == 1)
      myfile << mesh->top_cells[n][0] << " ";
    else if (tag == 2)
      myfile << mesh->left_cells[n][0] << " ";
    else 
      myfile << mesh->right_cells[n][0] << " ";

    if ((n + 1) % 3 == 0) myfile << "\n";
  }
  if (cellnum % 3 != 0) myfile << "\n";

  myfile << "        </DataArray>"<<endl;
  
  //PRESSURE
  myfile<<"        <DataArray type=\"Float64\" Name=\"Pressure\" NumberOfComponents=\"1\" format=\"ascii\">"<<endl;;
  
  for (int n=0;n<cellnum;n++){
    if (tag == 0)
      myfile << mesh->btm_cells[n][3] << " ";
    else if (tag == 1)
      myfile << mesh->top_cells[n][3] << " ";
    else if (tag == 2)
      myfile << mesh->left_cells[n][3] << " ";
    else 
      myfile << mesh->right_cells[n][3] << " ";

    if ((n + 1) % 3 == 0) myfile << "\n";
  }
  if (cellnum % 3 != 0) myfile << "\n";


  myfile << "        </DataArray>"<<endl;

  //VELOCITY
  myfile<<"        <DataArray type=\"Float64\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">"<<endl;;
  
  for (int n=0;n<cellnum;n++){
    if (tag == 0){
      myfile << mesh->btm_cells[n][1] << " " << mesh->btm_cells[n][2]<<" "<<0.0<<endl;
    }
    else if (tag == 1){
      myfile << mesh->top_cells[n][1] << " " << mesh->top_cells[n][2]<<" "<<0.0<<endl;
    }
    else if (tag == 2){
      myfile << mesh->left_cells[n][1] << " " << mesh->left_cells[n][2]<<" "<<0.0<<endl;
    }
    else {
      myfile << mesh->right_cells[n][1] << " " << mesh->right_cells[n][2]<<" "<<0.0<<endl;
    }

      if ((n + 1) % 3 == 0) myfile << "\n";
  }


  if (cellnum % 3 != 0) myfile << "\n";

  myfile << "        </DataArray>"<<endl;
  myfile << "      </CellData>"<<endl;

  //-----------------------------FOOTER-------------------------------
  myfile << "    </Piece>\n";
  myfile << "  </StructuredGrid>\n";
  myfile << "</VTKFile>\n";


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
    myfile<<"    <DataSet timestep=\""<<iter_visuals[n]<<"\""<<" part=\"0\""<<" file=\"Iteration_"<<iter_visuals[n]<<".vts\"/>"<<endl;

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
void Output::WriteGhostCellPVDFile(const char* &filename,vector<string> &iter_visuals){

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
  for (int n=0;n<(int)iter_visuals.size();n++){
    myfile<<"    <DataSet timestep=\""<<iter_visuals[n]<<"\""<<" group="<<"\"\""<<" part=\"0\""<<" file=\"BTM_Iteration_"<<iter_visuals[n]<<".vts\"/>"<<endl;

    myfile<<"    <DataSet timestep=\""<<iter_visuals[n]<<"\""<<" group="<<"\"\""<<" part=\"1\""<<" file=\"TOP_Iteration_"<<iter_visuals[n]<<".vts\"/>"<<endl;

    myfile<<"    <DataSet timestep=\""<<iter_visuals[n]<<"\""<<" group="<<"\"\""<<" part=\"2\""<<" file=\"LEFT_Iteration_"<<iter_visuals[n]<<".vts\"/>"<<endl;

    myfile<<"    <DataSet timestep=\""<<iter_visuals[n]<<"\""<<" group="<<"\"\""<<" part=\"3\""<<" file=\"RIGHT_Iteration_"<<iter_visuals[n]<<".vts\"/>"<<endl;


  }

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
void Output::WriteSolutions(int iter,vector<array<double,4>>* &field,vector<array<double,4>>* &resid,[[maybe_unused]] vector<array<double,4>>* &field_ms,[[maybe_unused]] vector<array<double,4>>* &field_ms_error,array<double,4> &ResidualNorms,EulerBASE* &euler,int scenario,bool &resid_stall,vector<string> &iter_visuals_primitive,vector<string> &iter_visuals_resid,[[maybe_unused]] vector<string> &iter_visuals_MMSerror){

  if (iter % iterout == 0) {
  
    //Outputting primitive variables in field
    string full_filename = results_prefix;
    string it = zeroPad(iter,4);

    full_filename += "Iteration_";
    full_filename += it;
    full_filename += ".vts";
    const char* full_filename_iter = full_filename.c_str();

    OutputPrimitiveVariables_VTS(full_filename_iter,field);
    iter_visuals_primitive.push_back(it); //saving iter value to iter_visuals list

    //Outputting ghost cell solutions
    full_filename = ghostcells_prefix;

    string full_filename_btm = full_filename + "BTM_Iteration_";
    string full_filename_top = full_filename + "TOP_Iteration_";
    string full_filename_left = full_filename + "LEFT_Iteration_";
    string full_filename_right = full_filename + "RIGHT_Iteration_";

    full_filename_btm += it; full_filename_top += it;
    full_filename_left += it; full_filename_right += it;

    full_filename_btm += ".vts"; full_filename_top += ".vts";
    full_filename_left += ".vts"; full_filename_right += ".vts";

    const char* full_filename_iter_btm = full_filename_btm.c_str();
    const char* full_filename_iter_top = full_filename_top.c_str();
    const char* full_filename_iter_left = full_filename_left.c_str();
    const char* full_filename_iter_right = full_filename_right.c_str();

    WriteAllGhostCellSolutions(full_filename_iter_btm,full_filename_iter_top,full_filename_iter_left,full_filename_iter_right);

  
    //Outputting residuals
    //error->OutputResidualNorms(resid_file,iter,ResidualNorms); //saving in norms file
    full_filename = resids_prefix;

    full_filename += "Iteration_";
    full_filename += it;
    full_filename += ".vts";
    full_filename_iter = full_filename.c_str();

    OutputPrimitiveVariables_VTS(full_filename_iter,resid);
    iter_visuals_resid.push_back(it);

    //Outputting MMS (if applicable)
    if (scenario == 3 || scenario == 4){
      full_filename = mms_error_prefix;

      full_filename += "Iteration_";
      full_filename += it;
      full_filename += ".vts";
      full_filename_iter = full_filename.c_str();

      euler->ComputeMSError(field_ms_error,field,field_ms);
      OutputPrimitiveVariables_VTS(full_filename_iter,field_ms_error);
      iter_visuals_MMSerror.push_back(it);
    }
 

    //! PRINTING RESIDUAL NORMS TO SCREEN
    Tools::print("------Iteration #: %d----------\n",iter);
    Tools::print("Continuity:%e\nX-Momentum:%e\nY-Momentum:%e\nEnergy:%e\n",ResidualNorms[0],ResidualNorms[1],ResidualNorms[2],ResidualNorms[3]);

    Tools::print("Epsilon: %e\n",euler->epsilon);

    if (resid_stall == true)
      Tools::print("Residuals are detected to be stalled!\n"); //printing message is temp. for now

    // Writing Residuals history to "SolResids.txt" file
    //myresids<<iter<<"  "<<ResidualNorms[0]<<"  "<<ResidualNorms[1]<<"  "<<ResidualNorms[2]<<endl;

  }

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
