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
void Output::DiscretizationErrorNorms(vector<array<double,3>> &Field,vector<array<double,3>> &ExactField,vector<array<double,3>> &Errors,SpaceVariables1D Sols){

  for (int n=0;n<(int)Field.size();n++){ //calculating errors
    for (int i=0;i<3;i++)
      Errors[n][i] = Field[n][i] - ExactField[n][i];
  }
  
  //L2 Norms of Error
  array<double,3> ErrorNorms = Sols.ComputeSolutionNorms(Errors);

  Tools::print("-------------------------\n");
  Tools::print("Discretization Error Norms\n");
  Tools::print("Density: %e\n",ErrorNorms[0]);
  Tools::print("Velocity: %e\n",ErrorNorms[1]);
  Tools::print("Pressure: %e\n",ErrorNorms[2]);

  return;
}

//-----------------------------------------------------------
void Output::CalculateOrderofAccuracy(const char *filename_read,const char *filename_write){

  ifstream myfileread(filename_read);
  ofstream myfilewrite(filename_write);

  if (!myfileread){ //Error Handling
    cerr<<"Error Opening Discretization Error Norms File!"<<endl;
    return; 
  }

  std::string line;

  vector<double> CellSize;
  vector<double> Density;
  vector<double> Velocity;
  vector<double> Pressure;

  // Reading Discreization Error File(.txt)
  while (std::getline(myfileread,line)){ //reading the line as a string

    std::stringstream ss(line); //reading the line as a string
    std::string label;
    double value;

    if (line.find("Cell Size:") != std::string::npos) { //found Cell Size

      ss >> label >> label >> value; 
      CellSize.push_back(static_cast<int>(value));
    }

    else if (line.find("Density:") != std::string::npos) { //found Cell Size
      ss >> label >> value; 
      Density.push_back(value);
    }

    else if (line.find("Velocity:") != std::string::npos) { //found Cell Size
      ss >> label >> value; 
      Velocity.push_back(value);
    }
    else if (line.find("Pressure:") != std::string::npos) { //found Cell Size
      ss >> label >> value; 
      Pressure.push_back(value);
    }

  }

  // Calculating Observed Order of Accuracy
  //using Section 4 Slide 31 Notes to calc. order of accuracy (p)
  // NOTE: arrangement of PHat lists start from coarsest and go to finest grids

  vector<double> PHat_density((int)CellSize.size()-1,0); //order of accuracy value
  vector<double> PHat_velocity((int)CellSize.size()-1,0); //order of accuracy value
  vector<double> PHat_pressure((int)CellSize.size()-1,0); //order of accuracy value

  double r = 2.0; //mesh refinement factor
  for (int n=0;n<(int)CellSize.size()-1;n++){ 
    PHat_density[n] = (log(Density[n]/Density[n+1])) / log(r);
    PHat_velocity[n] = (log(Velocity[n]/Velocity[n+1])) / log(r);
    PHat_pressure[n] = (log(Pressure[n]/Pressure[n+1])) / log(r);
 }

  vector<double> h; //grid spacing 
  h.push_back(1.0); //1st element is the finest grid
  for (int i=1;i<=(int)CellSize.size()-2;i++) //-2 b/c not evaluating coarsest mesh
    h.push_back(h[i-1]*r); //r times the previous mesh spacing


  reverse(h.begin(),h.end()); //reversing order to match with phat calc.


  // Outputting Observed Order of Accuracy in .dat format
  if (!myfilewrite){ //Error Handling
    cerr<<"Error Opening Output for Observed Order of Accuracy File!"<<endl;
    return; 
  }

  myfilewrite<<"variables= \"grid spacing(h)\" \"Phat(density)\" \"Phat(velocity)\"  \"Phat(Pressure)\""<<endl;

  myfilewrite<<"zone T= "<<"\""<<0<<"\""<<endl;
  myfilewrite<<"I="<<(int)h.size()<<endl;
  myfilewrite<<"DATAPACKING=POINT"<<endl;
  myfilewrite<<"DT=(DOUBLE DOUBLE DOUBLE DOUBLE )"<<endl;

  for (int n=0;n<(int)h.size();n++)
    myfilewrite<<h[n]<<"  "<<PHat_density[n]<<"  "<<PHat_velocity[n]<<"  "<<PHat_pressure[n]<<"  "<<endl;

  
  myfilewrite.close();

  return;

}

//-----------------------------------------------------------
void Output::ConvertToDatFile(const char*filename_read,const char *filename_write){


  return;
}

//-----------------------------------------------------------

Output::~Output(){}
