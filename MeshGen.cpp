//User-defined functions
#include "MeshGen.h" 

// MESHGEN1D DEFINITIONS

MeshGen1D::MeshGen1D(double &a,double &b,int &c)
  : xmin(a), xmax(b), cellnumber(c) {}

//-----------------------------------------------------------
double MeshGen1D::GetCellVolume(int &loc,double &dx,vector<double> &xcoords){

  //using trapezoidal rule since areas are stored at cell faces
  double area_leftface = Tools::AreaVal(xcoords[loc]); 
  double area_rightface = Tools::AreaVal(xcoords[loc+1]); 
  double area_characteristic = 0.5*(area_leftface+area_rightface); //avg. of left and right face cell area
  //double DArea = (dx/2.0) * abs(area_rightface - area_leftface);
  //previous:double DArea = (dx/2.0) * abs(area_rightface - area_leftface);

  double vol = area_characteristic * dx;
  return vol;

}

//-----------------------------------------------------------
void MeshGen1D::GenerateMesh(vector<double> &xcoords) {

  int facenum = cellnumber + 1; //number of faces
  xcoords = Tools::RetrievePoints(xmin,xmax,facenum);  

  return;

}
//-----------------------------------------------------------
void MeshGen1D::OutputNozzleAreas(vector<double> &xcoords,const char *filename){

  //Computing Areas
  vector<double> Areas((int)xcoords.size());
  for (int n=0;n<(int)xcoords.size();n++){
    Areas[n] = Tools::AreaVal(xcoords[n]);
  }


  //Printing out in filename
  ofstream myfile;
  myfile.open(filename);
  //ofstream myfile(filename);

  if (!myfile.is_open()){ //checking if file opened successfully
    cerr<<"Error: Could Not Open File "<<filename<<endl;
    return;
  }
    myfile<<"Areas"<<endl;
    myfile<<"XCoord"<<"  "<<"Area"<<endl;
    for (int i=0;i<(int)xcoords.size();i++){
  
      myfile<<xcoords[i]<<"  "<<Areas[i]<<endl;

    }

  
  myfile.close(); //closing file writing to it
  //myfile.flush();

  return;


}

//-----------------------------------------------------------
MeshGen1D::~MeshGen1D(){}

//-----------------------------------------------------------
MeshGen2D::MeshGen2D(const char* name) 
  : filename(name) {

  ReadMeshFile(); //extract nodal coords. of mesh
}
  
//-----------------------------------------------------------
void MeshGen2D::ReadMeshFile(){

  ifstream myfileread(filename);
  
  if (!myfileread) { //Error Handling
    cerr<<"No Mesh File Provided!"<<endl;
    //cerr<<"Error Opening Mesh File "<<filename<<" !"<<endl; 
    return;
  }

  std::string line;
  std::getline(myfileread,line); //!< skips 1st line

  if (std::getline(myfileread,line)){ //!< calling 2nd line
    std::istringstream iss(line); //converts string into stream
    iss >> imax >> jmax >> kmax; //setting values of 2nd line to imax,jmax,kmax, respectively
  }
  else{
    cerr<<"Mesh File does not contain a second line!"<<endl; 
    return;
  }

  double val; //value of i,j,&k indices
  int total_pts = imax*jmax;
  
  cellnumber = (imax-1) * (jmax-1); //1 more faces than each dir.
  //int pt_ct = 1;
  
  for (int j=0;j<3;j++){
    for (int i=0;i<total_pts;i++){
      myfileread >> val;
      if (j==0)
        xcoords.push_back(val);  
      else if (j==1)
        ycoords.push_back(val);
      else
        zcoords.push_back(val);
    }
  }
  //skipping first line for now
  //2nd line: 1st int refers to imax and 2nd int refers to jmax

}
//-----------------------------------------------------------
void MeshGen2D::OutputMesh(){

  ifstream myfileread(filename); 
  ofstream myfilewrite("Mesh.dat",ios::out);

  if (!myfilewrite){
    cerr<<"Error: Could not Open \"Mesh.dat File\"!"<<endl;
    return;
  }

  //Writing header to Mesh Output File
  myfilewrite<<"TITLE = \" 2D Structured Mesh \""<<endl;
  myfilewrite<<"VARIABLES = \"X\",\"Y\",\"Z\""<<endl;
  myfilewrite<<"ZONE I="<<imax<<", J="<<jmax<<", K="<<kmax<<" DATAPACKING=BLOCK"<<endl;
  
  
  //Reading vals. of input mesh file & writing to output mesh file
  double val;
  string line;
  std::getline(myfileread,line);//!< skipping 1st and 2nd line
  std::getline(myfileread,line);

  //std::istringstream iss(line);
  int count = 0;
  while(myfileread>>val){
    count++;
    myfilewrite<<std::setw(15)<<val;
    if (count % 4 == 0)
      myfilewrite<<endl;
 
  }
    
  //myfilewrite.close(); 

}

//-----------------------------------------------------------
MeshGen2D::~MeshGen2D(){}
//-----------------------------------------------------------
