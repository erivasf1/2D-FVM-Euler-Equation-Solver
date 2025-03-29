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

  //TODO: Need to add ghost cells here
  //Tools::print("1D EULER EQ. SOLVER\n");
  double facenum = cellnumber + 1; //number of faces
  int midface_loc = facenum / 2;
  xcoords = Tools::RetrievePoints(xmin,xmax,facenum);  
  //Tools::print("- Mesh Statistics:\n");
  //Tools::print("XCoords: [%e,%e]\n",xmin,xmax);
  //Tools::print("Number of Cells: %d\n",cellnumber);
  //Tools::print("Location of middle face: %f\n",xcoords[midface_loc]);
  
/*  for (int n=0;(int)n<xcoords.size();n++) {
    Tools::print("- Mesh Statistics:\n");
    Tools::print("XCoords: [%e,%e]\n",xmin,xmax);
    Tools::print("Number of Cells: %d\n",cellnumber);
  
  }
*/
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
