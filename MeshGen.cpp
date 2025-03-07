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
  double DArea = (dx/2.0) * abs(area_rightface - area_leftface);

  double vol = DArea * dx;
  return vol;

}

//-----------------------------------------------------------
void MeshGen1D::GenerateMesh(vector<double> &xcoords) {

  //TODO: Need to add ghost cells here
  Tools::print("1D EULER EQ. SOLVER?\n");
  double facenum = cellnumber + 1; //number of faces
  int midface_loc = facenum / 2;
  xcoords = Tools::RetrievePoints(xmin,xmax,facenum);  
  Tools::print("- Mesh Statistics:\n");
  Tools::print("XCoords: [%e,%e]\n",xmin,xmax);
  Tools::print("Number of Cells: %d\n",cellnumber);
  Tools::print("Location of middle face: %f\n",xcoords[midface_loc]);
  
/*  for (int n=0;(int)n<xcoords.size();n++) {
    Tools::print("- Mesh Statistics:\n");
    Tools::print("XCoords: [%e,%e]\n",xmin,xmax);
    Tools::print("Number of Cells: %d\n",cellnumber);
  
  }
*/
  return;

}
//-----------------------------------------------------------

MeshGen1D::~MeshGen1D(){}
