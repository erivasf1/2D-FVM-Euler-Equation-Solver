//User-defined functions
#include "MeshGen.h" 

// MESHGEN1D DEFINITIONS

MeshGen1D::MeshGen1D(double &a,double &b,int &c)
  : xmin(a), xmax(b), cellnumber(c) {}

void MeshGen1D::GenerateMesh(vector<double> &xcoords) {

  Tools::print("1D EULER EQ. SOLVER\n");
  xcoords = Tools::RetrievePoints(xmin,xmax,cellnumber);  
  Tools::print("- Mesh Statistics:\n");
  Tools::print("XCoords: [%e,%e]\n",xmin,xmax);
  Tools::print("Number of Cells: %d\n",cellnumber);
  
/*  for (int n=0;(int)n<xcoords.size();n++) {
    Tools::print("- Mesh Statistics:\n");
    Tools::print("XCoords: [%e,%e]\n",xmin,xmax);
    Tools::print("Number of Cells: %d\n",cellnumber);
  
  }
*/
  return;

}

MeshGen1D::~MeshGen1D(){}
