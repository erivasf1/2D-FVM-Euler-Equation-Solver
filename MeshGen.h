// Responsible for creating a Mesh (1D for now)
#ifndef _MESHGEN_H_
#define _MESHGEN_H_
#include "ExactNozzle.h"
#include <fstream>
#include <iostream>

class MeshGenBASE {

  public:
  vector<double> xcoords,ycoords;
  int cellnumber;
  int imax,jmax,kmax;

  MeshGenBASE();

  virtual double GetCellVolume(int cell_id);
  virtual void GenerateMesh();
  virtual void ReadMeshFile();
  virtual void OutputMesh();
  
  virtual ~MeshGenBASE();

};


class MeshGenNozzle : public MeshGenBASE { //creates a uniform mesh (in x)
  double xmin,xmax;
  double dx;
  

  public:

  MeshGenNozzle(double &a, double &b, int &c);

  double GetCellVolume(int cell_id) override;
  //static double GetCellVolume(int cell_id) override;
 
  void GenerateMesh() override;

  void OutputNozzleAreas(vector<double> &xcoords,const char *filename);

  ~MeshGenNozzle();


};

class MeshGen2D : public MeshGenBASE { //reads in a non-uniform 2D mesh
  //double xmin,xmax;
  //double ymin,ymax;
  const char* filename;

  public:
  MeshGen2D(const char* name);


  void ReadMeshFile() override;

  void OutputMesh();

  //double GetCellVolume(int cell_id) override;

  ~MeshGen2D();

};

#endif
