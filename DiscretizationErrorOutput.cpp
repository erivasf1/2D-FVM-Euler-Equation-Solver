// File to Output Order of Accuracy into Tecplot (.dat) file
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <mpi.h>

#include "Output.h"

using namespace std;

int main(){

  Output Out;
  const char* filename_read = "DiscretizationErrorNorms.txt";
  const char* filename_write = "OOA_VanLeer_1stOrder_Super.txt";

  Out.CalculateOrderofAccuracy(filename_read,filename_write);


  return 0;
}
