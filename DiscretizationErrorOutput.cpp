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
  const char* filename_read = "DiscretizationError.txt";
  const char* filename_write = "ObservedOrderofAccuracy.dat";

  Out.CalculateOrderofAccuracy(filename_read,filename_write);


  return 0;
}
