#pragma once
#include <array>
#include <vector>
//For usage of i,j indexing of 2D Field variable that is stored in the memory as a single vector

inline std::array<double, 4>& fieldij(std::vector<std::array<double, 4>>* &field, int i, int j, int Nx) {
    return (*field)[i + j * Nx];
}

//Indexing example of Cell index 4 or pos. (1,1) in total cellnum of 9:
// array<double,4> cell = fieldij(field,1,1,3);
