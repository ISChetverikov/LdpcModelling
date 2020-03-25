#pragma once
#include <vector>

// Read usual csv matrix from filename, size of matrix - out parameters MxN
std::vector<std::vector<int>> readAsRowSparseMatrix(std::string filename);