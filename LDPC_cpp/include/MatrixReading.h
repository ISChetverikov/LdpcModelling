#pragma once
#include <vector>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>

// Read usual csv matrix from filename, size of matrix - out parameters MxN
std::vector<std::vector<int>> readAsRowSparseMatrix(std::string filename);