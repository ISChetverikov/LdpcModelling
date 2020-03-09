#pragma once

#include <vector>
#include <fstream>
#include <sstream>
#include <math.h>

#include "../include/ONMS_decoder.h"

class base_simulate {
public:
    std::vector<double> _snr_array;
    std::vector<double> _fer_array;
    std::vector<double> _sigma_array;
    std::vector<int> _tests_count_array;
    std::vector<int> _tests_times;
protected:
    int _n, _m;
    ONMS_decoder* _decoder;
    std::vector<std::vector<int>> _H;
    std::vector<std::vector<int>> readAsRowSparseMatrix(std::string filename);
    double get_sigma(double snr);
    double EbN0(double snr);
};