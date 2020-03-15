#pragma once

#include <vector>
#include <string>
#include <random>
#include <stdio.h>
#include <math.h>
#include <chrono>
#include <iostream>
#include <algorithm>
#include <execution>
#include <fstream>
#include <sstream>

#include "../include/base_simulate.h"
#include "../include/ONMS_decoder.h"

class FFH_simulate : public base_simulate {
private:
    int _num_iterations;
    double _eps, _perc;
    double normal_pdf(double x, double m, double s);
    double loss_func(const std::vector<double>& z, const std::vector<int>& codeword);
    bool hist_is_flat(std::vector<std::vector<int>>& H, int iter);
    std::pair<double, double> find_opt_V(int L, int snr, const std::vector<int>& codeword,
                                         double sigma, double f);
public:
    FFH_simulate(std::string matrix_name,
                 std::string decoder,
                 int num_iterations,
	             const std::vector<double>& snr_array);
    void simulate();
};