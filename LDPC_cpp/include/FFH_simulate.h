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

#include "../include/ONMS_decoder.h"

class FFH_simulate {
private:
    std::vector<double> _snr_array;
    std::vector<double> _fer_array;
    std::vector<double> _sigma_array;
    std::vector<int> _tests_count_array;
    std::vector<int> _tests_times;
    size_t _n, _m;
    ONMS_decoder* _decoder;
    std::vector<std::vector<int>> _H;
    int _num_iterations, _skip_iterations, _rejections_num;
    double _eps, _perc;
    double normal_pdf(double x, double m, double s);
    double loss_func(const std::vector<double>& z, const std::vector<int>& codeword);
    bool hist_is_flat(std::vector<std::vector<int>>& H, int iter);
    std::pair<double, double> find_opt_V(int L, int snr, const std::vector<int>& codeword,
                                         double sigma, double f);
    double get_sigma(double snr);
    double EbN0(double snr, double R);
public:
    FFH_simulate(std::string matrix_name,
                 std::string decoder,
                 int num_iterations,
                 int rejections_num,
	             const std::vector<double>& snr_array,
                 std::vector<double>& fer_array);
    void simulate();
    std::vector<double> get_fer();
};