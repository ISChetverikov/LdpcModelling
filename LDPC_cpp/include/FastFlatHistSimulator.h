#pragma once

#include <string>
#include <chrono>
#include "Base_decoder.h"
#include "BaseSimulator.h"

class FastFlatHistSimulator : public BaseSimulator {
protected:
	int _skip_iterations = 0;
	int _minIterationsCount = 0;
	double _epsilon = 0;
	double _percent = 0;

	double normal_pdf(double x, double m, double s);
	double loss_func(const std::vector<double>& z, const std::vector<int>& codeword);
	bool hist_is_flat(std::vector<std::vector<int>>& H, size_t iter);
	std::tuple<double, double, std::vector<double>> find_opt_V(int L, double snr, const std::vector<int>& codeword,
		double sigma, double f);

public:
	FastFlatHistSimulator(int maxTests, int minTests, int rejectionsCount, Base_decoder * decoder, int skipInterations, double eps, double percent);
	~FastFlatHistSimulator() {};
	SimulationIterationResults Run(double snr) override;
};