#pragma once

#include <string>
#include <chrono>
#include "Base_decoder.h"
#include "BaseSimulator.h"

class FastFlatHistSimulator : public BaseSimulator {
protected:
	int _skip_iterations = 0;
	double _epsilon = 0;
	double _percent = 0;
	
	double normal_pdf(double x, double m, double s);
	double loss_func(const std::vector<double>& z, const std::vector<int>& codeword);
	bool hist_is_flat(std::vector<std::vector<int>>& H, int iter);
	std::pair<double, double> find_opt_V(int L, int snr, const std::vector<int>& codeword,
		double sigma, double f);

public:
	FastFlatHistSimulator(int maxTests, int rejectionsCount, Base_decoder * decoder, int skipInterations, double eps, double percent);
	~FastFlatHistSimulator() {};
	void Run(std::vector<double> snrArray,
		std::vector<double>& ebn0Array,
		std::vector<double>& ferArray,
		std::vector<double>& sigmaArray,
		std::vector<int>& testsCountArray,
		std::vector<std::chrono::milliseconds>& elapsedTimeArray) override;
};