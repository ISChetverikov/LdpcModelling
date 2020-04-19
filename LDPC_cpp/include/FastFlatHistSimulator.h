#pragma once

#include <string>
#include <chrono>
#include "Base_decoder.h"
#include "BaseSimulator.h"

#define FFH_DEBUG_FILE

class FastFlatHistSimulator : public BaseSimulator {
protected:
	int _skip_iterations = 0;
	int _minIterationsCount = 0;
	double _epsilon = 0;
	double _percent = 0;

	double normalPdf(double x, double m, double s);
	double lossFunc(const std::vector<double>& z, const std::vector<int>& codeword);
	bool isHistFlat(std::vector<std::vector<int>>& H, size_t iter);
	std::tuple<double, double, std::vector<double>> findStartCondition(int L, const std::vector<int>& codeword, double sigma, double f);

public:
	FastFlatHistSimulator(int maxTests, int minTests, int rejectionsCount, Base_decoder * decoder, int skipInterations, double eps, double percent);
	~FastFlatHistSimulator() {};
	SimulationIterationResults Run(double snr) override;
};