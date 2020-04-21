#pragma once

#include <string>
#include <chrono>
#include "Base_decoder.h"
#include "BaseSimulator.h"

#define FFH_DEBUG_FILE
#define FFH_DEBUG

class FastFlatHistSimulator : public BaseSimulator {
protected:
	int _skip_iterations = 0;
	int _iterationsCount = 0;
	double _epsilon = 0;
	double _percent = 0;

	double normalPdf(double x, double m, double s);
	double lossFunc(const std::vector<double>& z, const std::vector<int>& codeword);
	bool isHistFlat(std::vector<std::vector<int>>& H, size_t iter);
	std::tuple<double, double, std::vector<double>> findStartCondition(int L, const std::vector<int>& codeword, double sigma, double f);

public:
	FastFlatHistSimulator(Base_decoder * decoder, int iterationsCount, int skipInterations, double eps, double percent);
	~FastFlatHistSimulator() {};
	SimulationIterationResults Run(double snr) override;
};