#pragma once

#include <string>
#include <chrono>
#include "Base_decoder.h"

class MonteCarloSimulator {
private:
	int _maxTestsCount = 10000;
	int _maxRejectionsCount = 20;
	Base_decoder * _decoderPtr;
	size_t _n;
public:
	MonteCarloSimulator(int maxTests, int rejectionsCount, Base_decoder * decoder);
	void Run(std::vector<double> snrArray,
		std::vector<double>& ferArray,
		std::vector<double>& sigmaArray,
		std::vector<int>& testsCountArray,
		std::vector<std::chrono::milliseconds>& elapsedTimeArray);
};