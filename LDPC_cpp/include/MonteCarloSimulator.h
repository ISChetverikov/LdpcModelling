#pragma once

#include <string>
#include <chrono>
#include "Base_decoder.h"
#include "BaseSimulator.h"

class MonteCarloSimulator : public BaseSimulator {

public:
	MonteCarloSimulator(int maxTests, int rejectionsCount, Base_decoder * decoder);
	void Run(std::vector<double> snrArray,
		std::vector<double>& ebn0Array,
		std::vector<double>& ferArray,
		std::vector<double>& sigmaArray,
		std::vector<int>& testsCountArray,
		std::vector<std::chrono::milliseconds>& elapsedTimeArray) override;
};