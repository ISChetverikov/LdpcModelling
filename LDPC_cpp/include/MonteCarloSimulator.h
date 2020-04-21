#pragma once

#include <string>
#include <chrono>
#include "Base_decoder.h"
#include "BaseSimulator.h"

class MonteCarloSimulator : public BaseSimulator {
protected:
	int _maxTestsCount = 0;
	int _maxRejectionsCount = 0;
public:
	MonteCarloSimulator(int maxTests, int rejectionsCount, Base_decoder * decoder);
	~MonteCarloSimulator() {};
	SimulationIterationResults Run(double snr) override;
};