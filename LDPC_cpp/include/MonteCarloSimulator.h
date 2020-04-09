#pragma once

#include <string>
#include <chrono>
#include "Base_decoder.h"
#include "BaseSimulator.h"

class MonteCarloSimulator : public BaseSimulator {

public:
	MonteCarloSimulator(int maxTests, int rejectionsCount, Base_decoder * decoder);
	~MonteCarloSimulator() {};
	SimulationIterationResults Run(double snr) override;
};