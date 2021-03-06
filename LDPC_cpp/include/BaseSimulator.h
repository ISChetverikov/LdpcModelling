#pragma once

#include <string>
#include <chrono>
#include <math.h>
#include "Base_decoder.h"
#include "SimulationParameters.h"

class BaseSimulator {
protected:
	Base_decoder * _decoderPtr;
	size_t _n;
	size_t _m;

	double GetSigma(double snr);
	double GetEbN0(double snr, size_t m, size_t n);
public:
	BaseSimulator(Base_decoder * decoder);
	virtual ~BaseSimulator() {};
	virtual SimulationIterationResults Run(double snr) = 0;
	/*virtual void Run(std::vector<double> snrArray,
		std::vector<double>& ebn0Array,
		std::vector<double>& ferArray,
		std::vector<double>& sigmaArray,
		std::vector<int>& testsCountArray,
		std::vector<std::chrono::milliseconds>& elapsedTimeArray)=0;*/
};