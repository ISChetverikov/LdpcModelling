#include <random>
#include <chrono>
#include "../include/SimulationParameters.h"
#include "../include/MonteCarloSimulator.h"

MonteCarloSimulator::MonteCarloSimulator(int maxTests, int maxRejectionsCount, Base_decoder * decoderPtr)
	: BaseSimulator(maxTests, maxRejectionsCount, decoderPtr) {

}

SimulationIterationResults MonteCarloSimulator::Run(double snr)
{	
	SimulationIterationResults result;

	std::random_device randomDevice;

	std::vector<int> codeword(_n, 0);
	std::vector<double> llrs(_n, 0);
	std::vector<int> decoded(_n, 0);

	auto t1 = std::chrono::steady_clock::now();

	double sigma = GetSigma(snr);
	int tests = 0;
	int wrong_dec = 0;
	std::normal_distribution<double> distribution(0, sigma);

	while (((tests < _maxTestsCount) && (wrong_dec < _maxRejectionsCount)) || ((_maxTestsCount == -1) && (wrong_dec < _maxRejectionsCount))) {
		tests++;

		for (size_t i = 0; i < _n; i++) {
			llrs[i] = -2 * (2 * codeword[i] - 1 + distribution(randomDevice)) / (sigma * sigma);
		}

		decoded = _decoderPtr->Decode(llrs);
		if (decoded != codeword)
			wrong_dec += 1;
	}

	auto t2 = std::chrono::steady_clock::now();

	result.snr = snr;
	result.ebn0 = GetEbN0(snr, _m, _n);
	result.sigma = sigma;

	result.fer = (double)wrong_dec / tests;

	result.elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
	result.testsCount = tests;
	result.rejectionsCount = wrong_dec;
	
	return result;
}
