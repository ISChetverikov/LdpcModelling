#include <random>
#include <chrono>
#include <iostream>
#include "../include/MonteCarloSimulator.h"

MonteCarloSimulator::MonteCarloSimulator(int maxTests, int maxRejectionsCount, Base_decoder * decoderPtr)
	: BaseSimulator(maxTests, maxRejectionsCount, decoderPtr) {

}

void MonteCarloSimulator::Run(std::vector<double> snrArray,
	std::vector<double>& ebn0Array,
	std::vector<double>& ferArray,
	std::vector<double>& sigmaArray,
	std::vector<int>& testsCountArray,
	std::vector<std::chrono::milliseconds>& elapsedTimeArray)
{	
	std::random_device randomDevice;

	std::vector<int> codeword(_n, 0);
	std::vector<double> llrs(_n, 0);
	std::vector<int> decoded(_n, 0);

	for (size_t ii = 0; ii < snrArray.size(); ii++) {
		auto t1 = std::chrono::steady_clock::now();

		double snr = snrArray[ii];
		double sigma = GetSigma(snr);
		int tests = 0;
		int wrong_dec = 0;
		bool isFailed = false;
		std::normal_distribution<double> distribution(0, sigma);

		while (((tests < _maxTestsCount) && (wrong_dec < _maxRejectionsCount)) || ((_maxTestsCount == -1) && (wrong_dec < _maxRejectionsCount))) {
			tests++;

			for (size_t i = 0; i < _n; i++) {
				llrs[i] = -2 * (2 * codeword[i] - 1 + distribution(randomDevice)) / (sigma * sigma);
			}

			decoded = _decoderPtr->Decode(llrs, &isFailed);
            if (decoded != codeword) {
				wrong_dec += 1;
            }
		}
		auto t2 = std::chrono::steady_clock::now();

		elapsedTimeArray[ii] = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
		sigmaArray[ii] = sigma;
		ebn0Array[ii] = GetEbN0(snr, _m, _n);
		ferArray[ii] = (double)wrong_dec / tests;
		testsCountArray[ii] = tests;
	}
}
