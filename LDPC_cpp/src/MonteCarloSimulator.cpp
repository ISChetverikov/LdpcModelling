#include <random>
#include <chrono>
#include "../include/MonteCarloSimulator.h"

MonteCarloSimulator::MonteCarloSimulator(int maxTests, int maxRejectionsCount, Base_decoder * decoderPtr) {
	_decoderPtr = decoderPtr;
	_maxTestsCount = maxTests;
	_maxRejectionsCount = maxRejectionsCount;
	_n = decoderPtr->GetCodewordLegth();
}

void MonteCarloSimulator::Run(std::vector<double> snrArray,
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
		double snr = snrArray[ii];
		double sigma = sqrt(pow(10, -snr / 10) / 2);
		int tests = 0;
		int wrong_dec = 0;
		bool isFailed = false;
		std::normal_distribution<double> distribution(0, sigma);

		auto t1 = std::chrono::steady_clock::now();
		while ((tests < _maxTestsCount) && (wrong_dec < _maxRejectionsCount)) {
			tests = ++tests;

			for (size_t i = 0; i < _n; i++) {
				llrs[i] = -2 * (2 * codeword[i] - 1 + distribution(randomDevice)) / (sigma * sigma);
			}

			decoded = _decoderPtr->Decode(llrs, &isFailed);
			if (decoded != codeword)
				wrong_dec += 1;
		}
		auto t2 = std::chrono::steady_clock::now();

		elapsedTimeArray[ii] = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
		sigmaArray[ii] = sigma;
		ferArray[ii] = (double)wrong_dec / tests;
		testsCountArray[ii] = tests;
	}
}
