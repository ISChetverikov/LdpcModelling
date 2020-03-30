#include "../include/BaseSimulator.h"

BaseSimulator::BaseSimulator(int maxTests, int maxRejectionsCount, Base_decoder * decoderPtr) {
	_decoderPtr = decoderPtr;
	_maxTestsCount = maxTests;
	_maxRejectionsCount = maxRejectionsCount;
	_n = decoderPtr->GetCodewordLegth();
	_m = decoderPtr->GetChecksSymbolsCount();
}

double BaseSimulator::GetSigma(double snr) {
	return sqrt(pow(10, -snr / 10) / 2);
}

double BaseSimulator::GetEbN0(double snr, size_t m, size_t n) {
	return snr - 10 * log10(1 - (double)m / n);
}