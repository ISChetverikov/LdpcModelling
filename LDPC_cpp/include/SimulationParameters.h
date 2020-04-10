#pragma once

#include <string>
#include <unordered_map>
#include <chrono>
#include <sstream>
#include <iostream>
#include "DecoderType.h"
#include "SimulationType.h"

struct SimulationParams {
    simulationType type;
    std::unordered_map < std::string, std::string > simulationTypeParams;
    std::vector<double> snrArray;

	decoderType decoder;
	std::unordered_map < std::string, std::string > decoderParams;
	std::string H_MatrixFilename;
	std::string G_MatrixFilename;

	std::string resultsFilename;
};

struct SimulationIterationResults {
	double snr;
	double ebn0;
	double sigma;

	double fer;

	int rejectionsCount;
	int testsCount;
	std::chrono::milliseconds elapsedTime;

	static std::string GetHeader() {
		return "SNR, EbN0, sigma, FER, rejectionsCount, testsCount, time(ms)";
	}

	std::string ToString() {
		std::stringstream ss;
		
		ss << snr << ", " << ebn0 << ", " << sigma << ", " << fer << ", " << rejectionsCount
			<< ", " << testsCount << ", " << elapsedTime.count();
		
		return ss.str();
	}
};