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

	std::string ToString() {
		std::stringstream ss;

		ss << "Prameters of simulation run\n";
		ss << "Simulation type: " + simulationTypeToString(type) + "\n";
		ss << "Simulation parameters:\n";
		for (auto simulationParam : simulationTypeParams)
		{
			ss << "\t" << simulationParam.first + ": " + simulationParam.second + "\n";
		}

		ss << "\n";
		ss << "SNR:\n";
		for (size_t i = 0; i < snrArray.size(); i++)
		{
			ss << "\t" << snrArray[i] << "\n";
		}

		ss << "Decoder type: " + decoderTypeToString(decoder) + "\n";
		ss << "Decoder params:\n";
		for (auto decoderParam : decoderParams)
		{
			ss << "\t" << decoderParam.first + ": " + decoderParam.second + "\n";
		}

		ss << "Matrix filename:\n";
		ss << H_MatrixFilename + "\n";

		return ss.str();
	}
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