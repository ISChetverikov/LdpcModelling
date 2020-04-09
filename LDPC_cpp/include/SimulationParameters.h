#pragma once

#include <string>
#include <unordered_map>
#include <chrono>
#include <sstream>
#include <iostream>

enum decoderType { ONMS, MS, BF, SP };
enum simulationType { MC, FFH, LFH };

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
//
//struct SimulationResult {
//    std::vector<double> snrArray;
//    std::vector<double> ebn0Array;
//	std::vector<double> sigmaArray;
//
//    std::vector<double> ferArray;
//
//	std::vector<int> rejectionsCountArray;
//    std::vector<int> testsCountArray;
//    std::vector<> elapsedTimeArray;
//    
//    bool isError = false;
//    std::string errorText = "";
//    
//    std::string ToString() {
//        std::stringstream ss;
//        
//        if (isError)
//        {
//            ss << "Error during simulation:" << std::endl;
//            ss << errorText << std::endl;
//            return ss.str();
//        }
//        
//        ss << "SNR, EbN0, FER, sigma, testsCount, time(ms)" << std::endl;
//        for (size_t i = 0; i < snrArray.size(); i++)
//        {
//            ss << snrArray[i] << ", " << ebn0Array[i] << ", " << ferArray[i] << ", " << sigmaArray[i]
//            << ", " << testsCountArray[i] << ", " << elapsedTimeArray[i].count() << std::endl;
//        }
//        return ss.str();
//    }
//};
//
