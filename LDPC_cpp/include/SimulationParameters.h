#pragma once

#include <string>
#include <unordered_map>
#include <chrono>
#include <sstream>
#include <iostream>

enum decoderType { ONMS, MS, BF };
enum simulationType { MC, FFH, LFH };

struct CodeParams{
    decoderType decoder;
    std::unordered_map < std::string, std::string > decoderParams;
    std::string H_MatrixFilename;
    std::string G_MatrixFilename;
};

struct SimulationParams {
    simulationType type;
    std::unordered_map < std::string, std::string > simulationTypeParams;
    std::vector<double> snrArray;
};

struct SimulationResult {
    std::vector<double> snrArray;
    std::vector<double> ebn0Array;
    std::vector<double> ferArray;
    std::vector<double> sigmaArray;
    std::vector<int> testsCountArray;
    std::vector<std::chrono::milliseconds> elapsedTimeArray;
    
    bool isError = false;
    std::string errorText = "";
    
    std::string ToString() {
        std::stringstream ss;
        
        if (isError)
        {
            ss << "Error during simulation:" << std::endl;
            ss << errorText << std::endl;
            return ss.str();
        }
        
        ss << "SNR, EbN0, FER, sigma, testsCount, time(ms)" << std::endl;
        for (size_t i = 0; i < snrArray.size(); i++)
        {
            ss << snrArray[i] << ", " << ebn0Array[i] << ", " << ferArray[i] << ", " << sigmaArray[i]
            << ", " << testsCountArray[i] << ", " << elapsedTimeArray[i].count() << std::endl;
        }
        return ss.str();
    }
};

