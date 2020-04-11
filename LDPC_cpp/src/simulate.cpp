#include <exception>
#include <chrono>
#include <string>

#include "../include/simulate.h"
#include "../include/SimulationParameters.h"
#include "../include/Base_decoder.h"
#include "../include/ONMS_decoder.h"
#include "../include/BF_decoder.h"
#include "../include/SP_decoder.h"
#include "../include/MatrixReading.h"
#include "../include/MonteCarloSimulator.h"
#include "../include/FastFlatHistSimulator.h"
#include "../include/BaseSimulator.h"
#include "../include/ConfigReading.h"
#include "../include/Exceptions.h"

int ExtractInt(std::unordered_map<std::string, std::string> map, std::string key, std::string section) {
	if (map.count(key) <= 0)
		throw MissedParamException("Missed parameters \"" + key + "\" in section " + section);

	int result;
	try {
		result = std::stoi(map[key]);
	}
	catch (std::exception& e) {
		throw ParseParamException("Parse error of parameter \"" + key + "\": " + e.what());
	}
	
	return result;
}

double ExtractDouble(std::unordered_map<std::string, std::string> map, std::string key, std::string section) {
	if (map.count(key) <= 0)
		throw MissedParamException("Missed parameters \"" + key + "\" in section " + section);

	double result;
	try {
		result = std::stod(map[key]);
	}
	catch (std::exception& e) {
		throw ParseParamException("Parse error of parameter \"" + key + "\": " + e.what());
	}

	return result;
}


Base_decoder * BuildDecoder(
                            decoderType decoderType,
                            std::unordered_map<std::string, std::string> decoderParams,
                            std::vector<std::vector<int>> H_matrix) {
    Base_decoder * decoderPtr = NULL;
    
    switch (decoderType)
    {
        case decoderType::ONMS: {
            int interationsCount = ExtractInt(decoderParams, "iterationsCount", "ONMS decoder");
            double scale = ExtractDouble(decoderParams, "scale", "ONMS decoder");
            double offset = ExtractDouble(decoderParams, "offset", "ONMS decoder");
            
            decoderPtr = new ONMS_decoder(H_matrix, interationsCount, scale, offset);
        }
            break;
        case decoderType::MS: {
			int interationsCount = ExtractInt(decoderParams, "iterationsCount", "ONMS decoder");
            
            decoderPtr = new ONMS_decoder(H_matrix, interationsCount, 1, 0);
        }
            break;
        case decoderType::BF: {
			int interationsCount = ExtractInt(decoderParams, "iterationsCount", "BF decoder");
            
            decoderPtr = new BF_decoder(H_matrix, interationsCount);
        }
            break;
        case decoderType::SP: {
			int interationsCount = ExtractInt(decoderParams, "iterationsCount", "SP decoder");
            
            decoderPtr = new SP_decoder(H_matrix, interationsCount);
        }
            break;
        default:
            break;
    }
    
    return decoderPtr;
}

BaseSimulator * BuildSimulator(
                               simulationType simulationType,
                               std::unordered_map<std::string, std::string> simulationTypeParams,
                               Base_decoder * decoderPtr)
{
    BaseSimulator * simulator = NULL;
    
    switch (simulationType)
    {
        case simulationType::MC: {
            int maxTestsCount = ExtractInt(simulationTypeParams, "maxTestsCount", "MC simulator");
            int maxRejectionsCount = ExtractInt(simulationTypeParams, "maxRejectionsCount", "MC simulator");
            
            simulator = new MonteCarloSimulator(maxTestsCount, maxRejectionsCount, decoderPtr);
        }
            break;
        case simulationType::FFH: {
            int maxIterationssCount = ExtractInt(simulationTypeParams, "maxIterationsCount", "FFH simulator");
			int minIterationsCount = ExtractInt(simulationTypeParams, "minIterationsCount", "FFH simulator");
            int maxRejectionsCount = ExtractInt(simulationTypeParams, "maxRejectionsCount", "FFH simulator");
            int skipIterations = ExtractInt(simulationTypeParams, "skipIterations", "FFH simulator");
            double epsilon = ExtractDouble(simulationTypeParams, "epsilon", "FFH simulator");
            double percent = ExtractDouble(simulationTypeParams, "percent", "FFH simulator");
            
            simulator = new FastFlatHistSimulator(maxIterationssCount, minIterationsCount, maxRejectionsCount, decoderPtr, skipIterations, epsilon, percent);
        }
            break;
        case simulationType::LFH: {
            double epsilon = ExtractDouble(simulationTypeParams, "epsilon", "LFH simulator");
            int l = ExtractInt(simulationTypeParams, "l", "LFH simulator");
            int kMin = ExtractInt(simulationTypeParams, "kMin", "LFH simulator");
            int alpha = ExtractInt(simulationTypeParams, "alpha", "LFH simulator");
            int beta = ExtractInt(simulationTypeParams, "beta", "LFH simulator");
            int unconWithoutAB = ExtractInt(simulationTypeParams, "unconWithoutAB", "LFH simulator");
            int unconWithAB = ExtractInt(simulationTypeParams, "unconWithAB", "LFH simulator");
            int conWithoutAB = ExtractInt(simulationTypeParams, "conWithoutAB", "LFH simulator");
            int conWithAB = ExtractInt(simulationTypeParams, "conWithAB", "LFH simulator");
            int numIterForFindingV = ExtractInt(simulationTypeParams, "numIterForFindingV", "LFH simulator");
            
            simulator = new LocalFlatHistSimulator(decoderPtr, epsilon, l, kMin, alpha, beta, unconWithoutAB,
                                                   unconWithAB, conWithoutAB, conWithAB, numIterForFindingV);
        }
            break;
        default:
            break;
    }
    return simulator;
}

void LogIntoFile(std::string filename, std::string message, std::string stringPrefix="") {
	
	std::ofstream resultsFileStream;
	resultsFileStream.open(filename, std::fstream::out | std::fstream::app);

	std::string sentense;
	std::istringstream splitStream(message);
	while (std::getline(splitStream, sentense, '\n'))
	{
		if (!stringPrefix.empty())
			resultsFileStream << stringPrefix;

		resultsFileStream << sentense << std::endl;
	}
	
	resultsFileStream.close();
}

void LogIntoConsole(std::string message) {
	std::cout << message;
}

void simulate(std::string configFilename) {
    
    Base_decoder * decoderPtr = NULL;
    BaseSimulator * simulatorPtr = NULL;
	SimulationParams simulationParams;
	
	try {
		LogIntoConsole("Simulation initializing is starting....\n");

		simulationParams = ReadConfig(configFilename);
	}
	catch (const std::exception& err) {
		std::string message = "Error was ocurred:\n" + std::string(err.what()) + "\n";
		LogIntoConsole(message);
		return;
	}

	LogIntoConsole("\tConfiguration file has been read succesfully.\n");
	
	try {
        auto H_matrix = readAsRowSparseMatrix(simulationParams.H_MatrixFilename);
		LogIntoConsole("\tMatrix file has been read succesfully.\n");

        decoderPtr = BuildDecoder(simulationParams.decoder, simulationParams.decoderParams, H_matrix);
		LogIntoConsole("\tDecoder has been built succesfully.\n");

        simulatorPtr = BuildSimulator(simulationParams.type, simulationParams.simulationTypeParams, decoderPtr);
		LogIntoConsole("\tSimulator has been built succesfully.\n");

		LogIntoConsole("Simulation has been started.\n\n");
		LogIntoFile(simulationParams.resultsFilename, simulationParams.ToString(), "# ");
		LogIntoConsole(simulationParams.ToString());

		LogIntoFile(simulationParams.resultsFilename, SimulationIterationResults::GetHeader() + "\n");
		

		for (size_t i = 0; i < simulationParams.snrArray.size(); i++)
		{
			LogIntoConsole("Iteration has been started. SNR: " +  std::to_string(simulationParams.snrArray[i]) + "\n");

			auto result = simulatorPtr->Run(simulationParams.snrArray[i]);
			auto message = result.ToString() + "\n";

			LogIntoFile(simulationParams.resultsFilename, message);
			LogIntoConsole("Iteration has been ended with result:\n\t" + message);
		}
    }
    catch (const std::exception& err) {
		std::string message = "Error was ocurred:\n" + std::string(err.what()) + "\n";
		LogIntoConsole(message);
		LogIntoFile(simulationParams.resultsFilename, message);
    }
    
    delete decoderPtr;
    delete simulatorPtr;    
}
