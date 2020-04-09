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



Base_decoder * BuildDecoder(
                            decoderType decoderType,
                            std::unordered_map<std::string, std::string> decoderParams,
                            std::vector<std::vector<int>> H_matrix) {
    Base_decoder * decoderPtr = NULL;
    
    switch (decoderType)
    {
        case decoderType::ONMS: {
            int interationsCount = std::stoi(decoderParams["iterationsCount"]);
            double scale = std::stod(decoderParams["scale"]);
            double offset = std::stod(decoderParams["offset"]);
            
            decoderPtr = new ONMS_decoder(H_matrix, interationsCount, scale, offset);
            
        }
            break;
        case decoderType::MS: {
            int interationsCount = std::stoi(decoderParams["iterationsCount"]);
            
            decoderPtr = new ONMS_decoder(H_matrix, interationsCount, 1, 0);
        }
            break;
        case decoderType::BF: {
            int interationsCount = std::stoi(decoderParams["iterationsCount"]);
            
            decoderPtr = new BF_decoder(H_matrix, interationsCount);
        }
            break;
        case decoderType::SP: {
            int interationsCount = std::stoi(decoderParams["iterationsCount"]);
            
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
            int maxTestsCount = std::stoi(simulationTypeParams["maxTestsCount"]);
            int maxRjectionsCount = std::stoi(simulationTypeParams["maxRjectionsCount"]);
            
            simulator = new MonteCarloSimulator(maxTestsCount, maxRjectionsCount, decoderPtr);
        }
            break;
        case simulationType::FFH: {
            int maxIterationssCount = std::stoi(simulationTypeParams["maxIterationsCount"]);
			int minIterationsCount = std::stoi(simulationTypeParams["minIterationsCount"]);
            int maxRjectionsCount = std::stoi(simulationTypeParams["maxRjectionsCount"]);
            int skipIterations = std::stoi(simulationTypeParams["skipIterations"]);
            double epsilon = std::stod(simulationTypeParams["epsilon"]);
            double percent = std::stod(simulationTypeParams["percent"]);
            
            simulator = new FastFlatHistSimulator(maxIterationssCount, minIterationsCount, maxRjectionsCount, decoderPtr, skipIterations, epsilon, percent);
        }
            break;
        case simulationType::LFH: {
            double epsilon = std::stod(simulationTypeParams["epsilon"]);
            int l = std::stoi(simulationTypeParams["l"]);
            int kMin = std::stoi(simulationTypeParams["kMin"]);
            int alpha = std::stoi(simulationTypeParams["alpha"]);
            int beta = std::stoi(simulationTypeParams["beta"]);
            int unconWithoutAB = std::stoi(simulationTypeParams["unconWithoutAB"]);
            int unconWithAB = std::stoi(simulationTypeParams["unconWithAB"]);
            int conWithoutAB = std::stoi(simulationTypeParams["conWithoutAB"]);
            int conWithAB = std::stoi(simulationTypeParams["conWithAB"]);
            int numIterForFindingV = std::stoi(simulationTypeParams["numIterForFindingV"]);
            
            simulator = new LocalFlatHistSimulator(decoderPtr, epsilon, l, kMin, alpha, beta, unconWithoutAB,
                                                   unconWithAB, conWithoutAB, conWithAB, numIterForFindingV);
        }
            break;
        default:
            break;
    }
    return simulator;
}

void LogIntoFile(std::string filename, std::string message) {
	std::cout << "Log into file:" << message;
}

void LogIntoConsole(std::string message) {
	std::cout << message;
}

void simulate(std::string configFilename) {
    
    Base_decoder * decoderPtr = NULL;
    BaseSimulator * simulatorPtr = NULL;
	SimulationParams simulationParams;
	
    try {
		simulationParams = ReadConfig(configFilename);
        auto H_matrix = readAsRowSparseMatrix(simulationParams.H_MatrixFilename);
        decoderPtr = BuildDecoder(simulationParams.decoder, simulationParams.decoderParams, H_matrix);
        simulatorPtr = BuildSimulator(simulationParams.type, simulationParams.simulationTypeParams, decoderPtr);

		LogIntoFile(simulationParams.resultsFilename, SimulationIterationResults::GetHeader() + "\n");
		LogIntoConsole("Simulation has been started.");

		for (size_t i = 0; i < simulationParams.snrArray.size(); i++)
		{
			LogIntoConsole("Iteration has been started. SNR: " +  std::to_string(simulationParams.snrArray[i]) + "\n");

			auto result = simulatorPtr->Run(simulationParams.snrArray[i]);
			auto message = result.ToString() + "\n";

			LogIntoFile(simulationParams.resultsFilename, message);
			LogIntoConsole("Iteration has been ended with result:\n" + message);
		}
    }
    catch (const std::exception& err) {
		std::string message = "Error was ocurred:\n" + std::string(err.what());
		LogIntoConsole(message);
		LogIntoFile(simulationParams.resultsFilename, err.what());
    }
    
    delete decoderPtr;
    delete simulatorPtr;    
}
