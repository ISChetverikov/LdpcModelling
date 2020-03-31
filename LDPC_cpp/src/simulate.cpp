#include <exception>
#include <chrono>

#include "../include/simulate.h"
#include "../include/SimulationParameters.h"
#include "../include/Base_decoder.h"
#include "../include/ONMS_decoder.h"
#include "../include/BF_decoder.h"
#include "../include/MatrixReading.h"
#include "../include/MonteCarloSimulator.h"
#include "../include/FastFlatHistSimulator.h"
#include "../include/BaseSimulator.h"


Base_decoder * BuildDecoder(
                            decoderType decoderType,
                            std::unordered_map<std::string, std::string> decoderParams,
                            std::vector<std::vector<int>> H_matrix)
{
    Base_decoder * decoderPtr;
    
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
    BaseSimulator * simulator;
    
    switch (simulationType)
    {
        case simulationType::MC: {
            int maxTestsCount = std::stoi(simulationTypeParams["maxTestsCount"]);
            int maxRjectionsCount = std::stoi(simulationTypeParams["maxRjectionsCount"]);
            
            simulator = new MonteCarloSimulator(maxTestsCount, maxRjectionsCount, decoderPtr);
        }
            break;
        case simulationType::FFH: {
            int maxTestsCount = std::stoi(simulationTypeParams["maxTestsCount"]);
            int maxRjectionsCount = std::stoi(simulationTypeParams["maxRjectionsCount"]);
            int skipIterations = std::stoi(simulationTypeParams["skipIterations"]);
            double epsilon = std::stod(simulationTypeParams["epsilon"]);
            double percent = std::stod(simulationTypeParams["percent"]);
            
            simulator = new FastFlatHistSimulator(maxTestsCount, maxRjectionsCount, decoderPtr, skipIterations, epsilon, percent);
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

SimulationResult simulate(SimulationParams simulationParams, CodeParams codeParams) {
    
    Base_decoder * decoderPtr = NULL;
    BaseSimulator * simulatorPtr = NULL;
    
    SimulationResult results;
    results.snrArray = simulationParams.snrArray;
    
    auto snrCount = results.snrArray.size();
    results.ebn0Array = std::vector<double>(snrCount, 0);
    results.ferArray = std::vector<double>(snrCount, 0);
    results.sigmaArray = std::vector<double>(snrCount, 0);
    results.testsCountArray = std::vector<int>(snrCount, 0);
    results.elapsedTimeArray = std::vector<std::chrono::milliseconds>(snrCount);
    
    try {
        auto H_matrix = readAsRowSparseMatrix(codeParams.H_MatrixFilename);
        decoderPtr = BuildDecoder(codeParams.decoder, codeParams.decoderParams, H_matrix);
        simulatorPtr = BuildSimulator(simulationParams.type, simulationParams.simulationTypeParams, decoderPtr);
        simulatorPtr->Run(simulationParams.snrArray,
                          results.ebn0Array, results.ferArray, results.sigmaArray, results.testsCountArray, results.elapsedTimeArray);
    }
    catch (const std::exception& err) {
        results.isError = true;
        results.errorText = err.what();
    }
    
    delete decoderPtr;
    delete simulatorPtr;
    
    return results;
    
}
