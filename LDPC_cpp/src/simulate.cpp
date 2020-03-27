#include <exception>
#include <chrono>

#include "../include/simulate.h"
#include "../include/SimulationParameters.h"
#include "../include/Base_decoder.h"
#include "../include/ONMS_decoder.h"
#include "../include/MatrixReading.h"
#include "../include/MonteCarloSimulator.h"
#include "../include/BaseSimulator.h"


Base_decoder * BuildDecoder(
	decoderType decoderType,
	std::unordered_map<std::string, std::string> decoderParams,
	std::vector<std::vector<int>> H_matrix)
{
	Base_decoder * decoderPtr;

	switch (decoderType)
	{
	case decoderType::ONMS:
	{
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
		auto H_matrix = readCsvAsRowSparseMatrix(codeParams.H_MatrixFilename);
		decoderPtr = BuildDecoder(codeParams.decoderType, codeParams.decoderParams, H_matrix);
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
