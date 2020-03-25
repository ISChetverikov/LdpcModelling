#include <fstream>
#include <iostream>
#include <string>

#include "../include/simulate.h"
#include "../include/SimulationParameters.h"

int main() {
	SimulationParams simulationParams;
	simulationParams.type = simulationType::MC;
	simulationParams.simulationTypeParams = std::unordered_map<std::string, std::string>(
		{
			{ "maxTestsCount", "10" },
			{ "maxRjectionsCount", "20" }
		});
	simulationParams.snrArray = std::vector<double>{ -5, -4, -3 };

	CodeParams codeParams;
	codeParams.decoderType = decoderType::ONMS;
	codeParams.decoderParams = std::unordered_map<std::string, std::string>(
		{
			{ "iterationsCount", "20" },
			{ "scale", "0.72" },
			{ "offset", "0.0" }
		});
	codeParams.H_MatrixFilename = "../Matrices/H_R1f6K76.csv";
	codeParams.G_MatrixFilename = "../Matrices/G_R1f6K76.csv";

	auto results = simulate(simulationParams, codeParams);

	std::cout << results.ToString() << std::endl;

	std::string resultsFilename = "../Results/H_R1f6K76.csv";
	std::ofstream resultsFile;
	resultsFile.open(resultsFilename, std::fstream::out);
	resultsFile << results.ToString();
	resultsFile.close();

}