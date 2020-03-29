#include <fstream>
#include <iostream>
#include <string>
#include <tuple>

#include "../include/simulate.h"
#include "../include/SimulationParameters.h"
#include "../include/MatrixReading.h"

void McRun() {
	// Simulation Params
	SimulationParams simulationParams;
	simulationParams.type = simulationType::MC;
	simulationParams.simulationTypeParams = std::unordered_map<std::string, std::string>(
		{
			{ "maxTestsCount", "100000" },
			{ "maxRjectionsCount", "20" }
		});
	simulationParams.snrArray = std::vector<double>();
	for (double i = -3; i <= -2.8; i += 0.2)
	{
		simulationParams.snrArray.push_back(i);
	}

	// Codes Params
	std::string matricesFolder = "../Matrices/FromMatlabScript/";
	std::string resultsFolder = "../Results/FromMatlabScript/";

	std::vector<std::string> filenameArr = {
		"h_3_4_128.sprs",
		//"h_3_4_512.sprs"
		/*"h_3_4_2048.sprs",
		"h_3_6_128.sprs",
		"h_3_6_512.sprs",
		"h_3_6_2048.sprs",
		"h_3_15_128.sprs",
		"h_3_15_512.sprs",
		"h_3_15_2048.sprs"*/
	};

	for (auto& filename : filenameArr) {

		CodeParams codeParams;
		codeParams.decoder = decoderType::ONMS;
		codeParams.decoderParams = std::unordered_map<std::string, std::string>(
			{
				{ "iterationsCount", "20" },
				{ "scale", "0.72" },
				{ "offset", "0.0" }
			});
		codeParams.H_MatrixFilename = matricesFolder + filename;

		auto results = simulate(simulationParams, codeParams);

		std::cout << "\n" + filename + "\n";
		std::cout << "============================================================\n";
		std::cout << results.ToString() << std::endl;
		std::cout << "============================================================\n";

		std::string resultsFilename = resultsFolder + filename + ".mc.results";
		std::ofstream resultsFile;
		resultsFile.open(resultsFilename, std::fstream::out);
		resultsFile << results.ToString();
		resultsFile.close();
	}
}

void FfhRun() {
	
	std::string matricesFolder = "../Matrices/FromMatlabScript/";
	std::string resultsFolder = "../Results/FromMatlabScript/";

	std::vector<std::tuple<std::string, double, double>> tests = {
		{"h_3_4_128.sprs", -3.0, -3.0}/*,
		{"h_3_4_512.sprs", -6.0, -1.0},
		{"h_3_4_2048.sprs", -6.0, -1.0},
		{"h_3_6_128.sprs", -3.0, 1.0},
		{"h_3_6_512.sprs", -3.0, 1.0},
		{"h_3_6_2048.sprs", -3.0, 1.0},
		{"h_3_15_128.sprs", 1.0, 3.6},
		{"h_3_15_512.sprs", 1.0, 3.6},
		{"h_3_15_2048.sprs", 1.0, 3.6}*/
	};

	for (auto& test : tests) {

		std::string filename;
		double minSnr;
		double maxSnr;
		std::tie(filename, minSnr, maxSnr) = test;

		SimulationParams simulationParams;
		simulationParams.type = simulationType::FFH;
		simulationParams.simulationTypeParams = std::unordered_map<std::string, std::string>(
			{
				{ "maxTestsCount", "10" },
				{ "maxRjectionsCount", "500" },
				{ "skipIterations", "2000" },
				{ "epsilon", "0.5" },
				{ "percent", "0.95" }
			});
		simulationParams.snrArray = std::vector<double>();
		for (double i = minSnr; i <= maxSnr; i += 0.2)
		{
			simulationParams.snrArray.push_back(i);
		}

		CodeParams codeParams;
		codeParams.decoder = decoderType::ONMS;
		codeParams.decoderParams = std::unordered_map<std::string, std::string>(
			{
				{ "iterationsCount", "20" },
				{ "scale", "0.72" },
				{ "offset", "0.0" }
			});
		codeParams.H_MatrixFilename = matricesFolder + filename;

		auto results = simulate(simulationParams, codeParams);

		std::cout << "\n" + filename + "\n";
		std::cout << "============================================================\n";
		std::cout << results.ToString() << std::endl;
		std::cout << "============================================================\n";

		std::string resultsFilename = resultsFolder + filename + ".ffh.results";
		std::ofstream resultsFile;
		resultsFile.open(resultsFilename, std::fstream::out);
		resultsFile << results.ToString();
		resultsFile.close();
	}
}

void LfhRun() {
	
	std::string matricesFolder = "../Matrices/FromMatlabScript/";
	std::string resultsFolder = "../Results/FromMatlabScript/";

	std::vector<std::tuple<std::string, double, double>> tests = {
		{"h_3_4_128.sprs", -3.0, -3.0}/*
		"h_3_4_512.sprs"
		"h_3_4_2048.sprs",
		"h_3_6_128.sprs",
		"h_3_6_512.sprs",
		"h_3_6_2048.sprs",
		"h_3_15_128.sprs",
		"h_3_15_512.sprs",
		"h_3_15_2048.sprs",
		"h_8_16_128.sprs",
		"h_8_16_512.sprs",
		"h_8_16_2048.sprs"*/
	};

	for (auto& test : tests) {

		std::string filename;
		double minSnr;
		double maxSnr;
		std::tie(filename, minSnr, maxSnr) = test;

		SimulationParams simulationParams;
		simulationParams.type = simulationType::LFH;
		simulationParams.simulationTypeParams = std::unordered_map<std::string, std::string>(
			{
				{ "epsilon", "0.2" },
				{ "l", "300" },
				{ "kMin", "150" },
				{ "alpha", "2" },
				{ "beta", "2" },
				{ "unconWithoutAB", "2000" },
				{ "unconWithAB", "10" },
				{ "conWithoutAB", "10" },
				{ "conWithAB", "10" }
			});
		simulationParams.snrArray = std::vector<double>();
		for (double i = minSnr; i <= maxSnr; i += 0.2)
		{
			simulationParams.snrArray.push_back(i);
		}

		CodeParams codeParams;
		codeParams.decoder = decoderType::ONMS;
		codeParams.decoderParams = std::unordered_map<std::string, std::string>(
			{
				{ "iterationsCount", "50" },
				{ "scale", "0.72" },
				{ "offset", "0.0" }
			});
		codeParams.H_MatrixFilename = matricesFolder + filename;

		auto results = simulate(simulationParams, codeParams);

		std::cout << "\n" + filename + "\n";
		std::cout << "============================================================\n";
		std::cout << results.ToString() << std::endl;
		std::cout << "============================================================\n";

		std::string resultsFilename = resultsFolder + filename + ".ffh.results";
		std::ofstream resultsFile;
		resultsFile.open(resultsFilename, std::fstream::out);
		resultsFile << results.ToString();
		resultsFile.close();
	}
}


int main() {

	FfhRun();
}