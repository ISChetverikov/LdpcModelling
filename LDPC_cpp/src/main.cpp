#include <fstream>
#include <iostream>
#include <string>
#include <tuple>

#include "../include/simulate.h"
#include "../include/SimulationParameters.h"
#include "../include/MatrixReading.h"

void McRun() {
	
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
		simulationParams.type = simulationType::MC;
		simulationParams.simulationTypeParams = std::unordered_map<std::string, std::string>(
			{
				{ "maxTestsCount", "1000000" },
				{ "maxRjectionsCount", "2000" }
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

		std::string resultsFilename = resultsFolder + filename + ".mc.onms.results";
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
		//{"h_3_4_128.sprs", -6.0, -0.5},
		{"h_3_4_512.sprs", -6.0, -0.5},
		//{"h_3_4_2048.sprs", -6.0, -1.5},
		//{"h_3_6_128.sprs", -3.0, 1.5},
		{"h_3_6_512.sprs", -3.0, 1.5},
		//{"h_3_6_2048.sprs", -3.0, 1.5},
		//{"h_3_15_128.sprs", 1.0, 4.0},
		{"h_3_15_512.sprs", 1.0, 4.0},
		//{"h_3_15_2048.sprs", 1.0, 4.0}
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
				{ "maxIterationsCount", "10" },
				{ "minIterationsCount", "5" },
				{ "maxRjectionsCount", "500" },
				{ "skipIterations", "2000" },
				{ "epsilon", "0.1" },
				{ "percent", "0.9" }
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

		std::string resultsFilename = resultsFolder + filename + ".ffh.onms.results";
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
        {"h_3_4_128.sprs", -2.0, -1.8}/*
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
                { "epsilon", "1" },
                { "l", "200" },
                { "kMin", "100" },
                { "alpha", "2" },
                { "beta", "2" },
                { "unconWithoutAB", "2000" },
                { "unconWithAB", "4000" },
                { "conWithoutAB", "4000" },
                { "conWithAB", "5000" },
                { "numIterForFindingV", "500000" }
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

        std::string resultsFilename = resultsFolder + filename + ".lfh.results";
        std::ofstream resultsFile;
        resultsFile.open(resultsFilename, std::fstream::out);
        resultsFile << results.ToString();
        resultsFile.close();
    }
}

void test() {
	std::vector<double> v = { 265.068,262.068,264.068,267.068,267.068,271.068,276.068,289.068,283.068,288.068,290.068,290.068,295.068,294.068,302.068,304.068,307.068,310.068,311.068,311.068,311.068,308.068,311.068,314.068,311.068,313.068,312.068,312.068,312.068,310.068,309.068,308.068,311.068,306.068,310.068,303.068,305.068,307.068,307.068,301.068,298.068,300.068,295.068,284.068,296.068,279.068,289.068,288.068,280.068,275.068,285.068};

	double sum = 0;
	double ch = 0;
	for (size_t i = 0; i < v.size(); i++)
	{
		sum += std::exp(v[i]);
	}
	for (size_t i = 0; i < v.size(); i++)
	{
		ch += std::exp(v[i]) / sum;
		std::cout << std::exp(v[i])/sum << std::endl;
	}
	std::cout << sum << std::endl;
}

int main() {

	LfhRun();
}