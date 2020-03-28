#include <fstream>
#include <iostream>
#include <string>

#include "../include/simulate.h"
#include "../include/SimulationParameters.h"
#include "../include/FFH_simulate.h"
#include "../include/MatrixReading.h"

std::vector<CodeParams> MakeCodeParamsArray() {
	std::vector<CodeParams> codeParamsArray;

	

	return codeParamsArray;
}

int main() {

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

		std::string resultsFilename = resultsFolder + filename + ".res";
		std::ofstream resultsFile;
		resultsFile.open(resultsFilename, std::fstream::out);
		resultsFile << results.ToString();
		resultsFile.close();
	}
	

	// �� ����� FFH � Simulate ����
	/*auto snr_array = simulationParams.snrArray;
	auto ffh_fer_array = std::vector<double>{ 0, 0, 0, 0 };
	FFH_simulate ffh = FFH_simulate("../Matrices/FromMatlabScript/h_3_4_512.csv", "ONMS", 50, 20, snr_array, ffh_fer_array);
	ffh.simulate();
	ffh_fer_array = ffh.get_fer();

	std::string resultsFilenameFFH = "../Results/FromMatlabScript/h_3_4_512_FFH.csv";
	std::ofstream resultsFileFFH;
	resultsFileFFH.open(resultsFilenameFFH, std::fstream::out);
	std::cout << "FFH\n";
	resultsFileFFH << "SNR, FER" << std::endl;
	for (size_t i = 0; i < snr_array.size(); i++)
	{
		resultsFileFFH << snr_array[i] << ", " << ffh_fer_array[i] << std::endl;
		std::cout << snr_array[i] << ", " << ffh_fer_array[i] << std::endl;
	}
	resultsFileFFH.close();*/
	//--------------------------------------------
}