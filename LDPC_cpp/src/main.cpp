#include <fstream>
#include <iostream>
#include <string>

#include "../include/simulate.h"
#include "../include/SimulationParameters.h"
#include "../include/FFH_simulate.h"
#include "../include/MatrixReading.h"

int main() {
	SimulationParams simulationParams;
	simulationParams.type = simulationType::MC;
	simulationParams.simulationTypeParams = std::unordered_map<std::string, std::string>(
		{
			{ "maxTestsCount", "100000" },
			{ "maxRjectionsCount", "20" }
		});
	simulationParams.snrArray = std::vector<double>();
	for (double i = -6; i <= -1; i+=0.2)
	{
		simulationParams.snrArray.push_back(i);
	}
	
	CodeParams codeParams;
	codeParams.decoderType = decoderType::ONMS;
	codeParams.decoderParams = std::unordered_map<std::string, std::string>(
		{
			{ "iterationsCount", "20" },
			{ "scale", "0.72" },
			{ "offset", "0.0" }
		});
	codeParams.H_MatrixFilename = "../Matrices/FromMatlabScript/h_3_4_512.csv";
	//codeParams.G_MatrixFilename = "../Matrices/FromMatlabScript/g_3_4_512.csv";

	/*auto results = simulate(simulationParams, codeParams);

	std::cout << "\nMC\n";
	std::cout << results.ToString() << std::endl;

	std::string resultsFilename = "../Results/FromMatlabScript/h_3_4_512.csv";
	std::ofstream resultsFile;
	resultsFile.open(resultsFilename, std::fstream::out);
	resultsFile << results.ToString();
	resultsFile.close();*/
	
	// Не убрал FFH в Simulate пока
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

	size_t m;
	size_t n;
	auto T = readSprsAsRowSparseMatrix("../Matrices/FromMatlabScript/h_3_4_512.sprs", &m, &n);

}