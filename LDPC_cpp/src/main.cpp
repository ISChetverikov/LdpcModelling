#include <fstream>
#include <iostream>
#include <string>
#include <tuple>

#include "../include/simulate.h"
#include "../include/SimulationParameters.h"
#include "../include/MatrixReading.h"

void PrintUsage() {
#ifdef __linux__ 
	std::cout << "Usage: LDPC config_path" << std::endl;
#elif _WIN32
	std::cout << "Usage: LDPC.exe config_path" << std::endl;
#else

#endif
}

int main(int argc, char* argv[]) {
	
	if (argc != 2)
		PrintUsage();

	simulate(argv[1]);
}