#include <vector>
#include <string>
#include <random>
#include <stdio.h>
#include <math.h>
#include <chrono>
#include <iostream>
#include <algorithm>
#include <fstream>

#include "../include/ONMS_decoder.h"
#include "../include/BF_decoder.h"
#include "../include/SP_decoder.h"
#include "../include/MatrixReading.h"

using nano_s = std::chrono::nanoseconds;
using micro_s = std::chrono::microseconds;
using milli_s = std::chrono::milliseconds;
using seconds = std::chrono::seconds;
using minutes = std::chrono::minutes;
using hours = std::chrono::hours;

void simulate(int maxTests,
	std::vector<double> snr_array,
	int rejection_count,
	std::vector<double>& fer_array,
	std::vector<double>& sigma_array,
	std::vector<int>& tests_count_array) {
	
	size_t n = 0;
	size_t m = 0;
	std::vector<std::vector<int>> H = readAsRowSparseMatrix("../Matrices/H_R1f6K76.csv", &m, &n);
	ONMS_decoder * decoder = new ONMS_decoder(H, 20, 0.72, 0);

}


int main() {
	std::vector<double> snr_array{ -4};
	std::vector<double> fer_array(snr_array.size(), 0);
	std::vector<double> sigma_array(snr_array.size(), 0);
	std::vector<int> tests_count_array(snr_array.size(), 0);

	auto t1 = std::chrono::steady_clock::now();
	simulate(10000, snr_array, 200, fer_array, sigma_array, tests_count_array);
	auto t2 = std::chrono::steady_clock::now();

	auto d_s = std::chrono::duration_cast<seconds>(t2 - t1).count();
	std::ofstream output_file;
	output_file.open("sim_res.txt", std::fstream::out);
	for (size_t i = 0; i < snr_array.size(); i++) {
		printf("\nsigma:%f,\tfer: %f,\ttests numbers: %d\n", sigma_array[i], fer_array[i], tests_count_array[i]);
		output_file << snr_array[i] << ' ' << fer_array[i] << '\n';
	}
	
	std::cout << "Seconds: " << d_s << std::endl;
	return 0;
}
