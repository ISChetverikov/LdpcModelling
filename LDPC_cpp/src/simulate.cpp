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
#include "../include/FFH_simulate.h"

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
	std::vector<double>& fer_array2,
	std::vector<double>& sigma_array,
	std::vector<int>& tests_count_array) {
	
	size_t n = 0;
	size_t m = 0;
    std::vector<std::vector<int>> H = readAsRowSparseMatrix("../Matrices/H_R1f6K76.csv", &m, &n);
	ONMS_decoder * decoder = new ONMS_decoder(H, 20);

	std::random_device randomDevice;

	std::vector<int> codeword(n, 0);
	std::vector<double> llrs(n, 0);
	std::vector<int> decoded(n, 0);

    for (size_t ii = 0; ii < snr_array.size(); ii++) {
    	double snr = snr_array[ii];
    	double sigma = sqrt(pow(10, -snr / 10) / 2); 
    	int tests = 0;
        int my_wrong_dec = 0;
		int fedor_wrong_dec = 0;
        int errors = 0;
		bool isFailed = false;
		std::normal_distribution<double> distribution(0, sigma);

        while ((tests < maxTests) && (my_wrong_dec < rejection_count)) {
            tests = ++tests;
            
            for (size_t i = 0; i < n; i++) {
				llrs[i] = -2 * (2 * codeword[i] - 1 + distribution(randomDevice)) / (sigma * sigma);
            }
			   
			decoded = decoder->Decode(llrs, &isFailed);
			if (decoded != codeword)
				my_wrong_dec += 1;
        }

		sigma_array[ii] = sigma;
		fer_array[ii] = (double)my_wrong_dec / tests;
		fer_array2[ii] = (double)fedor_wrong_dec / tests; //
		tests_count_array[ii] = tests;	
    }
}


int main() {
	std::vector<double> snr_array{ -4};
	std::vector<double> fer_array(snr_array.size(), 0);
	std::vector<double> ffh_fer_array(snr_array.size(), 0);
	std::vector<double> fer_array2(snr_array.size(), 0);
	std::vector<double> sigma_array(snr_array.size(), 0);
	std::vector<int> tests_count_array(snr_array.size(), 0);

	auto t1 = std::chrono::steady_clock::now();
	simulate(100000, snr_array, 1000, fer_array, fer_array2, sigma_array, tests_count_array);
	auto t2 = std::chrono::steady_clock::now();

	auto d_s = std::chrono::duration_cast<seconds>(t2 - t1).count();
	std::ofstream output_file;
	output_file.open("sim_res.txt", std::fstream::out);
	for (size_t i = 0; i < snr_array.size(); i++) {
		printf("\nsigma:%f,\tfer: %f,\ttests numbers: %d\n", sigma_array[i], fer_array[i], tests_count_array[i]);
		output_file << snr_array[i] << ' ' << fer_array[i] << '\n';
	}

    std::cout << "Seconds: " << d_s << std::endl;

	t1 = std::chrono::steady_clock::now();
	FFH_simulate ffh = FFH_simulate("../Matrices/H_R1f6K76.csv", "ONMS", 6, 1000, snr_array, ffh_fer_array);
	ffh.simulate();
	ffh_fer_array = ffh.get_fer();
	t2 = std::chrono::steady_clock::now();

	d_s = std::chrono::duration_cast<seconds>(t2 - t1).count();
	std::ofstream output_file2;
	output_file2.open("sim_res.txt", std::fstream::out);
	for (size_t i = 0; i < snr_array.size(); i++) {
		printf("\nsigma:%f,\tfer: %f\n", sigma_array[i], ffh_fer_array[i]);
		output_file2 << snr_array[i] << ' ' << fer_array[i] << '\n';
	}
	
	std::cout << "Seconds: " << d_s << std::endl;
	return 0;
}
