#include <vector>
#include <string>
#include <random>
#include <stdio.h>
#include <math.h>
#include <chrono>
#include <iostream>
#include <algorithm>
#include <execution>
#include "Decoder.h"

using nano_s = std::chrono::nanoseconds;
using micro_s = std::chrono::microseconds;
using milli_s = std::chrono::milliseconds;
using seconds = std::chrono::seconds;
using minutes = std::chrono::minutes;
using hours = std::chrono::hours;

std::vector<std::vector<int>> read_ldpc_mat() {
	//TODO: read H from ldpc_filename
    int K = 12, M = 6;
    std::vector<std::vector<int>> H(M, std::vector<int>(K, 0));
    H[0][0] = H[0][1] = H[0][2] = H[0][5] = H[0][6] = H[0][10] = 1;
    H[1][0] = H[1][1] = H[1][2] = H[1][3] = H[1][4] = H[1][11] = 1;
    H[2][5] = H[2][6] = H[2][7] = H[2][9] = H[2][10] = H[2][11] = 1;
    H[3][0] = H[3][3] = H[3][7] = H[3][8] = H[3][9] = H[3][11] = 1;
    H[4][1] = H[4][3] = H[4][4] = H[4][6] = H[4][7] = H[4][8] = 1;
    H[5][2] = H[5][4] = H[5][5] = H[5][8] = H[5][9] = H[5][10] = 1;
    return H;
}

void simulate(int maxTests,
	std::vector<double> snr_array,
	double rejection_count,
	std::vector<double>& fer_array,
	std::vector<double>& sigma_array,
	std::vector<int>& tests_count_array) {
	
    std::vector<std::vector<int>> H = read_ldpc_mat();
	Decoder decoder = Decoder(H, 20);
	int m = H.size(), n = H[0].size();
	int k = n - m;

	std::random_device randomDevice;

	std::vector<int> codeword(n, 0);
	std::vector<double> llrs(n, 0);

    for (size_t ii = 0; ii < snr_array.size(); ii++) {
    	double snr = snr_array[ii];
    	double sigma = sqrt(pow(10, -snr / 10) / 2); 
    	int tests = 0;
        int wrong_dec = 0;
        int errors = 0;
		bool isFailed = false;
		std::normal_distribution<double> distribution(0, sigma);

        while ((tests < maxTests) && (wrong_dec < rejection_count)) {
            tests = ++tests;
            
            for (int i = 0; i < n; i++) {
				llrs[i] = -2 * (2 * codeword[i] - 1 + distribution(randomDevice)) / (sigma * sigma);
            }
			            
			decoder.Decode(llrs, &isFailed);
			
            wrong_dec += isFailed;
        }

		sigma_array[ii] = sigma;
		fer_array[ii] = (double)wrong_dec / tests;
		tests_count_array[ii] = tests;	
    }
}


int main() {
	std::vector<double> snr_array{ 3 };
	std::vector<double> fer_array(snr_array.size(), 0);
	std::vector<double> sigma_array(snr_array.size(), 0);
	std::vector<int> tests_count_array(snr_array.size(), 0);

	auto t1 = std::chrono::steady_clock::now();
	simulate(10000, snr_array, 30, fer_array, sigma_array, tests_count_array);
	auto t2 = std::chrono::steady_clock::now();

	auto d_s = std::chrono::duration_cast<seconds>(t2 - t1).count();
	for (size_t i = 0; i < snr_array.size(); i++) {
		printf("sigma:%f,\tfer: %f,\ttests numbers: %d\n", sigma_array[i], fer_array[i], tests_count_array[i]);
	}
	
	std::cout << d_s << std::endl;
	return 0;
}