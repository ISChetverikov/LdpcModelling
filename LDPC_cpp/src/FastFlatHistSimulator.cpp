#include <random>
#include <chrono>
#include <iostream>
#include <algorithm>
#include <tuple>
#include <numeric>
#include <fstream>
#include "../include/FastFlatHistSimulator.h"

FastFlatHistSimulator::FastFlatHistSimulator(
	int maxTests, int maxRejectionsCount, Base_decoder * decoderPtr, int skipInterations, double epsilon, double percent)
	: BaseSimulator(maxTests, maxRejectionsCount, decoderPtr) {

	_skip_iterations = skipInterations;
	_epsilon = epsilon;
	_percent = percent;
}

void FastFlatHistSimulator::Run(std::vector<double> snrArray,
	std::vector<double>& ebn0Array,
	std::vector<double>& ferArray,
	std::vector<double>& sigmaArray,
	std::vector<int>& testsCountArray,
	std::vector<std::chrono::milliseconds>& elapsedTimeArray)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::vector<int> codeword(_n, 0);
	std::vector<double> llrs(_n, 0);
	std::vector<int> decoded(_n, 0);
	bool isFailed = false;

	for (size_t ii = 0; ii < snrArray.size(); ii++) {

		auto t1 = std::chrono::steady_clock::now();

		// Hard code
		int L = 100; //(int)std::round(pow(10, 1 / sigma) * _n / 50);
		const int MaxFlatnessCheck = 50;

		double sigma = GetSigma(snrArray[ii]);
	
		std::vector<std::vector<int>> H, G;
		for (size_t itn = 0; itn < L; itn++) {
			H.push_back(std::vector<int>(_maxTestsCount, 0));
			G.push_back(std::vector<int>(_maxTestsCount, 0));
		}
		double f = std::exp(1);
		int check_const = 50 * L;

		std::tuple<double, double, std::vector<double>> V = find_opt_V(L, snrArray[ii], codeword, sigma, f);
		double Vmin, Vmax;
		std::vector<double> z(_n, 0);
		std::tie(Vmin, Vmax, z) = V;
		std::cout << Vmin << "-" << Vmax << std::endl;
		std::vector<double> prob(L, log(1.0 / L));

		
		std::normal_distribution<double> norm_distr(0, sigma);
		std::uniform_real_distribution<double> un_distr(0.0, 1.0);
		int cur_bin = 0;
		int wrong_dec = 0;
		int ar_total = 0;
		int ar_ac = 0;
		bool isReached = false;
		int iterations_count = 0;
		bool is_accepted = true; // first decoding is required
		size_t current_iteration = 0;
		std::vector<bool> isFlatArr(_maxTestsCount, false);

		std::cout << "SNR: " << snrArray[ii] << " ===== sigma: " << sigma << " =========" << std::endl;
		for (current_iteration = 0; current_iteration < _maxTestsCount; ++current_iteration) {
			bool is_flat = false;

			std::cout << "--New iteration---" << current_iteration << "--------------------------\n";

			int flatnessCheck = 0;
			while (!is_flat && flatnessCheck < MaxFlatnessCheck) {
				
				std::vector<double> new_z = z;

				for (size_t l = 0; l < _n; ++l) {
					double dz = norm_distr(gen);
					double new_z_l_bit = z[l] + _epsilon * dz;
					double prob_bit_acceptance = std::min(normal_pdf(new_z_l_bit, 0, sigma) /
						normal_pdf(z[l], 0, sigma),
						1.0);
					if (un_distr(gen) <= prob_bit_acceptance)
						new_z[l] = new_z_l_bit;
				}

				double new_loss = loss_func(new_z, codeword);
				int new_bin = (int)floor((new_loss - Vmin) / (Vmax - Vmin) * (L - 1));
				if (new_bin < 0 || new_bin >(L - 1))
					continue;

				//new_bin = std::min(new_bin, L - 1);
				//new_bin = std::max(new_bin, 0);

				double prob_stat_acceptance = std::min(std::exp(prob[cur_bin] - prob[new_bin]), 1.0);
				ar_total++;
				if (un_distr(gen) <= prob_stat_acceptance) {
					z = new_z;
					cur_bin = new_bin;
					ar_ac++;
					is_accepted = true;
				}
				prob[cur_bin] += log(f);
				if (is_accepted) { // economize on not accepted z - it does not change
					for (size_t i = 0; i < _n; i++) {
						llrs[i] = -2 * (2 * codeword[i] - 1 + z[i]) / (sigma * sigma);
					}
					decoded = _decoderPtr->Decode(llrs, &isFailed);
					is_accepted = false;
				}
				
				H[cur_bin][current_iteration] += 1;
				
				if (decoded != codeword) {
					G[cur_bin][current_iteration] += 1;
					if (current_iteration != 0)      // skip one iteration
						wrong_dec++;
				}
				//if (!isReached && wrong_dec >= _maxRejectionsCount) {
					//isReached = true;
					//std::cout << "----Found rejections: " << wrong_dec << " --------------\n";
				//}
				iterations_count += 1;
				if (iterations_count == check_const) {
					flatnessCheck++;
					if (hist_is_flat(H, current_iteration)) {
						is_flat = true;
						isFlatArr[current_iteration] = true;
					}
					std::cout << "------Is flat: " << is_flat << "----------\n";
					iterations_count = 0;
				}
			}
			f = sqrt(f);
			if (isReached)
				break;
		}

		if (current_iteration == _maxTestsCount)
			current_iteration--;
		size_t iteration_counts = current_iteration + 1;
//#define DEBUG
#ifdef DEBUG // Save info to file
		
		std::string resultsFilename = "ffh.prob.H.G.debug";
		std::ofstream resultsFile;
		resultsFile.open(resultsFilename, std::fstream::out);

		for (size_t i = 0; i < prob.size(); i++)
		{
			resultsFile << prob[i] << ",";
		}
		resultsFile << std::endl;
		std::cout << (double)ar_ac/ar_total << std::endl;
		for (size_t i = 0; i < iteration_counts; i++)
		{
			std::cout << isFlatArr[i] << ",";
			resultsFile << isFlatArr[i] << ",";
		}
		resultsFile << std::endl;
		for (size_t it = 0; it < iteration_counts; it++)
		{
			for (size_t i = 0; i < prob.size(); i++)
			{
				resultsFile << H[i][it] << ",";
			}
			resultsFile << std::endl;
			for (size_t i = 0; i < prob.size(); i++)
			{
				resultsFile << G[i][it] << ",";
			}
			resultsFile << std::endl;
		}
		resultsFile.close();

#endif // DEBUG

		double prob_error = 0;
		double prob_error2 = 0;
		double prob_sum = 0;

		double prob_mean = std::accumulate(prob.begin(), prob.end(), 0.0) / L;
		for (auto p : prob) {
			
			prob_sum += std::exp(p - prob_mean);
		}
		
		int iterationsCountsTotal = 0;
		for (size_t bin_num = 0; bin_num < L; ++bin_num) {
			int el_num = 0, er_num = 0;
			int el_num2 = 0, er_num2 = 0;
			for (size_t iter = 0; iter < iteration_counts; ++iter) {
				el_num += H[bin_num][iter];
				er_num += G[bin_num][iter];
			}
			if (el_num != 0) {
				iterationsCountsTotal += el_num;
				prob_error += std::exp(prob[bin_num] - prob_mean) / prob_sum * er_num / el_num;
			}
		}
		auto t2 = std::chrono::steady_clock::now();

		sigmaArray[ii] = sigma;
		ebn0Array[ii] = GetEbN0(snrArray[ii], _m, _n);
		ferArray[ii] = prob_error;
		elapsedTimeArray[ii] = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
		testsCountArray [ii] = iterationsCountsTotal;
	}
}

double FastFlatHistSimulator::normal_pdf(double x, double m, double s) {
	static const double inv_sqrt_2pi = 0.3989422804014327;
	double a = (x - m) / s;
	return inv_sqrt_2pi / s * std::exp(-0.5 * a * a);
}


double FastFlatHistSimulator::loss_func(const std::vector<double>& z, const std::vector<int>& codeword) {
	size_t n = codeword.size();
	double loss_sum = 0;
	for (size_t l = 0; l < n; l++) {
		int q = 1 - 2 * codeword[l]; 
		loss_sum += pow((q * z[l] > 0) * z[l], 2);
	}
	return sqrt(loss_sum / n);
}


bool FastFlatHistSimulator::hist_is_flat(std::vector<std::vector<int>>& H, size_t iter) {
	int sum = 0;
	for (size_t bin_num = 0; bin_num < H.size(); ++bin_num) {
		sum += H[bin_num][iter];
	}
	double floor_num = (double)sum / H.size() * _percent;
	for (size_t bin_num = 0; bin_num < H.size(); ++bin_num) {
		if (H[bin_num][iter] < floor_num)
			return false;
	}
	return true;
}


std::tuple<double, double, std::vector<double>> FastFlatHistSimulator::find_opt_V(int L, double snr, const std::vector<int>& codeword,
	double sigma, double f) {
	std::random_device rd;
	std::mt19937 gen(rd());
	double Vmin = 0, Vmax = 1;
	std::vector<double> prob(L, log(1.0 / L));
	std::vector<double> z(_n, 0);
	std::normal_distribution<double> norm_distr(0, sigma);
	std::uniform_real_distribution<double> un_distr(0.0, 1.0);
	int cur_bin = 0;
	std::vector<int> H(L, 0);
	auto e = GetEbN0(snr, _m, _n) * L * 200;
	for (size_t it = 0; it < e + _skip_iterations; ++it) {
		std::vector<double> new_z = z;
		for (size_t l = 0; l < _n; ++l) {
			double dz = norm_distr(gen);
			double new_z_l_bit = z[l] + _epsilon * dz;
			double prob_bit_acceptance = std::min(normal_pdf(new_z_l_bit, 0, sigma) /
				normal_pdf(z[l], 0, sigma),
				1.0);
			if (un_distr(gen) <= prob_bit_acceptance)
				new_z[l] = new_z_l_bit;
		}
		double new_loss = loss_func(new_z, codeword);
		int new_bin = (int)floor((new_loss - Vmin) / (Vmax - Vmin) * (L - 1));
		new_bin = std::min(new_bin, L - 1);
		new_bin = std::max(new_bin, 0);

		double prob_stat_acceptance = std::min(std::exp(prob[cur_bin] - prob[new_bin]), 1.0);
		if (un_distr(gen) <= prob_stat_acceptance) {
			z = new_z;
			cur_bin = new_bin;
		}
		prob[cur_bin] += log(f);
		if (it >= _skip_iterations) {
			H[cur_bin] += 1;
		}
	}
	size_t min_bin = L, max_bin = L;
	for (size_t bin_num = 0; bin_num < L; ++bin_num) {
		if (H[bin_num] > 0) {
			max_bin = bin_num;
			if (min_bin == L) {
				min_bin = bin_num;
			}
		}
	}

	return std::make_tuple((double)min_bin / L, ((double)max_bin + 1) / L, z);
}
