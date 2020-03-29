#include <random>
#include <chrono>
#include <iostream>
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

		double sigma = GetSigma(snrArray[ii]);
		
		int iterations_count = 0;

		int L = 20;// (int)std::round(pow(10, 1 / sigma) * _n / 25);
		std::vector<std::vector<int>> H, G;
		for (size_t itn = 0; itn < L; itn++) {
			H.push_back(std::vector<int>(_maxTestsCount, 0));
			G.push_back(std::vector<int>(_maxTestsCount, 0));
		}
		double f = std::exp(1);
		int check_const = 10 * L;

		std::pair<double, double> V = find_opt_V(L, snrArray[ii], codeword, sigma, f);
		double Vmin = V.first, Vmax = V.second;

		std::vector<double> prob(L, log(1.0 / L));
		std::vector<double> z(_n, 0);
		std::normal_distribution<double> norm_distr(0, sigma);
		std::uniform_real_distribution<double> un_distr(0.0, 1.0);
		int cur_bin = 0;
		int wrong_dec = 0;

		for (size_t it = 0; it < _maxTestsCount; ++it) {
			bool is_flat = false;

			std::cout << "---------New iteration-------------------------------" << it << "---\n";

			while (!is_flat) {
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

				double prob_stat_acceptance = std::min(std::exp(prob[cur_bin]) / std::exp(prob[new_bin]), 1.0);
				if (un_distr(gen) <= prob_stat_acceptance) {
					z = new_z;
					cur_bin = new_bin;
				}
				prob[cur_bin] += log(f);
				for (size_t i = 0; i < _n; i++) {
					llrs[i] = -2 * (2 * codeword[i] - 1 + z[i]) / (sigma * sigma);
				}
				decoded = _decoderPtr->Decode(llrs, &isFailed);
				H[cur_bin][it] += 1;
				
				if (decoded != codeword) {
					G[cur_bin][it] += 1;
					wrong_dec++;
				}
				if (wrong_dec >= _maxRejectionsCount) {
					it = _maxTestsCount;
					break;
				}
				iterations_count += 1;
				if (iterations_count == check_const) {
					if (hist_is_flat(H, it)) {
						is_flat = true;
					}
					std::cout << "Is flat: " << is_flat << "\n";
					iterations_count = 0;
				}
			}
			f = sqrt(f);
		}

		double prob_error = 0;
		double prob_sum = 0;
		for (auto p : prob) {
			prob_sum += std::exp(p);
		}
		for (size_t bin_num = 0; bin_num < L; ++bin_num) {
			int el_num = 0, er_num = 0;
			for (size_t it_num = 0; it_num < _maxTestsCount; ++it_num) {
				el_num += H[bin_num][it_num];
				er_num += G[bin_num][it_num];
			}
			if (el_num != 0) {
				prob_error += std::exp(prob[bin_num]) / prob_sum * er_num / el_num;
			}
		}
		auto t2 = std::chrono::steady_clock::now();

		sigmaArray[ii] = sigma;
		ebn0Array[ii] = GetEbN0(snrArray[ii], _m, _n);
		ferArray[ii] = prob_error;
		elapsedTimeArray[ii] = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
		testsCountArray [ii] = iterations_count;
	}
}


double FastFlatHistSimulator::normal_pdf(double x, double m, double s) {
	static const double inv_sqrt_2pi = 0.3989422804014327;
	double a = (x - m) / s;
	return inv_sqrt_2pi / s * std::exp(-0.5 * a * a);
}


double FastFlatHistSimulator::loss_func(const std::vector<double>& z, const std::vector<int>& codeword) {
	int n = codeword.size();
	double loss_sum = 0;
	for (size_t l = 0; l < n; l++) {
		int q = pow(-1, codeword[l]);
		loss_sum += pow((q * z[l] < 0) * z[l], 2);
	}
	return sqrt(loss_sum / n);
}


bool FastFlatHistSimulator::hist_is_flat(std::vector<std::vector<int>>& H, int iter) {
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


std::pair<double, double> FastFlatHistSimulator::find_opt_V(int L, int snr, const std::vector<int>& codeword,
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
	for (size_t it = 0; it < GetEbN0(snr, _m, _n) * L + _skip_iterations; ++it) {
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

		double prob_stat_acceptance = std::min(std::exp(prob[cur_bin]) / std::exp(prob[new_bin]), 1.0);
		if (un_distr(gen) <= prob_stat_acceptance) {
			z = new_z;
			cur_bin = new_bin;
		}
		prob[cur_bin] += log(f);
		if (it >= _skip_iterations) {
			H[cur_bin] += 1;
		}
	}
	double min_bin = -1.0, max_bin = -1.0;
	for (size_t bin_num = 0; bin_num < L; ++bin_num) {
		if (H[bin_num] != 0) {
			max_bin = bin_num;
			if (min_bin == -1.0) {
				min_bin = bin_num;
			}
		}
	}
	return std::make_pair(min_bin / L, (max_bin + 1) / L);
}
