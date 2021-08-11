#include <random>
#include <chrono>
#include <iostream>
#include <algorithm>
#include <tuple>
#include <numeric>
#include <ctime>
#include <fstream>
#include "../include/FastFlatHistSimulator.h"

FastFlatHistSimulator::FastFlatHistSimulator(Base_decoder * decoderPtr, int iterationsCount, int binCount, int maxFlatnessChecks,
	int testsForFlatnessCheck, int skipInterations, double epsilon, double percent)
	: BaseSimulator( decoderPtr) {

	_skip_iterations = skipInterations;
	_epsilon = epsilon;
	_percent = percent;
	_iterationsCount = iterationsCount;

	_L = binCount;
	_maxFlatnessChecks = maxFlatnessChecks;
	_testsForFlatnessCheck = testsForFlatnessCheck;

}

SimulationIterationResults FastFlatHistSimulator::Run(double snr)
{
	SimulationIterationResults result;

	auto t1 = std::chrono::steady_clock::now();

	double sigma = GetSigma(snr);

	std::random_device rd;
	std::mt19937 gen(rd());
	std::normal_distribution<double> norm_distr(0, sigma);
	std::uniform_real_distribution<double> un_distr(0.0, 1.0);

	std::vector<int> codeword(_n, 0);
	std::vector<double> llrs(_n, 0);
	std::vector<int> decoded(_n, 0);

	double f = 1.0;

	std::vector<std::vector<int>> H, G;
	for (size_t itn = 0; itn < _L; itn++) {
		H.push_back(std::vector<int>(_iterationsCount, 0));
		G.push_back(std::vector<int>(_iterationsCount, 0));
	}

	double Vmin, Vmax;
	std::vector<double> z(_n, 0);
	std::tie(Vmin, Vmax, z) = findStartCondition(_L, codeword, sigma, f);
	std::vector<double> start_z = z;

	std::vector<double> prob(_L, log(1.0 / _L));

#ifdef FFH_DEBUG_FILE
	std::cout << "== FFH debug file is turned on ==\n";
	
	std::vector<bool> isFlatArr(_iterationsCount, false);
	std::vector<std::vector<double>> probArr;
	for (size_t binNum = 0; binNum < _L; binNum++) {
		probArr.push_back(std::vector<double>(_iterationsCount, 0));
	}
#endif // FFH_DEBUG_FILE
	
	int ar_total = 0;
	int testsCount = 0;
	int currentBin = 0;
	int errorsCount = 0;
	int testsCountFromCheck = 0;
	size_t currentIteration = 0;

#ifdef FFH_DEBUG
	std::cout << "== FFH debug info is turned on ==\n";
#endif

	for (currentIteration = 0; currentIteration < _iterationsCount; ++currentIteration) {
		
#ifdef FFH_DEBUG
		std::cout << "--New ffh iteration---" << currentIteration << "--------------------------\n";
#endif
		z = start_z;
		bool isFlat = false;
		int flatnessCheck = 0;
		while (!isFlat && flatnessCheck < _maxFlatnessChecks) {
				
			std::vector<double> new_z = z;

			for (size_t l = 0; l < _n; ++l) {
				double dz = norm_distr(gen);
				double new_z_l_bit = z[l] + _epsilon * dz;
				double prob_bit_acceptance = std::min(normalPdf(new_z_l_bit, 0, sigma) /
					normalPdf(z[l], 0, sigma),
					1.0);
				if (un_distr(gen) <= prob_bit_acceptance)
					new_z[l] = new_z_l_bit;
			}

			double new_loss = lossFunc(new_z, codeword);
			int new_bin = (int)floor((new_loss - Vmin) / (Vmax - Vmin) * (_L - 1));
			if (new_bin < 0 || new_bin >(_L - 1)) {
				continue;
			}
			
			double prob_stat_acceptance = std::min(std::exp(prob[currentBin] - prob[new_bin]), 1.0);

			ar_total++;

			if (un_distr(gen) <= prob_stat_acceptance) {
				z = new_z;
				currentBin = new_bin;

				for (size_t i = 0; i < _n; i++) {
					llrs[i] = -2 * (2 * codeword[i] - 1 + z[i]) / (sigma * sigma);
				}
				decoded = _decoderPtr->Decode(llrs);

				testsCount++;
			}
			prob[currentBin] += f;
			H[currentBin][currentIteration] += 1;
				
			if (decoded != codeword) {
				G[currentBin][currentIteration] += 1;
				errorsCount++;
			}
			
			testsCountFromCheck += 1;
			if (testsCountFromCheck == _testsForFlatnessCheck) {
				flatnessCheck++;
				isFlat = isHistFlat(H, currentIteration);
#ifdef FFH_DEBUG_FILE
				isFlatArr[currentIteration] = isFlat;
#endif
#ifdef FFH_DEBUG
				std::cout << "--------Is flat: " << isFlat << "----------\n";
#endif
				testsCountFromCheck = 0;
			}
		}

		f = f / sqrt(2);
#ifdef FFH_DEBUG_FILE
		for (size_t binNum = 0; binNum < _L; binNum++) {
			probArr[binNum][currentIteration] = prob[binNum];
		}
#endif
	}

	if (currentIteration == _iterationsCount)
		currentIteration--;
	size_t iterationCounts = currentIteration + 1;

#ifdef FFH_DEBUG_FILE // Save info into a debug file
	std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
	std::time_t time = std::chrono::system_clock::to_time_t(now);
	auto date = std::localtime(&time);
	
	std::string resultsFilename = "ffh.prob.H.G." 
		+ std::to_string(date->tm_mon+1) + "-"
		+ std::to_string(date->tm_mday) + "_"
		+ std::to_string(date->tm_hour) + "-"
		+ std::to_string(date->tm_min) + '-'
		+ std::to_string(date->tm_sec)
		+ ".debug";

	std::ofstream resultsFile;
	resultsFile.open(resultsFilename, std::fstream::out);

	resultsFile << "# Vmin: " + std::to_string(Vmin) + ", Vmax: "
		+ std::to_string(Vmax) + ", start V: " + std::to_string(lossFunc(start_z, codeword)) << std::endl;

	resultsFile << "# Acceptance rate:" << (double)testsCount / ar_total << std::endl;

	resultsFile << "# Is flatness achieved: ";
	for (size_t i = 0; i < iterationCounts; i++)
	{
		resultsFile << isFlatArr[i] << ",";
	}
	resultsFile << std::endl;

	resultsFile << "# Histograms:" << std::endl;
	for (size_t it = 0; it < iterationCounts; it++)
	{
		resultsFile << "# Probability density:" << std::endl;
		for (size_t i = 0; i < _L; i++)
		{
			resultsFile << probArr[i][it] << ",";
		}
		resultsFile << std::endl;

		resultsFile << "# \tH" << it+1 << ":" << std::endl;
		for (size_t i = 0; i < _L; i++)
		{
			resultsFile << H[i][it] << ",";
		}
		resultsFile << std::endl;

		resultsFile << "# \tG" << it + 1 << ":" << std::endl;
		for (size_t i = 0; i < _L; i++)
		{
			resultsFile << G[i][it] << ",";
		}
		resultsFile << std::endl;
	}
	resultsFile.close();

#endif // FFH_DEBUG_FILE

	double prob_error = 0;
	double prob_sum = 0;

	double prob_mean = std::accumulate(prob.begin(), prob.end(), 0.0) / _L;
	for (auto p : prob) {
			
		prob_sum += std::exp(p - prob_mean);
	}
		
	for (size_t bin_num = 0; bin_num < _L; ++bin_num) {
		int el_num = 0, er_num = 0;
		for (size_t iter = 0; iter < iterationCounts; ++iter) {
			el_num += H[bin_num][iter];
			er_num += G[bin_num][iter];
		}
		if (el_num != 0) {
			prob_error += std::exp(prob[bin_num] - prob_mean) / prob_sum * er_num / el_num;
		}
	}
	auto t2 = std::chrono::steady_clock::now();

	result.snr = snr;
	result.ebn0 = GetEbN0(snr, _m, _n);
	result.sigma = sigma;

	result.fer = prob_error;

	result.elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
	result.testsCount = testsCount;
	result.rejectionsCount = errorsCount;	

	return result;
}

double FastFlatHistSimulator::normalPdf(double x, double m, double s) {
	static const double inv_sqrt_2pi = 0.3989422804014327;
	double a = (x - m) / s;
	return inv_sqrt_2pi / s * std::exp(-0.5 * a * a);
}


double FastFlatHistSimulator::lossFunc(const std::vector<double>& z, const std::vector<int>& codeword) {
	size_t n = codeword.size();
	double loss_sum = 0;
	for (size_t l = 0; l < n; l++) {
		int q = 1 - 2 * codeword[l]; 
		loss_sum += pow((q * z[l] > 0) * z[l], 2);
	}
	return sqrt(loss_sum / n);
}


bool FastFlatHistSimulator::isHistFlat(std::vector<std::vector<int>>& H, size_t iter) {
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


std::tuple<double, double, std::vector<double>> FastFlatHistSimulator::findStartCondition(int L, const std::vector<int>& codeword, double sigma, double f) {
	std::random_device rd;
	std::mt19937 gen(rd());
	double Vmin = 0, Vmax = 2;
	std::vector<double> prob(L, log(1.0 / L));
	std::vector<double> z(_n, 0);
	std::normal_distribution<double> norm_distr(0, sigma);
	std::uniform_real_distribution<double> un_distr(0.0, 1.0);
	int cur_bin = 0;
	std::vector<int> H(L, 0);
	auto e = L * 1000;
	for (size_t it = 0; it < e + _skip_iterations; ++it) {
		std::vector<double> new_z = z;
		for (size_t l = 0; l < _n; ++l) {
			double dz = norm_distr(gen);
			double new_z_l_bit = z[l] + _epsilon * dz;
			double prob_bit_acceptance = std::min(normalPdf(new_z_l_bit, 0, sigma) /
				normalPdf(z[l], 0, sigma),
				1.0);
			if (un_distr(gen) <= prob_bit_acceptance)
				new_z[l] = new_z_l_bit;
		}
		double new_loss = lossFunc(new_z, codeword);
		int new_bin = (int)floor((new_loss - Vmin) / (Vmax - Vmin) * (L - 1));
		new_bin = std::min(new_bin, L - 1);
		new_bin = std::max(new_bin, 0);

		double prob_stat_acceptance = std::min(std::exp(prob[cur_bin] - prob[new_bin]), 1.0);
		if (un_distr(gen) <= prob_stat_acceptance) {
			z = new_z;
			cur_bin = new_bin;
		}
		prob[cur_bin] += f;
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

	return std::make_tuple((double)min_bin * (Vmax - Vmin) / L, ((double)max_bin + 1) * (Vmax - Vmin) / L, z);
}
