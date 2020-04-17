#include "../include/SP_decoder.h"
#include "../include/Exceptions.h"
#include "../include/MathOperations.h"

#include <vector>
#include <map>
#include <cmath>
#include <string>
#include <iostream>
#include <algorithm>

using namespace std;

static float phi(float x)
{
	static const float lim = 31.;
	static const float inv_lim =
		log((exp((double)lim) + 1) / (exp((double)lim) - 1));

	if (x > lim)
	{
		return(0);
	}
	else if (x < inv_lim)
		return(lim);
	else
	{
		double t = exp((double)x);
		return (float)log((t + 1) / (t - 1));
	}
}

// Constructor
SP_decoder::SP_decoder(std::vector<std::vector<int>> H, int iterationsCount)
	: Base_decoder(H, iterationsCount) { }

//SP_decoder::SP_decoder(vector<vector<int>> H, int iterationsCount) {
//	
	//if ((_iterationsCount = iterationsCount) <= 0)
//	//	throw invalid_argument("Number of iterations is incorrect: " + to_string(_iterationsCount));
//	_iterationsCount = iterationsCount;
//	size_t m = H.size();
	//if (m <= 0)
	// {
	// 	printf("what");
	// 	throw 1;
	// }
	//size_t n = H[0].size();
//	_checks = H;
//	size_t n = 0;
//	for (size_t i = 0; i < m; i++)
//	{
//		int max = *max_element(H[i].begin(), H[i].end());
//		if (max > n)
//			n = max;
//	}
//	n++;
///	if (n <= 0)
//		throw 1;
///
//	// for (size_t i = 0; i < m; i++)
////	// {
//	// 	if (n != H[0].size())
//	// 		throw 1;
//	// }
//
//
//	//_H = H;
//	//_checks.resize(m, vector<int>());
//	_bits.resize(n, vector<int>());
//	for (size_t j = 0; j < m; j++)
//	{
//		for (size_t i = 0; i < _checks[j].size(); i++)
//		{
//			_bits[_checks[j][i]].push_back(j);
//		}
//	}
//	
//	_m = _checks.size();
//	_n = _bits.size();
//
//	return ;
//}/
///
int sp_sign(double x) {
	return (x < 0) ? 1 : 0;
}

//void sp_decoder::HorizontalStep(vector<map<int, int>> alpha, vector<map<int, double>> beta, vector<map<int, double>> &gamma) {
//	int sign = 0;
//	for (size_t j = 0; j < _m; j++)
//	{
//		int sum_sign = 0; // May be with count of sign it will be faster?
//		double values = 0;
//		for (auto &i : _checks[j])
//		{
//			sum_sign ^= alpha[j][i];
//			values += beta[j][i];
//		}
//		for (auto &i : _checks[j])
//		{
//			sign = sum_sign ^ alpha[j][i];
//			gamma[j][i] = (1 - 2*sign) * phi(values - (beta[j][i]));
//		}
//	}
//}

vector<int> SP_decoder::Decode(vector<double> llr) {

	size_t n = llr.size();
	// if (n != _n)
	// 	throw IncorrectCodewordException("The codeword is not from a code with given check matrix");

	vector<int> result(n, 0);
	vector<double> bits_values(n, 0);
	vector<vector<double>> alpha_beta(_m);
	vector<vector<double>> gamma(_m);
	//vector<map<int, int>> alpha(_m, map<int, int>());
	//vector<map<int, double>> beta(_m, map<int, double>());
	//vector<map<int, double>> gamma(_m, map<int, double>());

	for (int j = 0; j < _m; j++)
	{
		alpha_beta[j].resize(_checks[j].size());
		gamma[j].resize(_checks[j].size());
	}
	// Init
	int counter = 0;
	for (size_t j = 0; j < _m; j++)
	{
		counter = 0;
		for (auto &i : _checks[j])
		{
			alpha_beta[j][counter++] = llr[i];
			//alpha[j][i] = sign(llr[i]);
			//beta[j][i] = phi(fabs(llr[i]));
		}
	}
	size_t iteration = 0;

	while (iteration < _iterationsCount)
	{
		iteration++;
		//HorizontalStep(alpha, beta, gamma);
		int signn = 0;
		for (size_t j = 0; j < _m; j++)
		{
			counter = 0;
			int sum_sign = 0; // May be with count of sign it will be faster?
			double values = 0;
			for (auto &i : _checks[j])
			{
				if (alpha_beta[j][counter] < 0)
					sum_sign ^= 1;
				values += phi(fabs(alpha_beta[j][counter]));
				//sum_sign ^= alpha[j][i];
				//values += beta[j][i];
				counter++;
			}
			counter = 0;
			for (auto &i : _checks[j])
			{
				signn = sum_sign;
				if (alpha_beta[j][counter] < 0)
					signn ^= 1;
				//signn = sum_sign ^ alpha[j][i];
				//gamma[j][i] = (1 - 2*sign) * phi(values - (beta[j][i]));
				gamma[j][counter] = (1 - 2 * signn) * phi(values - phi(fabs(alpha_beta[j][counter])));
				counter++;
			}
		}
		// Result of iteration
		std::vector<size_t> counter_arr1(_m);
		for (size_t i = 0; i < n; i++)
		{
			bits_values[i] = llr[i];
			for (auto &j : _bits[i])
			{
				//bits_values[i] += gamma[j][i];
				bits_values[i] += gamma[j][counter_arr1[j]++];
			}
		}

		for (size_t i = 0; i < n; i++)
			result[i] = (bits_values[i] <= 0) ? 1 : 0;

		bool isFailed = false;
		for (size_t j = 0; j < _m; j++)
		{
			int sum = 0;
			for (auto &i : _checks[j])
			{
				sum ^= result[i];
			}
			if ((bool)sum) {
				isFailed = true;
				break;
			}
		}

		if (!isFailed)
			break;

		if (iteration >= _iterationsCount)
			break;

		//Vertical Step
		std::vector<size_t> counters_arr2(_m);
		for (size_t i = 0; i < n; i++)
		{
			double value = bits_values[i];
			for (auto &j : _bits[i])
			{
				alpha_beta[j][counters_arr2[j]] = value - gamma[j][counters_arr2[j]];
				counters_arr2[j]++;
				//double new_value = value - gamma[j][i];
				//alpha[j][i] = sign(new_value);
				//beta[j][i] = phi(fabs(new_value));
			}
		}
	}
	return (result);
}
