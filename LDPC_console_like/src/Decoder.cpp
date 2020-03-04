#include <vector>
#include <map>
#include <cmath>
#include <string>
#include <iostream>
#include <algorithm>

#include "../include/Decoder.h"

// Constructor
Decoder::Decoder(std::vector<std::vector<int>> H_row_sparse, int iterationsCount) {
	
	if ((_iterationsCount = iterationsCount) <= 0)
		throw std::invalid_argument("Number of iterations is incorrect: " + std::to_string(_iterationsCount));

	size_t m = H_row_sparse.size();
	if (m <= 0)
		throw IncorrectMatrixDimensionsException("Check matrix has incorrect row size");
	

	_checks = H_row_sparse;

	int n = 0;
	for (size_t i = 0; i < m; i++)
	{
		int max = *max_element(H_row_sparse[i].begin(), H_row_sparse[i].end());
		if (max > n)
			n = max;
	}
	n++;

	_bits.resize(n, std::vector<int>());
	for (size_t j = 0; j < m; j++)
	{
		for (size_t i = 0; i < _checks[j].size(); i++)
		{
			_bits[_checks[j][i]].push_back(j);
		}
	}
		
	_m = _checks.size();
	_n = _bits.size();

	/*for (size_t i = 0; i < n; i++)
	{
		cout << "-----------" << endl;
		for (size_t j = 0; j < _bits[i].size(); j++)
		{
			cout << _bits[i][j] << " ";
		}
		cout << endl;
	}*/

	return;
}

int sign(double x) {
	return (x < 0) ? -1 : 1;
}

double Decoder::MinSumFunction(std::vector<double> values) {
	double value = MinSumNorm * *min_element(values.begin(), values.end()) - MinSumOffset;
	return  value > 0 ? value : 0;
}

void Decoder::HorizontalStep(std::vector<std::map<int, int>> alpha,
                             std::vector<std::map<int, double>> beta,
							 std::vector<std::map<int, double>> &gamma) {
	for (size_t j = 0; j < _m; j++)
	{
		for (auto &i : _checks[j])
		{
			int sign = 1; // May be with count of sign it will be faster?
			std::vector<double> values;

			// TODO: Here we can get rid of redundant cycle ?
			for (auto &k : _checks[j])
			{
				if (k == i)
					continue;

				sign *= alpha[j][k];
				values.push_back(beta[j][k]);
			}

			gamma[j][i] = sign * MinSumFunction(values);
		}
	}
}

std::vector<int> Decoder::Decode(std::vector<double> llr, bool * isFailed) {
	
	size_t n = llr.size();
	if (n != _n)
		throw IncorrectCodewordException("The codeword is not from a code with given check matrix");

	std::vector<int> result(n);

	std::vector<int> alpha0(n);
	std::vector<double> beta0(n);
	std::vector<double> bits_values(n);

	for (size_t i = 0; i < n; i++)
	{
		alpha0[i] = sign(llr[i]);
		beta0[i] = abs(llr[i]);
	}

	std::vector<std::map<int, int>> alpha(_m, std::map<int, int>());
	std::vector<std::map<int, double>> beta(_m, std::map<int, double>());
	std::vector<std::map<int, double>> gamma(_m, std::map<int, double>());

	// Init
	for (size_t j = 0; j < _m; j++)
	{
		for (auto &i : _checks[j])
		{
			alpha[j][i] = sign(llr[i]);
			beta[j][i] = abs(llr[i]);
		}	
	}

	size_t iteration = 0;
	while (true)
	{
		iteration++;
		HorizontalStep(alpha, beta, gamma);

		// Result of iteration
		for (size_t i = 0; i < n; i++)
		{
			bits_values[i] = alpha0[i] * beta0[i];

			for (auto &j : _bits[i])
			{
				bits_values[i] += gamma[j][i];
			}

			alpha0[i] = sign(bits_values[i]);
			beta0[i] = abs(bits_values[i]);
		}
		
		for (size_t i = 0; i < n; i++)
		{
			result[i] = (alpha0[i] == -1) ? 1 : 0;
		}
/*
		if (H * result == vector<int>(_m, 0)) {
			*isFailed = false;
			break;
		}*/

		*isFailed = false;
		for (size_t j = 0; j < _m; j++)
		{
			int sum = 0;
			for (auto &i : _checks[j])
			{
				sum ^= result[i];
			}
			if ((bool)sum) {
				*isFailed = true;
				break;
			}
		}

		if (!*isFailed)
			break;

		if (iteration >= _iterationsCount) {
			*isFailed = true;
			break;
		}	

		// Vertical Step
		for (size_t i = 0; i < n; i++)
		{
			double value = bits_values[i];

			for (auto &j : _bits[i])
			{
				double new_value = value - gamma[j][i];
				alpha[j][i] = sign(new_value);
				beta[j][i] = abs(new_value);
			}
		}
	}

	return result;
}