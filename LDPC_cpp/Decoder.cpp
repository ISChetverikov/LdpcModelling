#include "Decoder.h"
#include "Exceptions.h"
#include "MathOperations.h"

#include <vector>
#include <map>
#include <cmath>
#include <string>
#include <iostream>
#include <algorithm>

using namespace std;
// Constructor
Decoder::Decoder(vector<vector<int>> H, int iterationsCount) {
	
	if ((_iterationsCount = iterationsCount) <= 0)
		throw invalid_argument("Number of iterations is incorrect: " + to_string(_iterationsCount));

	size_t m = H.size();
	if (m <= 0)
		throw IncorrectMatrixDimensionsException("Check matrix has incorrect row size");
	size_t n = H[0].size();

	if (n <= 0)
		throw IncorrectMatrixDimensionsException("Check matrix has incorrect column size");

	for (size_t i = 0; i < m; i++)
	{
		if (n != H[0].size())
			throw IncorrectMatrixDimensionsException("Check matrix has different column size in 0 and " + to_string(i) + " row");
	}

	_checks.resize(m, vector<int>());
	_bits.resize(n, vector<int>());
	for (size_t i = 0; i < m; i++)
	{
		for (size_t j = 0; j < n; j++)
		{
			if (H[i][j] == 0)
				continue;

			_checks[i].push_back(j);
			_bits[j].push_back(i);
		}
	}
	
	_m = _checks.size();
	_n = _bits.size();

	return;
}

int sign(double x) {
	return (x < 0) ? -1 : 1;
}

double Decoder::MinSumFucntion(vector<double> vector) {
	double value = MinSumNorm * *min_element(vector.begin(), vector.end()) - MinSumOffset;
	return  value > 0 ? value : 0;
}

void Decoder::HorizontalStep(vector<map<int, int>> alpha, vector<map<int, double>> beta, vector<map<int, double>> &gamma) {
	for (size_t j = 0; j < _m; j++)
	{
		for (auto &i : _checks[j])
		{
			int sign = 1; // May be with count of sign it will be faster?
			vector<double> values;

			// TODO: Here we can get rid of redundant cycle ?
			for (auto &k : _checks[j])
			{
				if (k == i)
					continue;

				sign *= alpha[j][k];
				values.push_back(beta[j][k]);
			}

			gamma[j][i] = sign * MinSumFucntion(values);
		}
	}
}

vector<int> Decoder::Decode(vector<double> llr, bool * isFailed) {
	
	size_t n = llr.size();
	if (n != _n)
		throw IncorrectCodewordException("The codeword is not from a code with given check matrix");

	vector<int> result(n);

	vector<int> alpha0(n);
	vector<double> beta0(n);
	vector<double> bits_values(n);

	for (size_t i = 0; i < n; i++)
	{
		alpha0[i] = sign(llr[i]);
		beta0[i] = abs(llr[i]);
	}

	vector<map<int, int>> alpha(_m, map<int, int>());
	vector<map<int, double>> beta(_m, map<int, double>());
	vector<map<int, double>> gamma(_m, map<int, double>());

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

		if (_H * result == vector<int>(_m, 0)) {
			*isFailed = false;
			break;
		}

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