#include <vector>
#include <map>
#include <cmath>
#include <string>
#include <iostream>
#include <algorithm>

#include "../include/SP_decoder.h"

static double phi(double x)
{
	static const double lim = 31.;
	static const double inv_lim = 
		log( (exp((double)lim) + 1)/(exp((double)lim) - 1) );
	
	if (x > lim)
	{
		return( 0 );
	}
	else if (x < inv_lim)
		return( lim );
	else 
	{
		double t = exp( (double) x );
		return (float) log( (t+1)/(t-1) );
	}
}

// Constructor
SP_decoder::SP_decoder(std::vector<std::vector<int>> H, int iterationsCount) {
	
	//if ((_iterationsCount = iterationsCount) <= 0)
	//	throw invalid_argument("Number of iterations is incorrect: " + to_string(_iterationsCount));

	size_t m = H.size();
	//if (m <= 0)
	//	throw IncorrectMatrixDimensionsException("Check matrix has incorrect row size");
	//size_t n = H[0].size();
	_checks = H;
	int n = 0;
	for (size_t i = 0; i < m; i++)
	{
		int max = *max_element(H[i].begin(), H[i].end());
		if (max > n)
			n = max;
	}
	n++;
	//if (n <= 0)
	//	throw IncorrectMatrixDimensionsException("Check matrix has incorrect column size");

	// for (size_t i = 0; i < m; i++)
	// {
	// 	//if (n != H[0].size())
	// 		//throw IncorrectMatrixDimensionsException("Check matrix has different column size in 0 and " + to_string(i) + " row");
	// }

	//_H = H;
	//_checks.resize(m, vector<int>());
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

	return ;
}

void SP_decoder::HorizontalStep(std::vector<std::map<int, int>> alpha,
                                std::vector<std::map<int, double>> beta,
								std::vector<std::map<int, double>> &gamma) {
	for (size_t j = 0; j < _m; j++)
	{
		for (auto &i : _checks[j])
		{
			int sign = 1; // May be with count of sign it will be faster?
			double values = 0;

			// TODO: Here we can get rid of redundant cycle ?
			for (auto &k : _checks[j])
			{
				if (k == i)
					continue;

				sign *= alpha[j][k];
				values += phi(beta[j][k]);
			}

			gamma[j][i] =sign * phi(values);
		}
	}
}

std::vector<int> SP_decoder::Decode(std::vector<double> llr, bool *isFailed) {
	
	size_t n = llr.size();
	// if (n != _n)
	// 	throw IncorrectCodewordException("The codeword is not from a code with given check matrix");

	std::vector<int> result(n, 0);

	std::vector<int> alpha0(n, 0);
	std::vector<double> beta0(n, 0);
	std::vector<double> bits_values(n, 0);
	for (size_t i = 0; i < n; i++)
	{
		alpha0[i] = sign(llr[i]);
		beta0[i] = fabs(llr[i]);
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
			beta[j][i] = fabs(llr[i]);
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
			beta0[i] = fabs(bits_values[i]);
		}	

		for (size_t i = 0; i < n; i++)
		{
			result[i] = (alpha0[i] == -1) ? 1 : 0;
		}

		//if (_H * result == vector<int>(_m, 0))
		//	break;
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


		if (iteration >= _iterationsCount)
			break;

		// Vertical Step
		for (size_t i = 0; i < n; i++)
		{
			double value = bits_values[i];

			for (auto &j : _bits[i])
			{
				double new_value = value - gamma[j][i];
				alpha[j][i] = sign(new_value);
				beta[j][i] = fabs(new_value);
			}
		}
	}

	return result;
}
