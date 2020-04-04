#include <vector>
#include <map>
#include <cmath>
#include <string>
#include <iostream>
#include <algorithm>

#include "../include/SP_decoder.h"

#include <vector>
#include <map>
#include <cmath>
#include <string>
#include <iostream>
#include <algorithm>

int spSign(double x) {
	return (x < 0) ? 1 : 0;
}

static float phi(float x)
{
	static const float lim = 31.;
	static const float inv_lim = 
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
SP_decoder::SP_decoder(std::vector<std::vector<int>> H, int iterationsCount) 
                 : Base_decoder(H, iterationsCount) { }


void SP_decoder::HorizontalStep(std::vector<std::map<int, int>> alpha, std::vector<std::map<int, double>> beta, std::vector<std::map<int, double>> &gamma) {
	int sign = 0;
	for (size_t j = 0; j < _m; j++)
	{
		int sum_sign = 0; // May be with count of sign it will be faster?
		double values = 0;
		for (auto &i : _checks[j])
		{
			sum_sign ^= alpha[j][i];
			values += beta[j][i];
		}
		for (auto &i : _checks[j])
		{
			sign = sum_sign ^ alpha[j][i];
			gamma[j][i] = (1 - 2*sign) * phi(values - (beta[j][i]));
		}
	}
}

std::vector<int> SP_decoder::Decode(std::vector<double> llr, bool *isFailed) {
	size_t n = llr.size();
	// if (n != _n)
	// 	throw IncorrectCodewordException("The codeword is not from a code with given check matrix");

	std::vector<int> result(n, 0);
	std::vector<double> bits_values(n, 0);
	std::vector<std::map<int, int>> alpha(_m, std::map<int, int>());
	std::vector<std::map<int, double>> beta(_m, std::map<int, double>());
	std::vector<std::map<int, double>> gamma(_m, std::map<int, double>());

	// Init
	for (size_t j = 0; j < _m; j++)
	{
		for (auto &i : _checks[j])
		{
			alpha[j][i] = spSign(llr[i]);
			beta[j][i] = phi(fabs(llr[i]));
		}
	}
	size_t iteration = 0;

	while (iteration < _iterationsCount)
	{
		iteration++;
		HorizontalStep(alpha, beta, gamma);

		 // Result of iteration
		for (size_t i = 0; i < n; i++)
		{
			bits_values[i] = llr[i];
			for (auto &j : _bits[i])
			{
				bits_values[i] += gamma[j][i];
			}
		}

		for (size_t i = 0; i < n; i++)
			result[i] = (bits_values[i] <= 0) ? 1 : 0;

		for (size_t j = 0; j < _m; j++)
        {
            int sum = 0;
            for (auto &i : _checks[j])
            {
                sum ^= result[i];
            }
            /*if ((bool)sum) {
                *isFailed = true;
                break;
            }*/
        }


		if (iteration >= _iterationsCount)
			break;

		//Vertical Step
		for (size_t i = 0; i < n; i++)
		{
			double value = bits_values[i];
			for (auto &j : _bits[i])
			{
				double new_value = value - gamma[j][i];
				alpha[j][i] = spSign(new_value);
				beta[j][i] = phi(fabs(new_value));
			}
		}
	}
	return (result);
}
