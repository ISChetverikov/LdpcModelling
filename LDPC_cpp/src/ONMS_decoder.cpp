#include <vector>
#include <map>
#include <unordered_map>
#include <cmath>
#include <string>
#include <iostream>
#include <algorithm>

#include "../include/ONMS_decoder.h"
#include "../include/MathOperations.h"
#include "../include/Exceptions.h"

#define DBL_MAX 1.7976931348623158e+308 

std::vector<int> ONMS_decoder::Decode(std::vector<double> llr, bool * isFailed) {
	
	size_t n = llr.size();
	if (n != _n)
		throw IncorrectCodewordException("The codeword is not from a code with given check matrix");

	std::vector<int> result(n);
	std::vector<double> bits_values(n);
	
	std::vector<std::map<int, double>> alpha_beta(_m, std::map<int, double>());
	std::vector<std::map<int, double>> gamma(_m, std::map<int, double>());

	// Init
	for (size_t j = 0; j < _m; j++)
	{
		for (auto &i : _checks[j])
		{
			alpha_beta[j][i] = llr[i];
		}	
	}

	size_t iteration = 0;
	while (true)
	{
		iteration++;

		// Horizontal step
		for (size_t j = 0; j < _m; j++)
		{
			double first_min = DBL_MAX;
			double second_min = DBL_MAX;
			int min_index = 0;
			int sum_sign = 0;

			for (auto &i : _checks[j])
			{
				double abs_value = abs(alpha_beta[j][i]);
				if (abs_value < first_min)
				{
					second_min = first_min;
					first_min = abs_value;
					min_index = i;
				}
				else if (abs_value < second_min)
				{
					second_min = abs_value;
				}

				if (alpha_beta[j][i] < 0)
				{
					sum_sign ^= 1;
				}
			}

			for (auto &i : _checks[j])
			{
				int sign = sum_sign;
				double abs_values = 0;

				if (alpha_beta[j][i] < 0)
				{
					sign ^= 1;
				}
				if (i == min_index)
				{
					abs_values = second_min;
				}
				else
				{
					abs_values = first_min;
				}

				double offseted_value = abs_values - MinSumOffset;
				offseted_value = offseted_value > 0 ? offseted_value : 0;
				gamma[j][i] = (1 - 2 * sign) * MinSumNorm * offseted_value;
			}
		}

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
		{
			result[i] = (bits_values[i] <= 0) ? 1 : 0;
		}

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
				alpha_beta[j][i] = value - gamma[j][i];
			}
		}
	}

	return result;
}