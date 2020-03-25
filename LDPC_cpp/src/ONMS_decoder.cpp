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

ONMS_decoder::ONMS_decoder(std::vector<std::vector<int>> H_row_sparse, int iterationsCount, double scale, double offset) 
	: Base_decoder(H_row_sparse, iterationsCount)
{
	_min_sum_scale = scale;
	_min_sum_offset = offset;

	_alpha_beta.resize(_m);
	_gamma.resize(_m);
	_result.resize(_n);
	_bits_values.resize(_n);

	for (size_t j = 0; j < _m; j++)
	{
		_alpha_beta[j].resize(_checks[j].size());
		_gamma[j].resize(_checks[j].size());
	}
}

std::vector<int> ONMS_decoder::Decode(std::vector<double> llr, bool * isFailed) {
	
	size_t n = llr.size();
	if (n != _n)
		throw IncorrectCodewordException("The codeword is not from a code with given check matrix");

	auto counter = 0;
	// Init
	for (size_t j = 0; j < _m; j++)
	{
		counter = 0;
		for (auto &i : _checks[j])
		{
			_alpha_beta[j][counter++] = llr[i];
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

			counter = 0;
			for (auto &i : _checks[j])
			{
				double abs_value = abs(_alpha_beta[j][counter]);
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

				if (_alpha_beta[j][counter] < 0)
				{
					sum_sign ^= 1;
				}
				counter++;
			}

			counter = 0;
			for (auto &i : _checks[j])
			{
				int sign = sum_sign;
				double abs_values = 0;

				if (_alpha_beta[j][counter] < 0)
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

				double offseted_value = abs_values - _min_sum_offset;
				offseted_value = offseted_value > 0 ? offseted_value : 0;
				_gamma[j][counter] = (1 - 2 * sign) * _min_sum_scale * offseted_value;
				counter++;
			}
		}

		// Result of iteration
		std::vector<size_t> counters_arr1(_m);
		for (size_t i = 0; i < n; i++)
		{
			_bits_values[i] = llr[i];

			for (auto &j : _bits[i])
			{
				_bits_values[i] += _gamma[j][counters_arr1[j]++];
			}
		}
		
		for (size_t i = 0; i < n; i++)
		{
			_result[i] = (_bits_values[i] <= 0) ? 1 : 0;
		}

		*isFailed = false;
		for (size_t j = 0; j < _m; j++)
		{
			int sum = 0;
			for (auto &i : _checks[j])
			{
				sum ^= _result[i];
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
		std::vector<size_t> counters_arr2(_m);
		for (size_t i = 0; i < n; i++)
		{
			double value = _bits_values[i];

			for (auto &j : _bits[i])
			{
				_alpha_beta[j][counters_arr2[j]] = value - _gamma[j][counters_arr2[j]];
				counters_arr2[j]++;
			}
		}
	}

	return _result;
}