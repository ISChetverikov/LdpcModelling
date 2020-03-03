#pragma once

#include <vector>
#include <map>

#include "../include/Exceptions.h"
#include "../include/MathOperations.h"

class Decoder {

private:
	// each check contain a vector of indices of value bits, which connected to the check
	std::vector<std::vector<int>> _checks;
	std::vector<std::vector<int>> _bits;
	
	const double MinSumNorm = 1;// 0.72;
	const double MinSumOffset = 0.0;

	size_t _m = 0;
	size_t _n = 0;

	size_t _iterationsCount = 0;

	double MinSumFunction(std::vector<double> values);
	void HorizontalStep(std::vector<std::map<int, int>> alpha,
	                    std::vector<std::map<int, double>> beta,
						std::vector<std::map<int, double>> &gamma);
public:
	Decoder(std::vector<std::vector<int>> H_row_sparse, int iterationsCount);
	std::vector<int> Decode(std::vector<double> llr, bool * isFailed);
};