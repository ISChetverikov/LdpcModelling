#pragma once

#include <vector>
#include <map>

#include "../include/Exceptions.h"
#include "../include/MathOperations.h"

class SP_decoder {
private:
	// each check contain a vector of indices of value bits, which connected to the check
	std::vector<std::vector<int>> _H;
	std::vector<std::vector<int>> _checks;
	std::vector<std::vector<int>> _bits;
	

	size_t _m = 0;
	size_t _n = 0;

	size_t _iterationsCount = 0;

	double MinSumFucntion(std::vector<double> vector);
	void HorizontalStep(std::vector<std::map<int, int>> alpha,
	                    std::vector<std::map<int, double>> beta,
						std::vector<std::map<int, double>> &gamma);
public:
	SP_decoder(std::vector<std::vector<int>> H, int iterationsCount);
	std::vector<int> Decode(std::vector<double> llr, bool *is_success);
};
