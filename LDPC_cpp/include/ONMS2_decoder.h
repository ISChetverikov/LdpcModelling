#pragma once

#include <vector>
#include <map>

#include "Base_decoder.h"

class ONMS2_decoder : Base_decoder {

private:
	const double MinSumNorm = 0.72;
	const double MinSumOffset = 0.0;

	double MinSumFunction(std::vector<double> values);
	void HorizontalStep(std::vector<std::map<int, int>> alpha,
		std::vector<std::map<int, double>> beta,
		std::vector<std::map<int, double>> &gamma);
public:
	ONMS2_decoder(std::vector<std::vector<int>> H_row_sparse, int iterationsCount)
		: Base_decoder(H_row_sparse, iterationsCount) {};
	std::vector<int> Decode(std::vector<double> llr, bool * isFailed) override;
};