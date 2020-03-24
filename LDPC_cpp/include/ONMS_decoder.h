#pragma once

#include <vector>
#include <map>

#include "Base_decoder.h"

class ONMS_decoder : Base_decoder {

private:
	const double MinSumNorm = 0.72;
	const double MinSumOffset = 0.0;

	// This vectors were put here in order to get better performance
	std::vector<int> _result;
	std::vector<double> _bits_values;
	std::vector<std::vector<double>> _gamma;
	std::vector<std::vector<double>> _alpha_beta;
	
public:
	ONMS_decoder(std::vector<std::vector<int>> H_row_sparse, int iterationsCount);
	std::vector<int> Decode(std::vector<double> llr, bool * isFailed) override;
};