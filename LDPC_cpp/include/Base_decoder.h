#pragma once

#include <vector>
#include <algorithm>
#include <stdexcept>

class Base_decoder {

protected:
	// each check contain a vector of indices of value bits, which connected to the check
	std::vector<std::vector<int>> _checks;
	std::vector<std::vector<int>> _bits;

	size_t _m = 0;
	size_t _n = 0;
	size_t _iterationsCount = 0;

public:
	// Input number of ones' positions in each row
	Base_decoder();
	Base_decoder(std::vector<std::vector<int>> H, int iterationsCount);

	size_t GetCodewordLegth();
	size_t GetChecksSymbolsCount();

	virtual std::vector<int> Decode(std::vector<double> llr) = 0;
	virtual ~Base_decoder() {};
};