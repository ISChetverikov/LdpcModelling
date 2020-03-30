#include "../include/Exceptions.h"
#include "../include/Base_decoder.h"

Base_decoder::Base_decoder() {
}

Base_decoder::Base_decoder(std::vector<std::vector<int>> H_row_sparse, int iterationsCount) {

	if ((_iterationsCount = iterationsCount) <= 0)
		throw std::invalid_argument("Number of iterations is incorrect: " + std::to_string(_iterationsCount));

	size_t m = H_row_sparse.size();
	if (m <= 0)
		throw IncorrectMatrixDimensionsException("Check matrix has incorrect row size");


	_checks = H_row_sparse;

	int n = 0;
	for (size_t i = 0; i < m; i++)
	{
		int max = *max_element(H_row_sparse[i].begin(), H_row_sparse[i].end());
		if (max > n)
			n = max;
	}
	n++;

	_bits.resize(n, std::vector<int>());
	for (size_t j = 0; j < m; j++)
	{
		for (size_t i = 0; i < _checks[j].size(); i++)
		{
			_bits[_checks[j][i]].push_back((int)j);
		}
	}

	_m = _checks.size();
	_n = _bits.size();

	return;
}

size_t Base_decoder::GetCodewordLegth() {
	return _n;
}

size_t Base_decoder::GetChecksSymbolsCount() {
	return _m;
}