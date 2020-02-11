#include "MathOperations.h"
#include "Exceptions.h"

#include <vector>
using namespace std;

vector<int> operator * (const vector<vector<int>> &A, const vector<int> &x) {
	size_t m = A.size();
	if (m <= 0)
		throw IncorrectMatrixDimensionsException("Matrix has incorrect row size");
	size_t n = A[0].size();

	if (n <= 0)
		throw IncorrectMatrixDimensionsException("Matrix has incorrect column size");

	for (size_t i = 0; i < m; i++)
	{
		if (n != A[0].size())
			throw IncorrectMatrixDimensionsException("Matrix has different column size in 0 and " + to_string(i) + " row");
	}

	if (x.size() != n)
		throw IncorrectDimensionsException("Matrix has incorrect column size");

	vector<int> result(m, 0);
	for (size_t i = 0; i < m; i++)
	{
		for (size_t j = 0; j < n; j++) {
			result[i] += A[i][j] * x[j];
		}
	}

	for (size_t i = 0; i < m; i++)
	{
		result[i] %= 2;
	}

	return result;
}