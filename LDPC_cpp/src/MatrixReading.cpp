#include "../include/MatrixReading.h"

std::vector<std::vector<int>> readAsRowSparseMatrix(std::string filename) {
	size_t m;
	size_t n;

	return readAsRowSparseMatrix(filename, &m, &n);
}

std::vector<std::vector<int>> readAsRowSparseMatrix(std::string filename, size_t * m, size_t * n) {
	std::vector<std::vector<int>> matrix;
	std::string line;
	std::ifstream myFile(filename);

	int val;

	while (std::getline(myFile, line))
	{
		std::stringstream ss(line);

		std::vector<int> temp;
		int columnIndex = 0;
		while (ss >> val) {
			if (val)
				temp.push_back(columnIndex);

			if (ss.peek() == ',')
				ss.ignore();

			columnIndex++;
		}

		matrix.push_back(temp);
	}

	size_t M_temp = matrix.size();
	size_t N_temp = 0;
	for (size_t i = 0; i < M_temp; i++)
	{
		size_t max = *max_element(matrix[i].begin(), matrix[i].end());
		if (max > N_temp)
			N_temp = max;
	}
	N_temp++;

	*m = M_temp;
	*n = N_temp;
	return matrix;
}