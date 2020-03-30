#include "../include/MatrixReading.h"
#include "../include/Exceptions.h"

std::vector<std::vector<int>> readAsRowSparseMatrix(std::string filename) {

	auto idx = filename.rfind('.');

	if (idx == std::string::npos)
		throw ExtensionException("File of matrix has no extension");

	std::string extension = filename.substr(idx + 1);
	if (extension == "sprs")
		return readSprsAsRowSparseMatrix(filename);
	else if (extension == "csv")
		return readCsvAsRowSparseMatrix(filename);
	else
		throw ExtensionException("Extension \"" + extension + "\n is not supported");
}


std::vector<std::vector<int>> readCsvAsRowSparseMatrix(std::string filename) {
	size_t m;
	size_t n;

	return readCsvAsRowSparseMatrix(filename, &m, &n);
}

std::vector<std::vector<int>> readCsvAsRowSparseMatrix(std::string filename, size_t * m, size_t * n) {
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

	getMatrixShape(matrix, *m, *n);

	return matrix;
}

std::vector<std::vector<int>> readSprsAsRowSparseMatrix(std::string filename) {
	size_t m;
	size_t n;

	return readSprsAsRowSparseMatrix(filename, &m, &n);
}

// sprs - just triples "row col val", numeration from 1
std::vector<std::vector<int>> readSprsAsRowSparseMatrix(std::string filename, size_t * m, size_t * n) {
	std::vector<std::vector<int>> matrix;
	std::string line;
	std::ifstream myFile(filename);

	int val;
	int previousRowIndex = -1;
	int rowIndex = 0;
	int columnIndex = 0;
	std::vector<int> temp;

	// First line
	std::getline(myFile, line);
	std::stringstream ss_first(line);
	ss_first >> rowIndex;
	rowIndex--; // matlab numeration
	ss_first >> columnIndex;
	columnIndex--;

	ss_first >> val;
	if (val != 1)
		throw new NotBinaryMatrixException("Sparse matrix elements not from GF(2)");

	temp.push_back(columnIndex);
	previousRowIndex = rowIndex;

	while (std::getline(myFile, line))
	{
		std::stringstream ss(line);
		
		ss >> rowIndex;
		rowIndex--;

		ss >> columnIndex;
		columnIndex--;

		if (rowIndex == previousRowIndex) {
			temp.push_back(columnIndex);
		}
		else {
			matrix.push_back(temp);
			temp.clear();
			temp.push_back(columnIndex);
		}

		ss >> val;
		if (val != 1)
			throw new NotBinaryMatrixException("Sprarse matrix elements not from GF(2)");

		if (previousRowIndex + 1 < rowIndex)
			throw new MatrixRowSkippedException("Sparse matrix has skipped row");

		previousRowIndex = rowIndex;
	}

	matrix.push_back(temp);
	temp.clear();
	temp.push_back(columnIndex);

	getMatrixShape(matrix, *m, *n);

	return matrix;
}

void getMatrixShape(std::vector<std::vector<int>> rowSparseMatrix, size_t& m, size_t& n) {

	size_t M_temp = rowSparseMatrix.size();
	size_t N_temp = 0;
	for (size_t i = 0; i < M_temp; i++)
	{
		size_t max = *max_element(rowSparseMatrix[i].begin(), rowSparseMatrix[i].end());
		if (max > N_temp)
			N_temp = max;
	}
	N_temp++;

	m = M_temp;
	n = N_temp;
}