#include <iostream>
#include "Decoder.h"

std::vector<std::vector<int>> get_matrix() {
	
	int K = 12, M = 6;
	std::vector<std::vector<int>> H(M, std::vector<int>(K, 0));
	H[0][0] = H[0][1] = H[0][2] = H[0][5] = H[0][6] = H[0][10] = 1;
	H[1][0] = H[1][1] = H[1][2] = H[1][3] = H[1][4] = H[1][11] = 1;
	H[2][5] = H[2][6] = H[2][7] = H[2][9] = H[2][10] = H[2][11] = 1;
	H[3][0] = H[3][3] = H[3][7] = H[3][8] = H[3][9] = H[3][11] = 1;
	H[4][1] = H[4][3] = H[4][4] = H[4][6] = H[4][7] = H[4][8] = 1;
	H[5][2] = H[5][4] = H[5][5] = H[5][8] = H[5][9] = H[5][10] = 1;
	return H;
}

int main() {

	auto H = get_matrix();
	auto decoder = Decoder(H, 1);

	auto codeword_llr = vector<double>({ -200.0, 50, -150, -160, -190, -160, -200.0, -198, -150, -160, -190, -160 });
	decoder.Decode(codeword_llr);
	std::cout << "Hello, LDPC!" << std::endl;
	return 0;
}