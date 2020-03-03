#include <iostream>
#include <vector>
#include "../include/Decoder.h"
#include "../include/MathOperations.h"


using namespace std;
std::vector<std::vector<int>> get_matrix() {
	
	int K = 12, M = 6;
	std::vector<std::vector<int>> H(M);
	H[0] = { 0,1,2,5,6,10 };
	H[1] = { 0,1,2,3,4,11 };
	H[2] = { 5,6,7,9,10,11 };
	H[3] = { 0,3,7,8,9,11};
	H[4] = { 1,3,4,6,7,8 };
	H[5] = { 2,4,5,8,9,10 };

	/*H[0][0] = H[0][1] = H[0][2] = H[0][5] = H[0][6] = H[0][10] = 1;
	H[1][0] = H[1][1] = H[1][2] = H[1][3] = H[1][4] = H[1][11] = 1;
	H[2][5] = H[2][6] = H[2][7] = H[2][9] = H[2][10] = H[2][11] = 1;
	H[3][0] = H[3][3] = H[3][7] = H[3][8] = H[3][9] = H[3][11] = 1;
	H[4][1] = H[4][3] = H[4][4] = H[4][6] = H[4][7] = H[4][8] = 1;
	H[5][2] = H[5][4] = H[5][5] = H[5][8] = H[5][9] = H[5][10] = 1;*/
	return H;
}

int main1() {
	bool result;
	auto H = get_matrix();
	auto decoder = new Decoder(H, 20);

	auto codeword_llr = vector<double>({ -20.0, -20, -15, -16, -19, -16, -20, -19, -15, -16, -19, -16 });
	auto decoded = decoder->Decode(codeword_llr, &result);

	for (size_t i = 0; i < decoded.size(); i++)
	{
		cout << decoded[i] << " ";
	}

	return 0;
}