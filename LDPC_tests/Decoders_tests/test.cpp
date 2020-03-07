#include "pch.h"
#include "../../LDPC_cpp/include/ONMS_decoder.h"
#include "../../LDPC_cpp/include/Base_decoder.h"

TEST(TestCaseName, TestName) {
	int K = 12, M = 6;
	std::vector<std::vector<int>> H(M);
	H[0] = { 0,1,2,5,6,10 };
	H[1] = { 0,1,2,3,4,11 };
	H[2] = { 5,6,7,9,10,11 };
	H[3] = { 0,3,7,8,9,11 };
	H[4] = { 1,3,4,6,7,8 };
	H[5] = { 2,4,5,8,9,10 };

	bool isFailed = false;
	auto codeword_llr = std::vector<double>({ -20.0, -20, -15, -16, -19, -16, -20, -19, -15, -16, -19, -16 });
	auto decoder = new ONMS_decoder(H, 20);
	auto result = decoder->Decode(codeword_llr, &isFailed);

	std::cerr << "[          ] random seed = " << "Hello!" << std::endl;
	EXPECT_EQ(1, 1);
	EXPECT_TRUE(true);
}