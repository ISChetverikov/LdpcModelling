#include "pch.h"
#include "../../LDPC_cpp/include/ONMS_decoder.h"
#include "../../LDPC_cpp/include/Base_decoder.h"

class OnmsDecoder_SmallMatrix : public ::testing::Test
{
protected:
	void SetUp()
	{
		std::vector<std::vector<int>> H(M);
		H[0] = { 0,1,2,5,6,10 };
		H[1] = { 0,1,2,3,4,11 };
		H[2] = { 5,6,7,9,10,11 };
		H[3] = { 0,3,7,8,9,11 };
		H[4] = { 1,3,4,6,7,8 };
		H[5] = { 2,4,5,8,9,10 };

		decoder = new ONMS_decoder(H, 20);
	}
	void TearDown()
	{
		return;
	}

	std::vector<double> ToLlr(std::vector<int> bits) {
		std::vector<double> result = std::vector<double>(bits.size());

		for (size_t i = 0; i < bits.size(); i++)
		{
			result[i] = 1 - 2 * bits[i];
		}

		return result;
	}
	
	int N = 12;
	int M = 6;
	ONMS_decoder* decoder;
	std::vector<int> ones_vector = std::vector<int>(N, 1);
	std::vector<int> zeros_vector = std::vector<int>(N, 0);
};

TEST_F(OnmsDecoder_SmallMatrix, OnesVectorWithoutError) {

	bool isFailed = false;
	auto codeword_llr = std::vector<double>({ -20, -20, -15, -16, -19, -16, -20, -19, -15, -16, -19, -16 });
	auto result = decoder->Decode(codeword_llr, &isFailed);
	
	ASSERT_TRUE(result == ones_vector);
}

TEST_F(OnmsDecoder_SmallMatrix, ZerosVectorWithoutError) {

	bool isFailed = false;
	auto codeword_llr = std::vector<double>(N, 1);
	auto result = decoder->Decode(codeword_llr, &isFailed);

	ASSERT_TRUE(result == zeros_vector);
}

TEST_F(OnmsDecoder_SmallMatrix, ZerosVectorWithOneError) {

	bool isFailed = false;

	for (size_t i = 0; i < N; i++)
	{
		auto codeword_llr = std::vector<double>(N, 1);
		codeword_llr[i] = -1;
		auto result = decoder->Decode(codeword_llr, &isFailed);

		ASSERT_TRUE(result == zeros_vector);
	}
	
}