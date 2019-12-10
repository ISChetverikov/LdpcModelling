#include <vector>
#include <string>
#include <random>
#include <stdio.h>
#include <math.h>

#include "Decoder.h"

std::vector<std::vector<int>> read_ldpc_mat() {
	//TODO: read H from ldpc_filename
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

void simulate(int maxIterationsCount, std::vector<double> snr_array, double fe, double min_fer, int min_tests) {
	
	// What is M?
	const int M = 2;

    std::vector<std::vector<int>> H = read_ldpc_mat();
	Decoder decoder = Decoder(H, maxIterationsCount);
	int m = H.size(), n = H[0].size();
	int k = n - m;
    double R = (1.0 * k) / n;

    std::vector<double> in_ber(snr_array.size(), 0);
    std::vector<double> ber(snr_array.size(), 0);
    std::vector<double> fer(snr_array.size(), 0);

    for (size_t ii = 0; ii < snr_array.size(); ii++) {
    	double snr = snr_array[ii];
    	double sigma = sqrt(pow(10, -snr / 10) / 2); 
    	printf("\n\n========== SNR = %f, EbN0 = %f, sigma = %f ===========\n", snr, snr - 10 * log10(R * log2(M)), sigma);
    	int tests = 0;
        int wrong_dec = 0;
        int in_errors = 0;
        int errors = 0;

        while ((tests <= min_tests) && (wrong_dec < fe)) {
            tests = ++tests;
            
            // generate inf. word 
            //std::vector<int> iwd(k, 0); //randi(q, 1, K) - 1;

            // encode by outer code
			// Use yet only null-vector
            std::vector<int> codeword(n, 0); //ldpc_encode(iwd, G, 2);

            // add noise
            // TODO: complex 

            std::vector<int> codeword_bpsk(n, 0);
            for (int i = 0; i < n; i++) {
            	codeword_bpsk[i] = 2 * codeword[i] - 1;
            }

            std::default_random_engine generator;
            generator.seed(100);
            std::normal_distribution<double> distribution(0, sigma);
            std::vector<double> noise_vector(n, 0);
            for (int i = 0; i < n; i++) {
            	noise_vector[i] = distribution(generator);
            }
            
            std::vector<double> received_vector(n, 0);
            for (int i = 0; i < n; i++) {
            	received_vector[i] =  codeword_bpsk[i] + noise_vector[i];
            }

            std::vector<double> in_llrs(n, 0);
            std::vector<int> est_cwd(n, 0);
            for (int i = 0; i < n; i++) {
            	in_llrs[i] = -2 * received_vector[i] / (sigma * sigma);
            	est_cwd[i] = (in_llrs[i] < 0);
            }
            for (int i = 0; i < n; i++) {
            	if (codeword[i] != est_cwd[i]) {
            		in_errors++;
            	}
            }
            
            //std::vector<double> scale_array(K + N, scale);
            //std::vector<double> offset_array(K + N, 0);
            
            // puncture
            //for (auto node: punc_nodes) {
            //	in_llrs[node] = 0;
            //}

            // decode
            //if (decoder_number_stub == 1) { //ONMS
            //	continue; //TODO
            //}
            /*switch decoder
                case 1 % Sum-Product
                    [result, number_of_iter, est_cwd, out_llrs] = decode_soft(decoder, ldpc, in_llrs, iter);
                case 2 % Layered Sum-Product
                    [result, number_of_iter, est_cwd, out_llrs] = decode_soft(decoder, ldpc, in_llrs, iter, row_sequence);
                case 3 % Min-Sum
                    [result, number_of_iter, est_cwd, out_llrs] = decode_soft(decoder, ldpc, in_llrs, iter, scale_array, offset_array);
                case 4 % Layered Min-Sum
                    [result, number_of_iter, est_cwd, out_llrs] = decode_soft(decoder, ldpc, in_llrs, iter, row_sequence, scale_array, offset_array);
                case 5 % Adjusted Min-Sum
                    [result, number_of_iter, est_cwd, out_llrs] = decode_soft(decoder, ldpc, in_llrs, iter);
                case 6 % Layered Adjusted Min-Sum
                    [result, number_of_iter, est_cwd, out_llrs] = decode_soft(decoder, ldpc, in_llrs, iter, row_sequence);
                case 7 % Adjusted Min-Sum (Enhanced)
                    [result, number_of_iter, est_cwd, out_llrs] = decode_soft(decoder, ldpc, in_llrs, iter);
                case 8 % Layered Adjusted Min-Sum (Enhanced)
                    [result, number_of_iter, est_cwd, out_llrs] = decode_soft(decoder, ldpc, in_llrs, iter, row_sequence);
                case 9 % Bit-flipping
                    in_est_cwd = est_cwd;
                    [result, number_of_iter, est_cwd] = decode_soft(decoder, ldpc, in_est_cwd, iter);
            end*/
			std::vector<int> decoded_vector = decoder.Decode(in_llrs);

            //std::vector<int> est_iwd;
            //for (int i = N + 1; i < K + N; i++) {
            //	est_iwd.push_back(est_cwd[i]); // information bits are stored in the last positions of cwd
            //}

            std::vector<int> error_indexes; 
            for (int i = 0; i < n; i++) {
            	if (codeword[i] != decoded_vector[i]) {
            		error_indexes.push_back(i);
            	}
            }

			if (error_indexes.empty())
				continue;

            wrong_dec = ++wrong_dec;
            errors = errors + error_indexes.size();
            in_ber[ii] = in_errors / n / tests;
            fer[ii] = wrong_dec / tests;
            ber[ii] = errors / k / tests;
            printf("in_ber = %f, fer = %f, ber = %f\n", in_ber[ii], fer[ii], ber[ii]);
            //save(sprintf('BF_result_M=%d_ldpc=%s_dectype=%d_scale=%f_iter=%d.mat', M, ldpc_filename, decoder, 0.75, iter));
            
			if (fer[ii] < min_fer)
				break;
        }
		printf("fer: %f, tests numbers: %d", fer[ii], tests);
        //save(sprintf('BF_result_M=%d_ldpc=%s_dectype=%d_scale=%f_iter=%d.mat', M, ldpc_filename, decoder, 0.75, iter)); 
    }
}


int main() {
	std::vector<double> snr_array(1, -1.8127);
	simulate(10, snr_array, 30, 1, 1000);
	return 0;
}