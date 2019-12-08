#include <vector>
#include <string>
#include <random>
#include <stdio.h>
#include <math.h>

std::vector<std::vector<int>> read_ldpc_mat(const std::string ldpc_filename) {
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

void simulate(int M, std::string ldpc_filename, int decoder, int iter, std::vector<int> punc_nodes,
	          std::vector<double> snr_array, double fe, double scale, double min_fer, int min_tests) {
	int K = 12, N = 6;
    std::vector<std::vector<int>> H = read_ldpc_mat(ldpc_filename);
    double R = (1.0 * K) / (K + N);

    std::vector<double> in_ber(snr_array.size(), 0);
    std::vector<double> ber(snr_array.size(), 0);
    std::vector<double> fer(snr_array.size(), 0);

    for (int ii = 0; ii < snr_array.size(); ii++) {
    	double snr = snr_array[ii];
    	double sigma = sqrt(pow(10, -snr / 10) / 2); 
    	printf("\n\n========== SNR = %f, EbN0 = %f, sigma = %f ===========\n", snr, snr - 10 * log10(R * log2(M)), sigma);
    	int tests = 0;
        int wrong_dec = 0;
        int in_errors = 0;
        int errors = 0;

        while ((tests <= min_tests) || (wrong_dec < fe)) {
            tests = tests + 1;
            
            // generate inf. word 
            std::vector<int> iwd(K, 0); //randi(q, 1, K) - 1;

            // encode by outer code
            std::vector<int> cwd(K + N, 0); //ldpc_encode(iwd, G, 2);

            // add noise
            // TODO: complex 

            std::vector<int> tx(K + N, 0);
            for (int i = 0; i < K + N; i++) {
            	tx[i] = 2 * cwd[i] - 1;
            }

            std::default_random_engine generator;
            generator.seed(100);
            std::normal_distribution<double> distribution(0, sigma);
            std::vector<double> noise_vector(K + N, 0);
            for (int i = 0; i < K + N; i++) {
            	noise_vector[i] = distribution(generator);
            }
            
            std::vector<double> rx(K + N, 0);
            for (int i = 0; i < K + N; i++) {
            	rx[i] =  tx[i] + noise_vector[i];
            }

            // demodulate
            std::vector<double> in_llrs(K + N, 0);
            std::vector<int> est_cwd(K + N, 0);
            for (int i = 0; i < K + N; i++) {
            	in_llrs[i] = -2 * rx[i] / (sigma * sigma);
            	est_cwd[i] = (in_llrs[i] < 0);
            }
            for (int i = 0; i < K + N; i++) {
            	if (cwd[i] != est_cwd[i]) {
            		in_errors++;
            	}
            }
            
            
            std::vector<double> scale_array(K + N, scale);
            std::vector<double> offset_array(K + N, 0);
            
            // puncture
            //for (auto node: punc_nodes) {
            //	in_llrs[node] = 0;
            //}

            // decode
            if (decoder == 1) { //ONMS
            	continue; //TODO
            }
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

            std::vector<int> est_iwd;
            for (int i = N + 1; i < K + N; i++) {
            	est_iwd.push_back(est_cwd[i]); // information bits are stored in the last positions of cwd
            }
            std::vector<int> error_indexes; 
            for (int i = 0; i < K; i++) {
            	if (iwd[i] != est_iwd[i]) {
            		error_indexes.push_back(i);
            	}
            }
            if (!error_indexes.empty()) {
                wrong_dec = wrong_dec + 1;
                errors = errors + error_indexes.size();
                in_ber[ii] = in_errors / (N + K) / tests;
                fer[ii] = wrong_dec / tests;
                ber[ii] = errors / K / tests;
                printf("in_ber = %f, fer = %f, ber = %f\n", in_ber[ii], fer[ii], ber[ii]);
                //save(sprintf('BF_result_M=%d_ldpc=%s_dectype=%d_scale=%f_iter=%d.mat', M, ldpc_filename, decoder, 0.75, iter));
            }
        }
        if (fer[ii] < min_fer)
            break;
        //save(sprintf('BF_result_M=%d_ldpc=%s_dectype=%d_scale=%f_iter=%d.mat', M, ldpc_filename, decoder, 0.75, iter)); 
    }
}


int main() {
	std::vector<double> snr_array(1, 1.5);
	simulate(2, "ll", 1, 1, std::vector<int>(3, 0), snr_array, 10, 10, 10, 10);
	return 0;
}