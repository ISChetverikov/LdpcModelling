#include "../include/FFH_simulate.h"


double FFH_simulate::normal_pdf(double x, double m, double s) {
    static const double inv_sqrt_2pi = 0.3989422804014327;
    double a = (x - m) / s;
    return inv_sqrt_2pi / s * std::exp(-0.5 * a * a);
}


double FFH_simulate::loss_func(const std::vector<double>& z, const std::vector<int>& codeword) {
    int n = codeword.size();
    double loss_sum = 0;
    for (size_t l = 0; l < n; l++) {
        int q = pow(-1, codeword[l]);
        loss_sum += pow((q * z[l] < 0) * z[l], 2);
    }
    return sqrt(loss_sum / n);
}


bool FFH_simulate::hist_is_flat(std::vector<std::vector<int>>& H, int iter) {
    for (size_t i = 0; i < H.size(); i++) {
        printf("%d ", H[i][iter]);
    }
    printf("\n");
    int sum = 0;
    for (size_t bin_num = 0; bin_num < H[iter].size(); ++bin_num) {
        sum += H[bin_num][iter];
    }
    double floor_num = (double) sum / H[iter].size() * _perc;
    for (size_t bin_num = 0; bin_num < H[iter].size(); ++bin_num) {
        if (H[bin_num][iter] < floor_num)
            return false;
    }
    return true;
}


std::pair<double, double> FFH_simulate::find_opt_V(int L, int snr, const std::vector<int>& codeword,
                                                   double sigma, double f) {
    double Vmin = 1, Vmax = 0;
    std::vector<double> prob(L, log(1.0 / L));
    std::vector<double> z(_n, 0);
    int cur_bin = 0;
    std::vector<int> H(L, 0);
    /*for (size_t it = 0; it < EbN0(snr, 1 - (double) n / m) * L; ++it) {
        dz = np.random.normal(0, sigma, n)
        for l in range(n):
            new_z_l_bit = z[l] + eps * dz[l]
            prob_bit_acceptance = np.min([norm(0, sigma).pdf(new_z_l_bit) / norm(0, sigma).pdf(z[l]), 1])
            if uniform.rvs() <= prob_bit_acceptance:
                new_z[l] = new_z_l_bit
        new_loss = loss_func(new_z, codeword)
        new_bin = int(np.floor((new_loss - Vmin) * L))
        new_bin = np.min([new_bin, L - 1])

        prob_stat_acceptance = np.min([np.exp(prob[cur_bin]) / np.exp(prob[new_bin]), 1])
        if uniform.rvs() <= prob_stat_acceptance:
            z = new_z
            cur_bin = new_bin

        prob[cur_bin] += np.log(f)
        if it > 10:
            H[cur_bin] += 1
        if it % 10 == 0:
            print(it)
    min_bin, max_bin = -1, -1
    for bin_num in range(L):
        if H[bin_num] != 0:
            max_bin = bin_num
            if min_bin == -1:
                min_bin = bin_num
    }
    return min_bin / L, (max_bin + 1) / L*/
    return std::make_pair(0.0, 1.0);
}


void FFH_simulate::simulate() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::vector<int> codeword(_n, 0);
	std::vector<double> llrs(_n, 0);
    bool isFailed = false;

    for (size_t ii = 0; ii < _snr_array.size(); ii++) {
        double sigma = get_sigma(_snr_array[ii]);
        printf("%f\n", sigma);
        int iterations_count = 0;
        
        int L = (int) std::round(pow(10, 1 / sigma) * _n / 10);
        printf("%d\n", L);
        std::vector<std::vector<int>> H, G;
        for (size_t itn = 0; itn < L; itn++) {
            H.push_back(std::vector<int>(_num_iterations, 0));
            G.push_back(std::vector<int>(_num_iterations, 0));
        }
        double f = std::exp(1);
        
        std::pair<double, double> V = find_opt_V(L, _snr_array[ii], codeword, sigma, f);
        double Vmin = V.first, Vmax = V.second;
        printf("%f %f\n", Vmin, Vmax);
        
        std::vector<double> prob(L, log(1.0 / L));
        std::vector<double> z(_n, 0);
        std::normal_distribution<double> norm_distr(0, sigma);
        std::uniform_real_distribution<double> un_distr(0.0, 1.0);
        int cur_bin = 0;
        for (size_t it = 0; it < _num_iterations; ++it) {
            bool is_flat = false;
            while (!is_flat) {
                /*for (size_t i = 0; i < z.size(); ++i) {
                    printf("%f ", z[i]);
                }
                printf("\n");*/
                std::vector<double> new_z = z;
                for (size_t l = 0; l < _n; ++l) {
                    double dz = norm_distr(gen);
                    double new_z_l_bit = z[l] + _eps * dz;
                    double prob_bit_acceptance = std::min(normal_pdf(new_z_l_bit, 0, sigma) / 
                                                          normal_pdf(z[l], 0, sigma), 
                                                          1.0);
                    if (un_distr(gen) <= prob_bit_acceptance)
                        new_z[l] = new_z_l_bit;
                }
                double new_loss = loss_func(new_z, codeword);
                int new_bin = (int) floor((new_loss - Vmin) / (Vmax - Vmin) * (L - 1));
                new_bin = std::min(new_bin, L - 1);
                new_bin = std::max(new_bin, 0);

                double prob_stat_acceptance = std::min(std::exp(prob[cur_bin]) / std::exp(prob[new_bin]), 1.0);
                //print("prob"+str(np.exp(prob[cur_bin]) / np.exp(prob[new_bin])))
                if (un_distr(gen) <= prob_stat_acceptance) {
                    z = new_z;
                    cur_bin = new_bin;
                }
                prob[cur_bin] += log(f);

                for (size_t i = 0; i < _n; i++) {
				    llrs[i] = -2 * (2 * codeword[i] - 1 + z[i]) / (sigma * sigma);
                }
                _decoder->Decode(llrs, &isFailed);
                H[cur_bin][it] += 1;
                if (isFailed) {
                    printf("new error!\n");
                    G[cur_bin][it] += 1;
                }
                iterations_count += 1;
                int check_const = L;
                if (iterations_count == check_const) {
                    if (hist_is_flat(H, it))
                        is_flat = true;
                    iterations_count = 0;
                }
            }
            f = sqrt(f);
        }
        double prob_error = 0;
        double prob_sum = 0;
        for (auto p: prob){
            prob_sum += std::exp(p);
        }
        for (size_t bin_num = 0; bin_num < L; ++bin_num) {
            int el_num = 0, er_num = 0;
            for (size_t it_num = 0; it_num < _num_iterations; ++it_num) {
                el_num += H[bin_num][it_num];
                er_num += G[bin_num][it_num];
            }
            prob_error += std::exp(prob[bin_num]) / prob_sum * er_num / el_num;
        }
        _fer_array.push_back(prob_error);
    }
}


FFH_simulate::FFH_simulate(std::string matrix_name,
                           std::string decoder,
                           int num_iterations,
	                       const std::vector<double>& snr_array) {
    
    _H = readAsRowSparseMatrix(matrix_name);
    _decoder = new ONMS_decoder(_H, 20);

	// size of sparse Matrix
	_n = 0;
	for (size_t i = 0; i < _H.size(); i++)
	{
		int max = *max_element(_H[i].begin(), _H[i].end());
		if (max > _n)
			_n = max;
	}
	_n++;
    _m = _H.size();

	_snr_array = snr_array;
    _num_iterations = num_iterations;

    _eps = 4.0;
    _perc = 0.03;
}