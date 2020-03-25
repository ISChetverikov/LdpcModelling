 #include "../include/FFH_simulate.h"
 #include "../include/MatrixReading.h"


std::vector<double> FFH_simulate::get_fer() {
    return _fer_array;
}

double FFH_simulate::get_sigma(double snr) {
    return sqrt(pow(10, (-snr / 10)) / 2);
}


double FFH_simulate::EbN0(double snr, double R) {
    return snr - 10 * log10(R);
}

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
    std::random_device rd;
    std::mt19937 gen(rd());
    double Vmin = 0, Vmax = 1;
    std::vector<double> prob(L, log(1.0 / L));
    std::vector<double> z(_n, 0);
    std::normal_distribution<double> norm_distr(0, sigma);
    std::uniform_real_distribution<double> un_distr(0.0, 1.0);
    int cur_bin = 0;
    std::vector<int> H(L, 0);
    for (size_t it = 0; it < EbN0(snr, 1 - (double) _m / _n) * L + _skip_iterations; ++it) {
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
        if (un_distr(gen) <= prob_stat_acceptance) {
            z = new_z;
            cur_bin = new_bin;
        }
        prob[cur_bin] += log(f);
        if (it >= _skip_iterations) {
            H[cur_bin] += 1;
        }
    }
    double min_bin = -1.0, max_bin = -1.0;
    for (size_t bin_num = 0; bin_num < L; ++bin_num) {
        if (H[bin_num] != 0) {
            max_bin = bin_num;
            if (min_bin == -1.0) {
                min_bin = bin_num;
            }
        }
    }
    return std::make_pair(min_bin / L, (max_bin + 1) / L);
}


void FFH_simulate::simulate() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::vector<int> codeword(_n, 0);
	std::vector<double> llrs(_n, 0);
    std::vector<int> decoded(_n, 0);
    bool isFailed = false;

    for (size_t ii = 0; ii < _snr_array.size(); ii++) {
        double sigma = get_sigma(_snr_array[ii]);
        //std::cout << sigma << "\n";

        int iterations_count = 0;
        
        int L = (int) std::round(pow(10, 1 / sigma) * _n / 10);
        std::vector<std::vector<int>> H, G;
        for (size_t itn = 0; itn < L; itn++) {
            H.push_back(std::vector<int>(_num_iterations, 0));
            G.push_back(std::vector<int>(_num_iterations, 0));
        }
        double f = std::exp(1);
        int check_const = 10 * L;
        
        std::pair<double, double> V = find_opt_V(L, _snr_array[ii], codeword, sigma, f);
        double Vmin = V.first, Vmax = V.second;
        std::vector<double> prob(L, log(1.0 / L));
        std::vector<double> z(_n, 0);
        std::normal_distribution<double> norm_distr(0, sigma);
        std::uniform_real_distribution<double> un_distr(0.0, 1.0);
        int cur_bin = 0;
        int wrong_dec = 0;
        for (size_t it = 0; it < _num_iterations; ++it) {
            bool is_flat = false;
            while (!is_flat) {
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
                if (un_distr(gen) <= prob_stat_acceptance) {
                    z = new_z;
                    cur_bin = new_bin;
                }
                prob[cur_bin] += log(f);
                for (size_t i = 0; i < _n; i++) {
				    llrs[i] = -2 * (2 * codeword[i] - 1 + z[i]) / (sigma * sigma);
                }
                decoded = _decoder->Decode(llrs, &isFailed);
                H[cur_bin][it] += 1;
                //if (isFailed) {
                if (decoded != codeword) {
                    G[cur_bin][it] += 1;
                    wrong_dec++;
                }
                if (wrong_dec >= _rejections_num) {
                    it = _num_iterations;
                    break;
                }
                iterations_count += 1;
                if (iterations_count == check_const) {
                    if (hist_is_flat(H, it)) {
                        is_flat = true;
                    }
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
            if (el_num != 0) {
                prob_error += std::exp(prob[bin_num]) / prob_sum * er_num / el_num;
            }
        }
        _fer_array.push_back(prob_error);
    }
}


FFH_simulate::FFH_simulate(std::string matrix_name,
                           std::string decoder,
                           int num_iterations,
                           int rejections_num,
	                       const std::vector<double>& snr_array,
                           std::vector<double>& _fer_array) {
    
    _H = readAsRowSparseMatrix(matrix_name, &_m, &_n);
    _decoder = new ONMS_decoder(_H, 20, 0.72, 0);

	_snr_array = snr_array;
    _num_iterations = num_iterations;
    _rejections_num = rejections_num;
    _skip_iterations = 2000;

    _eps = 4.45;
    _perc = 0.07;
}