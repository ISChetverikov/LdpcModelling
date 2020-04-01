#include <vector>
#include <map>
#include <cmath>
#include <string>
#include <iostream>
#include <algorithm>

#include "../include/Exceptions.h"
#include "../include/BF_decoder.h"


// Constructor
BF_decoder::BF_decoder(std::vector<std::vector<int>> H, int iterationsCount) : Base_decoder(H, iterationsCount) {
    
}

std::vector<int> BF_decoder::Decode(std::vector<double> llr, bool * isFailed) {
    std::vector<int> _llr;
    for (auto i : llr) {
        if (i >= 0) {
            _llr.push_back(0);
        }
        else {
            _llr.push_back(1);
        }
    }
    size_t n = _llr.size();
    if (n != _n)
        throw IncorrectCodewordException("The codeword is not from a code with given check matrix");
    std::vector<int> result(n);
    std::vector<int> check_values(_m, 0);
    result = _llr;
    for (size_t i = 0; i < n; i++) { // initial checks compute
        if (_llr[i] == 0) continue;
        for (size_t j = 0; j < _bits[i].size(); j++) {
            check_values[_bits[i][j]] = (check_values[_bits[i][j]] + 1) % 2;
        }
    }
    int iter_num = 0;
    bool full_cycle_flg = true;
    while (iter_num < _iterationsCount and full_cycle_flg) {
        iter_num++;
        full_cycle_flg = false;
        for (size_t i = 0; i < n; i++) {
            int sum_incorrect = 0;
            for (size_t j = 0; j < _bits[i].size(); j++) {
                sum_incorrect += check_values[_bits[i][j]];
            }
            if (sum_incorrect > int(_bits[i].size() / 2)) {
                result[i] = (result[i] + 1) % 2;
                for (size_t k = 0; k < _bits[i].size(); k++) {
                    check_values[_bits[i][k]] = (check_values[_bits[i][k]] + 1) % 2;
                    full_cycle_flg = true;
                }
            }
        }
    }
    return result;
}
