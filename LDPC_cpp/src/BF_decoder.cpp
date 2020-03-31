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
    result = _llr;
    /*for(size_t i=0; i < n; i++){
     if(_llr[i] == 0) continue;
     for(size_t j=0; j<_bits[i].size(); j++){
     std::cout << _bits[i][j] << std::endl;
     _checks[_bits[i][j]] = (_checks[_bits[i][j]] + 1) % 2;
     }
     }
     int iter_num = 0;
     int sum_incorect = 0;
     std::cout << "tttt";
     int full_cycle_flg = 1;
     while(iter_num<50 and full_cycle_flg==1){
     full_cycle_flg=0;
     for(size_t i=0; i<n; i++){
     iter_num++;
     sum_incorect = 0;
     for(size_t j=0; j<_bits[i].size(); j++)
     sum_incorect+= _checks[_bits[i][j]];
     if (sum_incorect > int(_bits[i].size()/2)){
     result[i] = (result[i] + 1) % 2;
     for(size_t k=0; k<_bits[i].size(); k++)
     _checks[_bits[i][k]] = (_checks[_bits[i][k]] + 1) % 2;
     full_cycle_flg = 1;
     }
     }
     }
     */
    return result;
}
