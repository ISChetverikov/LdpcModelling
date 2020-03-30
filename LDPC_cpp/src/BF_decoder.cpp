#include <vector>
#include <map>
#include <cmath>
#include <string>
#include <iostream>
#include <algorithm>

#include "../include/Exceptions.h"
#include "../include/BF_decoder.h"


// Constructor
BF_decoder::BF_decoder(std::vector<std::vector<int>> H, int iterationsCount) {
    
    if ((_iterationsCount = iterationsCount) <= 0)
        throw std::invalid_argument("Number of iterations is incorrect: " + std::to_string(_iterationsCount));

    size_t m = H.size();
    if (m <= 0)
        throw IncorrectMatrixDimensionsException("Check matrix has incorrect row size");
    size_t n = H[0].size();

    if (n <= 0)
        throw IncorrectMatrixDimensionsException("Check matrix has incorrect column size");

    for (size_t i = 0; i < m; i++)
    {
        if (n != H[0].size())
            throw IncorrectMatrixDimensionsException("Check matrix has different column size in 0 and " + std::to_string(i) + " row");
    }

    _bits.resize(n, std::vector<int>());
    _checks.resize(m, 0);
    for (size_t i = 0; i < m; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            if (H[i][j] == 0)
                continue;
            
            _bits[j].push_back(i);
        }
    }
    
    _n = _bits.size();
    
    return;
}

std::vector<int> BF_decoder::Decode(std::vector<int> llr, int max_iter){
    size_t n = llr.size();
    if (n != _n)
        throw IncorrectCodewordException("The codeword is not from a code with given check matrix");

    std::vector<int> result(n);
    result = llr;
    for(size_t i=0; i < n; i++){
        if(llr[i] == 0) continue;
        for(size_t j=0; j<_bits[i].size(); j++){
            _checks[_bits[i][j]] = (_checks[_bits[i][j]] + 1) % 2;
        }
    }
    int iter_num = 0;
    int sum_incorect = 0;
    
    int full_cycle_flg = 1;
    while(iter_num<max_iter and full_cycle_flg==1){
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
    
    return result;
}
