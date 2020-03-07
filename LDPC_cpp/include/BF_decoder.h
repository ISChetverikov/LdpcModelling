#pragma once

#include <vector>

class BF_decoder {

private:
    // each check contain a vector of indices of value bits, which connected to the check
    std::vector<std::vector<int>> _H;
    std::vector<int> _checks;
    std::vector<std::vector<int>> _bits;

    size_t _m = 0;
    size_t _n = 0;

    size_t _iterationsCount = 0;

public:
    BF_decoder(std::vector<std::vector<int>> H, int iterationsCount);
    std::vector<int> Decode(std::vector<int> llr, int max_iter);
};