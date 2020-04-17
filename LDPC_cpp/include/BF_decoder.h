#pragma once

#include <vector>

#include "Base_decoder.h"

class BF_decoder : public Base_decoder {
    
private:
    // each check contain a vector of indices of value bits, which connected to the check
    std::vector<std::vector<int>> _H;
    
public:
    BF_decoder(std::vector<std::vector<int>> H, int iterationsCount);
    std::vector<int> Decode(std::vector<double> llr) override;
};
