#pragma once

#include <vector>
#include <map>

#include "Base_decoder.h"
#include "MathOperations.h"

class SP_decoder : public Base_decoder {
    
private:
    
    void HorizontalStep(std::vector<std::map<int, int>> alpha, std::vector<std::map<int, double>> beta, std::vector<std::map<int, double>> &gamma);
public:
    SP_decoder(std::vector<std::vector<int>> H, int iterationsCount);
    std::vector<int> Decode(std::vector<double> llr, bool *isFailed) override;
};