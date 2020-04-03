#include <vector>
#include <map>
#include "Base_decoder.h"
using namespace std;

class SP_decoder : public Base_decoder {
    
private:
    
    void HorizontalStep(vector<map<int, int>> alpha, vector<map<int, double>> beta, vector<map<int, double>> &gamma);
public:
    SP_decoder(vector<vector<int>> H, int iterationsCount);
    vector<int> Decode(vector<double> llr, bool *isFailed) override;
};

