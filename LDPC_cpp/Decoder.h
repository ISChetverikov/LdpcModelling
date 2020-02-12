#include <vector>
#include <map>

using namespace std;

class Decoder {

private:
	// each check contain a vector of indices of value bits, which connected to the check
	vector<vector<int>> _H;
	vector<vector<int>> _checks;
	vector<vector<int>> _bits;
	
	const double MinSumNorm = 0.72;
	const double MinSumOffset = 0.5;

	size_t _m = 0;
	size_t _n = 0;

	size_t _iterationsCount = 0;

	double MinSumFucntion(vector<double> vector);
	void HorizontalStep(vector<map<int, int>> alpha, vector<map<int, double>> beta, vector<map<int, double>> &gamma);
public:
	Decoder(vector<vector<int>> H, int iterationsCount);
	vector<int> Decode(vector<double> llr, bool * isSuccess);
};