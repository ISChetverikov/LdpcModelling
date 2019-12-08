#include <vector>

using namespace std;

class Decoder {

private:
	vector<vector<int>> Check;
	vector<vector<int>> Values;

public:
	Decoder(vector<vector<int>> H);
	vector<int> Decode(vector<double> llr);
};