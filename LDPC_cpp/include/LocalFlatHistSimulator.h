#pragma once

#include <string>
#include <chrono>
#include <string> //del
#include "Base_decoder.h"
#include "BaseSimulator.h"

class LocalFlatHistSimulator : public BaseSimulator {
protected:
	double _epsilon;
	int _l;
	int _kMin;
	int _alpha;
	int _beta;
	int _unconWithoutAB;
	int _unconWithAB;
	int _conWithoutAB;
	int _conWithAB;
	
	double normalPDF(double x, double m, double s);
	double lossFunc(const std::vector<double>& z, const std::vector<int>& codeword);
	std::vector<double> findStartPoint();

public:
	LocalFlatHistSimulator(Base_decoder * decoderPtr, double epsilon, int l, int kMin, int alpha, int beta,
	int unconWithoutAB, int unconWithAB, int conWithoutAB, int conWithAB);
	~LocalFlatHistSimulator() {};
	void Run(std::vector<double> snrArray,
		std::vector<double>& ebn0Array,
		std::vector<double>& ferArray,
		std::vector<double>& sigmaArray,
		std::vector<int>& testsCountArray,
		std::vector<std::chrono::milliseconds>& elapsedTimeArray) override;
};