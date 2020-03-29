#include <random>
#include <chrono>
#include <iostream> //delete
#include "../include/LocalFlatHistSimulator.h"

LocalFlatHistSimulator::LocalFlatHistSimulator(
	Base_decoder * decoderPtr, double epsilon, int l, int kMin, int alpha, int beta,
	int unconWithoutAB, int unconWithAB, int conWithoutAB, int conWithAB)
	: BaseSimulator(0, 0, decoderPtr) {
	_epsilon = epsilon;
	_l = l;
	_kMin = kMin;
	_alpha = alpha;
	_beta = beta;
	_unconWithoutAB = unconWithoutAB;
	_unconWithAB = unconWithAB;
	_conWithoutAB = conWithoutAB;
	_conWithAB = conWithAB;
}

void LocalFlatHistSimulator::Run(std::vector<double> snrArray,
	std::vector<double>& ebn0Array,
	std::vector<double>& ferArray,
	std::vector<double>& sigmaArray,
	std::vector<int>& testsCountArray,
	std::vector<std::chrono::milliseconds>& elapsedTimeArray)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::vector<int> codeword(_n, 0);
	std::vector<double> llrs(_n, 0);
	std::vector<int> decoded(_n, 0);
	std::vector<double> F(_l);
	for (size_t i = 0; i < _l; i++) {
		if (i < _kMin) F[i] = _alpha * pow(_l - i, _beta);
		else F[i] = 1;
	}
	bool isFailed = false;
	double vMin = 0, vMax = 2;

	for (size_t ii = 0; ii < snrArray.size(); ii++) {
		double sigma = GetSigma(snrArray[ii]);

		std::vector<int> H(_l, 0), G(_l, 0), E(_l, 0);
		std::vector<double> prob(_l, log(1.0 / _l)), probNew(_l), probPrev(_l),
		                    probCond(_l, log(1.0 / _l)), probCondNew(_l);
		std::vector<double> z, newZ, g(_l);
		std::normal_distribution<double> normDistr(0, sigma);
		std::uniform_real_distribution<double> unDistr(0.0, 1.0);
		double newLoss, newBit, dz, probBitAcceptance, probStateAcceptance, probOmega, gSum, 
		       gHat, probSumAB, normCoef;
		int curBin, newBin;
		bool isFlat = false, isError;

		z = findStartPoint();
		newLoss = lossFunc(z, codeword);
		curBin = (int) floor((newLoss - vMin) / (vMax - vMin) * (_l - 1));
		curBin = std::min(curBin, _l - 1);
		curBin = std::max(curBin, 0);
		H[curBin]++;

		while (!isFlat) {
			probPrev = prob;

			for (size_t i = 0; i < _unconWithoutAB; i++) {
				for (size_t i = 0; i < _n; i++) {
					llrs[i] = -2 * (2 * codeword[i] - 1 + z[i]) / (sigma * sigma);
				}
				decoded = _decoderPtr->Decode(llrs, &isFailed);
				if (decoded != codeword) {
					G[curBin] += 1;
					std::cout << "error!\n";
				}

				newZ = z;
				for (size_t i = 0; i < _n; i++) {
					dz = normDistr(gen);
					newBit = z[i] + _epsilon * dz;
					probBitAcceptance = std::min(normalPDF(newBit, 0, sigma) /
						normalPDF(z[i], 0, sigma),
						1.0);
					if (unDistr(gen) <= probBitAcceptance)
						newZ[i] = newBit;
				}

				newLoss = lossFunc(newZ, codeword);
				std::cout << "loss " << newLoss << "\n";
				newBin = (int) floor((newLoss - vMin) / (vMax - vMin) * (_l - 1));
				newBin = std::min(newBin, _l - 1);
				newBin = std::max(newBin, 0);
				
				probStateAcceptance = std::min(prob[curBin] - prob[newBin], 0.0);
				if (log(unDistr(gen)) <= probStateAcceptance) {
					z = newZ;
					curBin = newBin;
				}
				H[curBin]++;
			}

			for (size_t i = 0; i < _unconWithAB; i++) {
				for (size_t i = 0; i < _n; i++) {
					llrs[i] = -2 * (2 * codeword[i] - 1 + z[i]) / (sigma * sigma);
				}
				decoded = _decoderPtr->Decode(llrs, &isFailed);
				if (decoded != codeword) {
					G[curBin] += 1;
				}

				newZ = z;
				for (size_t i = 0; i < _n; i++) {
					dz = normDistr(gen);
					newBit = z[i] + _epsilon * dz;
					probBitAcceptance = std::min(normalPDF(newBit, 0, sigma) /
						normalPDF(z[i], 0, sigma),
						1.0);
					if (unDistr(gen) <= probBitAcceptance)
						newZ[i] = newBit;
				}

				newLoss = lossFunc(newZ, codeword);
				newBin = (int) floor((newLoss - vMin) / (vMax - vMin) * (_l - 1));
				newBin = std::min(newBin, _l - 1);
				newBin = std::max(newBin, 0);
				probStateAcceptance = std::min(prob[curBin] * F[curBin] - prob[newBin] * F[newBin], 0.0);
				if (log(unDistr(gen)) <= probStateAcceptance) {
					z = newZ;
					curBin = newBin;
				}
				H[curBin]++;
			}

			probOmega = 1;
			for (size_t i = 0; i < _kMin; i++) {
				probOmega -= std::exp(prob[i]);
			}

			gSum = 0;
			for (size_t i = 1; i < _l; i++) {
				if (H[i] + H[i-1] == 0) g[i] = 0;
				else g[i] = (double) H[i] * H[i-1] / (H[i] + H[i-1]);
				gSum += g[i];
			}

			probNew[0] = 1;
			for (size_t i = 1; i < _l; i++) {
				if (g[i] == 0) gHat = 0;
				else gHat = g[i] / gSum;
				if (gHat == 0) probNew[i] = probNew[i-1] + prob[i] - prob[i-1];
				else probNew[i] = probNew[i-1] + prob[i] - prob[i-1] + gHat * log((double)H[i] / H[i-1]);
			}

			probSumAB = 0;
			for (size_t i = _kMin; i < _l; i++) {
				probSumAB += probOmega * std::exp(probNew[i]);
			}

			normCoef = 0;
			for (size_t i = 0; i < _l; i++) {
				probNew[i] -= probSumAB;
				normCoef += std::exp(probNew[i]);
			}
			for (size_t i = 0; i < _l; i++) {
				probNew[i] -= log(normCoef);
			}
			std::cout << "Prob:\n";
            for (int i = 0; i < _l; i++) {
				std::cout << std::exp(probNew[i]) << " ";
			}
			std::cout << "\n";

			isFlat = true;
			for (size_t i = 0; i < _l; i++) {
				if (fabs((std::exp(probNew[i]) - std::exp(probPrev[i])) / std::exp(probNew[i])) >= 0.1) {
					isFlat = false;
					break;
				}
			}
			prob = probNew;

			z = findStartPoint();
			newLoss = lossFunc(z, codeword);
			curBin = (int) floor((newLoss - vMin) / (vMax - vMin) * (_l - 1));
			curBin = std::min(curBin, _l - 1);
			curBin = std::max(curBin, 0);
			E[curBin]++;

			for (size_t i = 0; i < _conWithoutAB; i++) {
				for (size_t i = 0; i < _n; i++) {
					llrs[i] = -2 * (2 * codeword[i] - 1 + z[i]) / (sigma * sigma);
				}
				isError = false;
				decoded = _decoderPtr->Decode(llrs, &isFailed);
				if (decoded != codeword) {
					isError = true;
				}

				newZ = z;
				for (size_t i = 0; i < _n; i++) {
					dz = normDistr(gen);
					newBit = z[i] + _epsilon * dz;
					probBitAcceptance = std::min(normalPDF(newBit, 0, sigma) /
						normalPDF(z[i], 0, sigma),
						1.0);
					if (unDistr(gen) <= probBitAcceptance)
						newZ[i] = newBit;
				}

				newLoss = lossFunc(newZ, codeword);
				newBin = (int) floor((newLoss - vMin) / (vMax - vMin) * (_l - 1));
				newBin = std::min(newBin, _l - 1);
				newBin = std::max(newBin, 0);
				
				if (!isError) probBitAcceptance = 1.0; // illegal value (ln(P) <= 0)
				else probStateAcceptance = std::min(probCond[curBin] - probCond[newBin], 0.0);
				if (probBitAcceptance != 1.0 && log(unDistr(gen)) <= probStateAcceptance) {
					z = newZ;
					curBin = newBin;
				}
				E[curBin]++;
			}

			for (size_t i = 0; i < _conWithAB; i++) {
				for (size_t i = 0; i < _n; i++) {
					llrs[i] = -2 * (2 * codeword[i] - 1 + z[i]) / (sigma * sigma);
				}
				isError = false;
				decoded = _decoderPtr->Decode(llrs, &isFailed);
				if (decoded != codeword) {
					isError = true;
				}

				newZ = z;
				for (size_t i = 0; i < _n; i++) {
					dz = normDistr(gen);
					newBit = z[i] + _epsilon * dz;
					probBitAcceptance = std::min(normalPDF(newBit, 0, sigma) /
						normalPDF(z[i], 0, sigma),
						1.0);
					if (unDistr(gen) <= probBitAcceptance)
						newZ[i] = newBit;
				}

				newLoss = lossFunc(newZ, codeword);
				newBin = (int) floor((newLoss - vMin) / (vMax - vMin) * (_l - 1));
				newBin = std::min(newBin, _l - 1);
				newBin = std::max(newBin, 0);
				
				if (!isError) probBitAcceptance = 1.0; // illegal value (ln(P) <= 0)
				else probStateAcceptance = std::min(probCond[curBin] * F[curBin] - probCond[newBin] * F[newBin], 0.0);
				if (probBitAcceptance != 1.0 && log(unDistr(gen)) <= probStateAcceptance) {
					z = newZ;
					curBin = newBin;
				}
				E[curBin]++;
			}

			probOmega = 1;
			for (size_t i = 0; i < _kMin; i++) {
				probOmega -= std::exp(probCond[i]);
			}

			gSum = 0;
			for (size_t i = 1; i < _l; i++) {
				if (E[i] + E[i-1] == 0) g[i] = 0;
				else g[i] = (double) E[i] * E[i-1] / (E[i] + E[i-1]);
				gSum += g[i];
			}

			probCondNew[0] = 1;
			for (size_t i = 1; i < _l; i++) {
				if (g[i] == 0) gHat = 0;
				else gHat = g[i] / gSum;
				if (gHat == 0) probCondNew[i] = probCondNew[i-1] + probCond[i] - probCond[i-1];
				else probCondNew[i] = probCondNew[i-1] + probCond[i] - probCond[i-1] + gHat * log((double)E[i] / E[i-1]);
			}

			probSumAB = 0;
			for (size_t i = _kMin; i < _l; i++) {
				probSumAB += probOmega * std::exp(probCondNew[i]);
			}
			for (size_t i = 0; i < _l; i++) {
				probCondNew[i] -= probSumAB;
			}

			normCoef = 0;
			for (size_t i = 0; i < _l; i++) {
				probCondNew[i] -= probSumAB;
				normCoef += std::exp(probCondNew[i]);
			}
			for (size_t i = 0; i <  _l; i++) {
				probCondNew[i] -= log(normCoef);
			}

			probCond = probCondNew;
		}

		int maxInd = -1, maxEr = -1;
		for (size_t i = _l - _kMin; i < _l; i++) {
			if (maxEr < H[i] && E[i] > 0) {
				maxInd = i;
				maxEr = H[i];
			}
		}

		std::cout << "H:\n";
		for(int i = 0; i < _l; i++) {
			std::cout << H[i] << " ";
		}
		std::cout << "\n";
		std::cout << "G:\n";
		for(int i = 0; i < _l; i++) {
			std::cout << G[i] << " ";
		}
		std::cout << "\n";
		std::cout << "E:\n";
		for(int i = 0; i < _l; i++) {
			std::cout << E[i] << " ";
		}
		std::cout << "\n";
		std::cout << maxInd << " " << G[maxInd] << " " << H[maxInd] << " " << 
		             std::exp(prob[maxInd]) << " " << std::exp(probCond[maxInd]) << "\n";

		double PErr = (double) G[maxInd] / H[maxInd] * std::exp(prob[maxInd]) / std::exp(probCond[maxInd]);
		ferArray[ii] = PErr;
	}
}


double LocalFlatHistSimulator::normalPDF(double x, double m, double s) {
	static const double inv_sqrt_2pi = 0.3989422804014327;
	double a = (x - m) / s;
	return inv_sqrt_2pi / s * std::exp(-0.5 * a * a);
}


double LocalFlatHistSimulator::lossFunc(const std::vector<double>& z, const std::vector<int>& codeword) {
	double lossSum = 0;
	for (size_t l = 0; l < _n; l++) {
		int q = pow(-1, codeword[l]);
		lossSum += pow(fabs(((q * z[l]) > 0) * z[l]), 3);
	}
	return pow(lossSum / _n, 1.0 / 3);
}

std::vector<double> LocalFlatHistSimulator::findStartPoint() {
	// TODO: bsearch
	std::vector<double> startPoint(_n, 2.001);
	startPoint[0] = 2.001;
	startPoint[1] = 2.001;
	startPoint[2] = 2.001;
	//startPoint[3] = 2.001;
	return startPoint;
}
