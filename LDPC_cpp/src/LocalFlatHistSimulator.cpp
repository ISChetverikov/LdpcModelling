#include <random>
#include <chrono>
#include <iostream>
#include "../include/LocalFlatHistSimulator.h"

LocalFlatHistSimulator::LocalFlatHistSimulator(
    Base_decoder * decoderPtr, double epsilon, int l, int kMin, int alpha, int beta,
    int unconWithoutAB, int unconWithAB, int conWithoutAB, int conWithAB, int numIterForFindingV)
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
    _numIterForFindingV = numIterForFindingV;
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
    std::vector<double> llrs(_n, 0), startPoint, startPointCon;
    std::vector<int> decoded(_n, 0);
    std::vector<double> F(_l);
    for (size_t i = 0; i < _l; i++) {
        if (i < _kMin) F[i] = log(_alpha * pow(_l - i, _beta));
        else F[i] = 0;
    }
    startPoint = findStartPoint(codeword, 0);
    startPointCon = findStartPoint(codeword, 0);

    for (size_t ii = 0; ii < snrArray.size(); ii++) {
        auto t1 = std::chrono::steady_clock::now();
        double sigma = GetSigma(snrArray[ii]);
        int iterationsCount = 0;
        std::pair<double, double> V = findRangeV(startPoint, codeword, sigma);
        double vMin = 0.98 * V.first, vMax = 1.02 * V.second;
        std::cout << "Vmin: " << vMin << ", Vmax: " << vMax << "\n";
        std::vector<int> H(_l, 0), G(_l, 0), E(_l, 0);
        std::vector<double> prob(_l, log(1.0 / _l)), probNew(_l), 
                            probCond(_l, log(1.0 / _l)), probCondNew(_l), gSum(_l, 0);
        std::vector<double> z, newZ, g(_l);
        std::normal_distribution<double> normDistr(0, sigma);
        std::uniform_real_distribution<double> unDistr(0.0, 1.0);
        double newLoss, newBit, dz, probBitAcceptance, probStateAcceptance, probOmega,
               gHat, probSumAB, normCoef;
        int curBin, newBin, numPoints = 0, numAccepted = 0;
        bool isFlat = false, isError, isFailed, inV;

        z = startPoint;
        newLoss = lossFunc(z, codeword);
        curBin = (int) floor((newLoss - vMin) / (vMax - vMin) * (_l - 1));
        //curBin = std::min(curBin, _l - 1);
        //curBin = std::max(curBin, 0);
        if (curBin > _l - 1 || curBin < 0) inV = false;
        else inV = true;
        if (inV) H[curBin]++;

        while (!isFlat) {
            std::cout << "<<<=====NewIteration=====>>>\nSNR:" << snrArray[ii] << "\n";
            for (size_t it = 0; it < _unconWithoutAB; it++) {
                if (inV) {
                    for (size_t i = 0; i < _n; i++) {
                        llrs[i] = -2 * (2 * codeword[i] - 1 + z[i]) / (sigma * sigma);
                    }
                    iterationsCount++;
                    decoded = _decoderPtr->Decode(llrs, &isFailed);
                    if (decoded != codeword) {
                        G[curBin] += 1;
                    }
                }

                newZ = z;
                for (size_t i = 0; i < _n; i++) {
                    dz = normDistr(gen);
                    newBit = z[i] + _epsilon * dz;
                    probBitAcceptance = std::min(normalPDF(newBit, 0, sigma) /
                        normalPDF(z[i], 0, sigma),
                        1.0);
                    if (unDistr(gen) <= probBitAcceptance) {
                        newZ[i] = newBit;
                    }
                }

                newLoss = lossFunc(newZ, codeword);
                newBin = (int) floor((newLoss - vMin) / (vMax - vMin) * (_l - 1));
                //newBin = std::min(newBin, _l - 1);
                //newBin = std::max(newBin, 0);
                if (newBin > _l - 1 || newBin < 0) inV = false;
                else inV = true;
                
                if (inV) {
                    numPoints++;
                    probStateAcceptance = std::min(prob[curBin] - prob[newBin], 0.0);
                    if (log(unDistr(gen)) <= probStateAcceptance) {
                        z = newZ;
                        curBin = newBin;
                        numAccepted++;
                    }
                    H[curBin]++;
                }
            }

            for (size_t it = 0; it < _unconWithAB; it++) {
                if (inV) {
                    for (size_t i = 0; i < _n; i++) {
                        llrs[i] = -2 * (2 * codeword[i] - 1 + z[i]) / (sigma * sigma);
                    }
                    iterationsCount++;
                    decoded = _decoderPtr->Decode(llrs, &isFailed);
                    if (decoded != codeword) {
                        G[curBin] += 1;
                    }
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
                //newBin = std::min(newBin, _l - 1);
                //newBin = std::max(newBin, 0);
                if (newBin > _l - 1 || newBin < 0) inV = false;
                else inV = true;

                if (inV) {
                    probStateAcceptance = std::min(prob[curBin] + F[curBin] - prob[newBin] - F[newBin], 0.0);
                    if (log(unDistr(gen)) <= probStateAcceptance) {
                        z = newZ;
                        curBin = newBin;
                    }
                    H[curBin]++;
                }
            }

            probOmega = 1;
            for (size_t i = 0; i < _kMin; i++) {
                probOmega -= std::exp(prob[i]);
            }

            probNew[0] = -1;
            for (size_t k = 0; k < _l - 1; k++) {
                if (H[k] + H[k+1] == 0) g[k] = 0;
                else g[k] = (double) H[k] * H[k+1] / (H[k] + H[k+1]);
                gSum[k] += g[k];
                if (g[k] == 0) gHat = 0;
                else gHat = g[k] / gSum[k];
                if (gHat == 0) probNew[k+1] = probNew[k] + prob[k+1] - prob[k];
                else probNew[k+1] = probNew[k] + prob[k+1] - prob[k] + gHat * log((double) H[k+1]/H[k]);
            }

            probSumAB = 0;
            for (size_t i = _kMin; i < _l; i++) {
                probSumAB += probOmega * std::exp(probNew[i]);
            }

            normCoef = 0;
            for (size_t i = 0; i < _l; i++) {
                probNew[i] -= log(probSumAB);
                normCoef += std::exp(probNew[i]);
            }
            for (size_t i = 0; i < _l; i++) {
                probNew[i] -= log(normCoef);
            }

            isFlat = true;
            for (size_t i = 0; i < _l; i++) {
                if (log(fabs(std::exp(probNew[i]) - std::exp(prob[i]))) - probNew[i] >= log(0.1)) {
                    isFlat = false;
                    break;
                }
            }

            prob = probNew;

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

        }

        std::cout << "Acceptance rate: " << (double) numAccepted / numPoints << "\n";

        isFlat = false;
        for (int i = 0; i < gSum.size(); i++) {
            gSum[i] = 0;
        }

        z = startPointCon;
        newLoss = lossFunc(z, codeword);
        curBin = (int) floor((newLoss - vMin) / (vMax - vMin) * (_l - 1));
        //curBin = std::min(curBin, _l - 1);
        //curBin = std::max(curBin, 0);
        if (curBin > _l - 1 || curBin < 0) inV = false;
        else inV = true;
        if (inV) E[curBin]++;

        while (!isFlat) {
            std::cout << "<<<=====NewIteration=====>>>\nSNR:" << snrArray[ii] << "\n";
            for (size_t it = 0; it < _conWithoutAB; it++) {
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
                //newBin = std::min(newBin, _l - 1);
                //newBin = std::max(newBin, 0);
                if (newBin > _l - 1 || newBin < 0) inV = false;
                else inV = true;

                if (inV) {
                    for (size_t i = 0; i < _n; i++) {
                        llrs[i] = -2 * (2 * codeword[i] - 1 + newZ[i]) / (sigma * sigma);
                    }
                    isError = false;
                    iterationsCount++;
                    decoded = _decoderPtr->Decode(llrs, &isFailed);
                    if (decoded != codeword) {
                        isError = true;
                    }
                    
                    if (!isError) probStateAcceptance = 1.0; // illegal value (ln(P) <= 0)
                    else probStateAcceptance = std::min(probCond[curBin] - probCond[newBin], 0.0);
                    if (probStateAcceptance != 1.0 && log(unDistr(gen)) <= probStateAcceptance) {
                        z = newZ;
                        curBin = newBin;
                    }
                    E[curBin]++;
                }
            }

            for (size_t it = 0; it < _conWithAB; it++) {
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
                //newBin = std::min(newBin, _l - 1);
                //newBin = std::max(newBin, 0);

                if (newBin > _l - 1 || newBin < 0) inV = false;
                else inV = true;
                if (inV) {
                    for (size_t i = 0; i < _n; i++) {
                        llrs[i] = -2 * (2 * codeword[i] - 1 + newZ[i]) / (sigma * sigma);
                    }
                    isError = false;
                    iterationsCount++;
                    decoded = _decoderPtr->Decode(llrs, &isFailed);
                    if (decoded != codeword) {
                        isError = true;
                    }
                    
                    if (!isError) probStateAcceptance = 1.0; // illegal value (ln(P) <= 0)
                    else probStateAcceptance = std::min(probCond[curBin] + F[curBin] - probCond[newBin] - F[newBin], 0.0);
                    if (probStateAcceptance != 1.0 && log(unDistr(gen)) <= probStateAcceptance) {
                        z = newZ;
                        curBin = newBin;
                    }
                    E[curBin]++;
                }
            }

            probOmega = 1;
            for (size_t i = 0; i < _kMin; i++) {
                probOmega -= std::exp(probCond[i]);
            }

            probCondNew[0] = -1;
            for (size_t k = 0; k < _l - 1; k++) {
                if (E[k] + E[k+1] == 0) g[k] = 0;
                else g[k] = (double) E[k] * E[k+1] / (E[k] + E[k+1]);
                gSum[k] += g[k];
                if (g[k] == 0) gHat = 0;
                else gHat = g[k] / gSum[k];
                if (gHat == 0) probCondNew[k+1] = probCondNew[k] + probCond[k+1] - probCond[k];
                else probCondNew[k+1] = probCondNew[k] + probCond[k+1] - probCond[k] + gHat * log((double) E[k+1] / E[k]);
            }

            probSumAB = 0;
            for (size_t i = _kMin; i < _l; i++) {
                probSumAB += probOmega * std::exp(probCondNew[i]);
            }

            normCoef = 0;
            for (size_t i = 0; i < _l; i++) {
                probCondNew[i] -= log(probSumAB);
                normCoef += std::exp(probCondNew[i]);
            }

            for (size_t i = 0; i <  _l; i++) {
                probCondNew[i] -= log(normCoef);
            }


            isFlat = true;
            for (size_t i = 0; i < _l; i++) {
                if (log(fabs(std::exp(probCondNew[i]) - std::exp(probCond[i]))) - probCond[i] >= log(0.1)) {
                    isFlat = false;
                    break;
                }
            }

            probCond = probCondNew;

            std::cout << "E:\n";
            for(int i = 0; i < _l; i++) {
                std::cout << E[i] << " ";
            }
            std::cout << "\n";
        }

        double probPartSum = 0;
        int cntNonZeroBins = 0;
        for (size_t i = _kMin; i < _l; i++) {
            if (H[i] == 0) {
                continue;
            }
            probPartSum += (double) G[i] / H[i] * std::exp(prob[i]) / std::exp(probCond[i]);
            cntNonZeroBins++;
        }

        double PErr = probPartSum / cntNonZeroBins;
        auto t2 = std::chrono::steady_clock::now();
        sigmaArray[ii] = sigma;
        ebn0Array[ii] = GetEbN0(snrArray[ii], _m, _n);
        ferArray[ii] = PErr;
        elapsedTimeArray[ii] = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        testsCountArray[ii] = iterationsCount;
    }
}


double LocalFlatHistSimulator::normalPDF(double x, double m, double s) {
    static const double invSqrt2Pi = 0.3989422804014327;
    double a = (x - m) / s;
    return invSqrt2Pi / s * std::exp(-0.5 * a * a);
}


double LocalFlatHistSimulator::lossFunc(const std::vector<double>& z, const std::vector<int>& codeword) {
    double lossSum = 0;
    for (size_t l = 0; l < _n; l++) {
        int q = pow(-1, codeword[l]);
        lossSum += pow(fabs(((q * z[l]) > 0) * z[l]), 3);
    }
    return pow(lossSum / _n, 1.0 / 3);
}


std::vector<double> LocalFlatHistSimulator::findStartPoint(const std::vector<int>& codeword, int minN) {
    std::vector<double> startPoint(_n, 0), llrs(_n);
    std::vector<int> decoded;
    int l = minN - 1, r = _n + 1;
    bool isError, isFailed;
    int mPrev = -1, m;
    while (r - l > 1) {
        m = (l + r) / 2;

        // Invariant: [0, mPrev] inversed
        if (m > mPrev) {
            for (int i = mPrev + 1; i <= m; i++) {
                if (!codeword[i]) startPoint[i] = 1.001;
                else startPoint[i] = -1.001;
            }
        }
        else {
            for (size_t i = m + 1; i <= mPrev; i++) {
                startPoint[i] = 0;
            }
        }
        mPrev = m;

        for (size_t i = 0; i < _n; i++) {
            llrs[i] = -1 * (2 * codeword[i] - 1 + startPoint[i]);
        }
        isError = false;
        decoded = _decoderPtr->Decode(llrs, &isFailed);
        if (decoded != codeword) {
            isError = true;
        }

        if (!isError) {
            l = m;
        }
        else {
            r = m;
        }
    }
    if (!codeword[r]) startPoint[r] = 1.001;
    else startPoint[r] = -1.001;
    return startPoint;
}


std::pair<double, double> LocalFlatHistSimulator::findRangeV(const std::vector<double>& startPoint,
                                                             const std::vector<int>& codeword, double sigma) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> normDistr(0, sigma);
    std::uniform_real_distribution<double> unDistr(0.0, 1.0);
    std::vector<double> z, newZ;
    double minLoss = 10, maxLoss = 0, newLoss, newBit, probBitAcceptance, dz;
    z = startPoint;
    newLoss = lossFunc(z, codeword);
    minLoss = std::min(newLoss, minLoss);
    maxLoss = std::max(newLoss, maxLoss);
    for (int it = 0; it < _numIterForFindingV; it++) {
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
        minLoss = std::min(newLoss, minLoss);
        maxLoss = std::max(newLoss, maxLoss);
        z = newZ;
    }
    return std::make_pair(minLoss, maxLoss);
} 