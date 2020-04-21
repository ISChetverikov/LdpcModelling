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

SimulationIterationResults LocalFlatHistSimulator::Run(double snr) {
    SimulationIterationResults result;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::vector<int> codeword(_n, 0);
    std::vector<double> llrs(_n, 0), startPoint, startPointCon;
    std::vector<int> decoded(_n, 0);
    std::vector<double> F(_l);
    for (size_t i = 0; i < _l; i++) {
        if (i < _kMin) F[i] = log(_alpha * pow(_kMin - i, _beta));
        else F[i] = 0;
    }
    startPoint = findStartPoint(codeword, 0);
    startPointCon = findStartPoint(codeword, 0);

    auto t1 = std::chrono::steady_clock::now();
    double sigma = GetSigma(snr);
    int iterationsCount = 0, startPointBin, numCycles;
    std::pair<double, double> V = findRangeV(startPoint, codeword, sigma);
    double vMin = V.first, vMax = 1.02 * V.second;
    std::cout << "Vmin: " << vMin << ", Vmax: " << vMax << ", SNR: " << snr << "\n";
    std::vector<int> H(_l, 0), G(_l, 0), E(_l, 0), curH(_l, 0), curE(_l, 0);
    std::vector<double> prob(_l, log(1.0 / _l)), probNew(_l), cStartPoint(_l),
                        probCond(_l, log(1.0 / _l)), probCondNew(_l), gSum(_l, 0);
    std::vector<std::string> realProb;
    std::vector<double> z, newZ, g(_l);
    std::normal_distribution<double> normDistr(0, sigma);
    std::uniform_real_distribution<double> unDistr(0.0, 1.0);
    double newLoss, newBit, dz, probBitAcceptance, probStateAcceptance, probOmega,
            gHat, probSumAB, normCoef;
    int curBin, newBin, numPoints = 0, numAccepted = 0;
    bool isFlat = false, isError, isFailed, inV, startToOne = false, isCStart = false;

    realProb = MCHist(codeword, sigma, vMin, vMax, 400000);
    z = startPoint;
    newLoss = lossFunc(z, codeword);
    curBin = (int) floor((newLoss - vMin) / (vMax - vMin) * (_l - 1));
    startPointBin = curBin;
    if (curBin > _l - 1 || curBin < 0) inV = false;
    else inV = true;
    if (inV) curH[curBin]++;

    numCycles = 0;
    while (!isFlat) {
        for (size_t it = 0; it < _unconWithoutAB; it++) {
            if (inV) {
                for (size_t i = 0; i < _n; i++) {
                    llrs[i] = -2 * (2 * codeword[i] - 1 + z[i]) / (sigma * sigma);
                }
                iterationsCount++;
                decoded = _decoderPtr->Decode(llrs, &isFailed);
                if (decoded != codeword) {
                    G[curBin] += 1;
                    if (curBin != startPointBin && !isCStart) {
                        isCStart = true;
                        cStartPoint = z;
                    }
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
                curH[curBin]++;
            }
        }

        probNew[0] = -1;
        for (size_t k = 0; k < _l - 1; k++) {
            if (curH[k] + curH[k+1] == 0) g[k] = 0;
            else g[k] = (double) curH[k] * curH[k+1] / (curH[k] + curH[k+1]);
            gSum[k] += g[k];
            if (g[k] == 0) gHat = 0;
            else gHat = g[k] / gSum[k];
            if (gHat == 0) probNew[k+1] = probNew[k] + prob[k+1] - prob[k];
            else probNew[k+1] = probNew[k] + prob[k+1] - prob[k] + gHat * log((double) curH[k+1]/curH[k]);
        }

        normCoef = 0;
        for (int i = 0; i < _l; i++) {
            normCoef += std::exp(probNew[i]);
        }
        for (size_t i = 0; i < _l; i++) {
            probNew[i] -= log(normCoef);
        }

        isFlat = true;
        for (size_t i = 0; i < _kMin; i++) {
            if (log(fabs(std::exp(probNew[i]) - std::exp(prob[i]))) - probNew[i] >= log(0.1)) {
                isFlat = false;
                break;
            }
        }

        prob = probNew;

        for (int i = 0; i < _l; i++) {
            H[i] += curH[i];
            curH[i] = 0;
        }

        numCycles++;
    }

    std::cout << "[0, kMin) is flat!\nnumber of cycles: " << numCycles << "\n";

    probOmega = 1;
    for (size_t i = 0; i < _kMin; i++) {
        probOmega -= std::exp(prob[i]);
    }

    isFlat = false;
    for (int i = 0; i < gSum.size(); i++) {
        gSum[i] = 0;
    }

    numCycles = 0;
    while (!isFlat) {
        for (size_t it = 0; it < _unconWithAB; it++) {
            if (inV) {
                for (size_t i = 0; i < _n; i++) {
                    llrs[i] = -2 * (2 * codeword[i] - 1 + z[i]) / (sigma * sigma);
                }
                iterationsCount++;
                decoded = _decoderPtr->Decode(llrs, &isFailed);
                if (decoded != codeword) {
                    G[curBin] += 1;
                    if (curBin != startPointBin && !isCStart) {
                        isCStart = true;
                        cStartPoint = z;
                    }
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
            if (newBin > _l - 1 || newBin < 0) inV = false;
            else inV = true;

            if (inV) {
                probStateAcceptance = std::min(prob[curBin] + F[curBin] - prob[newBin] - F[newBin], 0.0);
                if (log(unDistr(gen)) <= probStateAcceptance) {
                    z = newZ;
                    curBin = newBin;
                }
                curH[curBin]++;
            }
        }

        probNew[0] = -1;
        for (size_t k = 0; k < _l - 1; k++) {
            if (curH[k] + curH[k+1] == 0) g[k] = 0;
            else g[k] = (double) curH[k] * curH[k+1] / (curH[k] + curH[k+1]);
            gSum[k] += g[k];
            if (g[k] == 0) gHat = 0;
            else gHat = g[k] / gSum[k];
            if (gHat == 0) probNew[k+1] = probNew[k] + prob[k+1] - prob[k];
            else probNew[k+1] = probNew[k] + prob[k+1] - prob[k] + gHat * log((double) curH[k+1]/curH[k]);
        }

        normCoef = 0;
        for (int i = 0; i < _l; i++) {
            normCoef += std::exp(probNew[i]);
        }
        for (size_t i = 0; i < _l; i++) {
            probNew[i] -= log(normCoef);
        }
        /*std::cout << "Prob before: \n";
        for (int i = 0; i < _l; i++) {
            std::cout << probNew[i] << " ";
        }
        std::cout << "\n";*/

        probSumAB = 0;

        for (size_t i = _kMin; i < _l; i++) {
            probSumAB += std::exp(probNew[i]);
        }

        for (size_t i = _kMin; i < _l; i++) {
            probNew[i] = probNew[i] - log(probSumAB) + log(probOmega);
        }

        for (int i = 0; i < _kMin; i++) {
            probNew[i] = probNew[i] - log(1 - probSumAB) + log(1 - probOmega);
        }

        normCoef = 0;
        for (int i = 0; i < _l; i++) {
            normCoef += std::exp(probNew[i]);
        }
        for (size_t i = 0; i < _l; i++) {
            probNew[i] -= log(normCoef);
        }

        /*std::cout << "Prob after:\n";
        for (int i = 0; i < _l; i++) {
            std::cout << prob[i] << " ";
        }
        std::cout << "\n";*/

        isFlat = true;
        for (size_t i = _kMin; i < _l; i++) {
            if (log(fabs(std::exp(probNew[i]) - std::exp(prob[i]))) - probNew[i] >= log(0.1)) {
                isFlat = false;
                break;
            }
        }

        if (isFlat) {
            int sumG = 0;
            for (int i = 0; i < _l; i++) {
                sumG += G[i];
            }
            if (sumG < 2) isFlat = false;
        }

        prob = probNew;

        /*std::cout << "CurH:\n";
        for (int i = 0; i < _l; i++) {
            std::cout << i << " " << curH[i] << " " << prob[i] << "\n";
        }*/

        for (int i = 0; i < _l; i++) {
            H[i] += curH[i];
            curH[i] = 0;
        }

        std::cout << "\nprob omega, prom omega est: " << probOmega << " " << probSumAB << "\n";
        numCycles++;
    }

    if (startPointBin >= 0 && startPointBin <= _l - 1) G[startPointBin]--;

    std::cout << "Real prob, observed prob:\n";
    for (int i = 0; i < _l; i++) {
        std::cout << i << " " << realProb[i] << " " << prob[i] << "\n";
    }
    std::cout << "number of cycles: " << numCycles << "\n";
    std::cout << "Acceptance rate: " << (double) numAccepted / numPoints << "\n";

    isFlat = false;
    for (int i = 0; i < gSum.size(); i++) {
        gSum[i] = 0;
    }

    z = cStartPoint;
    newLoss = lossFunc(z, codeword);
    curBin = (int) floor((newLoss - vMin) / (vMax - vMin) * (_l - 1));
    startPointBin = curBin;
    if (curBin > _l - 1 || curBin < 0) inV = false;
    else inV = true;
    if (inV) curE[curBin]++;

    numCycles = 0;
    while (!isFlat) {
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
                curE[curBin]++;
                if (curBin != startPointBin && !startToOne) {
                    startToOne = true;
                    curE[startPointBin] = 1;
                }
            }
        }

        probCondNew[0] = -1;
        for (size_t k = 0; k < _l - 1; k++) {
            if (curE[k] + curE[k+1] == 0) g[k] = 0;
            else g[k] = (double) curE[k] * curE[k+1] / (curE[k] + curE[k+1]);
            gSum[k] += g[k];
            if (g[k] == 0) gHat = 0;
            else gHat = g[k] / gSum[k];
            if (gHat == 0) probCondNew[k+1] = probCondNew[k] + probCond[k+1] - probCond[k];
            else probCondNew[k+1] = probCondNew[k] + probCond[k+1] - probCond[k] + gHat * log((double) curE[k+1] / curE[k]);
        }

        normCoef = 0;
        for (int i = 0; i < _l; i++) {
            normCoef += std::exp(probCondNew[i]);
        }
        for (size_t i = 0; i < _l; i++) {
            probCondNew[i] -= log(normCoef);
        }

        isFlat = true;
        for (size_t i = 0; i < _kMin; i++) {
            if (log(fabs(std::exp(probCondNew[i]) - std::exp(probCond[i]))) - probCondNew[i] >= log(0.1)) {
                isFlat = false;
                break;
            }
        }

        probCond = probCondNew;

        /*for (int i = 0; i < _l; i++) {
            std::cout << i << " " << curE[i] << "\n";
        }
        for (int i = 0; i < _l; i++) {
            std::cout << i << " " << probCond[i] << "\n";
        }*/

        for (int i = 0; i < _l; i++) {
            E[i] += curE[i];
            curE[i] = 0;
        }
        numCycles++;
    }

    std::cout << "[0, kMin) is flat!\n number of cycles: " << numCycles << "\n";

    probOmega = 1;
    for (size_t i = 0; i < _kMin; i++) {
        probOmega -= std::exp(probCond[i]);
    }

    isFlat = false;
    for (int i = 0; i < gSum.size(); i++) {
        gSum[i] = 0;
    }

    numCycles = 0;
    while (!isFlat) {
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
                curE[curBin]++;
                if (curBin != startPointBin && !startToOne) {
                    startToOne = true;
                    curE[startPointBin] = 1;
                }
            }
        }

        isFlat = false;
        for (int i = 1; i < _l; i++) {
            if (E[i] + curE[i] > 0 && E[i - 1] + curE[i - 1] > 0) {
                isFlat = true;
                break;
            }
        }

        probCondNew[0] = -1;
        for (size_t k = 0; k < _l - 1; k++) {
            if (curE[k] + curE[k+1] == 0) g[k] = 0;
            else g[k] = (double) curE[k] * curE[k+1] / (curE[k] + curE[k+1]);
            gSum[k] += g[k];
            if (g[k] == 0) gHat = 0;
            else gHat = g[k] / gSum[k];
            if (gHat == 0) probCondNew[k+1] = probCondNew[k] + probCond[k+1] - probCond[k];
            else probCondNew[k+1] = probCondNew[k] + probCond[k+1] - probCond[k] + gHat * log((double) curE[k+1] / curE[k]);
        }

        normCoef = 0;
        for (int i = 0; i < _l; i++) {
            normCoef += std::exp(probCondNew[i]);
        }
        for (size_t i = 0; i < _l; i++) {
            probCondNew[i] -= log(normCoef);
        }

        probSumAB = 0;

        for (size_t i = _kMin; i < _l; i++) {
            probSumAB += std::exp(probCondNew[i]);
        }

        for (size_t i = _kMin; i < _l; i++) {
            probCondNew[i] = probCondNew[i] - log(probSumAB) + log(probOmega);
        }

        for (int i = 0; i < _kMin; i++) {
            probCondNew[i] = probCondNew[i] - log(1 - probSumAB) + log(1 - probOmega);
        }

        normCoef = 0;
        for (int i = 0; i < _l; i++) {
            normCoef += std::exp(probCondNew[i]);
        }
        for (size_t i = 0; i < _l; i++) {
            probCondNew[i] -= log(normCoef);
        }


        isFlat = true;
        for (size_t i = _kMin; i < _l; i++) {
            if (log(fabs(std::exp(probCondNew[i]) - std::exp(probCond[i]))) - probCondNew[i] >= log(0.1)) {
                isFlat = false;
                break;
            }
        }

        probCond = probCondNew;

        /*for (int i = 0; i < _l; i++) {
            std::cout << i << " " << curE[i] << "\n";
        }
        for (int i = 0; i < _l; i++) {
            std::cout << i << " " << probCond[i] << "\n";
        }*/

        for (int i = 0; i < _l; i++) {
            E[i] += curE[i];
            curE[i] = 0;
        }

        std::cout << "\nprob omega, prom omega est: " << probOmega << " " << probSumAB << "\n";
        numCycles++;
    }

    std::cout << "number of cycles: " << numCycles << "\n";

    double probPartSum = 0;
    int cntNonZeroBins = 30;
    for (size_t i = 310; i < 310 + 30; i++) {
        if (H[i] == 0) {
            std::cout << "0\n";
            continue;
        }
        probPartSum += (double) G[i] / H[i] * std::exp(prob[i] - probCond[i]);
        std::cout << (double) G[i] / H[i] * std::exp(prob[i] - probCond[i]) << "\n";
    }
    /*std::cout << "Prob:\n";
    for (int i = 0; i < _l; i++){
        std::cout << prob[i] << " ";
    }
    std::cout << "\nProbCond:\n";
    for (int i = 0; i < _l; i++) {
        std::cout << probCond[i] << " ";
    }
    std::cout << "\n";*/
    std::cout << "bin, H, G, prob, probCond\n";
    for (int i = 0; i < _l; i++) {
        std::cout << i << " " << H[i] << " " << G[i] << " " << prob[i] << " " << probCond[i] << "\n";
    }

    double PErr = probPartSum / cntNonZeroBins;
    auto t2 = std::chrono::steady_clock::now();
    result.snr = snr;
    result.ebn0 = GetEbN0(snr, _m, _n);
    result.sigma = sigma;

    result.fer = std::min(1.0, PErr);

    result.elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    result.testsCount = iterationsCount;
    result.rejectionsCount = 0; // idk where

    std::cout << PErr << " " << iterationsCount << "\n";
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


std::vector<double> LocalFlatHistSimulator::findEvristicStartPoint(const std::vector<int>& codeword, double sigma, double minV, double maxV, int& iterationsCount) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> normDistr(0, sigma);
    std::uniform_real_distribution<double> unDistr(0.0, 1.0);
    std::vector<double> z(_n, 0), newZ, llrs(_n);
    std::vector<int> decoded(_n);
    double newLoss, newBit, probBitAcceptance, dz;
    bool isFailed;
    newLoss = lossFunc(z, codeword);
    while (true) {
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
        for (size_t i = 0; i < _n; i++) {
            llrs[i] = -2 * (2 * codeword[i] - 1 + newZ[i]) / (sigma * sigma);
        }
        iterationsCount++;
        decoded = _decoderPtr->Decode(llrs, &isFailed);
        if (decoded != codeword && newLoss >= minV && newLoss <= maxV) {
            return newZ;
        }
        z = newZ;
    }
}

std::vector<std::string> LocalFlatHistSimulator::MCHist(const std::vector<int>& codeword, double sigma, double minV, double maxV, int iterationsCount) {
    std::vector<std::string> res;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> normDistr(0, sigma);
    std::vector<double> z(_n, 0), newZ, llrs(_n);
    std::vector<int> decoded(_n), H(_l, 0);
    double newLoss, newBit, probBitAcceptance, dz;
    int curBin;
    for (int it = 0; it < iterationsCount; it++) {
        for (size_t i = 0; i < _n; i++) {
            llrs[i] = z[i] + normDistr(rd);
        }
        newLoss = lossFunc(llrs, codeword);
        if (newLoss >= minV && newLoss <= maxV) {
            curBin = (int) floor((newLoss - minV) / (maxV - minV) * (_l - 1));
            H[curBin]++;
        }
    }
    int sumH = 0;
    for (int i = 0; i < _l; i++) {
        sumH += H[i];
    }
    for (int i = 0; i < _l; i++) {
        if (!H[i]) res.push_back("inf");
        else res.push_back(std::to_string(log((double) H[i] / sumH)));
    }
    return res;
}