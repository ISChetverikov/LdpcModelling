

#include "../include/SP_decoder.h"
#include "../include/Exceptions.h"
#include "../include/MathOperations.h"

#include <vector>
#include <map>
#include <cmath>
#include <string>
#include <iostream>
#include <algorithm>

using namespace std;

static float phi(float x)
{
    static const float lim = 31.;
    static const float inv_lim =
    log( (exp((double)lim) + 1)/(exp((double)lim) - 1) );
    
    if (x > lim)
    {
        return( 0 );
    }
    else if (x < inv_lim)
        return( lim );
    else
    {
        double t = exp( (double) x );
        return (float) log( (t+1)/(t-1) );
    }
}

// Constructor
SP_decoder::SP_decoder(vector<vector<int>> H, int iterationsCount) : Base_decoder(H, iterationsCount) {
    
}

void SP_decoder::HorizontalStep(vector<map<int, int>> alpha, vector<map<int, double>> beta, vector<map<int, double>> &gamma) {
    int sign = 0;
    for (size_t j = 0; j < _m; j++)
    {
        int sum_sign = 1; // May be with count of sign it will be faster?
        double values = 0;
        for (auto &i : _checks[j])
        {
            sum_sign ^= alpha[j][i];
            values += phi(beta[j][i]);
            // TODO: Here we can get rid of redundant cycle ?
        }
        for (auto &i : _checks[j])
        {
            sign = sum_sign ^ alpha[j][i];
            gamma[j][i] = sign * phi(values - phi(beta[j][i]));
        }
    }
}

vector<int> SP_decoder::Decode(vector<double> llr, bool *isFailed) {
    
    size_t n = llr.size();
    // if (n != _n)
    //     throw IncorrectCodewordException("The codeword is not from a code with given check matrix");
    
    vector<int> result(n, 0);
    
    //vector<int> alpha0(n, 0);
    //vector<double> beta0(n, 0);
    vector<double> bits_values(n, 0);
    //for (size_t i = 0; i < n; i++)
    //{
    //    alpha0[i] = sign(llr[i]);
    //    beta0[i] = fabs(llr[i]);
    //}
    
    vector<map<int, int>> alpha(_m, map<int, int>());
    vector<map<int, double>> beta(_m, map<int, double>());
    vector<map<int, double>> gamma(_m, map<int, double>());
    
    // Init
    for (size_t j = 0; j < _m; j++)
    {
        for (auto &i : _checks[j])
        {
            alpha[j][i] = sign(llr[i]);
            beta[j][i] = fabs(llr[i]);
        }
    }
    
    size_t iteration = 0;
    
    while (true)
    {
        iteration++;
        HorizontalStep(alpha, beta, gamma);
        
        // Result of iteration
        for (size_t i = 0; i < n; i++)
        {
            bits_values[i] = llr[i];
            
            for (auto &j : _bits[i])
            {
                bits_values[i] += gamma[j][i];
            }
            
            //alpha0[i] = sign(bits_values[i]);
            //beta0[i] = fabs(bits_values[i]);
        }
        
        for (size_t i = 0; i < n; i++)
        {
            result[i] = (bits_values[i] <= 0) ? 1 : 0;
        }
        
        //if (_H * result == vector<int>(_m, 0))
        //    break;
        for (size_t j = 0; j < _m; j++)
        {
            int sum = 0;
            for (auto &i : _checks[j])
            {
                sum ^= result[i];
            }
            if ((bool)sum) {
                *isFailed = true;
                break;
            }
        }
        
        
        if (iteration >= _iterationsCount)
            break;
        
        // Vertical Step
        for (size_t i = 0; i < n; i++)
        {
            double value = bits_values[i];
            
            for (auto &j : _bits[i])
            {
                double new_value = value - gamma[j][i];
                alpha[j][i] = sign(new_value);
                beta[j][i] = fabs(new_value);
            }
        }
    }
    
    return result;
}

