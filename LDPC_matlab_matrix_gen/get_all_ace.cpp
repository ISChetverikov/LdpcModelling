#include <mex.h>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <algorithm> // for min
#include <functional>
#ifdef PROFILE
#include <callgrind.h>
#endif

using namespace std;

const double *col_weight, *offsets, *nb_values;
vector<int> H_exp, H_q;
size_t m, max_ace, tests;
double* current_ace;
int factor, q, base_n;

vector<int> blocked;
size_t vertex;
vector<vector<size_t>> A;
vector<vector<size_t>> Bidx;
vector<vector<int>> Bmask;
vector<size_t> stack;
vector<int> gcd_table;
vector<int> gcd_offsets;
size_t col;
size_t max_depth;

extern void list_paths(int u, int t, int md, vector<size_t>& path, vector<int> &blocked, const vector<vector<size_t>> &G, const function<void(const vector<size_t>&)>& output);
extern void list_paths(int s, int t, int md, const vector<vector<size_t>> &G, const function<void(const vector<size_t>&)>& output);

inline size_t gcd(size_t a, size_t b) {
    while (b!=0) {
        size_t c = a % b;
        a = b;
        b = c;
    }
    return a;
}

void calc_ace(const vector<size_t>& cycle_p);

bool circuit(size_t v, size_t col, int max_depth, vector<size_t>& stack, vector<int>& blocked, const vector<vector<size_t>>& A, vector<vector<size_t>>& Bidx, vector<vector<int>> &Bmask, const function<void(const vector<size_t>&)>& calc_ace);

template<typename T>
void fillA(const T* H, size_t m, size_t n) {
    for (size_t i=0; i<m; i++)
        for (size_t j=0; j<n; j++)
            if (H[i + j*m]) {
                A[j].push_back(i+n);
                A[i+n].push_back(j);
            }
}

// Calculated ace for the specified cycle and all offsets/nb_values
void calc_ace(const vector<size_t>& cycle_p)
{
    static int counter = 0;
    if (++counter % 100000 == 0)
        if (mexEvalString("drawnow"))
            throw runtime_error("interrupted");
    size_t cycle_len = cycle_p.size()-1;
    static vector<size_t> cycle;
    cycle.resize(cycle_len+1);
    for (size_t i=0; i<=cycle_len; i++)
        cycle[i] = cycle_p[i];
    for (size_t i=1; i<=cycle_len; i+=2)
        cycle[i] = cycle[i]-base_n;

    H_exp[col*m + m-1] = 0;
    H_q[col*m + m-1] = 0;
    int log_det_off = 0, log_det_nb = 0;
    int cycle_delta = 0;
    for (size_t i=0; i<cycle_len-1; i+=2) {
        if (cycle[ i ] == col && cycle[i+1] == m-1) cycle_delta++;
        if (cycle[i+2] == col && cycle[i+1] == m-1) cycle_delta--;
        log_det_off +=
            H_exp[cycle[ i ]*m + cycle[i+1]] -
            H_exp[cycle[i+2]*m + cycle[i+1]];
        log_det_nb +=
            H_q[cycle[ i ]*m + cycle[i+1]] -
            H_q[cycle[i+2]*m + cycle[i+1]];
    }

    if (cycle_delta < 0) {
        log_det_off = -log_det_off;
        log_det_nb = -log_det_nb;
        cycle_delta = -cycle_delta;
    }

    int weight = 0;
    for (size_t l=0; l<cycle_len; l+=2)
        weight += (int)col_weight[cycle[l]];
    //assert(cycle_len >= 4);

    int max_mult = min(factor, (int)floor(max_depth/cycle_len));
    if (cycle_delta == 0) {
        size_t multiplicity = gcd_table[(abs(log_det_off)%factor)*2*(q-1)+(abs(log_det_nb)%(q-1))];
        if (multiplicity <= max_mult) {
            // we have a cycle => calculate ACE
            int ace = (int)multiplicity*weight - (int)multiplicity*(int)cycle_len;
            for (size_t tt=0; tt<tests; tt++) {
                size_t idx = max_ace*tt + (size_t)multiplicity*cycle_len/2 - 2;
                // FIXME: can idx point outside of the matrix?
                if (current_ace[idx] < 0 || ace < current_ace[idx])
                    current_ace[idx] = ace;
            }
        }
    } else {
        mxAssert(cycle_delta == 1, "Cycle proved to be not simple");
        log_det_off = log_det_off % factor;
        if (log_det_off < 0) log_det_off += factor;
        log_det_nb = log_det_nb % (q-1);
        if (log_det_nb < 0) log_det_nb += q-1;
        int* pGCD = gcd_table.data() + log_det_off*2*(q-1) + log_det_nb;
        for (int tt=0; tt<tests; tt++) {
            int multiplicity = pGCD[gcd_offsets[tt]];
            if (multiplicity <= max_mult) {
                // we have a cycle => calculate ACE
                int ace = multiplicity*weight - multiplicity*(int)cycle_len;
                int idx = max_ace*tt + multiplicity*cycle_len/2 - 2;
                // FIXME: can idx point outside of the matrix?
                if (current_ace[idx] < 0 || ace < current_ace[idx])
                    current_ace[idx] = ace;
            }
        }
    }
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray * prhs[])
{
    if (nrhs != 11) mexErrMsgTxt("Usage: get_ace(H_exp, H_q, factor, max_depth, H_base, q, col_weight, col, offsets, nb_values, compat)");
    const double*H_exp_p = mxGetPr(prhs[0]);
    const double*H_q_p = mxGetPr(prhs[1]);
    size_t mat_size = mxGetNumberOfElements(prhs[0]);
    m = mxGetM(prhs[0]);
    H_exp = vector<int>(H_exp_p, H_exp_p + mat_size);
    H_q = vector<int>(H_q_p, H_q_p + mat_size);
    factor = (int)mxGetScalar(prhs[2]);
    max_depth = (int)mxGetScalar(prhs[3]);
    size_t base_m = mxGetM(prhs[4]);
    base_n = (int)mxGetN(prhs[4]);
    q = (int)mxGetScalar(prhs[5]);
    col_weight = mxGetPr(prhs[6]);
    col = (size_t)mxGetScalar(prhs[7])-1;
    tests = mxGetNumberOfElements(prhs[8]);
    offsets = mxGetPr(prhs[8]);
    nb_values = mxGetPr(prhs[9]);
    int compat = (int)mxGetScalar(prhs[10]);
    max_ace = max_depth/2-1;
    plhs[0] = mxCreateDoubleMatrix((int)max_ace, (int)tests, mxREAL);
    current_ace = mxGetPr(plhs[0]);
    for (size_t i=0; i<max_ace*tests; i++)
        current_ace[i] = -1;

    gcd_table.resize(4*factor*(q-1));
    for (size_t off=0; off<2*factor; off++)
        for (size_t nb=0; nb<2*(q-1); nb++) {
            size_t m1 = factor / gcd(off%factor, factor);
            size_t m2 = (q-1) / gcd(nb%factor, q-1);
            gcd_table[off*2*(q-1)+nb] = m1*m2/gcd(m1,m2);
        }
    gcd_offsets.resize(tests);
    for (int tt=0; tt<tests; tt++)
        gcd_offsets[tt] = (int)offsets[tt]*2*(q-1) + (int)nb_values[tt];

    vertex = base_m + base_n;
    A.assign(vertex, vector<size_t>());
    if (mxIsSparse(prhs[4]))
        mexErrMsgIdAndTxt("get_all_ace:sparse_not_supported", "Sparse matices are not supported");
    if (mxIsLogical(prhs[4])) fillA(mxGetLogicals(prhs[4]), base_m, base_n);
    else if (mxIsDouble(prhs[4]) && !mxIsComplex(prhs[4])) fillA(mxGetPr(prhs[4]), base_m, base_n);
    else mexErrMsgIdAndTxt("getCycles:type_not_supported", "Only double and logical are supported for matrix H_base");
    blocked.assign(vertex, 0);
    Bmask.assign(vertex, vector<int>(vertex, 0));
    Bidx.assign(vertex, vector<size_t>());
    stack.clear();
    int s = col;
    if (compat != 1) {
        stack.push_back(col);
        s = base_n+m-1;
        int next_v = col;
        int Asz = A[s].size();
        for (int i=0; i<A[s].size(); i++) if (A[s][i] == next_v) {
            swap(A[s][i], A[s].back());
            A[s].pop_back();
            break;
        }
        if (A[s].size() != Asz - 1) {
            mexPrintf("%d-variable and %d-check node are not adjacent\n", col+1, s-base_n);
            return;
        }
        for (int i=0; i<A[next_v].size(); i++) if (A[next_v][i] == s) {
            swap(A[next_v][i], A[next_v].back());
            A[next_v].pop_back();
            break;
        }
    }
#ifdef PROFILE
    CALLGRIND_START_INSTRUMENTATION;
#endif
    if (compat < 0) {
        list_paths(s, col, max_depth, stack, blocked, A, calc_ace);
    } else circuit(s, col, max_depth, stack, blocked, A, Bidx, Bmask, calc_ace);
#ifdef PROFILE
    CALLGRIND_STOP_INSTRUMENTATION;
#endif
}
