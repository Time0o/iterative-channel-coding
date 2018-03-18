#include <vector>

#include "ctrlmat.h"
#include "blockdecoder.h"

int main()
{
    // control matrices
    std::vector<std::vector<int>> H_bf_vals;
    H_bf_vals.push_back({1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1});
    H_bf_vals.push_back({1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0});
    H_bf_vals.push_back({1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0});
    H_bf_vals.push_back({0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0});
    H_bf_vals.push_back({0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1});
    H_bf_vals.push_back({0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1});
    CtrlMat H_bf(12, 6, 0, H_bf_vals);

    CtrlMat H_mlg(15, 7, 5, "11010001", true);

    std::vector<std::vector<int>> H_min_sum_vals;
    H_min_sum_vals.push_back({1, 1, 1, 1, 0, 0, 0});
    H_min_sum_vals.push_back({1, 1, 0, 0, 1, 1, 0});
    H_min_sum_vals.push_back({1, 0, 1, 0, 1, 0, 1});
    CtrlMat H_min_sum(7, 4, 3, H_min_sum_vals);

    // input vectors
    std::vector<double> b_bf {
        -0.9, -0.9, 1.0, 0.1, 0.3, 1.0, 0.8, -0.7, 1.0, 0.9, -1.0, 0.4
    };

    std::vector<double> b_mlg {
        -1.0, -1.0, -1.0, 1.0, -1.0, 1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0,
         1.0, 1.0
    };

    std::vector<double> b_min_sum {0.4, 0.2, -1.0, 0.3, 0.8, -0.8, 0.1};

    // BF
    BF bf(H_bf);
    bf.decode(b_bf);
    std::cout << '\n';

    // WBF
    WBF wbf(H_bf);
    wbf.decode(b_bf);
    std::cout << '\n';

    // MWBF
    MWBF mwbf(H_bf);
    mwbf.set_alpha(0.7);
    mwbf.decode(b_bf);
    std::cout << '\n';

    // One Step MLG
    OneStepMLG one_step_mlg(H_mlg);
    one_step_mlg.decode(b_mlg);
    std::cout << '\n';

    // Min Sum
    MinSum min_sum(H_min_sum);
    min_sum.decode(b_min_sum);
}
