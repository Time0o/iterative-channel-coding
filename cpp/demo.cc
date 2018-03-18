#include <vector>

#include "ctrlmat.h"
#include "blockdecoder.h"

int main()
{
    // BF
    std::vector<std::vector<int>> H_bf_vals;
    H_bf_vals.push_back({1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1});
    H_bf_vals.push_back({1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0});
    H_bf_vals.push_back({1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0});
    H_bf_vals.push_back({0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0});
    H_bf_vals.push_back({0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1});
    H_bf_vals.push_back({0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1});
    CtrlMat H_bf(12, 6, 0, H_bf_vals);

    std::vector<double> b_bf {-0.9, -0.9, 1.0, 0.1, 0.3, 1.0, 0.8, -0.7, 1.0,
                               0.9, -1.0, 0.4};

    BF bf(H_bf);
    bf.decode(b_bf);

    std::cout << '\n';

    // WBF
    WBF wbf(H_bf);
    wbf.decode(b_bf);

    std::cout << '\n';

    // MWBF
    MWBF mwbf(H_bf);
    mwbf.decode(b_bf);

    std::cout << '\n';

    // One Step MLG
    CtrlMat H1(15, 7, 5, "11010001", true);

    std::vector<double> b1 {-1.0, -1.0, -1.0, 1.0, -1.0, 1.0, -1.0, -1.0, -1.0,
                            1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

    OneStepMLG one_step_mlg(H1);
    one_step_mlg.decode(b1);

    std::cout << '\n';

    // Min Sum
    std::vector<std::vector<int>> H;
    H.push_back({1, 1, 1, 1, 0, 0, 0});
    H.push_back({1, 1, 0, 0, 1, 1, 0});
    H.push_back({1, 0, 1, 0, 1, 0, 1});

    CtrlMat H_hamming(7, 4, 3, H);

    MinSum min_sum(H_hamming);

    std::vector<double> b {0.4, 0.2, -1.0, 0.3, 0.8, -0.8, 0.1};
    min_sum.decode(b);
}
