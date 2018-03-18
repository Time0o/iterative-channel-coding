#include <vector>

#include "ctrlmat.h"
#include "blockdecoder.h"

int main()
{
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
