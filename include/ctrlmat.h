#pragma once

#include <iostream>
#include <string>
#include <vector>

class CtrlMat
{
friend std::ostream &operator<<(std::ostream &os, const CtrlMat& ctrl_mat);

public:
    CtrlMat(int n, int l, int dmin, std::vector<std::vector<int>> H,
            bool orthogonal = false);

    CtrlMat(int n, int l, int dmin, std::string h, bool orthogonal = false);

    int n, l, k, dmin;
    bool orthogonal;
    std::vector<std::vector<int>> K, N;
};
