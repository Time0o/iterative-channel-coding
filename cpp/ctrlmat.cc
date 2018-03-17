#include <algorithm>
#include <cassert>
#include <vector>

#include "ctrlmat.h"

std::ostream &operator<<(std::ostream &os, const CtrlMat& H)
{
    for (int i = 0; i < (H.orthogonal ? H.n : H.k); ++i) {
        os << '|';
        for (int j = 0; j < H.n; ++j) {
            std::vector<int> K = H.K[i];
            os << (std::find(K.begin(), K.end(), j) != K.end());
        }
        os << "|\n";
    }
    return os;
}

CtrlMat::CtrlMat(int n, int l, int dmin, std::vector<std::vector<int>> H,
    bool orthogonal) : n(n), l(l), k(n -l), dmin(dmin), orthogonal(orthogonal)
{
    for (int i = 0; i < (orthogonal ? n : k); ++i)
        K.push_back(std::vector<int>());

    for (int i = 0; i < n; ++i)
        N.push_back(std::vector<int>());

    for (int i = 0; i < (orthogonal ? n : k); ++i) {
        for (int j = 0; j < n; ++j) {
            if (H[i][j]) {
                K[i].push_back(j);
                N[j].push_back(i);
            }
        }
    }
}

CtrlMat::CtrlMat(int n, int l, int dmin, std::string h, bool orthogonal) :
    n(n), l(l), k(n -l), dmin(dmin), orthogonal(orthogonal)
{
    std::vector<int> hn(n);
    for (size_t j = 0u; j < h.size(); ++j)
        hn[j] = h[j] == '1' ? 1 : 0;

    std::vector<std::vector<int>> H;
    for (int j = 0; j < (orthogonal ? n : k); ++j) {
        H.push_back(hn);
        std::rotate(hn.begin(), hn.end() - 1, hn.end());
    }

    for (int i = 0; i < (orthogonal ? n : k); ++i)
        K.push_back(std::vector<int>());

    for (int i = 0; i < n; ++i)
        N.push_back(std::vector<int>());

    for (int i = 0; i < (orthogonal ? n : k); ++i) {
        for (int j = 0; j < n; ++j) {
            if (H[i][j]) {
                K[i].push_back(j);
                N[j].push_back(i);
            }
        }
    }
}
