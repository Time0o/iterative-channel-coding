#pragma once

#include <iostream>
#include <string>
#include <vector>

#include "ctrlmat.h"

#define DEFAULT_ALPHA 1.0
#define DEFAULT_MAX_ITER 10

class IterativeBlockdecoder
{
public:
    IterativeBlockdecoder(CtrlMat const &H) : H(H) {}
    IterativeBlockdecoder(CtrlMat const &H, double alpha) : H(H), alpha(alpha) {}

    bool decode(std::vector<double> const &y) {
        b = std::vector<int>(H.n);

        clock_t start = clock();
        bool success = _decode(y, b);
        exec_time = (clock() - start) / (double) CLOCKS_PER_SEC;

        return success;
    }

    double get_exec_time() const { return exec_time; }
    CtrlMat get_H() const { return H; }
    std::vector<int> get_result() const { return b; }

    void set_alpha(double alpha) { this->alpha = alpha; }

    virtual double get_alpha_low() const { return 0.0; }
    virtual double get_alpha_high() const { return 0.0; }
    virtual bool has_adaptable_alpha() const { return false; }

    virtual ~IterativeBlockdecoder() {}

protected:
    virtual bool _decode(std::vector<double> const &in, std::vector<int> &out) = 0;

    CtrlMat const &H;
    double alpha = DEFAULT_ALPHA;
    int max_iter = DEFAULT_MAX_ITER;

    std::vector<int> b;
    double exec_time;
};

class BF : public IterativeBlockdecoder
{
public:
    BF(CtrlMat const &H) : IterativeBlockdecoder(H) {}
private:
    bool _decode(std::vector<double> const &in, std::vector<int> &out);
};

class WBF : public IterativeBlockdecoder
{
public:
    WBF(CtrlMat const &H) : IterativeBlockdecoder(H) {}
private:
    bool _decode(std::vector<double> const &in, std::vector<int> &out);
};

class MWBF : public IterativeBlockdecoder
{
public:
    MWBF(CtrlMat const &H, double alpha = 0.3) : IterativeBlockdecoder(H, alpha) {}

    double get_alpha_low() const { return 0.3; }
    double get_alpha_high() const { return 2.0; }
    bool has_adaptable_alpha() const { return true; }
private:
    bool _decode(std::vector<double> const &in, std::vector<int> &out);
};

class IMWBF : public IterativeBlockdecoder
{
public:
    IMWBF(CtrlMat const &H, double alpha = 0.3) : IterativeBlockdecoder(H, alpha) {}

    double get_alpha_low() const { return 0.3; }
    double get_alpha_high() const { return 2.0; }
    bool has_adaptable_alpha() const { return true; }
private:
    bool _decode(std::vector<double> const &in, std::vector<int> &out);
};

class HardMLG : public IterativeBlockdecoder
{
public:
    HardMLG(CtrlMat const &H) : IterativeBlockdecoder(H) {}
private:
    bool _decode(std::vector<double> const &in, std::vector<int> &out);
};

class SoftMLG : public IterativeBlockdecoder
{
public:
    SoftMLG(CtrlMat const &H) : IterativeBlockdecoder(H) {}
private:
    bool _decode(std::vector<double> const &in, std::vector<int> &out);
};

class AdaptiveSoftMLG : public IterativeBlockdecoder
{
public:
    AdaptiveSoftMLG(CtrlMat const &H) : IterativeBlockdecoder(H) {}

    double get_alpha_low() const { return 0.0; }
    double get_alpha_high() const { return 1.0; }
    bool has_adaptable_alpha() const { return true; }
private:
    bool _decode(std::vector<double> const &in, std::vector<int> &out);
};

class MinSum : public IterativeBlockdecoder
{
public:
    MinSum(CtrlMat const &H) : IterativeBlockdecoder(H) {}
private:
    bool _decode(std::vector<double> const &in, std::vector<int> &out);
};

class NormalizedMinSum : public IterativeBlockdecoder
{
public:
    NormalizedMinSum(CtrlMat const &H) : IterativeBlockdecoder(H) {}

    double get_alpha_low() const { return 1.25; }
    double get_alpha_high() const { return 2.0; }
    bool has_adaptable_alpha() const { return true; }
private:
    bool _decode(std::vector<double> const &in, std::vector<int> &out);
};

class OffsetMinSum : public IterativeBlockdecoder
{
public:
    OffsetMinSum(CtrlMat const &H) : IterativeBlockdecoder(H) {}

    double get_alpha_low() const { return 0.15; }
    double get_alpha_high() const { return 0.20; }
    bool has_adaptable_alpha() const { return true; }
private:
    bool _decode(std::vector<double> const &in, std::vector<int> &out);
};
