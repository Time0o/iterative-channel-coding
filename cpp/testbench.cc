#include <algorithm>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <vector>

#include "ctrlmat.h"
#include "blockdecoder.h"

#define EBN0_LOW 0.5
#define EBN0_HIGH 8
#define EBN0_STEP 0.1
#define N_ADAPT_ALPHA 20
#define N_AVG 500

static std::vector<double> awgn_vector(int n, double R, double EbN0)
{
    static std::default_random_engine awgn_re(time(0));

    double sigma = 1.0 / std::sqrt(2.0 * R * std::pow(10.0, EbN0 / 10.0));
    std::normal_distribution<> normal(0, sigma);

    std::vector<double> err(n, 1.0);

    for (double &e : err)
        e += normal(awgn_re);

    return err;
}

struct IterativeBenchmarkResult {
    std::vector<double> ebn0, alpha, time;
    std::vector<int> correct, incorrect, failures;
};

static std::ostream &operator<<(std::ostream &os,
    const IterativeBenchmarkResult &res)
{
    os << "EbN0,alpha,time,correct,incorrect,failures\n";

    for (auto i = 0u; i < res.ebn0.size(); ++i) {
        os << res.ebn0[i] << ','
           << res.alpha[i] << ','
           << res.time[i] << ','
           << res.correct[i] << ','
           << res.incorrect[i] << ','
           << res.failures[i] << '\n';
    }

    return os;
}

static IterativeBenchmarkResult benchmark_iterative(IterativeBlockdecoder &dec)
{
    IterativeBenchmarkResult res;

    CtrlMat H = dec.get_H();
    double R = (double) (H.n - H.k) / (double) H.n;

    double alpha_low = dec.get_alpha_low();
    double alpha_high = dec.get_alpha_high();
    double alpha_step = (alpha_high - alpha_low) / N_ADAPT_ALPHA;

    for (double ebn0 = EBN0_LOW; ebn0 <= EBN0_HIGH; ebn0 += EBN0_STEP) {

        double optimal_alpha = 0.0;
        int optimal_alpha_correct = -1;
        int optimal_alpha_incorrect = 0;
        int optimal_alpha_failures = 0;
        double optimal_alpha_time = 0.0;

        for (double alpha = alpha_low; alpha <= alpha_high; alpha += alpha_step) {
            dec.set_alpha(alpha);

            int correct = 0;
            int incorrect = 0;
            int failures = 0;
            double time = 0.0;

            for (int j = 0; j < N_AVG; ++j) {
                std::vector<double> y = awgn_vector(H.n, R, ebn0);

                if (dec.decode(y)) {
                    std::vector<int> b = dec.get_result();
                    if (std::find(b.begin(), b.end(), 1) == b.end()) {
                        ++correct;
                        time += dec.get_exec_time();
                    } else {
                        ++incorrect;
                    }
                } else {
                    ++failures;
                }
            }

            if (correct > optimal_alpha_correct) {
                optimal_alpha = alpha;
                optimal_alpha_correct = correct;
                optimal_alpha_incorrect = incorrect;
                optimal_alpha_failures = failures;
                optimal_alpha_time = time;
            }

            if (!dec.has_adaptable_alpha())
                break;
        }

        res.ebn0.push_back(ebn0);
        res.alpha.push_back(optimal_alpha);
        if (optimal_alpha_correct == 0)
            res.time.push_back(std::nan(""));
        else
            res.time.push_back(optimal_alpha_time / optimal_alpha_correct);

        res.correct.push_back(optimal_alpha_correct);
        res.incorrect.push_back(optimal_alpha_incorrect);
        res.failures.push_back(optimal_alpha_failures);
    }

    return res;
}

int main(int argc, char **argv)
{
    if (argc != 7) {
        std::cout << "Usage: " << argv[0]
                  << " n l dmin h (ortho|nonortho) csvdir\n";
        return -1;
    }

    int n = std::stoi(argv[1]);
    int l = std::stoi(argv[2]);
    int dmin = std::stoi(argv[3]);

    char const *h = argv[4];

    bool is_ortho;
    std::string ortho(argv[5]);
    if (ortho == "ortho")
        is_ortho = true;
    else if (ortho == "nonortho")
        is_ortho = false;
    else
        throw std::invalid_argument("did not recognize: '" + ortho + "'");

    CtrlMat H(n, l, dmin, h, is_ortho);

    std::string csvdir(argv[6]);
    if (csvdir.back() != '/')
        csvdir += '/';

    std::vector<std::unique_ptr<IterativeBlockdecoder>> decoders;
    std::vector<char const *> decoder_names;

    if (is_ortho) {
        decoders.push_back(std::make_unique<OneStepMLG>(H));
        decoder_names.push_back("one-step MLG");
        decoders.push_back(std::make_unique<HardMLG>(H));
        decoder_names.push_back("hard MLG");
        decoders.push_back(std::make_unique<SoftMLG>(H));
        decoder_names.push_back("soft MLG");
        decoders.push_back(std::make_unique<AdaptiveSoftMLG>(H));
        decoder_names.push_back("adaptive soft MLG");
    } else {
        decoders.push_back(std::make_unique<BF>(H));
        decoder_names.push_back("BF");
        decoders.push_back(std::make_unique<WBF>(H));
        decoder_names.push_back("WBF");
        decoders.push_back(std::make_unique<MWBF>(H));
        decoder_names.push_back("MWBF");
        decoders.push_back(std::make_unique<IMWBF>(H));
        decoder_names.push_back("IMWBF");
        decoders.push_back(std::make_unique<MinSum>(H));
        decoder_names.push_back("Min Sum");
        decoders.push_back(std::make_unique<NormalizedMinSum>(H));
        decoder_names.push_back("normalized Min Sum");
        decoders.push_back(std::make_unique<OffsetMinSum>(H));
        decoder_names.push_back("offset Min Sum");
    }

    for (size_t i = 0u; i < decoders.size(); ++i) {
        IterativeBlockdecoder &decoder = *decoders[i];
        char const *decoder_name = decoder_names[i];

        std::string csvfile(
            csvdir + std::to_string(i + 1) + ". " + decoder_name + ".csv");

        std::ifstream is(csvfile);
        if (is.good())
            continue;

        std::ofstream os(csvfile);
        std::cout << "creating: " << csvfile << "...\n";
        os << benchmark_iterative(decoder);
    }
}
