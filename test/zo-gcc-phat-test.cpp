#include "zo-gcc-phat.hpp"
#include <vector>
#include <cmath>
#include <stdint.h>
#include <iostream>

void rotate(std::vector<int16_t>& out, const std::vector<int16_t>& in, int rot) {
    int n = in.size();
    for (unsigned int i=0; i<n ;i++) {
        unsigned int j = (i+rot)%n;
        out[j] = in[i];
    }
}

int main(int argc, char* argv[]) {
    zo::GccPhat* gcc_phat = zo::GccPhat::create();

    int    true_delay = -23;  /* The "true" delay value. This is the answer we are looking for. */
    double amp = 1.0;       /* amplitude of the sine wave - arbitrary */
    int    samp_freq = 512;   /* frequency of the sine wave - arbitrary */

    double samp_per = 1.0/(double)samp_freq; /* sample period of the test signals */
    int    nsamps = 8192;                  /* number of samples in the test signals */
    std::vector<int16_t> siga;                   /* storage for the test signal (sine wave) */
    siga.resize(nsamps);
    std::vector<int16_t> sigb;                   /* shifted version of test signal */
    sigb.resize(nsamps);

    /* Generate a test signal - a sine wave */
    for (unsigned i=0;i<nsamps;i++) {
        double v = amp * std::sin((2.0 * M_PI * samp_per * i));
        v = float(INT16_MAX) * v;
        int16_t vi = (int16_t)v;
        siga[i] = vi;

    }

    /* Generate a delayed version of the sine wave */
    rotate(sigb, siga, true_delay);

    int margin = 50;
    int distance = gcc_phat->execute(siga, sigb, 50);

    std::cout << "Delay: " << true_delay << "\n"
            << "Calculated Delay: " << distance << "\n";

    delete gcc_phat;
    
    return 1;
}