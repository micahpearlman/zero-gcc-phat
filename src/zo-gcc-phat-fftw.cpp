#include "zo-gcc-phat.hpp"
#include <fftw3.h>
#include <complex>
#include <algorithm>
#include <iostream>
#include <cfloat>
#include <string.h>

namespace zo {
    class GccPhatFFTWImpl : public GccPhat {
    public:
        void init(const int sample_cnt) override {
            
        }

        void terminate() override {

        }

        float execute(const std::vector<int16_t>& siga, const std::vector<int16_t>& sigb) override {
            std::vector<std::complex<float>> siga_fft;
            //make sure the length for the FFT is larger or equal than len(sig) + len(refsig)
            const size_t n = siga.size() + sigb.size();
            calculate_real_fft(siga_fft, siga, n);

            std::vector<std::complex<float>> sigb_fft;
            calculate_real_fft(sigb_fft, sigb, n);
            // for (auto& v : sigb_dft_conj) { // calculate conjugate
            //     v = std::conj(v);
            // }

            // R = SIG * REFSIG_CONJ
            std::vector<std::complex<float>> R;
            R.resize(siga_fft.size());            
            for (int i = 0; i < siga_fft.size(); i++) {
                std::complex<float> v = siga_fft[i] * std::conj(sigb_fft[i]);
                v = v / (std::abs(v) + FLT_MIN);
                R[i] = v;
            }

            // // IR = R / np.abs(R)
            // std::vector<std::complex<float>> IR;
            // for (int i = 0; i < R.size(); i++) {
            //     std::complex<float> v = R[i] / std::abs(R[i]);
            //     IR.push_back(v);
            // }

            // Inverse
            std::vector<float> cc;
            calculate_inverse_fft(cc, R, n);

            // calculate argmax
            int arg_max = std::distance(cc.begin(), std::max_element(cc.begin(), cc.end()));
            for (auto v : cc) {
                std::cout << v << "\n";
            }
            std::cout << "ArgMax: " << arg_max << "\n";

            return (float)arg_max;
        }

    protected:

        /**
         * @brief Calculate real valued signal fft
         * See: http://www.fftw.org/fftw3_doc/Complex-One_002dDimensional-DFTs.html
         * 
         * @param result DFT result 
         * @param signal audio signal
         * @param n size of input vector.  If n is smaller than the length of the input, 
         * the input is cropped. If it is larger, the input is padded with zeros.
         */
        void calculate_real_fft(std::vector<std::complex<float>>& result, const std::vector<int16_t>& signal, const int n) {
             
            // initialize input vector
            double* in = (double*)fftw_malloc(sizeof(double) * n);
            memset(in, 0, n);
            for(int i = 0; i < signal.size(); i++) {
                in[i] = (double)signal[i];
            }
            int nc = (n / 2) + 1;
            fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nc);

            fftw_plan plan = fftw_plan_dft_r2c_1d(n, in, out, FFTW_ESTIMATE);

            fftw_execute(plan);

            // copy result
            result.clear();
            result.resize(nc);
            for (int i = 0; i < nc; i++) {
                result[i] = std::complex<float>((float)out[i][0], (float)out[i][1]);
            }

            // cleanup
            fftw_destroy_plan(plan);
            fftw_free(in);
            fftw_free(out);
        }

        void calculate_inverse_fft(std::vector<float>& result, const std::vector<std::complex<float>>& signal, const int n) {

            // initialize input vector
            fftw_complex* in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * signal.size());
            memset(in, 0, signal.size());
            for(int i = 0; i < signal.size(); i++) {
                in[i][0] = (double)signal[i].real();
                in[i][1] = (double)signal[i].imag();
            }
            double* out = (double*)fftw_malloc(sizeof(double) * n);

            fftw_plan plan = fftw_plan_dft_c2r_1d(n, in, out, FFTW_ESTIMATE);

            fftw_execute(plan);

            // copy result
            result.clear();
            result.resize(n);
            for (int i = 0; i < n; i++) {
                result[i] = (float)out[i];
            }

            // cleanup
            fftw_destroy_plan(plan);
            fftw_free(in);
            fftw_free(out);
        }

    };

    GccPhat* GccPhat::create() {
        return new GccPhatFFTWImpl();
    }
} // namespace zo