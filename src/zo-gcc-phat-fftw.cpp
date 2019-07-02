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

        int execute(const std::vector<int16_t>& siga, const std::vector<int16_t>& sigb, int margin) override {
            std::vector<std::complex<double>> siga_fft;
            //make sure the length for the FFT is larger or equal than len(sig) + len(refsig)
            const size_t n = siga.size();// + sigb.size();
            calculate_real_fft(siga_fft, siga, n);

            std::vector<std::complex<double>> sigb_fft;
            calculate_real_fft(sigb_fft, sigb, n);
            // for (auto& v : sigb_dft_conj) { // calculate conjugate
            //     v = std::conj(v);
            // }

            // R = SIG * REFSIG_CONJ
            std::vector<std::complex<double>> R;
            R.resize(siga_fft.size());            
            for (int i = 0; i < siga_fft.size(); i++) {
                std::complex<double> v = sigb_fft[i] * std::conj(siga_fft[i]);
                v = v / (std::abs(v) + FLT_MIN);
                R[i] = v;
            }


            // Inverse
            std::vector<double> cross_correlation;
            calculate_inverse_fft(cross_correlation, R, n);

            /*
            * Shift the values in xcorr[] so that the 0th lag is at the center of
            * the output array. 
            * [Note: the index of the center value in the output will be: ceil(_N/2) ]
            */
            std::vector<double> shifted;
            shift<double>(shifted, cross_correlation);

            // First, make sure the margin is within the bounds of the computed lags 
            
            unsigned center_i = ceil(n/2.0); 
            unsigned newmargin=margin;
            if (((int)(center_i - newmargin)) < 0) {
                newmargin = center_i;
            }
            if ((center_i + newmargin) >= n) {
                newmargin = (n-1) - center_i;
            }

            /* Compute the begin index and length of the lags_loc[] array */
            unsigned start_i = center_i-newmargin;
            unsigned len = 2*newmargin+1;

            // calculate argmax
            int max_index = std::distance(shifted.begin() + start_i, std::max_element(shifted.begin() + start_i, shifted.begin() + start_i + len));
            int arg_max = max_index - newmargin;

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
        void calculate_real_fft(std::vector<std::complex<double>>& result, const std::vector<int16_t>& signal, const int n) {
             
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
                result[i] = std::complex<double>((double)out[i][0], (double)out[i][1]);
            }

            // cleanup
            fftw_destroy_plan(plan);
            fftw_free(in);
            fftw_free(out);
        }

        void calculate_inverse_fft(std::vector<double>& result, const std::vector<std::complex<double>>& signal, const int n) {

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
                result[i] = (double)out[i];
            }

            // cleanup
            fftw_destroy_plan(plan);
            fftw_free(in);
            fftw_free(out);
        }

        /**
         * Shift the output of an FFT.
         *
         * The index of the mid-point in the output will be located at: ceil(_N/2)
         * @ingroup GCC
         */
        template<typename SCALAR=double>
        void shift(std::vector<SCALAR>& out, const std::vector<SCALAR>& in) {
            const size_t N = in.size();
            
            // mid-point of out[] will be located at index ceil(N/2) 
            const size_t xx = (size_t) std::floor((SCALAR) N/2.0);

            auto in_xx_iter = in.begin() + xx;

            out.clear();

            // Copy last half of in[] to first half of out[] 
            std::copy(in_xx_iter, in.end(), std::back_inserter(out));

            // Copy first half of in[] to end of out[] 
            std::copy(in.begin(), in_xx_iter, std::back_inserter(out));

            // Copy last half of in[] to first half of out[] 
            // memcpy(&out[0],&in[xx],sizeof(float)*(N-xx));

            // Copy first half of in[] to end of out[] 
            // memcpy(&out[N-xx],&in[0],sizeof(float)*xx);
        }


    };

    GccPhat* GccPhat::create() {
        return new GccPhatFFTWImpl();
    }
} // namespace zo