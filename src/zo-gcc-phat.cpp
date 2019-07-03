#include "zo-gcc-phat.hpp"
#include "zo-fft.hpp"
#include <complex>
#include <algorithm>
#include <iostream>
#include <cfloat>
#include <string.h>

namespace zo {
    class GccPhatImpl : public GccPhat {
    public:
        void init(const int sample_cnt) override {
            _forward_FFT = FFT_forward::create();
            _forward_FFT->init(sample_cnt);

            _inverse_FFT = FFT_inverse::create();
            _inverse_FFT->init(sample_cnt);
        }

        void terminate() override {
            _forward_FFT->terminate();
            delete _forward_FFT;
            _inverse_FFT->terminate();
            delete _inverse_FFT;
        }

        int execute(const std::vector<int16_t>& siga, const std::vector<int16_t>& sigb, int margin) override {
            std::vector<std::complex<double>> siga_fft;
            const size_t n = siga.size();
            //calculate_real_fft(siga_fft, siga, n);
            _forward_FFT->execute(siga_fft, siga);

            std::vector<std::complex<double>> sigb_fft;
            // calculate_real_fft(sigb_fft, sigb, n);
            _forward_FFT->execute(sigb_fft, sigb);

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
            // calculate_inverse_fft(cross_correlation, R, n);
            _inverse_FFT->execute(cross_correlation, R);

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

    protected:
        FFT_forward* _forward_FFT = nullptr;
        FFT_inverse* _inverse_FFT = nullptr;

    };

    GccPhat* GccPhat::create() {
        return new GccPhatImpl();
    }
} // namespace zo