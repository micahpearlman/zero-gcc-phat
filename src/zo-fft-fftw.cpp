#include "zo-fft.hpp"
#include <fftw3.h>
#include <stdint.h>

namespace zo {
    class FFT_forward_fftw : public FFT_forward {
    public:
        void init(const int n) override {
            _insize = n;
            _in = (double*)fftw_malloc(sizeof(double) * n);

            const int nc = (n / 2) + 1;
            _outsize = nc;
            _out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nc);

            _plan = fftw_plan_dft_r2c_1d(n, _in, _out, FFTW_ESTIMATE);

        }

        void terminate() override {
            fftw_destroy_plan(_plan);
            fftw_free(_in);
            fftw_free(_out);
        }

        void execute(std::vector<std::complex<double>>& out, const std::vector<int16_t>& in) {
            for (int i = 0; i < in.size(); i++) {
                _in[i] = (double)in[i];
            }

            fftw_execute(_plan);

            out.clear();
            out.resize(_outsize);
            for (int i = 0; i < _outsize; i++) {
                out[i] = std::complex<double>((double)_out[i][0], (double)_out[i][1]);
            }
        }
    protected:
        int _insize = 0;
        int _outsize = 0;
        double* _in = nullptr;
        fftw_complex* _out = nullptr;
        fftw_plan _plan;
    };

    FFT_forward* FFT_forward::create() {
        return new FFT_forward_fftw();
    }

    class FFT_inverse_fftw : public FFT_inverse {
    public:
        void init(const int n) override {
            _outsize = n;
            _out = (double*)fftw_malloc(sizeof(double) * n);

            const int nc = (n / 2) + 1;
            _insize = nc;
            _in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nc);

            _plan = fftw_plan_dft_c2r_1d(n, _in, _out, FFTW_ESTIMATE);

        }

        void terminate() override {
            fftw_destroy_plan(_plan);
            fftw_free(_in);
            fftw_free(_out);
        }

        void execute(std::vector<double>& out, const std::vector<std::complex<double>>& in) {
            for(int i = 0; i < in.size(); i++) {
                _in[i][0] = (double)in[i].real();
                _in[i][1] = (double)in[i].imag();
            }

            fftw_execute(_plan);

            out.clear();
            out.resize(_outsize);
            for (int i = 0; i < _outsize; i++) {
                out[i] = _out[i];
            }
        }

    protected:
        int _insize = 0;
        int _outsize = 0;
        double* _out = nullptr;
        fftw_complex* _in = nullptr;
        fftw_plan _plan;

        
    };

    FFT_inverse* FFT_inverse::create() {
        return new FFT_inverse_fftw();
    }

}