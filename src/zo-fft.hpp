#ifndef __ZO_FFT_HPP__
#define __ZO_FFT_HPP__

#include <vector>
#include <complex>
#include <cstdint>

namespace zo {
    class FFT_forward {
    public:


        /**
         * @brief Creates a forward FFT instance
         * 
         * @return FFT_forward* forward FFT instance
         */
        static FFT_forward* create();

        /**
         * @brief 
         * 
         * @param sample_cnt 
         */
        virtual void init(const int sample_cnt) = 0;

        /**
         * @brief 
         * 
         */
        virtual void terminate() = 0;

        /**
         * @brief 
         * 
         * @param out 
         * @param in 
         */
        virtual void execute(std::vector<std::complex<double>>& out, const std::vector<int16_t>& in) = 0;
    };

    class FFT_inverse {
    public:

        /**
         * @brief 
         * 
         * @return FFT_inverse* 
         */
        static FFT_inverse* create();

        /**
         * @brief 
         * 
         * @param sample_cnt 
         */
        virtual void init(const int sample_cnt) = 0;

        /**
         * @brief 
         * 
         */
        virtual void terminate() = 0;

        /**
         * @brief 
         * 
         * @param out 
         * @param in 
         */
        virtual void execute(std::vector<double>& out, const std::vector<std::complex<double>>& in) = 0;
    };

}

#endif // __ZO_FFT_HPP__


