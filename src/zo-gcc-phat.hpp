#ifndef __ZO_GCC_PHAT_HPP__
#define __ZO_GCC_PHAT_HPP__
#include <stdint.h>
#include <vector>

namespace zo {
    class GccPhat {
    public:

        /**
         * @brief Creates an instance of GccPhat
         * 
         * @return GccPhat* 
         */
        static GccPhat* create();

        /**
         * @brief Initialize GccPhat
         * 
         * @param sample_cnt The number of samples.
         */
        virtual void init(const int sample_cnt) = 0;

        /**
         * @brief Terminate and cleanup instance
         * 
         */
        virtual void terminate() = 0;

        /**
         * @brief Execute GCC PHAT algorithm.
         * See: http://www.xavieranguera.com/phdthesis/node92.html
         * See: https://github.com/respeaker/mic_array/blob/master/gcc_phat.py
         * 
         * @param signal First signal, say from microphone 1
         * @param refsignal Second signal, say from microphone 2
         * @return float "Tau" or distance signal is from refsignal
         */
        virtual float execute(const std::vector<int16_t>& signal, const std::vector<int16_t>& refsignal) = 0;
    };
} // namespace zo

#endif // __ZO_GCC_PHAT_HPP__
