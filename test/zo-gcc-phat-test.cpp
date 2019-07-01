#include "zo-gcc-phat.hpp"

int main(int argc, char* argv[]) {
    zo::GccPhat* gcc_phat = zo::GccPhat::create();

    std::vector<int16_t> signal = {0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    std::vector<int16_t> refsignal = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    refsignal.resize(refsignal.size() * 2);
    gcc_phat->execute(signal, refsignal);

    delete gcc_phat;
    
    return 1;
}