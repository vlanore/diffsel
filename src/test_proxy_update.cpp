// to compile (the reste of the project should have been compiled already):
//    clang++ --std=c++11 -I../utils test_proxy_update.cpp ../_build/libdiffsel_lib.a

#define DOCTEST_CONFIG_IMPLEMENT

#include "MSCodonSubMatrix.hpp"

using namespace std;

int main() {
    Eigen::VectorXd fitness0(5);
    fitness0 << 4, 6, 2, 12, 8;
    Eigen::VectorXd fitness(5);
    fitness << 1111, 1117, 1119, 1111, 1115;
    Eigen::Matrix<bool, Eigen::Dynamic, 1> ind_conv(5);
    ind_conv << true, false, true, true, false;

    SparseFitness myproxy(fitness0, fitness, ind_conv);

    cout << "==============================\n";
    for (int i = 0; i < 5; i++) {
        cout << "proxy[" << i << "] = " << myproxy.GetFitness(i) << '\n';
    }

    ind_conv(4) = true;
    ind_conv(0) = false;
    fitness0(1) = 33;
    fitness(2) = 1111111111;

    cout << "==============================\n";
    for (int i = 0; i < 5; i++) {
        cout << "proxy[" << i << "] = " << myproxy.GetFitness(i) << '\n';
    }
}
