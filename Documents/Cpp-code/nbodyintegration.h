//
#ifndef NBODYINTEGRATION_H
#define NBODYINTEGRATION_H

std::array<std::array<long double, 3>, N*2> evolve(const int N, std::array<long double,N> M, std::array<std::array<long double, 3>, N*2> X, long double T, int n=1, long double eta = 0.02);

#endif