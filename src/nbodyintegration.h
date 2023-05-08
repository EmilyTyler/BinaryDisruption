//
#ifndef NBODYINTEGRATION_H
#define NBODYINTEGRATION_H
#include <array>
#include <vector>

void initialise_arrays(int N);

std::vector<std::array<long double, 3>> evolve(int N, std::vector<long double> M, std::vector<std::array<long double, 3>> X, long double T, int n=1, long double eta = 0.02, bool ini_arrays = true);

#endif