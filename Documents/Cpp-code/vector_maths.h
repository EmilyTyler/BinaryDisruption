#ifndef VECTOR_MATHS_H
#define VECTOR_MATHS_H
#include <array>

std::array<double, 3> cross(std::array<double,3> x, std::array<double,3> y);

double norm(std::array<double,3> x);

std::array<double, 3> normalise(std::array<double, 3> x);

double dot(std::array<double,3> x, std::array<double,3> y);

#endif