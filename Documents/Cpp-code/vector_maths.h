#ifndef VECTOR_MATHS_H
#define VECTOR_MATHS_H
#include <array>
using namespace std;

array<double, 3> cross(array<double,3> x, array<double,3> y);

double norm(array<double,3> x);

array<double, 3> normalise(array<double, 3> x);

#endif