#ifndef VECTOR_MATHS_H
#define VECTOR_MATHS_H
#include <array>
#include <vector>

//Tested
std::array<double, 3> cross(std::array<double,3> x, std::array<double,3> y);
std::array<long double, 3> cross(std::array<long double,3> x, std::array<long double,3> y);

//Tested
double norm(std::array<double,3> x);
long double norm(std::array<long double,3> x);

//Tested
std::array<double, 3> normalise(std::array<double, 3> x);
std::array<long double, 3> normalise(std::array<long double, 3> x);

//Tested
double dot(std::array<double,3> x, std::array<double,3> y);
long double dot(std::array<long double,3> x, std::array<long double,3> y);

std::vector<long double> where_positive(std::vector<long double> x);

std::vector<long double> rotate(std::vector<long double> vec, long double x_angle, long double y_angle, long double z_angle);
std::array<long double, 3> rotate(std::array<long double, 3> vec, long double x_angle, long double y_angle, long double z_angle);

std::vector<std::array<long double,3>> x_yPlaneAndAlignX(std::vector<std::array<long double,3>> X, long double x_align, long double m1, long double m2);

#endif