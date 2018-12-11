// Set of functions to generate binaries and to find their orbital elements
#ifndef BINARY_H
#define BINARY_H
#include <array>
#include <tuple>

//Tested
// Find the eccentric anomaly of a binary given its eccentricity e and mean anomaly M
double eccentricAnomaly(double e, double M);

//Tested with setupRandomBinary
//Return the semi-major axis and eccentricity of a binary and whether or not it is bound from the positions and velocities of the stars
std::tuple<double, double, bool> orbitalElements(std::array<std::array<double, 3>, 4> X, double m1, double m2);

//Tested with orbitalElements
// Open a binary: find the position and velocity vectors given the semi-major axis and eccentricity
std::array<std::array<double, 3>, 4> setupRandomBinary(double a, double e, double m1, double m2);

#endif
