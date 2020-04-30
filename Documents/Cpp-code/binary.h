// Set of functions to generate binaries and to find their orbital elements
#ifndef BINARY_H
#define BINARY_H
#include <array>
#include <tuple>

//Tested
// Find the eccentric anomaly of a binary given its eccentricity e and mean anomaly M
long double eccentricAnomaly(long double e, long double M, bool &non_converged_binary);

long double eccentricAnomalyIonised(long double e, long double M, bool notBound, bool &non_converged_binary);

//Tested with setupRandomBinary
//Return the semi-major axis and eccentricity of a binary and whether or not it is bound from the positions and velocities of the stars
std::tuple<long double, long double, bool> orbitalElements(std::array<std::array<long double, 3>, 4> X, long double m1, long double m2);
std::tuple<long double, long double, bool> orbitalElements(std::vector<std::array<long double, 3>> X, long double m1, long double m2);

std::tuple<long double, long double, long double, bool, long double> orbitalElementsIonised(std::array<std::array<long double, 3>, 4> X, long double m1, long double m2);
std::tuple<long double, long double, long double, bool, long double> orbitalElementsIonised(std::vector<std::array<long double, 3>> X, long double m1, long double m2);

//Tested with orbitalElements
// Open a binary: find the position and velocity vectors given the semi-major axis and eccentricity
std::array<std::array<long double, 3>, 4> setupRandomBinary(long double a, long double e, long double m1, long double m2);
std::vector<std::array<long double, 3>> setupRandomBinaryVector(long double a, long double e, long double m1, long double m2);

std::array<std::array<long double, 3>, 4> setupRandomBinaryIonised(long double a, long double e, long double m1, long double m2, long double E, bool notBound, bool &non_converged_binary);

#endif
