// Set of functions to generate binaries and to find their orbital elements
#ifndef BINARY_H
#define BINARY_H
#include <array>
#include <tuple>

// Find the eccentric anomaly of a binary given its eccentricity e and mean anomaly M
double eccentricAnomaly(double e, double M);

//Return the semi-major axis and eccentricity of a binary and whether or not it is bound from the positions and velocities of the stars
tuple<double, double, bool> orbitalElements(array<array<double, 3>, 4> X, double m1, double m2);

// Open a binary: find the position and velocity vectors given the semi-major axis and eccentricity
array<array<double, 3>, 4> setupRandomBinary(double a, double e, double m1, double m2);

#endif
