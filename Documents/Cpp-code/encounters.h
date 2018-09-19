// Set of functions required to simulate encounters between binary stars and a sea of perturbers
#ifndef ENCOUNTERS_H
#define ENCOUNTERS_H
#include <array>
#include <tuple>
using namespace std;

// Rate of encounters between impact parameters b0 and b1 and between relative speeds of v0 and v1 in a sea of particles with number density n_p and relative velocity dispersion v_rel
double encounterRate(double n_p, double v_rel, double b0, double b1, double v0, double v1);

// Calcultates the impact parameter at which the fractional change in semi-major axis of the binary will be equal to delta, for perturber mass M_p, relative velocity dispersion v_rel, semi-major axis a and binary star masses m1 and m2.
 double calc_b_max(double M_p, double v_rel, double a, double m1, double m2, double delta = pow(10.0, -3.0));

 // Finds the impact parameter and velocity vectors given the magnitudes of both and that they should be randomly distributed and perpendicular
 tuple<array<double,3>, array<double,3>> impactAndVelocityVectors(double b, double v);

#endif