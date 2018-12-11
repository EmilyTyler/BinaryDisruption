// Set of functions required to simulate encounters between binary stars and a sea of perturbers
#ifndef ENCOUNTERS_H
#define ENCOUNTERS_H
#include <array>
#include <tuple>
#include <vector>

//Tested
// Rate of encounters between impact parameters b0 and b1 and between relative speeds of v0 and v1 in a sea of particles with number density n_p and relative velocity dispersion v_rel
double encounterRate(double n_p, double v_rel, double b0, double b1, double v0, double v1);

//Tested
// Calculates the impact parameter at which the fractional change in semi-major axis of the binary will be equal to delta, for perturber mass M_p, relative velocity dispersion v_rel, semi-major axis a and binary star masses m1 and m2.
double calcBMax(double M_p, double v_rel, double a, double m1, double m2, double delta = pow(10.0, -3.0));

//Tested magnitude, direction assumed to be correct from testing randomDirection
// Finds the impact parameter and velocity vectors given the magnitudes of both and that they should be randomly distributed and perpendicular
std::tuple<std::array<double,3>, std::array<double,3>> impactAndVelocityVectors(double b, double v);

//Tested
//Finds the impact parameter for a star in a binary given the impact parameter b_vec, velocity of perturber v_vec, and star position x
std::array<double,3> calcBStar(std::array<double, 3> x, std::array<double, 3> v_vec, std::array<double, 3> b_vec);

// Implements an encounter at impact parameter b and relative speed v
std::tuple<double, double, bool> impulseEncounter(double m1, double m2, double M_p, double a, double e, double b, double v);

//Draw an impact parameter from a distribution linear in b up to b_max
double drawB(double b_max);

//Evolves a population of binaries (a,e) forwards by time T
std::tuple<std::vector<double>, std::vector<double>> MCEncounters(double v_rel, double n_p, double T, double m1, double m2, double M_p, std::vector<double> a, std::vector<double> e);

#endif