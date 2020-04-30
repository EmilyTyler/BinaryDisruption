// Set of functions required to simulate encounters between binary stars and a sea of perturbers
#ifndef ENCOUNTERS_H
#define ENCOUNTERS_H
#include <array>
#include <tuple>
#include <vector>

//Tested
// Rate of encounters between impact parameters b0 and b1 and between relative speeds of v0 and v1 in a sea of particles with number density n_p and relative velocity dispersion v_rel
long double encounterRate(long double n_p, long double v_rel, long double b0, long double b1, long double v0, long double v1);

//Tested
// Calculates the impact parameter at which the fractional change in semi-major axis of the binary will be equal to delta, for perturber mass M_p, relative velocity dispersion v_rel, semi-major axis a and binary star masses m1 and m2.
long double calcBMax(long double M_p, long double v_rel, long double a, long double m1, long double m2, long double delta = pow(10.0, -3.0));

//Tested magnitude, direction assumed to be correct from testing randomDirection
// Finds the impact parameter and velocity vectors given the magnitudes of both and that they should be randomly distributed and perpendicular
std::tuple<std::array<long double,3>, std::array<long double,3>> impactAndVelocityVectors(long double b, long double v);

//Tested
//Finds the impact parameter for a star in a binary given the impact parameter b_vec, velocity of perturber v_vec, and star position x
std::array<long double,3> calcBStar(std::array<long double, 3> x, std::array<long double, 3> v_vec, long double v_norm, std::array<long double, 3> b_vec);

//Tested
// Implements an encounter at impact parameter b and relative speed v
std::tuple<long double, long double, bool> impulseEncounter(long double m1, long double m2, long double M_p, long double a, long double e, long double b, long double v);

//Tested
//Draw an impact parameter from a distribution linear in b up to b_max
long double drawB(long double b_max, long double b_min=0.0);

std::tuple<std::vector<long double>, std::vector<long double>, std::vector<long double>, std::vector<long double>> MCEncountersIonised(long double v_rel, long double n_p, long double T, long double m1, long double m2, long double M_p, std::vector<long double> a, std::vector<long double> e);


#endif