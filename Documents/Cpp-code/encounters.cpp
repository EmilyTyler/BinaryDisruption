// Set of functions required to simulate encounters between binary stars and a sea of perturbers
#include <math.h>
#include <array>
#include <tuple>
#include "constants.h"
#include "random_direction.h"
#include "vector_maths.h"
using namespace std;

// Rate of encounters between impact parameters b0 and b1 and between relative speeds of v0 and v1 in a sea of particles with number density n_p and relative velocity dispersion v_rel
 double encounterRate(double n_p, double v_rel, double b0, double b1, double v0, double v1)
 {
 	return sqrt(2.0*pi)*n_p/v_rel*(pow(b1, 2.0)-pow(b0, 2.0))*((pow(v0, 2.0)+2.0*pow(v_rel, 2.0))*exp(-pow(v0, 2.0)/(2.0*pow(v_rel, 2.0)))-(pow(v1, 2.0)+2.0*pow(v_rel, 2.0))*exp(-pow(v1, 2.0)/(2.0*pow(v_rel, 2.0))));
 }

// Calculates the impact parameter at which the fractional change in semi-major axis of the binary will be equal to delta, for perturber mass M_p, relative velocity dispersion v_rel, semi-major axis a and binary star masses m1 and m2.
 double calc_b_max(double M_p, double v_rel, double a, double m1, double m2, double delta = pow(10.0, -3.0))
 {
 	return sqrt(2.0)* pow(G, 0.25)* sqrt(M_p)* pow(m1+m2, -0.25)* pow(a, 0.75)* pow(v_rel, -0.5)* pow(delta, -0.25);
 }

  // Finds the impact parameter and velocity vectors given the magnitudes of both and that they should be randomly distributed and perpendicular
 tuple<array<double,3>, array<double,3>> impactAndVelocityVectors(double b, double v)
 {
 	// Velocity vector
 	array<double, 3> v_vec = v* randomDirection();
 	// Other random vector
 	array<double, 3> n = randomDirection();
 	// Impact parameter vector
 	array<double, 3> b_vec = cross(v_vec, n);
 	// Correct magnitude of impact parameter vector
 	double b_vec_norm = norm(b_vec);
 	for (int i=0; i<3; i++)
 	{
 		b_vec[i] *= b/b_vec_norm;
 	}
 }