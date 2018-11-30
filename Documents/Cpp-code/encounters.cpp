// Set of functions required to simulate encounters between binary stars and a sea of perturbers
#include <cmath>
#include <array>
#include <tuple>
#include "constants.h"
#include "random_direction.h"
#include "vector_maths.h"
#include "binary.h"
using namespace std;

// Rate of encounters between impact parameters b0 and b1 and between relative speeds of v0 and v1 in a sea of particles with number density n_p and relative velocity dispersion v_rel
double encounterRate(double n_p, double v_rel, double b0, double b1, double v0, double v1)
{
	return sqrt(2.0*pi)*n_p/v_rel*(pow(b1, 2.0)-pow(b0, 2.0))*((pow(v0, 2.0)+2.0*pow(v_rel, 2.0))*exp(-pow(v0, 2.0)/(2.0*pow(v_rel, 2.0)))-(pow(v1, 2.0)+2.0*pow(v_rel, 2.0))*exp(-pow(v1, 2.0)/(2.0*pow(v_rel, 2.0))));
}

// Calculates the impact parameter at which the fractional change in semi-major axis of the binary will be equal to delta, for perturber mass M_p, relative velocity dispersion v_rel, semi-major axis a and binary star masses m1 and m2.
// Should we change this to the definition from Weinberg et al.?
double calc_b_max(double M_p, double v_rel, double a, double m1, double m2, double delta = pow(10.0, -3.0))
{
	return sqrt(2.0)* pow(G, 0.25)* sqrt(M_p)* pow(m1+m2, -0.25)* pow(a, 0.75)* pow(v_rel, -0.5)* pow(delta, -0.25);
}

// Finds the impact parameter and velocity vectors given the magnitudes of both and that they should be randomly distributed and perpendicular
tuple<array<double,3>, array<double,3>> impactAndVelocityVectors(double b, double v)
{
	// Velocity vector
	array<double, 3> v_vec = randomDirection();
	for (int i = 0; i < 3; i++)
    	v_vec[i] *= v;
	// Other random vector
	array<double, 3> n = randomDirection();
	// Impact parameter vector
	array<double, 3> b_vec = cross(v_vec, n);
	// Correct magnitude of impact parameter vector
	b_vec = normalise(b_vec);
	for (int i = 0; i < 3; i++)
    	b_vec[i] *= b;
return make_tuple(b_vec, v_vec);
}

//Finds the impact parameter for a star in a binary given the impact parameter b_vec, velocity of perturber v_vec, and star position x
array<double,3> calcBStar(array<double, 3> x, array<double, 3> v_vec, array<double, 3> b_vec)
{
	array<double,3> b_star;
	double v = norm(v_vec);
	for (int i=0; i<3; i++){
		b_star[i] = dot(x,v_vec)/(v*v) * v_vec[i] + b_vec[i] - x[i];
	}
	return b_star;
}

// Implements an encounter at impact parameter b and relative speed v
tuple<double, double, bool> impulseEncounter(double m1, double m2, double M_p, double a, double e, double b, double v)
{
	//Star masses
	array<double, 2> m = {{m1, m2}};
	//Open binary
	array<array<double, 3>, 4> X = setupRandomBinary(a, e, m1, m2);
	//Find impact parameter and velocity vectors
	tuple<array<double,3>, array<double,3>> bvvectors = impactAndVelocityVectors(b, v);
	array<double,3> b_vec = get<0>(bvvectors);
	array<double,3> v_vec = get<1>(bvvectors);
	//Declare variables
	double b_90, b_star_norm, v_perp, v_para;
	array<double,3> b_star;
	for (int i=0; i<2; i++){
		//90 degree deflection radius
		b_90 = G*(M_p + m[i])/(v*v);
		//Calculate impact parameter for this star
		b_star = calcBStar(X[i], v_vec, b_vec);
		//Calculate norm of b_star
		b_star_norm = norm(b_star);
		//Calculate speed change in b_star direction
		v_perp = 2.0*M_p*v/(m[i]+M_p) * (b_star_norm/b_90)/(1.0 + b_star_norm*b_star_norm/b_90*b_90);
		//Calculate speed change in -v_vec direction
		v_para = 2.0*M_p*v/(m[i]+M_p) * 1.0/(1.0 + b_star_norm*b_star_norm/b_90*b_90);
		//Change star velocity
		for (int j=0; j<3; j++){
			X[i+2][j] += v_perp * b_star[j]/b_star_norm - v_para * v_vec[j]/v;
		}
	}
	//Close binary
	return orbitalElements(X, m1, m2);
}