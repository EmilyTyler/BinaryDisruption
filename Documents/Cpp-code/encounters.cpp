// Set of functions required to simulate encounters between binary stars and a sea of perturbers
#include <cmath>
#include <array>
#include <tuple>
#include <vector>

#include <iostream>

#include "constants.h"
#include "random_direction.h"
#include "vector_maths.h"
#include "binary.h"
#include "random_numbers.h"
#include "MC_velocity.h"
using namespace std;


//Tested
// Rate of encounters between impact parameters b0 and b1 and between relative speeds of v0 and v1 in a sea of particles with number density n_p and relative velocity dispersion v_rel
double encounterRate(double n_p, double v_rel, double b0, double b1, double v0, double v1)
{
	return sqrt(2.0*pi)*n_p/v_rel*(pow(b1, 2.0)-pow(b0, 2.0))*((pow(v0, 2.0)+2.0*pow(v_rel, 2.0))*exp(-pow(v0, 2.0)/(2.0*pow(v_rel, 2.0)))-(pow(v1, 2.0)+2.0*pow(v_rel, 2.0))*exp(-pow(v1, 2.0)/(2.0*pow(v_rel, 2.0))));
}

//Tested
// Calculates the impact parameter at which the fractional change in semi-major axis of the binary will be equal to delta, for perturber mass M_p, relative velocity dispersion v_rel, semi-major axis a and binary star masses m1 and m2.
double calcBMax(double M_p, double v_rel, double a, double m1, double m2, double delta = pow(10.0, -3.0))
{
	return pow((64.0*G*M_p*M_p*pow(a,3.0)/((m1+m2)*v_rel*v_rel*delta*delta)), 0.25);
}

//Tested magnitude, direction assumed to be correct from testing randomDirection
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

//Tested
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

	cout << "X = ";
	for(int i=0; i<4; i++){
		for(int j=0; j<3; j++){
			cout << " " << X[i][j];
		}
		cout << endl;
	}

	//Find impact parameter and velocity vectors
	tuple<array<double,3>, array<double,3>> bvvectors = impactAndVelocityVectors(b, v);
	array<double,3> b_vec = get<0>(bvvectors);
	array<double,3> v_vec = get<1>(bvvectors);

	cout << "b_vec = " << b_vec[0] << ", " << b_vec[1] << ", " << b_vec[2] << endl;
	cout << "v_vec = " << v_vec[0] << ", " << v_vec[1] << ", " << v_vec[2] << endl;

	//Declare variables
	double b_90, b_star_norm, v_perp, v_para;
	array<double,3> b_star;
	for (int i=0; i<2; i++){
		//90 degree deflection radius
		b_90 = G*(M_p + m[i])/(v*v);
		cout << "b_90 = " << b_90 << endl;
		//Calculate impact parameter for this star
		b_star = calcBStar(X[i], v_vec, b_vec);
		cout << "b_star = " << b_star[0] << ", " << b_star[1] << ", " << b_star[2] << endl;
		//Calculate norm of b_star
		b_star_norm = norm(b_star);
		cout << "norm b_star = " << b_star_norm << endl;
		//Calculate speed change in b_star direction
		v_perp = 2.0*M_p*v/(m[i]+M_p) * (b_star_norm/b_90)/(1.0 + b_star_norm*b_star_norm/(b_90*b_90));
		cout << "v_perp = " << v_perp << endl;
		//Calculate speed change in -v_vec direction
		v_para = 2.0*M_p*v/(m[i]+M_p) * 1.0/(1.0 + b_star_norm*b_star_norm/(b_90*b_90));
		cout << "v_para = " << v_para << endl;
		//Change star velocity
		for (int j=0; j<3; j++){
			X[i+2][j] += v_perp * b_star[j]/b_star_norm - v_para * v_vec[j]/v;
		}
	}

	cout << "X = ";
	for(int i=0; i<4; i++){
		for(int j=0; j<3; j++){
			cout << " " << X[i][j];
		}
		cout << endl;
	}

	//Close binary
	return orbitalElements(X, m1, m2);
}

double drawB(double b_max)
{
	return b_max*sqrt(randomUniformDoubleOpen());
}

tuple<vector<double>, vector<double>> MCEncounters(double v_rel, double n_p, double T, double m1, double m2, double M_p, vector<double> a, vector<double> e)
{
	//Minimum impact parameter
	double b_min = 0.0;
	//Minimum relative velocity of encounter
	double v_min = 0.01 * v_rel;
	//Maximum relative velocity of encounter
	double v_max = 100.0 * v_rel;
	//Maximum semimajor axis
	double a_T = 1000.0 * parsec/length_scale;
	//Number of binaries
	int N_bin = static_cast<int>(a.size());
	//Declare variable types
	double t, b_max, rate, v, b;
	tuple<double, double, bool> result;
	bool notBound;
	//Iterate over binaries
	for (int i=0; i<N_bin; i++){
		//Time passed
		t = 0.0;
		//Implement encounters
		while (t<=T){
			//Maximum impact parameter
			b_max = calcBMax(M_p, v_rel, a[i], m1, m2);
			//Encounter rate
			rate = encounterRate(n_p, v_rel, b_min, b_max, v_min, v_max);
			//Increment time passed
			t += randomExponential(rate);
			//Draw velocity from distribution
			v = drawVMaxwellian(v_rel, v_min, v_max);
			//Draw impact parameter from distribution
			b = drawB(b_max);
			//Encounter
			result = impulseEncounter(m1, m2, M_p, a[i], e[i], b, v);
			a[i] = get<0>(result);
			e[i] = get<1>(result);
			notBound = get<2>(result);
			if(notBound or (a[i]>=a_T)){
				a[i] = -1.0;
				e[i] = -1.0;
				break;
			}
		}
	}
	return make_tuple(a, e);
}