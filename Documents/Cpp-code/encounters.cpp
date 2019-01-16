// Set of functions required to simulate encounters between binary stars and a sea of perturbers
#include <cmath>
#include <array>
#include <tuple>
#include <vector>

#include <iostream>
#include <iomanip>
#include <ctime>

#include "constants.h"
#include "random_direction.h"
#include "vector_maths.h"
#include "binary.h"
#include "random_numbers.h"
#include "MC_velocity.h"
using namespace std;


//Tested
// Rate of encounters between impact parameters b0 and b1 and between relative speeds of v0 and v1 in a sea of particles with number density n_p and relative velocity dispersion v_rel
long double encounterRate(long double n_p, long double v_rel, long double b0, long double b1, long double v0, long double v1)
{
	return sqrt(2.0*pi)*n_p/v_rel*(pow(b1, 2.0)-pow(b0, 2.0))*((pow(v0, 2.0)+2.0*pow(v_rel, 2.0))*exp(-pow(v0, 2.0)/(2.0*pow(v_rel, 2.0)))-(pow(v1, 2.0)+2.0*pow(v_rel, 2.0))*exp(-pow(v1, 2.0)/(2.0*pow(v_rel, 2.0))));
}

//Tested
// Calculates the impact parameter at which the fractional change in semi-major axis of the binary will be equal to delta, for perturber mass M_p, relative velocity dispersion v_rel, semi-major axis a and binary star masses m1 and m2.
long double calcBMax(long double M_p, long double v_rel, long double a, long double m1, long double m2, long double delta = pow(10.0, -3.0))
{
	return pow((64.0*G*M_p*M_p*pow(a,3.0)/((m1+m2)*v_rel*v_rel*delta*delta)), 0.25);
}

//Tested magnitude, direction assumed to be correct from testing randomDirection
// Finds the impact parameter and velocity vectors given the magnitudes of both and that they should be randomly distributed and perpendicular
tuple<array<long double,3>, array<long double,3>> impactAndVelocityVectors(long double b, long double v)
{
	// Velocity vector
	array<long double, 3> v_vec = randomDirection();
	for (int i = 0; i < 3; ++i)
    	v_vec[i] *= v;
	// Other random vector
	array<long double, 3> n = randomDirection();
	// Impact parameter vector
	array<long double, 3> b_vec = cross(v_vec, n);
	// Correct magnitude of impact parameter vector
	b_vec = normalise(b_vec);
	for (int i = 0; i < 3; ++i)
    	b_vec[i] *= b;
return make_tuple(b_vec, v_vec);
}

//Tested
//Finds the impact parameter for a star in a binary given the impact parameter b_vec, velocity of perturber v_vec, and star position x
array<long double,3> calcBStar(array<long double, 3> x, array<long double, 3> v_vec, long double v_norm, array<long double, 3> b_vec)
{
	array<long double,3> b_star;
	for (int i=0; i<3; ++i){
		b_star[i] = dot(x,v_vec)/(v_norm*v_norm) * v_vec[i] + b_vec[i] - x[i];
	}
	return b_star;
}

//Tested
// Implements an encounter at impact parameter b and relative speed v
tuple<long double, long double, bool> impulseEncounter(long double m1, long double m2, long double M_p, long double a, long double e, long double b, long double v)
{
	//Star masses
	array<long double, 2> m = {m1, m2};
	//Open binary
	array<array<long double, 3>, 4> X = setupRandomBinary(a, e, m1, m2);
	//Find impact parameter and velocity vectors
	tuple<array<long double,3>, array<long double,3>> bvvectors = impactAndVelocityVectors(b, v);
	array<long double,3> b_vec = get<0>(bvvectors);
	array<long double,3> v_vec = get<1>(bvvectors);
	//Declare variables
	long double b_90, b_star_norm, v_perp, v_para;
	array<long double,3> b_star;
	for (int i=0; i<2; ++i){
		//90 degree deflection radius
		b_90 = G*(M_p + m[i])/(v*v);
		//Calculate impact parameter for this star
		b_star = calcBStar(X[i], v_vec, v, b_vec);
		//Calculate norm of b_star
		b_star_norm = norm(b_star);
		//Calculate speed change in b_star direction
		v_perp = 2.0*M_p*v/(m[i]+M_p) * (b_star_norm/b_90)/(1.0 + b_star_norm*b_star_norm/(b_90*b_90));
		//Calculate speed change in -v_vec direction
		v_para = 2.0*M_p*v/(m[i]+M_p) * 1.0/(1.0 + b_star_norm*b_star_norm/(b_90*b_90));
		//Change star velocity
		for (int j=0; j<3; ++j){
			X[i+2][j] += v_perp * b_star[j]/b_star_norm - v_para * v_vec[j]/v;
		}
	}
	//Close binary
	return orbitalElements(X, m1, m2);
}

//Tested
//Draw an impact parameter from a distribution linear in b up to b_max
long double drawB(long double b_max)
{
	return b_max*sqrt(randomUniformDoubleClosed(0.0, 1.0));
}

tuple<vector<long double>, vector<long double>> MCEncounters(long double v_rel, long double n_p, long double T, long double m1, long double m2, long double M_p, vector<long double> a, vector<long double> e)
{
	//Minimum impact parameter
	long double b_min = 0.0;
	//Minimum relative velocity of encounter
	long double v_min = 0.0;
	//Maximum relative velocity of encounter
	long double v_max = 100.0 * v_rel;
	//Maximum semimajor axis
	long double a_T = 1000.0 * parsec/length_scale;
	//Number of binaries
	int N_bin = static_cast<int>(a.size());
	//Declare variable types
	long double t, b_max, rate, v, b;
	tuple<long double, long double, bool> result;
	bool notBound;
	double t_start;
	//Iterate over binaries
	for (int i=0; i<N_bin; ++i){
		cout << "Binary " << i+1 << " of " << N_bin << endl;
		//Start time
		clock_t t_start;
		t_start = clock();
		//Time passed
		t = 0.0;
		//Implement encounters
		while (t<T){
			//Maximum impact parameter
			b_max = calcBMax(M_p, v_rel, a[i], m1, m2);
			//Encounter rate
			rate = encounterRate(n_p, v_rel, b_min, b_max, v_min, v_max);
			//Increment time passed
			t += randomExponential(rate);
			//Draw velocity from distribution
			v = drawVMaxwellian(v_rel, v_max);
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
		//Print how long it took
		cout << "Time taken = " << (clock() - t_start)/(double)(CLOCKS_PER_SEC) << " s" << endl;
	}
	return make_tuple(a, e);
}

// Test impulse encounter
tuple<long double, long double, long double> testImpulseEncounter(long double m1, long double m2, long double M_p, long double a, long double e, long double b, long double v)
{
	//Star masses
	array<long double, 2> m = {m1, m2};
	//Open binary
	array<array<long double, 3>, 4> X = setupRandomBinary(a, e, m1, m2);
	//Initial relative velocity of stars
	//array<long double, 3> v_initial;
	//for (int i=0; i<3; i++){
		//v_initial[i] = X[2][i] - X[3][i];
	//}
	//Find impact parameter and velocity vectors
	tuple<array<long double,3>, array<long double,3>> bvvectors = impactAndVelocityVectors(b, v);
	array<long double,3> b_vec = get<0>(bvvectors);
	array<long double,3> v_vec = get<1>(bvvectors);
	//Declare variables
	long double b_90, v_perp, v_para;
	array<long double,3> b_star;
	array<long double,2> b_star_norm;
	//array<array<long double, 3>, 2> delta_v_i; 
	for (int i=0; i<2; ++i){
		//90 degree deflection radius
		b_90 = G*(M_p + m[i])/(v*v);
		//Calculate impact parameter for this star
		b_star = calcBStar(X[i], v_vec, v, b_vec);
		//Calculate norm of b_star
		b_star_norm[i] = norm(b_star);
		//Calculate speed change in b_star direction
		//v_perp = 2.0*M_p*v/(m[i]+M_p) * (b_star_norm[i]/b_90)/(1.0 + b_star_norm[i]*b_star_norm[i]/(b_90*b_90));
		v_perp = 2.0*G*M_p/(v*b_star_norm[i]);
		//Calculate speed change in -v_vec direction
		//v_para = 2.0*M_p*v/(m[i]+M_p) * 1.0/(1.0 + b_star_norm[i]*b_star_norm[i]/(b_90*b_90));
		v_para = 0.0;
		//Change star velocity
		for (int j=0; j<3; ++j){
			//delta_v_i[i][j] = v_perp * b_star[j]/b_star_norm[i] - v_para * v_vec[j]/v;
			X[i+2][j] += v_perp * b_star[j]/b_star_norm[i] - v_para * v_vec[j]/v;
		}
	}
	//Relative velocity change
	//array<long double,3> delta_v;
	//for(int i=0; i<3; i++){
		//delta_v[i] = delta_v_i[0][i] - delta_v_i[1][i];
	//}
	//Initial energy
	long double E_ini = -G*m1*m2/(2.0*a);
	//Close binary
	tuple<long double, long double, bool> result = orbitalElements(X, m1, m2);
	//Final energy
	long double E_fin = -G*m1*m2/(2.0*get<0>(result));
	//Effective impact parameter
	long double b_star_norm_min = min(b_star_norm[0], b_star_norm[1]);
	//Print accuracy testing terms
	//cout << setprecision(16) << "v dot delta v term = " << m1*m2/(m1+m2)*dot(v_initial, delta_v) << endl;
	//cout << setprecision(16) << "delta v squared term = " << 0.5*m1*m2/(m1+m2)*dot(delta_v, delta_v) << endl;
	//cout << setprecision(16) << "two above added = " << m1*m2/(m1+m2)*dot(v_initial, delta_v) + 0.5*m1*m2/(m1+m2)*dot(delta_v, delta_v) << endl;
	//cout << setprecision(16) << "delta E = " << E_fin - E_ini << endl;
	//cout << setprecision(16) << "difference between energy changes = " << m1*m2/(m1+m2)*dot(v_initial, delta_v) + 0.5*m1*m2/(m1+m2)*dot(delta_v, delta_v) - (E_fin - E_ini) << endl;
	//cout << "b_star_norm = " << b_star_norm[0] << " , " << b_star_norm[1] << endl;
	//cout << "b_star_norm_min = " << b_star_norm_min << endl;
	return make_tuple(E_ini, E_fin, b_star_norm_min);
}