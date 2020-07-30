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
#include "nbodyintegration.h"
using namespace std;


//Tested
// Rate of encounters between impact parameters b0 and b1 and between relative speeds of v0 and v1 in a sea of particles with number density n_p and relative velocity dispersion v_rel
long double encounterRate(long double n_p, long double v_rel, long double b0, long double b1, long double v0, long double v1)
{
	return sqrt(2.0L*pi)*n_p/v_rel*(pow(b1, 2.0L)-pow(b0, 2.0L))*((pow(v0, 2.0L)+2.0L*pow(v_rel, 2.0L))*exp(-pow(v0, 2.0L)/(2.0L*pow(v_rel, 2.0L)))-(pow(v1, 2.0L)+2.0L*pow(v_rel, 2.0L))*exp(-pow(v1, 2.0L)/(2.0L*pow(v_rel, 2.0L))));
}

//Tested
// Calculates the impact parameter at which the fractional change in semi-major axis of the binary will be equal to delta, for perturber mass M_p, relative velocity dispersion v_rel, semi-major axis a and binary star masses m1 and m2.
long double calcBMax(long double M_p, long double v_rel, long double a, long double m1, long double m2, long double delta = pow(10.0L, -3.0L))
{
	return pow((64.0L*G*M_p*M_p*pow(a,3.0L)/((m1+m2)*v_rel*v_rel*delta*delta)), 0.25L);
}

long double BHTBMax(long double M_p, long double v_rel, long double a, long double m1, long double m2, long double e)
{
	long double C = 8.0;
	long double D = 0.6*(1.0 + e);
	long double v_c = G*m1*m2*(m1 + m2 + M_p)/(M_p*(m1 + m2)*a);
	return (C*v_c/v_rel + D)/a;
}

long double YCGBMax(long double a, long double M_p, long double n_p, long double v_rel, long double T)
{
	long double b_min = sqrt(1.0/(pi*n_p*v_rel*T));
	return max(10.0*b_min, 2.0*a);
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
long double drawB(long double b_max, long double b_min=0.0L)
{
	return b_max*sqrt(randomUniformDoubleClosed(b_min*b_min/(b_max*b_max), 1.0L));
}

tuple<long double, long double, long double, bool, long double> impulseEncounterIonised(long double m1, long double m2, long double M_p, long double a, long double e, long double b, long double v, long double E, long double r, bool notBound, bool &non_converged_binary, bool linear)
{
	//Star masses
	array<long double, 2> m = {m1, m2};
	//Open binary
	array<array<long double, 3>, 4> X = setupRandomBinaryIonised(a, e, m1, m2, E, notBound, non_converged_binary);

	//Find impact parameter and velocity vectors
	tuple<array<long double,3>, array<long double,3>> bvvectors = impactAndVelocityVectors(b, v);
	array<long double,3> b_vec = get<0>(bvvectors);
	//cout << "b = " << b_vec[0] << ", " << b_vec[1] << ", " << b_vec[2] << endl;
	array<long double,3> v_vec = get<1>(bvvectors);
	//cout << "v = " << v_vec[0] << ", " << v_vec[1] << ", " << v_vec[2] << endl;
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
		//cout << "Velocity change = " << v_perp * b_star[0]/b_star_norm - v_para * v_vec[0]/v << ", " << v_perp * b_star[1]/b_star_norm - v_para * v_vec[1]/v << ", " << v_perp * b_star[2]/b_star_norm - v_para * v_vec[2]/v << endl;
		for (int j=0; j<3; ++j){
			X[i+2][j] += v_perp * b_star[j]/b_star_norm - v_para * v_vec[j]/v;
		}
	}
	//Close binary
	return orbitalElementsIonised(X, m1, m2);
}

tuple<vector<long double>, vector<long double>, vector<long double>, vector<long double>> MCEncountersIonised(long double v_rel, long double n_p, long double T, long double m1, long double m2, long double M_p, vector<long double> a, vector<long double> e)
{
	//Minimum impact parameter
	long double b_min = 0.0L;
	//Minimum relative velocity of encounter
	long double v_min = 0.0L;
	//Maximum relative velocity of encounter
	long double v_max = 100.0L * v_rel;
	//Maximum semimajor axis
	long double a_T =10000.0L * parsec/length_scale;
	//Number of binaries
	int N_bin = static_cast<int>(a.size());
	cout << "Number of binaries = " << N_bin << endl;
	//Declare variable types
	long double b_max, rate, v, b;
	tuple<long double, long double, long double, bool, long double> result;
	bool notBound;
	//Number of binaries broken
	int N_broken = 0;

	vector<long double> r_inis, r;
	r_inis.resize(N_bin);
	r.resize(N_bin);
	r_inis.shrink_to_fit();
	r.shrink_to_fit();

	int N_enc;
	long double E, M, n, f, dt, E_before_evolution;
	bool hasBroken;
	bool rebound;
	int N_rebound = 0;
	int N_close = 0;
	int N_nonconverged = 0;

	vector<long double> ms = {m1, m2};
	bool ini_arrays=true;

	vector<array<long double, 3>> X_nbody;

	tuple<array<long double,3>, array<long double,3>> bvvectors;
	array<long double,3> b_vec;
	array<long double,3> v_vec;
	long double b_90, b_star_norm, v_perp, v_para;
	array<long double,3> b_star;
	array<array<long double, 3>, 4> X;

	bool non_converged_binary;

	//Iterate over binaries
	for (int i=0; i<N_bin; ++i){
		cout << '\r' << "Binary " << i+1 << " of " << N_bin << flush;

		hasBroken = false;
		rebound = false;
		notBound = false;
		non_converged_binary = false;

		b_max = calcBMax(M_p, v_rel, a[i], m1, m2);
		rate = encounterRate(n_p, v_rel, b_min, b_max, v_min, v_max);
		N_enc = randomPoisson(rate*T);

		// Randomise mean anomaly
		M = randomUniformDoubleOpen(0.0L, 2.0L*pi);
		// Find eccentric anomaly
		E = eccentricAnomaly(e[i], M, non_converged_binary);

		//Find initial separation
		//True anomaly
		f = 2.0*atan(sqrt((1.0+e[i])/(1.0-e[i]))*tan(E/2.0));
		// Initial separation
		r_inis[i] = a[i]*(1.0 - e[i]*e[i])/(1.0 + e[i]*cos(f));

		if (N_enc == 0){
			// Randomise mean anomaly
			M = randomUniformDoubleOpen(0.0L, 2.0L*pi);
			// Find eccentric anomaly
			E = eccentricAnomaly(e[i], M, non_converged_binary);
			//True anomaly
			f = 2.0*atan(sqrt((1.0+e[i])/(1.0-e[i]))*tan(E/2.0));
			// Final separation
			r[i] = a[i]*(1.0 - e[i]*e[i])/(1.0 + e[i]*cos(f));
		}

		//Implement encounters
		for (int j=0; j<N_enc; j++){

			//Draw impact parameter
			b = drawB(b_max);
			//Draw relative encounter velocity
			v = drawVMaxwellian(v_rel, v_max);
			//Generate impact parameter and PBH velocity vectors
			bvvectors = impactAndVelocityVectors(b, v);
			b_vec = get<0>(bvvectors);
			v_vec = get<1>(bvvectors);
	
			//Evolve E in time
			//Mean motion
			n = sqrt(G*(m1+m2)/(pow(a[i],3)));
			//Time to evolve for
			dt = randomExponential(rate);
			//Calculate mean anomaly
			if (e[i] < 1.0L){
				M = E - e[i]*sin(E);
			} else {
				M = e[i]*sinh(E) - E;
			}
			//Evolve mean anomaly forward in time
			M += n*dt;
			if (e[i] < 1){
				M = fmod(M, 2.0L*pi);
			}
			E_before_evolution = E;
			//Calculate new eccentric anomaly
			E = eccentricAnomalyIonised(e[i], M, notBound, non_converged_binary);

			//If calculating eccentric anomaly doesn't work use n body integration to evolve instead
			if (non_converged_binary){
				E = E_before_evolution;
				n = sqrt(G*(m1+m2)/(pow(a[i],3.0L)));
				if (notBound){
					X_nbody = { {
						{0.0L, 0.0L, 0.0L},
						{a[i]*(cosh(E) - e[i]), a[i]*sqrt(e[i]*e[i]-1.0L)*sinh(E), 0.0L},
						{0.0L, 0.0L, 0.0L}, 
						{n*a[i]*sinh(E)/(e[i]*cosh(E) - 1.0L), n*a[i]*sqrt(e[i]*e[i]-1.0L)*cosh(E)/(e[i]*cosh(E) - 1.0L), 0.0L}} };
				} else {
					X_nbody = { {
						{0.0L, 0.0L, 0.0L},
						{a[i]*(cos(E)-e[i]), a[i]*sqrt(1.0L-e[i]*e[i])*sin(E), 0.0L},
						{0.0L, 0.0L, 0.0L}, 
						{-n*a[i]/(1.0L-e[i]*cos(E))*sin(E), n*a[i]/(1.0L-e[i]*cos(E))*sqrt(1.0L-e[i]*e[i])*cos(E), 0.0L}} };
				}
				X_nbody = evolve(2, ms, X_nbody, dt, ini_arrays = ini_arrays);
				result = orbitalElementsIonised(X_nbody, m1, m2);
				a[i]= get<0>(result);
				e[i] = get<1>(result);
				E = get<2>(result);
				notBound = get<3>(result);
				r[i] = get<4>(result);
			}
			X = setupRandomBinaryIonised(a[i], e[i], m1, m2, E, notBound, non_converged_binary);
			result = orbitalElementsIonised(X, m1, m2);
			a[i] = get<0>(result);
			e[i] = get<1>(result);
			E = get<2>(result);
			notBound = get<3>(result);
			r[i] = get<4>(result);
			X = setupRandomBinaryIonised(a[i], e[i], m1, m2, E, notBound, non_converged_binary);
			//Implement encounter
			for (int k=0; k<2; ++k){
				//90 degree deflection radius
				b_90 = G*(M_p + ms[k])/(v*v);
				//Calculate impact parameter for this star
				b_star = calcBStar(X[k], v_vec, v, b_vec);
				//Calculate norm of b_star
				b_star_norm = norm(b_star);
				//Calculate speed change in b_star direction
				v_perp = 2.0L*M_p*v/(ms[k]+M_p) * (b_star_norm/b_90)/(1.0L + b_star_norm*b_star_norm/(b_90*b_90));
				//Calculate speed change in -v_vec direction
				v_para = 2.0L*M_p*v/(ms[k]+M_p) * 1.0L/(1.0L + b_star_norm*b_star_norm/(b_90*b_90));
				//Change star velocity
				for (int l=0; l<3; ++l){
					X[k+2][l] += v_perp * b_star[l]/b_star_norm - v_para * v_vec[l]/v;
				}
			}
			result = orbitalElementsIonised(X, m1, m2);
			a[i] = get<0>(result);
			e[i] = get<1>(result);
			E = get<2>(result);
			notBound = get<3>(result);
			r[i] = get<4>(result);
			//Check if separation exceeds cut off radius
			if(r[i]>=a_T){
				a[i] = -1.0L;
				e[i] = -1.0L;
				r[i] = -1.0L;
				notBound = true;
				break;
			}
			//Check if rebound
			if (notBound){
				hasBroken = true;
			} else {
				if (hasBroken){
					rebound = true;
					hasBroken = false;
				}
			}
		}

		
		if (rebound && (notBound == false)) {
			N_rebound ++;
		}
		if ((a[i]>0.0L) && (a[i]<100.0L*parsec/length_scale) && notBound){
			N_close ++;
		}
		if (notBound){
			N_broken ++;
		}
	}
	cout << endl;
	cout << "Number of binaries rebound = " << N_rebound << endl;
	cout << "Number of binaries unbound within 100pc = " << N_close << endl;
	cout << "Number of binaries not converged = " << N_nonconverged << endl;
	cout << "Number of binaries broken = " << N_broken << endl;
	return make_tuple(a, e, r_inis, r);
}