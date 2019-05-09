// Set of functions required to simulate encounters between binary stars and a sea of perturbers
#include <cmath>
#include <array>
#include <tuple>
#include <vector>
#include <algorithm>

#include <iostream>
#include <iomanip>
#include <ctime>
#include <fstream>

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
	return 0.5*pow((64.0*G*M_p*M_p*pow(a,3.0)/((m1+m2)*v_rel*v_rel*delta*delta)), 0.25);
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

tuple<long double, long double, long double, bool> impulseEncounterIonised(long double m1, long double m2, long double M_p, long double a, long double e, long double b, long double v, long double E, bool notBound)
{
	//Star masses
	array<long double, 2> m = {m1, m2};
	//Open binary
	array<array<long double, 3>, 4> X = setupRandomBinaryIonised(a, e, m1, m2, E, notBound);

	//array<long double,3> V;
	//for(int i=0; i<3; ++i){
	//	V[i] = X[2][i] - X[3][i];
	//}
	//array<long double,3> dV;

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
		//cout << "Velocity change = " << v_perp * b_star[0]/b_star_norm - v_para * v_vec[0]/v << ", " << v_perp * b_star[1]/b_star_norm - v_para * v_vec[1]/v << ", " << v_perp * b_star[2]/b_star_norm - v_para * v_vec[2]/v << endl;
		for (int j=0; j<3; ++j){
			X[i+2][j] += v_perp * b_star[j]/b_star_norm - v_para * v_vec[j]/v;
		}
	}

	//for(int i=0; i<3; ++i){
	//	dV[i] = X[2][i] - X[3][i] - V[i];
	//}

	//cout << "Energy change = " << m1*m2/(m1+m2)*(dot(V,dV) + 0.5*dot(dV, dV)) << endl;

	//Close binary
	return orbitalElementsIonised(X, m1, m2);
}

//Tested
//Draw an impact parameter from a distribution linear in b up to b_max
long double drawB(long double b_max, long double b_min=0.0)
{
	return b_max*sqrt(randomUniformDoubleClosed(b_min*b_min/(b_max*b_max), 1.0));
}

tuple<vector<long double>, vector<long double>, int, int, int, int, int> MCEncounters(long double v_rel, long double n_p, long double T, long double m1, long double m2, long double M_p, vector<long double> a, vector<long double> e)
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
	//Number of binaries broken
	int N_broken = 0;
	//double t_start;

	int N_encounters = 0;
	int N_encounters_close = 0;
	int N_encounters_far = 0;
	int N_encounters_mid = 0;
	//long double a_0 = a[0];
	//long double e_0 = e[0];
	//b_max = calcBMax(M_p, v_rel, a[0], m1, m2);
	//rate = encounterRate(n_p, v_rel, b_min, b_max, v_min, v_max);
	//long double N_mean = T*rate;
	//cout << "Mean number of encounters for constant a = " << N_mean << endl;

	//Iterate over binaries
	for (int i=0; i<N_bin; ++i){
		//cout << endl;
		cout << "Binary " << i+1 << " of " << N_bin << endl;
		//Start time
		//clock_t t_start;
		//t_start = clock();
		//Time passed
		t = 0.0;

		//ofstream myfile;
		//string filename = "binary_history_" + to_string(i) + ".csv";
		//myfile.open(filename);

		//Implement encounters
		while (true){
			//Maximum impact parameter
			b_max = calcBMax(M_p, v_rel, a[i], m1, m2);
			//b_max = BHTBMax(M_p, v_rel, a_0, m1, m2, e_0);
			//Encounter rate
			rate = encounterRate(n_p, v_rel, b_min, b_max, v_min, v_max);
			//Increment time passed
			t += randomExponential(rate);
			if (t>T){
				break;
			}
			//Draw velocity from distribution
			v = drawVMaxwellian(v_rel, v_max);
			//v = drawMaxwellian(v_rel);
			//Draw impact parameter from distribution
			b = drawB(b_max);
			//Encounter
			if (b<a[i]){
				N_encounters_close += 1;
			}
			if (b>a[i]){
				N_encounters_far += 1;
			}
			if ((b>0.01*a[i]) && (b<100.0*a[i])){
				N_encounters_mid += 1;
			}
			N_encounters += 1;
			result = impulseEncounter(m1, m2, M_p, a[i], e[i], b, v);
			a[i] = get<0>(result);
			e[i] = get<1>(result);
			notBound = get<2>(result);

			/*
			cout << "t, Gyr = " << t*time_scale/(giga*year) << endl;
			cout << "v_rel, km/s = " << v*length_scale/time_scale/1000.0 << endl;
			cout << "b, pc = " << b*length_scale/parsec << endl;
			cout << "a_fin, pc = " << a[i]*length_scale/parsec << endl;
			cout << "e_fin = " << e[i] << endl;
			cout << "N_enc = " << N_encounters << endl;
			cout << "Binary broken? = " << (notBound or (a[i]>=a_T)) << endl;
			cout << endl;
			cin.ignore();
			*/
			//myfile << setprecision(16) << t*time_scale << " ," << v*length_scale/time_scale << " , " << b*length_scale << " , " << a[i]*length_scale << " , " << e[i] << " , " << N_encounters << endl;

			if(notBound or (a[i]>=a_T)){
				N_broken += 1;
				a[i] = -1.0;
				e[i] = -1.0;
				break;
			}
		}
		//myfile.close();
		//Print how long it took
		//cout << "Time taken = " << (clock() - t_start)/(double)(CLOCKS_PER_SEC) << " s" << endl;
	}
	return make_tuple(a, e, N_broken, N_encounters, N_encounters_close, N_encounters_far, N_encounters_mid);
}

// Test impulse encounter
tuple<long double, long double, long double, long double, long double, array<long double,3>, array<long double,3>, long double> testImpulseEncounter(long double m1, long double m2, long double M_p, long double a, long double e, long double b, long double v)
{
	//Star masses
	array<long double, 2> m = {m1, m2};
	//Open binary
	array<array<long double, 3>, 4> X = setupRandomBinary(a, e, m1, m2);
	//Initial relative velocity of stars
	array<long double, 3> v_initial;
	for (int i=0; i<3; i++){
		v_initial[i] = X[2][i] - X[3][i];
	}
	//Find impact parameter and velocity vectors
	tuple<array<long double,3>, array<long double,3>> bvvectors = impactAndVelocityVectors(b, v);
	array<long double,3> b_vec = get<0>(bvvectors);
	array<long double,3> v_vec = get<1>(bvvectors);
	//Declare variables
	long double b_90, v_perp, v_para;
	array<long double,3> b_star;
	array<long double,2> b_star_norm;
	array<array<long double, 3>, 2> delta_v_i; 
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
			delta_v_i[i][j] = v_perp * b_star[j]/b_star_norm[i] - v_para * v_vec[j]/v;
			X[i+2][j] += v_perp * b_star[j]/b_star_norm[i] - v_para * v_vec[j]/v;
		}
	}
	//Relative velocity change
	array<long double,3> delta_v;
	for(int i=0; i<3; i++){
		delta_v[i] = delta_v_i[0][i] - delta_v_i[1][i];
	}
	//Initial energy
	long double E_ini = -G*m1*m2/(2.0*a);
	//Close binary
	tuple<long double, long double, bool> result = orbitalElements(X, m1, m2);
	//Final energy
	long double E_fin = -G*m1*m2/(2.0*get<0>(result));
	//Effective impact parameter
	long double b_star_norm_min = min(b_star_norm[0], b_star_norm[1]);
	//v dot dv term
	//long double v_dv;
	long double v_dv = m1*m2/(m1+m2)*dot(v_initial, delta_v);
	//dv^2 term
	//long double dv_dv;
	long double dv_dv = 0.5*m1*m2/(m1+m2)*dot(delta_v, delta_v);

	array<long double,3> r;
	for (int i=0; i<3; i++){
		r[i] = X[0][i] - X[1][i];
	}
	long double phi = acos(dot(r, v_vec)/(norm(r)*norm(v_vec)));
	long double theta = acos(dot(v_initial, delta_v)/(norm(v_initial)*norm(delta_v)));
	//Print accuracy testing terms
	//cout << setprecision(16) << "v dot delta v term = " << m1*m2/(m1+m2)*dot(v_initial, delta_v) << endl;
	//cout << setprecision(16) << "delta v squared term = " << 0.5*m1*m2/(m1+m2)*dot(delta_v, delta_v) << endl;
	//cout << setprecision(16) << "two above added = " << m1*m2/(m1+m2)*dot(v_initial, delta_v) + 0.5*m1*m2/(m1+m2)*dot(delta_v, delta_v) << endl;
	//cout << setprecision(16) << "delta E = " << E_fin - E_ini << endl;
	//cout << setprecision(16) << "difference between energy changes = " << m1*m2/(m1+m2)*dot(v_initial, delta_v) + 0.5*m1*m2/(m1+m2)*dot(delta_v, delta_v) - (E_fin - E_ini) << endl;
	//cout << "b_star_norm = " << b_star_norm[0] << " , " << b_star_norm[1] << endl;
	//cout << "b_star_norm_min = " << b_star_norm_min << endl;
	return make_tuple(E_ini, E_fin, b_star_norm_min, v_dv, dv_dv, v_initial, delta_v, theta);
}

tuple<vector<long double>, vector<long double>, int, int, int, int, int> MCEncountersNClosest(int N_closest, long double v_rel, long double n_p, long double T, long double m1, long double m2, long double M_p, vector<long double> a, vector<long double> e)
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
	long double b_max, rate, v;
	tuple<long double, long double, bool> result;
	bool notBound;
	//Number of binaries broken
	int N_broken = 0;
	//double t_start;

	int N_encounters = 0;
	int N_encounters_close = 0;
	int N_encounters_far = 0;
	int N_encounters_mid = 0;

	int N_enc;
	vector<long double> bs;
	vector<long double> bs_sorted;
	vector<long double> bs_closest;
	long double b_limit;
	long double a_initial;

	//string filename = "N_enc_broken_dist_1000Msol_a10e4au.csv";
	//ofstream myfile;
	//myfile.open(filename);

	//Iterate over binaries
	for (int i=0; i<N_bin; ++i){
		cout << "Binary " << i+1 << " of " << N_bin << endl;
		N_encounters = 0;
		a_initial = a[i];

		N_enc = 0;
		bs.resize(0);
		bs_sorted.resize(0);
		bs_closest.resize(0);
		//b_max = YCGBMax(a[i], M_p, n_p, v_rel, T);
		b_max = calcBMax(M_p, v_rel, a[i], m1, m2);
		rate = encounterRate(n_p, v_rel, b_min, b_max, v_min, v_max);
		N_enc = randomPoisson(rate*T);

		/*
		if(M_p < 700.0*msol/mass_scale){
			N_closest = N_enc;
		}

		//cout << "Total number of encounters = " << N_enc << endl;

		bs.resize(N_enc);
		bs_sorted.resize(N_enc);
		for (int j=0; j<N_enc; ++j){
			bs[j] = drawB(b_max);
			bs_sorted[j] = bs[j];
		}
		sort(bs_sorted.begin(), bs_sorted.end());

		if (N_enc > N_closest){			
			b_limit = bs_sorted[N_closest];
			for (int j=0; j<N_enc; ++j){
				if (bs[j] < b_limit){
					bs_closest.push_back(bs[j]);
				}
			}
		} else {
			bs_closest.resize(N_enc);
			for (int j=0; j<N_enc; ++j){
				bs_closest[j] = bs[j];
			}
		}
		*/
		bs_closest.resize(N_enc);
		for (int j=0; j<N_enc; ++j){
			bs_closest[j] = drawB(b_max);
		}


		//cout << "Number of encounters = " << static_cast<int>(bs_closest.size()) << endl;

		//Implement encounters
		for (int j=0; j<static_cast<int>(bs_closest.size()); j++){

			//cout << "Encounter " << j+1 << " of " << static_cast<int>(bs_closest.size()) << endl;

			//Draw velocity from distribution
			v = drawVMaxwellian(v_rel, v_max);
			//v = drawMaxwellian(v_rel);

			if (bs_closest[j]<a[i]){
				N_encounters_close += 1;
			}
			if (bs_closest[j]>a[i]){
				N_encounters_far += 1;
			}
			if ((bs_closest[j]>0.01*a[i]) && (bs_closest[j]<100.0*a[i])){
				N_encounters_mid += 1;
			}
			N_encounters += 1;
			result = impulseEncounter(m1, m2, M_p, a[i], e[i], bs_closest[j], v);
			a[i] = get<0>(result);
			e[i] = get<1>(result);
			notBound = get<2>(result);

			/*
			cout << setprecision(16) << "b, pc = " << bs_closest[j]*length_scale/parsec << endl;
			cout << "a_fin, pc = " << a[i]*length_scale/parsec << endl;
			cout << "e_fin = " << e[i] << endl;
			cout << "N_enc = " << N_encounters << endl;
			cout << "Binary broken? = " << (notBound or (a[i]>=a_T)) << endl;
			cout << endl;
			*/
			//myfile << setprecision(16) << t*time_scale << " ," << v*length_scale/time_scale << " , " << b*length_scale << " , " << a[i]*length_scale << " , " << e[i] << " , " << N_encounters << endl;

			if(notBound or (a[i]>=a_T)){
				N_broken += 1;
				a[i] = -1.0;
				e[i] = -1.0;
				//myfile << setprecision(16) << N_encounters << ", " << a_initial*length_scale << endl;
				break;
			}
		}
	}
	//myfile.close();
	return make_tuple(a, e, N_broken, N_encounters, N_encounters_close, N_encounters_far, N_encounters_mid);
}



tuple<vector<long double>, vector<long double>, int, int, int, int, int> MCEncountersIonised(long double v_rel, long double n_p, long double T, long double m1, long double m2, long double M_p, vector<long double> a, vector<long double> e)
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
	long double b_max, rate, v;
	tuple<long double, long double, long double, bool> result;
	bool notBound;
	//Number of binaries broken
	int N_broken = 0;
	//double t_start;

	int N_encounters = 0;
	int N_encounters_close = 0;
	int N_encounters_far = 0;
	int N_encounters_mid = 0;
	int N_enc;
	vector<long double> bs;
	long double E, M, n;
	bool hasBroken;
	bool rebound;
	int N_rebound = 0;
	int N_close = 0;

	//Iterate over binaries
	for (int i=0; i<N_bin; ++i){
		cout << '\r' << "Binary " << i+1 << " of " << N_bin;

		N_encounters = 0;
		N_enc = 0;
		b_max = calcBMax(M_p, v_rel, a[i], m1, m2);
		rate = encounterRate(n_p, v_rel, b_min, b_max, v_min, v_max);
		N_enc = randomPoisson(rate*T);
		bs.resize(N_enc);
		for (int j=0; j<N_enc; ++j){
			bs[j] = drawB(b_max);
		}

		hasBroken = false;
		rebound = false;
		notBound = false;

		//cout << "a_0 = " << a[i]*length_scale/au << " , N_enc = " << N_enc << endl;

		//Implement encounters
		for (int j=0; j<N_enc; j++){
			//cout << '\r' << "Encounter " << j+1 << " of " << N_enc;

			//Draw velocity from distribution
			v = drawVMaxwellian(v_rel, v_max);
			//v = drawMaxwellian(v_rel);

			if (notBound){
				//Evolve E in time
				//Mean motion
				n = sqrt(G*(m1+m2)/(pow(a[i],3)));
				//Mean anomaly
				M = e[i]*sinh(E) - E;
				M += n*randomExponential(rate);
				M = fmod(M, 2.0*pi);
				//New eccentric anomaly
				E = eccentricAnomalyIonised(e[i], M, notBound);
			} else {
				if (hasBroken){
					rebound = true;
				}
			}

			if (bs[j]<a[i]){
				N_encounters_close += 1;
			}
			if (bs[j]>a[i]){
				N_encounters_far += 1;
			}
			if ((bs[j]>0.01*a[i]) && (bs[j]<100.0*a[i])){
				N_encounters_mid += 1;
			}
			N_encounters += 1;
			result = impulseEncounterIonised(m1, m2, M_p, a[i], e[i], bs[j], v, E, notBound);
			a[i] = get<0>(result);
			e[i] = get<1>(result);
			E = get<2>(result);
			notBound = get<3>(result);
			if (notBound){
				hasBroken = true;
			}

			//cout << "Not bound? " << notBound << ", a/au = " << a[i]*length_scale/au << ", e = " << e[i] << ", E = " << E << endl;

			if(a[i]>=a_T){
				hasBroken = true;
				//cout << "Binary broken!" << endl;
				N_broken += 1;
				a[i] = -1.0;
				e[i] = -1.0;
				break;
			}
		}
		//if (notBound){
			//cout << "Unbound binary at end time. a/au = " << a[i]*length_scale/au << ", e = " << e[i] << ", E = " << E << ", N_enc = " << N_enc << endl;
		//}
		if (rebound && (notBound == false)) {
			//cout << "Rebound binary!" << endl;
			N_rebound ++;
		}
		if ((a[i]>0.0) && (a[i]<100.0*parsec/length_scale) && notBound){
			N_close ++;
		}
		//cout << endl;
	}
	cout << endl;
	cout << "Number of binaries rebound = " << N_rebound << endl;
	cout << "Number of binaries unbound within 100pc = " << N_close << endl;
	return make_tuple(a, e, N_broken, N_encounters, N_encounters_close, N_encounters_far, N_encounters_mid);
}