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
#include "nbodyintegration.h"
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

tuple<vector<long double>, vector<long double>> MCEncountersNClosest(int N_closest, long double v_rel, long double n_p, long double T, long double m1, long double m2, long double M_p, vector<long double> a, vector<long double> e)
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


	int N_enc;
	vector<long double> bs;
	vector<long double> bs_sorted;
	vector<long double> bs_closest;
	long double b_limit;


	//Iterate over binaries
	for (int i=0; i<N_bin; ++i){
		cout << '\r' << "Binary " << i+1 << " of " << N_bin << flush;

		//b_max = YCGBMax(a[i], M_p, n_p, v_rel, T);
		b_max = calcBMax(M_p, v_rel, a[i], m1, m2);
		rate = encounterRate(n_p, v_rel, b_min, b_max, v_min, v_max);
		N_enc = randomPoisson(rate*T);
		

		if (N_closest == 0){
			bs_closest.resize(N_enc);
			for (int j=0; j<N_enc; ++j){
				bs_closest[j] = drawB(b_max);
			}
		} else {
			bs.resize(N_enc);
			bs_sorted.resize(N_enc);
			for (int j=0; j<N_enc; ++j){
				bs[j] = drawB(b_max);
				bs_sorted[j] = bs[j];
			}
			sort(bs_sorted.begin(), bs_sorted.end());

			if (N_enc > N_closest){			
				b_limit = bs_sorted[N_closest];
				bs_closest.resize(0);
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

		}


		//cout << "Number of encounters = " << static_cast<int>(bs_closest.size()) << endl;

		//Implement encounters
		for (int j=0; j<static_cast<int>(bs_closest.size()); j++){

			//cout << "Encounter " << j+1 << " of " << static_cast<int>(bs_closest.size()) << endl;

			//Draw velocity from distribution
			v = drawVMaxwellian(v_rel, v_max);
			//v = drawMaxwellian(v_rel);

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
				a[i] = -1.0;
				e[i] = -1.0;
				break;
			}
		}
	}
	return make_tuple(a, e);
}

tuple<long double, long double, long double, bool, long double> impulseEncounterIonised(long double m1, long double m2, long double M_p, long double a, long double e, long double b, long double v, long double E, long double r, bool notBound, bool &non_converged_binary, bool linear)
{
	//Star masses
	array<long double, 2> m = {m1, m2};
	//Open binary
	array<array<long double, 3>, 4> X = setupRandomBinaryIonised(a, e, m1, m2, E, r, notBound, non_converged_binary, linear);

	//cout << "a = " << a << endl;
	//cout << "e = " << e << endl;

	//cout << "X = " << X[0][0] << ", " << X[0][1] << ", " << X[0][2] << endl;
	//cout << X[1][0] << ", " << X[1][1] << ", " << X[1][2] << endl;
	//cout << X[2][0] << ", " << X[2][1] << ", " << X[2][2] << endl;
	//cout << X[3][0] << ", " << X[3][1] << ", " << X[3][2] << endl;



	//array<long double,3> V;
	//for(int i=0; i<3; ++i){
	//	V[i] = X[2][i] - X[3][i];
	//}
	//array<long double,3> dV;

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
	//cout << "ENCOUNTER" << endl;

	//for(int i=0; i<3; ++i){
	//	dV[i] = X[2][i] - X[3][i] - V[i];
	//}
	//cout << "Energy change = " << m1*m2/(m1+m2)*(dot(V,dV) + 0.5*dot(dV, dV)) << endl;

	//Close binary
	return orbitalElementsIonised(X, m1, m2);
}

vector<array<long double, 3>> impulseEncounterXV(vector<array<long double, 3>> X, long double M_p, long double m1, long double m2, long double b, long double v)
{
	//Star masses
	array<long double, 2> m = {m1, m2};

	// Centre of mass position vector
	array<long double, 3> R;
	// Centre of mass velocity vector
	array<long double,3> V;
	for (int i=0; i<3; ++i){
		R[i] = (m1*X[0][i] + m2*X[1][i])/(m1 + m2);
		V[i] = (m1*X[2][i] + m2*X[3][i])/(m1 + m2);
	}
	// Move into centre of mass rest frame
	for (int i=0; i<3; ++i){
		X[0][i] -= R[i];
		X[1][i] -= R[i];
		X[2][i] -= V[i];
		X[3][i] -= V[i];
	}

	//cout << "X = " << X[0][0] << ", " << X[0][1] << ", " << X[0][2] << endl;
	//cout << X[1][0] << ", " << X[1][1] << ", " << X[1][2] << endl;
	//cout << X[2][0] << ", " << X[2][1] << ", " << X[2][2] << endl;
	//cout << X[3][0] << ", " << X[3][1] << ", " << X[3][2] << endl;

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
	return X;
}

tuple<vector<long double>, vector<long double>> MCEncountersIonised(long double v_rel, long double n_p, long double T, long double m1, long double m2, long double M_p, vector<long double> a, vector<long double> e)
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
	cout << "Number of binaries = " << N_bin << endl;
	//Declare variable types
	long double b_max, rate, v, b;
	tuple<long double, long double, long double, bool, long double> result;
	bool notBound;
	//Number of binaries broken
	int N_broken = 0;
	//double t_start;


	int N_enc;
	long double E, M, n, En, t, r;
	bool hasBroken;
	bool rebound;
	int N_rebound = 0;
	int N_close = 0;
	int N_nonconverged = 0;
	bool linear;
	vector<array<long double, 3>> X;

	vector<long double> ms = {m1, m2};
	bool ini_arrays=true;
	array<long double,3> R;
	tuple<long double, long double, bool> result2;
	long double a_break;

	long double r_rebound_max = 0.0;
	long double r_rebound_max_total = 0.0;
	long double r_nonconverged_min = a_T;

	//long double En_previous;
	//long double En;
	//bool notBound_previous;
	//ofstream myfile;
	//myfile.open("dE_break_binary_1000Msol.csv");

	bool non_converged_binary;

	//Iterate over binaries
	for (int i=0; i<N_bin; ++i){
		cout << '\r' << "Binary " << i+1 << " of " << N_bin << flush;

		hasBroken = false;
		rebound = false;
		notBound = false;
		non_converged_binary = false;
		linear = false;
		r_rebound_max = 0.0;

		b_max = calcBMax(M_p, v_rel, a[i], m1, m2);
		rate = encounterRate(n_p, v_rel, b_min, b_max, v_min, v_max);
		N_enc = randomPoisson(rate*T);

		//En_previous = -G*m1*m2/(2.0*a[i]);

		//cout << "a_0 = " << a[i]*length_scale/au << " , N_enc = " << N_enc << endl;

		// Randomise mean anomaly
		M = randomUniformDoubleOpen(0.0, 2.0*pi);
		//cout << "M_0 = " << M << endl;
		// Find eccentric anomaly
		E = eccentricAnomaly(e[i], M, non_converged_binary);
		//cout << "E_0 = " << E << endl;
		r = a[i];

		//Implement encounters
		for (int j=0; j<N_enc; j++){
			//cout << '\r' << "Encounter " << j+1 << " of " << N_enc;
			if(r>=a_T){
				//cout << "Binary broken!" << endl;
				a[i] = -1.0;
				e[i] = -1.0;
				notBound = true;
				break;
			}

			//Draw impact parameter
			b = drawB(b_max);
			//cout << "b = " << b*length_scale << endl;
			//Draw relative encounter velocity
			v = drawVMaxwellian(v_rel, v_max);
			//cout << "v = " << v*length_scale/time_scale << endl;

			//Evolve E in time
			//Mean motion
			n = sqrt(G*(m1+m2)/(pow(a[i],3)));
			//cout << "n = " << n << endl;

			t = randomExponential(rate);

			if (linear){
				r += sqrt(G*(m1+m2)/a[i])* randomExponential(rate);
			} else {
				if (e[i] < 1){
					//Mean anomaly
					M = E - e[i]*sin(E);

					M += n*t;
					if (e[i] < 1){
						M = fmod(M, 2.0*pi);
					}
					E = eccentricAnomalyIonised(e[i], M, notBound, non_converged_binary);

					result = impulseEncounterIonised(m1, m2, M_p, a[i], e[i], b, v, E, r, notBound, non_converged_binary, linear);
					a[i] = get<0>(result);
					e[i] = get<1>(result);
					E = get<2>(result);
					notBound = get<3>(result);
					r = get<4>(result);
					a_break = a[i];
				} else if (e[i] > 1){
					//Hyperbolic equations

					//Mean anomaly
					M = e[i]*sinh(E) - E;
					M += n*t;
					E = eccentricAnomalyIonised(e[i], M, notBound, non_converged_binary);
					result = impulseEncounterIonised(m1, m2, M_p, a[i], e[i], b, v, E, r, notBound, non_converged_binary, linear);
					a[i] = get<0>(result);
					e[i] = get<1>(result);
					E = get<2>(result);
					notBound = get<3>(result);
					r = get<4>(result);
					//cout << "r = " << r*length_scale << endl;


					/*
					//linear
					r += sqrt(G*(m1+m2)/a_break)* randomExponential(rate);
					result = impulseEncounterIonised(m1, m2, M_p, a[i], e[i], b, v, E, r, notBound, non_converged_binary, true);
					a[i] = get<0>(result);
					e[i] = get<1>(result);
					E = get<2>(result);
					notBound = get<3>(result);
					r = get<4>(result);
					*/

					/*
					//nbody integration
					X = setupRandomBinaryVector(a[i], e[i], m1, m2);
					X = evolve(2, ms, X, t, ini_arrays = ini_arrays);
					X = impulseEncounterXV(X, M_p, m1, m2, b, v);

					array<long double,3> R;
					for (int k=0; k<3; ++k){
						R[k] = X[0][k] - X[1][k];
					}
					r = norm(R);
					result2 = orbitalElements(X, m1, m2);
					a[i] = get<0>(result2);
					e[i] = get<1>(result2);
					notBound = get<2>(result2);
					*/

				} else {
					cout << endl << "e = 1" << endl;
					break;
				}
				//cout << "dt = " << t << endl;
				//M += n*t;
				//if (e[i] < 1){
				//	M = fmod(M, 2.0*pi);
				//}
				//E = eccentricAnomalyIonised(e[i], M, notBound, non_converged_binary);
			}
			//cout << "M = " << M << endl;
			//New eccentric anomaly
			//cout << "Separation = " << a[i]*(1.0 - e[i]*cos(E)) << endl;
			//cout << "Ecc pre-encounter = " << E << endl;
			//cout << "M pre-encounter = " << E - e[i]*sin(E) << endl;

			//cout << "Energy = " << En << endl;
			//result = impulseEncounterIonised(m1, m2, M_p, a[i], e[i], b, v, E, r, notBound, non_converged_binary, linear);
			//En_previous = En;
			//notBound_previous = notBound;
			//a[i] = get<0>(result);
			//cout << "a = " << a[i] << endl;
			//e[i] = get<1>(result);
			//cout << "e = " << e[i] << endl;
			//E = get<2>(result);
			
			//notBound = get<3>(result);
			//r = get<4>(result);
			//cout << "Energy = " << En << endl;

			//cout << "Separation = " << a[i]*(1.0 - e[i]*cos(E)) << endl;
			//cout << endl << "Ecc = " << E << endl;
			//cout << "M = " << E - e[i]*sin(E) << endl;

			if (non_converged_binary){
				r_nonconverged_min = min(r_nonconverged_min, r);
				//cout << "Non-converged binary!" << endl;
				//cout << "Separation, au = " << r*length_scale/au << endl;
				notBound = false;
				N_nonconverged ++;
				break;
			}

			if(r>=a_T){
				//cout << "Binary broken!" << endl;
				a[i] = -1.0;
				e[i] = -1.0;
				notBound = true;
				break;
			}
			

			//if (r>2.0*a[i]/(1.0-pow(1.0-pow(10.0, -10.0), 2.0))){
			if ((E>9.21) && (e[i]>20000)){
				linear = false;
			} else {
				linear = false;
			}

			if (notBound){
				//cout << endl << "Binary broken!" << endl;
				//cout << endl;
				hasBroken = true;
			} else {
				if (hasBroken){
					rebound = true;
					r_rebound_max = max(r_rebound_max, r);
					//hasBroken = false;
					//cout << "Rebound binary!" << endl;
					//cout << "Separation, pc = " << r*length_scale/parsec << endl;
					//cout << "r_rebound_max, pc = " << r_rebound_max*length_scale/parsec << endl;
				}
			}

			//if (notBound && !(notBound_previous)){
				//myfile << setprecision(16) << (En - En_previous)/En_previous << endl;
			//}

			//cout << endl << "Not bound? " << notBound << ", a/au = " << a[i]*length_scale/au << ", e = " << e[i] << ", E = " << E << endl;
			//cout << endl;
			//cin.ignore();

		}
		//if (notBound){
			//cout << "Unbound binary at end time. a/au = " << a[i]*length_scale/au << ", e = " << e[i] << ", E = " << E << ", N_enc = " << N_enc << endl;
		//}
		if (rebound && (notBound == false)) {
			//cout << "Rebound binary bound at the end!!!!!!!!!!!!!!!!!!!!!" << endl;
			//cout << endl << "Rebound binary bound at the end!!" << endl;
			r_rebound_max_total = max(r_rebound_max_total, r_rebound_max);
			N_rebound ++;
		}
		if ((a[i]>0.0) && (a[i]<100.0*parsec/length_scale) && notBound){
			N_close ++;
		}
		if (notBound){
			//cout << "Binary considered broken at the end" << endl;
			N_broken ++;
		}
		//cout << endl;
		//cout << "a = " << a[i] << endl;
		//cout << "e = " << e[i] << endl;
		//cout << endl << "Ecc = " << E << endl;
		//cin.ignore();
	}
	//myfile.close();
	cout << endl;
	cout << "Number of binaries rebound = " << N_rebound << endl;
	cout << "Number of binaries unbound within 100pc = " << N_close << endl;
	cout << "Number of binaries not converged = " << N_nonconverged << endl;
	cout << "Number of binaries broken = " << N_broken << endl;
	cout << "Maximum Separation of rebound binaries, pc = " << r_rebound_max_total*length_scale/parsec << endl;
	cout << "Minimum Separation of non-converged binaries, pc = " << r_nonconverged_min*length_scale/parsec << endl;
	return make_tuple(a, e);
}



tuple<vector<long double>, vector<long double>> MCEncountersXV(long double v_rel, long double n_p, long double T, long double m1, long double m2, long double M_p, vector<long double> a, vector<long double> e)
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
	cout << "Number of binaries = " << N_bin << endl;
	//Declare variable types
	long double b_max, rate, v, b, t, Ecc, r_mag, n, r_dot_mag, f;
	int N_enc;
	//Position and velocity vector
	vector<array<long double, 3>> X;
	X.resize(4);
	X.shrink_to_fit();
	tuple<long double, long double, bool> result;
	bool notBound, hasBroken, rebound;
	int N_rebound = 0;

	vector<long double> M = {m1, m2};
	M.shrink_to_fit();

	initialise_arrays(2);
	bool ini_arrays = false;

	array<long double, 3> r, r_dot;

	//Iterate over binaries
	for (int i=0; i<N_bin; ++i){
		cout << '\r' << "Binary " << i+1 << " of " << N_bin << flush;

		//cout << endl << setprecision(16) << "Initial energy = " << -G*m1*m2/(2.0*a[i]) << endl;

		notBound = false;
		hasBroken = false;
		rebound = false;

		b_max = calcBMax(M_p, v_rel, a[i], m1, m2);
		rate = encounterRate(n_p, v_rel, b_min, b_max, v_min, v_max);
		N_enc = randomPoisson(rate*T);
		//Set-up binary
		X = setupRandomBinaryVector(a[i], e[i], m1, m2);

		//Implement encounters
		for (int j=0; j<N_enc; j++){
			//cout << '\r' << "Encounter " << j+1 << " of " << N_enc << flush;
			//cout << setprecision(16) << "x1 = " << X[0][0]*length_scale << ", " << X[0][1]*length_scale << ", " << X[0][2]*length_scale << endl;
			//cout << "x2 = " << X[1][0]*length_scale << ", " << X[1][1]*length_scale << ", " << X[1][2]*length_scale << endl;
			//cout << "v1 = " << X[2][0]*length_scale << ", " << X[2][1]*length_scale << ", " << X[2][2]*length_scale << endl;
			//cout << "v2 = " << X[3][0]*length_scale << ", " << X[3][1]*length_scale << ", " << X[3][2]*length_scale << endl;
			//Draw impact parameter
			b = drawB(b_max);
			//cout << "b = " << b*length_scale << endl;
			//Draw relative encounter velocity
			v = drawVMaxwellian(v_rel, v_max);
			//cout << "v = " << v*length_scale/time_scale << endl;

			//result = orbitalElements(X, m1, m2);
			//cout << "Energy = " << -G*m1*m2/(2.0*get<0>(result)) << endl;
			//n = sqrt(G*(m1+m2)/(pow(a[i],3)));
			//cout << "n = " << n << endl;

			//Evolve binary forward in time
			t = randomExponential(rate);
			//cout << "dt = " << t << endl;
			X = evolve(2, M, X, t, ini_arrays = ini_arrays);

			//CALCULATE Ecc
			// Separation vector
			//r = {X[0][0] - X[1][0], X[0][1] - X[1][1], X[0][2] - X[1][2]};
			//r_dot = {X[2][0] - X[3][0], X[2][1] - X[3][1], X[2][2] - X[3][2]};
			// Magnitudes of the above vectors
			//r_mag = norm(r);
			//r_dot_mag = norm(r_dot);
			//cout << "Separation = " << r_mag << endl;
			//if (notBound){
				//Hyperbolic orbit
			//	Ecc = acosh((r_mag/a[i] + 1.0)/e[i]);
			//} else{
				//Elliptical orbit
			//	Ecc = acos((1.0 - r_mag/a[i])/e[i]);
			//	f = asin(sqrt(1.0-e[i]*e[i])*r_dot_mag/(n*a[i]*e[i]));
			//	if (((0<=Ecc<pi) && (pi<=f<2*pi)) || ((pi<=Ecc<2.0*pi) && (0<=f<pi))){
			//		Ecc = 2.0*pi - Ecc;
			//	} 
			//}
			//cout << "Ecc pre-encounter = " << Ecc << endl;
			//cout << "M pre-encounter = " << Ecc - e[i]*sin(Ecc) << endl;


			X = impulseEncounterXV(X, M_p, m1, m2, b, v);

			//cout << "Encounter" << endl;
			//cout << "x1 = " << X[0][0]*length_scale << ", " << X[0][1]*length_scale << ", " << X[0][2]*length_scale << endl;
			//cout << "x2 = " << X[1][0]*length_scale << ", " << X[1][1]*length_scale << ", " << X[1][2]*length_scale << endl;
			//cout << "v1 = " << X[2][0]*length_scale << ", " << X[2][1]*length_scale << ", " << X[2][2]*length_scale << endl;
			//cout << "v2 = " << X[3][0]*length_scale << ", " << X[3][1]*length_scale << ", " << X[3][2]*length_scale << endl;

			//Separation vector
			array<long double,3> R;
			for (int k=0; k<3; ++k){
				R[k] = X[0][k] - X[1][k];
			}
			if(norm(R)>=a_T){
				a[i] = -1.0;
				e[i] = -1.0;
				notBound = true;
				break;
			}

			//Check if unbound
			result = orbitalElements(X, m1, m2);
			//a[i] = get<0>(result);
			//cout << "a = " << a[i] << endl;
			//e[i] = get<1>(result);
			//cout << "e = " << e[i] << endl;
			//cout << "Energy = " << -G*m1*m2/(2.0*get<0>(result)) << endl;
			notBound = get<2>(result);

			/*
			//CALCULATE Ecc
			// Separation vector
			r = {X[0][0] - X[1][0], X[0][1] - X[1][1], X[0][2] - X[1][2]};
			r_dot = {X[2][0] - X[3][0], X[2][1] - X[3][1], X[2][2] - X[3][2]};
			// Magnitudes of the above vectors
			r_mag = norm(r);
			r_dot_mag = norm(r_dot);
			//cout << "Separation = " << r_mag << endl;
			if (notBound){
				//Hyperbolic orbit
				Ecc = acosh((r_mag/a[i] + 1.0)/e[i]);
			} else{
				//Elliptical orbit
				Ecc = acos((1.0 - r_mag/a[i])/e[i]);
				f = asin(sqrt(1.0-e[i]*e[i])*r_dot_mag/(n*a[i]*e[i]));
				if (((0<=Ecc<pi) && (pi<=f<2*pi)) || ((pi<=Ecc<2.0*pi) && (0<=f<pi))){
					Ecc = 2.0*pi - Ecc;
				} 
			}
			//cout << "Ecc post-encounter = " << Ecc << endl;
			//cout << "M post-encounter = " << Ecc - e[i]*sin(Ecc) << endl;
			*/

			//cout << "notBound = " << notBound << endl;
			if (notBound){
				hasBroken = true;
			} else {
				if (hasBroken){
					rebound = true;
				}
			}

			

			//cin.ignore();
		}

		//Close binary
		result = orbitalElements(X, m1, m2);
		a[i] = get<0>(result);
		e[i] = get<1>(result);
		notBound = get<2>(result);

		if (rebound && !(notBound)){
			//cout << endl << "Rebound binary!" << endl;
			N_rebound += 1;
		}
		//cout << endl;
		//cout << "a = " << a[i] << endl;
		//cout << "e = " << e[i] << endl;
		//cout << "Ecc = " << Ecc << endl;
		//cin.ignore();
	}
	cout << endl;
	cout << "Number of rebound binaries = " << N_rebound << endl;
	return make_tuple(a, e);
}



tuple<long double, long double, long double> impulseEncounterPEF(long double m1, long double m2, long double M_p, long double p, long double e, long double b, long double v, long double f)
{
	//Star masses
	array<long double, 2> m = {m1, m2};
	//Open binary
	array<array<long double, 3>, 4> X = setupBinaryPEF(p, e, m1, m2, f);
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
	return orbitalElementsPEF(X, m1, m2);
}

tuple<vector<long double>, vector<long double>> MCEncountersPEF(long double v_rel, long double n_p, long double T, long double m1, long double m2, long double M_p, vector<long double> a, vector<long double> e)
{
	cout << "Running MCEncountersPEF" << endl;
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
	cout << "Number of binaries = " << N_bin << endl;
	//Declare variable types
	long double b_max, rate, v, b;
	tuple<long double, long double, long double> result;
	bool notBound;
	//Number of binaries broken
	int N_broken = 0;
	//double t_start;


	int N_enc;
	long double p, f, dt;
	bool hasBroken;
	bool rebound;
	int N_rebound = 0;
	int N_close = 0;

	//long double En_previous;
	//long double En;
	//bool notBound_previous;
	//ofstream myfile;
	//myfile.open("dE_break_binary_1000Msol.csv");

	//Iterate over binaries
	for (int i=0; i<N_bin; ++i){
		cout << '\r' << "Binary " << i+1 << " of " << N_bin << flush;

		hasBroken = false;
		rebound = false;
		notBound = false;

		b_max = calcBMax(M_p, v_rel, a[i], m1, m2);
		rate = encounterRate(n_p, v_rel, b_min, b_max, v_min, v_max);
		N_enc = randomPoisson(rate*T);

		p = a[i]*(1.0 - e[i]*e[i]);
		f = randomF(e[i]);

		//Implement encounters
		for (int j=0; j<N_enc; j++){
			//cout << '\r' << "Encounter " << j+1 << " of " << N_enc;

			//Draw velocity from distribution
			b = drawB(b_max);
			v = drawVMaxwellian(v_rel, v_max);

			result = impulseEncounterPEF(m1, m2, M_p, p, e[i], b, v, f);
			p = get<0>(result);
			e[i] = get<1>(result);
			f = get<2>(result);
			notBound = e[i]>= 1.0;

			if(p >= a_T){
				//cout << "Binary broken!" << endl;
				a[i] = -1.0;
				e[i] = -1.0;
				break;
			}



			//Evolve in time
			dt = randomExponential(rate);
			f = evolveF(f, dt, e[i], m1, m2, p);

			if (notBound){
				//cout << endl << "Binary broken!" << endl;
				//cout << endl;
				hasBroken = true;
			} else {
				if (hasBroken){
					rebound = true;
				}
			}


			
		}
		if (e[i] < 1.0){
			a[i] = p/(1.0 - e[i]*e[i]);
		} else {
			a[i] = -1.0;
			e[i] = -1.0;
		}

		if (rebound && (notBound == false)) {
			//cout << "Rebound binary!" << endl;
			N_rebound ++;
		}
		if ((a[i]>0.0) && (a[i]<100.0*parsec/length_scale) && notBound){
			N_close ++;
		}
		//cout << endl;
	}
	//myfile.close();
	cout << endl;
	cout << "Number of binaries rebound = " << N_rebound << endl;
	cout << "Number of binaries unbound within 100pc = " << N_close << endl;
	return make_tuple(a, e);
}