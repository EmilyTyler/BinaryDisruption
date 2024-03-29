// Set of functions to generate binaries and to find their orbital elements
#include <array>
#include <cmath>
#include <iostream>
#include <tuple>
#include "constants.h"
#include "vector_maths.h"
#include "random_numbers.h"
using namespace std;

// Find the eccentric anomaly of a binary given its eccentricity e and mean anomaly M
long double eccentricAnomaly(long double e, long double M, bool &non_converged_binary){
	// Solves Kepler's equation E-esinE=M to find eccentric anomaly E given eccentricity e and mean anomaly M.
	// From page 36 of Solar System Dynamics, or Danby 1988.
	// Initial value for eccentric anomaly
	long double E = fmod(M + copysign(1.0L, sin(M))*0.85L*e, 2.0L*pi);
	// Initialise loop counter
	int count = 0;
	// Define other variables
	long double f, f_p, f_pp, f_ppp, d_1, d_2, d_3;
	while (abs(E - e*sin(E) - M) > pow(10.0L, -16.0L)){
		f = E - e*sin(E) - M;
		f_p = 1.0L - e*cos(E);
		f_pp = e*sin(E);
		f_ppp = e*cos(E);

		d_1 = -f/f_p;
		d_2 = -f/(f_p + 0.5L*d_1*f_pp);
		d_3 = -f/(f_p + 0.5L*d_2*f_pp + pow(d_2, 2.0L)*f_ppp/6.0L);

		E += d_3;

		count += 1;
		if (count > 100){
			cout << '\n' << "eccentricAnomaly did not converge" << ", E = " << E << ", M = " << M << ", e = " << e << endl;
			non_converged_binary = true;
			break;
		}
	}
	return fmod(E, 2.0L*pi);
}

long double eccentricAnomalyIonised(long double e, long double M, bool notBound, bool &non_converged_binary){
	// Solves Kepler's equation E-esinE=M to find eccentric anomaly E given eccentricity e and mean anomaly M.
	// From page 36 of Solar System Dynamics, or Danby 1988.
	if (notBound){
		// Initial value for eccentric anomaly
		long double E;
		if (M > 10000){
			E = log(M/e + sqrt(1.0L+M*M/(e*e)));
		} else {
			E = pi;
		}

		// Initialise loop counter
		int count = 0;
		// Define other variables
		long double f, f_p, f_pp, f_ppp, d_1, d_2, d_3;
		while (abs(E - e*sinh(E) + M) > pow(10.0L, -16.0L)){
				f = E - e*sinh(E) + M;
				f_p = 1.0L - e*cosh(E);
				f_pp = -e*sinh(E);
				f_ppp = -e*cosh(E);

				d_1 = -f/f_p;
				d_2 = -f/(f_p + 0.5L*d_1*f_pp);
				d_3 = -f/(f_p + 0.5L*d_2*f_pp + pow(d_2, 2.0L)*f_ppp/6.0L);

				E += d_3;

			count += 1;
			if (count > 100){
				non_converged_binary = true;
				break;
			}
		}
		return E;
	} else {
		return eccentricAnomaly(e, M, non_converged_binary);
	}
}

//Return the semi-major axis and eccentricity of a binary and whether or not it is bound from the positions and velocities of the stars
tuple<long double, long double, bool> orbitalElements(array<array<long double,3>, 4> X, long double m1, long double m2){
	// Separation vector
	array<long double, 3> r = {X[0][0] - X[1][0], X[0][1] - X[1][1], X[0][2] - X[1][2]};
	// Relative velocity vector
	array<long double, 3> v = {X[2][0] - X[3][0], X[2][1] - X[3][1], X[2][2] - X[3][2]};
	// Magnitudes of the above vectors
	long double R = norm(r);
	long double V = norm(v);
	// Total energy
	long double E = m1*m2*(V*V/(2.0L*(m1+m2)) - G/R);
	// Total angular momentum
	array<long double, 3> L = cross(r, v);
	L[0] *= m1*m2/(m1+m2);
	L[1] *= m1*m2/(m1+m2);
	L[2] *= m1*m2/(m1+m2);
	long double L_norm = norm(L);
	// Semi-major axis
	long double a = G*m1*m2/(2.0L*abs(E));
	// Eccentricity
	long double e = sqrt(1.0L + 2.0L*(m1+m2)*L_norm*L_norm*E/(G*G*pow(m1,3.0L)*pow(m2,3.0L)));
	// Not bound?
	bool notBound = E >= 0.0L;
	return make_tuple(a, e, notBound);
}
tuple<long double, long double, bool> orbitalElements(vector<array<long double,3>> X, long double m1, long double m2){
	// Separation vector
	array<long double, 3> r = {X[0][0] - X[1][0], X[0][1] - X[1][1], X[0][2] - X[1][2]};
	// Relative velocity vector
	array<long double, 3> v = {X[2][0] - X[3][0], X[2][1] - X[3][1], X[2][2] - X[3][2]};
	// Magnitudes of the above vectors
	long double R = norm(r);
	long double V = norm(v);
	// Total energy
	long double E = m1*m2*(V*V/(2.0L*(m1+m2)) - G/R);
	// Total angular momentum
	array<long double, 3> L = cross(r, v);
	L[0] *= m1*m2/(m1+m2);
	L[1] *= m1*m2/(m1+m2);
	L[2] *= m1*m2/(m1+m2);
	long double L_norm = norm(L);
	// Semi-major axis
	long double a = G*m1*m2/(2.0L*abs(E));
	// Eccentricity
	long double e = sqrt(1.0L + 2.0L*(m1+m2)*L_norm*L_norm*E/(G*G*pow(m1,3.0L)*pow(m2,3.0L)));
	// Not bound?
	bool notBound = (E >= 0.0L);
	return make_tuple(a, e, notBound);
}

tuple<long double, long double, long double, bool, long double> orbitalElementsIonised(array<array<long double,3>, 4> X, long double m1, long double m2){
	// Separation vector
	array<long double, 3> r = {X[0][0] - X[1][0], X[0][1] - X[1][1], X[0][2] - X[1][2]};
	// Relative velocity vector
	array<long double, 3> v = {X[2][0] - X[3][0], X[2][1] - X[3][1], X[2][2] - X[3][2]};
	// Magnitudes of the above vectors
	long double R = norm(r);
	long double V = norm(v);
	// Total energy
	long double E = m1*m2*(V*V/(2.0L*(m1+m2)) - G/R);
	// Total angular momentum
	array<long double, 3> L = cross(r, v);
	L[0] *= m1*m2/(m1+m2);
	L[1] *= m1*m2/(m1+m2);
	L[2] *= m1*m2/(m1+m2);
	long double L_norm = norm(L);
	// Semi-major axis
	long double a = G*m1*m2/(2.0L*abs(E));
	// Eccentricity
	long double e = sqrt(1.0L + 2.0L*(m1+m2)*L_norm*L_norm*E/(G*G*pow(m1,3.0L)*pow(m2,3.0L)));
	// Not bound?
	bool notBound = (E >= 0.0L);
	//Eccentric anomaly
	long double Ecc = 0.0L;
	if (E>0.0L){
		//Hyperbolic orbit
		Ecc = acosh((R/a + 1.0L)/e);
	} else if(E == 0.0L){
		//Parabolic orbit
		cout << "Energy = 0" << endl;
		Ecc = 0L;
	} else{
		//Elliptical orbit
		Ecc = acos((1.0L - R/a)/e);
		long double n = sqrt(G*(m1+m2)/(pow(a,3L)));
		long double f = asin(sqrt(1.0L-e*e)*V/(n*a*e));
		if (((0L<=Ecc<pi) && (pi<=f<2L*pi)) || ((pi<=Ecc<2.0L*pi) && (0L<=f<pi))){
			Ecc = 2.0L*pi - Ecc;
		} 
	}
	return make_tuple(a, e, Ecc, notBound, R);
}
tuple<long double, long double, long double, bool, long double> orbitalElementsIonised(vector<array<long double,3>> X, long double m1, long double m2){
	// Separation vector
	array<long double, 3> r = {X[0][0] - X[1][0], X[0][1] - X[1][1], X[0][2] - X[1][2]};
	// Relative velocity vector
	array<long double, 3> v = {X[2][0] - X[3][0], X[2][1] - X[3][1], X[2][2] - X[3][2]};
	// Magnitudes of the above vectors
	long double R = norm(r);
	long double V = norm(v);
	// Total energy
	long double E = m1*m2*(V*V/(2.0L*(m1+m2)) - G/R);
	// Total angular momentum
	array<long double, 3> L = cross(r, v);
	L[0] *= m1*m2/(m1+m2);
	L[1] *= m1*m2/(m1+m2);
	L[2] *= m1*m2/(m1+m2);
	long double L_norm = norm(L);
	// Semi-major axis
	long double a = G*m1*m2/(2.0L*abs(E));
	// Eccentricity
	long double e = sqrt(1.0L + 2.0L*(m1+m2)*L_norm*L_norm*E/(G*G*pow(m1,3.0L)*pow(m2,3.0L)));
	// Not bound?
	bool notBound = (E >= 0.0L);
	//Eccentric anomaly
	long double Ecc = 0.0L;
	if (E>0.0L){
		//Hyperbolic orbit
		Ecc = acosh((R/a + 1.0L)/e);
	} else if(E == 0.0L){
		//Parabolic orbit
		cout << "Energy = 0" << endl;
		Ecc = 0L;
	} else{
		//Elliptical orbit
		Ecc = acos((1.0L - R/a)/e);
		long double n = sqrt(G*(m1+m2)/(pow(a,3L)));
		long double f = asin(sqrt(1.0L-e*e)*V/(n*a*e));
		if (((0L<=Ecc<pi) && (pi<=f<2L*pi)) || ((pi<=Ecc<2.0L*pi) && (0L<=f<pi))){
			Ecc = 2.0L*pi - Ecc;
		} 
	}
	return make_tuple(a, e, Ecc, notBound, R);
}

// Open a binary: find the position and velocity vectors given the semi-major axis and eccentricity
array<array<long double, 3>, 4> setupRandomBinary(long double a, long double e, long double m1, long double m2){
	// Randomise mean anomaly
	long double M = randomUniformDoubleOpen(0.0L, 2.0L*pi);
	// Find eccentric anomaly
	bool arg3 = false;
	long double E = eccentricAnomaly(e, M, arg3);
	// Find true anomaly
	long double f = 2.0L*atan(sqrt((1.0L+e)/(1.0L-e))*tan(E/2.0L));
	// Separation of stars
	long double r = a*(1.0L - e*e)/(1.0L + e*cos(f));
	// Mean motion
	long double n = sqrt(G*(m1+m2)/(pow(a,3L)));
	// Position and velocity vectors
	array<array<long double, 3>, 4> X = { {
		{0.0L, 0.0L, 0.0L},
		{r*cos(f), r*sin(f), 0.0L},
		{0.0L, 0.0L, 0.0L}, 
		{-n*a/(sqrt(1.0L-e*e))*sin(f), n*a/(sqrt(1.0L-e*e))*(e+cos(f)), 0.0L}} };
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
	return X;
}

vector<array<long double, 3>> setupRandomBinaryVector(long double a, long double e, long double m1, long double m2){
	// Randomise mean anomaly
	long double M = randomUniformDoubleOpen(0.0L, 2.0L*pi);
	// Find eccentric anomaly
	bool arg3 = false;
	long double E = eccentricAnomaly(e, M, arg3);
	// Mean motion
	long double n = sqrt(G*(m1+m2)/(pow(a,3L)));
	// Position and velocity vectors
	vector<array<long double, 3>> X = { {
		{0.0L, 0.0L, 0.0L},
		{a*(cos(E)-e), a*sqrt(1.0L-e*e)*sin(E), 0.0L},
		{0.0L, 0.0L, 0.0L}, 
		{-n*a/(1.0L-e*cos(E))*sin(E), n*a/(1.0L-e*cos(E))*sqrt(1.0L-e*e)*cos(E), 0.0L}} };
	X.shrink_to_fit();
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
	return X;
}

array<array<long double, 3>, 4> setupRandomBinaryIonised(long double a, long double e, long double m1, long double m2, long double E, bool notBound, bool &non_converged_binary){
	long double M, f, n, r;
	array<array<long double, 3>, 4> X;
	//Mean motion
	n = sqrt(G*(m1+m2)/(pow(a,3.0L)));
	if (notBound){
		// Position and velocity vectors
		X = { {
			{0.0L, 0.0L, 0.0L},
			{a*(cosh(E) - e), a*sqrt(e*e-1.0L)*sinh(E), 0.0L},
			{0.0L, 0.0L, 0.0L}, 
			{n*a*sinh(E)/(e*cosh(E) - 1.0L), n*a*sqrt(e*e-1.0L)*cosh(E)/(e*cosh(E) - 1.0L), 0.0L}} };
	} else {
		// Randomise mean anomaly
		M = randomUniformDoubleOpen(0.0, 2.0*pi);
		//Find eccentric anomaly
		E = eccentricAnomaly(e, M, non_converged_binary);
		// Find true anomaly
		f = 2.0*atan(sqrt((1.0+e)/(1.0-e))*tan(E/2.0));
		// Separation of stars
		r = a*(1.0 - e*e)/(1.0 + e*cos(f));
		// Position and velocity vectors
		X = { {
			{0.0, 0.0, 0.0},
			{r*cos(f), r*sin(f), 0.0},
			{0.0, 0.0, 0.0}, 
			{-n*a/(sqrt(1.0-e*e))*sin(f), n*a/(sqrt(1.0-e*e))*(e+cos(f)), 0.0}} };
	}
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
	return X;
}
