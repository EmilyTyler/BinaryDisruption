// Set of functions to generate binaries and to find their orbital elements
#include <array>
#include <math.h>
#include <iostream>
#include "constants.h"
#include "vector_maths.h"
#include "random_numbers.h"
using namespace std;


// Find the eccentric anomaly of a binary given its eccentricity e and mean anomaly M
double eccentricAnomaly(double e, double M){
	// Solves Kepler's equation E-esinE=M to find eccentric anomaly E given eccentricity e and mean anomaly M.
	// From page 36 of Solar System Dynamics, or Danby 1988.
	// Initial value for eccentric anomaly
	double E = M + copysign(1.0, sin(M))*0.85*e;
	// Initialise loop counter
	int count = 0;
	// Define other variables
	double f, f_p, f_pp, f_ppp, d_1, d_2, d_3;
	while (abs(E - e*sin(E) - M) > pow(10.0, -8.0)){
		f = E - e*sin(E) - M;
		f_p = 1.0 - e*cos(E);
		f_pp = e*sin(E);
		f_ppp = e*cos(E);

		d_1 = -f/f_p;
		d_2 = -f/(f_p + 0.5*d_1*f_pp);
		d_3 = -f/(f_p + 0.5*d_2*f_pp + pow(d_2, 2.0)*f_ppp/6.0);

		E += d_3;

		if (count > 100){
			cout << "eccentricAnomaly did not converge" << endl;
			break;
		}
	}
	return E;
}

//Return the semi-major axis and eccentricity of a binary and whether or not it is bound from the positions and velocities of the stars
tuple<double, double, bool> orbitalElements(array<double,(4,3)> X, double m1, double m2){
	// Separation vector
	array<double, 3> r = {X[0,0] - X[1,0], X[0,1] - X[1,1], X[0,2] - X[1,2]};
	// Relative velocity vector
	array<double, 3> v = {X[2,0] - X[3,0], X[2,1] - X[3,1], X[2,2] - X[3,2]};
	// Magnitudes of the above vectors
	double R = norm(r);
	double V = norm(v);
	// Total energy
	double E = m1*m2*(V*V/(2.0*(m1+m2)) - G/R);
	// Total angular momentum
	array<double, 3> L = cross(r, v);
	L[0] *= m1*m2/(m1+m2);
	L[0] *= m1*m2/(m1+m2);
	L[0] *= m1*m2/(m1+m2);
	double L_norm = norm(L);
	// Semi-major axis
	double a = G*m1*m2/(2.0*abs(E));
	// Eccentricity
	double e = sqrt(1.0 + 2.0*(m1+m2)*L_norm*L_norm*E/(G*G*pow(m1,3.0)*pow(m2,3.0)));
	// Not bound?
	bool notBound = E > 0.0;
	return make_tuple(a, e, notBound);
}

// Open a binary: find the position and velocity vectors given the semi-major axis and eccentricity
array<double,(4,3)> setupRandomBinary(double a, double e, double m1, double m2){
	// Randomise mean anomaly
	double M = randomUniformDoubleOpen(0.0, 2.0*pi);
	// Find eccentric anomaly
	double E = eccentricAnomaly(e, M);
	// Find true anomaly
	double f = 2.0*atan(sqrt((1.0+e)/(1.0-e))*tan(E/2.0));
	// Separation of stars
	double r = a*(1.0 - e*e)/(1.0 + e*cos(f));
	// Mean motion
	double n = sqrt(G*(m1+m2)/(pow(a,3.0)));
	// Position and velocity vectors
	array<double, (4,3)> = {
		{0.0, 0.0, 0.0},
		{r*cos(f), r*sin(f), 0.0},
		{0.0, 0.0, 0.0}, 
		{-n*a/(sqrt(1.0-e*e))*sin(f), n*a/(sqrt(1.0-e*e))*(e+cos(f)), 0.0}};
	// Centre of mass position vector
	array<double, 3> R;
	// Centre of mass velocity vector
	array<double,3> V;
	for (i=0, i<3, i++){
		R[i] = (m1*X[0,i] + m2*X[1,i])/(m1 + m2);
		V[i] = (m1*X[2,i] + m2*X[3,i])/(m1 + m2);
	}
	// Move into centre of mass rest frame
	for (i=0, i<3, i++){
		X[0,i] -= R[i];
		X[1,i] -= R[i];
		X[2,i] -= V[i];
		X[3,i] -= V[i];
	}
	return X;
}