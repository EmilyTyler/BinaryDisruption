#include <iostream>
#include <iomanip>
#include <array>
#include <fstream>
#include <cmath>
#include <vector>
#include <tuple>
#include <string>
#include "constants.h"
#include "random_numbers.h"
#include "encounters.h"
#include "vector_maths.h"
#include "random_direction.h"
#include "binary.h"
#include "nbodyintegration.h"
using namespace std;

//Draw a from distribution dN/dloga\propto a^{1-alpha}
//Draw e from distribution uniform in e^2 between 0 and 1
tuple<vector<long double>, vector<long double>> initialDistributions(int N_bin, long double a_min, long double a_max, long double alpha){
	// Semimajor axis array
	vector<long double> a;
	//Set size
	a.resize(N_bin);
	//Reduce capacity
	a.shrink_to_fit();
	//Declare variables
	long double c;
	if (alpha == 1.0L){
		c = log(a_min)/log(a_max/a_min);
		for (int i=0; i<N_bin; ++i){
			a[i] = pow(a_max/a_min, (randomUniformDoubleClosed(0.0L, 1.0L) + c));
		}
	} else {
		for (int i=0; i<N_bin; ++i){
			a[i] = pow((randomUniformDoubleClosed(0.0L, 1.0L)*(pow(a_max, 1.0L-alpha) - pow(a_min, 1.0L-alpha)) + pow(a_min, 1.0L-alpha)), 1.0L/(1.0L - alpha));
		}
	}
	//Eccentricity array
	vector<long double> e;
	//Set size
	e.resize(N_bin);
	//Reduce capacity
	e.shrink_to_fit();
	for (int i=0; i<N_bin; ++i){
		e[i] = pow(randomUniformDoubleClosed(0.0L, 1.0L), 1.0L/3.0L);
	}
	return make_tuple(a, e);
}

void evolvePopulation(string filename, int N_bin, long double a_min, long double a_max, long double alpha, long double v_rel, long double n_p, long double T, long double m1, long double m2, long double M_p){
	//Initial semimajor axis and eccentricity distributions
	cout << "Generating initial binaries" << endl;
	tuple<vector<long double>, vector<long double>> initial_dists = initialDistributions(N_bin, a_min, a_max, alpha);
	vector<long double> a_ini = get<0>(initial_dists);
	vector<long double> e_ini = get<1>(initial_dists);

	cout << setprecision(16) << "M_p, M_sol = " << M_p*mass_scale/msol << endl;

	//Final semimajor axis and eccentricity distributions
	ofstream myfile;
	myfile.open(filename);
	cout << "Evolving binaries" << endl;
	tuple<vector<long double>, vector<long double>, vector<long double>, vector<long double>> final_dists = MCEncountersIonised(v_rel, n_p, T, m1, m2, M_p, a_ini, e_ini);
	//Extract results
	vector<long double> a_fin = get<0>(final_dists);
	vector<long double> e_fin = get<1>(final_dists);
	vector<long double> r_ini = get<2>(final_dists);
	vector<long double> r_fin = get<3>(final_dists);
	//Save results to file
	cout << "Saving" << endl;
	for (int i=0; i<N_bin; ++i){
		myfile << setprecision(16) << a_ini[i]*length_scale << "," << r_ini[i]*length_scale << "," << e_ini[i] << "," << a_fin[i]*length_scale << "," << r_fin[i]*length_scale << "," << e_fin[i] << endl; 
	}
    myfile.close();
    cout << "Finished" << endl;
}

int main() {

	long double m1 = 0.5*msol/mass_scale;
	long double m2 = 0.5*msol/mass_scale;
	long double M_p = 10.0L*msol/mass_scale;
	long double rho = 0.009L * msol/pow(parsec, 3.0L) * (pow(length_scale, 3.0L)/mass_scale);
	long double n_p = rho/M_p;
	long double v_rel = 220.0L*1000.0L *(time_scale/length_scale);
	long double T = 10.0L * giga * year /time_scale;

	long double alpha = 1.0L;
	long double a_min = pow(10.0L, 1.0L) * au/length_scale;
	long double a_max = pow(10.0L, 5.5L) * au/length_scale;

	int N_bin = pow(10,5);

	string filename = "final_r_and_a_distributions_rho0_009_Mp10_vrel_220_Nbin10e5_format_ai_ri_ei_af_rf_ef.csv";

	//Run simulation
	evolvePopulation(filename, N_bin, a_min, a_max, alpha, v_rel, n_p, T, m1, m2, M_p);
	
	return 1;
}

