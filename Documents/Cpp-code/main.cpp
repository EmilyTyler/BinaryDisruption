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
	if (alpha == 2.0){
		c = log(a_min)/log(a_max/a_min);
		for (int i=0; i<N_bin; ++i){
			a[i] = pow(a_max/a_min, (randomUniformDoubleClosed() + c));
		}
	} else {
		for (int i=0; i<N_bin; ++i){
			a[i] = pow((randomUniformDoubleClosed()*(pow(a_max, 2.0-alpha) - pow(a_min, 2.0-alpha)) + pow(a_min, 2.0-alpha)), 1.0/(2.0 - alpha));
		}
	}
	//Eccentricity array
	vector<long double> e;
	//Set size
	e.resize(N_bin);
	//Reduce capacity
	e.shrink_to_fit();
	for (int i=0; i<N_bin; ++i){
		e[i] = pow(randomUniformDoubleClosed(), 1.0/3.0);
	}
	return make_tuple(a, e);
}

void evolvePopulation(string filename, int N_bin, long double a_min, long double a_max, long double alpha, long double v_rel, long double n_p, long double T, long double m1, long double m2, long double M_p){
	//Initial semimajor axis and eccentricity distributions
	tuple<vector<long double>, vector<long double>> initial_dists = initialDistributions(N_bin, a_min, a_max, alpha);
	//Final semimajor axis and eccentricity distributions
	tuple<vector<long double>, vector<long double>> final_dists = MCEncounters(v_rel, n_p, T, m1, m2, M_p, get<0>(initial_dists), get<1>(initial_dists));
	//Extract results
	vector<long double> a_fin = get<0>(final_dists);
	vector<long double> e_fin = get<1>(final_dists);
	//Save results to file
	ofstream myfile;
	myfile.open(filename);
	for (int i; i<N_bin; ++i){
		myfile << a_fin[i] << ", " << e_fin[i] << endl; 
	}
    myfile.close();
}

void WSWEncounterTest(string filename, long double m1, long double m2, long double M_p, long double a, long double e, long double v){
	//Number of encounters for each b
	int N_enc = pow(10.0, 4.0);
	//b's to run encounters
	int N_b = 11;
	array<long double,11> b = {3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0};
	for(int i=0; i<N_b; ++i){
		b[i] = pow(10.0,b[i])*au/length_scale;
	}
	//Initial energies
	array<long double, 110000> E_ini;
	//Final energies
	array<long double, 110000> E_fin;
	//Actual impact parameters
	array<long double, 110000> b_star;
	//Run encounters
	tuple<long double, long double, long double> result;
	cout << "Simulating encounters" << endl;
	for(int i=0; i<N_b; ++i){
		for(int j=0; j<N_enc; ++j){
			result = testImpulseEncounter(m1, m2, M_p, a, e, b[i], v);
			E_ini[i*N_enc + j] = get<0>(result) * mass_scale*(length_scale*length_scale/(time_scale*time_scale));
			E_fin[i*N_enc + j] = get<1>(result) * mass_scale*(length_scale*length_scale/(time_scale*time_scale));
			b_star[i*N_enc + j] = get<2>(result) * length_scale;
			if (get<2>(result) == 0.0){
				cout << "here" << endl;
			}

		}
	}
	//Write to file
	cout << "Writing to file" << endl;
	ofstream myfile;
	myfile.open(filename);
	for (int i=0; i<N_b*N_enc; ++i){
		myfile << setprecision(16) << E_ini[i] << ", " << E_fin[i] << ", " << b_star[i] << endl; 
	}
    myfile.close();
    cout << "Finished" << endl;
}


int main() {
	/*
	long double m1 = msol/mass_scale;
	long double m2 = msol/mass_scale;
	long double M_p = msol/mass_scale;
	long double rho = 0.009 * msol/pow(parsec, 3.0) * (pow(length_scale, 3.0)/mass_scale);
	long double n_p = rho/M_p;
	long double v_rel = 2.2 * pow(10.0, 5.0) *(time_scale/length_scale);
	long double T = 10.0 * giga * year /time_scale;

	long double alpha = 1.0;
	long double a_min = pow(10.0, 3.0) * au/length_scale;
	long double a_max = pow(10.0, 6.0) * au/length_scale;

	int N_bin = 1;

	string filename = "binary_pop.csv";

	//Run simulation
	evolvePopulation(filename, N_bin, a_min, a_max, alpha, v_rel, n_p, T, m1, m2, M_p);
	*/

	//Test MCEncounters


	//Test impulse approx against WSW
	string filename = "WSW_encounters.csv";

	long double m1 = 2.0*msol/mass_scale;
	long double m2 = 2.0*msol/mass_scale;
	long double M_p = 3.0*msol/mass_scale;
	long double a = pow(10.0, 5.0) * au/length_scale;
	long double e = 0.7;
	long double v = 2.2 * pow(10.0, 5.0) *(time_scale/length_scale);

	WSWEncounterTest(filename, m1, m2, M_p, a, e, v);


}

