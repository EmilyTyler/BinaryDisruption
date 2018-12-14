#include <iostream>
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
tuple<vector<double>, vector<double>> initialDistributions(int N_bin, double a_min, double a_max, double alpha){
	// Semimajor axis array
	vector<double> a;
	//Set size
	a.resize(N_bin);
	//Reduce capacity
	a.shrink_to_fit();
	//Declare variables
	double c;
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
	vector<double> e;
	//Set size
	e.resize(N_bin);
	//Reduce capacity
	e.shrink_to_fit();
	for (int i=0; i<N_bin; ++i){
		e[i] = pow(randomUniformDoubleClosed(), 1.0/3.0);
	}
	return make_tuple(a, e);
}

void evolvePopulation(string filename, int N_bin, double a_min, double a_max, double alpha, double v_rel, double n_p, double T, double m1, double m2, double M_p){
	//Initial semimajor axis and eccentricity distributions
	tuple<vector<double>, vector<double>> initial_dists = initialDistributions(N_bin, a_min, a_max, alpha);
	//Final semimajor axis and eccentricity distributions
	tuple<vector<double>, vector<double>> final_dists = MCEncounters(v_rel, n_p, T, m1, m2, M_p, get<0>(initial_dists), get<1>(initial_dists));
	//Extract results
	vector<double> a_fin = get<0>(final_dists);
	vector<double> e_fin = get<1>(final_dists);
	//Save results to file
	ofstream myfile;
	myfile.open(filename);
	for (int i; i<N_bin; ++i){
		myfile << a_fin[i] << ", " << e_fin[i] << endl; 
	}
    myfile.close();
}


int main() {

	double m1 = msol/mass_scale;
	double m2 = msol/mass_scale;
	double M_p = msol/mass_scale;
	double rho = 0.009 * msol/pow(parsec, 3.0) * (pow(length_scale, 3.0)/mass_scale);
	double n_p = rho/M_p;
	double v_rel = 2.2 * pow(10.0, 5.0) *(time_scale/length_scale);
	double T = 10.0 * giga * year /time_scale;

	double alpha = 1.0;
	double a_min = pow(10.0, 3.0) * au/length_scale;
	double a_max = pow(10.0, 6.0) * au/length_scale;

	int N_bin = 1;

	string filename = "binary_pop.csv";

	//Run simulation
	evolvePopulation(filename, N_bin, a_min, a_max, alpha, v_rel, n_p, T, m1, m2, M_p);


	//Test MCEncounters


}

