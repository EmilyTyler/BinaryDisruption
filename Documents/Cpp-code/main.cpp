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

void testBAndVVectors(){
	long double b = pow(10.0, 5.0)*au/length_scale;
	long double v = pow(10.0, 5.0)*2.2/length_scale*time_scale;
	tuple<array<long double,3>, array<long double,3>> result;
	array<long double,3> b_vec, v_vec;
	long double b_norm, v_norm;
	ofstream myfile;
	myfile.open("test_data.csv");
	for (int i=0; i<pow(10, 6); ++i){
		result = impactAndVelocityVectors(b, v);
		b_vec = get<0>(result);
		v_vec = get<1>(result);
		b_norm = norm(b_vec);
		v_norm = norm(v_vec);
		if (abs((b_norm - b)/b) > pow(10.0, -16.0)){
			cout << "Wrong magnitude! b = " << b << " , b_norm = " << b_norm << endl; 
		}
		if (abs((v_norm - v)/v) > pow(10.0, -16.0)){
			cout << "Wrong magnitude! v = " << v << " , v_norm = " << v_norm << endl; 
		}
		if (abs(dot(b_vec, v_vec)) > pow(10.0,-14.0)){
			cout << "Dot product not equal to zero! b dot v = " << dot(b_vec, v_vec) << endl;
		}
		b_vec = normalise(b_vec);
		v_vec = normalise(v_vec);
		myfile << setprecision(16) << v_vec[0] << " , " << v_vec[1] << " , " << v_vec[2] << endl;
	} 
	myfile.close();
}

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
	if (alpha == 1.0){
		c = log(a_min)/log(a_max/a_min);
		for (int i=0; i<N_bin; ++i){
			a[i] = pow(a_max/a_min, (randomUniformDoubleClosed(0.0, 1.0) + c));
		}
	} else {
		for (int i=0; i<N_bin; ++i){
			a[i] = pow((randomUniformDoubleClosed(0.0, 1.0)*(pow(a_max, 1.0-alpha) - pow(a_min, 1.0-alpha)) + pow(a_min, 1.0-alpha)), 1.0/(1.0 - alpha));
		}
	}
	//Eccentricity array
	vector<long double> e;
	//Set size
	e.resize(N_bin);
	//Reduce capacity
	e.shrink_to_fit();
	for (int i=0; i<N_bin; ++i){
		e[i] = pow(randomUniformDoubleClosed(0.0, 1.0), 1.0/3.0);
	}
	return make_tuple(a, e);
}

void evolvePopulation(string filename, int N_bin, long double a_min, long double a_max, long double alpha, long double v_rel, long double n_p, long double T, long double m1, long double m2, long double M_p){
	//Initial semimajor axis and eccentricity distributions
	cout << "Generating initial binaries" << endl;
	tuple<vector<long double>, vector<long double>> initial_dists = initialDistributions(N_bin, a_min, a_max, alpha);
	vector<long double> a_ini = get<0>(initial_dists);
	vector<long double> e_ini = get<1>(initial_dists);

	for(int i=0; i<N_bin; i++){
		a_ini[i] = parsec/length_scale;
	}
	cout << setprecision(16) << "M_p, M_sol = " << M_p*mass_scale/msol << endl;
	cout << "a_ini, au = " << a_ini[0]*length_scale/au << endl;

	//Final semimajor axis and eccentricity distributions
	//ofstream myfile;
	//myfile.open(filename);
	cout << "Evolving binaries" << endl;
	//tuple<vector<long double>, vector<long double>> final_dists = MCEncountersXV(v_rel, n_p, T, m1, m2, M_p, a_ini, e_ini);
	tuple<vector<long double>, vector<long double>> final_dists = MCEncountersIonised(v_rel, n_p, T, m1, m2, M_p, a_ini, e_ini);
	//tuple<vector<long double>, vector<long double>> final_dists = MCEncountersNClosest(100, v_rel, n_p, T, m1, m2, M_p, a_ini, e_ini);
	//Extract results
	vector<long double> a_fin = get<0>(final_dists);
	vector<long double> e_fin = get<1>(final_dists);
	//Save results to file
	//cout << "Saving" << endl;
	//for (int i=0; i<N_bin; ++i){
	//	myfile << setprecision(16) << a_ini[i]*length_scale << " , " << e_ini[i] << " , " << a_fin[i]*length_scale << ", " << e_fin[i] << endl; 
	//}
    //myfile.close();
    cout << "Finished" << endl;
}

void WSWEncounterTest(string filename, long double m1, long double m2, long double M_p, long double a, long double e, long double v){
	//Number of encounters for each b
	const unsigned int N_enc = pow(10, 8);
	//b's to run encounters
	const int N_b = 1;
	array<long double, N_b> b = {6.0};
	//array<long double, N_b> b = {3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0};
	for(int i=0; i<N_b; ++i){
		b[i] = pow(10.0,b[i])*au/length_scale;
	}
	//Declare variables
	tuple<long double, long double, long double, long double, long double, array<long double, 3>, array<long double, 3>, long double> result;
	long double E_ini, E_fin, b_star, dE_v_dv, dE_dv_dv, b_input, delta_v_norm, v_initial_norm, theta, phi;
	array<long double,3> v_initial, delta_v;
	cout << "Simulating encounters" << endl;	
	ofstream myfile;
	myfile.open(filename);
	//Theoretical average energy change
	long double dE_avg_analytic;
	//Theoretical standard deviation
	long double std_dev_analytic;
	if (b[0] < a){
		dE_avg_analytic = 2.0*(G*M_p/(b[0]*v))*(G*M_p/(b[0]*v));
		std_dev_analytic = sqrt(4.0/3.0*G*(m1+m2)/a*(G*M_p/(b[0]*v))*(G*M_p/(b[0]*v)));
	} else{
		dE_avg_analytic = 4.0/3.0 * (G*M_p/(b[0]*v))*(G*M_p/(b[0]*v)) * (a/b[0])*(a/b[0]) * (1.0 + 3.0*e*e/2.0);
		std_dev_analytic = sqrt(4.0/5.0*G*(m1+m2)/a*(G*M_p/(b[0]*v))*(G*M_p/(b[0]*v))*(a/b[0])*(a/b[0])*(1.0 - e*e/3.0) + 16.0/45.0*pow(G*M_p/(b[0]*v), 4.0)*pow(a/b[0], 4.0)*(1.0 +15.0*e*e));
	}
	dE_avg_analytic *= m1*m2/(m1+m2);
	//Maximum energy change
	long double dE_max = m1*m2/(m1+m2)*(sqrt(G*(m1+m2)*(1.0+e)/(a*(1.0-e)))*(2.0*G*M_p*a*(1.0+e)/(b[0]*b[0]*v)) + 0.5*(2.0*G*M_p*a*(1.0+e)/(b[0]*b[0]*v))*(2.0*G*M_p*a*(1.0+e)/(b[0]*b[0]*v)));
	long double dE_min = m1*m2/(m1+m2)*(-sqrt(G*(m1+m2)*(1.0+e)/(a*(1.0-e)))*(2.0*G*M_p*a*(1.0+e)/(b[0]*b[0]*v)) + 0.5*(2.0*G*M_p*a*(1.0+e)/(b[0]*b[0]*v))*(2.0*G*M_p*a*(1.0+e)/(b[0]*b[0]*v)));
	unsigned int N_enc_so_far = 0;
	int counter = 0;
	long double dE_mean = 0.0;
	long double dE2_mean = 0.0;
	long double dE_v_dv_mean = 0.0;
	long double dE_dv_dv_mean = 0.0;
	long double dE_mean_old = 0.0;
	long double std_dev; 
	long double b_star_min = 100.0*b[0]*length_scale;
	long double b_star_max = 0.0; 
	long double b_max = calcBMax(M_p, v, a, m1, m2);

	tuple<array<long double,3>, array<long double,3>> result2;
	while(N_enc_so_far < N_enc){
		for(int i=0; i<N_b; ++i){
			
			b_input = b[0];
			result = testImpulseEncounter(m1, m2, M_p, a, e, b_input, v);
			//Convert to SI units
			E_ini = get<0>(result) * mass_scale*(length_scale*length_scale/(time_scale*time_scale));
			E_fin = get<1>(result) * mass_scale*(length_scale*length_scale/(time_scale*time_scale));
			
			b_star = get<2>(result) * length_scale;

			dE_v_dv = get<3>(result) * mass_scale*(length_scale*length_scale/(time_scale*time_scale));
			dE_dv_dv = get<4>(result) * mass_scale*(length_scale*length_scale/(time_scale*time_scale));

			v_initial = get<5>(result);
			for (int j=0; j<3; j++){
				v_initial[j] *= length_scale/time_scale;
			}
			v_initial_norm = norm(v_initial);

			delta_v = get<6>(result);
			for (int j=0; j<3; j++){
				delta_v[j] *= length_scale/time_scale;
			}
			delta_v_norm = norm(delta_v);

			theta = get<7>(result);

			//cout << "b_star = " << b_star/au << endl;

			//cout << "Minimum impact parameter, au = " << b_star_min/au << endl;
			//cout << "Maximum impact parameter, au = " << b_star_max/au << endl;
			//Write to file
			//myfile << setprecision(16) << E_ini << ", " << E_fin << ", " << b_star << endl;
			
			if ((0.9*b[0] < b_star/length_scale) && (b_star/length_scale < 1.1*b[0])){
				N_enc_so_far += 1;
				//myfile << setprecision(16) << E_fin-E_ini << " , " << dE_v_dv << " , " << dE_dv_dv << endl;

				
				dE_mean = dE_mean*(N_enc_so_far-1)/N_enc_so_far + (E_fin-E_ini)/N_enc_so_far;
				dE2_mean = dE2_mean*(N_enc_so_far-1)/N_enc_so_far + (E_fin-E_ini)*(E_fin-E_ini)/N_enc_so_far;
				dE_v_dv_mean = dE_v_dv_mean*(N_enc_so_far-1)/N_enc_so_far + dE_v_dv/N_enc_so_far;
				dE_dv_dv_mean = dE_dv_dv_mean*(N_enc_so_far-1)/N_enc_so_far + dE_dv_dv/N_enc_so_far;
				
				if (N_enc_so_far > pow(10.0, counter*0.1)-1){
				//if (N_enc_so_far % static_cast<int>(pow(10, 6)) == 0){
					std_dev = sqrt(dE2_mean - dE_mean*dE_mean);
					//cout << setprecision(16) << dE_mean << " , " << std_dev << " , " << N_enc_so_far << endl;
					myfile << setprecision(16) << dE_mean << " , " << std_dev << " , " << N_enc_so_far << endl;
					//cout << setprecision(16) << E_fin-E_ini<< " , " << dE_v_dv << " , " << dE_dv_dv << endl;
					//myfile << setprecision(16) << E_fin-E_ini<< " , " << dE_v_dv << " , " << dE_dv_dv << endl;
					//v_initial = get<1>(impactAndVelocityVectors(b[0], v));
					//v_initial_norm = norm(v_initial);
					//myfile << setprecision(16) << dE_v_dv << " , " << asin(abs(v_initial[2])/v_initial_norm) << endl;
					//cout << setprecision(16) << v_initial_norm << " , " << delta_v_norm << " , " << cos(theta) << " , " << dE_v_dv << endl;
					//myfile << setprecision(16) << v_initial_norm << " , " << delta_v_norm << " , " << cos(theta) << " , " << dE_v_dv << endl;
					counter += 1;
					
				}
				
				
				b_star_min = min(b_star_min, b_star);
				b_star_max = max(b_star_max, b_star);
				
							
				if (copysign(1, dE_mean) != copysign(1, dE_mean_old)){
					cout << endl;
					cout << "Number of encounters so far = " << N_enc_so_far << endl;
					cout << "Old mean energy change = " << dE_mean_old << endl;
					cout << "Energy change = " << E_fin-E_ini << endl;
					cout << "Maximum energy change = " << dE_max* mass_scale*(length_scale*length_scale/(time_scale*time_scale)) << endl;
					cout << "Minimum energy change = " << dE_min* mass_scale*(length_scale*length_scale/(time_scale*time_scale)) << endl;
					cout << "New mean energy change = " << dE_mean << endl;
					cout << "New standard deviation = " << sqrt(dE2_mean - dE_mean*dE_mean) << endl;
					cout << "Analytical average energy change = " << dE_avg_analytic * mass_scale*(length_scale*length_scale/(time_scale*time_scale)) << endl;
					cout << "Analytical standard deviation = " << std_dev_analytic * mass_scale*(length_scale*length_scale/(time_scale*time_scale)) << endl;
					cout << "Number required for convergence = " << pow(std_dev_analytic/(0.1*dE_avg_analytic) ,2.0) << endl;
					cout << endl;
				}
				
				dE_mean_old = dE_mean;
				/*
				if (abs(E_fin-E_ini) < pow(10.0, 25.0)){
					cout << E_fin-E_ini << " , " << dE_v_dv + dE_dv_dv -(E_fin-E_ini)<< endl;
				}
				*/
			}
			

		}
	}
	myfile.close();
	cout << "Analytical average energy change = " << dE_avg_analytic * mass_scale*(length_scale*length_scale/(time_scale*time_scale)) << endl;
	cout << "Average of v dv term = " << dE_v_dv_mean << endl;
	cout << "Average of dv dv term = " << dE_dv_dv_mean << endl;
	cout << "Average energy change = " << dE_mean << endl;
	cout << endl;
	cout << "Analytical standard deviation = " << std_dev_analytic * mass_scale*(length_scale*length_scale/(time_scale*time_scale)) << endl;
	cout << "Number required for convergence = " << pow(std_dev_analytic/(0.1*dE_avg_analytic) ,2.0) << endl;
	cout << "Minimum impact parameter, au = " << b_star_min/au << endl;
	cout << "Maximum impact parameter, au = " << b_star_max/au << endl;
    cout << "Finished" << endl;
}

void BHT_survival_probability(){
	cout << "Initialising" << endl;
	//Input parameters
	//Mass of binary
	long double m1 = msol / mass_scale;
	long double m2 = msol / mass_scale;
	long double M_b = m1 + m2;
	//Perturber mass
	long double M_p = 3.0 * msol / mass_scale;
	//Relative velocity dispersion
	long double v_rel = sqrt(1.0/3.0) * pow(10.0, 5.0) /length_scale*time_scale;
	//Density of dark matter halo
	long double rho = 0.1 *msol/(pow(parsec, 3.0)) /mass_scale*pow(length_scale, 3.0);
	//Number density of perturbers
	long double n_p = rho/M_p;
	//Total simulation time
	long double T = 10.0*giga*year / time_scale;
	//Initial semi-major axis
	long double a_0 = 0.1*parsec / length_scale;
	//Eccentricity
	long double e_0 = 0.7;
	//Number of binaries per simulation
	int N_bin = 1000;
	//Number of simulations
	int N_sim = 1;
	//Starting index in file names
	int i_start = 0;

	//Time steps
	//Minimum impact parameter
	long double b_min = 0.0;
	//Maximum impact parameter
	long double b_max = calcBMax(M_p, v_rel, a_0, m1, m2);
	//Minimum velocity
	long double v_min = 0.0;
	//Maximum velocity
	long double v_max = 100.0 * v_rel;
	//Encounter rate
	long double rate = encounterRate(n_p, v_rel, b_min, b_max, v_min, v_max);
	//Time step
	long double dt = 0.5/rate;
	//Number of timesteps
	int N_t = static_cast<int>(floor(T/dt)) + 1;
	//Adjust timestep
	dt = T/(N_t-1);
	//Time array
	vector<long double> t;
	t.resize(N_t);
	t.shrink_to_fit();
	//Number of binaries broken array
	vector<long double> N_broken;
	N_broken.resize(N_t);
	N_broken.shrink_to_fit();

	//Semi-major axis and eccentricity arrays
	vector<long double> a;
	a.resize(N_bin);
	a.shrink_to_fit();
	vector<long double> e;
	e.resize(N_bin);
	e.shrink_to_fit();

	tuple<vector<long double>, vector<long double>, int, int, int, int, int> result;
	bool previous_number_zero;

	int N_encounters = 0;
	int N_encounters_close = 0;
	int N_encounters_far = 0;
	int N_encounters_mid = 0;

	//Run simulations
	for (int i=0; i<N_sim; i++){
		cout << "Simulation " << i+1 << " of " << N_sim << endl;
		//Initialise semi-major axis and eccentricity arrays
		a.resize(N_bin);
		e.resize(N_bin);
		for (int j=0; j<N_bin; j++){
			a[j] = a_0;
			e[j] = e_0;
		}
		//Initialise time and number broken arrays
		t.resize(N_t);
		N_broken.resize(N_t);
		for (int j=0; j<N_t; j++){
			t[j] = j*dt;
			N_broken[j] = 0;
		}

		for(int j=1; j<N_t; j++){
			result = MCEncounters(v_rel, n_p, t[j]-t[j-1], m1, m2, M_p, a, e);
			a = get<0>(result);
			e = get<1>(result);
			N_broken[j] = static_cast<long double>(get<2>(result));
			a = where_positive(a);
			e = where_positive(e);
			N_encounters += get<3>(result);
			N_encounters_close += get<4>(result);
			N_encounters_far += get<5>(result);
			N_encounters_mid += get<6>(result);
		}


		cout << "Filtering" << endl;
		//Filter out zeros
		previous_number_zero = false;
		for (int j=N_t-1; j>0; j--){
			if (N_broken[j] < 1){
				if (previous_number_zero) {
					N_broken[j] = -1;
					t[j] = -1.0;
				}
				previous_number_zero = true;
			} else {
				previous_number_zero = false;
			}
		}
		N_broken = where_positive(N_broken);
		t = where_positive(t);
		int N_N_broken = static_cast<int>(N_broken.size());
		//Make number broken cumulative
		for (int j=2; j<N_N_broken; j++){
			N_broken[j] += N_broken[j-1];
		}
		//Normalise to find fraction broken
		for (int j=0; j<N_N_broken; j++){
			N_broken[j] /= N_bin;
		}
		//Convert fraction broken to survival fraction
		for (int j=0; j<N_N_broken; j++){
			N_broken[j] = 1 - N_broken[j];
		}
		
		//Save data
		cout << "Saving" << endl;
		int file_index = i_start+i;
		string filename = "BHTfig2_mysim_new_5bmax_" + to_string(N_bin) + "bin_" + to_string(file_index) + ".csv";
		ofstream myfile;
		myfile.open(filename);
		for (int j=0; j<N_N_broken; j++){
			myfile << setprecision(16) << t[j]*time_scale << " , " << N_broken[j] << endl;
		}
		myfile.close();
		
	}
	cout << "Total number of encounters = " << N_encounters << endl;
	cout << "Average number of encounters = " << N_encounters/N_bin/N_sim << endl;
	cout << "Number of encounters at b<a = " << N_encounters_close << endl;
	cout << "Number of encounters at b>a = " << N_encounters_far << endl;
	cout << "Number of encounters between close and far regimes = " << N_encounters_mid << endl;
	cout << endl;
	cout << "Finished" << endl;	
}

void WSWEncounterTest_MeanvB(string filename, long double m1, long double m2, long double M_p, long double a, long double e, long double v){
	//Number of encounters for each b
	const unsigned int N_enc = pow(10, 6);
	//b's to run encounters 
	const int N_b = 41;
	//array<long double, N_b> b = {3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0};
	array<long double, N_b> b = {3.0, 3.125, 3.25, 3.375, 3.5, 3.625, 3.75, 3.875, 4.0, 4.125, 4.25, 4.375, 4.5, 4.625, 4.75, 4.875, 5.0, 5.125, 5.25, 5.375, 5.5, 5.625, 5.75, 5.875, 6.0, 6.125, 6.25, 6.375, 6.5, 6.625, 6.75, 6.875, 7.0, 7.125, 7.25, 7.375, 7.5, 7.625, 7.75, 7.875, 8.0};
	for(int i=0; i<N_b; ++i){
		b[i] = pow(10.0,b[i])*au/length_scale;
	}
	//Declare variables
	tuple<long double, long double, long double, long double, long double, array<long double,3>, array<long double,3>, long double> result;
	long double E_ini, E_fin, b_star, dE_v_dv, dE_dv_dv;
	cout << "Simulating encounters" << endl;	
	ofstream myfile;
	myfile.open(filename);
	unsigned int N_enc_so_far = 0;
	long double dE_mean = 0.0;
	long double dE2_mean = 0.0;
	long double dE_v_dv_mean = 0.0;
	long double dE_dv_dv_mean = 0.0;
	long double std_dev; 
	//Theoretical average energy change
	long double dE_avg_analytic;
	//Theoretical standard deviation
	long double std_dev_analytic;
	for(int i=0; i<N_b; ++i){
		N_enc_so_far = 0;
		dE_mean = 0.0;
		dE2_mean = 0.0;
		dE_v_dv_mean = 0.0;
		dE_dv_dv_mean = 0.0;
		while(N_enc_so_far < N_enc){
			result = testImpulseEncounter(m1, m2, M_p, a, e, b[i], v);
			//Convert to SI units
			E_ini = get<0>(result) * mass_scale*(length_scale*length_scale/(time_scale*time_scale));
			E_fin = get<1>(result) * mass_scale*(length_scale*length_scale/(time_scale*time_scale));
			b_star = get<2>(result) * length_scale;
			dE_v_dv = get<3>(result) * mass_scale*(length_scale*length_scale/(time_scale*time_scale));
			dE_dv_dv = get<4>(result) * mass_scale*(length_scale*length_scale/(time_scale*time_scale));
			
			if ((0.9*b[i] < b_star/length_scale) && (b_star/length_scale < 1.1*b[i])){
				N_enc_so_far += 1;
				dE_mean += E_fin-E_ini;
				dE2_mean += (E_fin-E_ini)*(E_fin-E_ini);
				dE_v_dv_mean += dE_v_dv;
				dE_dv_dv_mean += dE_dv_dv;
			}
		}
		//Normalise
		dE_mean /= N_enc;
		dE2_mean /= N_enc;
		dE_v_dv_mean /= N_enc;
		dE_dv_dv_mean /= N_enc;
		std_dev = sqrt(dE2_mean - dE_mean*dE_mean);
		//cout << setprecision(16) << dE_mean << " , " << std_dev << " , " << b[i]*length_scale << endl;
		//myfile << setprecision(16) << dE_mean << " , " << std_dev << " , " << b[i]*length_scale << endl;

		if (b[i] < a){
			dE_avg_analytic = 2.0*(G*M_p/(b[i]*v))*(G*M_p/(b[i]*v));
			std_dev_analytic = sqrt(4.0/3.0*G*(m1+m2)/a*(G*M_p/(b[i]*v))*(G*M_p/(b[i]*v)));
		} else{
			dE_avg_analytic = 4.0/3.0 * (G*M_p/(b[i]*v))*(G*M_p/(b[i]*v)) * (a/b[i])*(a/b[i]) * (1.0 + 3.0*e*e/2.0);
			std_dev_analytic = sqrt(4.0/5.0*G*(m1+m2)/a*(G*M_p/(b[i]*v))*(G*M_p/(b[i]*v))*(a/b[i])*(a/b[i])*(1.0 - e*e/3.0) + 16.0/45.0*pow(G*M_p/(b[i]*v), 4.0)*pow(a/b[i], 4.0)*(1.0 +15.0*e*e));
		}
		dE_avg_analytic *= m1*m2/(m1+m2);
		cout << setprecision(16) << dE_dv_dv_mean / (dE_avg_analytic* mass_scale*(length_scale*length_scale/(time_scale*time_scale))) << " , " << b[i]*length_scale << endl;
		myfile << setprecision(16) << dE_dv_dv_mean / (dE_avg_analytic* mass_scale*(length_scale*length_scale/(time_scale*time_scale))) << " , " << b[i]*length_scale << endl;
	}
	myfile.close();
    cout << "Finished" << endl;
}


//Solve recurrence relation for max N_enc as a function of a_0
void recurrenceSolve(){
	array<long double,4> a_0 = {3.0, 4.0, 5.0, 6.0};
	long double delta = pow(10.0, -3.0);
	long double m1 = msol/mass_scale;
	long double m2 = msol/mass_scale;

	long double M_p = 1000.0*msol/mass_scale;
	long double v_rel = 2.0*pow(10.0, 5.0)*time_scale/length_scale;

	long double M_b = m1+m2;
	M_b *= 1.0;
	const int N_a = static_cast<int>(a_0.size());
	array<long double,N_a> E_0;
	for(int i=0; i<N_a; i++){
		a_0[i] = pow(10.0, a_0[i]) * au/length_scale;
		E_0[i] = -G*M_b/(2.0*a_0[i]);
	}
	int N_enc;
	long double E;
	long double dE;
	for(int i=0; i<N_a; i++){
		N_enc = 0;
		E = E_0[i];
		dE = 5.0/2.0*pow(G, 3.0/2.0)*delta*M_p/v_rel*sqrt(M_b)*pow(a_0[i], -3.0/2.0);
		while (E < 0.0){
			//E = E - 189.0/512.0*delta*delta*pow(E_0[i],3.0)/(pow(E,2.0));
			E = E + dE;
			N_enc += 1;
			//cout << "E = " << E << ", " << "N_enc = " << N_enc << endl;
		}
		cout << "a_0, au = " << a_0[i]*length_scale/au << ", " << "N_enc = " << N_enc << endl;
	}
}

void testEvolve(){
	vector<long double> M = {msol/mass_scale, msol/mass_scale};
	long double a = pow(10.0, 1.0)*au/length_scale;
	long double e = 1.5;
	long double T = 1.0*year/time_scale;
	bool arg3 = false;
	array<array<long double, 3>, 4> X = setupRandomBinaryIonised(a, e, M[0], M[1], 0.0, 0.0, true, arg3, false);
	vector<array<long double, 3>> X_0;
	X_0.resize(4);
	for (int i=0; i<4; ++i){
		for (int j=0; j<3; ++j){
			X_0[i][j] = X[i][j];
		}
	}
	X_0.shrink_to_fit();
	cout << setprecision(16) << "Initial energy, internal units = " << -G*M[0]*M[1]/(2.0*a) << endl;
	X_0 = evolve(2, M, X_0, T);
	return;
}

void testImpulseApprox(){
	tuple<long double, long double, bool> result;
	long double a_imp, e_imp, a_thr, e_thr, w;
	bool notBound;
	array<long double, 3> x_pbh, v_vec, b_vec;
	array<array<long double, 3>, 4> X;
	vector<array<long double, 3>> X_thr, X_imp;
	tuple<array<long double,3>, array<long double,3>> bvvectors;
	bool ini_arrays;

	long double b_90, b_star_norm, v_perp, v_para;
	array<long double,3> b_star;

	long double m1 = msol/mass_scale;
	long double m2 = msol/mass_scale;
	long double M_p = 1000.0*msol/mass_scale;
	long double rho = 0.009 * msol/pow(parsec, 3.0) * (pow(length_scale, 3.0)/mass_scale);
	long double n_p = rho/M_p;
	long double v_rel = 200.0*1000.0 *(time_scale/length_scale);
	vector<long double> M = {m1, m2, M_p};
	vector<long double> m = {m1, m2};
	long double T = 10.0 * giga * year /time_scale;

	long double v = v_rel;
	long double a = pow(10.0, 4.0) * au;
	long double e = 0.7;
	int N_enc = 100;



	long double b_min = sqrt(M_p/(pi*rho*v_rel*T));
	long double b_max = calcBMax(M_p, v_rel, a, m1, m2);
	int N_b = 10;
	long double db = (log(b_max) - log(b_min))/(N_b - 1);
	vector<long double> bs;
	bs.resize(N_b);
	for (int i=0; i<N_b; i++){
		bs[i] = b_min*exp(i*db);
	}
	vector<long double> a_frac_avg;
	a_frac_avg.resize(N_b);

	ofstream myfile;
	myfile.open("impulse_test_a10e4au_Mp1000Msol.csv");

	for (int l=0; l<N_b; l++){
		cout << '\r' << "Impact parameter " << l+1 << " of " << N_b << flush;
		a_frac_avg[l] = 0.0;
		for(int i=0; i<N_enc; i++){
			//cout << '\r' << "Encounter " << i+1 << " of " << N_enc << flush;

			w = sqrt(pow(10.0, 6.0)*M_p*a*a/(min(m1, m2)) - bs[l]*bs[l])/v;
			X = setupRandomBinary(a, e, m1, m2);
			bvvectors = impactAndVelocityVectors(bs[l], v);
			b_vec = get<0>(bvvectors);
			v_vec = get<1>(bvvectors);

			X_imp = {X[0], X[1], X[2], X[3]};
			X_imp = evolve(2, m, X_imp, w, ini_arrays=true);
			for (int j=0; j<2; ++j){
				//90 degree deflection radius
				b_90 = G*(M_p + m[j])/(v*v);
				//Calculate impact parameter for this star
				b_star = calcBStar(X_imp[j], v_vec, v, b_vec);
				//Calculate norm of b_star
				b_star_norm = norm(b_star);
				//Calculate speed change in b_star direction
				v_perp = 2.0*M_p*v/(m[j]+M_p) * (b_star_norm/b_90)/(1.0 + b_star_norm*b_star_norm/(b_90*b_90));
				//Calculate speed change in -v_vec direction
				v_para = 2.0*M_p*v/(m[j]+M_p) * 1.0/(1.0 + b_star_norm*b_star_norm/(b_90*b_90));
				//Change star velocity
				for (int k=0; k<3; ++k){
					X_imp[j+2][k] += v_perp * b_star[k]/b_star_norm - v_para * v_vec[k]/v;
				}
			}
			X_imp = evolve(2, m, X_imp, w, ini_arrays=true);
			result = orbitalElements(X_imp, m1, m2);
			a_imp = get<0>(result);
			e_imp = get<1>(result);
			notBound = get<2>(result);


			x_pbh = {b_vec[0]-w*v_vec[0], b_vec[1]-w*v_vec[1], b_vec[2]-w*v_vec[2]};
			X_thr = {X[0], X[1], x_pbh, X[2], X[3], v_vec};
			X_thr = evolve(3, M, X_thr, 2.0*w, ini_arrays = true);
			X_thr.erase(X_thr.begin()+2);
			X_thr.erase(X_thr.begin()+4);
			result = orbitalElements(X_thr, m1, m2);
			a_thr = get<0>(result);
			e_thr = get<1>(result);

			a_frac_avg[l] += (a_imp-a_thr)/a_thr/N_b;

			//cout << setprecision(16) << "a_imp = " << a_imp << endl;
			//cout << setprecision(16) << "a_thr = " << a_thr << endl;
			//cin.ignore();
		}
		cout << endl;
	}
	cout << "Saving" << endl;
	for (int i=0; i< N_b; i++){
		myfile << setprecision(16) << a_frac_avg[i] << ", " << bs[i]*length_scale << endl;
	}
	myfile.close();
	cout << "Finished" << endl;
}



int main() {

	//recurrenceSolve();
	

	long double m1 = msol/mass_scale;
	long double m2 = msol/mass_scale;
	long double M_p = 1.0*msol/mass_scale;
	long double rho = 0.009 * msol/pow(parsec, 3.0) * (pow(length_scale, 3.0)/mass_scale);
	long double n_p = rho/M_p;
	long double v_rel = 200.0*1000.0 *(time_scale/length_scale);
	long double T = 10.0 * giga * year /time_scale;

	long double alpha = 1.0;
	long double a_min = pow(10.0, 1.0) * au/length_scale;
	long double a_max = pow(10.0, 5.5) * au/length_scale;

	int N_bin = pow(10, 4);

	string filename = "";

	//Test evolve
	//testEvolve();

	//Run simulation
	evolvePopulation(filename, N_bin, a_min, a_max, alpha, v_rel, n_p, T, m1, m2, M_p);
	
	//testImpulseApprox();

	//To do
	//Test MCEncounters
	//Re-test MC_velocity

	//testBAndVVectors();
		
	
	//Test impulse approx against WSW
	
	//string filename = "WSW_encounters_N_enc_b10e6au_log.csv";
	/*
	long double m1 = msol/mass_scale;
	long double m2 = msol/mass_scale;
	long double M_p = 3.0*msol/mass_scale;
	long double a = pow(10.0, 5.0) * au/length_scale;
	long double e = 0.0;
	long double v = 2.2 * pow(10.0, 5.0) *(time_scale/length_scale);
	*/
	//WSWEncounterTest(filename, m1, m2, M_p, a, e, v);
	
	//WSWEncounterTest_MeanvB(filename, m1, m2, M_p, a, e, v);

	//BHT_survival_probability();

	/*
	long double M;
	long double E;
	long double f;
	long double r;
	long double n;
	array<long double,3> v_vec;
	long double v_norm;
	ofstream myfile;
	myfile.open("test_data.csv");
	for (int i=0; i<10000000; i++){
		M = randomUniformDoubleOpen(0.0, 2.0*pi);
		E = eccentricAnomaly(e, M);
		f = 2.0*atan(sqrt((1.0+e)/(1.0-e))*tan(E/2.0));
		r = a*(1.0 - e*e)/(1.0 + e*cos(f));
		n = sqrt(G*(m1+m2)/(pow(a,3)));
		v_vec = {-n*a/(sqrt(1.0-e*e))*sin(f), n*a/(sqrt(1.0-e*e))*(e+cos(f)), 0.0};
		for (int j=0; j<3; j++){
			v_vec[j] *= length_scale/time_scale;
		}
		v_norm = norm(v_vec);
		myfile << setprecision(16) << v_norm << endl;
	}
	myfile.close();
	*/

	/*

	int N_bin=1000;
	//Mass of binary
	long double m1 = msol / mass_scale;
	long double m2 = msol / mass_scale;
	//Perturber mass
	long double M_p = 3.0 * msol / mass_scale;
	//Relative velocity dispersion
	long double v_rel = sqrt(2.0/3.0) * pow(10.0, 5.0) /length_scale*time_scale;
	//Density of dark matter halo
	long double rho = 0.1 *msol/(pow(parsec, 3.0)) /mass_scale*pow(length_scale, 3.0);
	//Number density of perturbers
	long double n_p = rho/M_p;
	//Total simulation time
	long double T = 10.0*giga*year / time_scale;
	//Initial semi-major axis
	long double a_0 = 0.1*parsec / length_scale;
	//Eccentricity
	long double e_0 = 0.7;
	tuple<vector<long double>, vector<long double>, int, int, int, int, int> result;

	//Semi-major axis and eccentricity arrays
	vector<long double> a;
	a.resize(N_bin);
	a.shrink_to_fit();
	vector<long double> e;
	e.resize(N_bin);
	e.shrink_to_fit();
	for (int j=0; j<N_bin; j++){
		a[j] = a_0;
		e[j] = e_0;
	}
	result = MCEncounters(v_rel, n_p, T, m1, m2, M_p, a, e);
	*/
	return 1;
}

