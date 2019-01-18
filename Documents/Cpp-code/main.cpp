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
	if (alpha == 2.0){
		c = log(a_min)/log(a_max/a_min);
		for (int i=0; i<N_bin; ++i){
			a[i] = pow(a_max/a_min, (randomUniformDoubleClosed(0.0, 1.0) + c));
		}
	} else {
		for (int i=0; i<N_bin; ++i){
			a[i] = pow((randomUniformDoubleClosed(0.0, 1.0)*(pow(a_max, 2.0-alpha) - pow(a_min, 2.0-alpha)) + pow(a_min, 2.0-alpha)), 1.0/(2.0 - alpha));
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
	tuple<vector<long double>, vector<long double>> initial_dists = initialDistributions(N_bin, a_min, a_max, alpha);
	//Final semimajor axis and eccentricity distributions
	tuple<vector<long double>, vector<long double>> final_dists = MCEncounters(v_rel, n_p, T, m1, m2, M_p, get<0>(initial_dists), get<1>(initial_dists));
	//Extract results
	vector<long double> a_fin = get<0>(final_dists);
	vector<long double> e_fin = get<1>(final_dists);
	//Save results to file
	ofstream myfile;
	myfile.open(filename);
	for (int i=0; i<N_bin; ++i){
		myfile << setprecision(16) << a_fin[i]*length_scale << ", " << e_fin[i] << endl; 
	}
    myfile.close();
}

void WSWEncounterTest(string filename, long double m1, long double m2, long double M_p, long double a, long double e, long double v){
	//Number of encounters for each b
	const unsigned int N_enc = pow(10, 10);
	//b's to run encounters
	const int N_b = 1;
	array<long double, N_b> b = {4.0};
	//array<long double, N_b> b = {3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0};
	for(int i=0; i<N_b; ++i){
		b[i] = pow(10.0,b[i])*au/length_scale;
	}
	//Declare variables
	tuple<long double, long double, long double, long double, long double> result;
	long double E_ini, E_fin, b_star, dE_v_dv, dE_dv_dv, b_input;
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
	//Maximum energy change
	long double dE_max = m1*m2/(m1+m2)*(sqrt(G*(m1+m2)*(1.0+e)/(a*(1.0-e)))*(2.0*G*M_p*a*(1.0+e)/(b[0]*b[0]*v)) + 0.5*(2.0*G*M_p*a*(1.0+e)/(b[0]*b[0]*v))*(2.0*G*M_p*a*(1.0+e)/(b[0]*b[0]*v)));
	unsigned int N_enc_so_far = 0;
	int counter = 0;
	long double dE_mean = 0.0;
	long double dE2_mean = 0.0;
	long double dE_mean_old = 0.0;
	long double std_dev; 
	long double b_star_min = 100.0*b[0]*length_scale;
	long double b_star_max = 0.0; 
	long double b_max = calcBMax(M_p, v, a, m1, m2);
	while(N_enc_so_far < N_enc){
		for(int i=0; i<N_b; ++i){
			
			b_input = drawB(b_max);
			result = testImpulseEncounter(m1, m2, M_p, a, e, b_input, v);
			//Convert to SI units
			E_ini = get<0>(result) * mass_scale*(length_scale*length_scale/(time_scale*time_scale));
			E_fin = get<1>(result) * mass_scale*(length_scale*length_scale/(time_scale*time_scale));
			
			b_star = get<2>(result) * length_scale;

			//dE_v_dv = get<3>(result) * mass_scale*(length_scale*length_scale/(time_scale*time_scale));
			//dE_dv_dv = get<4>(result) * mass_scale*(length_scale*length_scale/(time_scale*time_scale));

			//cout << "b_star = " << b_star/au << endl;

			//cout << "Minimum impact parameter, au = " << b_star_min/au << endl;
			//cout << "Maximum impact parameter, au = " << b_star_max/au << endl;
			//Write to file
			//myfile << setprecision(16) << E_ini << ", " << E_fin << ", " << b_star << endl;
			
			//if ((0.9*b[0] < b_star/length_scale) && (b_star/length_scale < 1.1*b[0])){
				N_enc_so_far += 1;
				//myfile << setprecision(16) << E_fin-E_ini << " , " << dE_v_dv << " , " << dE_dv_dv << endl;

				
				dE_mean = dE_mean*(N_enc_so_far-1)/N_enc_so_far + (E_fin-E_ini)/N_enc_so_far;
				dE2_mean = dE2_mean*(N_enc_so_far-1)/N_enc_so_far + (E_fin-E_ini)*(E_fin-E_ini)/N_enc_so_far;
				
				if (N_enc_so_far > pow(10.0, counter*0.1)-1){
					std_dev = sqrt(dE2_mean - dE_mean*dE_mean);
					cout << setprecision(16) << dE_mean << " , " << std_dev << " , " << N_enc_so_far << endl;
					myfile << setprecision(16) << dE_mean << " , " << std_dev << " , " << N_enc_so_far << endl;
					counter += 1;
					
				}
				
				
				b_star_min = min(b_star_min, b_star);
				b_star_max = max(b_star_max, b_star);
				/*
				dE_mean_old = dE_mean;
				
				if (copysign(1, dE_mean) != copysign(1, dE_mean_old)){
					cout << endl;
					cout << "Number of encounters so far = " << N_enc_so_far << endl;
					cout << "Old mean energy change = " << dE_mean_old << endl;
					cout << "Energy change = " << E_fin-E_ini << endl;
					cout << "Maximum energy change = " << dE_max* mass_scale*(length_scale*length_scale/(time_scale*time_scale)) << endl;
					cout << "New mean energy change = " << dE_mean << endl;
					cout << "New standard deviation = " << sqrt(dE2_mean - dE_mean*dE_mean) << endl;
					cout << "Analytical average energy change = " << dE_avg_analytic * mass_scale*(length_scale*length_scale/(time_scale*time_scale)) << endl;
					cout << "Analytical standard deviation = " << std_dev_analytic * mass_scale*(length_scale*length_scale/(time_scale*time_scale)) << endl;
					cout << "Number required for convergence = " << pow(std_dev_analytic/(0.1*dE_avg_analytic) ,2.0) << endl;
					cout << endl;
				}
				*/

			//}
			

		}
	}
	myfile.close();
	cout << "Analytical average energy change = " << dE_avg_analytic * mass_scale*(length_scale*length_scale/(time_scale*time_scale)) << endl;
	cout << "Analytical standard deviation = " << std_dev_analytic * mass_scale*(length_scale*length_scale/(time_scale*time_scale)) << endl;
	cout << "Number required for convergence = " << pow(std_dev_analytic/(0.1*dE_avg_analytic) ,2.0) << endl;
	cout << "Minimum impact parameter, au = " << b_star_min/au << endl;
	cout << "Maximum impact parameter, au = " << b_star_max/au << endl;
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
	long double T = 1.0 * giga * year /time_scale;

	long double alpha = 1.0;
	long double a_min = pow(10.0, 3.0) * au/length_scale;
	long double a_max = pow(10.0, 6.0) * au/length_scale;

	int N_bin = 1;

	string filename = "binary_pop.csv";

	//Run simulation
	evolvePopulation(filename, N_bin, a_min, a_max, alpha, v_rel, n_p, T, m1, m2, M_p);
	*/

	//To do
	//Test MCEncounters
	//Re-test MC_velocity

	//testBAndVVectors();
		
	
	//Test impulse approx against WSW
	string filename = "WSW_encounters_N_enc_log.csv";

	long double m1 = 2.0*msol/mass_scale;
	long double m2 = 2.0*msol/mass_scale;
	long double M_p = 3.0*msol/mass_scale;
	long double a = pow(10.0, 5.0) * au/length_scale;
	long double e = 0.7;
	long double v = 2.2 * pow(10.0, 5.0) *(time_scale/length_scale);

	WSWEncounterTest(filename, m1, m2, M_p, a, e, v);
	
}

