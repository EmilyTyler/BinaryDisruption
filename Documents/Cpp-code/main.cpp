#include <iostream>
#include <array>
#include <fstream>
#include <cmath>
#include "constants.h"
#include "random_direction.h"
#include "random_numbers.h"
#include "vector_maths.h"
#include "MC_velocity.h"
#include "binary.h"
#include "encounters.h"
using namespace std;

int main() {

	double a = 0.1*parsec;
	double e = 0.7;
	double m1 = msol;
	double m2 = msol;
	double M_p = msol;
	double rho = 0.009 * msol/pow(parsec, 3.0);
	double n_p = rho/M_p;
	cout << "n_p = " << n_p << endl;
	double b_min = 0.0;
	double v_rel = 2.2 * pow(10.0, 5.0);

	double b_max = calcBMax(M_p, v_rel, a, m1, m2);
	cout << "b_max = " << b_max << endl;

	double b = 0.05*b_max;
	cout << "b = " << b << endl;

	//Test impulseEncounter
	tuple<double, double, bool> result = impulseEncounter(m1, m2, M_p, a, e, b, v_rel);

	//Test drawB

	//Test MCEncounters



	/*
	double number;
	ofstream myfile;
	myfile.open("test_data.csv");
	for (int i; i<1000000; i++){
		number = drawVMaxwellian(50.0, 0.01*50.0, 100.0*50.0);
		myfile << number << endl; 
	}
    myfile.close();
	*/

}

