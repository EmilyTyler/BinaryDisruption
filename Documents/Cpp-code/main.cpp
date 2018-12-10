#include <iostream>
#include <array>
#include <fstream>
#include "random_direction.h"
using namespace std;

int main() {
	array<double, 3> arr;

	ofstream myfile;
	myfile.open("test_data.csv");
	for (int i; i<1000; i++){
		arr = randomDirection();
		myfile << arr[0] << "," << arr[1] << "," << arr[2] << endl; 
	}

}

