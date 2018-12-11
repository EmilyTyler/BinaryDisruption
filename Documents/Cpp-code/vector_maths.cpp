//Functions for basic vector operations
#include <array>
#include <cmath>
using namespace std;

//Tested
array<double, 3> cross(array<double,3> x, array<double,3> y)
{
	array<double, 3> z = {x[1]*y[2] - x[2]*y[1], x[2]*y[0] - x[0]*y[2], x[0]*y[1] - x[1]*y[0]};
	return z;
}

//Tested
double norm(array<double,3> x)
{
	return sqrt(pow(x[0], 2.0) + pow(x[1], 2.0) + pow(x[2], 2.0));
}

//Tested
array<double, 3> normalise(array<double, 3> x)
{
	double x_norm = norm(x);
	for (int i=0; i<3; i++){
		x[i] /= x_norm;
	}
	return x;
}

//Tested
double dot(array<double,3> x, array<double,3> y)
{
	return (x[0]*y[0] + x[1]*y[1] + x[2]*y[2]);
}