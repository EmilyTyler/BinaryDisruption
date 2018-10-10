// Generates a random direction vector
#include "random_numbers.h"
#include "constants.h"
#include <math.h>
#include <array>
using namespace std;

array<double, 3> randomDirection()
{
	double u = randomUniformDoubleClosed(-1.0, 1.0);
	double theta = randomUniformDoubleOpen(0.0, 2.0*pi);
	array<double, 3> vec = {sqrt(1.0 - pow(u,2.0))*cos(theta), sqrt(1.0 - pow(u,2.0))*sin(theta), u};
	// Normalise
	double vec_norm = sqrt(pow(vec[0],2.0) + pow(vec[1],2.0) + pow(vec[2],2.0));
	for (int i=0; i<3; i++)
	{
		vec[i] /= vec_norm;
	}
	return vec;
}