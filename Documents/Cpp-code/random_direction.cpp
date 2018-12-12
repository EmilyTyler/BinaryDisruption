// Generates a random direction vector
#include "random_numbers.h"
#include "constants.h"
#include "vector_maths.h"
#include <cmath>
#include <array>
using namespace std;

//Tested
array<double, 3> randomDirection()
{
	double u = randomUniformDoubleClosed()*2.0 -1.0;
	double theta = randomUniformDoubleOpen()*2.0*pi;
	array<double, 3> vec = {sqrt(1.0 - pow(u,2.0))*cos(theta), sqrt(1.0 - pow(u,2.0))*sin(theta), u};
	// Normalise
	vec = normalise(vec);
	return vec;
}