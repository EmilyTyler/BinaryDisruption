// Generates a random direction vector
#include "random_numbers.h"
#include "constants.h"
#include "vector_maths.h"
#include <cmath>
#include <array>
using namespace std;

//Tested
array<long double, 3> randomDirection()
{
	long double u = randomUniformDoubleClosed(-1.0, 1.0);
	long double theta = randomUniformDoubleOpen(0.0, 2.0*pi);
	array<long double, 3> vec = {sqrt(1.0 - pow(u,2.0))*cos(theta), sqrt(1.0 - pow(u,2.0))*sin(theta), u};
	// Normalise
	vec = normalise(vec);
	return vec;
}