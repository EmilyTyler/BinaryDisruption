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
	long double u = randomUniformDoubleClosed(-1.0L, 1.0L);
	long double theta = randomUniformDoubleOpen(0.0L, 2.0L*pi);
	array<long double, 3> vec = {sqrt(1.0L - pow(u,2.0L))*cos(theta), sqrt(1.0L - pow(u,2.0L))*sin(theta), u};
	// Normalise
	vec = normalise(vec);
	return vec;
}