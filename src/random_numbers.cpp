// Generates random numbers
#include <random>
#include <limits>
#include "constants.h"
using namespace std;

const long double DBL_MAX = numeric_limits<long double>::max();

random_device rd;
mt19937 gen(seed);

long double randomUniformDoubleClosed(long double min, long double max)
{
	uniform_real_distribution<long double> closed_distribution(min, nextafter(max, DBL_MAX));
	return closed_distribution(gen);
}


long double randomUniformDoubleOpen(long double min, long double max)
{
	uniform_real_distribution<long double> open_distribution(min, max);
	return open_distribution(gen);
}

long double randomExponential(long double rate)
{
	exponential_distribution<long double> exp_distribution(rate);
	return exp_distribution(gen);
}

long double randomNormal(long double mean, long double stddev)
{
	normal_distribution<long double> norm_distribution(mean, stddev);
	return norm_distribution(gen);
}

int randomPoisson(long double mean)
{
	poisson_distribution<int> pois_distribution(mean);
	return pois_distribution(gen);
}