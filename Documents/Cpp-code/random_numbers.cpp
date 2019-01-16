// Generates random numbers
#include <random>
#include <limits>
using namespace std;

const long double DBL_MAX = numeric_limits<long double>::max();

random_device rd;
mt19937 gen(rd());


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