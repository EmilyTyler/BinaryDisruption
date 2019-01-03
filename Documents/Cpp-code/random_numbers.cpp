// Generates random numbers
#include <random>
#include <limits>
using namespace std;

const long double DBL_MAX = numeric_limits<long double>::max();

random_device rd;
mt19937 gen(rd());
uniform_real_distribution<long double> closed_distribution(0.0, nextafter(1.0, DBL_MAX));
uniform_real_distribution<long double> open_distribution(0.0, 1.0);


long double randomUniformDoubleClosed()
{
	return closed_distribution(gen);
}


long double randomUniformDoubleOpen()
{
	return open_distribution(gen);
}

long double randomExponential(long double rate)
{
	exponential_distribution<long double> exp_distribution(rate);
	return exp_distribution(gen);
}