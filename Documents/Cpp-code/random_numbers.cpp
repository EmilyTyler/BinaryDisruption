// Generates random numbers
#include <random>
#include <limits>
using namespace std;

const double DBL_MAX = numeric_limits<double>::max();

random_device rd;
mt19937 gen(rd());
uniform_real_distribution<double> closed_distribution(0.0, nextafter(1.0, DBL_MAX));
uniform_real_distribution<double> open_distribution(0.0, 1.0);


double randomUniformDoubleClosed()
{
	return closed_distribution(gen);
}


double randomUniformDoubleOpen()
{
	return open_distribution(gen);
}

double randomExponential(double rate)
{
	exponential_distribution<double> exp_distribution(rate);
	return exp_distribution(gen);
}