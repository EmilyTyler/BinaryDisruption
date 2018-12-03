// Generates random numbers
#include <random>
#include <limits>
using namespace std;

const double DBL_MAX = numeric_limits<double>::max();

double randomUniformDoubleClosed(double min=0.0, double max=1.0)
{
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<double> distribution(min, nextafter(max, DBL_MAX));
	return distribution(gen);
}

double randomUniformDoubleOpen(double min=0.0, double max=1.0)
{
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<double> distribution(min, max);
	return distribution(gen);
}

double randomExponential(double rate)
{
	random_device rd;
	mt19937 gen(rd());
	exponential_distribution<double> distribution(rate);
	return distribution(gen);
}