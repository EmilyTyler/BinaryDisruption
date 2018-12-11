// Generates random numbers
#include <random>
#include <limits>
using namespace std;

const double DBL_MAX = numeric_limits<double>::max();

//Indirectly tested through random_direction
double randomUniformDoubleClosed(double min=0.0, double max=1.0)
{
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<double> distribution(min, nextafter(max, DBL_MAX));
	return distribution(gen);
}

//Indirectly tested through random_direction
double randomUniformDoubleOpen(double min=0.0, double max=1.0)
{
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<double> distribution(min, max);
	return distribution(gen);
}

//Tested
double randomExponential(double rate = 1.0)
{
	random_device rd;
	mt19937 gen(rd());
	exponential_distribution<double> distribution(rate);
	return distribution(gen);
}