//Functions for basic vector operations
#include <array>
using namespace std;

array<double, 3> cross(array<double,3> x, array<double,3> y)
{
	array<double, 3> z = {x[1]*y[2] - x[2]*y[1], x[2]*y[0] - x[0]*y[2], x[0]*y[1] - x[1]*y[0]};
	return z;
}

inline double norm(array<double,3> x)
{
	return sqrt(pow(x[0], 2.0) + pow(x[1], 2.0) + pow(x[2], 2.0))
}