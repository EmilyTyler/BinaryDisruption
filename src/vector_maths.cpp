//Functions for basic vector operations
#include <array>
#include <vector>
#include <cmath>
using namespace std;

array<double, 3> cross(array<double,3> x, array<double,3> y)
{
	array<double, 3> z = {x[1]*y[2] - x[2]*y[1], x[2]*y[0] - x[0]*y[2], x[0]*y[1] - x[1]*y[0]};
	return z;
}
array<long double, 3> cross(array<long double,3> x, array<long double,3> y)
{
	array<long double, 3> z = {x[1]*y[2] - x[2]*y[1], x[2]*y[0] - x[0]*y[2], x[0]*y[1] - x[1]*y[0]};
	return z;
}

double norm(array<double,3> x)
{
	return sqrt(pow(x[0], 2.0) + pow(x[1], 2.0) + pow(x[2], 2.0));
}
long double norm(array<long double,3> x)
{
	return sqrt(pow(x[0], 2.0L) + pow(x[1], 2.0L) + pow(x[2], 2.0L));
}

array<double, 3> normalise(array<double, 3> x)
{
	double x_norm = norm(x);
	for (int i=0; i<3; ++i){
		x[i] /= x_norm;
	}
	return x;
}
array<long double, 3> normalise(array<long double, 3> x)
{
	long double x_norm = norm(x);
	for (int i=0; i<3; ++i){
		x[i] /= x_norm;
	}
	return x;
}

double dot(array<double,3> x, array<double,3> y)
{
	return (x[0]*y[0] + x[1]*y[1] + x[2]*y[2]);
}
long double dot(array<long double,3> x, array<long double,3> y)
{
	return (x[0]*y[0] + x[1]*y[1] + x[2]*y[2]);
}

vector<long double> where_positive(vector<long double> x)
{
	int N_x;
	bool reached_the_end = false;
	while (reached_the_end == false){
		N_x = static_cast<int>(x.size());
		for (int j=0; j<N_x; j++){
			if (j==N_x - 1){
				reached_the_end = true;
			}
			if (x[j] < 0){
				x.erase(x.begin() + j);
				break;
			}
		}
	}
	return x;
}

vector<long double> rotate(vector<long double> vec, long double x_angle, long double y_angle, long double z_angle){
	long double x_new = cos(y_angle)*cos(z_angle)*vec[0] + cos(y_angle)*sin(z_angle)*vec[1] - sin(y_angle)*vec[2];
	long double y_new = (sin(x_angle)*sin(y_angle)*cos(z_angle) - cos(x_angle)*sin(z_angle))*vec[0] + (sin(x_angle)*sin(y_angle)*sin(z_angle) + cos(x_angle)*cos(z_angle))*vec[1] + sin(x_angle)*cos(y_angle)*vec[2];
	long double z_new = (cos(x_angle)*sin(y_angle)*cos(z_angle) + sin(x_angle)*sin(z_angle))*vec[0] +(cos(x_angle)*sin(y_angle)*sin(z_angle) - sin(x_angle)*cos(z_angle))*vec[1] + cos(x_angle)*cos(y_angle)*vec[2];
	vec = {x_new, y_new, z_new};
	return vec; 
}
array<long double, 3> rotate(array<long double, 3> vec, long double x_angle, long double y_angle, long double z_angle){
	long double x_new = cos(y_angle)*cos(z_angle)*vec[0] + cos(y_angle)*sin(z_angle)*vec[1] - sin(y_angle)*vec[2];
	long double y_new = (sin(x_angle)*sin(y_angle)*cos(z_angle) - cos(x_angle)*sin(z_angle))*vec[0] + (sin(x_angle)*sin(y_angle)*sin(z_angle) + cos(x_angle)*cos(z_angle))*vec[1] + sin(x_angle)*cos(y_angle)*vec[2];
	long double z_new = (cos(x_angle)*sin(y_angle)*cos(z_angle) + sin(x_angle)*sin(z_angle))*vec[0] +(cos(x_angle)*sin(y_angle)*sin(z_angle) - sin(x_angle)*cos(z_angle))*vec[1] + cos(x_angle)*cos(y_angle)*vec[2];
	vec = {x_new, y_new, z_new};
	return vec; 
}