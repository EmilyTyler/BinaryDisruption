//Functions to draw a random number from a v times maxwellian distribution
#include <cmath>
#include "random_numbers.h"
using namespace std;

//Indirectly tested through drawVMaxwellian
double VMaxwellianPdf(double x, double v_rel)
{
	return pow(x,3.0)/(2.0*pow(v_rel,4.0))*exp(-x*x/(2.0*v_rel*v_rel));
}

//Indirectly tested through drawVMaxwellian
double VMaxwellianComparison(double x, double v_rel, double v_min, double v_max, double x_max)
{
	if (v_min<x<v_max){
		return VMaxwellianPdf(x_max, v_rel);
	} else{
		return 0.0;
	}
}

//Indirectly tested through drawVMaxwellian
double VMaxwellianXFromArea(double A, double v_rel, double v_min, double x_max){
	return v_min + A/VMaxwellianPdf(x_max, v_rel);
}

//Tested
double drawVMaxwellian(double v_rel, double v_min, double v_max)
{
	//Value of x for which the pdf is maximum
	double x_max = sqrt(3.0)*v_rel;
	//Total area under comparison function
	double area = (v_max - v_min) * VMaxwellianPdf(x_max, v_rel);
	while (true){
		double u = randomUniformDoubleOpen()*area;
		double x = VMaxwellianXFromArea(u, v_rel, v_min, x_max);
		double y = randomUniformDoubleOpen() *VMaxwellianComparison(x, v_rel, v_min, v_max, x_max);
		if (y < VMaxwellianPdf(x, v_rel)){
			return x;
		}
	}
}