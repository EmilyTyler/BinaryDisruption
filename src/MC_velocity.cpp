//Functions to draw a random number from a v times maxwellian distribution
#include <cmath>
#include "constants.h"
#include "random_numbers.h"
using namespace std;


long double VMaxwellianPdf(long double x, long double v_rel)
{
	return pow(x,3.0L)/(2.0L*pow(v_rel,4.0L))*exp(-x*x/(2.0L*v_rel*v_rel));
}


long double VMaxwellianComparison(long double x, long double v_rel,long double v_max, long double f_max)
{
	if(x < v_rel){
		return pow(x/v_rel, 3.0L)*f_max;
	} else{
		if(x < 2.0L*v_rel) {
			return f_max;
		} else {
			return f_max * exp(-(x-2.0L*v_rel)/(sqrt(2.0L)*v_rel));
		}
	}
}


long double VMaxwellianXFromArea(long double A, long double v_rel, long double f_max){
	if (A < f_max*v_rel/4.0L){
		return pow((4.0L*pow(v_rel, 3.0L)*A)/f_max, 0.25L);
	} else {
		if (A < f_max*v_rel*5.0L/4.0L){
			return A/f_max + 3.0L*v_rel/4.0L;
		} else {
			return 2.0L*v_rel + sqrt(2.0L)*v_rel*log(sqrt(2.0L)) - sqrt(2.0L)*v_rel*log(-A/(f_max*v_rel) + 5.0L/4.0L + sqrt(2.0L));
		}
	}
}


long double drawVMaxwellian(long double v_rel, long double v_max)
{
	//Maximum of pdf
	long double f_max = VMaxwellianPdf(sqrt(3.0L)*v_rel, v_rel);
	//Total area under comparison function
	long double area = f_max*v_rel*(5.0L/4.0L + sqrt(2.0L) - sqrt(2.0L)*exp(sqrt(2.0L) - v_max/v_rel/(sqrt(2.0L))));
	while (true){
		long double u = randomUniformDoubleOpen(0.0L, area);
		long double x = VMaxwellianXFromArea(u, v_rel, f_max);
		long double y = randomUniformDoubleOpen(0.0L, VMaxwellianComparison(x, v_rel, v_max, f_max));
		if (y < VMaxwellianPdf(x, v_rel)){
			return x;
		}
	}
}

long double drawMaxwellian(long double v_rel)
{
	long double v_x = randomNormal(0.0L, v_rel);
	long double v_y = randomNormal(0.0L, v_rel);
	long double v_z = randomNormal(0.0L, v_rel);
	return sqrt(v_x*v_x + v_y*v_y + v_z*v_z);
}