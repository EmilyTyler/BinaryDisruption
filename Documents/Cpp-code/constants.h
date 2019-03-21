// Constants
#ifndef CONSTANTS_H
#define CONSTANTS_H
#include <cmath>

const long double pi = 3.14159265359;
const long double G_SI = 6.67408 * pow(10.0, -11.0);
const long double giga = pow(10.0, 9.0);

const long double year = 365.25*24*60*60;
const long double au = 149597870700;
const long double parsec = 3.0856775813057292 * pow(10.0, 16.0);
const long double msol = 2.0 * pow(10.0, 30.0);

//Internal units
const long double time_scale = giga*year;
//const long double time_scale = 1.0;
const long double mass_scale = msol;
//const long double mass_scale = 1.0;
const long double length_scale = pow((G_SI * mass_scale * time_scale*time_scale), 1.0/3.0);
//const long double length_scale = 1.0;
const long double G = 1.0;
//const long double G = G_SI;

//Random number generator seed
//const long double seed = 3263214825;
const long double seed = 2137850808;


#endif