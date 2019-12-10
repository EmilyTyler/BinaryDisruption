// Constants
#ifndef CONSTANTS_H
#define CONSTANTS_H
#include <cmath>

const long double pi = 3.14159265358979323846L;
const long double G_SI = 6.67430L * pow(10.0L, -11.0L);
const long double giga = pow(10.0L, 9.0L);

const long double year = 365.25L*24L*60L*60L;
const long double au = 149597870700L;
const long double parsec = 3.0856775813057292L * pow(10.0L, 16.0L);
const long double msol = 2.0L * pow(10.0L, 30.0L);

//Internal units
const long double time_scale = giga*year;
//const long double time_scale = 1.0;
const long double mass_scale = msol;
//const long double mass_scale = 1.0;
const long double length_scale = pow((G_SI * mass_scale * time_scale*time_scale), 1.0L/3.0L);
//const long double length_scale = 1.0;
const long double G = 1.0L;
//const long double G = G_SI;

//Random number generator seed
//const long double seed = 3263214825;
const long double seed = 2137850808;


#endif