//Functions to draw a random number from a v times maxwellian distribution
#ifndef MC_VELOCITY_H
#define MC_VELOCITY_H

//Tested
long double drawVMaxwellian(long double v_rel, long double v_min, long double v_max);

//Indirectly tested through drawVMaxwellian
long double VMaxwellianPdf(long double x, long double v_rel);

//Indirectly tested through drawVMaxwellian
long double VMaxwellianComparison(long double x, long double v_rel, long double v_min, long double v_max);

//Indirectly tested through drawVMaxwellian
long double VMaxwellianXFromArea(long double A, long double v_rel, long double v_min);

#endif