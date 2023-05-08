//Functions to draw a random number from a v times maxwellian distribution
#ifndef MC_VELOCITY_H
#define MC_VELOCITY_H


long double drawVMaxwellian(long double v_rel, long double v_max);


long double VMaxwellianPdf(long double x, long double v_rel);


long double VMaxwellianComparison(long double x, long double v_rel, long double v_max, long double f_max);


long double VMaxwellianXFromArea(long double A, long double v_rel, long double f_max);


long double drawMaxwellian(long double v_rel);

#endif