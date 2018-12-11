//Functions to draw a random number from a v times maxwellian distribution
#ifndef MC_VELOCITY_H
#define MC_VELOCITY_H

//Tested
double drawVMaxwellian(double v_rel, double v_min, double v_max);

//Indirectly tested through drawVMaxwellian
double VMaxwellianPdf(double x, double v_rel);

//Indirectly tested through drawVMaxwellian
double VMaxwellianComparison(double x, double v_rel, double v_min, double v_max);

//Indirectly tested through drawVMaxwellian
double VMaxwellianXFromArea(double A, double v_rel, double v_min);

#endif