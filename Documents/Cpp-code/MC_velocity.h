//Functions to draw a random number from a v times maxwellian distribution
#ifndef MC_VELOCITY_H
#define MC_VELOCITY_H

double drawVMaxwellian(double v_rel, double v_min, double v_max);

double VMaxwellianPdf(double x, double v_rel);

double VMaxwellianComparison(double x, double v_rel, double v_min, double v_max);

double VMaxwellianXFromArea(double A, double v_rel, double v_min);

#endif