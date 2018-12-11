#ifndef RANDOM_NUMBERS_H
#define RANDOM_NUMBERS_H

//Indirectly tested through random_direction
double randomUniformDoubleClosed(double min=0.0, double max=1.0);

//Indirectly tested through random_direction
double randomUniformDoubleOpen(double min=0.0, double max=1.0);

//Tested
double randomExponential(double rate=1.0);

#endif