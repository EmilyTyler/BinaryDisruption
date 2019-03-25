#ifndef RANDOM_NUMBERS_H
#define RANDOM_NUMBERS_H

//Indirectly tested through random_direction
long double randomUniformDoubleClosed(long double min, long double max);

//Indirectly tested through random_direction
long double randomUniformDoubleOpen(long double min, long double max);

//Tested
long double randomExponential(long double rate=1.0);

long double randomNormal(long double mean, long double stddev);

int randomPoisson(long double mean);

#endif