#ifndef RANDOM_NUMBERS_H
#define RANDOM_NUMBERS_H

long double randomUniformDoubleClosed(long double min, long double max);

long double randomUniformDoubleOpen(long double min, long double max);

long double randomExponential(long double rate=1.0);

long double randomNormal(long double mean, long double stddev);

int randomPoisson(long double mean);

#endif