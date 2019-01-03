#ifndef RANDOM_NUMBERS_H
#define RANDOM_NUMBERS_H

//Indirectly tested through random_direction
long double randomUniformDoubleClosed();

//Indirectly tested through random_direction
long double randomUniformDoubleOpen();

//Tested
long double randomExponential(long double rate=1.0);

#endif