// written by Szymon Winczewski

#ifndef constants_h
#define constants_h

#include <cmath>


// liczba PI oraz jej rozne wielokrotnosci
const double HALF_OF_PI = 2.0 * atan(1.0);
const double PI = 4.0 * atan(1.0);
const double THREE_HALFS_OF_PI = 6.0 * atan(1.0);
const double TWO_PI = 8.0 * atan(1.0);
const double FOUR_PI = 16.0 * atan(1.0);
const double INV_OF_FOUR_PI = 1.0 / ( 16.0 * atan(1.0) );
const double SQRT_OF_PI = sqrt(4.0 * atan(1.0));

// konwersja radiany-stopnie oraz stopnie-radiany
const double RAD_TO_DEG = 180.0 / PI;
const double DEG_TO_RAD = PI / 180.0;

#endif
