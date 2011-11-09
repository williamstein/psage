#ifndef __MATHLIB_H_
#define __MATHLIB_H_

#include <complex>

std::complex<double> log_GAMMA(std::complex<double> z, int n = 0);
double sinc(double u, double Delta);
double sinc_square(double u, double Delta);
double triangle(double u, double Delta);

#endif
