#include <stdlib.h>
#include <stdio.h>
#include "ff.h"
#include "ffpoly.h"
#include "g2tor3poly.h"

#define MAX_POWER			19

struct poly_4coeff *g2tor3[41] = { g2tor3d0, g2tor3d1, g2tor3d2, g2tor3d3, g2tor3d4, g2tor3d5, g2tor3d6, g2tor3d7, g2tor3d8, g2tor3d9,
g2tor3d10, g2tor3d11, g2tor3d12, g2tor3d13, g2tor3d14, g2tor3d15, g2tor3d16, g2tor3d17, g2tor3d18, g2tor3d19,	
g2tor3d20, g2tor3d21, g2tor3d22, g2tor3d23, g2tor3d24, g2tor3d25, g2tor3d26, g2tor3d27, g2tor3d28, g2tor3d29,
g2tor3d30, g2tor3d31, g2tor3d32, g2tor3d33, g2tor3d34, g2tor3d35, g2tor3d36, g2tor3d37, g2tor3d38, g2tor3d39,g2tor3d40};

void ff_poly_g2tor3_modpoly (ff_t g[POLY_G2TOR3_DEGREE+1], ff_t f[6])
{
	register int i, j, k;
	ff_t powers[4][20];
	ff_t t0,c;
	
	// compute coefficient powers
	for ( i = 0 ; i < 4 ; i++ ) {
		_ff_set_one(powers[i][0]);
		_ff_set (powers[i][1], f[i]);
		for ( j = 2 ; j <= MAX_POWER ; j++ ) _ff_mult(powers[i][j],powers[i][j-1],f[i]);
	}
	
	/*
		The 3-torsion modular poly provided by Gaudry & Schost assumes f(x) = x^5+f0x^3+f1x^2+f2x + f3, while we use x^5+f3x^3+f2x^2+f1x+f0.
		Thus the need for replacing k with 3-k below.
	*/
	for ( i = 0 ; i <= POLY_G2TOR3_DEGREE ; i++ ) {
		_ff_set_zero(c);
		for ( j = 0 ; g2tor3[i][j].c ; j++ ) {
			_ff_set_i (t0,g2tor3[i][j].c);
			for ( k = 0 ; k < 4 ; k++ ) ff_mult(t0,t0,powers[3-k][g2tor3[i][j].f[k]]);
			_ff_addto(c,t0);
		}
		_ff_set(g[i],c);
	}
}
