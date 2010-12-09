#ifndef _POINTCOUNT_INCLUDE_
#define _POINTCOUNT_INCLUDE_

#include "ff.h"

/*
    Copyright 2007 Andrew V. Sutherland

    This file is part of smalljac.

    smalljac is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    smalljac is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with smalljac.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
	This module implements pointcounting over prime fields based on finite differences,
	as described in [KedlayaSutherland2007]
*/

#define POINTCOUNT_MULTI_X		32

// Allocates memory for map of residues - must call first
void pointcount_init (unsigned maxp);

// Note that precompute returns integer values, independent of p
// These need to be reduced mod p prior to calling any of the functions below
void pointcount_precompute (mpz_t D[], mpz_t f[], int degree);
void pointcount_precompute_long (long D[], long f[], int degree);

// handy inline for doing reduction (note sign handling, d[i] must have nonnegative values)
static inline void pointcount_reduce (unsigned long d[], long D[], int degree, long p)
	{ register long x; register int i; for ( i = 0 ; i <= degree ; i++ ) { x=D[i]%p;  d[i] = (x < 0 ? x+p : x); } }

// Standard hyperelliptic point counting functions for curves of the form y^2=f(x)  - default is degree 2g+1
unsigned pointcount_g1 (unsigned long D[4], unsigned p);
unsigned pointcount_g2 (unsigned long D[6], unsigned p);
unsigned pointcount_g2d6 (unsigned long D[7], unsigned p, unsigned long f6);
unsigned pointcount_g3 (unsigned long D[8], unsigned p);
unsigned pointcount_g3d8 (unsigned long D[9], unsigned p, unsigned long f8);
unsigned pointcount_g4 (unsigned long D[10], unsigned p);

// Point counting over Picard curves y^3 = f(x)
unsigned pointcount_pd4 (unsigned long D[5], unsigned p);

// Point counting in F_p^2
unsigned pointcount_g2_d2 (unsigned long f[6], unsigned p);

// This use half as much memory but are slower - use when p/8 > L2 cache
unsigned pointcount_big_g1 (unsigned long D[4], unsigned p);
unsigned pointcount_big_g2 (unsigned long D[6], unsigned p);
unsigned pointcount_big_g2d6 (unsigned long D[7], unsigned p, unsigned long f6);
unsigned pointcount_big_g3 (unsigned long D[8], unsigned p);
unsigned pointcount_big_g4 (unsigned long D[10], unsigned p);

// These routines return point counts on 32 curves f(x), f(x)+1, ..., f(x)+32
int pointcount_multi_g2 (unsigned pts[], unsigned long D[6], unsigned p);
int pointcount_multi_g2d6 (unsigned pts[], unsigned long D[6], unsigned p, unsigned long f6);
int pointcount_multi_g3 (unsigned pts[], unsigned long D[8], unsigned p);
int pointcount_multi_g3d8 (unsigned pts[], unsigned long D[6], unsigned p, unsigned long f8);

// naive polynomial evaluation, provided for comparison
unsigned pointcount_slow (ff_t f[], int d, unsigned p);
unsigned pointcount_tiny (unsigned long f[], int d, unsigned p);

// naive polynomial evaluation for pointcounting over F_p^2 - only used for small p
unsigned pointcount_slow_d2 (ff_t f[], int d, unsigned p);
unsigned pointcount_tiny_d2 (unsigned long f[], int d, unsigned p);

// naive polynomial evaluation for pointcounting over F_p^3 - only used for small p
unsigned pointcount_tiny_d3 (unsigned long f[], int d, unsigned p);

#endif
