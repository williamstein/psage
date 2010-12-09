#ifndef _SMALLJAC_INTERNAL_
#define _SMALLJAC_INTERNAL_

#include "gmp.h"
#include "smalljac.h"
#include "ffpoly.h"

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

// Internal smalljac functions not part of the public interface

#define SMALLJAC_MAX_BAD_PRIMES		32
#define SMALLJAC_MAX_FACTOR_BITS		64						// don't waste time trying to factor anything hard
#define SMALLJAC_SMALL_INTERVAL			1000					// don't bother factoring for small intervals

#define SMALLJAC_SPECIAL_X6PA			1						// special curve type 1: y^2=x^6+a

#define SMALLJAC_QCURVE_DELTA			0x1						// indicates Delta values have been precomputed
#define SMALLJAC_QCURVE_WS			0x2						// indicates genus 1 curve specified in Weierstrass form
#define SMALLJAC_QCURVE_2G				0x4						// indicates the 2g coefficient was non-zero - modified curve not defined for primes dividing 2g
															// however the original curve may have had good reduction at these primes and should be checked

#define SMALLJAC_GROUP_FLAG			0x1000


typedef struct smalljac_Qcurve_struct {
	mpz_t f[SMALLJAC_MAX_DEGREE+1];							// Integer poly f(x) representing numerator of rational poly over common denominator
	mpz_t Deltas[SMALLJAC_MAX_DEGREE+1];						// Deltas[k] = (\Delta^k f)(0), the k-th difference poly of f evaluated at 0
	mpz_t denom;												// common denominator - applys to both f and Deltaf0
	mpz_t disc;												// mostly D is used, but its handy to keep this around
	mpz_t D;													// discriminant * denom  - the divisors of D are bad primes (curve undefined or singular)
	int f_inits;												//# of f elements initiailized
	int Delta_inits;												// # of Deltaf0 elements initialized
	int degree;												// degree of the curve - currently always 2*genus + 1
	int genus;												// genus of the curve
	unsigned flags;											// mask of boolean flags indicating the state of initialization of other parameters
	unsigned ws2;											// for Weierstrass specified curves, holds the coefficients mod 2 in bottom 5 bits
	unsigned ws3;											// for Weierstrass specified curves, holds the coefficients mod 3 in bottom 10 bits (2 bits per)
	long a[SMALLJAC_MAX_GENUS];								// L_p(T) coefficients (or group structure coefficients) most recently computed
	int n;													// number of coefficients
	int special;												// flag for special curves, e.g. x^6+a
	unsigned long pts;										// pointcount over F_p, if performed, zero o.w.
	unsigned long p;											// prime for which a[] are the coefficients of L_p(T)
	char str[1024];											// pointer to original curve specification - null terminated
} smalljac_Qcurve;

int padic_charpoly(long a[], long f[], int n, unsigned long p);			// c -> c++ interface function for David Harvey's frobenius() code - used in genus 3

int smalljac_internal_Lpoly (long a[], smalljac_Qcurve *qc, unsigned long p, unsigned long flags);			// doesn't handle small primes, assumes good reduction
unsigned long smalljac_pointcount (smalljac_Qcurve *qc,  unsigned long p);
int smalljac_padic_Lpoly (long a[], smalljac_Qcurve *qc, unsigned long p, unsigned long flags);
int smalljac_generic_Lpoly (long a[], smalljac_Qcurve *qc, unsigned long p, unsigned long pts, unsigned long flags);
int smalljac_tiny_Lpoly (long a[], smalljac_Qcurve *qc, int p, unsigned long flags);
int smalljac_x6pa_Lpoly (long a[], smalljac_Qcurve *qc, long p, unsigned long flags);
int smalljac_Lpoly_extend (long a[], int n, unsigned long p, int h);		// extend coefficients for prime field p to extension field of size q = p^h			

// the following functions update an existing Qcurve structure to reflect new curve parameters without reallocating.
int smalljac_Qcurve_set_str (smalljac_Qcurve *qc, char *str, int *err);
int smalljac_Qcurve_set_mpz (smalljac_Qcurve *qc, mpz_t f[], int degree, mpz_t denom, mpz_t disc, char *str);	// note str is not validated
int smalljac_Qcurve_set_i (smalljac_Qcurve *qc, long f[], int degree, char *str);							// ditto

smalljac_Qcurve *smalljac_Qcurve_alloc ();							// simply allocates an unitialized Qcurve structure to be used above

static inline void smalljac_Qcurve_init_f (smalljac_Qcurve *qc, int degree)
	{ qc->degree = degree;  qc->genus = (degree-1)>>1;  while ( qc->f_inits <= qc->degree ) mpz_init (qc->f[qc->f_inits++]); }

static inline void smalljac_Qcurve_init_Deltas(smalljac_Qcurve *qc)
	{ while ( qc->Delta_inits <= qc->degree ) mpz_init (qc->Deltas[qc->Delta_inits++]); }

static inline int smalljac_Qcurve_reduce (ff_t f[], smalljac_Qcurve *qc)
	{ return ff_poly_set_rational_mpz (f, qc->f, qc->degree, qc->denom); }

static inline int smalljac_Qcurve_degree (smalljac_Qcurve *qc)
	{ return qc->degree; }

static inline int smalljac_Qcurve_special (smalljac_Qcurve *qc)
	{ return qc->special; }

	
void smalljac_init (void);

// computes L_p(1) - does not check for overflow
static inline unsigned long smalljac_Lp1_ui (long a[], int genus, unsigned long p)
{
	switch (genus) {
	case 1: return (unsigned long)((long)(p+1)+a[0]);
	case 2: return (unsigned long)((long)(p*p+1)+(long)(p+1)*a[0]+a[1]);
	case 3: return (unsigned long)((long)(p*p*p+1)+(long)(p*p+1)*a[0]+(long)(p+1)*a[1]+a[2]);
	default: printf ("unhandled genus %d\n", genus);  exit (0);
	}
}

static inline int smalljac_last_Lpoly (long a[], smalljac_Qcurve *qc)
	{ register int i;  for ( i = 0 ; i < qc->n ; i++ ) a[i] = qc->a[i];  return qc->n; }
static inline unsigned long smalljac_last_p (smalljac_Qcurve *qc)
	{ return qc->p; }
static inline unsigned long smalljac_last_pts (smalljac_Qcurve *qc)
	{ return qc->pts; }
	
#endif
