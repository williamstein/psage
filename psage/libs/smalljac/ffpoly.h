#ifndef _FFPOLY_INCLUDE
#define _FFPOLY_INCLUDE

#include "gmp.h"
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

// This module consists mostly of functions that are probably better implemented in NTL
// but are handy for a standalone C implementation.

// The discriminant and resultant code is woefully incomplete

#define POLY_MAX_DEGREE			80
#define POLY_G2TOR3_DEGREE		40

// handy macros - these require integer variables d_f be declared for each poly f
#define _poly_set(a,b)				ff_poly_copy(a,&d_ ## a,b,d_ ## b)
#define _poly_print(a)				ff_poly_print(a,d_ ## a)
#define _poly_neg(c,a)				ff_poly_neg(c,&d_ ## c,a,d_ ## a)
#define _poly_monic(c,a)				ff_poly_monic(c,&d_ ## c,a,d_ ## a)
#define _poly_addto(c,a)				ff_poly_addto(c,&d_ ## c,a,d_ ## a)
#define _poly_add(c,a,b)				ff_poly_add(c,&d_ ## c,a,d_ ## a,b,d_ ## b)
#define _poly_sub(c,a,b)				ff_poly_sub(c,&d_ ## c,a,d_ ## a,b,d_ ## b)
#define _poly_mult(c,a,b)			ff_poly_mult(c,&d_ ## c,a,d_ ## a,b,d_ ## b)
#define _poly_div(q,r,a,b)			ff_poly_div(q,&d_ ## q,r,&d_ ## r,a,d_ ## a,b,d_ ## b)
#define _poly_mod(c,a,b)			ff_poly_mod(c,&d_ ## c,a,d_ ## a,b,d_ ## b)
#define _poly_gcd(c,a,b)				ff_poly_gcd(c,&d_ ## c,a,d_ ## a,b,d_ ## b)
#define _poly_gcdext(c,u,v,a,b)		ff_poly_gcdext(c,&d_ ## c,u,&d_ ## u, v,&d_ ## v,a,d_ ## a,b,d_ ## b)
#define _poly_expmod(c,a,e,b)		ff_poly_exp_mod(c,&d_ ## c,a,&d_ ## a,e,b,&d_ ## b)

int poly_expr_degree (char *str, int *ws);
void mpq_poly_weierstrass (mpq_t f[4], mpz_t W[5]);
int mpq_poly_expr (mpq_t f[], int maxd, char *expr);				// handles Weierstrass or standard form
int mpq_poly_standardize (mpq_t f[], int d);						// make d-1 coefficient zero
int mpq_poly_jinv_i (mpq_t j, long a[5]);

int mpz_poly_expr (mpz_t f[], int maxd, char *expr);
int mpz_poly_sprint (char *s, mpz_t f[], int d_f);
void mpz_poly_print (mpz_t f[], int d_f);
void mpz_poly_discriminant (mpz_t disc, mpz_t f[], int d);
void mpq_poly_discriminant (mpq_t disc, mpq_t f[], int d);
int ui_poly_expr (unsigned long f[], int maxd, char *expr);
int ui_poly_expr_mod_p (unsigned long f[], int maxd, char *expr, unsigned long p);
int i_poly_expr (long f[], int maxd, char *expr);
int ui_poly_sprint (char *s, unsigned long f[], int d_f);
void ui_poly_print (unsigned long f[], int d_f);
int i_poly_sprint (char *s, long f[], int d_f);
void i_poly_print (long f[], int d_f);

unsigned long i_poly_set_mpq (long f[], mpq_t F[], int d);			// sets f[i] to numberator of F[i] with common denominator and returns denominator or zero if overflow
void mpz_poly_set_mpq (mpz_t denom, mpz_t f[], mpq_t F[], int d);	// computes integer denom and integer poly f s.t. f/denom = F
int ff_poly_set_rational_mpz (ff_t f[], mpz_t F[], int d, mpz_t denom);
int ui_poly_set_rational_mpz_mod_p (unsigned long f[], mpz_t F[], int d, mpz_t denom, unsigned long p);

int ff_poly_expr (ff_t f[], int maxd, char *expr);
void ff_poly_print (ff_t f[], int d_f);
int ff_poly_sprint (char *s, ff_t f[], int d_f);

void ff_poly_eval (ff_t y[1], ff_t f[], int d, ff_t x[1]);
void ff_poly_twist_mpz (ff_t g[], mpz_t f[], int d);
void ff_poly_twist (ff_t g[], ff_t f[], int d);
int ff_poly_irreducible (ff_t f[], int d_f, int *pnroots);	// pnroots is the number of distinct roots
int ff_poly_factors (ff_t f[], int d_f);					// number of (not-necessarily distinct) irreducible factors
int ff_poly_roots (ff_t f[], int d_f);
void ff_poly_derivative (ff_t g[], int *pd_g, ff_t f[], int d_f);
void ff_poly_discriminant (ff_t disc[1], ff_t f[], int d);
void ff_poly_exp_mod (ff_t g[], int *pd_g, ff_t a[], int d_a, mpz_t e, ff_t f[], int d_f);
void ff_poly_gcd (ff_t g[], int *pd_g, ff_t a[], int d_a, ff_t b[], int d_b);
void ff_poly_gcd_red (ff_t g[], int *pd_g, ff_t a[], int d_a, ff_t b[], int d_b);
void ff_poly_gcdext (ff_t D[], int *pd_D, ff_t u[], int *pd_u, ff_t v[], int *pd_v, ff_t a[], int d_a, ff_t b[], int d_b);
void ff_poly_mod (ff_t g[], int *pd_g, ff_t a[], int d_a, ff_t f[], int d_f);
void ff_poly_div (ff_t q[], int *pd_q, ff_t r[], int *pd_r, ff_t a[], int d_a, ff_t b[], int d_b);
void ff_poly_neg (ff_t c[], int *pd_c, ff_t a[], int d_a);
void ff_poly_addto (ff_t c[], int *pd_c, ff_t a[], int d_a);
void ff_poly_add (ff_t c[], int *pd_c, ff_t a[], int d_a, ff_t b[], int d_b);
void ff_poly_sub (ff_t c[], int *pd_c, ff_t a[], int d_a, ff_t b[], int d_b);
void ff_poly_mult (ff_t c[], int *pd_c, ff_t a[], int d_a, ff_t b[], int d_b);
void ff_poly_monic (ff_t c[], int *pd_c, ff_t b[], int d_b);

int mpz_poly_bad_primes (unsigned long bad[], int maxbad, mpz_t f[], int d);
int mpq_poly_bad_primes (unsigned long bad[], int maxbad, mpq_t f[], int d);

void ff_depress_cubic (ff_t t[1], ff_t f[4]);					// puts monic f into the form x^3+ax+b via translation f(x-t)
void ff_depressed_cubic_resolvent (ff_t t[1], ff_t g[4], ff_t f[5]);	// forms cubic resolvent of depressed quartic and depresses it via translation t
int ff_poly_roots_d2 (ff_t r[2], ff_t f[3], int d_f);				// f may be quadratic or linear, needn't be monic
int _ff_poly_roots_d3 (ff_t r[3], ff_t f[4], ff_t *pd);				// f must be of the form x^3+ax+b, r may be null if only a count is desired
static inline int ff_poly_roots_d3 (ff_t r[3], ff_t f[4]) { return _ff_poly_roots_d3 (r, f, 0); }
int ff_old_poly_roots_d3 (ff_t r[3], ff_t f[4]);
int _ff_poly_roots_d4 (ff_t r[4], ff_t f[5], ff_t *pd);
static inline int ff_poly_roots_d4 (ff_t r[3], ff_t f[4]) { return _ff_poly_roots_d4 (r, f, 0); }
int ff_old_poly_roots_d4 (ff_t r[4], ff_t f[5]);
int ff_poly_g1_3tor (ff_t f[4]);								// f must be of the form x^3+ax+b
int ff_poly_g1_4tor (int *o, ff_t f[4], int flag8);				// f must be of the form x^3+ax+b (flag8 set indicates that Z/8Z should also be checked)
int ff_poly_g2_3tor (ff_t f[5]);								// f must be of the form x^5+ax+b
void ff_poly_xn_mod_d3 (ff_t g[3], unsigned long n, ff_t f[4]);	// computes x^n mod f=x^3+ax+b
void ff_poly_xn_mod_d4 (ff_t g[3], unsigned long n, ff_t f[5]);
void ff_poly_xpan_mod_d3 (ff_t g[3], ff_t a, unsigned long n, ff_t f[4]);
void ff_poly_xpan_mod_d4 (ff_t g[4], ff_t a, unsigned long n, ff_t f[5]);
void ff_poly_g2tor3_modpoly (ff_t g[41], ff_t f[6]);

int ff_poly_g1_2Sylow (ff_t f[4]);							// returns size of the 2-Sylow subgroup of the elliptic curve y^2=f(x)=x^3+f1x+f0, provided it is cyclic, returns -1 otherwise.

static inline void ff_poly_copy (ff_t b[], int *pd_b, ff_t a[], int d_a)
{
	register int i;
	
	for ( i = 0 ; i <= d_a ; i++ ) _ff_set(b[i], a[i]);
	*pd_b = d_a;
}

static inline int ff_poly_equal (ff_t a[], int d_a, ff_t b[], int d_b)
{
	register int i;
	
	if ( d_a != d_b ) return 0;
	for ( i = 0 ; i <= d_a ; i++ ) if ( ! _ff_equal(a[i],b[i]) ) return 0;
	return 1;
}

static inline void ff_poly_set_i (ff_t f[], long F[], int d) { register int i;  for ( i = 0 ; i <= d ; i++ ) _ff_set_i(f[i],F[i]); }

#endif
