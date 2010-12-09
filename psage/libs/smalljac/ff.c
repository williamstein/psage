#include <stdlib.h>
#include <stdio.h>
#include "ff.h"

/*
    Copyright 2007-2008 Andrew V. Sutherland

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


ff_t _ff_t1;
ff_t _ff_p;
ff_t _ff_2g;		// generator of Sylow 2-subgroup and necessarily a quadratic non-residue
ff_t _ff_2gi;		// inverse of _ff_2g
ff_t _ff_3g;		// generator of Sylow 3-
ff_t _ff_negone;
ff_t _ff_half;
ff_t _ff_third;
static unsigned long mod256itab[256];
static ff_t __ff_mont_itab[FF_ITAB_SIZE];				// private copy for single modulus environment

static ff_t _ff_2exp;								// 2^(p+3)/4-1 set only for p=5mod8, used for fast sqrts
unsigned long _ff_p2_m;						// odd part of p-1, p=2^p2_e*p2_m+1
int _ff_p2_e;									// power of 2 dividing p-1
unsigned long _ff_p3_m;						// p=3^p3_e*p3_m+1 when p=1mod3, p=3^p3_3*p3_m-1 when p=2mod3
int _ff_p3_m1mod3;								// true if p3_m is 1 mod 3
int _ff_p3_e;									// power of 3 dividing  p-1 when p=1mod3, power of 3 dividing p+1 when p=2mod3
ff_t _ff_2Sylow_tab[64][2];						// 2Sylow_tab[i][0]=g^(2^i) for 0<=i<=e-1, where g generates the 2-Sylow subgroup.
											// 2Sylow_tab[i][1] = 2Sylow_tab[i][0]*2Sylow_tab[i+1][0]
static ff_t _ff_3Sylow_tab[42][2];					// 3Sylow_tab[i][0]=g^(3^i) for 0<=i<=e-1, where g generates the 3-Sylow subgroup (when it is non-trivial)
											// 3Sylow_tab[i][1] = 3Sylow_tab[i][0]^2
ff_t _ff_cbrt_unity;
int _ff_sqrt_setup;
int _ff_cbrt_setup;

int _ff_p1mod3;

ff_t *_ff_mont_itab = __ff_mont_itab;				// my be repointed to mange multiple moduli
ff_t _ff_mont_R;
ff_t _ff_mont_R2;
unsigned long _ff_mont_pni;


void ff_montgomery_setup (int ring);
int ff_ui_legendre (unsigned long a, unsigned long b);

void ff_ext_setup(void);

// nr35_tab[p%35] is either 0 or a non-residue mod p for p=1mod4 (if p is 2 or 3 mod 5, then 5 is a non-residue mod p; if p is 3, 5, or 6 mod 7, then 7 is a non-residue mod p).
// nric35_tab[p%35] is, the coefficient k s.t. kp+1=0 mod ff_nr35_tab[p%35]
int _ff_nr35_tab[35] =    { 0, 0, 5, 5, 0, 7, 7, 5, 5, 0, 7, 0, 5, 5, 0, 0, 0, 5, 5, 7, 7, 0, 5, 5, 7, 0, 7, 5, 5, 0, 0, 7, 5, 5, 7};
int _ff_nric35_tab[35] = { 0, 0, 2, 3, 0, 4, 1, 2, 3, 0, 2, 0, 2, 3, 0, 0, 0, 2, 3, 4, 1, 0, 2, 3, 2, 0, 4, 2, 3, 0, 0, 2, 2, 3, 1};


// Note that it is possible to use non-prime p, and we don't want to waste time checking primality anyway.
// We assume that either the caller has checked the primality of p, or is aware that it is working over a ring
void ff_setup_ui (unsigned long p)
{
	static int init = 0;
	unsigned long n;
	ff_t t;
	
#if  FF_WORDS != 1
	err_printf ("ff_setup_ui only supports single-word Montgomery or native representation\n");
	exit (0);
#endif
#if FF_MONTGOMERY
	if ( ! init ) { for ( n = 1 ; n < 256 ; n+= 2 ) mod256itab[n] = ff_ui_inverse(n,256); init = 1; }		// This is only done once
#endif

	if ( p && _ff_p == p ) return;
	if ( p > (1UL<<FF_BITS) ) { printf ("ff_setup: p=%lu is too large, FF_BITS = %d\n", p, FF_BITS);  exit (0); }
	
	// Don't deal with char 2
	if ( p < 3 ) { printf ("ff_setup: p=%lu must be an odd prime\n", p);  exit (0); }

	_ff_p = p;
#if FF_MONTGOMERY
	ff_montgomery_setup (0);
#endif
	_ff_p2_m = (_ff_p-1)/2;
	_ff_set_ui (_ff_half, _ff_p2_m+1);
	_ff_set_zero(_ff_negone);  _ff_dec(_ff_negone);
	for ( _ff_p2_e = 1 ; !(_ff_p2_m&0x1) ; _ff_p2_m>>=1, _ff_p2_e++ );
	_ff_sqrt_setup = _ff_cbrt_setup = 0;

	_ff_p1mod3 = ((p%3)==1);
	n = (_ff_p1mod3? (2*p+1UL)/3UL : (p==3?0:(p+1UL)/3UL));
	_ff_set_ui(_ff_third,n);				// note this is set to zero for p=3
	if ( _ff_p1mod3 ) {
		_ff_p3_m = (_ff_p-1)/3;
		for ( _ff_p3_e = 1 ; !(_ff_p3_m%3) ; _ff_p3_m /= 3, _ff_p3_e++ );
		_ff_p3_m1mod3 = ((_ff_p3_m%3)==1);
	} else {
		// this code is meaningless if p=3, but we won't bother checking this because these values won't be used in this case
		_ff_p3_m = (_ff_p+1)/3;
		for ( _ff_p3_e = 1 ; !(_ff_p3_m%3) ; _ff_p3_m /= 3, _ff_p3_e++ );
		_ff_p3_m1mod3 = ((_ff_p3_m%3)==1);
	}
	
	_ff_2exp = _ff_2g = _ff_3g = 0;			// will get set when needed
	ff_ext_setup();
	return;
}

// macros to handle setting ff_t types in raw mode -- needed to handle multiples of p and montgomery values
#if FF_WORDS == 1
#define _ff_raw_set_mpz(z,X)		((z) = mpz_get_ui(X))
#define _ff_raw_get_mpz(Z,x)		(mpz_set_ui(Z,x),Z)
#endif


#if FF_WORDS == 1 && FF_MONTGOMERY

void ff_montgomery_setup (int ring)
{
	register unsigned long p, t;
	register int i;
	
	// precompute -p^{-1} mod B (B is either 2^32 or 2^64), R = B mod p and R^2 mod p
	// Use Jebelean's trick for computing p^{-1} mod B (see HECHECC p. 190 remark 10.4 (ii))
	p = _ff_p;
	t = mod256itab[p&0xFF];
	t = (2*t + (-(p*t*t)))&0xFFFF;
	t = (2*t + (-(p*t*t)))&0xFFFFFFFF;
#if FF_HALF_WORDS == 1
	_ff_mont_pni = (-t)&0xFFFFFFFF;
	t = (1UL<<FF_MONTGOMERY_RBITS)%p;
	_ff_mont_R = (ff_t)t;
	t = (t*t)%p;
	_ff_mont_R2 = (ff_t)t;
#else
	t = (2*t + (-(p*t*t)))&0xFFFFFFFFFFFFFFFF;
	_ff_mont_pni = -t;
	_ff_mont_R = (1UL<<(FF_MONTGOMERY_RBITS-1))%p;
	_ff_mont_R += _ff_mont_R;
	if ( _ff_mont_R > p ) _ff_mont_R -= p;
		
	// to avoid depending on multi-precision arithmetic, we square R mod p by doubling it RBITS times (R = 2^RBITS)
	// probably no slower than a modular reduction anyway
	_ff_mont_R2 = _ff_mont_R;
	for ( i = 0 ; i < FF_MONTGOMERY_RBITS ; i++ ) {
		_ff_mont_R2 += _ff_mont_R2;
		if ( _ff_mont_R2 > p ) _ff_mont_R2 -= p;
	}
#endif
	if ( ring ) return;
	t = 1;
	for ( i = 0 ; i <= 3*FF_MONTGOMERY_RBITS ; i++ ) {
		_ff_mont_itab[i] = t;
		t += t;													// note that _ff_p < 2^(ULONG_BITS-1) so overflow can't occur
		if ( t >= _ff_p ) t -= _ff_p;
	}
//	printf ("Montgomery setup p = %lu, -p^{-1} mod b = %lu, R mod p = %lu, R^2 mod p = %lu\n", _ff_p, _ff_mont_pni, _ff_mont_R, _ff_mont_R2);
}

// This is a streamlined version of Algorithm 11.12 in HECHECC p. 208
// The input is in Montgomery representation: x = [a] = aR mod p
// This functions computes v = [a^{-1}] = a^{-1}R mod p = x^{-1}R^2 mod p
unsigned long ff_montgomery1_invert (unsigned long x)
{
	register unsigned long r, s, t, v;
	register int k;

	if ( _ff_zero(x) ) { printf ("ff_invert: attempt to invert 0,\n");exit(0);  return 0; }
	r = x;  s = 1;  t = _ff_p;  v = 0;  k = 0;
	while ( r > 0 ) {
		if ( !(t&0x1UL) ) {
			t >>= 1;  s <<= 1;
		} else if ( !(r&0x1UL) ) {
			r >>= 1;  v <<= 1;
		} else if ( t > r ) {
			t = (t-r)>>1;  v += s;  s <<= 1;
		} else {
			r = (r-t)>>1; s += v;  v <<= 1;
		}
		k++;
//		if ( k > 2*_ff_pbits ) { printf ("ff_invert: k = %d > 2*pbits!  raw input %lu\n", k, x);  exit (0); }
	}
	if ( v >= _ff_p ) v -= _ff_p;				
	v = _ff_p - v;														// v cannot be zero if x is invertible (which we assume to be the case)
	// This accomplishes steps 10-12 of Alg 11.12 in one multiplication
	v = ff_montgomery1_mult (v, _ff_mont_itab[3*FF_MONTGOMERY_RBITS-k]);		// Montgomery lookup table holds powers 2^k to save time
#if ! FF_FAST
	if ( ff_montgomery1_mult (v,x) != _ff_mont_R ) {printf ("ff_montgomery1_invert failed, %lu*%lu != 1, raw input: %lu  raw output %lu\n", _ff_get_ui(x), _ff_get_ui(v), x, v); exit (0); }
#endif
	return v;
}

#endif

// Algorithm 11.15 of [HECHECC], P. 209 - note that z and x can overlap
void ff_parallel_invert (ff_t z[], ff_t x[], unsigned n)
{
	ff_t c[FF_MAX_PARALLEL_INVERTS];
	register ff_t u, v;
	register unsigned i;

	if ( ! n ) return;
	if ( n > FF_MAX_PARALLEL_INVERTS ) { printf ("Exceeded FF_MAX_INVERTS, %d > %d\n", n, FF_MAX_PARALLEL_INVERTS);  exit (0); }
	
	_ff_set (c[0], x[0]);
	for ( i = 1 ; i < n ; i++ ) _ff_mult (c[i], c[i-1], x[i]);
	 _ff_invert (u, c[n-1]);
	for ( i = n-1 ; i > 0 ; i-- ) {
		_ff_mult (v, c[i-1], u);
		_ff_mult (u, u, x[i]);
		_ff_set (z[i], v);
	}
	_ff_set (z[0], u);
}

static inline void ff_setup_2exp(void) { ff_t t; if ( _ff_2exp ) return;  _ff_set_one(t);  _ff_x2(t);  ff_exp_ui (&_ff_2exp,&t,(_ff_p+3)/4-1); }


void ff_setup_2g (void)
{
	ff_t t;
	register int i,n;
	
	if ( _ff_2g ) return;

	// When p=3 mod 4 life is very simple, since -1 then generates the 2-Sylow subgroup
	if ( _ff_p2_e == 1 ) { _ff_set(_ff_2g,_ff_negone);  _ff_set(_ff_2gi,_ff_negone); _ff_set(_ff_2Sylow_tab[0][0],_ff_negone);  return; }
	
	// hardwire a few quick tests for non-residues to speed things up.  This should catch all but 1/32 of the primes, on average.
	_ff_set_zero(t);
	if ( ! _ff_p1mod3 ) { _ff_set_i (t,-3); goto SylowSetup; }					// if p is 2 mod 3 then -3 is a non-residue
	if ( (_ff_p&0x7UL) == 5 ) { _ff_set_ui (t, 2); goto SylowSetup; }			// if p is 5 mod 8 then 2 is a non-residue
	n = _ff_nr35_tab[_ff_p%35];
	if ( n ) { _ff_set_ui (t, n); goto SylowSetup; }
	// We only reach this prob 1/32, and the expected cost of the computation below is less than 32 modular reductions
	// Note that legendre is faster here than computing square roots because n is small
	for ( n = 11 ; ; n++ ) if ( ff_ui_legendre (n, _ff_p) < 0 ) break;
	_ff_set_ui(t,n);
	
SylowSetup:

	ff_exp_ui (&_ff_2g, &t, _ff_p2_m);
	_ff_set(_ff_2Sylow_tab[0][0], _ff_2g);
	_ff_set(_ff_2gi,_ff_2g);
	// compute the inverse of _ff_2g as we square up, this is faster than calling _ff_inverse
	for ( i = 1 ; i < _ff_p2_e-1 ; i++ ) { _ff_square(_ff_2Sylow_tab[i][0],_ff_2Sylow_tab[i-1][0]);  ff_mult(_ff_2gi,_ff_2gi,_ff_2Sylow_tab[i][0]); }
	_ff_set(_ff_2Sylow_tab[i][0],_ff_negone);
	ff_negate(_ff_2gi);
	for ( i = 0 ; i < _ff_p2_e-1 ; i++ ) _ff_mult(_ff_2Sylow_tab[i][1], _ff_2Sylow_tab[i][0], _ff_2Sylow_tab[i+1][0]);
	_ff_set(_ff_2Sylow_tab[i][1],_ff_negone);
}

/*
   Extended sqrt algorithm, based on Tonelli-Shanks, compute sthe square root in F_p^2 if necessary.
   The return value is 0 if the root is in F_p^2, in which case o=a^{-1/2}/z where F_p^2=F_p[z]/(z^2-_ff_2g).

   This algorithm is deterministic (given precomputed 2-Sylow generator) - caller should flip coin to randomize choice of sqrt.
*/
int ff_invsqrt (ff_t o[1], ff_t a[1], int ext)
{
	ff_t b, x, y;
	int sts;

	if ( _ff_zero(a[0]) ) { ff_setup_2g();  _ff_set_zero (o[0]);  return 1; }			// need to handle zero separately (make sure setup_2g gets called though)

	ff_exp_ui (&x, a, (_ff_p2_m-1)/2);			// x = a^{(m-1)/2}
	ff_mult (b, a[0], x);
	ff_mult (b, b, x);						// b = a^m is in the 2-Sylow subgroup
	sts = ff_2Sylow_invsqrt(&y,&b, ext);
	if ( ! ext && ! sts ) return 0;
	ff_mult (x, x, y);						// x = a^{(m-1)/2} * b^{-1/2} = a^{(m-1)/2} * a^{-m/2} = a^{-1/2}
#if ! FF_FAST
	_ff_square (y, x);
	if ( ! sts ) ff_mult(y,y,_ff_2g);
	ff_mult(y,y,a[0]);
	if ( ! _ff_one(y) ) { 
		printf ("p=%lu, a=%lu, b=%lu, 2g=%lu, 2gi=%lu\n", _ff_p, _ff_get_ui(a[0]), _ff_get_ui(b), _ff_get_ui(_ff_2g), _ff_get_ui(_ff_2gi));
		if ( sts ) {
			printf ("ff_sqrt failed, %lu^2 != 1/%lu\n", _ff_get_ui(x),  _ff_get_ui(a[0])); exit (0);
		} else {
			printf ("ff_sqrt failed, %lu^2*%lu != 1/%lu\n", _ff_get_ui(x), _ff_get_ui(_ff_2g), _ff_get_ui(a[0])); exit (0);
		}
	}
#endif
	_ff_set (o[0], x); 
	return sts;
}

/*
   Computes a^{-1/2} for a in the 2-Sylow subgroup, computing in F_p^2 if necessary.
   Returns 0 if the root is in F_p^2, in which case o=a^{-1/2}/z as above.

   Uses precomputed table of 2e powers of the 2-Sylow generator to reduce running time by a factor of 4 over standard Tonelli-Shanks (still O(e^2)).
   Assumes _ff_p2_e > 2 and a has order at least 4 (easy cases handle inline by ff.h)
*/
int _ff_2Sylow_invsqrt (ff_t o[1], ff_t a[1], int ext)
{
	register ff_t b, q, q1, q2, t, w1, w2;
	register int i, j, k, s;
	int sts;

	// set w1 and w2 to the two elements of order 4 in the 2-Sylow subgroup (note we know e > 2 and a has order at least 4 since it is not 1 or -1)
	_ff_set (w1, _ff_2Sylow_tab[_ff_p2_e-2][0]);
	_ff_set (w2, _ff_2Sylow_tab[_ff_p2_e-2][1]);
	_ff_set_one (t);
	_ff_set (b,a[0]);
	sts = 1;
	for(;;) {
		_ff_set (q, b);
		for ( s = 2 ; s<_ff_p2_e+1 ; s++ ) {					// s<_ff_p2_e+1 is an unnecessary safety check (in case a isn't in the 2-Sylow)
			j=1;
			if ( _ff_equal(q,w1) ) break;
			j=0;
			if ( _ff_equal(q,w2) ) break;
			ff_square (q,q);
		}
		k = _ff_p2_e-s;
#if ! FF_FAST
		if ( k < 0 ) { printf ("Unexpected result: k<0 in ff_2Sylow_invsqrt?!  a = %lu, p = %lu\n", _ff_get_ui(a[0]), _ff_p);  exit (1); }
#endif	
		if ( !k ) {											// this can only happen the first time through
			if ( ! ext ) return 0;
			sts=0;
			if ( j ) {
				ff_mult(b,b,_ff_2gi);							// clear the low bit rather than adding
				ff_mult(t,t,_ff_2gi);							// setting sts=0 implicitly multiplied t by s^{1/2}, so multiplying by s^{-1} gives s^{-1/2}
			} else {
				ff_mult(b,b,_ff_2g);							// sts=0 multiplied t by s^{1/2}, so we don't need to do anything to t
			}
		} else {
			ff_mult (b, b, _ff_2Sylow_tab[k][j]);					// the product of all the elements S[k] we multiply into b here is b^{-1}, since we terminate with b=1
			ff_mult (t, t, _ff_2Sylow_tab[k-1][j]);					// multiply t by S[k-1]=sqrt(S[k]), the product of all these will be b^{-1/2}
		}
		if ( _ff_equal(b,_ff_negone) ) { ff_mult(t, t, w1); break; }
		if ( _ff_one(b) ) break;
		ff_mult(q2,q2,_ff_2Sylow_tab[_ff_p2_e-4][j]);
		if ( _ff_one(q2) ) continue;
		if ( _ff_equal(q2,_ff_negone) ) {
			_ff_mult(b,b,_ff_2Sylow_tab[k+3][0]);
			_ff_mult(t,t,_ff_2Sylow_tab[k+2][0]);
		} else {
			if ( _ff_equal(q2,w1) ) j = 1; else j = 0;
			_ff_mult(b,b,_ff_2Sylow_tab[k+2][j]);
			_ff_mult(t,t,_ff_2Sylow_tab[k+1][j]);
		}
		if ( _ff_equal(b,_ff_negone) ) { ff_mult(t, t, w1); break; }
		if ( _ff_one(b) ) break;
	}
check:
#if ! FF_FAST
	_ff_square (b, t);
	if ( ! sts ) ff_mult(b,b,_ff_2g);
	ff_mult(b,b,a[0]);
	if ( ! _ff_one(b) ) {
		printf ("a=%lu, p=%lu\n", _ff_get_ui(a[0]), _ff_p);
		if ( sts ) {
			printf ("ff_2Sylow_invsqrt failed, %lu^2 *  %lu != 1\n", _ff_get_ui(t), _ff_get_ui(a[0]));
		} else {		
			printf ("ff_2Sylow_invsqrt failed, %lu^2 * %lu *  %lu != 1\n", _ff_get_ui(t), _ff_get_ui(_ff_2g), _ff_get_ui(a[0]));
		}
		exit (0); 
	}
#endif
	_ff_set(o[0],t);
	return sts;
}

// Recursive discrete log algorithm - not currently used as it is  slower than ff_2Sylow_invsqrt for p < 2^64
// Given a and b in the 2 Sylow subgroup, returns a nonnegative integer k < 2^e s.t. _a^k = b, or -1 if no such k exists
int ff_2Sylow_dlog (ff_t a, ff_t b, int e)
{
	register ff_t  x0, x1,y;
	register int d, i, k0, k1;
	ff_t t;
//printf ("2Sylow_dlog(%ld,%ld,%d)\n", _ff_get_ui(a), _ff_get_ui(b), e);
	if ( _ff_one(b) ) return 0;
	if ( _ff_equal(b,_ff_negone) ) return (1<<(e-1));			// non-generic optimization specific to fields
	if ( e < 2 ) return -1;
	if ( e == 2 ) {
		if ( _ff_equal (b,a) ) return 1;
		_ff_mult(y,a,b);
		if ( _ff_one(y) ) return 3;
		return -1;
	}
	d = e/2;
	// compute x0=a^(2^(e-d)), y=b^(2^(e-d)), x1 = a^(2^d)
	_ff_set(x0,a); _ff_set(y,b);
	for ( i = 0 ; i < e-d ; i++ ) {
		if ( i == d ) _ff_set(x1,x0);
		ff_square(x0,x0);
		ff_square(y,y);
	}
	if ( i == d ) _ff_set(x1,x0);
	k0 = ff_2Sylow_dlog (x0,y,d);
//printf ("k0=%d\n", k0);
	if ( k0 < 0 ) return -1;
	k1 = (1<<e)-k0;
	ff_exp_ui(&t, &a, k1);
	_ff_mult(y,b,t);
	k1 = ff_2Sylow_dlog(x1,y,e-d);
//printf("k1=%d\n", k1);
	if ( k1 < 0 ) return -1;
	return (1<<d)*k1+k0;
}

void _ff_setup_3g (void)
{
	register int i;
	register unsigned long n;
	ff_t r,s,t;

	if ( _ff_3g || ! _ff_p1mod3 ) return;
	_ff_cbrt_setup = 1;
	n = (_ff_p-1)/(3*_ff_p3_m);					// p-1 = 3*m*n

	// we could use cubic reciprocity here
	_ff_set_one(t);
	_ff_x2(t);
	for(;;) { 
		ff_exp_ui (&r,&t, _ff_p3_m);				// exponentiation into 3-Sylow group first
		if ( _ff_one(r) ) { _ff_inc(t); continue; }
		if ( _ff_p3_e == 1 ) break;
		ff_exp_ui (&s, &r, n);					// s = t^((p-1)/3)
		if ( ! _ff_one(s) ) break;					// if s is not 1 then r is a generator of the 3-Sylow (and a cubic nonresidue)
		_ff_inc(t);
	}
	_ff_set (_ff_3g, r);
	_ff_set (_ff_3Sylow_tab[0][0], r);
	_ff_square (_ff_3Sylow_tab[0][1],_ff_3Sylow_tab[0][0]);
	for ( i = 1 ; i < _ff_p3_e ; i++ ) {
		_ff_mult(_ff_3Sylow_tab[i][0],_ff_3Sylow_tab[i-1][0],_ff_3Sylow_tab[i-1][1]);		// tab[i][0] = tab[i-1][0]^3
		_ff_square(_ff_3Sylow_tab[i][1], _ff_3Sylow_tab[i][0]);
	}
	_ff_set(_ff_cbrt_unity,_ff_3Sylow_tab[_ff_p3_e-1][0]);
//printf("cbrt_unity = %d\n", _ff_get_ui(_ff_cbrt_unity));
}

/*
   Tonelli-Shanks cuberoot algorithm

   This algorithm is deterministic (given precomputed 3-Sylow generator) - caller should flip coin to randomize choice of sqrt.
*/
int ff_cbrt (ff_t o[1], ff_t a[1])
{
	ff_t b, x, y;

	ff_setup_3g();						// always do this first to make sure _ff_cbrt_unity is set for caller
	if ( _ff_zero(a[0]) ) { _ff_set_zero (o[0]);  return 1; }
	
	if ( ! _ff_p1mod3 ) {
		// every element of F_p is a cubic residue
		if ( _ff_p == 3 ) {_ff_set(o[0],a[0]);  return 1; }
		ff_exp_ui (o, a, (2*_ff_p-1)/3);			// o = a^((2p-1)/3) (note that (2p-1)/3 = 1/3 mod (p-1)
		return 1;
	}

	if ( _ff_p3_m1mod3 ) {
		ff_exp_ui (&x, a, (_ff_p3_m-1)/3);		// x = a^((m-1))/3)
		_ff_square (y, x); _ff_mult (b, x, y);		// b = x^3 = a^(m-1)
		ff_mult (b, b, a[0]);					// b = a^m is in the 3-Sylow subgroup
	} else {
		ff_exp_ui (&x, a, (_ff_p3_m-2)/3);		// x = a^(m-2)/3
		_ff_square (y, x); _ff_mult (b, x, y);		// b = x^3 = a^(m-2)
		_ff_square (y,a[0]); ff_mult (b, b, y);	// b = a^m is in the 3-Sylow subgroup
	}

	// Check if be is the 3-Sylow generator or its square--this happens with probability 2/3^e.  For 2/3 of the primes this is 2/3, which makes it worth avoiding the call to ff_2Sylow_invcbrt
	if ( _ff_equal(b,_ff_3Sylow_tab[0][0]) || _ff_equal(b,_ff_3Sylow_tab[0][1]) ) return 0;
	if ( ! ff_3Sylow_invcbrt(&y,&b) ) return 0;
	if ( _ff_p3_m1mod3 ) { ff_square(x,x); ff_square(y,y); }	// when m is 1mod3 we want to x=a^((2m-2)/3) and y=a^(-2m/3)
	ff_mult (x, x, a[0]);  ff_mult(x,x,y);					// either x = a^(2m-2)/3 *a* a^(-2m/3) = a^(1/3) or x=a^(m-2)/3*a*a^(-m/3) = a^(1/3)
#if ! FF_FAST
	_ff_square (y, x);  ff_mult(y,y,x);
	if ( ! _ff_equal(y,a[0]) ) { printf ("ff_cbrt failed, %lu^3 = %lu != %lu\n", _ff_get_ui(x), _ff_get_ui(y), _ff_get_ui(a[0])); exit (0); }
#endif
	_ff_set (o[0], x); 
	return 1;
}

// computes a^{-1/3} for a in the 3-Sylow subgroup, returns 0 if not a quadratic residue
// uses precomputed table of 2e powers of the 3-Sylow generator to reduce running time by a factor of 2 over standard Tonelli-Shanks (still O(e^2)).
int ff_3Sylow_invcbrt (ff_t o[1], ff_t a[1])
{
	register ff_t b, q, q1, t, w1, w2;
	register int i, j, k, s;
	
	ff_setup_3g();										// always do this first to make sure _ff_cbrt_unity is set for caller	
	// handle easy cases first
	if ( _ff_one(a[0]) ) { _ff_set_one(o[0]);  return 1; }		// use 1 as the cube root of 1
	if ( _ff_p3_e == 1 ) return 0;

	// set w1 and w2 to the two elements of order 3 in the 3-Sylow (i.e. the two non-trivial cube roots of unity)
	_ff_set (w1, _ff_3Sylow_tab[_ff_p3_e-1][0]);
	_ff_set (w2, _ff_3Sylow_tab[_ff_p3_e-1][1]);
	_ff_set_one (t);
	_ff_set (b,a[0]);
	do {
		_ff_set (q, b);
		for ( s = 1 ; s < _ff_p3_e+1 ; s++ ) {		// s<e+1 is just a safety check in case a isn't in the 3-Sylow, this could be removed
			j=1;
			if ( _ff_equal(q,w1) ) break;
			j=0;
			if ( _ff_equal(q,w2) ) break;
			_ff_set(q1,q);
			ff_square (q,q);  ff_mult(q,q,q1);
		}
		k = _ff_p3_e-s;
#if ! FF_FAST
		if ( k < 0 ) { printf ("Unexpected result: k<0 in ff_3Sylow_invsqrt?!  a = %lu, p = %lu\n", _ff_get_ui(a[0]), _ff_p);  exit (1); }
#endif
		if ( k <= 0 ) return 0;
		ff_mult (b, b, _ff_3Sylow_tab[k][j]);		// the product of all the elements S[k] we multiply into b here is b^{-1}, since we terminate with b=1
		ff_mult (t, t, _ff_3Sylow_tab[k-1][j]);		// multiply t by S[k-1]=cbrt(S[k]), the product of all these will be b^{-1/3}
	} while ( !_ff_one(b) );
#if ! FF_FAST
	_ff_square (b, t); ff_mult(b,b,t);
	ff_mult(b,b,a[0]);
	if ( ! _ff_one(b) ) { printf ("ff_3Sylow_invcbrt failed, %lu^3 *  %lu != 1\n", _ff_get_ui(t), _ff_get_ui(a[0])); exit (0); }
#endif
	_ff_set(o[0],t);
	return 1;
}


// standard 4-ary exponentiation (fixed 2-bit window)
void ff_exp_ui (ff_t o[1], ff_t a[1], unsigned long e)
{
	register int i, j;
	register ff_t c;
	register unsigned long m;
	ff_t b[4];
	
	if ( ! e ) { _ff_set_one (o[0]);  return; }
//	if ( _ff_one(e) ) { _ff_set (o[0], a[0]);  return; }
	// avoid tests to optimize for e <3 or a==0,1,-1
	i = _asm_highbit(e);
	if ( i&1 ) i--;
	m = 3UL<<i;
	_ff_set (b[1], a[0]);
	_ff_square (b[2],b[1]);
	_ff_mult(b[3],b[2],b[1]);
	_ff_set (c, b[(m&e)>>i]);
	for ( m>>=2,i-=2 ; m ; m>>=2, i-=2 ) {
		ff_square(c,c);  ff_square(c,c);
		j = (m&e)>>i;
		if ( j ) ff_mult(c,c,b[j]);
	}
//printf ("%ld^%ld=%ld mod %ld\n", _ff_get_ui(a[0]),e,_ff_get_ui(c),_ff_p);
	_ff_set (o[0], c);
}


// This function is copied from mpzutil.c to make ff.c/ff.h/asm.h self-contained
unsigned long ff_ui_inverse (unsigned long a, unsigned long m)
{
	register unsigned long q, r, r0, r1;
	register long t, t0, t1;

	if ( a >= m ) a = a%m;
	if ( a == 0 ) return 0;

	t0 = 0;  t1 = 1;
	r0 = m;  r1 = a;
	while ( r1 > 0 ) {
		q = r0 / r1;
		r = r0 - q*r1;
		r0 = r1;  r1 = r;
		t = t0 - q*t1;
		t0 = t1;  t1 = t;
	}
	if ( r0 > 1 ) return 0;
	if ( t0 < 0 ) return m - ((unsigned long)(-t0));
	return (unsigned long)t0;
}

/*
  Algorithm 1.4.10, p. 29 of [CANT] trimmed down to just handle a >= 0, b > 0 odd

  This algorithm is used to find a non-residue for use by ff_sqrt.

  NOTE: this is slower than both ff_fast_sqrt and ff_sqrt (after setup) due to the
  high cost of modular reductions (as opposed to Montgomery arithmetic).  It does have
  asymptotically superior "bit" complexity, O(lg^2(n)), but this doesn't mean much when
  the cost of multiplying two unsigned longs is effectively a constant.
*/
int ff_ui_legendre (unsigned long a, unsigned long b)
{
	register unsigned long r;
	register int k, v;
	
#if ! FF_FAST
	if ( !  (b&0x1) ) { printf ("ff_ui_Legendre, b = %lu must be odd!\n", b);  exit (0); }
#endif
	k = 1;
	while ( a ) {
		for ( v = 0 ; ! (a&0x1) ; v++ ) a >>= 1;
		if ( v&0x1 )  if ( (b&0x7) ==3 || (b&0x7) ==5 ) k = -k;
		if ( (a&0x3) == 3 && (b&0x3) == 3 ) k = -k;
		r = a;   a = b % r;  b = r;
	}
	return ( b == 1 ? k : 0 );
}

/*
  Fast square-root algorithm using a single exponentiation, requires p = 3 mod 4 or p = 5 mod 8.
  Relies on Thm 8.1 on p. 107 of Bressoud "Factorization and Primality Testing".
  Returns 1 if successful, 0 if a is not a quadratic residue.

  Not currently used, ff_sqrt is effectively just as fast.
	
  Currently only supported in single word implementations.  Overlap is OK.
*/

int ff_fast_sqrt (ff_t o[1], ff_t a[1])
{
	register ff_t  i;
	ff_t temp;
	register int k;
	
#if FF_WORDS > 1
	return 0;
#endif
	if ( _ff_zero(a[0]) ) { _ff_set_zero(o[0]);  return 1; }
	if ( (_ff_p&3)==3 ) {
		ff_exp_ui (&temp, a, (_ff_p+1)>>2);
	} else if ((_ff_p&7)==5 ) {
		ff_exp_ui (&temp, a, (_ff_p+3)>>3);
		_ff_square (i,temp);
		if ( _ff_equal (i, a[0]) ) {
			_ff_set(o[0], temp);
			return 1;
		}
		ff_setup_2exp();
		ff_mult(temp, temp, _ff_2exp);
	} else {
		return 0;
	}
	_ff_square(i,temp);
	if ( ! _ff_equal(i,a[0]) ) return 0;
	_ff_set(o[0],temp);
	return 1;
}
