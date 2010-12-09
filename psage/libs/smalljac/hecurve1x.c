#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include "gmp.h"
#include "mpzutil.h"
#include "hecurve.h"
#include "ffwrapper.h"
#include "ffpoly.h"
#include "cstd.h"

/*
    Copyright 2008 Andrew V. Sutherland

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
	This module implements operations for genus 1 curves, both low level
	point operations and functions for exponentiation computing the
	group order, group structure, etc...
	
	These functions are heavily optimized and responsible for a more than a 2x
	improvement in genus 1 performance from  smalljac version 3 over version 2.
	
	We begin with basic point arithmetic functions (mostly inlined), the higher
	level functions appear below.  	At the end of the module are a collection of functions
	that were implement and test but are not used because they were found to be
	sub-optimal in the current application.
*/


long hecurve_g1_JC_pp_order (hec_jc_t *a, long p, long q, ff_t f1);
int hecurve_g1_dlog (hec_jc_t *a, hec_jc_t *b, long o, ff_t f1);
int hecurve_g1_bsgs_search (long *ord, hec_jc_t *b, long low, long high, int m, int a, ff_t f1);
long hecurve_g1_fastorder (ff_t x, ff_t y, long k, ff_t f1);
int hecurve_g1_test_exponent (long e, ff_t f[4]);
void hecurve_g1_p_reduce (hec_jc_t *b1, long *q1, hec_jc_t *b2, long *q2, long p, ff_t f1);

unsigned long hecurve_expbits;
unsigned long hecurve_steps;
unsigned long hecurve_retries;


/*
	We use a reduced form of the Chudnovsky Jacobian representation (JC) which uses (x,y,z^2,z^3) to represent the affine point (x/z^2,y/z^3), but does not maintain z.
	Not computing z saves 1 field multiplication when compsing an affine point with JC point (usually called mixed addition).  Squaring (doubling) costs an extra multiplication
	but two fewer additions and the performance is comparable (in fact, when tested on an AMD Athlon 64, JC doubling was slighly faster, not clear why).

	In terms of the Mumford representation used elsewhere, the point point (x,y,z^2,z^3) is u[1]=z, u[0]=-x/z^2, v[0]=y/x^3. (notice the sign on x)
        When z is 1 this corresponds to u(t)=(t-x) and v(t)=y, so that x is a root of u and y=v(x).

	Negating y inverts an element, as in the affine rep, and 2-torsion elements have y=0.
	Any element with z=0 (or z2 for reduced Chudnovsky) is considered the identity, regardless of the value of x and y.

	As we use multiplicative notation for group operations elsewhere, we generally speak of the elliptic curve group operation multiplicatively,
        so we multiply, square, and exponentiate points, rather than adding, doubling, or multiplying by a scalar.

	In our operation counts we don't distinguish field multiplication and squaring.  For prime fields where p fits in a machine word, they are effectively the same.
        We do count additions, as these are not negligible--roughly speaking, 1M is about 2.5A.  Field inversions are horribly expensive relative to Montgomery multiplication,
        costing 40M or more (for p~2^30, say).  Whenever possible, we do n field inversions in parallel for a cost of I+3(n-1)M, effectively 3M per inverted element, for n large.
*/

// We trust the compiler to be smart about register allocation and don't hassle with macros (the speed difference appears to be negligible in testing)
// It would be nice to have C++ reference data types here...


static inline int hecurve_g1_JC_id (hec_jc_t *p) { return _ff_zero(p->z2); }
static inline int hecurve_g1_JC_2tor (hec_jc_t *p) { return _ff_zero(p->y); }
static inline void hecurve_g1_JC_invert (hec_jc_t *p) { ff_negate(p->y); }

static inline int hecurve_g1_JC_cmp (hec_jc_t *p1, hec_jc_t *p2)
{
	register ff_t t0,t1;
	
	_ff_mult(t0,p1->x,p2->z2); _ff_mult(t1,p2->x,p1->z2);
	if ( ! _ff_equal(t0,t1) ) return 0;
	_ff_mult(t0,p1->y,p2->z3); _ff_mult(t1,p2->y,p1->z3);
	return ( _ff_equal(t0,t1) ? 1 : -1 );
}

// converts reduced Chudnovsky Jacobian to affine coords, given the inverse of z3, overlap ok
static inline void hecurve_g1_JC_to_A (ff_t *px, ff_t *py, hec_jc_t *p, ff_t z3inv)
{
	register ff_t t0,t1;

	ff_mult(*py,p->y,z3inv);
	_ff_mult(t0,p->z2,z3inv);
	_ff_square(t1,t0);
	ff_mult(*px,p->x,t1);
	// 4M
}

// squares (doubles) an affine point into reduced Chudnovsky Jacobian coords
static inline void hecurve_g1_2AJC (hec_jc_t *p, ff_t x, ff_t y, ff_t f1)
{
	register ff_t t0,t1,t2,a,b;

	_ff_add(t0,y,y);  _ff_mult(t1,t0,y); _ff_add(p->z2,t1,t1); _ff_mult(a,x,p->z2);		// a=4xy^2, t1=2y^2, z2=4y^2
	_ff_mult(p->z3,t0,p->z2);												// z3=8y^3
	_ff_square(t1,x); _ff_add(b,t1,t1); _ff_addto(b,t1); _ff_addto(b,f1);				// b = 3x^2+f1*1^4
	_ff_square(t0,b); _ff_add(t2,a,a); _ff_sub(p->x,t0,t2);						// x = b^2-2a
	_ff_sub(t1,a,p->x); _ff_mult(t2,t1,b); _ff_mult(t0, y,p->z3); _ff_sub(p->y,t2,t0);	// y = b(a-x)-8y^4 (note we use the new x here and new z3=8y^3)
	// 7M+9A
}

// squares (doubles) an affine point into reduced Chudnovsky Jacobian coords and also sets az4 to f1*z3^4=16y^4 (cached value for squaring used below)
static inline void hecurve_g1_2AJC_cache (hec_jc_t *p, ff_t x, ff_t y, ff_t f1, ff_t *az4)
{
	register ff_t t0,t1,t2,a,b;

	_ff_add(t0,y,y);  _ff_mult(t1,t0,y); _ff_add(p->z2,t1,t1); _ff_mult(a,x,p->z2);		// a=4xy^2, t1=2y^2, z2=4y^2
	_ff_mult(p->z3,t0,p->z2);												// z3=8y^3
	_ff_square(t1,x); _ff_add(b,t1,t1); _ff_addto(b,t1); _ff_addto(b,f1);				// b = 3x^2+f1*1^4
	_ff_square(t0,b); _ff_add(t2,a,a); _ff_sub(p->x,t0,t2);						// x = b^2-2a
	_ff_sub(t1,a,p->x); _ff_mult(t2,t1,b); _ff_mult(t0, y,p->z3); _ff_sub(p->y,t2,t0);	// y = b(a-x)-8y^4 (note we use the new x here and new z3=8y^3)
	_ff_add(*az4,t0,t0); 
	// 7M+10A
}

// squares a point p1 in reduced Chudnovsky Jacobian coords (p3 is the output, may be equal to p1)
// This code requires one more multiplication than doubling in standard Jacobian coordinates (but two fewer additions, which makes it a close call)
// Surprisingly, it is actually slightly faster than doubling in Jacobian coords when tested on an AMD Athlon 64 (YMMV).
static inline void hecurve_g1_2JC (hec_jc_t *p3, hec_jc_t *p1, ff_t f1)
{
	register ff_t a, b, c, t0, t1, t2;
	
	_ff_square(t0,p1->x); _ff_add(t1,t0,t0); _ff_addto(t1,t0); 					// t1 = 3x^2
	_ff_square(t2,p1->z2);  _ff_mult(t0,f1,t2); _ff_add(b,t1,t0);					// b = 3x^2+f1*z^4	(note that f1=a4 in 13.2.1.c  of HECHECC p.282)
	_ff_add(c,p1->y,p1->y); _ff_square(t2,c); _ff_mult(a,t2,p1->x);				// a=4xy^2, c=2y, t2=4y^2
	_ff_add(t0,a,a); _ff_square(t1,b); _ff_sub(p3->x,t1,t0);						// x = b^2-2a
	_ff_mult(p3->z2,p1->z2,t2); _ff_mult(t1,c,t2); _ff_mult(p3->z3,p1->z3,t1);		// z2=4y^2z2, z3=8y^3z3, t1=8y^3
	_ff_mult(c,t1,p1->y);												// c = 8y^4
	_ff_sub(t0,a,p3->x); _ff_mult(t2,t0,b); _ff_sub(p3->y,t2,c);					// y = b(a-x)-c   -- note we use the new x value here
	// 11M+8A
}

// same as above except the parameter az4 contains a_4*z1^4 = f1*z1^4 and is updated to hold a_4*z3^4 (cost 1M less, 1A more)
static inline void hecurve_g1_2JC_cache (hec_jc_t *p3, hec_jc_t *p1, ff_t *az4)
{
	register ff_t a, b, c, t0, t1, t2;
	
	_ff_square(t0,p1->x); _ff_add(t1,t0,t0); _ff_addto(t1,t0);  _ff_add(b,t1,*az4);	// b = 3x^2+f1*z^4
	_ff_add(c,p1->y,p1->y); _ff_square(t2,c); _ff_mult(a,t2,p1->x);				// a=4xy^2, c=2y, t2=4y^2
	_ff_add(t0,a,a); _ff_square(t1,b); _ff_sub(p3->x,t1,t0);						// x = b^2-2a
	_ff_mult(p3->z2,p1->z2,t2); _ff_mult(t1,c,t2); _ff_mult(p3->z3,p1->z3,t1);		// z2=4y^2z2, z3=8y^3z3, t1=8y^3
	_ff_mult(c,t1,p1->y); _ff_add(t0,c,c); ff_mult(*az4,*az4,t0);					// c = 8y^4, az4 = az4*16y^4
	_ff_sub(t0,a,p3->x); _ff_mult(t2,t0,b); _ff_sub(p3->y,t2,c);					// y = b(a-x)-c   -- note we use the new x value here
	// 10M+9A
}

// multiplies a non-identity point p1 in reduced Chudnovsky Jacobians coords by a (non-identity) affine point and puts result in p
// using the reduced Chudnosky form (no z, only z^2 and z^3) saves one field mult (this is faster than any of the alternatives in table 13.3 of HECHECC on p.284)
static inline void hecurve_g1_AJC (hec_jc_t *p, hec_jc_t *p1, ff_t x0, ff_t y0, ff_t f1)
{
	register ff_t a, c, e, f, t0, t1, t2;
	
	if ( _ff_zero(p1->z2) ) { _ff_set(p->x,x0); _ff_set(p->y,y0); _ff_set_one(p->z2); _ff_set_one(p->z3); return; }
	_ff_mult(a,x0,p1->z2);											// a = x0z^2, b = x*1^2=x (z0=1)
	_ff_sub(e,p1->x,a);												// e = a-b
	_ff_mult(c,y0,p1->z3);											// c = y0z^3, d=y*1^3=y (z0=1)
	_ff_sub(f,p1->y,c);												// f = d-c
	if ( _ff_zero(e) && _ff_zero(f) ) { hecurve_g1_2AJC (p, x0, y0, f1); return; }	// must use doubling code here, but at least it's an affine double (7M+9A) (inverses are ok)
	_ff_square(t0,e); _ff_mult(t1,t0,e); _ff_mult(t2,t0,a);					// e2=e^2, t1=e^3, t2=ae^2
	ff_mult(p->z2,p1->z2,t0); ff_mult(p->z3,p1->z3,t1);					// z2 = 1*z2*e^2, z3=1*z3*e^3
	_ff_square(t0,f); _ff_sub(p->x,t0,t1); _ff_add(t0,t2,t2); _ff_subfrom(p->x,t0);	// x = f^2-e^3-2ae^2
	_ff_sub(t0,t2,p->x); _ff_mult(t2,f,t0); _ff_mult(t0,t1,c); _ff_sub(p->y,t2,t0);	// y = f(ae^2-x) - ce^3
	// 10M+6A
}

// multiplies p1 in reduced Chudnovsky Jacobian coords by p2 in reduced Chudnovsky Jacobian coords  and puts the result in p1.
// Handles identity and will square (double) if needed.  Cost is 13M+7A, one less than with standard Chudnovsky Jacobian coords.
static inline void hecurve_g1_JCJC (hec_jc_t *p1, hec_jc_t *p2, ff_t f1)
{
	register ff_t a, b, c, d, e, f, t0, t1, t2;

	// we need to handle identity cases in general (although in some cases we might be able to rule this out)
	if ( _ff_zero(p2->z2) ) return;
	if ( _ff_zero(p1->z2) ) { *p1=*p2; return; }
	_ff_mult(a,p1->x,p2->z2);  _ff_mult(b,p2->x,p1->z2);  _ff_sub(e,b,a);			// a = x1z2^2, b = x2z1^2, e = a-b
	_ff_mult(c,p1->y,p2->z3);  _ff_mult(d,p2->y,p1->z3);  _ff_sub(f,d,c);			// c = y1z2^2, d = y2z1^3, f = d-c
	if ( _ff_zero(e) && _ff_zero(f) ) { hecurve_g1_2JC (p1,p1,f1); return; }			// must double if pts are equal (11M+8A), inverses will end up with z2=0, so we let them through
	_ff_square(t1,e); ff_mult(t0,p1->z2,p2->z2);  _ff_mult(p1->z2,t0,t1);			// z^2 = z1^2z2^2e^2, t1=e^2
	_ff_mult(t2,t1,e); ff_mult(t0,p1->z3,p2->z3);  _ff_mult(p1->z3,t0,t2);			// z^2 = z1^2z2^2e^2, t2=e^3
	ff_mult(t1,t1,a);													// t1=ae^2
	_ff_square(t0,f); _ff_sub(p1->x,t0,t2); _ff_add(t0,t1,t1); _ff_subfrom(p1->x,t0);	// x = f^2-e^3-2ae^2
	_ff_sub(t0,t1,p1->x); _ff_mult(t1,f,t0); _ff_mult(t0,t2,c); _ff_sub(p1->y,t1,t0);	// y = f(ae^2-x) - ce^3 -- note we use the new x here
	// 13M+7A
}

/*
	Computes p=(x,y)^n where (x,y) !=1 is in affine coordinates and p is in Jacobian coordinates
	It takes advantage of the fact that all additions are of the form J+A, requiring only 11M, doubling is done in J+J form, using 10M
*/
void hecurve_g1_AJC_exp_ui (hec_jc_t *p, ff_t x0, ff_t y0, unsigned long n, ff_t f1)
{
	register unsigned long m;
	unsigned long pbits, nbits;
	register ff_t negy0;
	int i;

	if ( n == 0 ) { _ff_set_zero(p->z2); return; }
	if ( n == 1 ) { _ff_set(p->x,x0); _ff_set(p->y,y0); _ff_set_one(p->z2); _ff_set_one(p->z3);  return; }
	hecurve_g1_2AJC(p,x0,y0,f1);
	if ( n == 2 ) return;
	ui_NAF(&pbits,&nbits,n);
	i = _asm_highbit(pbits)-2;						// we know the top two bits of the NAF are 10
	_ff_neg(negy0,y0);
	m = (1UL<<i);
	for ( ; m ; m >>= 1 ) {
		hecurve_g1_2JC(p,p,f1);	 					// 11M+8A
		if ( m&pbits ) hecurve_g1_AJC(p,p,x0,y0,f1);		// 10M+6A
		if ( m&nbits ) hecurve_g1_AJC(p,p,x0,negy0,f1);	// 10M+6A
	}
}

/*
	Computes p=a^e where a!=1 is in reduced Jacobian Chudnovsky coords, and so is p. overlap is ok.
*/
void hecurve_g1_JC_exp_ui (hec_jc_t *p, hec_jc_t *a, unsigned long n, ff_t f1)
{
	register unsigned long m;
	unsigned long pbits, nbits;
	hec_jc_t t, ai;
	int i;
	
	if ( n == 0 ) { _ff_set_zero(p->z2); return; }
	if ( n == 1 ) { *p=*a; return; }
	hecurve_g1_2JC(&t,a,f1);							// 11M+8A

	if ( n == 2 ) { *p = t; return; }
	ui_NAF(&pbits,&nbits,n);
	i = _asm_highbit(pbits)-2;
	ai=*a;  hecurve_g1_JC_invert(&ai);
	m = (1UL<<i);
	for ( ; m ; m >>= 1 ) {
		hecurve_g1_2JC(&t,&t,f1);					// 11M+8A
		if ( m&pbits ) hecurve_g1_JCJC(&t,a,f1);			// 13M+7A
		if ( m&nbits ) hecurve_g1_JCJC(&t,&ai,f1);		// 13M+7A
	}
	*p = t;
}

// Computes b=a^e where a and b are both in affine coordinates
void hecurve_g1_exp_ui (ff_t u[HECURVE_GENUS+1], ff_t v[HECURVE_GENUS], ff_t u1[HECURVE_GENUS+1], ff_t v1[HECURVE_GENUS], unsigned long n, ff_t f[HECURVE_DEGREE+1])
{
	ff_t zinv,x0;
	hec_jc_t o[1];
	
	if ( _hecurve_is_identity(u1,v1) ) { _hecurve_set_identity(u,v);  return; }
	if ( ! _ff_one(u1[1]) ) { printf ("p=%ld, input to hecurve_g1_exp_ui most be in affine coords!\n",_ff_p); hecurve_print(u,v); exit(0); }
	_ff_neg(x0,u1[0]);
	hecurve_g1_AJC_exp_ui (o, x0, v1[0], n, f[1]);
	if ( hecurve_g1_JC_id(o) ) { _hecurve_set_identity(u,v);  return; }
	ff_invert(zinv,o->z3);
	hecurve_g1_JC_to_A(&x0,v,o,zinv);
	_ff_neg(u[0],x0);
	_ff_set_one(u[1]);
}

// Combines precomputed values p[i]=p[0]^(2^i) to compute p[0]^n.  Assumes all required powers are present: NAF potentially requires one more bit!
// CAUTION: if this is false, the resulting behavior is very strange and unpredictable.  You have been warned.
void hecurve_g1_JC_exp_powers(hec_jc_t *o, hec_jc_t *p, unsigned long n, ff_t f1)
{
	register unsigned long m;
	unsigned long pbits, nbits;
	hec_jc_t t;
	int i;

	// handle small n quickly
	switch (n) {
	case 0:	_ff_set_zero(o->z2); return; 
	case 1:	*o=p[0]; return;
	case 2:	*o=p[1]; return;
	case 3:	*o=p[0]; hecurve_g1_JCJC(o,p+1,f1); return;
	case 4:	*o=p[2]; return;
	case 5:	*o=p[2]; hecurve_g1_JCJC(o,p,f1); return;
	case 6:	*o=p[2]; hecurve_g1_JCJC(o,p+1,f1); return;
	case 7:	*o=p[0]; hecurve_g1_JC_invert(o); hecurve_g1_JCJC(o,p+3,f1); return;
	case 8:	*o=p[3]; return;
	}
	ui_NAF(&pbits,&nbits,n);
	i = _asm_highbit(pbits);
	*o = p[i];
	i-=2;
	m = (1UL<<i);
	for ( ; m ; m >>= 1, i-- ) {
		if ( (m&pbits) ) hecurve_g1_JCJC(o,p+i,f1);								// 13M+8A
		if ( (m&nbits) ) {													// 13M+10A (slightly more expensive)
			hecurve_g1_JC_invert(o); hecurve_g1_JCJC(o,p+i,f1); hecurve_g1_JC_invert(o);
		}
	}
}

// Sets p[i] = p[0]^(2^i) for i in [0,k]. Input is an affine pt not equal to the identity, output is in reduced Jacobian Chudnovsky coords.
void hecurve_g1_AJC_powers (hec_jc_t p[], ff_t x0, ff_t y0, int k, ff_t f1)
{
	register hec_jc_t *e;
	ff_t az4;
	
	_ff_set(p->x,x0);  _ff_set(p->y,y0);  _ff_set_one(p->z2); _ff_set_one(p->z3);
	hecurve_g1_2AJC_cache (p+1, x0, y0, f1, &az4);						// 7M + 10A
	for ( e=p+k,p+=2 ; p <= e ; p++ ) hecurve_g1_2JC_cache(p,p-1,&az4);	// 10M + 9A
}

// Sets p[i] = p[0]^(2^i) for i in [1,k]. Input is an pt not equal to the identity, input and output is in reduced Jacobian Chudnovsky coords.
void hecurve_g1_JC_powers (hec_jc_t p[], int k, ff_t f1)
{
	register hec_jc_t *e;

	// don't bother caching
	for ( e=p+k,p++ ; p <= e ; p++ ) hecurve_g1_2JC(p,p-1,f1);		 	// 11M + 8A
}

// Sets s[i] = s[0]*(x0,y0)^i for i < n.  Returns k<n if s[k] is the identity (in which case no further elements are computed), otherwise returns n.
int hecurve_g1_AJC_steps_1 (hec_jc_t s[], ff_t x0, ff_t y0, int n, ff_t f1)
{
	register hec_jc_t *end;
	
	if ( hecurve_g1_JC_id(s) ) return 0;
	end = s+n;
	for ( s++ ; s < end ; s++ ) {
		hecurve_g1_AJC (s,s-1,x0,y0,f1);								// 10M + 6A
		if ( hecurve_g1_JC_id(s) ) break;
	}
	hecurve_steps += n-(end-s);
	return n-(end-s);
}

// Sets s[-i] = s[0]*(x0,-y0)^i for i < n.  Returns k<n if s[k] is the identity (in which case no further elements are computed), otherwise returns n.
// Used for downward stepping - NOTE THAT s SHOULD POINT TO THE LAST ENTRY IN THE ARRAY
int hecurve_g1_AJC_dsteps_1 (hec_jc_t *s, ff_t x0, ff_t y0, int n, ff_t f1)
{
	register hec_jc_t *end;
	register ff_t negy0;
	
	if ( hecurve_g1_JC_id(s) ) return 0;
	end = s-n;
	_ff_neg(negy0,y0);
	for ( s-- ; s > end ; s-- ) {
		hecurve_g1_AJC (s,s+1,x0,negy0,f1);							// 10M + 6A
		if ( hecurve_g1_JC_id(s) ) break;
	}
	hecurve_steps += n-(s-end);
	return n-(s-end);
}

// sets s[2*i+j] = s[0]*(x0,y0)^(i+j)*(x1,y1)^i for 2*i+j < n with j=0 or 1.  Returns k<n if s[k] is the identity (in which case no further elements are computed), otherwise returns n.
int hecurve_g1_AJC_steps_2 (hec_jc_t s[], ff_t x0, ff_t y0, ff_t x1, ff_t y1, int n, ff_t f1)
{
	register hec_jc_t *end;
	register int i;
	
	if ( hecurve_g1_JC_id(s) ) return 0;
	end = s+n;
	for ( i=1,s++ ; s < end ; i++,s++ ) {
		if ( i&1) {
			hecurve_g1_AJC (s,s-1,x0,y0,f1);		// 10M + 6A
		} else {
			hecurve_g1_AJC (s,s-1,x1,y1,f1);		// 10M + 6A
		}
		if ( hecurve_g1_JC_id(s) ) break;
	}
	hecurve_steps += n-(end-s);
	return n-(end-s);
}

static inline int inv_mod3 (int e) { return ((e%3)==1?1:2); }

/*
	Fast group order computation for elliptic curves over F_p, designed for p<2^40, optimized for p~2^30.
	Returns 1 if successful and sets *o to the group order and if pd is non-null, sets *pd to gcd(m,6), where the group structure is Z/mZ x Z/nZ with m dividing n (possibly m=1).
	This information can be used to speed up group structure computations.

	Failure can only happen for p < 229, in which case 0 is returned and *o is a divisor of the group order highly likely equal to the group exponent (*pd is as also set, as above).
*/
long hecurve_g1_order (long *pd, ff_t f[4])
{
	long r, d, e, E, min, max, low, high, exp, M;
	hec_jc_t t[1];
	ff_t g[4],x,y,*h;
	int a,i,j,m,twist,s2known,tor3,tor4,flag8;

	/*
		We compute the 3-torsion subgroup first, because this information may lead us to work in the twist.
		The function ff_poly_g1_3tor() computes the 3-torsion subgroup of y^2=f(x), but if it is trivial, it will try to
		determine the 3-torsion subgroup of the twist, using the factorization pattern of the 3-division polynomial (a quartic).
	*/
	h = f;  twist = 0;
	tor3 = ff_poly_g1_3tor(f);	
	if ( tor3<0 ) { ff_poly_twist(g,f,3); h = g; twist = 1; tor3=-tor3; }		// negative return value indicates 3-torsion in the twist (but not in the primary)
	d = ( tor3==9 ? 3 : 1);										// d is a known divisor of |G|/lambda(|G|), i.e. it divides the group exponent (lambda(G)) and its square divides |G|.

	if ( _ff_p < HECURVE_G1_4TOR_MINP ) {							// when p is small, only compute 2-torsion but not 4-torsion (the marginal benefit does not justify the cost)
		 i = ff_poly_roots_d3(0,h);
		if ( i == 0 ) {											// no roots means there is no 2-torsion and the group order must be odd
			s2known = 1;  tor4 = 1;								// the flag s2known indicates we know the entire 2-Sylow subgroup
		} else {												// this means that once we factor it out, the remaining exponent is known to be odd
			s2known = 0;  tor4 = 2;								// if there is one root, we have 2-torsion Z/2Z, otherwise there are 3 roots and 2-torsion Z/2Z x Z/2Z
			if ( i > 1 ) { tor4 = 4;  d *= 2; }			
		}		
	} else {
		s2known = 0;
		flag8 = ( _ff_p < HECURVE_G1_8TOR_MINP ? 0 : 1);				// flag8 set if we also want to check whether the 2-Sylow subgroup is a cyclic group containing Z/8Z
		d *= ff_poly_g1_4tor(&tor4,h,flag8);						// compute the 4-torsion subgroup (and, optionally, check Z/8Z as well)
		if ( flag8 ) {
			if ( tor4 <= 4 ) s2known = 1;							// In these cases we know the entire 2-Sylow subgoup
		} else {
			if ( tor4 < 4 || (tor4==4 && d==2) ) s2known = 1;			// if the 4-torsion subgroup contains no elements order 4, it must be equal to the 2-Sylow subgroup.
		}
	}
	e = tor4*tor3/d;											// e is our known divisor of lambda(|G|)
	E = d*e;													// E is our known divisor of |G|
	m = (s2known==1?2:1)*(tor3==1?3:1);							// we know that gcd(|G|/tor4,m) = 1, which we can use to speed up the BSGS search
	
	a = ( tor3==1 && _ff_p1mod3 ? -1 : 0 );							// if neither the curve or its twist has 3-torsion and p=1mod3, we must have a_p=0 and group order 2 mod 3
															// For any divisor x of |G| this means |G|/x must be equal to -1/x mod 3 (and also mod 6 if m=6)
	r = (long)(2.0*sqrt(_ff_p));
	exp = _ff_p+1;
	min = exp-r;
	max = exp+r;
	low = _ui_ceil_ratio(min,E);
	high = max/E;
	if ( low > high ) { printf ("Error, no multiple of torsion derived e=%ld exists in interval [%ld,%ld] for p=%ld\n", e, min, max, _ff_p); exit(0); }
	if ( pd ) { *pd = ( twist && tor3==9 ? d/3 : d );  if ( !((*pd)&3) ) *pd/=2; }	// set *pd to reflect 2-torsion and 3-torsion information in y^2=f(x) (but not the twist)

//printf("%ld: tor3=%d, tor4=%d, e=%d, d=%d ", _ff_p, tor3, tor4, e, d); ff_poly_print(h,3);
	
	// Note that for a non-CM curve the probability that hecurve_g1_bsgs succeeds on the first try is close to 1 and grows with p (say 99.9% for p~2^20)
	for ( i = 0 ; i < HECURVE_G1_ORDER_RETRIES ; i++ ) {
		if ( low==high ) {exp=low; break; }
		if ( ! hecurve_random_point(&x,&y,h) ) continue;				// this can fail, e.g.  y^2=x^3+2x+2 over F_3 has no finite rational points
//printf ("%ld: e=%d, Random pt (%ld,%ld) on curve ", _ff_p, e, _ff_get_ui(x), _ff_get_ui(y));  ff_poly_print(h,3);
		hecurve_g1_AJC_exp_ui (t, x, y, e, h[1]);  if ( hecurve_g1_JC_id(t) ) continue;
		// handle order 2 elements here so that bsgs search can assume |t| > 2
		if ( hecurve_g1_JC_2tor(t) ) { 
			exp = 2;
 		//On paper, the code below should speed things up (it uses no inversions and fewer operations for small intervals), but in testing it actually slows things down slightly
		//} else if ( high-low < HECURVE_G1_SHORT_INTERVAL ) {
		//	if ( hecurve_g1_bsgs_short (&exp, &t, low, high, h[1]) ) break;
		} else {
			if ( a ) { a = -inv_mod3(E); if ( m==6 && a==-2 ) a = -5; }
			if ( hecurve_g1_bsgs_search (&exp, t, low, high, m, a, h[1]) ) break;
		}
		e *= exp;  E *= exp;
		low = _ui_ceil_ratio(min,E);
		high = max/E;
		// there are lot's of optimizations we could insert here, but none of them apply often enough to be worth doing, so keep it simple.
		hecurve_retries++;
	}
	if ( i < HECURVE_G1_ORDER_RETRIES ) {
		E *= exp;
		if ( twist ) E = 2*(_ff_p+1)-E;
		if ( E < min || E > max ) { printf ("Error, computed order %ld is not in Hasse-Weil interval [%ld,%ld] for p=%ld, h= ", E, min, max, _ff_p); ff_poly_print(h,3); exit(0); }
		return E;
	}
	// For non-CM curves we will almost never reach this point.
	
	if ( e > 2*r ) { printf("Error, exponent %ld for p=%ld has unique multiple in [%ld,%ld] that was not detected.\n", e, _ff_p, min, max);  exit(0); }
	if ( e < sqrt(min) ) { printf ("Error, exponent %ld is impossibly small for p=%ld\n", e, _ff_p); exit(0); }
	
	/*
		With very high probability (exponential in HECURVE_G1_ORDER_RETRIES), e is now the group exponent.   For all but 21 primes (the largest of which is 547)
		this uniquely determines the group order (this result is due to Cremona and Harley, or see Washington, "Elliptic Curves", Prop. 4.19).
		In fact, taking our knowledge of 2-torsion and 3-torsion into account, the group order is uniquely determined in every case (AVS: to be written up).
		Unfortunately, we can't prove that e is definitely the group exponent without computing the group structure.
	
		Instead, we apply results of Schoof and Mestre (see Washington, Prop. 4.18) which tells us that either the curve or its twist contains an element whose
		order has a unique multiple in the Hasse interval provided p > 229.  If one examines the exceptional cases for p<=229, one finds that if we also
		consider the information provided by 2-torsion, in every case but two (y^2=x^3+2 over F_13 and F_19) we obtain a known divisor of the group order
		with a unique multiple in the Hasse interval, and the two exceptional cases are addressed if we consider 3-torsion. (AVS: to be written up).

		In the two 3-torsion cases, we will have chosen h to be the twist with 3-torsion and should have succeeded above.  If not we mustn't have computed
		the full group exponent and we should try again (by recursively calling ourselves below)
	
		For the 2-torsion cases, the twist will also have 2-torsion, we just need to use this information when we try to find a unique multiple of the exponent
		in the twist below.  If we fail, it must mean we got unlucky and should try again.
	
		In summary, this yields a Las Vegas algorithm with provably correct output, for all odd p (for p=2 the curve y^2=f(x) is singular).
	*/
	
	m = (d%2?1:2);		// m is a known divisor of the twist group order based on the rank of the 2-torsion subgroup, note that we don't keep track of d anymore and reuse it below

	for ( M = e*_ui_ceil_ratio(min,e) ; M <= max ; M += e ) { d = M/e;  if ( !(e%d) && !((_ff_p-1)%d) ) break; }
	if ( M > max ) { printf ("Error: no multiple M of ambiguous e=%ld satisfies M/e | gcd(e,p-1) (p=%ld).\n", e, _ff_p); exit (0); }
	if ( ! twist ) { ff_poly_twist(g,f,3); h = g; } else { h = f; }	

	do {
		r = 2*(_ff_p+1)-M;
		// Attempt to prove that r is the order of the twist.  We don't bother trying to be efficient here, we expect to succeed on the first try.
		for ( i = 0 ; i < HECURVE_G1_ORDER_RETRIES ; i++ ) {
			if ( ! hecurve_random_point(&x,&y,h) ) { printf("hecurve_random_point failed!"); exit (0); }
			hecurve_g1_AJC_exp_ui (t, x, y, r, h[1]);
			if ( ! hecurve_g1_JC_id(t) ) break;											// can't be the right M
			exp = m*hecurve_g1_fastorder(x,y,r,h[1]);									// compute the order of our random point, times our 2-torsion derived info
			low = exp*_ui_ceil_ratio(min,exp);											// low is the first multiple of exp >= min
			if ( low <= max && low+exp > max ) { return ( twist ? 2*(_ff_p+1)-M : M ); }		// success
		}
		// try another candidate, if there is one (this can happen for small p)
		for ( M+= e; M <= max ; M += e ) { d = M/e;  if ( !(e%d) && !((_ff_p-1)%d) ) break; }			
	} while ( M <= max );
	// Otherwise, try again -- we must have not gotten the full group exponent (this essentially never happens unless HECURVE_G1_ORDER_RETRIES is set quite low).
	return hecurve_g1_order(pd,f);
}

/*
	Given the group order N, compute the isomorphism type of the group structure of an elliptic curve, of the form Z/n1Z x Z/n2Z with n1 dividing n2 (and also p-1), possibly n1=1.
	The parameter d specifies gcd(6,n1), which can be derived from 2-torsion and 3-torsion info by hecurve_g1_order above (specify 0 if not known).
	The rank of the group (1 or 2) is returned, and n[0]=n2 if n1=1, o.w. n[0]=n1 and n[1]=n2.  Note that n2 is the group exponent.

	It is possible to compute the group invariants (n1 and n2) in probabilistic (Las Vegas) polynomial time using the Weil pairing via Miller's algorithm
	(although this doesn't actually determine a basis, i.e. independent elements of order n1 and n2).

	We take a simpler generic approach which is very fast in the typical case.  We may in the future want to invoke (refinements of) Miller's algorithm for
	hard cases (large prime dividing p-1 whose square divides N).  As it stands, the algorithm below has complexity O(q^(1/2)) where q is the largest prime dividing p-1
	whose square divides N.  This is O(N^(1/4)) in the worst case, but this rarely happens and in practice computing the group structure takes about 10% the
	time it takes to compute the group order (for a non-CM curve).  (It might be interesting to do a performance comparison wth Miller's algorithm at some point)
	
	We just return the group invariants here, and optimize for computing these as quickly as possible, but the algorithm below can easily be modified to return a basis
	without changing the complexity (but the constant factors will be slightly worse).  This is a Las Vegas algorithm, i.e. provably correct output and bounded expected running time.
*/
int hecurve_g1_group_structure (long n[2], long N, long d, ff_t f[4])
{
	long m,n1,n2,e;
	hec_jc_t b1, b2;
	unsigned long p[MAX_UI_PP_FACTORS], q[MAX_UI_PP_FACTORS], h[MAX_UI_PP_FACTORS];
	ff_t x,y;
	long q0, q1, q2, r;
	int i, j, k;
	
	/*
		For each prime q|N, if q does not divide p-1 we know the q-Sylow subgroup is cyclic.  Even when q does divide p-1, if q^2 does not divide N, we
		again know the q-Sylow is cyclic (of prime order in this case).  Applying just these 2 facts will often suffice to determine the group structure without
		using any group operations, particularly if information in d is also applied (e.g. 2-torsion and 3-torsion info).
	*/
	k = ui_factor(p,h,ui_gcd(_ff_p-1,N));
	for ( i = j = 0 ; i < k ; i++ ) {												// make a list of primes p[i] dividing p-1 whose square divides N
		if ( p[i]==2 && d && d%2 ) continue;									// d=gcd(6,n1) not divisible by 2 implies 2-Sylow must be cyclic
		if ( p[i]==3 && d && d%3 ) continue;									// d=gcd(6,n1) not divisible by 3 implies 3-Sylow must be cyclic
		if ( !(N%(p[i]*p[i])) ) p[j++] = p[i];
	}
	n1 = 1;  n2 = N;
	for ( i = 0 ; i < j ; i++ ) {
		for ( q[i] = p[i]*p[i] ; !(N%(q[i]*p[i])) ; q[i] *= p[i] );						// determine power of p[i] dividing N
		n2 /= q[i];														// remove q[i] from n2
	}
//printf("%ld: N=%d ", _ff_p, N); ff_poly_print(f,3);
	// Now determine the p[i]-Sylow subgroups for p[i] dividing p-1 and p[i]^2 dividing N (if there are none, we have no work to do).
	for ( i = 0 ; i < j ; i++ ) {
		if ( d ) for ( q0 = 1 ; ! (d%(p[i]*q0)) ; q0 *= p[i] );	else q0 = 1;				// if d was specified, use d to compute q0, a known divisor of q1 (and q2), possibly 1
		e = N/(q[i]*n1);													// the p[i]-Sylow subgroup must lie in the image of the e-power map (we can take out the n1 we know)
		for ( q1 = q2 = q0, k=0 ; q1*q2 < q[i] && k < 1000 ; k++ ) {				// bound the number of retries as a safety valve, just in case the group order was invalid
			if ( ! hecurve_random_point(&x,&y,f) ) continue;						// get a random finite point
			hecurve_g1_AJC_exp_ui (&b1, x, y, e, f[1]);							// get an element of the p[i]-Sylow via the e-power map
			r = hecurve_g1_JC_pp_order (&b1,p[i],q[i],f[1]);						// compute its order, which will be a power of p[i]
			if ( r <= q1 ) continue;											// if the order of b1 is not strictly larger than q1, it isn't going to help us
			if ( q2==q0 ) { b2=b1; q2 = r; continue; }							// first time through we will just set q2
			hecurve_g1_p_reduce (&b1, &r, &b2, &q2, p[i], f[1]);					// reduce to a basis for <b1,b2>
			q1 = ( r > q0 ? r : q0);											// note that if r<q0, then q1 is not the order of b1, but we don't care
		}
		if ( k == 1000 ) { printf ("Group structure computation failed on %d-Sylow of size %d with group of order %d on curve over F_%d: ", p[i], q[i], N, _ff_p); ff_poly_print(f,3); exit (0); }
		n1 *= q1;  n2 *= q2;
	}
	if ( n1*n2 != N ) { printf ("bug in hecurve_g1_group_order, %d*%d != %d,  F_%ld curve ", n1, n2, N, _ff_p);  ff_poly_print(f,3); exit(0); }		// sanity check
	if ( n1 == 1 ) { n[0] = n2; return 1; }
	n[0] = n1;  n[1] = n2;
	return 2;
}

/*
	Slow Monte Carlo algorithm to test the group exponent, use for debugging/testing only.
*/
int hecurve_g1_test_exponent (long e, ff_t f[4])
{
	hec_jc_t t[1];
	ff_t x,y;
	long m, n;
	int i;
	
	n = 1;
	for ( i = 0 ; i < 100 ; i++ ) {
		if ( ! hecurve_random_point(&x,&y,f) ) continue;
		hecurve_g1_AJC_exp_ui(t,x,y,e,f[1]);
		if ( ! hecurve_g1_JC_id(t) ) return 0;
		m = hecurve_g1_fastorder (x,y,e,f[1]);
		if ( n%m ) n = m*(n/ui_gcd(m,n));
	}
	return (n==e);
}

/*
	Given non-trivial elements b1 and b2 of the p-Sylow subgroup with |b1|=q1 and |b2|=q2  powers of p (p here is a divisor of the group order, not the characteristic of the field),
	The function below computes a basis for <b1,b2> and updates b1,b2,q1, and q2 appropriately, with q1 < q2 (possibly q1=1 if <b1,b2> is cyclic).

	The complexity is O(log_p(q1)*(sqrt(p)+log(q2)).
	AVS: the code below is a special case of a new algorithm for group structure computation (a simpler and faster version of Algorithm 9.2 in my thesis which also incorporates
	a few ideas from Teske's Pohlig-Hellman paper).  I am in the process of writing a paper which I will eventually reference here.
*/
void hecurve_g1_p_reduce (hec_jc_t *b1, long *q1, hec_jc_t *b2, long *q2, long p, ff_t f1)
{
	hec_jc_t t[1],a1[1],a2[1];
	long k;
	int i;
	
	if ( *q1 > *q2 ) { t[0] = *b1;  *b1 = *b2; *b2 = t[0];  k = *q1; *q1 = *q2; *q2 = k; }	// swap inputs if needed so that q1 <= q2

	hecurve_g1_JC_exp_ui(a2,b2,*q2/p,f1);										// a2=b2^(q2/p), <a2> is the subgroup of <b2> of order p
	while (*q1>1) {
		hecurve_g1_JC_exp_ui(a1,b1,(*q1)/p,f1);									// a1=b1^(q1/p), <a1> is the subgroup of <b1> of order p
		k = hecurve_g1_dlog (a2,a1,p,f1);										// compute least k>0 s.t. a1=a2^k, if possible
		if ( ! k ) break;													// if we can't, then a1 and a2 are independent, hence so are b1 and b2, and we are done
		hecurve_g1_JC_exp_ui(t,b2,(*q2)/(*q1)*k,f1);  hecurve_g1_JC_invert(t);			// compute t=b2^-(q2/q1*k)
		hecurve_g1_JCJC (b1,t,f1);											// replace b1 with b1*b2^-(q2/q1*k), this reduces the order of b1 to at most q1/p since
																		// (b1*b2^-(q2/q1*k))^(q1/p) = a1*a2^{-k} = 1, and we do not change <b1,b2> by doing this.
		k = hecurve_g1_JC_pp_order(b1,p,*q1,f1);								// recompute q1 (it could be less than q1/p)
		if ( k >= *q1 ) { printf ("%ld: Error in hecurve_g1_p_reduce(%d,%d) element order %d not reduced\n", _ff_p, *q1, *q2, k); exit (0); }		// sanity check
		*q1 = k;
	}
}


long hecurve_g1_JC_pp_order (hec_jc_t *a, long p, long q, ff_t f1)
{
	register long r;
	hec_jc_t t[1];
	
	if ( hecurve_g1_JC_id(a) ) return 1;
	if ( p==q ) return p;
	hecurve_g1_JC_exp_ui (t,a,p,f1);
	if ( hecurve_g1_JC_id(t) ) return p;
	for ( r = p*p ; r < q ; r *= p ) {
		hecurve_g1_JC_exp_ui (t,t,p,f1);
		if ( hecurve_g1_JC_id(t) ) return r;
	}
	return q;
}

/*
	Fast and simple small hash table lookup used by hecurve_g1_bsgs_search and also by hecurve_g1_dlog.
*/

#if FF_HALF_WORDS == 1
#define BSGS_MAX_STEPS		256					// guaranteed to handle  p<2^31 (single word ff_t), assuming 2 torsion is known, enlarge if needed but it's nice to fit in L1 cache
#else
#define BSGS_MAX_STEPS		1024				// will BSGS up to 2^40, assuming 2 torsion is known, and nearly as high for dlogs
#endif
#define BSGS_TABSIZE			BSGS_MAX_STEPS		// don't make this too big, it takes time to initialize it.  A few collisions won't kill us.
#define BSGS_TABMASK			((unsigned long)(BSGS_TABSIZE-1))

static hec_jc_t babys[BSGS_MAX_STEPS];
static hec_jc_t giants[BSGS_MAX_STEPS];

static ff_t stepzs[2*BSGS_MAX_STEPS];
#if FF_HALF_WORDS == 1
static unsigned char hashtab[BSGS_TABSIZE];
#else
static unsigned short hashtab[BSGS_TABSIZE];
#endif
static struct tab_entry {
	ff_t x;
	short i;
	short next;
} entries[BSGS_MAX_STEPS+1];
short nexttabentry;

static inline void tab_clear() { memset(hashtab,0,sizeof(hashtab)); nexttabentry = 1; } // don't use entry 0

// we require inserts to have unique x values, return -1 if unique, otherwise return index value for existing entry
static inline int tab_insert(ff_t x, short i)
{
	register struct tab_entry *p;
	register int n,h;

	h = x&BSGS_TABMASK;
	n = hashtab[h];
	if ( n ) {
		p=entries+n;
		for(;;) {
			if ( _ff_equal(p->x,x) ) return p->i;
			if ( ! p->next ) break;
			p=entries+p->next;
		}
	}
	p = entries+nexttabentry;
	_ff_set(p->x,x);
	p->i = i;
	p->next = n;
	hashtab[h] = nexttabentry++;
	return -1;
}

static inline int tab_lookup(ff_t x)
{
	register struct tab_entry *p;
	register int n;

	n = hashtab[x&BSGS_TABMASK];
	if ( ! n ) return -1;
	p=entries+n;
	for(;;) {
		if ( _ff_equal(p->x,x) ) return p->i;
		if ( ! p->next ) return -1;
		p=entries+p->next;
	}
	return -1;
}

/*
	Computes the discrete log of b wrt a given the order of a using a BSGS search.  Optimized for small searches (shares table lookup code with hecurve_g1_bsgs_search).
	Returns 0 if b is not an element of <a> (this is not assumed).

	Note that we are only called for prime values whose square divides the group order, so we only use O(N^(1/4)) steps in the worst case, and even this is very rare
	since the prime must also divide p-1.
*/
int hecurve_g1_dlog (hec_jc_t *a, hec_jc_t *b, long o, ff_t f1)
{
	hec_jc_t c[1];
	long bsteps, giant, gstep, gsteps;
	ff_t zinv[2], baby_x[1], baby_y[1], giant_x[1], giant_y[1];
	ff_t t0, t1;
	register int i, j;
		
	// hard-wire small o cases
	if ( hecurve_g1_JC_id(a) ) return o;
	switch ( o ) {
	case 1:	return 0;
	case 2:	return ( hecurve_g1_JC_cmp (a,b) ? 1 : 0 );
	case 3: 	i = hecurve_g1_JC_cmp (a,b);  return ( i ? (i>0?1:2) : 0 );
	}

//printf ("dlog(%d) a=(%d,%d,%d,%d) b=(%d,%d,%d,%d)\n", o, _ff_get_ui(a->x), _ff_get_ui(a->y), _ff_get_ui(a->z2), _ff_get_ui(a->z3), _ff_get_ui(b->x), _ff_get_ui(b->y), _ff_get_ui(b->z2), _ff_get_ui(b->z3));
	
	// For small searches we will just search the whole range so we only match once (using just one inversion).   The only optimization we use is the usual one for inverses.
	bsteps = (long)sqrt(o/2);
	gstep = 2*bsteps+1;
	gsteps = _ui_ceil_ratio(o,gstep);
	if ( bsteps > BSGS_MAX_STEPS || gsteps > BSGS_MAX_STEPS ) { printf ("%ld: order %ld caused BSGS_MAX_STEPS to be exceeded in hecurve_g1_dlog\n", _ff_p, o); exit (0); }
//printf("bsteps=%d, gsteps=%d, giant step=%d\n", bsteps, gsteps, gstep);
	
	// Convert baby and giant step to affine coords for stepping
	hecurve_g1_JC_exp_ui (c,a,gstep,f1);
	hecurve_g1_JC_invert(c);
	_ff_set(zinv[0],a->z3);
	_ff_set(zinv[1],c->z3);
	ff_parallel_invert (zinv, zinv, 2);
	hecurve_g1_JC_to_A (baby_x, baby_y, a, zinv[0]);
	hecurve_g1_JC_to_A (giant_x, giant_y, c, zinv[1]);
//printf ("baby step(a) = (%d,%d)\n", _ff_get_ui(baby_x[0]), _ff_get_ui(baby_y[0]));
//printf ("giant step(a^-%d) = (%d,%d)\n", gstep, _ff_get_ui(giant_x[0]), _ff_get_ui(giant_y[0]));

	// Step
	babys[0] = *a;  giants[0] = *b;
	j = hecurve_g1_AJC_steps_1 (babys, baby_x[0], baby_y[0], bsteps, f1);							// baby steps
	if ( j < bsteps ) { printf ("%ld: baby step hit identity in dlog with o=%ld\n", _ff_p, o); exit(0); }
	j = hecurve_g1_AJC_steps_1 (giants, giant_x[0], giant_y[0], gsteps, f1);						// giant steps
	if ( j < gsteps ) return j*gstep;

	// Convert to affine to get unique x coords
	for ( i = j = 0 ; i < bsteps ; i++, j++ ) _ff_set(stepzs[j], babys[i].z2);
	for ( i = 0 ; i < gsteps ; i++, j++ ) _ff_set(stepzs[j], giants[i].z2);
	ff_parallel_invert(stepzs, stepzs, j);																	// we assume j < FF_MAX_PARALLEL_INVERTS here
	for ( i = j = 0 ; i < bsteps ; i++,j++ ) ff_mult(babys[i].x, babys[i].x, stepzs[j]);								// only invert x-coord here, we never need to invert y
	for ( i = 0 ; i < gsteps ; i++,j++ ) ff_mult(giants[i].x, giants[i].x, stepzs[j]);

/*
printf ("affine baby x's:\n", j);
for ( i = 0 ; i < bsteps ; i++ ) printf ("   %d: %ld\n", i+1, _ff_get_ui(babys[i].x));
printf ("affine giant x's:\n", j);
for ( i = 0 ; i < gsteps ; i++ ) printf ("   %ld: %ld\n", i*gstep, _ff_get_ui(giants[i].x));
*/
	// Populate the table with babys.  Inserts should never fail (baby steps can't be inverses provided bsteps < o/2)
	tab_clear();
	for ( i = 0 ; i < bsteps ; i++ ) if ( (j=tab_insert(babys[i].x,i)) >= 0 ) break;
	if ( i < bsteps ) { printf("%ld: baby step insert failed in dlog\n", _ff_p); exit (0); }

	// Now match giant steps
	for ( i = 0 ; i < gsteps ; i++ ) {
		if ( (j=tab_lookup(giants[i].x)) >= 0 ) {
			_ff_mult(t0,babys[j].y,giants[i].z3);	_ff_mult(t1,giants[i].y,babys[j].z3);								// normalize y values for comparison
			if ( _ff_equal(t0,t1) ) return i*gstep+j+1;
			return (i?i*gstep-j-1:o-j-1);
		}
	}
	return 0;
}


/*
	Support functions for hecurve_g1_bsgs_search...
*/

double baby_stretch[8] = {0.0,1.0,2.0,1.5,3.0,0.0,3.0,6.0};

// returns the (j+1)-th integer satisfying property determined by modcase.
static inline int baby_index(int j, int modcase)
{
	register int k;
	
	k=j&1;
	switch (modcase) {
	case 1: return j+1;				// all values
	case 2: return 2*j+1;				// odd values	(bspan,gsteps and giants even)
	case 3: return 3*(j-k)/2+1+k;		// prime to 3    (bspan, gstep, and giants multiples of 3)
	case 4: return 3*(j+1);				// 0 mod 3        (bspan and gstep multiples of 3, giants a mod 3)
	case 6: return 6*(j-k)/2+1+4*k;		// prime to 6    (bspan, gstep, and giants multiples of 6)
	case 7: return 6*(j+1);				// 0 mod 6        (bspan and gstep multiples of 6, giants a mod 6)
	}
}

/*
	sets exp to the unique multiple of k in [low,high] if there is one (and returns 1) or sets exp=k and returns 0
*/
static inline int set_exp (long *exp, long k, long low, long high)
{
	register int i;
	i = ui_multiples_in_range(k,low,high);
	if ( ! i ) { printf ("no multiple of element order %ld in interval [%ld,%ld]\n", k, low, high); exit(0); }
	if ( i == 1 ) { *exp = (high/k)*k;  return 1; } else { *exp = k; return 0; }
}


/*
	Uses BSGS (and fastorder) to compute the order k of the element b, given that some multiple of k lies in [low,high] (this is assumed, not verified).
	The values m and a specify modular constraints on k as follows:

		m=1		no constraint
		m=2		k is odd						(use odd baby steps, even giant steps)
		m=3		k is prime to 3					(use baby steps 1,2 mod 3, giant steps 0 mod 3)
		m=3 (a!=0)	k = a mod 3					(use baby steps 0 mod 3, giant steps congruent to a mod 3)
		m=6		k is prime to 6					(use baby steps 1,5 mod 6, giant steps 0 mod 6)
		m=6 (a!=0)	k = a mod 6					(use baby steps 0 mod 6, giant steps congruent to a mod 6)

	The value modcase = m+(a?1:0) is used internally to uniquely identify these 6 situations.
	Note that the giant steps are congruent to a mod m in every case.

	If the return value is 1, then exp is the unique multiuple of k in [low,high].
        If the return value is 0, exp=k is the order of b.  This almost always means there is more than one multiple of k in [low,high], but not necessarily.
*/
int hecurve_g1_bsgs_search (long *exp, hec_jc_t *b, long low, long high, int m,int a, ff_t f1)
{
	hec_jc_t p[64];
	hec_jc_t baby_steps[3], giant_steps[1];
	ff_t zinv[4], baby_x[3], baby_y[3], giant_x[1], giant_y[1];
	register ff_t t0, t1, t2;
	register int i,j,k,bsteps,gsteps,dsteps,usteps;
	int bspan,gstep,modcase;
	long o, o1, o2, gbase, range;
	
//printf ("%d: hecurve_g1_bsgs pt (%ld,%ld,%ld,%ld), low=%d, high=%d, m=%d, a=%d\n", _ff_p, _ff_get_ui(b->x), _ff_get_ui(b->y),_ff_get_ui(b->z2),_ff_get_ui(b->z3), low, high, m, a); 
	
	range = high-low+1;												// range is inclusive, we assume high>low, so range > 1
	bsteps = (long)sqrt((double)(range) / (2.0*baby_stretch[m]));
	if ( ! bsteps ) bsteps = 1;
	if ( m>2 && !a && (bsteps&1) ) bsteps++;								// make sure bsteps is even if we are covering steps prime to 3 or 6
	modcase = m + (a?1:0);
	bspan = bsteps * baby_stretch[modcase];								// effective baby coverage depends on modular contraints m and a
	if ( m > 1 ) bspan = m * _ui_ceil_ratio(bspan,m);							// make sure bspan is a multiple of modulus m (we can safely round up)
	gstep= 2*bspan;													// giant span is effectively 2*bspan+1 due to inverses, but we use 2*bspan so bspan divides gspan...
	if ( m==1 ) gstep++;												// ...except when m=1, where we can use 2*bspan+1.

	/*
		Pick first giant step to be a multiple of as large a value of 2 as possible (and still lie in the interval [low,high]).
		We then need to adjust to make it congruent to a mod m, but this changes at most 2 bits.
		(An optimal approach might search for the value in [low,high] congruent to a mod m with minimal hamming weight, but this would likely take more time than it is worth).
		
		The overall savings is only a few percent (it cuts the cost of computing the first giant step by 1/6, but computing the first giant step is less than 1/4 the total cost).
		However, it's an easy optimization, so there is no reason not to do it.
		
		This idea was suggested by Dan Bernstein.
	*/
	k = ui_lg_floor(range);
	for (;;) {
		o1 = (1L<<(k+1));  o = o1*_ui_ceil_ratio(low,o1);
		if ( o > high ) break;
		k++;
	}
	o1 = 1L<<k;
	gbase = o1*_ui_ceil_ratio(low,o1);										// gbase is now a multiple of a largish power of 2 in [low,high]
	if ( m > 1 ) {
		i = (gbase-a)%m;
		gbase -= i;													// gbase is now congruent to a mod m
		if ( gbase < low ) gbase += m;									// make sure we stay in [low,high]
		if ( gbase > high ) gbase -= m;
	}
	for ( dsteps = (gbase-low)/gstep ; gbase - (dsteps-1)*gstep - bspan > low ; dsteps++ );
	while ( gbase-(dsteps-1)*gstep <= 0 ) gbase += m;						// make sure we don't step on zero or negative values
	for ( usteps = (high-gbase)/gstep ; gbase + (usteps-1)*gstep + bspan < high ; usteps++ );
	if ( dsteps+usteps+1 > BSGS_MAX_STEPS )  { printf ("Exceeded BSGS_MAX_STEPS=%d! p=%ld, low=%ld, high=%ld, m=%d, a=%d\n", BSGS_MAX_STEPS, _ff_p, low, high, m, a);  exit (0); }
	gsteps = dsteps+usteps-1;											// dsteps and usteps both include starting step, so total is one less then the sum
	
//printf ("%d: bsteps=%d, gsteps=%d(%d,%d), bspan=%d, gstep=%d, first giant = %ld\n", _ff_p, bsteps, gsteps, dsteps, usteps, bspan, gstep, gbase);
	
	// compute 2^k powers of b sufficient to compute b^gstep, and b^gbase
	o = _ui_max(gstep,gbase);
	k = ui_lg_floor(o);
	p[0] = *b;
	hecurve_g1_JC_powers (p, k+1, f1);									// add 1 extra power for NAF
	baby_steps[0] = *b; k = 1;
	
	// compute baby step sizes depending on the value of modcase
	switch (modcase) {																				// if m=1, fist step=step size=1
	case 1: break;
	case 2:
	case 3: hecurve_g1_2JC (baby_steps+1, baby_steps, f1); k = 2; break;										// first step is 1, steps size is 1 or 1,2
	case 4: hecurve_g1_2JC (baby_steps+1,baby_steps,f1);  hecurve_g1_JCJC(baby_steps+1,baby_steps,f1); k = 2;  break;	// first step = step size is 3
	case 6: hecurve_g1_2JC (baby_steps+1, baby_steps, f1);											
		     hecurve_g1_2JC (baby_steps+2, baby_steps+1, f1); k = 3; break;										// first step is 1, step size is 2 or 4
	case 7: hecurve_g1_2JC (baby_steps+1,baby_steps,f1);  hecurve_g1_JCJC(baby_steps+1,baby_steps,f1);
		     hecurve_g1_2JC(baby_steps+1,baby_steps+1,f1); k = 2;  break;										// first step = step size is 6
	default:
		printf("Invalid modcase value %d in hecure_g1_bsgs\n", m); exit (0);
	}

	hecurve_g1_JC_exp_powers(giant_steps, p, gstep, f1);
	hecurve_g1_JC_exp_powers(giants+dsteps-1, p, gbase, f1);
	
	// If any of our baby steps is the identity, we won't get very far, so we need to check this situation (it can happen).   The giant steps are handled below.
	if ( k > 1 && hecurve_g1_JC_id(baby_steps+1) ) return set_exp(exp,(a?m:2),low,high);
	if ( k > 2 && hecurve_g1_JC_id(baby_steps+2) ) return set_exp(exp,4,low,high);
	
	// We need to convert the steps to affine coords before using them--do this with one field inversion.
	for ( i = 0 ; i < k ; i++ ) _ff_set(zinv[i], baby_steps[i].z3);
	j = k;
	if ( ! hecurve_g1_JC_id(giant_steps) ) { _ff_set(zinv[j], giant_steps[0].z3); j++; }						 			// giant step could be the identity
	// invert all the z's
	ff_parallel_invert (zinv, zinv, j);
	for ( i = 0 ; i < k ; i++ ) hecurve_g1_JC_to_A (baby_x+i, baby_y+i, baby_steps+i, zinv[i]);
	j=k;
	if ( ! hecurve_g1_JC_id(giant_steps) ) { hecurve_g1_JC_to_A (giant_x, giant_y, giant_steps, zinv[j]); j++; }

	// Now handle giant step = identity (we needed to get the baby in affine coords for the fastorder computation first)
	if ( hecurve_g1_JC_id(giant_steps) ) { o = hecurve_g1_fastorder (baby_x[0], baby_y[0], gstep, f1); return set_exp(exp,o,low,high); }
	if ( hecurve_g1_JC_id(giants+dsteps-1) ) { o = hecurve_g1_fastorder (baby_x[0], baby_y[0], gbase, f1); return set_exp(exp,o,low,high); }

/*
printf("affine baby step 0 (%ld,%ld)\n", _ff_get_ui(baby_x[0]),_ff_get_ui(baby_y[0]));
if (m>1) printf("affine baby step 1 (%ld,%ld)\n", _ff_get_ui(baby_x[1]),_ff_get_ui(baby_y[1]));
if (m==6) printf("affine baby step 2 (%ld,%ld)\n", _ff_get_ui(baby_x[2]),_ff_get_ui(baby_y[2]));
printf("affine giant step %d: (%ld,%ld)\n", gstep, _ff_get_ui(giant_x[0]),_ff_get_ui(giant_y[0])); 
*/
		
	// Baby steps
	babys[0] = ( a ? baby_steps[1] : baby_steps[0] );														// take first baby step
	switch (modcase) {
	case 1:  j=hecurve_g1_AJC_steps_1 (babys, baby_x[0], baby_y[0], bsteps, f1);  break;						// step on everything (use step size 1)
	case 2:  j=hecurve_g1_AJC_steps_1 (babys, baby_x[1], baby_y[1], bsteps, f1);  break;						// step on odd powers (use step size 2)
	case 3:  j=hecurve_g1_AJC_steps_2 (babys, baby_x[0], baby_y[0], baby_x[1], baby_y[1], bsteps, f1);  break;		// step on powers prime to 3 (alternate step sizes 1 and 2: 1,2,4,5,7,8,...)
	case 4:  j=hecurve_g1_AJC_steps_1 (babys, baby_x[1], baby_y[1], bsteps, f1);  break;						// step on multiples of 3
	case 6:  j=hecurve_g1_AJC_steps_2 (babys, baby_x[2], baby_y[2], baby_x[1], baby_y[1], bsteps, f1);  break;		// step on powers prime to 6 (alternate step sizes 4 and 2: 1,5,7,11,13,17,...)
	case 7:  j=hecurve_g1_AJC_steps_1 (babys, baby_x[1], baby_y[1], bsteps, f1);  break;						// step on multiples of 6
	}
	if ( j < bsteps ) { 																				// if a baby step is the identity, we know the order of the element
		k = baby_index(j,modcase);																	// note that we rely on the order being relatively prime to 2,3,6 (per m)
		if ( a ) k/=m;  																				// need to adjust for common divisor of baby steps when a!=0
		return set_exp(exp,k,low,high);
	}	
	
	// Giant steps
	j = hecurve_g1_AJC_dsteps_1 (giants+dsteps-1, giant_x[0], giant_y[0], dsteps, f1);							// downward steps first
	if ( j < dsteps ) {																				// if a giant step hit the identity, it could be a multiple of the element order
		// It might be faster to just continue in this situation, but it's rare in any event
		o = hecurve_g1_fastorder (baby_x[0], baby_y[0], gbase-j*gstep, f1);
		return set_exp(exp,o,low,high);
	}
	j = hecurve_g1_AJC_steps_1 (giants+dsteps-1, giant_x[0], giant_y[0], usteps, f1);							// now upward steps
	if ( j < usteps ) {																				// if a giant step hit the identity, it could be a multiple of the element order
		// It might be faster to just continue in this situation, but it's rare in any event
		o = hecurve_g1_fastorder (baby_x[0], baby_y[0], gbase+j*gstep, f1);
		return set_exp(exp,o,low,high);
	}
	gbase -= (dsteps-1)*gstep;																		// reset gbase to match giants[0]
	
	// Convert to affine to get unique x coords
	for ( i = j = 0 ; i < bsteps ; i++, j++ ) _ff_set(stepzs[j], babys[i].z2);
	for ( i = 0 ; i < gsteps ; i++, j++ ) _ff_set(stepzs[j], giants[i].z2);
	ff_parallel_invert(stepzs, stepzs, j);																	// we assume j < FF_MAX_PARALLEL_INVERTS here
	for ( i = j = 0 ; i < bsteps ; i++,j++ ) ff_mult(babys[i].x, babys[i].x, stepzs[j]);								// only invert x-coord here, we never need to invert y
	for ( i = 0 ; i < gsteps ; i++,j++ ) ff_mult(giants[i].x, giants[i].x, stepzs[j]);
/*
printf ("affine baby x's:\n", j);
for ( i = 0 ; i < bsteps ; i++ ) printf ("   %d(%d): %ld\n", i, baby_index(i,modcase), _ff_get_ui(babys[i].x));
printf ("affine giant x's:\n", j);
for ( i = 0 ; i < gsteps ; i++ ) printf ("   %ld: %ld\n", gbase+i*gstep, _ff_get_ui(giants[i].x));
*/
	// Populate the table with babys.  Insert will fail if we try to insert the same x value twice - this can only happen if two babys are inverses (they can't be equal because we didn't hit the identity).
	tab_clear();
	for ( i = 0 ; i < bsteps ; i++ ) if ( (j=tab_insert(babys[i].x,i)) >= 0 ) break;
	/*
	   If we encountered two baby steps with the same x value, they must be inverses, since otherwise we would have hit the identity in between.
	   If m=1,2,3,6 we have stepped on every possible element order between the two, and if m=4,5 implies we have stepped on every multiple of a common divisor of the two,
	   and hence on every possible value of the difference of their indices (which would be the identity if the steps were equal).
	*/
	if ( i < bsteps ) {
		k = baby_index(i,modcase)+baby_index(j,modcase); 
		if ( a ) k/=m; 																				// need to adjust for common divisor of baby steps when a != 0
		return set_exp(exp,k,low,high);
	}

	// Now match giant steps
	o1 = o2 = 0;
	for ( i = 0 ; i < gsteps ; i++ ) {
		if ( (j=tab_lookup(giants[i].x)) >= 0 ) {
			_ff_mult(t0,babys[j].y,giants[i].z3);	_ff_mult(t1,giants[i].y,babys[j].z3);								// normalize y values for comparison
			if ( _ff_zero(t0) ) { k = 2*baby_index(j,modcase);  return set_exp(exp,k,low,high); }					// handle 2-torsion case
			k = baby_index(j, modcase);
			if ( _ff_equal(t0,t1) ) k = -k;
			o = gbase+i*gstep+k;
			if ( o==o1 ) continue;																	// ignore duplicate matches, this can happen for a != 0
			if ( o1 )  { o2 = o; break; }
			o1 = o;
		}
	}
	if ( ! o1 ) { printf ("%ld: Error, BSGS search of entire Hasse-Weil interval failed, low=%ld, high=%ld, m=%d\n", _ff_p, low, high, m);  exit (0); }
	if ( ! o2 ) { *exp = o1; return 1; }									// we know the o1 is the unique multiple of the element order in [low,high]
	k=i_abs(o2-o1);												// otherwise, we know |o2-o1| is a multiple of the element order and there are at least two multiples in [low,high]
	*exp = (m==1 ? k : hecurve_g1_fastorder (baby_x[0], baby_y[0], k, f1) );	// if m is not 1 we don't know that o1 and o2 are adjacent multiples, but a fastorder computation will figure it out
	return 0;
}

/*
	Compute the order of affine point (x,y) given that (x,y)^e=1 using the classical algorithm.

	This algorithm is not particularly fast, but it doesn't need to be, as it is not invoked for most p.
	The average number of bits exponentiated per prime for a non-CM curve (fixed curve, varying p) is less than 1.
*/
long hecurve_g1_fastorder (ff_t x, ff_t y, long e, ff_t f1)
{
	hec_jc_t b[MAX_UI_PP_FACTORS];
	unsigned long q, p[MAX_UI_PP_FACTORS], h[MAX_UI_PP_FACTORS];
	register long o;
	register int i, j, k, w;

	w = ui_factor(p,h,e);				// takes about 2 microseconds for e~2^32 (2.5MHz AMD-64), but e is generally much smaller
	if ( w == 1 && h[0] == 1 ) return e;	// not hard when the exponent is prime
		
	for ( i = 0 ; i < w ; i++ ) {
		for ( q=p[i],j=1 ; j < h[i] ; j++ ) q *= p[i]; 
		o = e/q;
		hecurve_g1_AJC_exp_ui (b+i, x, y, o, f1);
	}
	o = 1;
	for ( i = 0 ; i < w ; i++ ) {
		for ( j = 0 ; j < h[i]-1 ; j++ ) {
			if ( hecurve_g1_JC_id(b+i) ) break;
			hecurve_g1_JC_exp_ui (b+i, b+i, p[i], f1);
		}
		if ( ! hecurve_g1_JC_id (b+i) ) j++;
		for ( k = 0 ; k < j ; k++ ) o *= p[i];
	}
 //printf ("fastorder computed order of (%ld,%ld) is %ld from exponent %ld\n", _ff_get_ui(x), _ff_get_ui(y), o, e);
	return o;
}

/*
	BSGS search hardwired for very short intervals (typically less than 50), which uses only one or two baby steps.
	These intervals can arise when the first order computation does not uniquely determine the exponent of the group (or when p is small).

	We assume b does not have 2 torsion (hence its order is at least 3) and high > low.
	We rely on the existence of an exponent of b in [low,high] and do not always verify this (e.g., if high=low+1 and b^low != id then we conclude b^high == id)
*/
int hecurve_g1_bsgs_short (long *exp, hec_jc_t *b, long low, long high, ff_t f1)
{
	long range;
	hec_jc_t b2[1],b5[1],g[1];
	register ff_t t0, t1;
	register long e;
	long o1,o2;
	int i;

	range = high-low+1;
	switch (range) {
	case 2:		// E no computation other than the g exponentiation, cost E
		hecurve_g1_JC_exp_ui(g, b, low, f1);
		*exp = ( hecurve_g1_JC_id(g) ? low : high );
		return 1;
	case 3:		// E+4/3M (on average)
		hecurve_g1_JC_exp_ui(g, b, low+1, f1);
		if ( hecurve_g1_JC_id(g) ) { *exp = low+1;  return 1; }
		// it must be the case that g = +/- b, we don't verify this but simply compare y values to distinguish
		_ff_mult(t0,b->y,g->z3);  _ff_mult(t1,g->y,b->z3);
		*exp =  ( _ff_equal (t0, t1) ? low : high );
		return 1;
	case 4:		// E+1/4D+2M = E+5M+2A (roughly) (D indicates doubling(squaring) a point on the curve, which costs 11M+8A field operations)
		hecurve_g1_JC_exp_ui(g, b, low+1, f1);
		if ( hecurve_g1_JC_id(g) ) { *exp = low+1;  return 1; }
		_ff_mult(t0,b->x,g->z2);  _ff_mult(t1,g->x,b->z2);
		if ( _ff_equal (t0, t1) ) {
			_ff_mult(t0,b->y,g->z3);  _ff_mult(t1,g->y,b->z3);
			if ( ! _ff_equal (t0, t1) ) { *exp = low+2;  return 1; }
			// we know b^low = id, but it could also be that b^high = id if b has order 3, so we need to check this
			hecurve_g1_2JC (b2,b, f1);
			_ff_mult(t0,b->x,b2->z2);  _ff_mult(t1,b2->x,b->z2);					// given b != id, b^3==id iff b and b^2 have the same x coord
			if ( _ff_equal(t0,t1) ) { *exp = 3; return 0; }
			*exp = low;
			return 1;
		}
		*exp = high;
		return 1;
	case 5:		// < E+4/5D+5M (on average)
		hecurve_g1_JC_exp_ui(g, b, low+2, f1);
		if ( hecurve_g1_JC_id(g) ) { *exp = low+2;  return 1; }
		hecurve_g1_2JC (b2,b, f1);
		if ( hecurve_g1_JC_id(b2) ) { *exp = 2; return 0; }
		if ( hecurve_g1_JC_2tor(b2) ) { *exp = 4; return 0; }
		_ff_mult(t0,b->x,b2->z2);  _ff_mult(t1,b2->x,b->z2);
		if ( _ff_equal(t0,t1) ) { *exp = 3; return 0; }
		// we now know |b|>4, so there is a unique multiple in the interval
		_ff_mult(t0,b->x,g->z2);  _ff_mult(t1,g->x,b->z2);
		if ( _ff_equal (t0, t1) ) {
			_ff_mult(t0,b->y,g->z3);  _ff_mult(t1,g->y,b->z3);
			*exp = ( _ff_equal(t0,t1) ? low+1 : low+3 );
			return 1;
		}
		_ff_mult(t0,b2->y,g->z3);  _ff_mult(t1,g->y,b2->z3);
		*exp =  ( _ff_equal (t0, t1) ? low : low+4 );
		return 1;
	}
	hecurve_g1_2JC (b2,b, f1);
	if ( hecurve_g1_JC_id(b2) ) { *exp = 2; return 0; }
	if ( hecurve_g1_JC_2tor (b2) ) { *exp = 4; return 0; }	
	_ff_mult(t0,b->x,b2->z2);  _ff_mult(t1,b2->x,b->z2);
	if ( _ff_equal(t0,t1) ) { *exp = 3; return 0; }
	hecurve_g1_2JC (b5,b2, f1);
	hecurve_g1_JCJC (b5,b, f1);
	if ( hecurve_g1_JC_id(b5) ) { *exp = 5; return 0; }
	// we know |b|>5
	e = low+2;
	hecurve_g1_JC_exp_ui(g, b, e, f1);
	o1 = o2 = 0;
	while ( e < high+3 ) {
		if ( hecurve_g1_JC_id(g) ) { o2 = o1; o1 = e; goto next; }
		_ff_mult(t0,b->x,g->z2);  _ff_mult(t1,g->x,b->z2);
		if ( _ff_equal (t0, t1) ) {
			_ff_mult(t0,b->y,g->z3);  _ff_mult(t1,g->y,b->z3);
			i = ( _ff_equal(t0,t1) ? -1 : 1 );
			o2 = o1;  o1 = e+i; goto next;
		}
		_ff_mult(t0,b2->x,g->z2);  _ff_mult(t1,g->x,b2->z2);
		if ( _ff_equal (t0, t1) ) {
			_ff_mult(t0,b2->y,g->z3);  _ff_mult(t1,g->y,b2->z3);
			i = ( _ff_equal(t0,t1) ? -2 : 2 );
			o2 = o1;  o1 = e+i; goto next;
		}
next:	if ( o2 ) break;
		hecurve_g1_JCJC (g,b5,f1);
		e += 5;
	}
	if ( ! o1 ) { printf ("%ld: No match found in short interval [%ld,%ld]\n", _ff_p, low, high);  exit (0); }
	if ( o2 ) { *exp = o1-o2; return 0; }
	*exp = o1;
	// total cost is roughly E + 2D + (1+(range/5-1))*A + (range/5)*4M (D indicates point doubling, A indicates point addition)
	return 1;
}

/*
	The rest of this file contains code that is not currently used because it was found to to be suboptimal in testing (on an AMD Athlon) but
	may be useful on other platforms and/or other applications.
*/
#if 0

/*
	Below are elliptic curve arithmetic functions for Jacobian coordinates.  The reduced Chudnovsky coordinates
        are faster at everything except 2JC versus 2J is (11M+8A) vs (10M+10A), but this difference is negligible (1M~2.5A),
        and in testing it was found to be faster to use JC everywhere.
*/

// converts Jacobian to affine coords, given the inverse of z
static inline void hecurve_g1_J_to_A (ff_t *px, ff_t *py, hec_j_t *p, ff_t zinv)
{
	register ff_t t0,t1;
	
	_ff_square(t0,zinv);
	_ff_mult(*px,p->x,t0);
	_ff_mult(t1,t0,zinv);
	_ff_mult(*py,p->y,t1);
	// 4M
}

// converts Jacobian to reduced Chudnovsky Jacobian coords
static inline void hecurve_g1_J_to_JC (hec_jc_t *o, hec_j_t *p)
{
	_ff_set(o->x,p->x);
	_ff_set(o->y,p->y);
	_ff_square(o->z2,p->z);
	_ff_mult(o->z3,o->z2,p->z);
	// 2M
}

// squares an affine point into Jacobian coords
static inline void hecurve_g1_2AJ (hec_j_t *p, ff_t x, ff_t y, ff_t f1)
{
	register ff_t t0,t1,t2,t3,a,b;

	_ff_add(p->z,y,y);  _ff_mult(t1,p->z,y); _ff_add(t2,t1,t1); _ff_mult(a,x,t2);		// z=2y, a=4xy^2, t1=2y^2, t2=4y^2
	_ff_mult(t3,t1,t2);													// t3=8y^4
	_ff_square(t1,x); _ff_add(b,t1,t1); _ff_addto(b,t1); _ff_addto(b,f1);				// b = 3x^2+f1*1^4
	_ff_square(t0,b); _ff_add(t2,a,a); _ff_sub(p->x,t0,t2);						// x = b^2-2a
	_ff_sub(t1,a,p->x); _ff_mult(t2,t1,b); _ff_sub(p->y,t2,t3);					// y = b(a-x)-8y^4 (note we use the new x here)
	// 6M+9A
}


// squares a point p1 in Jacobian coords (p3 is the output, may be equal to p1)
static inline void hecurve_g1_2J (hec_j_t *p3, hec_j_t *p1, ff_t f1)
{
	register ff_t a, b, c, t0, t1, t2;
	
	_ff_square(t0,p1->x); _ff_add(t1,t0,t0); _ff_addto(t1,t0); 					// t1 = 3x^2
	_ff_square(t0,p1->z); _ff_square(t2,t0);  _ff_mult(t0,f1,t2); _ff_add(b,t1,t0);		// b = 3x^2+f1*z^4	(note that f1=a4 in 13.2.1.c  of HECHECC p.282)
	_ff_add(c,p1->y,p1->y); ff_mult(p3->z,p1->z,c);							// c=2y, z = 2yz
	_ff_mult(t2,c,p1->y); _ff_add(t0,t2,t2);  _ff_mult(a,t0,p1->x);					// a=4xy^2,t2=2y^2
	_ff_add(t0,a,a); _ff_square(t1,b); _ff_sub(p3->x,t1,t0);						// x = b^2-2a
	_ff_square(t0,t2); _ff_add(c,t0,t0);										// c = 8y^4
	_ff_sub(t0,a,p3->x); _ff_mult(t2,t0,b); _ff_sub(p3->y,t2,c);					// y = b(a-x)-c   -- note we use the new x value here
	// 10M+10A
}


// multiplies a point p in Jacobian coords by a (non-identity) affine point (p is an input and an output)
static inline void hecurve_g1_AJ (hec_j_t *p, ff_t x0, ff_t y0, ff_t f1)
{
	register ff_t a, c, e, f, t0, t1, t2;
	
	if ( _ff_zero(p->z) ) { _ff_set(p->x,x0); _ff_set(p->y,y0); _ff_set_one(p->z); return; }
	_ff_square(t0,p->z);  _ff_mult(a,x0,t0);									// a = x0z^2, b = x*1^2=x (since z0=1), and t0=z^2
	_ff_sub(e,p->x,a);													// e = a-b
	_ff_mult(t1,t0,p->z); _ff_mult(c,y0,t1);									// c = y0z^3, d=y*1^3=y
	_ff_sub(f,p->y,c);													// f = d-c
	if ( _ff_zero(e) && _ff_zero(f) ) { hecurve_g1_2AJ (p, x0,y0,f1); return; }		// must use doubling code here, but at least its an affine double (6M+9A)
	_ff_square(t0,e); _ff_mult(t1,t0,e); _ff_mult(t2,t0,a);						// t1=e^3, t2=ae^2
	_ff_square(t0,f); _ff_sub(p->x,t0,t1); _ff_add(t0,t2,t2); _ff_subfrom(p->x,t0);		// x = f^2-e^3-2ae^2
	_ff_sub(t0,t2,p->x); _ff_mult(t2,f,t0); _ff_mult(t0,t1,c); _ff_sub(p->y,t2,t0);		// y = f(ae^2-x) - ce^3
	ff_mult(p->z,p->z,e);												// z = 1*z*e
	// 11M+6A
}


// multiplies p1 in Jacobian coords by p2 Jacobian coords  and puts the result in p1.  Handles identity and will square (double) if needed
// this is the slowest case and should be avoided when possible, AJ, AJC, and JCJC are all better
static inline void hecurve_g1_JJ (hec_j_t *p1, hec_j_t *p2, ff_t f1)
{
	register ff_t a, b, c, d, e, f, t0, t1, t2;
	
	// we need to handle identity cases in general (although in some cases we might be able to rule this out)
	if ( _ff_zero(p2->z) ) return;
	if ( _ff_zero(p1->z) ) { *p1=*p2; return; }
	_ff_square(t0,p2->z);  _ff_mult(a,p1->x,t0);								// a = x1z2^2, and t0=z2^2
	_ff_square(t1,p1->z);  _ff_mult(b,p2->x,t1);								// b = x2z1^2, and t1=z1^2
	_ff_sub(e,b,a);														// e = a-b
	_ff_mult(t2,t0,p2->z); _ff_mult(c,p1->y,t2);								// c = y1z2^3
	_ff_mult(t2,t1,p1->z); _ff_mult(d,p2->y,t2);								// d = y2z1^3
	_ff_sub(f,d,c);														// f = d-c
	if ( _ff_zero(e) && _ff_zero(f) ) { hecurve_g1_2J (p1,p1,f1); return; }			// must double if pts are equal (10M+10A), inverses will end up with z=0, so we let them through
	_ff_square(t0,e); _ff_mult(t1,t0,e); _ff_mult(t2,t0,a);						// t1=e^3, t2=ae^2
	_ff_square(t0,f); _ff_sub(p1->x,t0,t1); _ff_add(t0,t2,t2); _ff_subfrom(p1->x,t0);	// x = f^2-e^3-2ae^2
	_ff_sub(t0,t2,p1->x); _ff_mult(t2,f,t0); _ff_mult(t0,t1,c); _ff_sub(p1->y,t2,t0);	// y = f(ae^2-x) - ce^3	-- note we use the new x here
	ff_mult(t0,p1->z,p2->z);	_ff_mult(p1->z,t0,e);							// z = z1z2e
	// 16M+7A (ouch)
}


// Computes p=(x,y)^n where (x,y) !=1 is in affine coordinates and p is in Jacobian coordinates
// It takes advantage of the fact that all additions are of the form J+A, requiring only 11M, doubling is done in J+J form, using 10M
void hecurve_g1_AJ_exp_ui (hec_j_t *p, ff_t x0, ff_t y0, unsigned long n, ff_t f1)
{
	register unsigned long m;
	unsigned long pbits, nbits;
	register ff_t negy0;
	int i;

	if ( n == 0 ) { _ff_set_zero(p->z); return; }
	if ( n == 1 ) { _ff_set(p->x,x0); _ff_set(p->y,y0); _ff_set_one(p->z); return; }
	hecurve_g1_2AJ(p,x0,y0,f1);
	if ( n == 2 ) return;
	ui_NAF(&pbits,&nbits,n);
	i = _asm_highbit(pbits)-2;						// we know the top two bits of the NAF are 10
hecurve_expbits+=i+1;
	_ff_neg(negy0,y0);
	m = (1UL<<i);
	for ( ; m ; m >>= 1 ) {
		hecurve_g1_2J(p,p,f1);	 					// 10M+10A
		if ( m&pbits ) hecurve_g1_AJ(p,x0,y0,f1);		// 11M+6A
		if ( m&nbits ) hecurve_g1_AJ(p,x0,negy0,f1);		// 11M+6A
	}
}


// Computes p=a^e where a!=1 is in Jacobian coords, and so is p.
// This is slow and should be avoided.  Overlap is ok
void hecurve_g1_J_exp_ui (hec_j_t *p, hec_j_t *a, unsigned long n, ff_t f1)
{
	register unsigned long m;
	unsigned long pbits, nbits;
	hec_j_t t, ai;
	int i;
	
	if ( n == 0 ) { _ff_set_zero(p->z); return; }
	if ( n == 1 ) { *p=*a; return; }
	hecurve_g1_2J(&t,a,f1);							// 10M+10A
	if ( n == 2 ) { *p = t; return; }
	ui_NAF(&pbits,&nbits,n);
	i = _asm_highbit(pbits)-2;
hecurve_expbits+=i+1;
	ai=*a;  ff_negate(ai.y);
	m = (1UL<<i);
	for ( ; m ; m >>= 1 ) {
		hecurve_g1_2J(&t,&t,f1);						// 10M+10A
		if ( m&pbits ) hecurve_g1_JJ(&t,a,f1);			// 16M+7A
		if ( m&nbits ) hecurve_g1_JJ(&t,&ai,f1);			// 16M+7A
	}
	*p = t;
}


// Combines precomputed values p[i]=p[0]^(2^i) to compute p[0]^n.  Assumes all required powers are present: NAF potentially requires one more bit!
// CAUTION: if this is false, the resulting behavior is very strange and unpredictable.  You have been warned.
void hecurve_g1_J_exp_powers(hec_j_t *o, hec_j_t *p, unsigned long n, ff_t f1)
{
	register unsigned long m;
	unsigned long pbits, nbits;
	hec_j_t t;
	int i;

	// handle small n quickly
	switch (n) {
	case 0:	_ff_set_zero(o->z); return; 
	case 1:	*o=p[0]; return;
	case 2:	*o=p[1]; return;
	case 3:	*o=p[0]; hecurve_g1_JJ(o,p+1,f1); return;
	case 4:	*o=p[2]; return;
	case 5:	*o=p[2]; hecurve_g1_JJ(o,p,f1); return;
	case 6:	*o=p[2]; hecurve_g1_JJ(o,p+1,f1); return;
	case 7:	*o=p[0]; ff_negate(o->y); hecurve_g1_JJ(o,p+3,f1); return;
	case 8:	*o=p[3]; return;
	}
	ui_NAF(&pbits,&nbits,n);
	i = _asm_highbit(pbits);
	*o = p[i];
	i-=2;
hecurve_expbits+=i+1;
	m = (1UL<<i);
	for ( ; m ; m >>= 1, i-- ) {
		if ( (m&pbits) ) hecurve_g1_JJ(o,p+i,f1);						 	// 16M+7A (these are expensive)
		if ( (m&nbits) )
			{ ff_negate(o->y); hecurve_g1_JJ(o,p+i,f1); ff_negate(o->y); }	// 16M+9A (slightly more expensive)
	}
}

// Sets p[i] = p[0]^(2^i) for i in [0,k]. Input is an affine pt not equal to the identity, output is in Jacobian coords.
void hecurve_g1_AJ_powers (hec_j_t p[], ff_t x0, ff_t y0, int k, ff_t f1)
{
	ff_t a, b, c, x, y, z, t0, t1, t2;
	register hec_j_t *e;
	
	_ff_set(p->x,x0);  _ff_set(p->y,y0);  _ff_set_one(p->z);
	hecurve_g1_2AJ (p+1, x0, y0, f1);
	for ( e=p+k,p+=2 ; p <= e ; p++ ) hecurve_g1_2J(p,p-1,f1); 			// 10M + 10A
}
#endif
