#ifndef _FF_INCLUDE_
#define _FF_INCLUDE_

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

// Stripped-down, self-contained implementation if 32 or 64 bit single precision
// finite fields using a Montgomery representation.  All you need is ff.h, asm.h, and ff.c

// THIS CODE IS NOT PORTABLE, but it is fast.  It assumes an unsigned is 32 bits and that an
// unsigned long is 64 bits and (most critically) uses  inline assembly directives.  These
// should work on any AMD/Intel compatible instruction set.

// The FF_MONTGOMERY flag should generally be set, the only reason to unset it is for 
// testing or benchmarking.  It can be helpful to find bugs in code that breaks when
// the field representation changes (e.g. when the bit string "000...01" is not the identity element).

#include "asm.h"
#include "smalljac.h"								// only needed to pick word size

#define ULONG_BITS	(8*sizeof(unsigned long))		// must be 32 or 64

#define FF_FAST					1				// define to turn off certain error checking

#define FF_MAX_PARALLEL_INVERTS		2048			// make big enough to not force chunking in hecurve1

#define FF_WORDS					1				// must be 1
#define FF_HALF_WORDS				(SMALLJAC_FF_64BITS?2:1)
#define FF_MONTGOMERY			1				// if 0 FF_HALF_WORDS must be 1
#define FF_NAIL_BITS				1				// must be at least 1 (need to be able to add without overflow)
#define FF_BITS					((FF_HALF_WORDS*ULONG_BITS/2) - FF_NAIL_BITS)
#define FF_MONTGOMERY_RBITS		(FF_HALF_WORDS*ULONG_BITS/2)	// note that R=B are identical - only one word is ever used
#define FF_ITAB_SIZE				(3*FF_MONTGOMERY_RBITS+1)


#if FF_MONTGOMERY && FF_HALF_WORDS == 1
typedef unsigned ff_t;
#else
typedef unsigned long ff_t;
#endif

// temporary value used by macros - could remove by replacing macros with inlines, probably should
extern ff_t _ff_t1;

// The globals below all assume a single modulus environment.  Multiple moduli can be managed by
// copying and resetting these (provided this is reasonably infrequent).  We don't put them into a
// dynamically allocated context to save the cost of dereferencing a pointer every time we want to
// perform a field operation.

extern ff_t _ff_p;
extern ff_t _ff_2g;					// generator of the Sylow 2-subgroup, necessarily a quadratic non-residue
extern ff_t _ff_2gi;					// inverse of _ff_2g
extern ff_t _ff_2Sylow_tab[64][2];
extern ff_t _ff_3g;					// generator of the Sylow 3-subgroup, necessarily a cubic non-residue
extern ff_t _ff_half;					// 1/2 = (p+1)/2
extern ff_t _ff_third;				// 1/3 = (p+1)/3 or (2p+1)/3 for p=2mod3 or 1mod3 (resp.)
extern ff_t _ff_negone;
extern int _ff_p2_e;
extern unsigned long _ff_p2_m;		// p = 2^e*m+1
extern int _ff_p3_e;
extern unsigned long _ff_p3_m;		// p = 3^e*m+1 when p=1mod3
extern int _ff_p3_m1mod3;
extern ff_t *_ff_mont_itab;
extern ff_t _ff_mont_R;
extern ff_t _ff_mont_R2;
extern unsigned long _ff_mont_pni;
extern int _ff_p1mod3;

extern int _ff_cbrt_setup;
extern ff_t _ff_cbrt_unity;			// MUST CALL ff_cbrt_setup() or ff_cbrt() to initialize.  Set to 1 if p=2mod3

void ff_setup_ui (unsigned long p);	// note - DOES NOT CHECK PRIMALITY


// WARNING - several of the macros below are decidedly unsafe.  In general, if it starts with an underscore,
// don't use side effects, and be aware that many can't be used in expressions.
// These could be converted to inline functions, but many involve assignments and would require passing
// by reference and then dereferencing a pointer (this is handled by reference types in C++). 
// A smart compiler could produce inline function code that is just as fast as these macros, but we won't rely on this.

// Input/Output
#define _ff_sprint(s,x)				sprintf(s,"%lu", _ff_get_ui(x))
#define _ff_set_raw(z,x)				((z) = (x))								// useful for setting random field elements.  assumes 0<=x<p
#if FF_MONTGOMERY
#define _ff_set_ui(z,x)				((z) = ff_montgomery1_mult(((x)%_ff_p), _ff_mont_R2));
#define _ff_rset_ui(z,x)				((z) = ff_montgomery1_mult((x), _ff_mont_R2))	// assumes 0<=x<p
#define _ff_get_ui(x)				(ff_montgomery1_mult((x),1UL))
#else
#define _ff_set_ui(z,x)				((z) = ((x)%_ff_p))
#define _ff_rset_ui(z,x)				((z) = (x))			 					// assumes 0<=x<p
#define _ff_get_ui(x)				(x)
#endif
// signed conversions - you must use these for negative values
#define _ff_get_i(x)						((long)_ff_get_ui(x))
#define _ff_set_i(x,z)					if ( (z) < 0 ) { _ff_set_ui(x,-(z));  ff_negate(x); } else { _ff_set_ui(x,z); }

// Basic operations and 0/1 - IMPORTANT: the binary value 1 is not the field identity in a Montgomery representation (but 0 is zero)
#define _ff_set_zero(z)				((z)=0)
#define _ff_zero(x)					(!(x))
#define _ff_nonzero(x)				(x)
#define _ff_parity(x)					((x)&0x1)
#define _ff_set(z,x)					((z)=(x))
#define _ff_equal(x,y)				((x) == (y))
#define _ff_low_word(x)				(x)
#if FF_MONTGOMERY
#define _ff_set_one(z)				_ff_set(z,_ff_mont_R)
#define _ff_one(z)					_ff_equal(z,_ff_mont_R)
#else
#define _ff_set_one(z)				((z)=1UL)
#define _ff_one(z)					((z)==1UL)
#endif
// end basic operations and 0/1

// Core arithmetic operations - these may be applied to unreduced values that fit in the word limit
#define _ff_core_addto(z,x)			((z) += (x))
#define _ff_core_subfrom(z,x)			((z) -= (x))											// assumes z dominates x
#define _ff_core_shiftl(z)				((z) <<= 1)
#define _ff_core_shiftr(z)				((z) >>= 1)
#define _ff_core_ge(x,y)				((x) >= (y))
#define _ff_addto(z,x)				{register ff_t _ff_t; _ff_set(_ff_t,z); _ff_core_addto(_ff_t,x);_ff_core_red(_ff_t); _ff_set(z,_ff_t);}
#define _ff_add(z,x,y)				{_ff_set(z,x);_ff_addto(z,y);}
#define _ff_subfrom(z,x)				{register ff_t _ff_t; _ff_set(_ff_t,z); _ff_core_dom(_ff_t,x); _ff_core_subfrom(_ff_t,x); _ff_set(z,_ff_t);}
#define _ff_sub(z,x,y)				{_ff_set(z,x);_ff_subfrom(z,y);}
#if FF_MONTGOMERY
#define _ff_core_inc(z)				_ff_core_addto(z,_ff_mont_R)
#else
#define _ff_core_inc(z)				((z)++)
#endif

// derived arithmetic operations - note that inputs cannot overlap outputs! use higher level versions if needed
// Note that core_red requires z < 2p
#define _ff_core_red(z)				if (_ff_core_ge(z,_ff_p) ) _ff_core_subfrom (z,_ff_p)
#define _ff_core_dom(z,x)			if ( !_ff_core_ge(z,x) ) _ff_core_addto (z,_ff_p)
#define _ff_neg(z,x)					if (_ff_nonzero(x) ) {_ff_set(z,_ff_p); _ff_core_subfrom(z,x); } else { _ff_set_zero(z); }
#define _ff_x2(z)					_ff_addto(z,z);		// this is much faster than shifting
#define _ff_inc(z)					{_ff_core_inc(z);_ff_core_red(z);}
#if FF_MONTGOMERY
#define _ff_dec(z)					_ff_subfrom((z),_ff_mont_R)
#else
#define _ff_dec(z)					if (z) { (z)--; } else { (z) = _ff_p-1; }
#endif
// end core arithmetic operations

// Safer versions - overlap ok, but still need to watch side effects and can't use in expressions
#define ff_negate(z)				if (_ff_nonzero(z) ) {_ff_set(_ff_t1,_ff_p); _ff_core_subfrom(_ff_t1,z);  _ff_set(z,_ff_t1); } else { _ff_set_zero(z); }
#define ff_add(z,x,y)				{_ff_set(_ff_t1,x);_ff_addto(_ff_t1,y);_ff_set(z,_ff_t1);}
#define ff_sub(z,x,y)					{_ff_set(_ff_t1,x);_ff_subfrom(z,y);_ff_set(z,_ff_t1);}

// higher arithmetic operations
#if FF_MONTGOMERY
#define _ff_mult(z,x,y)				((z) = ff_montgomery1_mult(x,y))
#define _ff_square(z,x)				_ff_mult(z,x,x)
#define _ff_invert(z,x)				((z)=ff_montgomery1_invert(x))
#else
#define _ff_mult(z,x,y)				((z) = ((x)*(y)) % _ff_p)
#define _ff_square(z,x)				_ff_mult(z,x,x)
#define _ff_invert(z,x)				((z) = ff_ui_inverse(x,_ff_p))
#endif
#define _ff_div2(z,x)					_ff_mult(z,x,_ff_half);

#define _ff_incmult(z,x,w)				{ _ff_set(z,x);_ff_inc(z);_ff_mult(z,z,w); }	// could be optimized
#define _ff_multadd(z,x,y,a)			{_ff_mult(z,x,y); _ff_addto(z,a); }		// no optimization here
#define ff_multadd(z,x,y,a)			_ff_multadd(z,x,y,a)					// overlap of x and z ok but a can't overlap

unsigned long ff_montgomery1_invert (unsigned long x);
unsigned long ff_ui_inverse (unsigned long a, unsigned long m);
void ff_exp_ui (ff_t o[1], ff_t a[1], unsigned long e);
int ff_ui_legendre (unsigned long a, unsigned long b);
int ff_invsqrt (ff_t o[1], ff_t a[1], int ext);							// returns 1 if sqrt is rational, 0 if sqrt is in F_p^2.  In the later case, if ext is true, sqrt(a)=o*z where z^2=nr
void ff_setup_2g (void);	
int ff_cbrt (ff_t o[1], ff_t a[1]);										// returns 1 for success, 0 if a is a non-residue
int ff_3Sylow_invcbrt (ff_t o[1], ff_t a[1]);								// computes a^{-1/3} for a in the 3-Sylow subgroup, returns 0 if a is not a residue in the 3-Sylow
int ff_fast_sqrt (ff_t o[1], ff_t a[1]);									// fast sqrt via single exponentiation, returns 1 for success, 0 if no power of a is sqrt(a) (obsolete, not used)

// called automatically by ff_sqrt, but also used in ffext.c - note that for p=3mod4, e=1 and -1 generates the 2-Sylow subgroup
void ff_setup_2g(void);

// ditto for cbrt, only relevant when p=1mod3
void _ff_setup_3g(void);
static inline void ff_setup_3g(void) { if ( ! _ff_3g ) _ff_setup_3g(); }

void ff_parallel_invert (ff_t z[], ff_t x[], unsigned n);					// z and x may overlap

#define ff_mult(z,x,y)				_ff_mult(z,x,y)					// overlap is ok for these macros (but not side effects)
#define ff_square(z,x)				_ff_square(z,x)
#define ff_invert(z,x)				{ if ( ! _ff_one(x) ) _ff_invert(z,x); else _ff_set(z,x); }

// end higher arithmetic operations

// Inline Functions
static inline int ff_is_negation (ff_t x, ff_t y)
{
	register ff_t t;
	
	_ff_set(t,x);
	_ff_core_addto(t,y);
	return ( _ff_equal(t,_ff_p) || _ff_zero(t) );
}

#if FF_MONTGOMERY && FF_WORDS == 1 && FF_HALF_WORDS == 1
static  inline unsigned  ff_montgomery1_mult (unsigned  x, unsigned  y)
{
	register unsigned long z;
	register unsigned a;

	z = (unsigned long)x*y;
	a = z*_ff_mont_pni;
	z  += ((unsigned long)a * _ff_p);
	z >>= 32;
	if ( z >= _ff_p ) z -= _ff_p;
	return z;
}
#endif

#if FF_MONTGOMERY && FF_WORDS == 1 && FF_HALF_WORDS == 2
static  inline unsigned long ff_montgomery1_mult (unsigned long x, unsigned long y)
{
	register unsigned long x0, x1,a0, a1;

	_asm_mult_1_1 (x1,x0,x,y);
	a0 = x0*_ff_mont_pni;
	_asm_mult_1_1 (a1,a0,a0,_ff_p);
	_asm_addto_2_2 (a1,a0,x1,x0);
	if ( a1 >= _ff_p ) a1 -= _ff_p;
	return a1;
}
#endif

static inline int ff_residue (ff_t z) { ff_t t;  if ( _ff_zero(z) ) return 1;  ff_exp_ui(&t,&z,(_ff_p-1)/2);  return _ff_one(t); } // this is faster than using Legendre
static inline int ff_sqrt(ff_t o[1], ff_t a[1]) { ff_t t;  if ( ! ff_invsqrt(&t,a,0) ) return 0;  _ff_mult(o[0],t,a[0]);  return 1; }
static inline int ff_sqrt_ext(ff_t o[1], ff_t a[1]) { register int sts; ff_t t;  sts=ff_invsqrt(&t,a,1);  _ff_mult(o[0],t,a[0]);  return sts; }
int _ff_2Sylow_invsqrt (ff_t o[1], ff_t a[1], int ext);
static inline int ff_2Sylow_invsqrt (ff_t o[1], ff_t a[1], int ext)			// computes a^{-1/2} for a in the 2-Sylow subgroup, returns 0 if a^{-1/2} is not in F_p and returns a^{-1/2}/z
{
	// handle easy cases inline, otherwise call _ff_2Sylow_invsqrt (if e <= 2 this will never happen)
	ff_setup_2g();
	if ( _ff_one(a[0]) ) { _ff_set_one(o[0]);  return 1; }				// use 1 as the inverse square root of 1
	if ( _ff_equal(a[0],_ff_2g) ) { _ff_set(o[0],_ff_2gi);  return 0; }		// if e=1, this must happen
	if ( _ff_equal(a[0],_ff_negone) ) { _ff_set(o[0], _ff_2Sylow_tab[_ff_p2_e-2][1]); return 1;}
	if ( _ff_equal(a[0],_ff_2gi) ) { _ff_set_one (o[0]);  return 0; }			// if e=2, this must happen
	return _ff_2Sylow_invsqrt (o, a,ext);
}		

#endif
