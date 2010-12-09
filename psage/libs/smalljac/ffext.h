#ifndef _FFEXT_INCLUDE_
#define _FFEXT_INCLUDE_

#include "mpzutil.h"		// ui_randomm
#include "ffwrapper.h"

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

// subset of finite field arithmetic supported in degree 2 and 3 extension fields
// elements are represented as polys over F_p[x]/(g(x)) where g(x)=x^2-_ff_2g for degree 2
// and in degree 3, g(x)=x^3-_ff_3g if p=1mod3 else g(x)=x^3-x-_ff_3g

extern int _ff2_cbrt_setup;
extern ff_t _ff2_cbrt_unity[2];

void ff_ext_setup(void);

// for convenience
static inline void ff2_random(ff_t o[2]) { _ff_random(o[0]);  _ff_random(o[1]); }
static inline void ff3_random(ff_t o[3]) { _ff_random(o[0]);  _ff_random(o[1]);  _ff_random(o[2]); }
static inline void ff2_set_zero(ff_t o[2]) { _ff_set_zero(o[0]); _ff_set_zero(o[1]); }
static inline void ff3_set_zero(ff_t o[3]) { _ff_set_zero(o[0]); _ff_set_zero(o[1]); _ff_set_zero(o[2]); }
static inline void ff2_set_one(ff_t o[2]) { _ff_set_one(o[0]); _ff_set_zero(o[1]); }
static inline void ff3_set_one(ff_t o[3]) { _ff_set_one(o[0]); _ff_set_zero(o[1]); _ff_set_zero(o[2]); }
static inline int ff2_zero(ff_t o[3]) { return ( _ff_zero(o[0]) && _ff_zero(o[1]) ); }
static inline int ff3_zero(ff_t o[3]) { return ( _ff_zero(o[0]) && _ff_zero(o[1]) &&  _ff_zero(o[2]) ); }
static inline int ff2_one(ff_t o[3]) { return ( _ff_one(o[0]) && _ff_zero(o[1]) ); }
static inline int ff3_one(ff_t o[3]) { return ( _ff_one(o[0]) && _ff_zero(o[1]) &&  _ff_zero(o[2]) ); }
static inline int ff2_set(ff_t o[3], ff_t a[3]) { _ff_set(o[0],a[0]);  _ff_set(o[1],a[1]);  }
static inline int ff3_set(ff_t o[3], ff_t a[3]) { _ff_set(o[0],a[0]);  _ff_set(o[1],a[1]);  _ff_set(o[2],a[2]);  }
static inline int ff2_equal(ff_t o[3], ff_t a[3]) { return (_ff_equal(a[0],o[0])&&_ff_equal(a[1],o[1])); }
static inline int ff3_equal(ff_t o[3], ff_t a[3]) { return (_ff_equal(a[0],o[0])&&_ff_equal(a[1],o[1])&&_ff_equal(a[2],o[2])); }

// overlap ok
static inline void ff2_add(ff_t o[2], ff_t a[2], ff_t b[2]) { _ff_add(o[0],a[0],b[0]);   _ff_add(o[1],a[1],b[1]);  }
static inline void ff3_add(ff_t o[3], ff_t a[3], ff_t b[3]) { _ff_add(o[0],a[0],b[0]);   _ff_add(o[1],a[1],b[1]);   _ff_add(o[2],a[2],b[2]);  }
static inline void ff2_sub(ff_t o[2], ff_t a[2], ff_t b[2]) { _ff_sub(o[0],a[0],b[0]);   _ff_sub(o[1],a[1],b[1]);  }
static inline void ff3_sub(ff_t o[3], ff_t a[3], ff_t b[3]) { _ff_sub(o[0],a[0],b[0]);   _ff_sub(o[1],a[1],b[1]);   _ff_sub(o[2],a[2],b[2]);  }
static inline void ff2_neg(ff_t o[2], ff_t a[2]) { _ff_neg(o[0],a[0]); _ff_neg(o[1],a[1]); }
static inline void ff3_neg(ff_t o[3], ff_t a[3]) { _ff_neg(o[0],a[0]); _ff_neg(o[1],a[1]); _ff_neg(o[2],a[2]); }
static inline void ff2_negate(ff_t a[2]) { ff_negate(a[0]); ff_negate(a[1]); }
static inline void ff3_negate(ff_t a[3]) { ff_negate(a[0]); ff_negate(a[1]); ff_negate(a[2]); }

#define ff2_norm(o,a) 	ff2_norm_s(o,a,_ff_2g)
static inline void ff2_norm_s (ff_t o[1], ff_t a[2], ff_t s)
{
	register ff_t t1,t2;
	
	_ff_square(t1,a[1]);
	ff_square(o[0],a[0]);
	_ff_mult(t2,t1,s);
	_ff_subfrom(o[0],t2);
	// 3M+1A
}

#define ff2_trace(o,a)	ff_add(*(o),*(a),*(a))

static inline void ff2_scalar_mult (ff_t o[3], ff_t c, ff_t a[3]) { ff_mult(o[0],c,a[0]); ff_mult(o[1],c,a[1]); }

#define ff2_mult(o,a,b)	ff2_mult_s(o,a,b,_ff_2g)
static inline void ff2_mult_s (ff_t o[2], ff_t a[2], ff_t b[2], ff_t s)
{
	register ff_t t0, t1, t2;
	
	_ff_add(t0,a[0],a[1]);  _ff_add(t1,b[0],b[1]);  _ff_mult(t2,t0,t1);
	_ff_mult(t0, a[0], b[0]);  _ff_mult(t1, a[1], b[1]);
	_ff_subfrom(t2,t0);
	_ff_sub(o[1],t2,t1);
	_ff_mult(t2,t1,s);
	_ff_add(o[0],t0,t2);
	// 4M+5A  could use 5M+2A (essentially the same speed)
}

// multiplies (b[1]z+b[0])*(z+a) mod (z^2-s)
static inline void ff2_mult_zpa_s (ff_t o[2], ff_t b[2], ff_t a, ff_t s)
{
	register ff_t t0,t1,t2;
	
	_ff_mult(t1,a,b[1]);  _ff_mult(t0,a,b[0]);  _ff_mult(t2,b[1],s);
	_ff_add(o[1],b[0],t1);  _ff_add(o[0],t0,t2);
	/// 3M+2A
}

#define ff2_square(o,a)	ff2_square_s(o,a,_ff_2g)
static inline void ff2_square_s (ff_t o[2], ff_t a[2], ff_t s)
{
	register ff_t t0, t1, t2;
	
	_ff_square(t0, a[0]);  _ff_square(t1, a[1]);
	_ff_mult(t2,a[0],a[1]);
	_ff_add(o[1],t2,t2);
	_ff_mult(t2,t1,s);
	_ff_add(o[0],t0,t2);
	// 4M+2A
}

#define ff2_exp_ui(o,a,e)		ff2_exp_ui_s(o,a,e,_ff_2g)
void ff2_exp_ui_s (ff_t o[2], ff_t a[2], unsigned long e, ff_t s);	// computes in F_p[z]/(z^2-s)

// Cube roots in F_p^2 are only relevant when p=2mod3, since otherwise the 3-Sylow of F_p^2 is the 3-Sylow of F_p.
void _ff2_setup_cbrt(void);
static inline void ff2_setup_cbrt(void) { if ( ! _ff2_cbrt_setup ) _ff2_setup_cbrt(); }
int ff2_3Sylow_invcbrt (ff_t o[2], ff_t a[2]);

// Note ff3_setup_ncr() must be called before using any of the ff3 functions below
// ff3_poly_eval, ff3_sqrt, and ff3_exp_ui do this automatically

static inline void ff3_scalar_mult (ff_t o[3], ff_t c, ff_t a[3]) { ff_mult(o[0],c,a[0]); ff_mult(o[1],c,a[1]); ff_mult(o[2],c,a[2]); }
void ff3_mult (ff_t o[3], ff_t a[3], ff_t b[3]);			// too long to bother inlining
void ff3_square (ff_t o[3], ff_t a[3]);

// exponentiating by p (the Frobenius map) is very fast, potentially only 2M and at most 6M+4A, faster than ff3_mult by a lot
void ff3_exp_p (ff_t o[3], ff_t a[3]);

extern int _ff3_trace_z2;

// tr(a[0]+a[1]z+a[2]z^2 = 3a[0]+a[1]tr(z)+a[2]tr(z^2) = 3a[0]+a[2]tr(z^2) since tr(z)=0 for z^3-z-s=0 and z^3-s=0
static inline void ff3_trace (ff_t o[1], ff_t a[3])
{
	register ff_t t1,t2;
	
	_ff_add(t1,a[0],a[0]);
	_ff_add(o[0],t1,a[0]);
	if ( ! _ff3_trace_z2 ) return;
	_ff_add(t2,a[2],a[2]);			// we rely on the fact that tr(z^2) is either 0 or 2
	_ff_addto(o[0],t2);	
}

// norm(a)=a*a^p*a^(p^2), uses 16M+17A or 27M+23A, about 2 or 3 ff3_mults.
static inline void ff3_norm (ff_t o[1], ff_t a[3])
{
	ff_t x[3],y[3];
	register ff_t t1, t2;
	
	ff3_exp_p (x,a);				// x=a^p
	ff3_exp_p (y,x);				// y=a^(p^2)
	ff3_mult(x,x,y);					// x=a^(p^2+p)
	// we know the norm is in F_p, so we only compute the constant coefficient of a*x
	_ff_mult(t1,a[2],x[1]);			// t1=a2x1
	_ff_mult(t2,a[1],x[2]);			// t2=a1x2
	_ff_addto(t1,t2);				// t1=a2x1+a1x2
	_ff_mult(t2,t1,_ff_3g);			// t2=(a2x1+a1x2)s
	_ff_mult(t1,a[0],x[0]);			// t1=a0x0
	_ff_add(o[0],t1,t2);				// o=a0x0+(a2x2+a1x2)s		note that this works for both p=1mod3 or p=2mod3
}

void ff2_poly_eval (ff_t y[2], ff_t f[], int d, ff_t x[2]);
void ff3_poly_eval (ff_t y[3], ff_t f[], int d, ff_t x[3]);

int ff2_sqrt(ff_t o[2], ff_t a[2]);
int ff3_sqrt (ff_t o[3], ff_t a[3]);

int ff3_trsqrt_zpa_mod_rs (ff_t o[1], ff_t a, ff_t r, ff_t s);			// Computes tr(sqrt(z)) in F_p^3=F_p[z]/(z^3-rz-s).  This is a support function for factoring quartics.

void ff3_exp_ui (ff_t o[3], ff_t a[3], unsigned long e);
void ff3_exp_ui_rs (ff_t o[3], ff_t a[3], unsigned long e, ff_t r, ff_t s);
// inversion not currently supported/needed in degree 3

void ff3_zn_mod (ff_t o[3], unsigned long n, ff_t f[4]);			// only looks at f[0] and f[1], assumes f[2]=0 and f[3]=1

	
#endif
