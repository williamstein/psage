#ifndef _FFWRAPPER_INCLUDE_
#define _FFWRAPPER_INCLUDE_

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

// This header file contains declarations and wrappers for code that was
// written for a more general multi-precision version of ff.c and needs
// to deal with multi-precision and initialization issues.  It assume GMP is
// available

#define _ff_t_declare				ff_t
#define _ff_t_declare_static			static ff_t
#define _ff_t_declare_dynamic			ff_t
#if FF_WORDS == 1
#define _ff_t_declare_reg			register ff_t
#else
#define _ff_t_declare_reg			ff_t
#endif

#define _ff_init(x)
#define _ff_clear(x)

#define _ff_set_mpz(x,Z)				_ff_set_ui(x, mpz_fdiv_ui(Z,_ff_p))
#define _ff_get_mpz(Z,x)				(mpz_set_ui(Z,_ff_get_ui(x)),(Z))

#define FF_MPZ_WRAPS				16

// shared global temps - of course these are not thread safe and should probably have a home
int _ff_wrapper_initialized;
mpz_t __ff_mpz_p;
#define _ff_mpz_p					(mpz_set_ui(__ff_mpz_p,_ff_p), __ff_mpz_p)

mpz_t _ff_mpz_wraps[FF_MPZ_WRAPS];						// Wrappers used to convert to mpz for formatted output
#define _ff_wrap_mpz(x,i)				_ff_get_mpz(_ff_mpz_wraps[i],x)

#define _ff_random(x)				((x) = ui_randomm(_ff_p))
#define ff_random(x)				_ff_random(*(x))

// These macros are used for intermediate operations where reduction mod p may not necessary
// With single precision, its generally better to keep things reduced, so we don't attempt to optmize these
#define _ff_qneg(z,x)				_ff_neg(z,x);
#define _ff_qnegsum(z,x,y)			{ _ff_neg(z,x); _ff_subfrom(z,y); }
#define _ff_qsub(z,x,y)				_ff_sub(z,x,y);
#define _ff_qsubfrom(z,x)			_ff_subfrom(z,x);
#define _ff_qadd(z,x,y)				_ff_add(z,x,y);
#define _ff_qaddto(z,x)				_ff_addto(z,x);


// MUST BE CALLED AT LEAST ONCE
static inline void _ff_wrapper_init (void)
{
	int i;
	if ( _ff_wrapper_initialized ) return;
	mpz_init2(__ff_mpz_p, ULONG_BITS);
	for ( i = 0 ; i < FF_MPZ_WRAPS ; i++ ) mpz_init2 (_ff_mpz_wraps[i], ULONG_BITS);
	_ff_wrapper_initialized = 1;
}

static inline void ff_setup (mpz_t P)
	{  _ff_wrapper_init();  ff_setup_ui (mpz_get_ui(P)); }

#endif
