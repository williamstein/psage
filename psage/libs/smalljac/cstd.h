#ifndef _CSTD_INCLUDE
#define _CSTD_INCLUDE

#include "gmp.h"

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

// General utility items that potentially could be used by any module

#define DEBUG_LEVEL		2
#define INFO_LEVEL			1
#define WARN_LEVEL		-1
#define ERROR_LEVEL		-2
#define OUTPUT_LEVEL		0

int dbg_level;

#define dbg_setlevel(i)		(dbg_level = (i))
#define dbg_printf			if ( dbg_level >= DEBUG_LEVEL ) gmp_printf
#define info_printf			if ( dbg_level >= INFO_LEVEL ) gmp_printf
#define warn_printf			if ( dbg_level >= WARN_LEVEL ) gmp_printf
#define err_printf			if ( dbg_level >= ERROR_LEVEL ) gmp_printf
#define out_printf			if ( dbg_level >= OUTPUT_LEVEL ) gmp_printf

#define delta_msecs(s,t)		(1000UL*(t-s)/CLOCKS_PER_SEC)
#define delta_nsecs(s,t)		(1000000000UL*(t-s)/CLOCKS_PER_SEC)			// assumes 64 bit UL
#define delta_secs(s,t)		((double)(t-s)/CLOCKS_PER_SEC)
#define delta_wall_msecs(w1,w2)		(1000*((w2)->tv_sec - (w1)->tv_sec) + ((w2)->tv_usec - (w1)->tv_usec)/1000)

// Centralized memory allocation to aid debugging.
// mem_alloc will never fail (it aborts on errors) and zero-fills all allocated memory.
void *mem_alloc (unsigned long bytes);
void mem_free (void *ptr);

#define my_mpz_init(z)		{ puts ("mpz_init"); mpz_init(z); }
#define my_mpz_clear(z)		{ puts ("mpz_clear"); mpz_clear(z); }

#endif
