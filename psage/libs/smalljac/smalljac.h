#ifndef _SMALLJAC_INCLUDE_
#define _SMALLJAC_INCLUDE_

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

/*
	IMPORTANT NOTES
	
	The genus 1 algorithms now reside in hecurve1x.c - the general Jacobian algorithms are not used -
	but the same calling interface through smalljac_Lpolys is retained.  In addition to being more than
	twice as fast, version 3 corrects a bug in version 2 which could have resulted in some incorrect a_p
	values (this typically only impacted CM curves).
	
	Currently in genus > 1 computation of L-poly coefficients and group structure is
	limited to hyperelliptic curves of the form y^2=f(x) with f of odd degree (and also x^6+a), however
	pointcounting support (to compute a_1=-a_p) is provided for even degree f in genus 2 and 3
	(use the SMALLJAC_A1_ONLY flag for this purpose).

	The genus 3 code is still incomplete.  It will handle the typical case (minimal endomorphism ring)
	reliably, but for exceptional curves it may often fail for small values of p due to unresolved ambiguities.
	(when successful, the output is (barring a bug) provably correct).  Pointcounting is supported
	(within memory constraints - see SMALLJAC_INIT_COUNT_P) for all genus 3 curves.
*/

#define SMALLJAC_VERSION_STRING	"smalljac version 3.0"

//#define SMALLJAC_GENUS		1					// currently genus is hardwired at compile time - this will change

#define SMALLJAC_FF_64BITS		1					// set to 1 to use 64-bit finite field elements, default is 32 (this forces p < 2^31)

#define SMALLJAC_MAX_END		\
	(SMALLJAC_FF_64BITS ? (1UL<<63) : (1UL<<31))

#define SMALLJAC_MAX_GENUS	4					// genus 4 curves only supported for pointcounting
#define SMALLJAC_MAX_DEGREE	9					
#define SMALLJAC_RETRIES		40					// number of random elements to use to test group exponents 
#define SMALLJAC_FULL_P		(1UL<<33)			// determines when to do a full search of the Weil interval in genus 1(only)
#define SMALLJAC_MAX_COUNT_P	(1<<25)				// genus independent
#define SMALLJAC_BIG_COUNT_P	3600000				// experimentally determined on an AMD Athlon-64 4800+,  YMMV

#if SMALLJAC_GENUS == 1
#define SMALLJAC_MAX_P			(1UL<<40)			// this can easily be increased to 2^63 by increasing BSGS_STEPS (to approx 2^16) in hecurve1x.c
#define SMALLJAC_TINY_P			3					// must be at least 3
#define SMALLJAC_COUNT_P		(1UL<<10)			// no minimum requirement
#define SMALLJAC_PADIC_P		SMALLJAC_MAX_P		// never use p-adic
#define SMALLJAC_TABBITS		16					// no longer used in genus 1
#endif
#if SMALLJAC_GENUS == 2
#define SMALLJAC_MAX_P			(1UL<<31)			// group size needs to fit in a long
#define SMALLJAC_TINY_P			61					// should be >= 61.  p > 61 ensures a1 < p/2, only need p > 23 to ensures unique group order in Weil interval given a1.
#define SMALLJAC_COUNT_P		1500000				// crossover point determined on an AMD Athonl 2.5GHz - YMMV, should be at least 147 (p > 147 ensures unique group order in Weil interval)
#define SMALLJAC_PADIC_P		SMALLJAC_MAX_P		// p-adic is slower within the feasible range (crossover is >2^32)		
#define SMALLJAC_TABBITS		23
#endif
#if SMALLJAC_GENUS == 3
#define SMALLJAC_MAX_P			(1UL<<31)			// p^2+1 needs to fit in a long, its ok for the group size to be larger than 64 bits
#define SMALLJAC_TINY_P			47					// must be at least 47, needs to be 143 to ensure a1 < p/2.
#define SMALLJAC_COUNT_P		SMALLJAC_MAX_P		// always pointcount
#define SMALLJAC_PADIC_P		(1<<14)				// Should not be greater than 2^21 - group computation code assumes order < 2^63
#define SMALLJAC_TABBITS		19
#endif

#define SMALLJAC_INIT_COUNT_P	(1<<26)				// allocate extra space

// bit-masks for flags parameter
#define SMALLJAC_A1_ONLY    		0x1    				// only compute a1 coefficient
#define SMALLJAC_GOOD_ONLY  	0x2 			   		// only callback for good primes
#define SMALLJAC_SPLIT4			0x4					// only process half the primes mod 4, SMALLJAC_HIGH4 set => p mod 4 > 2, else p mod 4 <= 2
#define SMALLJAC_SPLIT8			0x8					// similarly for 8, 16, and 32 - allows primes to be conveniently partitioned among up to 16 jobs
#define SMALLJAC_SPLIT16		0x10
#define SMALLJAC_SPLIT32		0x20
#define SMALLJAC_HIGH4			0x40
#define SMALLJAC_HIGH8			0x80
#define SMALLJAC_HIGH16		0x100
#define SMALLJAC_HIGH32		0x200
#define SMALLJAC_FILTER			0x400				// flag to permit more general prime filtering - calls callback function with good=-1 to test whether p should be processed
#define SMALLJAC_GROUP		0x1000				// indicates group structure computation rather than L-poly coefficients
#define SMALLJAC_SGROUP		0x2000				// compute only the group structure actually generated by the black box, may be a subgroup.
												// only relevant when SMALLJAC_GROUP is set

#define SMALLJAC_SPLIT			0x3c					// all 2^k split flags
#define SMALLJAC_SPLIT_SHIFT		1
#define SMALLJAC_HIGH			0x3c0				// all 2^k high flags
#define SMALLJAC_HIGH_SHIFT		5

// error codes
#define SMALLJAC_INTERNAL_ERROR			-1			// general unspecified internal error - probably indicates a bug
#define SMALLJAC_INVALID_INTERVAL		-2			// end < start or end > SMALLJAC_MAX_END
#define SMALLJAC_PARSE_ERROR			-3			// couldn't parse curve string
#define SMALLJAC_UNSUPPORTED_CURVE	-4			// currently curve must have odd degree in [3,SMALLJAC_MAX_DEGREE]
#define SMALLJAC_SINGULAR_CURVE		-5			// curve is singular over Q
#define SMALLJAC_INVALID_PP				-6			// prime power (or prime) required
#define SMALLJAC_WRONG_GENUS			-7			// curve may be valid, but it doesn't have the compiled genus (this will go away eventually)
#define SMALLJAC_INVALID_FLAGS			-8			// invalid combination of flag bits set

// public interface for SAGE
typedef void *smalljac_Qcurve_t;						// private type - we don't want to spell out the details here - see smalljac_internal.h

smalljac_Qcurve_t smalljac_Qcurve_init (char *str, int *err);	// creates a curve over Q from a null-terminated curve specification, either [a1,a2,a3,a4,a6] (integers a1,a2,...,a6)
												// or a monic poly f(x) of odd degree with rational coefficients to define y^2 = f(x), err is an optional pointer to an error code
												
char *smalljac_Qcurve_str (smalljac_Qcurve_t c);			// returns a (readonly) pointer to curve specification used to create the curve
int smalljac_Qcurve_genus (smalljac_Qcurve_t c);			// returns the genus of the curve
void smalljac_Qcurve_clear (smalljac_Qcurve_t c);	  		// frees memory allocated for curve

long smalljac_Lpolys (smalljac_Qcurve_t c,				// curve over Q created via smalljac_Qcurve_init
                      unsigned long start,					// start and end specify closed search interval [start,end].  Set start=end=p to compute L_p(T)
                      unsigned long end,				
                      unsigned long flags,					// bitmask of options defined above
                      int (*callback)(smalljac_Qcurve_t c,			// pointer to callback function
                                      unsigned long p,				// prime in [start,end]
                                      int good,						// 1 if good reduction, 0 if bad
                                      long a[],						// n coefficients of L_p(T) with a[0] = a_1
                                      int n,							// either 1 or g for good p, 0 for bad p
                                      void *arg),					// forwarded arg from caller
                     void *arg);								// pass-through arg uninterpreted by smalljac

// smalljac_groups has the same interface as smalljac_Lpolys, except instead of an array a[] of L-poly coefficients,
// the array m[] contains the orders of the cyclic factors of J(C/F_p) ~ Z_m[0] x Z_m[1] x ... x Z_m[n-1]
// with m[0] | m[1] | ... | m[n-1] and 1 <= n <= 2g
static inline long smalljac_groups (smalljac_Qcurve_t curve, unsigned long start, unsigned long end, unsigned long flags,
					       int (*callback)(smalljac_Qcurve_t curve, unsigned long p, int good, long m[], int n, void *arg), void *arg)
    { return smalljac_Lpolys(curve, start, end, flags|SMALLJAC_GROUP, callback, arg); }

// One shot version for computing a single L-poly.  flags may be set to SMALLJAC_A1_ONLY
// Returns the number of coefficients computed (1 or g), 0 for bad reduction, or a negative error code
// Supports prime powers (but curve must still have rational coefficients)
int smalljac_Lpoly (long a[], char *curve, unsigned long q, unsigned long flags);

// One shot version for computing a single group structure.  Computes m[0],...,m[n-1] s.t.
// J(C/F_q) is isomorphic to Z_m[0] x ... x Z_m[r] with m[0] | m[1] | ... | m[n-1], where 1 <= n <= 2g.
// Returns n, the number of cyclic factors, 0 for bad reduction, or a negative error code
// Currently q must be prime and no flag values are used.
static inline int smalljac_group (long m[], char *curve, unsigned long q, unsigned long flags)
    { return smalljac_Lpoly (m, curve, q, flags|SMALLJAC_GROUP); }


#endif
