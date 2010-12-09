#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>
#include "gmp.h"
#include "mpzutil.h"
#include "cstd.h"

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

// This module is a grab-bag of stuff, not all of which has anything to do with
// large integer arithmetic (mpx).  This module should be split up and refined


#define MPZ_MAX_TINY_PRIME			251
#define MPZ_MAX_TINY_INTEGER		255
#define MPZ_TINY_PRIMES			54


#define _abs(a)					((a)<0?-(a):(a))

#define expr_valid_op(op)		((op) == 'E' || (op) == 'e' || (op) == '*' || (op) == '+' || (op) == '-')

static mpz_t mpz_util_primorial;
unsigned short mpz_util_primes[MPZ_SMALL_PRIMES+1];
unsigned short mpz_util_prime_index[MPZ_MAX_SMALL_INTEGER+1];
static mpz_t mpz_util_max_prime;
static mpz_t _mpz_temp;
static unsigned char mpz_small_factors[MPZ_MAX_SMALL_INTEGER+1];
unsigned short mpz_small_primorial_map[MPZ_SMALL_PRIMORIAL];		// i-th entry is 0 if gcd(i,MPZ_SMALL_PRIMORIAL)>1, otherwise it indicates that i is the k-th integer in F_P*, where 1 is the 1st.

int mpz_max_primorial_w;
unsigned long mpz_primorials[MPZ_MAX_PRIMORIAL_W+1];
unsigned long mpz_primorial_phis[MPZ_MAX_PRIMORIAL_W+1];
wheel_t mpz_small_wheels[MPZ_SMALL_PRIMORIAL_W+1];

static unsigned char _wheel_gaps0[1] = { 1 };		// The trivial wheel - every number is relatively prime to 1
static unsigned char _wheel_gaps1[2] = { 2 };		// The first wheel - odd numbers
static unsigned char _small_map[MPZ_MAX_SMALL_INTEGER+1];

#define MPZ_PP_TABSIZE						28

static struct {
	unsigned short i0, i1, i2;					// m0 contains product of primes with index in (i0,i1], m covers (i0,i2] and m = m0*m1
	unsigned short unsused;
	unsigned long b;							// any composite below b contains a prime factor with index <= i2
	mpz_t m0, m;
} _pptab[MPZ_PP_TABSIZE];

unsigned long _mpz_randomf (mpz_t o, mpz_t N, mpz_t factors[], mpz_ftree_t factor_trees[], unsigned long w);
int mpz_pollard_rho (mpz_t d, mpz_t n);

static gmp_randstate_t mpz_util_rands;
static int mpz_util_inited;

void mpz_util_init ()
{
	unsigned char *mp, *np;
	unsigned  i, j, k, n, p, maxp, w, gap, maxgap;
	unsigned long x;
	unsigned char *_small_map;
	
	if ( mpz_util_inited ) return;
	mpz_init (_mpz_temp);

	// i-th entry of small factors array holds a prime divisor of i for composite i, 0 o.w.
	for ( i = 2 ; i <= MPZ_MAX_TINY_PRIME ; ) {
		n = MPZ_MAX_SMALL_INTEGER / i;
		for ( j = 2 ; j <= n ; j++ ) {
			mpz_small_factors[i*j] = (unsigned char)i;
		}
		i++;
		while ( mpz_small_factors[i] ) i++;
	}
	
	// run through them in reverse order so that array holds the smallest prime divisor - more convenient for factoring
	for ( i = MPZ_MAX_TINY_PRIME ; i >= 2 ; ) {
		n = MPZ_MAX_SMALL_INTEGER / i;
		for ( j = 2 ; j <= n ; j++ ) {
			mpz_small_factors[i*j] = (unsigned char)i;
		}
		i--;
		while ( mpz_small_factors[i] ) i--;
	}
	
	mpz_primorials[0] = 1;
	mpz_primorial_phis[0] = 1;
	mpz_primorials[1] = 2;
	mpz_primorial_phis[1] = 1;
	mpz_small_wheels[0].n = 1;
	mpz_small_wheels[0].phi = 1;
	mpz_small_wheels[0].gaps = _wheel_gaps0;
	mpz_small_wheels[1].n = 2;
	mpz_small_wheels[1].phi = 1;
	mpz_small_wheels[1].gaps = _wheel_gaps1;
	_small_map = (unsigned char *) mem_alloc (MPZ_MAX_SMALL_INTEGER+1);
	_small_map[1] = 1;		// roll the first to prime the pump
	_small_map[3] = 1;
	mpz_util_primes[1] = 2;
	for ( i = 2 ; i <= MPZ_SMALL_PRIMORIAL_W ; i++ ) {
		for ( mp = _small_map+2 ; ! *mp ; mp++ );					// find first marked entry > 1 in previous map
		*mp;												// must be the next prime.  clear it.
		p = mp-_small_map;
		mpz_util_primes[i] = p;
		mpz_small_wheels[i].n = mpz_primorials[i] = mpz_primorials[i-1]*p;
		mpz_small_wheels[i].phi = mpz_primorial_phis[i] = mpz_primorial_phis[i-1]*(p-1);
		mp = _small_map+mpz_primorials[i-1] +1;
		for ( k = 1 ; k < p ; k++ ) {								// roll the previous wheel (p-1) more times
			for ( j = 0 ; j < mpz_primorial_phis[i-1] ; j++ ) { *mp = 1;  mp += mpz_small_wheels[i-1].gaps[j]; }
		}
		if ( mp != _small_map+mpz_primorials[i]+1 ) { err_printf ("mpz_util_init: primorial wheel alignment error\n");  exit (0); }
		*mp = 1;												// mark primorial+1 entry to bound last gap
		np = _small_map+mpz_primorials[i];
		for ( mp = _small_map+p ; mp < np ; mp += p ) *mp = 0;		// clear the multiples of p
		mpz_small_wheels[i].gaps = (unsigned char *) mem_alloc(mpz_primorial_phis[i]+4);	// make alloc big enough to avoid small alloc error
		maxgap = 0;
		for ( mp = _small_map+1, j = 0 ; j < mpz_primorial_phis[i] ; j++ ) {
			for ( np = mp+1 ; ! *np ; np++ );
			gap = np-mp;
			if ( gap > maxgap ) {
				maxgap = gap;
				if ( gap > UCHAR_MAX ) { err_printf ("gap %d too large to store in primorial wheel %d\n", gap, i);  exit (0); }
			}
			mpz_small_wheels[i].gaps[j] = gap;
			mp = np;
		}
		if ( np != _small_map+mpz_primorials[i]+1 ) { err_printf ("mpz_util_init: primorial wheel alignment error\n");  exit (0); }
		mpz_small_wheels[i].maxgap = maxgap;
		dbg_printf ("Wheel %d, maxgap %d\n", i, maxgap);
	}
	w = --i;
	// roll last wheel across rest of the small map
	mp = _small_map+mpz_primorials[w]+1;
	np = _small_map+MPZ_MAX_SMALL_INTEGER;
	while ( mp < np ) {
		for ( j = 0 ; j < mpz_primorial_phis[w] ; j++ ) {
			mp += mpz_small_wheels[w].gaps[j];
			if ( mp >= np ) break;
			*mp = 1;
		}
	}
	// now sieve for primes up to MPZ_MAX_SMALL_INTEGER
	maxp = (unsigned) floor(sqrt((double)MPZ_MAX_SMALL_INTEGER));
	mp = _small_map + mpz_util_primes[i];
	for (;;) {
		mp++;
		while ( ! *mp ) mp++;
		p = mp-_small_map;
		if ( p > maxp ) break;
		mpz_util_primes[++i] = p;
		for ( j = 0, k=1 ; ; j++ ) {
			if ( j == mpz_small_wheels[w].phi ) j = 0;
			k += mpz_small_wheels[w].gaps[j];
			if ( k*p > MPZ_MAX_SMALL_INTEGER ) break;
			_small_map[k*p] = 0;
		}
	}
	mp = _small_map+mpz_util_primes[i]+1;
	for (;;) {
		while ( ! *mp ) mp++;
		if ( mp >= np ) break;
		mpz_util_primes[++i] = mp-_small_map;
		mp++;
	}
	if ( i != MPZ_SMALL_PRIMES || mpz_util_primes[i] != MPZ_MAX_SMALL_PRIME ) { err_printf ("Error sieving small primes\n");  exit (0); }
	mem_free (_small_map);

	// index primes
	for ( i = 1 ; i <= MPZ_SMALL_PRIMES ; i++ ) mpz_util_prime_index[mpz_util_primes[i]] = i;
	// spread prime indexes to get pi(n) table
	for ( i = 0, j = 0 ; i <= MPZ_MAX_SMALL_INTEGER ; i++ ) {
		if ( mpz_util_prime_index[i] ) {
			j = mpz_util_prime_index[i];
		} else {
			mpz_util_prime_index[i] = j;
		}
	}
	
	// Compute various primorials and prime products
	mpz_init2 (mpz_util_primorial, 94027);					// the big Kahuna - the product of all small primes
	mpz_set_ui (mpz_util_primorial, 1);
	for ( i = 1 ; i <= MPZ_TINY_PRIMES ; i++ ) {
		if ( i > MPZ_SMALL_PRIMORIAL_W && i <= MPZ_MAX_PRIMORIAL_W ) {
			mpz_primorials[i] = mpz_primorials[i-1]*mpz_util_primes[i];
			mpz_primorial_phis[i] = mpz_primorial_phis[i-1]*(mpz_util_primes[i]-1);
		}
		mpz_mul_ui (mpz_util_primorial, mpz_util_primorial, mpz_util_primes[i]);
	}
	
	// early gcd gaps are smaller
	for ( j = 0 ; j < MPZ_PP_TABSIZE ; j++ ) {
		if ( j < 2 ) {
			n = 64;
		} else if ( j < 4 ) {
			n = 128;
		} else {
			n = 256;
		}
		mpz_init2 (_pptab[j].m0,n/2*16);   mpz_init2 (_pptab[j].m,n*16);
		mpz_set_ui (_pptab[j].m0, 1);  mpz_set_ui (_pptab[j].m, 1);
		_pptab[j].i0 = i-1;
		if ( i+256 < MPZ_SMALL_PRIMES ) {
			_pptab[j].i1 = i+n/2-1;			// note top of range is closed
			_pptab[j].i2 = i+n-1;
			_pptab[j].b = (unsigned long)mpz_util_primes[_pptab[j].i2+1]*(unsigned long)mpz_util_primes[_pptab[j].i2+1];
			for ( k = 0 ; k < n/8 ; k++ ) {
				x = (unsigned long)mpz_util_primes[i]*(unsigned long)mpz_util_primes[i+1]*(unsigned long)mpz_util_primes[i+2]*(unsigned long)mpz_util_primes[i+3];
				mpz_mul_ui (_pptab[j].m0, _pptab[j].m0, x);
				i += 4;
			}
			mpz_set (_pptab[j].m, _pptab[j].m0);
			for ( k = 0 ; k < n/8 ; k++ ) {
				x = (unsigned long)mpz_util_primes[i]*(unsigned long)mpz_util_primes[i+1]*(unsigned long)mpz_util_primes[i+2]*(unsigned long)mpz_util_primes[i+3];
				mpz_mul_ui (_pptab[j].m, _pptab[j].m, x);
				i += 4;
			}
		} else {
			_pptab[j].i2 = MPZ_SMALL_PRIMES;
			_pptab[j].i1 = _pptab[j].i0 + (_pptab[j].i2-_pptab[j].i0) / 2;
			_pptab[j].b = 0xFFFFFFFF;
			while ( i <= _pptab[j].i1 ) mpz_mul_ui (_pptab[j].m0, _pptab[j].m0, mpz_util_primes[i++]);
			mpz_set (_pptab[j].m, _pptab[j].m0);		
			while ( i <= _pptab[j].i2 ) mpz_mul_ui (_pptab[j].m, _pptab[j].m, mpz_util_primes[i++]);
		}
		mpz_mul (mpz_util_primorial, mpz_util_primorial, _pptab[j].m);
	}
	if ( i <= MPZ_SMALL_PRIMES ) { err_printf ("Initializion alignment problem computing prime products, only %d small primes used\n", i);  exit (0); }

	// roll small primorial wheel to create map of integers mod MPZ_MAX_SMALL_PRIMORIAL
	w = MPZ_SMALL_PRIMORIAL_W;
	for ( j = 0, k = 1; j < mpz_small_wheels[w].phi ; j++ ) {
		mpz_small_primorial_map[k] = j+1;
		k += mpz_small_wheels[w].gaps[j];
	}
	
	gmp_randinit_default (mpz_util_rands);
	x = (((unsigned long)gethostid())<<32) + getpid();
	x *= time(0);
//x=42;		//uncomment if you want to fix the seed for debugging
	gmp_randseed_ui (mpz_util_rands, x);
	mpz_util_inited = 1;
}


int ui_is_small_prime (unsigned long p)
{
	if ( ! p || p > MPZ_MAX_SMALL_PRIME ) return 0;
	if ( ! mpz_util_inited ) mpz_util_init();
	return ( mpz_util_primes[mpz_util_prime_index[p]] == p ? 1 : 0 );
}

int ui_is_prime (unsigned long p)
{
	static mpz_t P;
	static int init;
	
	if ( p <= MPZ_MAX_SMALL_PRIME ) return ui_is_small_prime(p);
	if ( ! init ) { mpz_init(P); init = 1; }
	mpz_set_ui(P,p);
	return mpz_probab_prime_p(P,5);
}

unsigned long ui_small_prime (unsigned long n)
{
	if ( ! n || n > MPZ_SMALL_PRIMES ) { err_printf ("Invalid request for small prime - %d > %d\n", n, MPZ_SMALL_PRIMES);  exit (0); }
	if ( ! mpz_util_inited ) mpz_util_init();
	return mpz_util_primes[n];
}


unsigned long ui_small_prime_index (unsigned long n)			// returns pi(n) for n <= MPZ_MAX_SMALL_INTEGER
{
	if ( n > MPZ_MAX_SMALL_PRIME ) { err_printf ("Invalid request for small prime index %d > %d\n", n, MPZ_MAX_SMALL_INTEGER);  exit (0); }
	if ( ! mpz_util_inited ) mpz_util_init();
	return mpz_util_prime_index[n];
}


unsigned long ui_primorial (int w)							// returns P_w = the w_th primorial or 0 if P_w > ULONG_MAX
{
	if ( w > MPZ_MAX_PRIMORIAL_W ) { err_printf ("Requested primorial P_%d is too large\n", w);  exit(0); }
	if ( ! mpz_util_inited ) mpz_util_init();
	return mpz_primorials[w];
}


unsigned long ui_primorial_phi (int w)							// returns the phi(P_w) or 0 if P_w > ULONG_MAX
{
	if ( w > MPZ_MAX_PRIMORIAL_W ) { err_printf ("Requested primorial P_%d is too large\n", w);  exit(0); }
	if ( ! mpz_util_inited ) mpz_util_init();
	return mpz_primorial_phis[w];
}

unsigned long ui_phi (unsigned long n)
{
	unsigned long p[MPZ_MAX_UI_PP_FACTORS];
	unsigned long h[MPZ_MAX_UI_PP_FACTORS];
	unsigned long m;
	int i, j, w;
	
	w = ui_factor(p,h,n);
	for ( m = 1, i = 0 ; i < w ; i++ ) {
		m *= (p[i]-1);
		for ( j = 1 ; j < h[i] ; j++ ) m *= p[i];
	}
	return m;
}

// returns explicit upper bound on pi(n) = (x/log x)(1+3/(2log x))  based on Shoup p.91 (from Rosser and Schoenfeld)
unsigned long ui_pi_bound (unsigned long n)
{
	double x,y;
	
	if ( n < 59 ) return 18;
	y = log((double)n);
	x = (double) n / y;
	x *= (1.0 + 3.0/(2.0*y));
	return ((unsigned long)ceil(x));
}

wheel_t *primorial_wheel_alloc (int w)
{
	wheel_t *wheel;
	char *map, *mp, *np;
	unsigned char *small_gaps, gap;
	unsigned long i, j, k, p, small_phi;
	
	if ( w <= MPZ_SMALL_PRIMORIAL_W ) return mpz_small_wheels+w;
	if ( w > MPZ_MAX_WHEEL_W ) {err_printf ("Requested wheel exceeds MPZ_MAX_WHEEL_W %d > %d\n", w, MPZ_MAX_WHEEL_W);  exit (0); }
	wheel = (wheel_t*)mem_alloc (sizeof(*wheel));
	wheel->n = mpz_primorials[w];
	wheel->phi = mpz_primorial_phis[w];
	wheel->gaps = (unsigned char *)mem_alloc (mpz_primorial_phis[w]);
	map = (char *)mem_alloc (wheel->n+1);
	np = map + wheel->n + 1;
	*np = 1;		// end marker to bound the last gap;
	small_phi = mpz_small_wheels[MPZ_SMALL_PRIMORIAL_W].phi;
	small_gaps = mpz_small_wheels[MPZ_SMALL_PRIMORIAL_W].gaps;
	for ( i = MPZ_SMALL_PRIMORIAL_W+1 ; i <= w ; i++ ) {
		p = mpz_util_primes[i];
		mp = map+p;
		j = 0;
		while ( mp < np ) {		// wheel the prime over the map
			*mp = 1;
			if ( j == small_phi ) j = 0;
			mp += p*small_gaps[j++];
		}
	}
	// compute the gaps by rolling the small wheel and skipping marked entries
	gap = 0;
	i = j = 0;
	mp = map+1;
	while ( mp < np ) {
		if ( j == small_phi  ) j = 0;
		mp += small_gaps[j];
		gap += small_gaps[j];
		if ( ! *mp ) {
			wheel->gaps[i++] = gap;
			if ( gap > wheel->maxgap ) wheel->maxgap = gap;
			gap = 0;
		}
		j++;
	}
	if ( mp != np) { err_printf ("wheel alignment error creating wheel for w = %d, %d != %d\n", w, mp-map, np-map);  exit (0); }
	if ( i != wheel->phi-1 ) { err_printf ("gap count error creating wheel for w = %d, %d != %d\n", w, i, wheel->phi);  exit (0); }
	wheel->gaps[i] = 2;	// the last gap is always 2
	return wheel;		
}


void primorial_wheel_free (wheel_t *wheel)
{
	if ( wheel->n <= MPZ_SMALL_PRIMORIAL ) return;
	mem_free (wheel->gaps);
	mem_free (wheel);
}

/*
	Fast prime enumeration for p <= L < 2^32 based on a wheeled sieve.
	The powers flag simply indicates that instead of returning p, return the largest power q=p^h bounded by L - enumeration is still ordered by p and only one power of p is returned.
*/
prime_enum_ctx_t *fast_prime_enum_start (unsigned start, unsigned L, int powers)
{
	prime_enum_ctx_t *ctx;
	unsigned char *mp, *np;
	unsigned long k, m;
	unsigned i, n, p, gap;

	if ( start > L ) { err_printf ("Invalid prime enumeration, start must be less than L\n");  exit (0); }
	if ( ! mpz_util_inited ) mpz_util_init();
	if ( L > MPZ_MAX_ENUM_PRIME ) { err_printf ("Attempted prime enumeration too large %d > %d = MPZ_MAX_ENUM_PRIME\n", L, MPZ_MAX_ENUM_PRIME);  exit (0); }
	ctx = (prime_enum_ctx_t *) mem_alloc (sizeof(*ctx));
	p = mpz_util_primes[MPZ_SMALL_PRIMORIAL_W+1];
//	if ( L < p*p ) L = p*p;									// make sure initial ctx->w is greater than MPZ_SMALL_PRIMORIAL_W
	ctx->L = L;
	ctx->powers = powers;									// indicates prime power enumeration
	ctx->h = ui_len(L);										// 2^h <= L <=2^{h+1}
	for ( i = 2 ; i <= ctx->h ; i++ ) {
		n = (unsigned) floor(pow((double)L,1.0/(double)i));		// n^i <= L < (n+1)^i
		ctx->r[i] = mpz_util_prime_index[n];						// r(i) = pi(L^{1/i})
		if ( ! powers ) break;
	}
	ctx->wheel = primorial_wheel_alloc (MPZ_SMALL_PRIMORIAL_W);
	ctx->w = ctx->r[2];										// p_w <= L^{1/2} < p_{w+1} so it suffices to sieve by the first w primes
	ctx->wv = (unsigned int*)mem_alloc ((ctx->w+1)*sizeof(*ctx->wv));			// the entry wv[i] stores the value of the i-th prime modulo MPZ_SMALL_PRIMORIAL = ctx->wheel->n
	ctx->wi = (unsigned int*)mem_alloc ((ctx->w+1)*sizeof(*ctx->wi));			// the entry wi[i] stores the index of the next wheel gap for the i-th prime - see below
	for ( i =1; i <= ctx->w ; i++ ) ctx->wv[i] = mpz_util_primes[i];	// these values aren't necessarily reduced mod MPZ_SMALL_PRIMORIAL, but that's ok, the code below doesn't assume this
	ctx->map = (unsigned char*)mem_alloc(ctx->wheel->n);
	np = ctx->map+ctx->wheel->n;
	/*
		The only relevant entries of our map are ones relatively prime to MPZ_SMALL_PRIMORIAL, since we enumerate the map via the wheel.
		Thus we only need to sieve multiples of primes whose index lies in [MPZ_SMALL_PRIMORIAL_W+1,ctx->w] and we only need to worry
	        about those multiples which are relatively prime to MPZ_SMALL_PRIMORIAL.  Thus for each prime, we effectively use the same wheel
	        to enumerate these multiples.  The value wi[i] holds the index into the wheel for the i-th prime.
	
		If the starting point is far from 0, we can skip ahead by starting at the first multiple of MPZ_SMALL_PRIMORIAL below start, call it n.
	        For simplicity, we avoid special code for the first ctx->w primes by always forcing enumeration of the first w primes
		For each prime p <= ctx->w we need to compute the least multiple kp > start where k is relatively prime to MPZ_SMALL_PRIMORIAL
	        and then set wv[i] = kp - n and set wi[i] = index of gap between k and next integer relatively prime to MPZ_SMALL_PRIMORIAL in wheel.
	
		To keep this simple, we set m to the first multiple of p*MPZ_SMALL_PRIMORIAL below start, set wv[i] = p, wi[i] = 0 and roll forward from there.
	*/
	if ( start > 2 ) {
		// This is a bit awkward but removes the need for any conditional code in fast_prime_enum to handle startup.
		// We want to enumerate just up to the greatest prime below start.  To facilitate this we backup start as required
		n = ui_len(start);
		for ( k = 2 ;; k += n ) {
			if ( k > start ) k = start;
			mpz_set_ui (_mpz_temp, start-k);
			mpz_nextprime (_mpz_temp, _mpz_temp);
			if ( mpz_cmp_ui(_mpz_temp,start) < 0 ) break;
		}
		do {
			m = mpz_get_ui(_mpz_temp);
			mpz_nextprime (_mpz_temp, _mpz_temp);
		} while ( mpz_cmp_ui(_mpz_temp,start) < 0 );
		start = m;
	} else {
		start = 0;
	}
	if ( start > MPZ_SMALL_PRIMORIAL && start > mpz_util_primes[ctx->w] ) {
		n = (start/MPZ_SMALL_PRIMORIAL) * MPZ_SMALL_PRIMORIAL;
		for ( i = MPZ_SMALL_PRIMORIAL_W+1 ; i <= ctx->w ; i++ ) {
			p = mpz_util_primes[i];
			ctx->wv[i] = p;
			m = _ui_ceil_ratio(n,p)*(unsigned long)p;								// m is the least multiple of p >= n - could be > 2^32
			m /= p;
			k = m % MPZ_SMALL_PRIMORIAL;
			while ( ! mpz_small_primorial_map[k] ) k++;								// find least k >= m%MPZ_SMALL_PRIMORIAL relatively prime to MPZ_SMALL_PRIMORIAL
			m = (unsigned long)p*((m/MPZ_SMALL_PRIMORIAL)*MPZ_SMALL_PRIMORIAL+k) - n;
			ctx->wv[i] = (unsigned) m;
			ctx->wi[i] = mpz_small_primorial_map[k]-1;								// gap index i gives gap between (i+1)st and (i+2)nd integers relatively prime to MPZ_SMALL_PRIMORIAL
		}
		ctx->pi = MPZ_SMALL_PRIMES;		// just need to make it bigger than ctx->w and MPZ_SMALL_PRIMORIAL_W
		ctx->p = n-1;
	} else {
		ctx->pi = 1;
		ctx->p = 0;
	}
	ctx->j = ctx->wheel->phi-1;
/*
	for ( i = MPZ_SMALL_PRIMORIAL_W+1 ; i <= ctx->w ; i++ ) {
		p = mpz_util_primes[i];
		mp = ctx->map+ctx->wv[i];
		while ( mp < np ) {
			*mp = 1;
			gap = p*ctx->wheel->gaps[ctx->wi[i]++];
			if ( ctx->wi[i] == ctx->wheel->phi ) ctx->wi[i] = 0;
			ctx->wv[i] += gap;
			mp += gap;
		}
		ctx->wv[i] -= ctx->wheel->n;
	}
	ctx->mp = ctx->map + 1;
*/
	if ( start ) {
		while ( (p=fast_prime_enum(ctx)) < start );
		if ( p != start ) { err_printf ("Unable to find start prime, p = %u, start = %u, L = %u\n", p, start, L);  exit (0); }
	}
	return ctx;
}


unsigned fast_prime_enum (prime_enum_ctx_t *ctx)
{
	register unsigned gap;

	// The first ctx->w primes need to be enumerated seperately, since they have all been sieved out of the map
	// Note that if ctx->w > MPZ_SMALL_PRIMORIAL this means the first map is effectively empty, but that's ok.
	// Prime powers are only relevant for the first ctx->w primes since ctx->w was chosen to ensure pi(ctx->w)^2 > L
	if ( ctx->pi <= ctx->w || ctx->pi <= MPZ_SMALL_PRIMORIAL_W ) {
		register unsigned i, p, q;
		q = mpz_util_primes[ctx->pi++];
		if ( q > ctx->L ) return 0;
		if ( ctx->powers ) {
			while ( ctx->h && ctx->r[ctx->h] < ctx->pi ) ctx->h--;
			for ( i = 1, p=q ; i < ctx->h ; i++ ) q *= p;
		}
		return q;
	}
	for (;;) {
		if ( ctx->j == ctx->wheel->phi-1 ) {
			register unsigned char *mp, *np;
			register unsigned i, p;
			memset (ctx->map, 0, ctx->wheel->n);
			np = ctx->map + ctx->wheel->n;
			for ( i = MPZ_SMALL_PRIMORIAL_W+1 ; i <= ctx->w ; i++ ) {
				p = mpz_util_primes[i];
				mp = ctx->map+ctx->wv[i];
				while ( mp < np ) {
					*mp = 1;
					gap = p*ctx->wheel->gaps[ctx->wi[i]++];
					if ( ctx->wi[i] == ctx->wheel->phi ) ctx->wi[i] = 0;
					ctx->wv[i] += gap;
					mp += gap;
				}
				ctx->wv[i] -= ctx->wheel->n;
			}
			ctx->mp = ctx->map+1;
			if ( ctx->p ) {
				ctx->p += 2;									// last gap is always 2
				ctx->j = 0;
			} else {
				ctx->p = 1+ctx->wheel->gaps[0]; 					// special case - skip 1
				ctx->mp += ctx->wheel->gaps[0]; 
				ctx->j = 1;
			}
		} else {
			gap = ctx->wheel->gaps[ctx->j++];
			ctx->mp += gap;
			ctx->p += gap;
		}
		if ( ctx->p > ctx->L ) return 0;
		if ( ! *ctx->mp ) return ctx->p;
	}
}


void fast_prime_enum_end (prime_enum_ctx_t *ctx)
{
	mem_free (ctx->map);
	mem_free (ctx->wv);
	mem_free (ctx->wi);
	primorial_wheel_free (ctx->wheel);
	mem_free (ctx);
}


// returns the least p for which the p-adic valuation of a and b differ
unsigned long ui_pdiff (unsigned long a, unsigned long b)
{
	unsigned long d, p;
	unsigned long q[MPZ_MAX_UI_PP_FACTORS];
	unsigned long h[MPZ_MAX_UI_PP_FACTORS];
	int w;
	int i;
	
	mpz_util_init();
	if ( a == b ) return 0;
	d = ui_gcd(a,b);
	a /= d;
	b /= d;
	// check small primes before factoring
	for ( i = 1 ; i < 25 ; i++ ) {
		if ( a%mpz_util_primes[i] ) {
			if ( ! (b%mpz_util_primes[i]) ) break;
		} else {
			if ( b%mpz_util_primes[i] ) break;
		}
	}
	if ( i < 25 ) return mpz_util_primes[i];
printf ("Factoring %lu and %lu\nb", a, b);
	w = ui_factor (q, h, a);
	if ( w ) p = q[0]; else p = 0;
	w = ui_factor (q, h, b);
	if ( w && p > q[0] ) p = q[0];
printf ("pdiff = %d\n", p);
	return p;
}


// note that o and a can overlap
void mpz_parallel_invert (mpz_t o[], mpz_t a[], unsigned n, mpz_t p)
{
	static int init;
	static mpz_t c[MPZ_MAX_INVERTS];
	static mpz_t u, v;
	unsigned i;
	
	if ( ! init ) {
		for ( i = 0 ; i < MPZ_MAX_INVERTS ; i++ ) mpz_init (c[i]);
		mpz_init (u);  mpz_init (v);
		init = 1;
	}
	if ( ! n ) return;
	if ( n > MPZ_MAX_INVERTS ) { err_printf ("Exceeded MPZ_MAX_INVERTS, %d > %d\n", n, MPZ_MAX_INVERTS);  exit (0); }
	mpz_set (c[0], a[0]);
	for ( i = 1 ; i < n ; i++ ) mpz_mulm (c[i], c[i-1], a[i], p);
	if ( ! mpz_invert (u, c[n-1], p) ) { err_printf ("Invert failed in mpz_parallel_invert!\n"); }
	for ( i = n-1 ; i > 0 ; i-- ) {
		mpz_mulm (v, c[i-1], u, p);
		mpz_mulm (u, u, a[i], p);
		mpz_set (o[i], v);
	}
	mpz_set (o[0], u);
}


void mpz_mul_set (mpz_t o, mpz_t *a, unsigned long n)
{
	unsigned long i, j;
	
	if ( ! n ) { mpz_set_ui(o, 1); return; }
	while ( n >= 2 ) {
		for ( i = 0 ; i < (n+1)/2 ; i++ ) {
			if ( 2*i+1 < n ) {
				mpz_mul (a[i], a[2*i], a[2*i+1]);
			} else {
				mpz_set (a[i], a[2*i]);
			}
		}
		n = i;
	}
	if ( n == 2 ) {
		mpz_mul (o, a[0], a[1]);
	} else {
		mpz_set (o, a[0]);
	}
}

#define MPZ_PM_CACHE_SIZE		100

static struct {
	mpz_t p;
	mpz_t maxp;
	mpz_t endp;
	mpz_t o;
	unsigned long maxbits;
	int n;
} _mpz_pm_cache[MPZ_PM_CACHE_SIZE];
unsigned long _mpz_pm_cache_count;


#define MPZ_TIER_SIZE		256

// multiplies o by the product of all primes in (p,maxp] up to maxbits and updates p to last prime used or > maxp if all used.  return number of primes multiplied.
int mpz_prime_mult (mpz_t o, mpz_t p, mpz_t maxp, unsigned long maxbits)
{
	static int init, warn;
	static mpz_t t1[MPZ_TIER_SIZE];
	static mpz_t t2[MPZ_TIER_SIZE];
	static mpz_t t3[MPZ_TIER_SIZE];
	unsigned long bits;
	register int i, i1, i2, i3, n;
	
	if ( ! init ) {
		for ( i = 0 ; i < MPZ_TIER_SIZE ; i++ ) { mpz_init (t1[i]); mpz_init (t2[i]); mpz_init (t3[i]); }
		for ( i = 0 ; i < MPZ_PM_CACHE_SIZE ; i++ ) {
			mpz_init (_mpz_pm_cache[i].p);
			mpz_init (_mpz_pm_cache[i].endp);
			mpz_init (_mpz_pm_cache[i].maxp);
			mpz_init (_mpz_pm_cache[i].o);
		}
		init = 1;
	}
	for ( i = 0 ; i < _mpz_pm_cache_count ; i++ ) {
		if ( mpz_cmp (p,_mpz_pm_cache[i].p) == 0 && mpz_cmp (maxp, _mpz_pm_cache[i].maxp) == 0 &&
		     maxbits == _mpz_pm_cache[i].maxbits ) {
			mpz_set (p, _mpz_pm_cache[i].endp);
			mpz_set (o, _mpz_pm_cache[i].o);
			return _mpz_pm_cache[i].n;
		}
	}
	if ( _mpz_pm_cache_count < MPZ_PM_CACHE_SIZE ) { 
		i = _mpz_pm_cache_count;
		mpz_set (_mpz_pm_cache[i].p, p);
		mpz_set (_mpz_pm_cache[i].maxp, maxp);
		_mpz_pm_cache[i].maxbits = maxbits;
	} else {
		if ( ! warn ) { err_printf ("Prime product cache full in mpz_prime_mult...\n");  warn = 1; }
	}
	i1 = i2 = i3 = 0;
	bits = mpz_sizeinbase (o, 2);
	for ( n = 0 ; bits < maxbits ; n++ ) {
		mpz_nextprime (p, p);
		if ( mpz_cmp (p, maxp) > 0 ) break;
		if ( i3 == MPZ_TIER_SIZE ) {
			if ( i2 == MPZ_TIER_SIZE ) {
				if ( i1 >= MPZ_TIER_SIZE-2 ) break;
				mpz_mul_set (t1[i1++], t2, i2);
				i2 = 0;
			}
			mpz_mul_set (t2[i2++], t3, i3);
			i3 = 0;
		}
		mpz_set (t3[i3++], p);
		bits += mpz_sizeinbase (p, 2);
	}
	if ( i2 == MPZ_TIER_SIZE ) {
		mpz_mul_set (t1[i1++], t2, i2);
		i2 = 0;
	}
	mpz_mul_set (t2[i2++], t3, i3);
	mpz_mul_set (t1[i1++], t2, i2);
	mpz_mul_set (t3[0], t1, i1);
	mpz_mul (o, o, t3[0]);
	if ( _mpz_pm_cache_count < MPZ_PM_CACHE_SIZE ) {
		i = _mpz_pm_cache_count++;
		mpz_set (_mpz_pm_cache[i].endp, p);
		mpz_set (_mpz_pm_cache[i].o, o);
		_mpz_pm_cache[i].n = n;
	}
	return n;
}

// computes the product of all prime powers <= L
// this function is designed to be called once - it dynamically allocates and frees all memory
void mpz_power_primorial (mpz_t o, unsigned long L)	
{
	mpz_t t1[MPZ_TIER_SIZE];
	mpz_t t2[MPZ_TIER_SIZE];
	mpz_t t3[MPZ_TIER_SIZE];
	mpz_t p, q, t;
	unsigned long roots[33];
	register int i, i1, i2, i3, j1, j2, n;
	
	n = ui_len(L);
	if ( n > 32 ) { err_printf ("primorial %u too large to compute - be reasonable.\n");  exit (0); }
	for ( i = n ; i >= 2 ; i-- ) roots[i] = (unsigned long) floor(pow((double)L,1.0/(double)i));
	roots[1] = L;
	roots[0] = -1;
	
	// just init tier 3 initially, init tier 1 and 2 variables as needed
	for ( i = 0 ; i < MPZ_TIER_SIZE ; i++ ) mpz_init2 (t3[i], MPZ_TIER_SIZE*n);
	j1 = j2 = 0;

	mpz_init (p);  mpz_init (q);  mpz_init (t);  mpz_set_ui (p, 1);
	i1 = i2 = i3 = 0;
	i = n;
	for ( ;; ) {
		mpz_nextprime (p, p);
		while ( i && mpz_cmp_ui (p,roots[i]) > 0 ) i--;
		if ( ! i ) break;
		mpz_pow_ui (q, p, i);
		if ( i3 == MPZ_TIER_SIZE ) {
			if ( i2 == MPZ_TIER_SIZE ) {
				if ( i1 >= MPZ_TIER_SIZE-2 ) { err_printf ("Insufficient memory in mpz_power_primorial on input %u - increase MPZ_TIER_SIZE or add a tier\n", L);  exit (0); }
				if ( i1 == j1 ) mpz_init (t1[j1++]);
				mpz_mul_set (t1[i1++], t2, i2);
				i2 = 0;
			}
			if ( i2 == j2 ) mpz_init (t2[j2++]);
			mpz_mul_set (t2[i2++], t3, i3);
			i3 = 0;
		}
		mpz_set (t3[i3++], q);
	}
	if ( i2 == MPZ_TIER_SIZE ) {
		if ( i1 == j1 ) mpz_init (t1[j1++]);
		mpz_mul_set (t1[i1++], t2, i2);
		i2 = 0;
	}
	// free memory as we go to keep peak usage down
	if ( i2 == j2 ) mpz_init (t2[j2++]);
	mpz_mul_set (t2[i2++], t3, i3);
	for ( i = 0 ; i < MPZ_TIER_SIZE ; i++ ) mpz_clear (t3[i]);
	if ( i1 == j1 ) mpz_init (t1[j1++]);
	mpz_mul_set (t1[i1++], t2, i2);
	for ( i = 0 ; i < j2 ; i++ ) mpz_clear (t2[i]);
	mpz_mul_set (o, t1, i1);
	for ( i = 0 ; i < j1 ; i++ ) mpz_clear (t1[i]);
	mpz_clear (p);  mpz_clear (q);  mpz_clear (t);
	return;
}


int mpz_remove_small_primes (mpz_t o, mpz_t n, unsigned long exps[], unsigned long maxprimes)
{
	static int init;
	static mpz_t d, x, t;
	int i, j, w;
	
	if ( ! init ) { mpz_util_init();  mpz_init (d);  mpz_init (x);  mpz_init (t);  init = 1; }
	if ( maxprimes > MPZ_SMALL_PRIMES ) { err_printf ("maxprimes value %d exceeded MAX_SMALL_PRIMES = %d in mpz_remove_small_primes\n", maxprimes, MPZ_SMALL_PRIMES); exit (0); }
	mpz_gcd (d, n, mpz_util_primorial);
	mpz_set (o, n);
	exps[0] = 0;
	w = 0;
	for ( i = 1 ; i <= maxprimes ; i++ ) {
		if ( mpz_divisible_ui_p (d, mpz_util_primes[i]) ) {
			mpz_set_ui (x, mpz_util_primes[i]);
			exps[i] = mpz_remove (o, o, x);
			w++;
		} else {
			exps[i] = 0;
		}
	}
	return w;
}

// Let p be the largest prime factor of n.  If n/p <= L, return n/p, otherwise return 0
unsigned long mpz_nearprime (mpz_t n, unsigned long L)
{
	static int init;
	static mpz_t d, x, t;
	int i;
	
	if ( ! init ) { mpz_util_init();  mpz_init (d);  mpz_init (x);  mpz_init (t);  init = 1; }
	if ( L > MPZ_MAX_SMALL_INTEGER ) { err_printf ("cofactor bound too large in mpz_nearprime\n"); exit (0); }
	mpz_gcd (d, n, mpz_util_primorial);
	if ( mpz_cmp_ui (d, L) > 0 ) return 0;
	mpz_set (t, n);
	for ( i = 1 ; i <= MPZ_SMALL_PRIMES ; i++ ) {
		if ( mpz_divisible_ui_p (d, mpz_util_primes[i]) ) {
			mpz_set_ui (x, mpz_util_primes[i]);
			mpz_remove (t, t, x);
			mpz_remove (d, d, x);
			if ( mpz_cmp_ui(d,1) == 0 ) break;
		}
	}
	mpz_divexact (x, n, t);
	if ( mpz_cmp_ui (x, L) > 0 ) return 0;
	if ( ! mpz_probab_prime_p (t, 10) ) return 0;
	return mpz_get_ui (x);
}


int mpz_remove_small_squares (mpz_t o, mpz_t n)
{
	static int init;
	static mpz_t d, x, m;
	int i, j, w;
	
	if ( ! init ) { mpz_util_init();  mpz_init (d);  mpz_init (x);   mpz_init (m);  init = 1; }
	mpz_set (o, n);
	if ( mpz_cmp_ui (n, 1) <= 0 ) return 0;
	mpz_gcd (d, n, mpz_util_primorial);
	mpz_set (m, n);
	for ( i = 1 ; i <= MPZ_SMALL_PRIMES ; i++ ) {
		if ( mpz_divisible_ui_p (d, mpz_util_primes[i]) ) {
			mpz_set_ui (x, mpz_util_primes[i]);
			j = mpz_remove (m, m, x);
			if ( j%2 ) mpz_mul (m, m, x);
		}
	}
	mpz_set (o, m);
	return 1;
}


int mpz_print_factors (mpz_t N)
{
	static int init;
	static mpz_t P;
	unsigned long p[MPZ_MAX_SMALL_FACTORS];
	unsigned long h[MPZ_MAX_SMALL_FACTORS];
	int i, w;
	
	if ( ! init ) { mpz_init (P);  init = 1; }
	w = mpz_factor_small (p, h, P, N, MPZ_MAX_SMALL_FACTORS, 128);
	if ( mpz_cmp_ui (P,1) > 0 ) out_printf ("%Zd", P);
	if ( w ) out_printf (" ");
	for ( i = 0 ; i < w ; i++ ) {
		if ( h[i] > 1 ) {
			out_printf ("%lu^%lu ", p[i], h[i]);
		} else if ( h[i] == 1 ) {
			out_printf ("%lu", p[i]);			
		}
	}
	puts ("");
}


#define _remove(n,d,c)			while ( !((n)%(d)) ) { (n)/=(d); (c)++; }
#define _tdiv(n,d,h,p,w)			{ _remove(n,d,h[w]); if ( (h)[(w)] ) { (p)[(w)++] = (d); (h)[(w)] = 0; } }

int _mod64res[64] = { 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, };
int _mod63res[63] = { 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 };
int _mod65res[65] = { 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, };
int _mod11res[11] = { 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0 };

static inline int _ui_sqrt(unsigned long n)
{
	register unsigned long x;
	
	if ( ! _mod64res[n&0x3F] ) return 0;
	if ( ! _mod63res[n%63] ) return 0;
	if ( ! _mod65res[n%65] ) return 0;
	if ( ! _mod11res[n%11] ) return 0;
	x = (unsigned long)(sqrt(n)+0.01);		// avoid potential rounding problems
	return (x*x == n ? x : 0);
}

int ui_sqrt(long n)
  { register int k; if (n<0) return -1;  if ( !n) return 0;  k = _ui_sqrt(n); if ( k ) return k; else return -1; }


int ui_factor (unsigned long p[MPZ_MAX_UI_PP_FACTORS], unsigned long h[MPZ_MAX_UI_PP_FACTORS], unsigned long n)
{
	static int init;
	static mpz_t N, D, D1;
	register unsigned long d, d0, d1;
	register int i, k, w;
	
	if ( ! init ) { mpz_init (N);  mpz_init (D);  mpz_init (D1); mpz_util_init(); init = 1; }
	if ( n == 0 ) { err_printf ("attempt to factor 0\n");  exit (0); }
	w = 0;
	h[w] = 0;
	while ( ! (n&0x1) ) { n >>= 1; h[w]++; }
	if ( h[w] ) { p[w++] = 2;  h[w] = 0; }
	if ( n == 1 ) return w;

	// factor small integers using factor table - but need to reverse the order of the factors
	if ( n <= MPZ_MAX_SMALL_INTEGER ) {
		register unsigned char tiny;
		unsigned long p2[16];
		
		k = 0;
		while ( (tiny = mpz_small_factors[n]) ) {
			if ( w && p[w-1] == tiny ) { h[w-1]++; } else { p[w] = tiny; h[w++] = 1; }
			n /= tiny;
		}
		if ( w && p[w-1] == n ) { h[w-1]++; } else { p[w] = n; h[w++] = 1; }
		return w;
	}
	// Clear out all tiny primes via trial division.  Yes, it is worth hardwiring all this - it's much faster.
	_tdiv(n,3,h,p,w); _tdiv(n,5,h,p,w); _tdiv(n,7,h,p,w); _tdiv(n,11,h,p,w); _tdiv(n,13,h,p,w); _tdiv(n,17,h,p,w); _tdiv(n,19,h,p,w); _tdiv(n,23,h,p,w); _tdiv(n,29,h,p,w);
	_tdiv(n,31,h,p,w); _tdiv(n,37,h,p,w); _tdiv(n,41,h,p,w); _tdiv(n,43,h,p,w); _tdiv(n,47,h,p,w); _tdiv(n,53,h,p,w); _tdiv(n,59,h,p,w); _tdiv(n,61,h,p,w); _tdiv(n,67,h,p,w); _tdiv(n,71,h,p,w);
	_tdiv(n,73,h,p,w); _tdiv(n,79,h,p,w); _tdiv(n,83,h,p,w); _tdiv(n,89,h,p,w); _tdiv(n,97,h,p,w); _tdiv(n,101,h,p,w); _tdiv(n,103,h,p,w); _tdiv(n,107,h,p,w); _tdiv(n,109,h,p,w); _tdiv(n,113,h,p,w);
	_tdiv(n,127,h,p,w); _tdiv(n,131,h,p,w); _tdiv(n,137,h,p,w); _tdiv(n,139,h,p,w); _tdiv(n,149,h,p,w); _tdiv(n,151,h,p,w); _tdiv(n,157,h,p,w); _tdiv(n,163,h,p,w); _tdiv(n,167,h,p,w); _tdiv(n,173,h,p,w);
	_tdiv(n,179,h,p,w); _tdiv(n,181,h,p,w); _tdiv(n,191,h,p,w); _tdiv(n,193,h,p,w); _tdiv(n,197,h,p,w); _tdiv(n,199,h,p,w); _tdiv(n,211,h,p,w); _tdiv(n,223,h,p,w); _tdiv(n,227,h,p,w); _tdiv(n,229,h,p,w);
	_tdiv(n,233,h,p,w); _tdiv(n,239,h,p,w); _tdiv(n,241,h,p,w); _tdiv(n,251,h,p,w);
	
	if ( n == 1 ) return w;
	if ( n < MPZ_MAX_SMALL_INTEGER ) {
		p[w] = n;
		h[w++] = 1;
		return w;
	}

	// at this point n > MPZ_MAX_SMALL_INTEGER and contains no primes <= MPZ_MAX_TINY_PRIME
	// time to start whacking it with some gcd's
	mpz_set_ui (N, n);
	for ( i = 0 ; i < MPZ_PP_TABSIZE ; i++ ) {
		mpz_gcd (D, N,_pptab[i].m);
		if ( mpz_cmp_ui(D,1) != 0 ) {
			d = mpz_get_ui(D);
			if ( d <= mpz_util_primes[_pptab[i].i2]  ) {		// sliced of just one prime - this is what we want
				p[w] = d;
				h[w] = 0;
				do { n /= d;  h[w]++; } while ( ! (n%d) );	// need to remove all occurrences of d
				w++;
				mpz_set_ui (N, n);
			} else {
//printf ("multiple primes in gcd span\n");
				mpz_gcd (D1, D, _pptab[i].m0);
				if ( mpz_cmp_ui(D1,1) ) {
					d1 = mpz_get_ui (D1);
					if ( d1 <= mpz_util_primes[_pptab[i].i1] ) {
						p[w] = d1;
						h[w] = 0;
						do { n /= d1;  h[w]++; } while ( ! (n%d1) );	// need to remove all occurrences of d1
						w++;
					} else {
//printf ("trial dividing low i=%d, d1=%lu, (%u,%u)\n", i, d1, mpz_util_primes[_pptab[i].i0],mpz_util_primes[_pptab[i].i1]);
						h[w] = 0;
						for ( k = _pptab[i].i0+1 ; k <= _pptab[i].i1 ; k++ ) _tdiv (n, mpz_util_primes[k], h, p, w);
					}
					d1 = d/d1;
				} else {
					d1 = d;
				}
				if ( d1 > 1 ) {
					if ( d1 < mpz_util_primes[_pptab[i].i2]) {
						p[w] = d1;
						h[w] = 0;
						do { n /= d1;  h[w]++; } while ( ! (n%d1) );	// need to remove all occurrences of d1
						w++;
					} else {
//printf ("trial dividing high i=%d, d1=%lu, (%u,%u)\n", i, d1, mpz_util_primes[_pptab[i].i1],mpz_util_primes[_pptab[i].i2]);
						h[w]  = 0;
						for ( k = _pptab[i].i1+1 ; k <= _pptab[i].i2 ; k++ ) _tdiv (n, mpz_util_primes[k], h, p, w);
					}
				}
			}
			if ( n == 1 ) return w;
		}
		if ( n < _pptab[i].b ) {
//printf ("added remainder n=%d\n", n);
			p[w] = n; h[w++] = 1;
			return w;
		}
		mpz_set_ui (N, n);
	}

	// If the original n's second largest prime was below MPZ_SMALL_PRIME, what remains is prime - this is likely for n in the 32-48 bit range, so check here.
	if  ( mpz_probab_prime_p (N, 10) ) {
		p[w] = n; h[w++] = 1;
		return w;
	}
	
	// resort to bigger guns - this can't happen unless N is a composite with more than 32 bits and no prime factors below MPZ_MAX_SMALL_PRIME
	// we could go direct to a single-precision version of pollard rho here, but we don't expect this to happen much anyway
	w += mpz_factor_small(p+w,h+w,D,N,MPZ_MAX_UI_PP_FACTORS-w, 64);
	return w;
}

// this code hasn't been fully tested and isn't currently used
int ui_lehman_factor (unsigned long p[MPZ_MAX_UI_PP_FACTORS], unsigned long h[MPZ_MAX_UI_PP_FACTORS], unsigned long n)
{
	unsigned long b;
	register unsigned long a, m, c, d, e, r, k, s, t, s1;
	register int w;

	// Lehman's O(n^{1/3}) deterministic algorithm - [CANT] p. 425
	b = (unsigned)pow(n,1.0/3.0)+1;	// add 1 to be safe
	t = b*b;
	s1 = n<<2;
	s = 0;
	for ( k = 1 ; k < b ; k++ ) {
		if ( k&0x1 ) {
			r = (k+n)&0x3;
			m = 4;
		} else {
			r = 1;
			m = 2;
		}
		s += s1;
		a = _ui_sqrt(s);
		while ( (a&(m-1)) != r ) a++;
		e = s+t;
		while ( (c=a*a) <= e ) {
			c -= s;
			if ( (c=_ui_sqrt(c)) ) goto lehman_done;
			a += m;
		}
	}
	p[w] = n;
	h[w++] = 1;
	return  w;
lehman_done:
	d = ui_gcd(a+c,n);
	if ( d == 1 || d == n ) { err_printf ("Lehman's algorithm failed for n = %lu with a = %lu, b = %lu, c = %lu, s = %lu, s+t = %lu, k = %lu, r = %lu, m = %lu\n", n, a, b, c, s, s+t, k, r, m);  exit (0); }
	p[w] = d;
	p[w+1] = n/d;
	h[w] = h[w+1] = 1;
	if ( p[w] > p[w+1] ) { p[w] = p[w+1];  p[w+1] = d; }		// keep primes in increasing order
	if ( p[w] == p[w+1] ) h[w++] = 2; else w += 2;			// check for square
	return w;
}


// returns primes in order
int mpz_factor_small (unsigned long p[], unsigned long h[], mpz_t bigp, mpz_t n, int max_factors, int max_hard_bits)
{
	static int init;
	static mpz_t d, x, m;
	unsigned long pt, ht;
	int i, j, w;
	
	if ( ! init ) { mpz_util_init();  mpz_init (d);  mpz_init (x);   mpz_init (m);  init = 1; }
	if ( ! mpz_sgn (n) ) { err_printf ("Attempt to factor 0");  exit (0); }
	mpz_set_ui (bigp, 1);
	if ( mpz_sgn(n) < 0 ) mpz_neg (m, n); else mpz_set (m,n);
	if ( mpz_cmp_ui(m,1) == 0 ) return 0;
	mpz_gcd (d, m, mpz_util_primorial);
	w = 0;
	for ( i = 1 ; i <= MPZ_SMALL_PRIMES ; i++ ) {
		if ( mpz_cmp_ui (d, 1) == 0 ) break;
		if ( mpz_divisible_ui_p (d, mpz_util_primes[i]) ) {
			p[w] = mpz_util_primes[i];
			mpz_set_ui (x, p[w]);
			h[w] = mpz_remove (m, m, x);
			w++;
			if ( w == max_factors ) { err_printf ("Exceeded max_factors %d in mpz_factor_small\n", max_factors);  exit (0); }
			mpz_remove (d, d, x);
		}
	}
	if ( mpz_cmp_ui (m, 1) == 0 ) return w;
	if ( mpz_sizeinbase (m, 2) > max_hard_bits ) return -1;
	while ( ! mpz_probab_prime_p (m, 5) ) {
		mpz_pollard_rho (d, m);
		if ( ! mpz_fits_ulong_p (d) ) { err_printf ("factor too big to fit in unsigned long!");   return -1; }
		pt = mpz_get_ui (d);
		ht = mpz_remove (m, m, d);
		for ( i = 0 ; i < w ; i++ ) if ( p[i] > pt ) break;
		for ( j = w ; j-1 >= i ; j-- ) { p[j] = p[j-1];  h[j] = h[j-1]; }
		p[j] = pt;
		h[j] = ht;			
		w++;
		if ( mpz_cmp_ui (m,1) == 0 ) return w;
		if ( w == max_factors ) { err_printf ("Exceeded max_factors %d in mpz_factor_small\n", max_factors);  exit (0); }
	}
	if ( ! mpz_fits_ulong_p (m) ) {
		mpz_set (bigp, m);
	} else {
		mpz_set_ui (bigp, 1);
		pt = mpz_get_ui (m);
		ht = 1;
		for ( i = 0 ; i < w ; i++ ) if ( p[i] > pt ) break;
		for ( j = w ; j-1 >= i ; j-- ) { p[j] = p[j-1];  h[j] = h[j-1]; }
		p[j] = pt;
		h[j] = ht;			
		w++;
	}
	return w;
}


int mpz_pollard_rho (mpz_t d, mpz_t n)
{
	static int init;
	static mpz_t x, y, z, t, P;
	register unsigned c, i, j, k, m;

	if ( ! init ) { mpz_util_init();  mpz_init (x);  mpz_init (y);  mpz_init(t);  mpz_init (P); mpz_init (z); init = 1; }
	m = 0;
	for ( c = 1 ;; c++ ) {
		mpz_set_ui (x, 2);
		mpz_set_ui (y, 2);
		i = 1;
		k = 2;
		for (;;) {
			mpz_set_ui (P, 1);
			mpz_set (t, x);
			for ( j = 0 ; i < k && j < 20 ; j++ ) {
				m++;
				i++;
				mpz_mul (x, x, x);
				mpz_add_ui (x, x, c);
				mpz_mod (x, x, n);
				mpz_sub (z, y, x);
				mpz_mulm (P, P, z, n);
			}
			mpz_gcd (d, P, n);
			if ( mpz_cmp_ui (d,1) != 0 ) {
				mpz_set (x,t);
				do {
					mpz_mul (x, x, x);
					mpz_add_ui (x, x, c);
					mpz_mod (x, x, n);
					mpz_sub (z, y, x);
					mpz_gcd (d, z, n);
				} while ( mpz_cmp_ui (d,1) == 0 );
				if ( mpz_cmp (d,n) != 0 ) return m;
				break;
			}
			if ( i == k ) {
				mpz_set (y,x);
				k *= 2;
			}
		}
	}
	return 1;
}


// computes the y-coarse part of x
int mpz_coarse_part (mpz_t o, mpz_t x, mpz_t y)
{
	static int init;
	static mpz_t d, z;
	int i, j, w;
	
	if ( ! init ) { mpz_util_init();  mpz_init (d);  mpz_init (z);   init = 1; }
	if ( mpz_cmp (x,y) <= 0 ) { mpz_set_ui (o,1);  return 1; }
	mpz_gcd (d, x, mpz_util_primorial);
	mpz_set (o, x);
	for ( i = 1 ; i <= MPZ_SMALL_PRIMES ; i++ ) {
		if ( mpz_cmp_ui (y, mpz_util_primes[i]) < 0 ) break;
		if ( mpz_divisible_ui_p (d, mpz_util_primes[i]) ) {
			mpz_set_ui (z, mpz_util_primes[i]);
			mpz_remove (o, o, z);
		}
	}
	if ( i <= MPZ_SMALL_PRIMES ) return 1;
	if ( mpz_cmp (o,y) <= 0 ) { mpz_set_ui (o, 1);  return 1; }
	mpz_set_ui (d, MPZ_MAX_SMALL_PRIME);
	while ( ! mpz_probab_prime_p (o, 5) ) {
		while ( ! mpz_divisible_p (o, d) ) {
			mpz_nextprime (d, d);
			if ( mpz_cmp (d, y) > 0 ) return 1;
		}
		mpz_remove (o, o, d);
	}
	return 1;
}

 
int mpz_eval_expr (mpz_t o, char *expr)
{
	static int init;
	static mpz_t p, x, y, z;
    char *s, *t, op, nextop;
    int i, digits, n;

 	if ( ! init ) { mpz_init (p);  mpz_init (x);  mpz_init (y);  mpz_init (z);  init = 1; }
    	mpz_util_init();	
    s = expr;
    if ( *s == 'D' || *s == 'd' ) {
        digits = atoi (s+1);
		mpz_ui_pow_ui (x, 10, digits);
		do {
			mpz_urandomm (o, mpz_util_rands, x);
		} while ( ! mpz_tstbit (o, 0) );
		if ( mpz_congruent_ui_p (o, 1, 4) ) mpz_mul_2exp (o, o, 2);
		info_printf ("Random Discriminant: %Zd\n", o);
		return 1;
    }
    if ( *s == 'R' || *s == 'r' ) {
        digits = atoi (s+1);
		mpz_ui_pow_ui (x, 10, digits);
		mpz_urandomm (o, mpz_util_rands, x);
		info_printf ("Random Number: %Zd\n", o);
		return 1;
    }
    if ( *s == 'P' || *s == 'p' ) {
        digits = atoi (s+1);
		mpz_ui_pow_ui (x, 10, digits);
		mpz_urandomm (o, mpz_util_rands, x);
		do {
			mpz_nextprime (o, o);
		} while ( ! mpz_probab_prime_p (o, 20) );			// We really want to be sure here so we don't screw up test cases
		info_printf ("Random Prime: %Zd\n", o);
		return 1;
    }
    if ( *s == 'B' || *s == 'b' ) {
        digits = atoi (s+1);
		mpz_ui_pow_ui (x, 10, digits);
		mpz_urandomm (o, mpz_util_rands, x);
		do {
			mpz_nextprime (o, o);
		} while ( ! mpz_congruent_ui_p (o, 3, 4) || ! mpz_probab_prime_p (o, 20) );			// We really want to be sure here so we don't screw up test cases
		info_printf ("Random Prime 3mod4: %Zd\n", o);
		return 1;
    }
    if ( *s == 'C' || *s == 'c' ) {
        digits = atoi (s+1);
		mpz_ui_pow_ui (x, 10, digits/2);
		mpz_urandomm (y, mpz_util_rands, x);
		mpz_nextprime (y, y);
		mpz_urandomm (z, mpz_util_rands, x);
		mpz_nextprime (z, z);
		mpz_mul (o, y, z);
		info_printf ("Random Composite of 2 Primes: %Zd\n", o);
		return 1;
    }
   mpz_set_ui (o, 1);
    mpz_set_ui (x, 1);
    op = '*';
    while ( expr_valid_op(op) ) {
	    t = s;
	    while ( isdigit (*s) ) s++;
	    if ( s == t ) return 0;
	    nextop = *s;
	    *s++ = '\0';
	    mpz_set (y, x);
	    if ( mpz_set_str (x, t, 0) != 0 ) return 0;
		if ( nextop == '!' || nextop == '#' || nextop == '$' ) {
			if ( nextop == '!' ) {
				n = atoi (t);
				mpz_fac_ui (x, n);
			} else if ( nextop == '#' ) {
				mpz_set_ui (p,1);
				mpz_set_ui (z,1);
				for(;;) {
					mpz_nextprime (p, p);
					if ( mpz_cmp (p, x) > 0 ) break;
					mpz_mul (z, z, p);
				}
				mpz_set (x, z); 
			} else if ( nextop == '$' ) {
				mpz_set_ui (p,1);
				mpz_set_ui (z,1);
				for ( i = 0 ; mpz_cmp_ui (x, i) > 0 ; i++ ) {
					mpz_nextprime (p, p);
					mpz_mul (z, z, p);
				}
				mpz_set (x, z); 
			}
			nextop = *s;
			*s++ = '\0';
		}
	    if ( op == 'E' || op == 'e' ) {
	    	n = atoi (t);
		  	if ( n > 0 ) {
		  		mpz_pow_ui (x, y, n-1);
		  		mpz_mul (o, o, x);
		  	}
	    } else if ( op == '*' ) {
	    	mpz_mul (o, o, x);
	    } else if ( op == '+' ) {
	    	mpz_add (o, o, x);
	    	return 1;
	    } else if ( op == '-' ) {
	    	mpz_sub (o, o, x);
	    	return 1;
	    }
	    op = nextop;
	}
    return 1;
}


static unsigned long mpz_mulm_counter, mpz_powm_counter, mpz_powm_tiny_counter;
static clock_t mpz_counter_reset_time;

void mpz_mulm (mpz_t o, mpz_t a, mpz_t b, mpz_t m)
    { mpz_mul (o, a, b);  mpz_mod (o, o, m); mpz_mulm_counter++; }

void mpz_addm (mpz_t o, mpz_t a, mpz_t b, mpz_t m)
    { mpz_add (o, a, b);  mpz_mod (o, o, m); }

void mpz_subm (mpz_t o, mpz_t a, mpz_t b, mpz_t m)
    { mpz_sub (o, a, b);  mpz_mod (o, o, m); }

// assumes input already mod m
void mpz_negm (mpz_t a, mpz_t m)
	{ if ( mpz_sgn(a) ) mpz_sub (a, m, a); }

void mpz_subm_ui (mpz_t o, mpz_t a, unsigned long b, mpz_t m)
    { mpz_sub_ui (o, a, b);  mpz_mod (o, o, m); }
    
void mpz_set_i (mpz_t o, long n)
{
	if ( n < 0 ) { mpz_set_ui (o, (unsigned long)(-n));  mpz_neg (o, o); } else { mpz_set_ui (o, (unsigned long)n); }
}

void mpq_set_i (mpq_t o, long n)
{
	if ( n < 0 ) { mpq_set_ui (o, (unsigned long)(-n), 1UL);  mpq_neg (o, o); } else { mpq_set_ui (o, (unsigned long)n, 1UL); }
}

// Optimize small powers - 0, 1, and 2 - assume b < m
void mpz_powm_tiny (mpz_t o, mpz_t b, unsigned e, mpz_t m)
{
    switch (e) {
    case 0: mpz_set_ui (o, 1);  return;
    case 1: if ( o != b ) mpz_set (o, b);  return;
    case 2: mpz_mulm (o, b, b, m);  return;
    default: mpz_powm_ui (o, b, e, m); mpz_powm_tiny_counter;
    }
    return;
}

void mpz_powm_big (mpz_t o, mpz_t b, mpz_t e, mpz_t m)
    { mpz_powm (o, b, e, m);  mpz_powm_counter++; }

void mpz_reset_counters ()
{
    mpz_mulm_counter = 0;
    mpz_powm_counter = 0;
    mpz_powm_tiny_counter = 0;
    mpz_counter_reset_time = clock();
}

void mpz_report_counters ()
{
    info_printf ("Counters: %d mult, %d tiny exp, %d exp, %d msec\n", mpz_mulm_counter, mpz_powm_tiny_counter, mpz_powm_counter,
    	    	 delta_msecs(mpz_counter_reset_time, clock()));
}

// Algorithm 1.5.1 in [CANT] - standard Tonelli-Shanks n^2 algorithm
int mpz_sqrt_modprime (mpz_t o, mpz_t a, mpz_t p)
{
	static mpz_t q, x, y, b;
	static int init;
	int i, r, m, e;

	if ( ! init ) { mpz_util_init();  mpz_init (q);  mpz_init (x);  mpz_init (y);  mpz_init (b);  init = 1; }
	if ( ! mpz_sgn (a) ) { mpz_set_ui (o, 0);  return 1; }
	if ( mpz_cmp_ui (a, 1) == 0 ) { mpz_set_ui (o, 1);  return 1; }
	mpz_sub_ui (q, p, 1);
	mpz_set_ui (x, 2);
	e = mpz_remove (q, q, x);
	do {
		mpz_urandomm (x, mpz_util_rands, p);
	} while ( mpz_jacobi (x, p) != -1 );
	mpz_powm_big (y, x, q, p);
	r = e;
	mpz_sub_ui (q, q, 1);
	mpz_divexact_ui (q, q, 2);
	mpz_powm_big (x, a, q, p);
	mpz_mulm (b, a, x, p);
	mpz_mulm (b, b, x, p);
	mpz_mulm (x, a, x, p);
	for (;;) {
		if ( mpz_cmp_ui (b, 1) == 0 ) { mpz_set (o, x);  break; }
		mpz_set (q, b);
		for ( m = 1 ; m <= r; m++ ) {
			mpz_mulm (q, q, q, p);
			if ( mpz_cmp_ui (q, 1) == 0 ) break;
		}
		if ( m > r ) { err_printf ("Unexpected result: m=%d > r=%d in mpz_sqrt_modprime?!  a = %Zd, p = %Zd\n", m, r, a, p);  exit (1); }
		if ( m == r ) return 0;
		mpz_set (q, y);
		for ( i = 0 ; i < r-m-1 ; i++ ) mpz_mulm (q, q, q, p);
		mpz_mulm (y, q, q, p);
		r = m;
		mpz_mulm (x, x, q, p);
		mpz_mulm (b, b, y, p);
	}
///	dbg_printf ("%Zd is a square root of %Zd mod %Zd\n", o, a, p);
	return 1;
}

// these are horribly inefficient
unsigned long mpz_get_bits_ui (mpz_t a, unsigned long i, unsigned long n)		// gets bits i thru i+n-1 and returns as ui
{
	register unsigned long j, k, x, y;

/*	y = 0;
	for ( j = i ; j < i+n ; j++ ) {
		if ( mpz_tstbit (a, j) ) y |= (1<<(j-i));
	}
	return y;
*/
	if ( mp_bits_per_limb != ULONG_BITS ) { err_printf ("mpz_get_bits_ui needs mp_bits_per_limb == ULONG_BITS!\n");  exit (0); }
	j = i  >> ULONG_BITS_LOG2;
	if ( j >= a->_mp_size ) return 0;
	x = a->_mp_d[j];
	k = i-(j<<ULONG_BITS_LOG2);		// k = i mod ULONG_BITS < ULONG_BITS
	x >>= k;
	k = ULONG_BITS-k;				// k = # bits in x > 0
	if ( k < n ) {
		if ( ++j >= a->_mp_size ) return x;
		x |= (a->_mp_d[j]&((1UL<<(n-k))-1))<<k;	// get n-k more bits
	} else if ( k > n ) {
		x &= (1UL<<n)-1;
	}

//	if ( x != y ) { err_printf ("bug in mpz_get_bits_ui(%Zx(hex),i,n) = %d != %d\n", a, i, n, x, y);  exit (0); }

	return x;
}

// Note bits is assumed to contain only the n bits to be set, the higher order bits must be zero
unsigned long mpz_set_bits_ui (mpz_t a, unsigned long i, unsigned long n, unsigned long bits)		// sets bits i thru i+n-1 and returns as ui
{
	register unsigned long j, k, k2, x;
	
/*
	for ( j = i ; j < i+n ; j++ ) {
		if ( bits&(1<<(j-i)) ) {
			mpz_setbit (a, j);
		} else {
			mpz_clrbit (a, j);
		}
	}
	return bits;
*/	
	if ( ! n ) return 0;
	if ( mp_bits_per_limb != ULONG_BITS ) { err_printf ("mpz_get_bits_ui needs mp_bits_per_limb == ULONG_BITS!\n");  exit (0); }
	j = i  >> ULONG_BITS_LOG2;
	if ( j >= a->_mp_size ) { mpz_setbit (a, i+n-1); }	// force size adjustment and realloc if required
	k = i-(j<<ULONG_BITS_LOG2);		// k = i mod ULONG_BITS < ULONG_BITS
	x = bits<<k;
	x |= a->_mp_d[j] & ((1UL<<k)-1);		// get k low bits that should remain the same
	k2 = ULONG_BITS-k;				// k2 = # bits in set so far > 0
	if ( k2 < n ) {
		a->_mp_d[j] = x;
		if ( ++j >= a->_mp_size ) { mpz_setbit (a, i+n-1); }	// force size adjustment and realloc if required
		x = bits >>k2;
		x |= (a->_mp_d[j]&~((1UL<<(n-k2))-1));	// get high bits that should remain the same
		a->_mp_d[j] = x;
	} else if ( k2 > n ) {
		x |= a->_mp_d[j] & ~((1UL<<(k+n))-1); // get high bits past k+n bits that should remain the same
		a->_mp_d[j] = x;
	} else {
		a->_mp_d[j] = x;
	}
	while ( a->_mp_size && ! a->_mp_d[a->_mp_size-1] ) a->_mp_size--;
	return bits;
}

unsigned long ui_randomm (unsigned long m)
{
	mpz_util_init(); 
	return gmp_urandomm_ui (mpz_util_rands, m);
}

unsigned long ui_randomb (unsigned long b)
{
	mpz_util_init(); 
	return gmp_urandomb_ui (mpz_util_rands, b);
}

unsigned long ui_pp_div (unsigned long n, unsigned long p)
{
	unsigned long q;
	
	if ( ! n ) return 1;
	q = 1;
	while ( (n%p) == 0 ) {
		n /= p;
		q *= p;
	}
	return q;
}


unsigned long ui_inverse (unsigned long a, unsigned long m)
{
	register unsigned long q, r, r0, r1;		// amazingly in this day and age, the register declaration gives a 10% improvement
	register long t, t0, t1;

	if ( a >= m ) a = a%m;
	if ( a == 0 ) return 0;

	t1 = 1;
	t0 = 0;
	r0 = m;
	r1 = a;
	while ( r1 > 0 ) {
		q = r0/r1;
		r = r0 - q*r1;
		r0 = r1;
		r1 = r;
		t = t0 - q*t1;
		t0 = t1;
		t1 = t;
	}
	if ( r0 > 1 ) return 0;
	if ( t0 < 0 ) return m - ((unsigned long)(-t0));
	return (unsigned long)t0;
}

// return 1 if the integer represented by the prime factorization (p,h,w) <ULONG_MAX contains a divisor in [min,max]
// slow implementation
int ui_divisor_in_interval (unsigned long p[], unsigned long h[], int w, unsigned long min, unsigned long max)
{
	register unsigned long n;
	unsigned long e[MPZ_MAX_FACTORS+1];
	register int i, j;
	
	if ( min > max ) return 0;
	if ( min <= 1 ) return 1;
	for ( i = 0 ; i < w ; i++ ) e[i] = 0;
	for (;;) {
		n = 1;
		for ( j = 0 ; j < w ; j++ ) {
			if ( e[j] ) {
				n *= p[j];
				for ( i = 1 ; i < e[j] ; i++ ) n *= p[j];
			}
		}
//		printf ("n=%d\n", n);
		if ( n >= min && n <= max ) return 1;
		e[0]++;
		for ( i = 0 ; e[i] > h[i] ; i++ ) { if ( i == w-1 ) return 0;  e[i] = 0;  e[i+1]++; }
	}
}

// a is compatible with b if it is composed entirely of primes dividing b
int ui_compatible  (unsigned long a, unsigned long b)
{
	unsigned long d;
	
	while ( a != 1 ) {
		d = ui_gcd(a,b);
		if ( d == 1 ) return 0;
		a /= d;
	}
	return 1;	
}

// a is compatible with b if it is composed entirely of primes dividing b
int mpz_compatible  (mpz_t a, mpz_t b)
{
	static int init;
	static mpz_t t, d;
	
	if ( ! init ) { mpz_init(d);  mpz_init(t);  init = 1; }
	mpz_set (t, a);
	while ( mpz_cmp_ui (t,1) != 0 ) {
		mpz_gcd (d, t, b);
		if ( mpz_cmp_ui(d,1) == 0 ) return 0;
		mpz_divexact(t,t,d);
	}
	return 1;	
}

unsigned long ui_crt (unsigned long a, unsigned long M, unsigned long b, unsigned long N)
{
	unsigned long x;
	
	if ( a > b ) {
		x = ui_inverse (N, M);
		return (x*N*(a-b)+b)%(M*N);
	} else {
		x = ui_inverse (M, N);
		return (x*M*(b-a)+a)%(M*N);
	}
}

unsigned long ui_gcd (unsigned long a, unsigned long b)
	{  return (ui_gcd_ext (a, b, NULL, NULL)); }

unsigned long ui_gcd_ext (unsigned long a, unsigned long b, long *x, long *y)
{
	register unsigned long q, r, r0, r1;
	register long s, t, s0, s1, t0, t1;
	
	if ( a < b ) return (ui_gcd_ext (b, a, y, x));
	if ( b == 0 ) {
		if ( x ) *x = 1;
		if ( y ) *y = 0;
		return a;
	}
	if ( x ) { s0 = 1;  s1 = 0; }
	if ( y ) { t1 = 1;  t0 = 0; }
	r0 = a;
	r1 = b;
	while ( r1 > 0 ) {
		q = r0/r1;
		r = r0 - q*r1;
		r0 = r1;
		r1 = r;
		if ( y ) {
			t = t0 - q*t1;
			t0 = t1;
			t1 = t;
		}
		if ( x ) {
			s = s0 - q*s1;
			s0 = s1;
			s1 = s;
		}
	}
	if ( x ) *x = s0;
	if ( y ) *y = t0;
	return r0;
}

// Assumes gcd(a,b,N) = 1 but does not verify this
unsigned long bach_gcd (long a, long b, unsigned long N)
{
	unsigned long c, g, h;
	unsigned long M, K;
	
	g = ui_gcd(abs(a),N);
	if ( g == 1 ) return 0;
	if ( ui_gcd(abs(a+b),N) == 1 ) return 1;
	M = N;
	for ( h = g ; h > 1 ; ) {
		do { M /= h; } while ( (M%h) == 0 );
		h = ui_gcd (M, g);			// could use h instead of g?
	}
	// M>1 should now be a factor of N that is relatively prime to both g and N
	if ( M <= 1 || (N%M)!=0 ) { out_printf ("Error M=%d in bach_gcd\n", M);  exit (1); }
	c = (M*ui_inverse(M,N/M)) % N;
	if ( c == 0 ) { out_printf ("Error c == 0 in bach_gcd\n", c);  exit (1); }
	return c;
}




int mpz_eval_term_ui (mpz_t o, unsigned long numvars, unsigned long vars[], unsigned long exps[])
{
	static int init;
	static mpz_t x;
	int i;
	
	if ( ! init ) { mpz_init (x);  init = 1; }
	
	mpz_set_ui (o, 1);
	for ( i = 0 ; i < numvars ; i++ ) {
		if ( exps[i] > 0 ) {
			mpz_ui_pow_ui (x, vars[i], exps[i]);
			mpz_mul (o, o, x);
		}
	}
	return 1;
}


char *ui_term_to_string (char *buf, unsigned long numvars, unsigned long vars[], unsigned long exps[])
{
	char *s;
	int i;
	
	if ( ! numvars ) { strcpy (buf, "1");  return buf; }
	s = buf;
	*s = '\0';
	for ( i = 0 ; i < numvars ; i++ ) {
		if ( exps[i] > 0 ) {
			if ( s > buf ) *s++ = '*';
			sprintf (s, "%d^%d", vars[i], exps[i]);
			while ( *s ) s++;
		}
	}
	return buf;
}

unsigned long ui_wt (unsigned long x)
{
	unsigned long i, n;
	
	n = 0;
	for ( i = 1 ; i <= x ; i <<= 1 ) {
		if ( i&x ) n++;
	}
	return n;
}

unsigned long ui_lg_ceil (unsigned long x)
{
	unsigned long i;

	i = _asm_highbit(x);
	if ( x==(1UL<<i) ) return i; else return i+1;
}

unsigned long ui_binexp_cost (unsigned long x)
	{ return ui_len(x) + ui_wt(x) - 2; }

unsigned long ui_get_bits (unsigned long x, unsigned long pos, unsigned long bits)
{
	unsigned long i, m;
	
	m = 0;
	for ( i = pos ; i < pos+bits ; i++ ) m |= (1<<i);
	m &= x;
	return (m >> pos);
}

char *ui_bit_string (char *buf, unsigned long x)
{
	int i;
	
	for ( i = ui_len(x)-1 ; i >= 0 ; i-- ) *buf++ = ( (x & (1<<i)) ? 1 : 0 );
	*buf = '\0';
	return buf;
}

void ui_invert_permutation (unsigned long p2[], unsigned long p1[], unsigned long n)
{
	unsigned long i;
	
	for ( i = 0 ; i < n ; i++ ) p2[p1[i]] = i;
}


double mpz_log2 (mpz_t a)
{
	unsigned long b, n;
	
	b = mpz_sizeinbase(a,2);
	if ( ! b ) return 0.0;			// define lg(0) = 0 for convenience
	if ( b < ULONG_BITS ) {
		n = mpz_get_ui (a);
		return log2((double)n);
	} else {
		n = mpz_get_bits_ui (a, b-(ULONG_BITS-1), ULONG_BITS-2);
		return (double)(b-1) + log2(1.0+(double)n/(double)(1UL<<(ULONG_BITS-2)));
	}
}


int ui_qsort_cmp (const void *a, const void *b)
{
	if ( ((unsigned long*)a) > ((unsigned long *)b) ) return 1;
	if ( ((unsigned long*)b) > ((unsigned long *)a) ) return -1;
	return 0;
}

int _qsort_mpz_cmp (const void *a, const void *b)
{
	return mpz_cmp (*((mpz_t *)a), *((mpz_t *)b));
}


void mpz_sort (mpz_t a[], unsigned long n)
{
	qsort (a, n, sizeof(a[0]), _qsort_mpz_cmp);
}


unsigned long mpz_randomf (mpz_t o, mpz_t N, mpz_t factors[], unsigned long w)
{
	return mpz_randomft (o, N, factors, NULL, w);
}


/*
	The algorithm below is an implementation of Bach's algorithm for computing random factored integers.
	See "Analytic Methods in the Analysis and Design of Number-Theoretic Algorithms", Eric Bach, MIT Press 1984
	Returns the factorization of a random integer in the range (N/2,N]
*/

unsigned long mpz_randomft (mpz_t o, mpz_t N, mpz_t factors[], mpz_ftree_t factor_trees[], unsigned long w)
{
	mpz_t d, Q, M;
	mpf_t x;
	mpz_ftree_t t;
	double b, u;
	unsigned long n;
	int i, j, k, m, q;
	
	mpz_init (d);
	mpz_init (M);
	mpz_util_init ();

	if ( w == 0 ) { err_printf ("Factor array too small - increase size!\n");  exit (1); }

	// Base case
	if ( mpz_cmp_ui (N, MPZ_MAX_SMALL_INTEGER) <= 0 ) {
		// pick random r in (N/2,N]
		mpz_fdiv_q_2exp (M, N, 1);
		mpz_sub (d, N, M);
		mpz_urandomm (o, mpz_util_rands, d);
		mpz_add (o, o, M);
		mpz_add_ui (o, o, 1);
		// factor r using small factors array
		n = mpz_get_ui (o);
		k = 0;
		while ( n > 1 ) {
			if ( k >= w ) { err_printf ("Factor array too small - increase size!\n");  exit (1); }
			q = mpz_small_factors[n];
			if ( ! q ) q = n;
			if ( n%q ) { err_printf ("Small factor array assert failed - %d not divisible by %d\n", n, q);  exit (1); }
			n /= q;
			for ( i = 0 ; i < k ; i++ ) {
				if ( mpz_divisible_ui_p (factors[i], q) ) {
					mpz_mul_ui (factors[i], factors[i], q);
					break;
				}
			}
			if ( i == k ) mpz_set_ui (factors[k++], q);
		}
		mpz_clear (d);
		mpz_clear (M);
		return k;
	}
	
	mpz_init (Q);
	mpf_init (x);
	for (;;) {
		if ( factor_trees ) {
			t = mpz_randomfpp (Q, N);
		} else {
			mpz_randompp (Q, N);
		}
		mpz_fdiv_q (M, N, Q);
		k = mpz_randomf (o, M, factors, w-1);
		for ( i = 0 ; i < k ; i++ ) {
			if ( mpz_divisible_p (factors[i], Q) ) {
				mpz_mul (factors[i], factors[i], Q);
				break;
			}
		}
		if ( i == k ) {
			if ( factor_trees ) factor_trees[k] = t;
			mpz_set (factors[k++], Q);
		} else {
			if ( factor_trees ) mpz_clear_ftree (t);
		}
		mpz_set_ui (o, 1);
		for ( i = 0 ; i < k ; i++ ) mpz_mul (o, o, factors[i]);
		b = (mpz_log2 (N) - 1.0) / mpz_log2 (o);
		mpf_urandomb (x, mpz_util_rands, 64);
		u = mpf_get_d (x);
		if ( u <= b ) break;
dbg_printf ("random reject %f > %f\n", u, b);
		if ( factor_trees ) for ( i = 0 ; i < k ; i++ ) mpz_clear_ftree(t);
	}
	// Need to clean-up factorizations of p-1 to make sure we deal with cases where we generated a prime power with
	// the factorization of p^k-1.  Fortunately (p-1)|p^k-1, so all the prime factors of p-1 are present in the factorization
	// of p^k-1, we just need to pick them out.
	if ( factor_trees ) {
dbg_printf ("cleaning up\n");
		for ( i = 0 ; i < k ; i++ ) {
			if ( ! mpz_perfect_power_p (factors[i]) ) continue;
			if ( ! mpz_pp_base (Q, factors[i]) ) { err_printf ("Unable to extract prime power base from a known prime power!");  exit (0); }
			mpz_sub_ui (Q, Q, 1);
			m = 0;
			t = factor_trees[i];
			for ( j = 0 ; j < t->w ; j++ ) {
				mpz_gcd (d, Q, t->factors[j]);
				if ( mpz_cmp_ui (d,1) > 0 ) mpz_set (t->factors[m++], d);
			}
			for ( j = m ; j < t->w ; j++ ) mpz_clear (t->factors[j]);
			t->w = m;
			// sanity check the results
			mpz_set_ui (d, 1);
			for ( j = 0 ; j < t->w ; j++ ) mpz_mul (d, d, t->factors[j]);
			if ( mpz_cmp (Q, d) != 0 ) { err_printf ("Invalid factorization of p-1 after prime power cleanup!  %Zd != %Zd\n", d, Q);  exit (0); }
		}
	}
	mpz_clear (d);
	mpz_clear (M);
	mpz_clear (Q);
	mpf_clear (x);
	return k;
}

// Returns a random prime power in the range [2,N]
void mpz_randompp (mpz_t Q, mpz_t N)
{
	static int init;
	static mpz_t M, M2, J, p, d, r;
	static mpf_t x, y;
	double b, u;
	unsigned long i, j, n;

	if ( ! init ) {
		mpz_init (M);  mpz_init (M2); mpz_init (J);  mpz_init (p);  mpz_init (d);  mpz_init (r);  mpf_init (x);  mpf_init (y); 
		mpz_util_init ();
		init = 1;
	}
	for (;;) {
		// select random j in [1,log2(N)]
		n = (unsigned long)ceil (mpz_log2 (N));
		j = gmp_urandomm_ui (mpz_util_rands, n) + 1;
		mpz_ui_pow_ui (J, 2, j);
		mpz_urandomm (r, mpz_util_rands, J);
		mpz_add (Q, J, r);
		if ( mpz_cmp (Q, N) > 0 ) continue;
		if ( mpz_perfect_power_p (Q) ) {
			if ( ! mpz_pp_base (p, Q) ) continue;
		} else {
			mpz_set (p, Q);
		}
		mpz_gcd (d, p, mpz_util_primorial);
		if ( mpz_cmp_ui (d,1) != 0 && mpz_cmp (d,p) != 0 ) continue;
		if ( ! mpz_probab_prime_p (p, 5) ) continue;
		mpf_urandomb (x, mpz_util_rands, 64);
		u = mpf_get_d (x);
		mpz_fdiv_q (M, N, Q);
		mpz_mul_ui (d, Q, 2);
		mpz_fdiv_q (M2, N, d);
		mpz_sub (M, M, M2);
//err_printf ("#(N/2q,N/q] = %Zd for N = %Zd, q = %Zd\n", M, N, Q);
		mpz_mul (M, M, J);
		mpf_set_z (x, M);
		mpf_set_z (y, N);
		mpf_div (x, x, y);
		b = mpf_get_d (x);
//err_printf ("b = %f, Q = %Zd, p = %Zd, log2(p) = %f, log2(N) = %f\n", b, Q, p, mpz_log2(p), mpz_log2(N));
		b *= mpz_log2(p) / mpz_log2(N);
		if ( u < b ) break;
//err_printf ("random reject %f >= %f\n", u, b);
	}
}

// optimize later
unsigned long ui_pp_base (unsigned long pp)
{
	static int init;
	static mpz_t x, b;
	
	if ( ! init ) { mpz_init (x);  mpz_init (b); }
	mpz_set_ui (x, pp);
	if ( mpz_probab_prime_p (x, 5) ) return mpz_get_ui (x);
	if ( ! mpz_pp_base (b, x) ) return 0;
	return mpz_get_ui (b);
}


/*
	Extracts base from a prime power and returns the exponent
*/
unsigned long mpz_pp_base (mpz_t b, mpz_t q)
{
	static int init;
	static mpz_t p, x, y, z;
	double d;
	int c, m, n;
	
	if ( ! init ) { mpz_util_init ();  mpz_init (p);  mpz_init (x);  mpz_init (y);  mpz_init (z);  init = 1; }
	mpz_gcd (p, q, mpz_util_primorial);
	if ( mpz_cmp_ui (p, 1) != 0 ) {
		if ( mpz_cmp_ui(p,MPZ_MAX_SMALL_PRIME) > 0 || ! ui_is_small_prime (mpz_get_ui(p)) ) return 0;		// q divisible by more than 1 small prime
		mpz_fdiv_q (x, q, p);
		n = 1;
		while ( mpz_cmp_ui (x,1) > 0 ) {
			mpz_fdiv_qr (x, z, x, p);
			if ( mpz_cmp_ui (z, 0) != 0 ) return 0;		// q is not a prime power
			n++;
		}
		mpz_set (b, p);
		return n;
	} else {
		// This is not particular efficient, but it rarely gets used.
		if ( mpz_probab_prime_p (q, 5) ) { mpz_set (b, q);  return 1; }
		n = mpz_sizeinbase (q,2)/ui_lg_floor(mpz_util_primes[MPZ_SMALL_PRIMES]);
		while ( n >= 0 ) {
			mpz_ui_pow_ui (z, 2, (mpz_sizeinbase (q,2) / n) + 1);		// upper bound
			mpz_ui_pow_ui (y, 2, (mpz_sizeinbase (q,2)-1) / n);			// lower bound
			do {
				mpz_add (x, y, z);
				mpz_tdiv_q_ui (x, x, 2);									// middle
				mpz_pow_ui (p, x, n);
				c = mpz_cmp (p, q);
				if ( ! c ) { if ( ! mpz_probab_prime_p (x, 5) ) return 0;  mpz_set (b, x);  return n; }
				if ( c < 0 ) {
					mpz_add_ui (y, x, 1);
				} else {
					mpz_sub_ui (z, x, 1);
				}
			} while ( mpz_cmp (y,z) <= 0 );
			n--;
		}
	}
	return 0;
}

mpz_ftree_t mpz_randomfp (mpz_t o, mpz_t N)
{
	static int init;
	static mpz_t N2;
	mpz_ftree_t t;
	int i;
	
	if ( ! init ) { mpz_init(N2);  init = 1; }
	mpz_sub_ui (N2, N, 1);
	mpz_fdiv_q_ui (N2, N2, 2);
	t = (struct factor_tree_item_struct *)mem_alloc (sizeof (*t));
	for ( i = 0 ; i < MPZ_MAX_TREE_FACTORS ; i++ ) mpz_init (t->factors[i]);
	for (;;) {
		t->w = mpz_randomf (o, N2, t->factors, MPZ_MAX_TREE_FACTORS);
		mpz_mul_ui (o, o, 2);
		mpz_add_ui (o, o, 1);
		if ( mpz_probab_prime_p (o, 5) ) break;
		out_printf (".");
	}
	mpz_sort (t->factors, t->w);
	for ( i = 0 ; i < t->w ; i++ ) {
		if ( ! mpz_tstbit (t->factors[i], 0) ) break;
	}
	if ( i < t->w  ) {
		mpz_mul_ui (t->factors[i], t->factors[i], 2);
	} else {
		if ( t->w >= MPZ_MAX_TREE_FACTORS ) { err_printf ("Exceeded MPZ_MAX_TREE_FACTORS!");  exit (0); }
		for ( i = t->w ; i > 0 ; i-- ) mpz_set (t->factors[i], t->factors[i-1]);
		mpz_set_ui (t->factors[0], 2);
		t->w++;
	}
	return t;	
}


mpz_ftree_t mpz_randomfpp (mpz_t o, mpz_t N)
{
	static int init;
	static mpz_t d, p, N2;
	mpz_ftree_t t;
	int i;
	
	if ( ! init ) { mpz_init (d);  mpz_init (p);  mpz_init(N2);  init = 1; }
	mpz_sub_ui (N2, N, 1);
	mpz_fdiv_q_ui (N2, N2, 2);
	t = (struct factor_tree_item_struct*)mem_alloc (sizeof (*t));
	for ( i = 0 ; i < MPZ_MAX_TREE_FACTORS ; i++ ) mpz_init (t->factors[i]);
	for (;;) {
		t->w = mpz_randomf (o, N2, t->factors, MPZ_MAX_TREE_FACTORS);
		out_printf (".");
		mpz_mul_ui (o, o, 2);
		mpz_add_ui (o, o, 1);
		if ( mpz_perfect_power_p (o) ) {
			if ( ! mpz_pp_base (p, o) ) continue;
		} else {
			mpz_set (p, o);
		}
		mpz_gcd (d, p, mpz_util_primorial);
		if ( mpz_cmp_ui (d,1) != 0 && mpz_cmp (d,p) != 0 ) continue;
		if ( ! mpz_probab_prime_p (p, 5) ) continue;
	}
	mpz_sort (t->factors, t->w);
	for ( i = 0 ; i < t->w ; i++ ) {
		if ( ! mpz_tstbit (t->factors[i], 0) ) break;
	}
	if ( i < t->w  ) {
		mpz_mul_ui (t->factors[i], t->factors[i], 2);
	} else {
		if ( t->w >= MPZ_MAX_TREE_FACTORS ) { err_printf ("Exceeded MPZ_MAX_TREE_FACTORS!");  exit (0); }
		for ( i = t->w ; i > 0 ; i-- ) mpz_set (t->factors[i], t->factors[i-1]);
		mpz_set_ui (t->factors[0], 2);
		t->w++;
	}
	return t;	
}


void mpz_clear_ftree (mpz_ftree_t t)
{
	int i;
	
	for ( i = 0 ; i < t->w ; i++ ) {
		if ( t->subtrees[i] ) mpz_clear_ftree (t->subtrees[i]);
	}
	for ( i = 0 ; i < MPZ_MAX_TREE_FACTORS ; i++ ) mpz_clear (t->factors[i]);
}


unsigned long mpz_store_bytes (mpz_t a)
{
	unsigned long size;
	
	size = _ui_ceil_ratio(mpz_sizeinbase(a,2),8);
	if ( size+2 >= (1<<15) ) { err_printf ("Integer %Zd too large to store!\n", a);  exit (0); }
	return size+2;
}


void mpz_store (char *s, mpz_t a)
{
	struct _group_storage_block *pStore;
	size_t size, count;
	short cnt;
	
	size = _ui_ceil_ratio(mpz_sizeinbase(a,2),8);
	if ( size+1 >= (1<<15) ) { err_printf ("Integer %Zd too large to store!\n", a);  exit (0); }
	mpz_export (s+sizeof(short), &count, 1, 1, 0, 0, a);
	if ( count >= (1<<15) ) { err_printf ("mpz_export wrote too many bytes for integer %Zd\n", a);  exit (0); }
	cnt = (short)count;
	if ( mpz_sgn(a) < 0 ) cnt = -cnt;
	*((short*)s) = cnt;
	return;	
}


void mpz_retrieve (mpz_t o, char *s)
{
	short cnt;
	int sign;
	
	cnt = *((short*)s);
	sign = 1;
	if ( cnt < 0 ) { sign = -1;  cnt = -cnt; }
	mpz_import (o, (unsigned long)cnt, 1, 1, 0, 0, s+2);
	if ( sign < 0 ) mpz_neg (o, o);
	return;
}


unsigned long ui_eval_expr (char *expr)
{
	char *s;
	unsigned long n;
	int i, j, k;

	for ( s = expr ; *s && *s != 'e' ; s++ );
	if ( *s ) {
		*s++ = '\0';
		k = atoi(expr);
		j = atoi(s);
		for ( i = 0, n = 1 ; i < j ; i++ ) n *= k;
	} else {
		n = atoll (expr);
	}
	return n;
}

NAF_entry_t ui_NAF_table[2][4] = { {{0,0}, {1,0}, {0,0}, {-1,1}}, {{1,0},{0,1},{-1,1},{0,1}} };
unsigned char NAF_pbits_tab[1024] = {0, 1, 2, 4, 4, 5, 8, 8, 8, 9, 10, 16, 16, 17, 16, 16, 16, 17, 18, 20, 20, 21, 32, 32, 32, 33, 34, 32, 32, 33, 32, 32, 32, 33, 34, 36, 36, 37, 40, 40, 40, 41, 42, 64, 64, 65, 64, 64, 64, 65, 66, 68, 68, 69, 64, 64, 64, 65, 66, 64, 64, 65, 64, 64, 64, 65, 66, 68, 68, 69, 72, 72, 72, 73, 74, 80, 80, 81, 80, 80, 80, 81, 82, 84, 84, 85, 128, 128, 128, 129, 130, 128, 128, 129, 128, 128, 128, 129, 130, 132, 132, 133, 136, 136, 136, 137, 138, 128, 128, 129, 128, 128, 128, 129, 130, 132, 132, 133, 128, 128, 128, 129, 130, 128, 128, 129, 128, 128, 128, 129, 130, 132, 132, 133, 136, 136, 136, 137, 138, 144, 144, 145, 144, 144, 144, 145, 146, 148, 148, 149, 160, 160, 160, 161, 162, 160, 160, 161, 160, 160, 160, 161, 162, 164, 164, 165, 168, 168, 168, 169, 170, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 8, 8, 8, 9, 10, 16, 16, 17, 16, 16, 16, 17, 18, 20, 20, 21, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 8, 8, 8, 9, 10, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 8, 8, 8, 9, 10, 16, 16, 17, 16, 16, 16, 17, 18, 20, 20, 21, 32, 32, 32, 33, 34, 32, 32, 33, 32, 32, 32, 33, 34, 36, 36, 37, 40, 40, 40, 41, 42, 64, 64, 65, 64, 64, 64, 65, 66, 68, 68, 69, 64, 64, 64, 65, 66, 64, 64, 65, 64, 64, 64, 65, 66, 68, 68, 69, 72, 72, 72, 73, 74, 80, 80, 81, 80, 80, 80, 81, 82, 84, 84, 85, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 8, 8, 8, 9, 10, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 8, 8, 8, 9, 10, 16, 16, 17, 16, 16, 16, 17, 18, 20, 20, 21, 32, 32, 32, 33, 34, 32, 32, 33, 32, 32, 32, 33, 34, 36, 36, 37, 40, 40, 40, 41, 42, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 8, 8, 8, 9, 10, 16, 16, 17, 16, 16, 16, 17, 18, 20, 20, 21, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 8, 8, 8, 9, 10, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 1, 2, 4, 4, 5, 8, 8, 8, 9, 10, 16, 16, 17, 16, 16, 16, 17, 18, 20, 20, 21, 32, 32, 32, 33, 34, 32, 32, 33, 32, 32, 32, 33, 34, 36, 36, 37, 40, 40, 40, 41, 42, 64, 64, 65, 64, 64, 64, 65, 66, 68, 68, 69, 64, 64, 64, 65, 66, 64, 64, 65, 64, 64, 64, 65, 66, 68, 68, 69, 72, 72, 72, 73, 74, 80, 80, 81, 80, 80, 80, 81, 82, 84, 84, 85, 128, 128, 128, 129, 130, 128, 128, 129, 128, 128, 128, 129, 130, 132, 132, 133, 136, 136, 136, 137, 138, 128, 128, 129, 128, 128, 128, 129, 130, 132, 132, 133, 128, 128, 128, 129, 130, 128, 128, 129, 128, 128, 128, 129, 130, 132, 132, 133, 136, 136, 136, 137, 138, 144, 144, 145, 144, 144, 144, 145, 146, 148, 148, 149, 160, 160, 160, 161, 162, 160, 160, 161, 160, 160, 160, 161, 162, 164, 164, 165, 168, 168, 168, 169, 170, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 8, 8, 8, 9, 10, 16, 16, 17, 16, 16, 16, 17, 18, 20, 20, 21, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 8, 8, 8, 9, 10, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 8, 8, 8, 9, 10, 16, 16, 17, 16, 16, 16, 17, 18, 20, 20, 21, 32, 32, 32, 33, 34, 32, 32, 33, 32, 32, 32, 33, 34, 36, 36, 37, 40, 40, 40, 41, 42, 64, 64, 65, 64, 64, 64, 65, 66, 68, 68, 69, 64, 64, 64, 65, 66, 64, 64, 65, 64, 64, 64, 65, 66, 68, 68, 69, 72, 72, 72, 73, 74, 80, 80, 81, 80, 80, 80, 81, 82, 84, 84, 85, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 8, 8, 8, 9, 10, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 8, 8, 8, 9, 10, 16, 16, 17, 16, 16, 16, 17, 18, 20, 20, 21, 32, 32, 32, 33, 34, 32, 32, 33, 32, 32, 32, 33, 34, 36, 36, 37, 40, 40, 40, 41, 42, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 8, 8, 8, 9, 10, 16, 16, 17, 16, 16, 16, 17, 18, 20, 20, 21, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 8, 8, 8, 9, 10, 0, 0, 1, 0, 0, 0, 1, 2, 4, 4, 5, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 0, };
unsigned char NAF_nbits_tab[1024] = {0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 10, 9, 8, 8, 8, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 21, 20, 20, 18, 17, 16, 16, 16, 17, 16, 16, 10, 9, 8, 8, 8, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 42, 41, 40, 40, 40, 37, 36, 36, 34, 33, 32, 32, 32, 33, 32, 32, 34, 33, 32, 32, 32, 21, 20, 20, 18, 17, 16, 16, 16, 17, 16, 16, 10, 9, 8, 8, 8, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 10, 9, 8, 8, 8, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 85, 84, 84, 82, 81, 80, 80, 80, 81, 80, 80, 74, 73, 72, 72, 72, 69, 68, 68, 66, 65, 64, 64, 64, 65, 64, 64, 66, 65, 64, 64, 64, 69, 68, 68, 66, 65, 64, 64, 64, 65, 64, 64, 42, 41, 40, 40, 40, 37, 36, 36, 34, 33, 32, 32, 32, 33, 32, 32, 34, 33, 32, 32, 32, 21, 20, 20, 18, 17, 16, 16, 16, 17, 16, 16, 10, 9, 8, 8, 8, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 10, 9, 8, 8, 8, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 21, 20, 20, 18, 17, 16, 16, 16, 17, 16, 16, 10, 9, 8, 8, 8, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 170, 169, 168, 168, 168, 165, 164, 164, 162, 161, 160, 160, 160, 161, 160, 160, 162, 161, 160, 160, 160, 149, 148, 148, 146, 145, 144, 144, 144, 145, 144, 144, 138, 137, 136, 136, 136, 133, 132, 132, 130, 129, 128, 128, 128, 129, 128, 128, 130, 129, 128, 128, 128, 133, 132, 132, 130, 129, 128, 128, 128, 129, 128, 128, 138, 137, 136, 136, 136, 133, 132, 132, 130, 129, 128, 128, 128, 129, 128, 128, 130, 129, 128, 128, 128, 85, 84, 84, 82, 81, 80, 80, 80, 81, 80, 80, 74, 73, 72, 72, 72, 69, 68, 68, 66, 65, 64, 64, 64, 65, 64, 64, 66, 65, 64, 64, 64, 69, 68, 68, 66, 65, 64, 64, 64, 65, 64, 64, 42, 41, 40, 40, 40, 37, 36, 36, 34, 33, 32, 32, 32, 33, 32, 32, 34, 33, 32, 32, 32, 21, 20, 20, 18, 17, 16, 16, 16, 17, 16, 16, 10, 9, 8, 8, 8, 5, 4, 4, 2, 1, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 10, 9, 8, 8, 8, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 21, 20, 20, 18, 17, 16, 16, 16, 17, 16, 16, 10, 9, 8, 8, 8, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 42, 41, 40, 40, 40, 37, 36, 36, 34, 33, 32, 32, 32, 33, 32, 32, 34, 33, 32, 32, 32, 21, 20, 20, 18, 17, 16, 16, 16, 17, 16, 16, 10, 9, 8, 8, 8, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 10, 9, 8, 8, 8, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 85, 84, 84, 82, 81, 80, 80, 80, 81, 80, 80, 74, 73, 72, 72, 72, 69, 68, 68, 66, 65, 64, 64, 64, 65, 64, 64, 66, 65, 64, 64, 64, 69, 68, 68, 66, 65, 64, 64, 64, 65, 64, 64, 42, 41, 40, 40, 40, 37, 36, 36, 34, 33, 32, 32, 32, 33, 32, 32, 34, 33, 32, 32, 32, 21, 20, 20, 18, 17, 16, 16, 16, 17, 16, 16, 10, 9, 8, 8, 8, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 10, 9, 8, 8, 8, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 21, 20, 20, 18, 17, 16, 16, 16, 17, 16, 16, 10, 9, 8, 8, 8, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 2, 1, 0, 0, 0, 5, 4, 4, 2, 1, 0, 0, 0, 1, 0, 0, 170, 169, 168, 168, 168, 165, 164, 164, 162, 161, 160, 160, 160, 161, 160, 160, 162, 161, 160, 160, 160, 149, 148, 148, 146, 145, 144, 144, 144, 145, 144, 144, 138, 137, 136, 136, 136, 133, 132, 132, 130, 129, 128, 128, 128, 129, 128, 128, 130, 129, 128, 128, 128, 133, 132, 132, 130, 129, 128, 128, 128, 129, 128, 128, 138, 137, 136, 136, 136, 133, 132, 132, 130, 129, 128, 128, 128, 129, 128, 128, 130, 129, 128, 128, 128, 85, 84, 84, 82, 81, 80, 80, 80, 81, 80, 80, 74, 73, 72, 72, 72, 69, 68, 68, 66, 65, 64, 64, 64, 65, 64, 64, 66, 65, 64, 64, 64, 69, 68, 68, 66, 65, 64, 64, 64, 65, 64, 64, 42, 41, 40, 40, 40, 37, 36, 36, 34, 33, 32, 32, 32, 33, 32, 32, 34, 33, 32, 32, 32, 21, 20, 20, 18, 17, 16, 16, 16, 17, 16, 16, 10, 9, 8, 8, 8, 5, 4, 4, 2, 1, 0, };
unsigned char NAF_c_tab[1024] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, };
	
// n must be less than 2^63, takes about 14 nsecs for n~2^30 and about 32 nsecs for n~2^60 (AMD 2.5GHz Athlon64)
void ui_NAF (unsigned long *pbits, unsigned long *nbits, unsigned long n)
{
	register int c, m;
	
	m=n&0x1FF;
	*pbits = NAF_pbits_tab[m];
	*nbits = NAF_nbits_tab[m];
	if ( n < 128 ) return;
	c = NAF_c_tab[m];
	n>>=8;
	m=n&0x1FF;
	*pbits |= (unsigned long)NAF_pbits_tab[m+(c<<9)] << 8;
	*nbits |= (unsigned long)NAF_nbits_tab[m+(c<<9)] << 8;
	if ( n < 128 ) return;
	c = NAF_c_tab[m+(c<<9)];
	n>>=8;
	m=n&0x1FF;
	*pbits |= (unsigned long)NAF_pbits_tab[m+(c<<9)] << 16;
	*nbits |= (unsigned long)NAF_nbits_tab[m+(c<<9)] << 16;
	if ( n < 128 ) return;
	c = NAF_c_tab[m+(c<<9)];
	n>>=8;
	m=n&0x1FF;
	*pbits |= (unsigned long)NAF_pbits_tab[m+(c<<9)] << 24;
	*nbits |= (unsigned long)NAF_nbits_tab[m+(c<<9)] << 24;
	if ( n < 128 ) return;
	c = NAF_c_tab[m+(c<<9)];
	n>>=8;
	m=n&0x1FF;
	*pbits |= (unsigned long)NAF_pbits_tab[m+(c<<9)] << 32;
	*nbits |= (unsigned long)NAF_nbits_tab[m+(c<<9)] << 32;
	if ( n < 128 ) return;
	c = NAF_c_tab[m+(c<<9)];
	n>>=8;
	m=n&0x1FF;
	*pbits |= (unsigned long)NAF_pbits_tab[m+(c<<9)] << 40;
	*nbits |= (unsigned long)NAF_nbits_tab[m+(c<<9)] << 40;
	if ( n < 128 ) return;
	c = NAF_c_tab[m+(c<<9)];
	n>>=8;
	m=n&0x1FF;
	*pbits |= (unsigned long)NAF_pbits_tab[m+(c<<9)] << 48;
	*nbits |= (unsigned long)NAF_nbits_tab[m+(c<<9)] << 48;
	if ( n < 128 ) return;
	c = NAF_c_tab[m+(c<<9)];
	n>>=8;
	m=n&0x1FF;
	*pbits |= (unsigned long)NAF_pbits_tab[m+(c<<9)] << 56;
	*nbits |= (unsigned long)NAF_nbits_tab[m+(c<<9)] << 56;
	return;
}

	
