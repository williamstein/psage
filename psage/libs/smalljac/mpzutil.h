#ifndef _MPZUTIL_INCLUDE_
#define _MPZUTIL_INCLUDE_

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

// This module is a grab-bag of stuff, a lot of which has nothing to do with
// GMP's multi-precision integer arithmetic (mpz).  This module should be split up and refined

// Much of this code is not used by smalljac

#include "asm.h"
#include "gmp.h"

#define MPZ_E						2.7182818284590452353602874713527
#define MPZ_SQRT2					1.41421356237309504880168872421
#define MPZ_LN2						0.69314718055994530941723212145818
#define MPZ_PI						3.1415926535897932384626433832795

#define ULONG_BITS					(8*sizeof(unsigned long))
#define ULONG_BITS_LOG2				(4+(sizeof(unsigned long)/4))		// this only works for 16, 32, and 64 bit word sizes

#define MPZ_SMALL_PRIMES				6542
#define MPZ_MAX_SMALL_PRIME			65521
#define MPZ_MAX_SMALL_INTEGER			0xFFFF
#define MPZ_MAX_SMALL_COMPOSITE		(1<<20)
#define MPZ_MIN_CUBE					0x100000000000UL				// lower bound on cube of a prime > MPZ_MAX_SMALL_PRIME
#define MPZ_MAX_ENUM_PRIME			(0xFFFFFFFF)
#define MPZ_MAX_ENUM_PRIME_LOG2		31							// floor(log2(MPZ_MAX_ENUM_PRIME))

#define MPZ_MAX_PRIMORIAL_W			15
#define MPZ_SMALL_PRIMORIAL_W			6							// P_6 = 2*3*5*7*11*13 = 30030
#define MPZ_SMALL_PRIMORIAL			(2*3*5*7*11*13)				// must be less than MPZ_MAX_SMALL_INTEGER
#define MPZ_SIEVE_LIMIT					(MPZ_MAX_SMALL_INTEGER*MPZ_MAX_SMALL_INTEGER)
#define MPZ_MAX_WHEEL_W				9							// needs 200 mb map to construct a wheel this large

#define MPZ_MAX_TREE_FACTORS			100
#define MPZ_MAX_INVERTS				200
#define MPZ_MAX_SMALL_FACTORS			256

#define MPZ_MAX_UI_PP_FACTORS			20
#define MAX_UI_PP_FACTORS				MPZ_MAX_UI_PP_FACTORS		// more sensible name
#define MPZ_MAX_FACTORS				64							// we could be a bit more generous here

#define log2(x)		(log(x)/MPZ_LN2)
static inline int mod (int n, int m) { n %= m;  return ( n<0 ? (m + n) : n ); }
static inline int modl (long n, long m) { n %= m;  return ( n<0 ? (m + n) : n ); }

struct wheel_struct {
	unsigned long n;
	unsigned long phi;
	unsigned char *gaps;
	unsigned maxgap;
};
typedef struct wheel_struct wheel_t;

struct prime_enum_ctx_struct {
	unsigned L;
	unsigned r[MPZ_MAX_ENUM_PRIME_LOG2+1];
	unsigned h;
	unsigned w;
	unsigned *wi;
	unsigned *wv;
	unsigned p;
	unsigned pi;
	unsigned j;
	unsigned start;
	wheel_t *wheel;
	unsigned char *map;
	unsigned char *mp;
	int powers;
};
typedef struct prime_enum_ctx_struct prime_enum_ctx_t;

struct factor_tree_item_struct {
	int w;
	mpz_t factors[MPZ_MAX_TREE_FACTORS];						// list of prime power factors
	struct factor_tree_item_struct *subtrees[MPZ_MAX_TREE_FACTORS];	// subtrees factor p_i - 1 where p_i is the base of the ith pp factor
};

typedef struct factor_tree_item_struct *mpz_ftree_t;

#define _ui_ceil_ratio(a,b)		(((a)%(b)) ? (a)/(b)+1 : (a)/(b))
#define _ui_max(a,b)			((a)>(b)?(a):(b))
#define _ui_min(a,b)			((a)>(b)?(b):(a))

#define i_sgn_c(x)		((x)<0?'-':'+')
#define i_abs(x)		((x)<0?-(x):(x))

#define mpz_get_i(z)			(((long)mpz_get_ui(z))*mpz_sgn(z))
#define mylround(x)			((long)floor((x)+0.5))				// assumes x is positive - workaround problem with builtin lround function

void mpz_util_init();

int mpz_eval_expr (mpz_t o, char *expr);
void mpz_mulm (mpz_t o, mpz_t a, mpz_t b, mpz_t m);
void mpz_powm_tiny (mpz_t o, mpz_t b, unsigned e, mpz_t m);
void mpz_powm_big (mpz_t o, mpz_t b, mpz_t e, mpz_t m);
void mpz_addm (mpz_t o, mpz_t a, mpz_t b, mpz_t m);
void mpz_subm (mpz_t o, mpz_t a, mpz_t b, mpz_t m);
void mpz_subm_ui (mpz_t o, mpz_t a, unsigned long b, mpz_t m);
void mpz_negm (mpz_t a, mpz_t m);
void mpz_set_i (mpz_t o, long n);
void mpq_set_i (mpq_t o, long n);
void mpz_parallel_invert (mpz_t o[], mpz_t a[], unsigned n, mpz_t p);
int mpz_eval_term_ui (mpz_t o, unsigned long numvars, unsigned long vars[], unsigned long exps[]);
char *ui_term_to_string (char *buf, unsigned long numvars, unsigned long vars[], unsigned long exps[]);

int mpz_sqrt_modprime (mpz_t o, mpz_t a, mpz_t p);
int mpz_factor_small (unsigned long p[], unsigned long h[], mpz_t bigp, mpz_t n, int max_factors, int max_hard_bits);
int mpz_remove_smallsquares (mpz_t o, mpz_t n);
int mpz_coarse_part (mpz_t o, mpz_t x, mpz_t y);
void mpz_mul_set (mpz_t o, mpz_t *a, unsigned long n);		// note overwrites array specified by a
int mpz_prime_mult (mpz_t o, mpz_t p, mpz_t maxp, unsigned long maxbits);
void mpz_power_primorial (mpz_t o, unsigned long L);
int mpz_compatible (mpz_t a, mpz_t b);

wheel_t *primorial_wheel_alloc (int w);
void primorial_wheel_free (wheel_t *pwheel);

prime_enum_ctx_t *fast_prime_enum_start (unsigned start, unsigned L, int powers);
unsigned fast_prime_enum (prime_enum_ctx_t *ctx);
unsigned fast_prime_power_enum (prime_enum_ctx_t *ctx);		// enumerates floor_p(L) by p
void fast_prime_enum_end (prime_enum_ctx_t *ctx);

unsigned long ui_randomm(unsigned long m);
unsigned long ui_randomb (unsigned long b);

unsigned long mpz_randomf (mpz_t o, mpz_t N, mpz_t factors[], unsigned long w);
unsigned long mpz_randomft (mpz_t o, mpz_t N,  mpz_t factors[], mpz_ftree_t factor_trees[], unsigned long w);
mpz_ftree_t mpz_randomfp (mpz_t o, mpz_t N);
mpz_ftree_t mpz_randomfpp (mpz_t o, mpz_t N);
void mpz_clear_ftree (mpz_ftree_t t);

void mpz_randompp (mpz_t Q, mpz_t N);
unsigned long mpz_pp_base (mpz_t b, mpz_t q);
unsigned long ui_pp_base (unsigned long pp);
unsigned long ui_pp_div (unsigned long n, unsigned long p);

unsigned long ui_pdiff (unsigned long a, unsigned long b);		// returns least p s.t. the power of p dividing a and b differ
int ui_factor (unsigned long p[], unsigned long h[], unsigned long n);
int ui_divisor_in_interval (unsigned long p[], unsigned long h[], int w, unsigned long min, unsigned long max);

int mpz_remove_small_primes (mpz_t o, mpz_t n, unsigned long exps[], unsigned long maxprimes);		// returns # distinct primes removed
unsigned long mpz_nearprime (mpz_t n, unsigned long L);

int ui_is_small_prime (unsigned long p);
int ui_is_prime (unsigned long p);
unsigned long ui_small_prime (unsigned long n);				// returns the nth prime for 0 < n <= MPZ_SMALL_PRIMES (2 is the 1st prime)
unsigned long ui_small_prime_index (unsigned long p);			// returns pi(n) for n <= MPZ_MAX_SMALL_INTEGER
unsigned long ui_primorial (int w);							// returns the P_w = the w_th primorial or 0 if P_w > ULONG_MAX
unsigned long ui_primorial_phi (int w);							// returns the phi(P_w) or 0 if P_w > ULONG_MAX
unsigned long ui_pi_bound (unsigned long n);					// returns a guaranteed upper bound on pi(n)
unsigned long ui_phi (unsigned long n);						// returns \phi(n) = |Z/nZ|^*
int ui_sqrt(long n);											// returns -1 for non square input

unsigned long mpz_get_bits_ui (mpz_t a, unsigned long i, unsigned long n);		// gets bits i thru i+n-1 and returns as ui
unsigned long mpz_set_bits_ui (mpz_t a, unsigned long i, unsigned long n, unsigned long bits);		// sets bits i thru i+n-1 and returns as ui
double mpz_log2 (mpz_t a);		// approximation guarunteed to be between floor(lg(a)) and lg(a)

void mpz_reset_counters ();
void mpz_report_counters ();

void mpz_sort (mpz_t a[], unsigned long n);

unsigned long ui_inverse (unsigned long a, unsigned long m);		// returns the multiplicative inverse of a mod m
unsigned long ui_gcd (unsigned long a, unsigned long b);
unsigned long ui_gcd_ext (unsigned long a, unsigned long b, long *x, long *y);
unsigned long ui_crt (unsigned long a, unsigned long M, unsigned long b, unsigned long N);	// returns x cong a%M and cong b%N with 0<=x<M*N, assumes gcd(M,N)=1
unsigned long bach_gcd (long a, long b, unsigned long N);
int ui_compatible  (unsigned long a, unsigned long b);

int ui_qsort_cmp (const void *a, const void *b);
unsigned long ui_wt (unsigned long x);								// number of bits in binary rep of x
unsigned long ui_lg_ceil (unsigned long x);							// ceil(lg(x))
#define ui_lg_floor(x)	_asm_highbit(x)								// floor(lg(x))
unsigned long ui_binexp_cost (unsigned long x);
unsigned long ui_get_bits (unsigned long x, unsigned long pos, unsigned long bits);
char *ui_bit_string (char *buf, unsigned long x);

unsigned long mpz_store_bytes (mpz_t a);
void mpz_store (char *s, mpz_t a);
void mpz_retrieve (mpz_t o, char *s);

// counts multiples of k in [Min,Max] - assumes that the least multiple of k greater than Max fits in an unsigned long (this will be true if k and Max are both less than U
static inline unsigned long ui_multiples_in_range (unsigned long k, unsigned long Min, unsigned long Max)
{
	register unsigned long a, b;
	
	a = _ui_ceil_ratio (Min, k)*k;			// least multiple of k >= Min
	b = (Max/k)*k;						// greatest multiple of k <= Max
	if ( b < a ) return 0;				// we could avoid this using signed arithmetic (and probably should, but we could have overflow problems if Max won't fit in a long)
	return (b-a)/k + 1;
}

static inline unsigned long ui_trim (unsigned long x) { while ( x & 0x1 ) x << 1; return x; }
static inline unsigned ui_len (unsigned long x) { register unsigned i;   for ( i = 0 ; x ; i++ ) x >>= 1;  return i; }
unsigned long ui_eval_expr (char *expr);

typedef struct NAF_entry_struct {
	char digit;
	char c;
} NAF_entry_t;

extern NAF_entry_t ui_NAF_table[2][4];

void ui_NAF (unsigned long *pbits, unsigned long *nbits, unsigned long n);

static inline int ui_NAF_byte (unsigned char *pbits, unsigned char *nbits, unsigned long n, int c)
{
	register int i,j,k;
	register unsigned char pb, nb;
	
	pb=nb=0;
	for ( i = 0 ; i < 8 ; i++ ) { j=n&3; k=ui_NAF_table[c][j].digit; if ( k>0 ) pb|=(1<<i); if ( k < 0 ) nb|=(1<<i); c = ui_NAF_table[c][j].c; n>>=1; }
	*pbits = pb; *nbits=nb;
	return c;
}

static inline long atol_expr(char *expr)
{
	register char *s;
	long i,b,n,N;
	
	for ( s = expr ; *s && *s != 'e' && *s != 'E' ; s++ );
	if ( *s ) { *s++; N = b = atol (expr); n = atol(s); for ( i = 1; i < n ; i++ )  N *= b; } else N = atol(expr);
	return N;
}

#endif
