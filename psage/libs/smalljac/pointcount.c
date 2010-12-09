#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include "ff.h"			// only used by slowcount
#include "cstd.h"
#include "pointcount.h"

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
    Fast pointcounting using finite differences as described in [KedlayaSutherland2007].
    Hardwired versions for each genus are provided for optimal performance.
    Characteristic 2 is not supported (for genus 1 p=2 is handled in smalljac.c).
*/


/*   We use the triangle T(j,k) are from OEIS A019538 - and add row 0 and col 0 which are just 1 0 0 0 ...
     Note that the square tables displayed by OEIS transpose the entries
     
      T(j,k) = k!S(j,k)    where S(j,k) is a Stirling number of the second kind
*/

static unsigned long T[10][10] = { {1,	0,	0,		0,		0,		0,		0,		0,		0,		0},
							   {0,	1,	0,		0,		0,		0,		0,		0,		0,		0},
							   {0,	1,	2,		0,		0,		0,		0,	    	0,		0,		0},
							   {0,	1,	6,		6,		0,		0,		0,	   	0,		0,		0},
							   {0,	1,	14,		36,		24,		0,		0,  		0,		0,		0},
							   {0,	1,	30,		150,		240,		120,		0,		0,		0,		0},
							   {0,	1,	62,		540,		1560,	1800,	720,		0,		0,		0},
							   {0,	1,	126,		1806,	8400,	16800,	15120,	5040,	0,		0},
							   {0,	1,	254,		5796,	40824,	126000,	191520,	141120,	40320,	0},
							   {0,	1,	510,		18150,	186480,	834120,	1905120,	2328480,	1451520,	362880}};

static unsigned long *map;

static unsigned long tab[64] = { 0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x80,
					      0x100, 0x200, 0x400, 0x800, 0x1000, 0x2000, 0x4000, 0x8000,
					      0x10000, 0x20000, 0x40000, 0x80000, 0x100000, 0x200000, 0x400000, 0x800000,
					      0x1000000, 0x2000000, 0x4000000, 0x8000000, 0x10000000, 0x20000000, 0x40000000, 0x80000000,
					      0x100000000, 0x200000000, 0x400000000, 0x800000000, 0x1000000000, 0x2000000000, 0x4000000000, 0x8000000000,
					      0x10000000000, 0x20000000000, 0x40000000000, 0x80000000000, 0x100000000000, 0x200000000000, 0x400000000000, 0x800000000000,
					      0x1000000000000, 0x2000000000000, 0x4000000000000, 0x8000000000000, 0x10000000000000, 0x20000000000000, 0x40000000000000, 0x80000000000000,
					      0x100000000000000, 0x200000000000000, 0x400000000000000, 0x800000000000000, 0x1000000000000000, 0x2000000000000000, 0x4000000000000000, 0x8000000000000000 };

#define _advance_3t_seq()			t0 -= t1; if ( t0 < 0 ) t0 += p;  t1 -= t2;  if ( t1 < 0 ) t1 += p; t2 -= t3; if ( t2 < 0 ) t2 += p; 
#define _advance_4t_seq()			t0 -= t1; if ( t0 < 0 ) t0 += p;  t1 -= t2;  if ( t1 < 0 ) t1 += p; t2 -= t3; if ( t2 < 0 ) t2 += p; t3 -= t4; if ( t3 < 0 ) t3 += p; 
#define _advance_5t_seq()			t0 -= t1; if ( t0 < 0 ) t0 += p;  t1 -= t2;  if ( t1 < 0 ) t1 += p; t2 -= t3; if ( t2 < 0 ) t2 += p; t3 -= t4; if ( t3 < 0 ) t3 += p; t4 -= t5; if ( t4 < 0 ) t4 += p;
#define _advance_6t_seq()			t0 -= t1; if ( t0 < 0 ) t0 += p;  t1 -= t2;  if ( t1 < 0 ) t1 += p; t2 -= t3; if ( t2 < 0 ) t2 += p; t3 -= t4; if ( t3 < 0 ) t3 += p; t4 -= t5; if ( t4 < 0 ) t4 += p; t5 -= t6; if ( t5 < 0 ) t5 += p; 
#define _advance_7t_seq()			t0 -= t1; if ( t0 < 0 ) t0 += p;  t1 -= t2;  if ( t1 < 0 ) t1 += p; t2 -= t3; if ( t2 < 0 ) t2 += p; t3 -= t4; if ( t3 < 0 ) t3 += p; t4 -= t5; if ( t4 < 0 ) t4 += p; \
								t5 -= t6; if ( t5 < 0 ) t5 += p; t6 -= t7; if ( t6 < 0 ) t6 += p; 
#define _advance_8t_seq()			t0 -= t1; if ( t0 < 0 ) t0 += p;  t1 -= t2;  if ( t1 < 0 ) t1 += p; t2 -= t3; if ( t2 < 0 ) t2 += p; t3 -= t4; if ( t3 < 0 ) t3 += p; t4 -= t5; if ( t4 < 0 ) t4 += p; \
								t5 -= t6; if ( t5 < 0 ) t5 += p; t6 -= t7; if ( t6 < 0 ) t6 += p; t7 -= t8; if ( t7 < 0 ) t7 += p;   
#define _advance_9t_seq()			t0 -= t1; if ( t0 < 0 ) t0 += p;  t1 -= t2;  if ( t1 < 0 ) t1 += p; t2 -= t3; if ( t2 < 0 ) t2 += p; t3 -= t4; if ( t3 < 0 ) t3 += p; t4 -= t5; if ( t4 < 0 ) t4 += p; \
								t5 -= t6; if ( t5 < 0 ) t5 += p; t6 -= t7; if ( t6 < 0 ) t6 += p;   t7 -= t8; if ( t7 < 0 ) t7 += p;   t8 -= t9; if ( t8 < 0 ) t8 += p; 

void pointcount_init (unsigned maxp)
{
	map = mem_alloc (maxp/8+64);			// budget extra space for wrapping
}

/*
	As described in KedlayaSutherland2007, we compute
		D[k] = (-1)^k\(Delta^k f)(0)
	where Delta is the difference operator: (Delta f)(x) = f(x+1) - f(x)
*/
void pointcount_precompute (mpz_t D[], mpz_t f[], int degree)
{
	static mpz_t x;
	static int init;
	int j, k;
	
	if ( ! init ) { mpz_init (x);  init = 1; }
	for ( k= 0 ; k <= degree ; k++ ) {
		mpz_set_ui (D[k], 0);
		for ( j = 0 ; j <= degree ; j++ ) {
			if ( mpz_sgn (f[j]) ) {
				mpz_mul_ui (x, f[j], T[j][k]);
				mpz_add (D[k], D[k], x);
			}
		}
		if ( k&1 ) mpz_neg (D[k], D[k]);
	}
}


void pointcount_precompute_long (long D[], long f[], int degree)
{
	int j, k;
	
	for ( k= 0 ; k <= degree ; k++ ) {
		D[k] = 0;
		for ( j = 0 ; j <= degree ; j++ ) D[k] += f[j] * T[j][k];
		if ( k&1 ) D[k] = -D[k];
	}
}


// Creates a bitmap of non-zero quadratic residues in F_p
void pointcount_map_residues (unsigned p)
{
	register signed t0, t1;
	unsigned long x;

	memset (map, 0, p/8+9);				// be sure to clear out 64 bits past the end
	t1 = p;
	x = (unsigned long)((p-1)/2);
	x *= x;
	t0 = (signed)(x%(unsigned long)p);		// the first square we check is [(p-1)/2]^2 - need to compute in long to handle overflow
										// work down toward zero but don't include zero since its handled explicitly
	while ( t0 ) {
		map[t0>>6] |= tab[t0&0x3F];
		t1 -= 2;
		t0 -= t1;  if ( t0 < 0 ) t0 += p;
	};
}	


// Creates a bitmap of non-zero quadratic residues in F_p which are <= (p-1)/2 (as integers mod p)
void pointcount_map_small_residues (unsigned p)
{
	register signed c, t0, t1;
	unsigned long x;

	memset (map, 0, p/16+9);				// be sure to clear out 64 bits past the end
	
	t1 = p;
	x = (unsigned long)((p-1)/2);
	x *= x;
	t0 = (signed)(x%(unsigned long)p);		// the first square we check is [(p-1)/2]^2 - need to compute in long to handle overflow
										// work down toward zero but don't include zero since its handled explicitly
	c = p>>1;
	while ( t0 ) {
		if ( t0 <= c ) map[t0>>6] |= tab[t0&0x3F];
		t1 -= 2;
		t0 -= t1;  if ( t0 < 0 ) t0 += p;
	};
}	


// Create a bitmap of non-zero quadratic residues in F_p^2 represented as az+b F[z]/(z^2-s) where s is the least
// non-residue mod p.  Note that z is not necessarily a primitive element of F_p^2 and we don't need it to be.
// The bit p*a+b is set to 1 if az+b is a residue.  The value of s is returned.

unsigned pointcount_map_residues_d2 (unsigned p)
{
	register signed i, k, s, t00, t01, t10, t11;
	unsigned long x;

	memset (map, 0, (p*p)/8+1);
	
	// This program is not as efficient as it could be, it maps every residue twice, but
	// efficiency is not so critical in extension fields and simplicity is a virtue
	// Map residues mod p to compute least non-residue s - this is overkill of course.
	pointcount_map_residues (p);
	for ( s = 2 ; (map[s>>6] & tab[s&0x3f]) ; s++ );
		
	// For simplicity, we use addition rather than subtraction
	for ( k = 1 ; k <= (p-1)/2 ; k++ ) {
		x = (unsigned long)k;  x *= k;  x *= s;  x %= p;					// x = (kz)^2 = k^2s
		t00 = (signed)x;  t01 = 0;									// t0 = 0z + k^2s
		t10 = 1;  t11 = k+k;  if ( t11 >= p ) t11 -= p;						// t1 = 2kz + 1
		do {
			i = p*t01+t00;
			map[i>>6] |= tab[i&0x3F];
			t00 += t10;  if ( t00 >= p ) t00 -= p;
			t01 += t11; if ( t01 >= p ) t01 -= p;
			t10 += 2; if ( t10 >= p ) t10 -= p;
		} while ( t01 );
	}
	map[0] &= ~(1UL);			// clear zero-bit, we only want non-zero residues
	return s;
}	

#define POINTCOUNT_MAX_D3_PRIME		47
static int d3stab[POINTCOUNT_MAX_D3_PRIME+1] = { 0,0,0,1,0,2,0,2,0,0,0,3,0,1,0,0,0,2,0,4,0,0,0,4,0,0,0,0,0,1,0,1,0,0,0,0,0,2,0,0,0,1,0,2,0,0,0,1};

// Create a bitmap of non-zero quadratic residues in F_p^3 represented as az^2+bz+c in F[z]/(z^3-z-s) where s the
// the least s for which z^3-z-s is irreducible mod p.  Note that z is not necessarily a primitive element of F_p^3
// and we don't need it to be.  The bit p^2*a+pb+c is set to 1 if az^2+bz+c is a residue.
// The value of s is returned - currently hardwired via tablelookup for odd primes up to 23
unsigned pointcount_map_residues_d3 (unsigned p)
{
	register signed i, j, a, b, s, t00, t01, t02, t10, t11, t12;

	if ( p > POINTCOUNT_MAX_D3_PRIME || ! (s=d3stab[p]) ) { err_printf ("invalid prime %u specified in pointcount_map_residues_d3\n", p);  exit (0); }
	memset (map, 0, (p*p)/8+1);	
		
	// For simplicity, we use addition rather than subtraction, assume p is small enough so that overflow is not a worry
	for ( a = 0 ; a <= (p-1)/2 ; a++ ) {									// only need to square half the elements to get all the residues
		for ( b = 0 ; b < p ; b++ ) {
			i=2*a*b;  j=a*a;
			t00 = (i*s)%p;  t01 = (j*s+i)%p;  t02 = (j+b*b)%p;				// t0 = (az^2+bz)^2 = (a^2+b^2)z^2 + (a^2s+2ab)z + 2abs  (applying z^3=z+s)
			t10 = 1;  t11 = (b+b)%p;  t12 = (a+a)%p;						// t1 = 2(az^2+bz)+1 = 2az^2 + 2bz + 1
			for ( j = 0  ; j < p ; j++ ) {
				i = p*p*t02+p*t01+t00;
				map[i>>6] |= tab[i&0x3F];
				t00 += t10;  if ( t00 >= p ) t00 -= p;
				t01 += t11; if ( t01 >= p ) t01 -= p;
				t02 += t12; if ( t02 >= p ) t02 -= p;
				t10 += 2; if ( t10 >= p ) t10 -= p;
			}
		}
	}
	map[0] &= ~(1UL);			// clear zero-bit, we only want non-zero residues
	return s;
}	


// Creates a bitmap of non-zero cubic residues in F_p
void pointcount_map_cubic_residues (unsigned p)
{
	register signed t0, t1, t2, t3;
	
	memset (map, 0, p/8+1);
	// enumerate non-zero values of f(x) = x^3 via finite differences starting at f(1)
	t0 = 1;					// f(1)
	t1 = 7%p;  if ( t1) t1 = p-t1;	// (-1)^1 (Delta^1 f)(1)
	t2 = 12%p;				// (-1)^2 (Delta^2 f)(1)
	t3 = 6%p;  if ( t3 ) t3 = p-t3;	// (-1)^3 (Delta^3 f)(1)

	while ( t0 ) {
		map[t0>>6] |= tab[t0&0x3F];
		_advance_3t_seq();
	};
}

// counts points on Picard curves y^3 = f(x) with f(x) of degree 4
unsigned pointcount_pd4 (unsigned long D[5], unsigned p)
{
	register signed  i, t0, t1, t2, t3, t4, m, n;
	
	pointcount_map_cubic_residues (p);
	t0 = (signed) D[0];  t1 = (signed) D[1];  t2 = (signed) D[2];  t3 = (signed) D[3];  t4 = (signed) D[4];
	m = n = 0;
	for ( i = 0 ; i < p ; i++ ) {
		if ( ! t0 ) m++;
		if ( (map[t0>>6] & tab[t0&0x3F]) ) n++;
		_advance_4t_seq();
	}
	return m+((p%3)==1?3:1)*n+1;		// add point at infinity
}



unsigned pointcount_g1 (unsigned long D[4], unsigned p)
{
	register signed  i, t0, t1, t2, t3, t, s0, m, n;

	pointcount_map_residues (p);
	t0 = (signed) D[0];  t1 = (signed) D[1];  t2 = (signed) D[2];  t3 = (signed) D[3];
	m = n = 0;
	for ( i = 0 ; i < p ; i++ ) {
		if ( ! t0 ) m++;
		if ( (map[t0>>6] & tab[t0&0x3F]) ) n++;
		_advance_3t_seq();
	}
	return m+2*n+1;		// add point at infinity
}


unsigned pointcount_g2 (unsigned long D[6], unsigned p)
{
	register signed  i, t0, t1, t2, t3, t4, t5, m, n;
	
	pointcount_map_residues (p);
	t0 = (signed) D[0];  t1 = (signed) D[1];  t2 = (signed) D[2];  t3 = (signed) D[3];  t4 = (signed) D[4];  t5 = (signed) D[5];
	m = n = 0;
	for ( i = 0 ; i < p ; i++ ) {
		if ( ! t0 ) m++;
		if ( (map[t0>>6] & tab[t0&0x3F]) ) n++;
		_advance_5t_seq();
	}
	return m+2*n+1;		// add point at infinity
}

unsigned pointcount_g2d6 (unsigned long D[7], unsigned p, unsigned long f6)
{
	register signed  i, t0, t1, t2, t3, t4, t5, t6, m, n;
	
	pointcount_map_residues (p);
	t0 = (signed) D[0];  t1 = (signed) D[1];  t2 = (signed) D[2];  t3 = (signed) D[3];  t4 = (signed) D[4];  t5 = (signed) D[5];  t6 = (signed) D[6];
	m = 1 + ((map[f6>>6] & tab[f6&0x3f])?1:-1);		// 2 or 0 points at infinity depending on whether leading coeff is square or not.
	n = 0;
	for ( i = 0 ; i < p ; i++ ) {
		if ( ! t0 ) m++;
		if ( (map[t0>>6] & tab[t0&0x3F]) ) n++;
		_advance_6t_seq();
	}
	return m+2*n;
}

unsigned pointcount_g3 (unsigned long D[8], unsigned p)
{
	register signed  i, t0, t1, t2, t3, t4, t5, t6, t7, m, n;
	unsigned long x;
	
	pointcount_map_residues (p);
	t0 = (signed) D[0];  t1 = (signed) D[1];  t2 = (signed) D[2];  t3 = (signed) D[3];  t4 = (signed) D[4];  t5 = (signed) D[5];  t6 = (signed) D[6];  t7 = (signed) D[7];
	m = n = 0;
	for ( i = 0 ; i < p ; i++ ) {
		if ( ! t0 ) m++;
		if ( (map[t0>>6] & tab[t0&0x3F]) ) n++;
		_advance_7t_seq();
	}
	return m+2*n+1;		// add point at infinity
}

unsigned pointcount_g3d8 (unsigned long D[9], unsigned p, unsigned long f8)
{
	register signed  i, t0, t1, t2, t3, t4, t5, t6, t7, t8, m, n;
	unsigned long x;
	
	pointcount_map_residues (p);
	t0 = (signed) D[0];  t1 = (signed) D[1];  t2 = (signed) D[2];  t3 = (signed) D[3];  t4 = (signed) D[4];  t5 = (signed) D[5];  t6 = (signed) D[6];  t7 = (signed) D[7];  t8 = (signed) D[8];
	m = 1 + ((map[f8>>6] & tab[f8&0x3F])?1:-1);		// 2 or 0 points at infinity depending on whether leading coeff is square or not.
	n = 0;
	for ( i = 0 ; i < p ; i++ ) {
		if ( ! t0 ) m++;
		if ( (map[t0>>6] & tab[t0&0x3F]) ) n++;
		_advance_8t_seq();
	}
	return m+2*n;
}

unsigned pointcount_g4 (unsigned long D[10], unsigned p)
{
	register signed  i, t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, m, n;
	unsigned long x;
	
	pointcount_map_residues(p);
	t0 = (signed) D[0];  t1 = (signed) D[1];  t2 = (signed) D[2];  t3 = (signed) D[3];  t4 = (signed) D[4];  t5 = (signed) D[5];
	t6 = (signed) D[6];  t7 = (signed) D[7];  t8 = (signed) D[8];  t9 = (signed) D[9];
	m = n = 0;
	for ( i = 0 ; i < p ; i++ ) {
		if ( ! t0 ) m++;
		if ( (map[t0>>6] & tab[t0&0x3F]) ) n++;
		_advance_9t_seq();
	}
	return m+2*n+1;		// add point at infinity
}

/*
	The "big" versions below are intended for use with "big" primes.  They use a "small" residue
	map (half the size) that may better fit in cache.
*/

unsigned pointcount_big_g1 (unsigned long D[4], unsigned p)
{
	register signed  i, j, k, t0, t1, t2, t3, t, c, m, n;

	pointcount_map_small_residues (p);
	t0 = (signed) D[0];  t1 = (signed) D[1];  t2 = (signed) D[2];  t3 = (signed) D[3];
	
	m = n = 0;
	c = p>>1;
	if ( (p&0x3) == 1 ) {							// -1 a quadratic residue mod p?
		for ( i = 0 ; i < p ; i++ ) {
			if ( ! t0 ) m++;
			j = ( t0>c ? p-t0 : t0 );
			if ( (map[j>>6] & tab[j&0x3F]) ) n++;
			_advance_3t_seq();
		}
	} else {										// this case is unfortunately nearly twice as slow - the conditional code kills it
		for ( i = 0 ; i < p ; i++ ) {
			if ( ! t0 ) m++;
			k = (t0>c);
			j = ( k ? p-t0 : t0 );
			n += ( (map[j>>6] & tab[j&0x3F]) ? !k : k );
			_advance_3t_seq();
		}
	}
	return m+2*n+1;		// add point at infinity
}

unsigned pointcount_big_g2 (unsigned long D[6], unsigned p)
{
	register signed  i, j, k, t0, t1, t2, t3, t4, t5, t, m, n, c;
	
	pointcount_map_small_residues (p);
	t0 = (signed) D[0];  t1 = (signed) D[1];  t2 = (signed) D[2];  t3 = (signed) D[3];  t4 = (signed) D[4];  t5 = (signed) D[5];

	m = n = 0;
	c = p>>1;

	if ( (p&0x3) == 1 ) {							// -1 a quadratic residue mod p?
		for ( i = 0 ; i < p ; i++ ) {
			if ( ! t0 ) m++;
			j = ( t0>c ? p-t0 : t0 );
			if ( (map[j>>6] & tab[j&0x3F]) ) n++;
			_advance_5t_seq();
		}
	} else {										// this case is unfortunately nearly twice as slow - the conditional code kills it
		for ( i = 0 ; i < p ; i++ ) {
			if ( ! t0 ) m++;
			k = (t0>c);
			j = ( k ? p-t0 : t0 );
			n += ( (map[j>>6] & tab[j&0x3F]) ? !k : k );
			_advance_5t_seq();
		}
	}
	return m+2*n+1;		// add point at infinity
}

unsigned pointcount_big_g2d6 (unsigned long D[7], unsigned p, unsigned long f6)
{
	register signed  i, j, k, t0, t1, t2, t3, t4, t5, t6, t, m, n, c;
	
	pointcount_map_small_residues (p);
	t0 = (signed) D[0];  t1 = (signed) D[1];  t2 = (signed) D[2];  t3 = (signed) D[3];  t4 = (signed) D[4];  t5 = (signed) D[5];  t6 = (signed) D[6];

	m = 1 + ((map[f6>>6] & tab[f6&0x3f])?1:-1);		// 2 or 0 points at infinity depending on whether leading coeff is square or not.
	n = 0;
	c = p>>1;

	if ( (p&0x3) == 1 ) {							// -1 a quadratic residue mod p?
		for ( i = 0 ; i < p ; i++ ) {
			if ( ! t0 ) m++;
			j = ( t0>c ? p-t0 : t0 );
			if ( (map[j>>6] & tab[j&0x3F]) ) n++;
			_advance_6t_seq();
		}
	} else {										// this case is unfortunately nearly twice as slow - the conditional code kills it
		for ( i = 0 ; i < p ; i++ ) {
			if ( ! t0 ) m++;
			k = (t0>c);
			j = ( k ? p-t0 : t0 );
			n += ( (map[j>>6] & tab[j&0x3F]) ? !k : k );
			_advance_6t_seq();
		}
	}
	return m+2*n;
}

unsigned pointcount_big_g3 (unsigned long D[8], unsigned p)
{
	register signed  i, j, k, t0, t1, t2, t3, t4, t5, t6, t7, m, n, c;
	unsigned long x;
	
	pointcount_map_small_residues (p);
	t0 = (signed) D[0];  t1 = (signed) D[1];  t2 = (signed) D[2];  t3 = (signed) D[3];  t4 = (signed) D[4];  t5 = (signed) D[5];  t6 = (signed) D[6];  t7 = (signed) D[7];
	
	c = p>>1;
	m = n = 0;
	if ( (p&0x3) == 1 ) {							// -1 a quadratic residue mod p?
		for ( i = 0 ; i < p  ; i++ ) {
			if ( ! t0 ) m++;
			j = ( t0>c ? p-t0 : t0 );
			if ( (map[j>>6] & tab[j&0x3F]) ) n++;
			_advance_7t_seq();
		}
	} else {										// this case is unfortunately nearly twice as slow - the conditional code kills it
		for ( i = 0 ; i < p ; i++ ) {
			if ( ! t0 ) m++;
			k = (t0>c);
			j = ( k ? p-t0 : t0 );
			n += ( (map[j>>6] & tab[j&0x3F]) ? !k : k );
			_advance_7t_seq();
		}
	}
	return m+2*n+1;		// add point at infinity
}

unsigned pointcount_big_g4 (unsigned long D[10], unsigned p)
{
	register signed  i, j, k, t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, m, n, c;
	unsigned long x;
	
	pointcount_map_small_residues (p);
	t0 = (signed) D[0];  t1 = (signed) D[1];  t2 = (signed) D[2];  t3 = (signed) D[3];  t4 = (signed) D[4];  t5 = (signed) D[5];
	t6 = (signed) D[6];  t7 = (signed) D[7];  t8 = (signed) D[8];  t9 = (signed) D[9];
	
	c = p>>1;
	m = n = 0;
	if ( (p&0x3) == 1 ) {							// -1 a quadratic residue mod p?
		for ( i = 0 ; i < p  ; i++ ) {
			if ( ! t0 ) m++;
			j = ( t0>c ? p-t0 : t0 );
			if ( (map[j>>6] & tab[j&0x3F]) ) n++;
			_advance_9t_seq();
		}
	} else {										// this case is unfortunately nearly twice as slow - the conditional code kills it
		for ( i = 0 ; i < p ; i++ ) {
			if ( ! t0 ) m++;
			k = (t0>c);
			j = ( k ? p-t0 : t0 );
			n += ( (map[j>>6] & tab[j&0x3F]) ? !k : k );
			_advance_9t_seq();
		}
	}
	return m+2*n+1;		// add point at infinity
}


#define _count32(z)		if ( (z & 1) ) p0 += 2;   z >>= 1;   if ( (z & 1) ) p1 += 2;   z >>= 1;   if ( (z & 1) ) p2 += 2;   z >>= 1;  if ( (z & 1) ) p3 += 2;  z >>= 1; \
						if ( (z & 1) ) p4 += 2;  z >>= 1;   if ( (z & 1) ) p5 += 2;   z >>= 1;   if ( (z & 1) ) p6 += 2;   z >>= 1;   if ( (z & 1) ) p7 += 2;   z >>= 1; \
						if ( (z & 1) ) p8 += 2;  z >>= 1;  if ( (z & 1) ) p9 += 2;  z >>= 1;  if ( (z & 1) ) p10 += 2;  z >>= 1;  if ( (z & 1) ) p11 += 2;   z >>= 1; \
						if ( (z & 1) ) p12 += 2;  z >>= 1;  if ( (z & 1) ) p13 += 2;  z >>= 1;  if ( (z & 1) ) p14 += 2;  z >>= 1;  if ( (z & 1) ) p15 += 2;  z >>= 1; \
						if ( (z & 1) ) q0 += 2;   z >>= 1;   if ( (z & 1) ) q1 += 2;   z >>= 1;   if ( (z & 1) ) q2 += 2;   z >>= 1;  if ( (z & 1) ) q3 += 2;  z >>= 1; \
						if ( (z & 1) ) q4 += 2;  z >>= 1;   if ( (z & 1) ) q5 += 2;   z >>= 1;   if ( (z & 1) ) q6 += 2;   z >>= 1;   if ( (z & 1) ) q7 += 2;   z >>= 1; \
						if ( (z & 1) ) q8 += 2;  z >>= 1;  if ( (z & 1) ) q9 += 2;  z >>= 1;  if ( (z & 1) ) q10 += 2;  z >>= 1;  if ( (z & 1) ) q11 += 2;   z >>= 1; \
						if ( (z & 1) ) q12 += 2;  z >>= 1;  if ( (z & 1) ) q13 += 2;  z >>= 1;  if ( (z & 1) ) q14 += 2;  z >>= 1;  if ( (z & 1) ) q15 += 2;

// hardwired for 32x
int pointcount_multi_g2 (unsigned pts[], unsigned long D[6], unsigned p)
{
	register signed  i, t0, t1, t2, t3, t4, t5;
	unsigned p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15;
	unsigned q0, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, q14, q15;
	register unsigned long z;

	if ( p < 32 ) { err_printf ("cannot multicount with p < 32\n");  exit (0); }
		
	pointcount_map_residues (p);
	// wrap table around top
	for ( t2 = 1 ; t2 < 32 ; t2++ )
		if ( map[t2>>6] & tab[t2&0x3f] ) map[(p+t2)>>6] |= tab[(p+t2)&0x3f];
		
	t0 = (signed) D[0];  t1 = (signed) D[1];  t2 = (signed) D[2];  t3 = (signed) D[3];  t4 = (signed) D[4];  t5 = (signed) D[5];

	for ( i = 0 ; i < 32 ; i++ ) pts[i] = 1;	// point at infinity
	p0 = p1 = p2  = p3 = p4 = p5 = p6 = p7 = p8 = p9 = p10 = p11 = p12 = p13 = p14 = p15 = 0;
	q0 = q1 = q2  = q3 = q4 = q5 = q6 = q7 = q8 = q9 = q10 = q11 = q12 = q13 = q14 = q15 = 0;
	for ( i = 0 ; i < p ; i++ ) {
		if ( ! t0 ) p0++;
		if ( t0+32 > p ) pts[p-t0]++; 		// assumes p > 32
		z = map[t0>>6];
		if ( (t0&0x3f) > 31 ) {
			z >>= 32;
			z |= map[(t0>>6)+1]<<32;
			z >>= (t0&0x3f) - 32;
		} else {
			z >>= (t0&0x3f);
		}
		_count32(z);
		_advance_5t_seq();
	}
	pts[0] += p0;  pts[1] += p1;  pts[2] += p2;  pts[3] += p3;  pts[4] += p4;  pts[5] += p5;  pts[6] += p6;  pts[7] += p7;
	pts[8] += p8;  pts[9] += p9;  pts[10] += p10;  pts[11] += p11;  pts[12] += p12;  pts[13] += p13;  pts[14] += p14;  pts[15] += p15;
	pts[16] += q0;  pts[17] += q1;  pts[18] += q2;  pts[19] += q3;  pts[20] += q4;  pts[21] += q5;  pts[22] += q6;  pts[23] += q7;
	pts[24] += q8;  pts[25] += q9;  pts[26] += q10;  pts[27] += q11;  pts[28] += q12;  pts[29] += q13;  pts[30] += q14;  pts[31] += q15;
	return 1;
}

// hardwired for 32x
int pointcount_multi_g2d6 (unsigned pts[], unsigned long D[7], unsigned p, unsigned long f6)
{
	register signed  i, t0, t1, t2, t3, t4, t5, t6;
	unsigned p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15;
	unsigned q0, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, q14, q15;
	register unsigned long z;

	if ( p < 32 ) { err_printf ("cannot multicount with p < 32\n");  exit (0); }
		
	pointcount_map_residues (p);
	// wrap table around top
	for ( t2 = 1 ; t2 < 32 ; t2++ )
		if ( map[t2>>6] & tab[t2&0x3f] ) map[(p+t2)>>6] |= tab[(p+t2)&0x3f];
		
	t0 = (signed) D[0];  t1 = (signed) D[1];  t2 = (signed) D[2];  t3 = (signed) D[3];  t4 = (signed) D[4];  t5 = (signed) D[5];  t6 = (signed) D[6];

	for ( i = 0 ; i < 32 ; i++ ) pts[i] = 1 + ((map[f6>>6] & tab[f6&0x3F])?1:-1);		// 2 or 0 points at infinity depending on whether leading coeff is square or not.
		
	p0 = p1 = p2  = p3 = p4 = p5 = p6 = p7 = p8 = p9 = p10 = p11 = p12 = p13 = p14 = p15 = 0;
	q0 = q1 = q2  = q3 = q4 = q5 = q6 = q7 = q8 = q9 = q10 = q11 = q12 = q13 = q14 = q15 = 0;
	for ( i = 0 ; i < p ; i++ ) {
		if ( ! t0 ) p0++;
		if ( t0+32 > p ) pts[p-t0]++; 		// assumes p > 32
		z = map[t0>>6];
		if ( (t0&0x3f) > 31 ) {
			z >>= 32;
			z |= map[(t0>>6)+1]<<32;
			z >>= (t0&0x3f) - 32;
		} else {
			z >>= (t0&0x3f);
		}
		_count32(z);
		_advance_6t_seq();
	}
	pts[0] += p0;  pts[1] += p1;  pts[2] += p2;  pts[3] += p3;  pts[4] += p4;  pts[5] += p5;  pts[6] += p6;  pts[7] += p7;
	pts[8] += p8;  pts[9] += p9;  pts[10] += p10;  pts[11] += p11;  pts[12] += p12;  pts[13] += p13;  pts[14] += p14;  pts[15] += p15;
	pts[16] += q0;  pts[17] += q1;  pts[18] += q2;  pts[19] += q3;  pts[20] += q4;  pts[21] += q5;  pts[22] += q6;  pts[23] += q7;
	pts[24] += q8;  pts[25] += q9;  pts[26] += q10;  pts[27] += q11;  pts[28] += q12;  pts[29] += q13;  pts[30] += q14;  pts[31] += q15;
	return 1;
}

// hardwired for 32x
int pointcount_multi_g3 (unsigned pts[], unsigned long D[8], unsigned p)
{
	register signed  i, t0, t1, t2, t3, t4, t5, t6, t7;
	unsigned p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15;
	unsigned q0, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, q14, q15;
	register unsigned long z;
	unsigned long x;

	if ( p < 32 ) return 0;
		
	pointcount_map_residues (p);
	
	// wrap table around top
	for ( t2 = 0 ; t2 < 32 ; t2++ )
		if ( map[t2>>6] & tab[t2&0x3f] ) map[(p+t2)>>6] |= tab[(p+t2)&0x3f];

	t0 = (signed) D[0];  t1 = (signed) D[1];  t2 = (signed) D[2];  t3 = (signed) D[3];  t4 = (signed) D[4];  t5 = (signed) D[5];  t6 = (signed) D[6];  t7 = (signed) D[7];

	for ( i = 0 ; i < 32 ; i++ ) pts[i] = 1;	// point at infinity

	p0 = p1 = p2  = p3 = p4 = p5 = p6 = p7 = p8 = p9 = p10 = p11 = p12 = p13 = p14 = p15 = 0;
	q0 = q1 = q2  = q3 = q4 = q5 = q6 = q7 = q8 = q9 = q10 = q11 = q12 = q13 = q14 = q15 = 0;
	for ( i = 0 ; i < p ; i++ ) {
		if ( ! t0 ) p0++;
		if ( t0+32 > p ) pts[p-t0]++; 		// assumes p > 32
		z = map[t0>>6];
		if ( (t0&0x3f) > 31 ) {
			z >>= 32;
			z |= map[(t0>>6)+1]<<32;
			z >>= (t0&0x3f) - 32;
		} else {
			z >>= (t0&0x3f);
		}
		_count32(z);
		_advance_7t_seq();
	}
	pts[0] += p0;  pts[1] += p1;  pts[2] += p2;  pts[3] += p3;  pts[4] += p4;  pts[5] += p5;  pts[6] += p6;  pts[7] += p7;
	pts[8] += p8;  pts[9] += p9;  pts[10] += p10;  pts[11] += p11;  pts[12] += p12;  pts[13] += p13;  pts[14] += p14;  pts[15] += p15;
	pts[16] += q0;  pts[17] += q1;  pts[18] += q2;  pts[19] += q3;  pts[20] += q4;  pts[21] += q5;  pts[22] += q6;  pts[23] += q7;
	pts[24] += q8;  pts[25] += q9;  pts[26] += q10;  pts[27] += q11;  pts[28] += q12;  pts[29] += q13;  pts[30] += q14;  pts[31] += q15;
	return 1;
}


// hardwired for 32x
int pointcount_multi_g3d8 (unsigned pts[], unsigned long D[9], unsigned p, unsigned long f8)
{
	register signed  i, t0, t1, t2, t3, t4, t5, t6, t7, t8;
	unsigned p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15;
	unsigned q0, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, q14, q15;
	register unsigned long z;
	unsigned long x;

	if ( p < 32 ) return 0;
		
	pointcount_map_residues (p);
	
	// wrap table around top
	for ( t2 = 0 ; t2 < 32 ; t2++ )
		if ( map[t2>>6] & tab[t2&0x3f] ) map[(p+t2)>>6] |= tab[(p+t2)&0x3f];

	t0 = (signed) D[0];  t1 = (signed) D[1];  t2 = (signed) D[2];  t3 = (signed) D[3];  t4 = (signed) D[4];  t5 = (signed) D[5];  t6 = (signed) D[6];  t7 = (signed) D[7];  t8 = (signed) D[8];

	for ( i = 0 ; i < 32 ; i++ )  1 + ((map[f8>>6] & tab[f8&0x3F])?1:-1);		// 2 or 0 points at infinity depending on whether leading coeff is square or not.

	p0 = p1 = p2  = p3 = p4 = p5 = p6 = p7 = p8 = p9 = p10 = p11 = p12 = p13 = p14 = p15 = 0;
	q0 = q1 = q2  = q3 = q4 = q5 = q6 = q7 = q8 = q9 = q10 = q11 = q12 = q13 = q14 = q15 = 0;
	for ( i = 0 ; i < p ; i++ ) {
		if ( ! t0 ) p0++;
		if ( t0+32 > p ) pts[p-t0]++; 		// assumes p > 32
		z = map[t0>>6];
		if ( (t0&0x3f) > 31 ) {
			z >>= 32;
			z |= map[(t0>>6)+1]<<32;
			z >>= (t0&0x3f) - 32;
		} else {
			z >>= (t0&0x3f);
		}
		_count32(z);
		_advance_8t_seq();
	}
	pts[0] += p0;  pts[1] += p1;  pts[2] += p2;  pts[3] += p3;  pts[4] += p4;  pts[5] += p5;  pts[6] += p6;  pts[7] += p7;
	pts[8] += p8;  pts[9] += p9;  pts[10] += p10;  pts[11] += p11;  pts[12] += p12;  pts[13] += p13;  pts[14] += p14;  pts[15] += p15;
	pts[16] += q0;  pts[17] += q1;  pts[18] += q2;  pts[19] += q3;  pts[20] += q4;  pts[21] += q5;  pts[22] += q6;  pts[23] += q7;
	pts[24] += q8;  pts[25] += q9;  pts[26] += q10;  pts[27] += q11;  pts[28] += q12;  pts[29] += q13;  pts[30] += q14;  pts[31] += q15;
	return 1;
}


unsigned pointcount_slow (ff_t f[], int d, unsigned p)
{
	register unsigned  m, n;
	register ff_t x, y, z;
	register int i;
	
	memset (map, 0, p/8+1);
	_ff_set_ui (z, (p/2)+1);
	_ff_set_one (x);
	do {					// note that we start checking squares at x=1, the square 0 is handled explicitly by the counter m
		_ff_square(y,x);
		map[y>>6] |= tab[y&0x3F];
		_ff_inc (x);
	} while ( ! _ff_equal(x,z) );
	
	m = 1;										// one point at infinity
	_ff_set(y,f[d]);
	if ( !(d&1) ) m += ((map[y>>6] & tab[y&0x3F])?1:-1);  	// add/sub point at infinity when f has even degree
	n = 0;
	
	_ff_set_zero(x);
	do {
		_ff_set (y, f[d]);
		for ( i = d-1 ; i >= 0 ; i-- ) {
			ff_mult (y, y, x);
			_ff_addto (y, f[i]);
		}
		if ( _ff_zero(y) ) m++;
		if ( (map[y>>6] & tab[y&0x3F]) ) n++;
		_ff_inc(x);
	} while ( _ff_nonzero(x) );
	return m+2*n;
}


// basically the same as pointcount_slow, only it doesn't require ff initialization and will handle the prime 3
unsigned pointcount_tiny (unsigned long f[], int d, unsigned p)
{
	register unsigned  m, n;
	register unsigned long x, y, z;
	register int i;
	
	if ( p == 2 ) { err_printf ("p==2 not supported in pointcount_tiny_d2\n");  exit (0); }

	// avoid worrying about overflow
	if ( p > (1<<20) ) { err_printf ("p=%d not tiny enough in pointcount_tiny_d2\n", p);  exit (0); }

	pointcount_map_residues (p);

	m = 1;											// one point at infinity
	if ( !(d&1) ) m += ((map[f[d]>>6] & tab[f[d]&0x3F])?1:-1);  	// add/sub point at infinity when f has even degree
	n = 0;
	
	x = 0;
	do {
		y = f[d];
		for ( i = d-1 ; i >= 0 ; i-- ) y = (x*y + f[i]) % p;
		if ( !y ) m++;
		if ( (map[y>>6] & tab[y&0x3F]) ) n++;
		x++;
	} while ( x < p );
	return m+2*n;		// add point at infinity
}


unsigned pointcount_g2_d2 (unsigned long f[6], unsigned p)
{
	register signed i, j, k, t00, t01, t10, t11, t20, t21, t30, t31, t40, t41, t50, t51;
	register unsigned m, n;
	register unsigned long s, x0, x1, y0, y1, t0, t1;
	
	// handle p <= d seperately, otherwise it interferes with our \Delta computations
	if ( p <= 5 ) return pointcount_tiny_d2 (f, 5, p);

	s = pointcount_map_residues_d2 (p);

	m = 1;				// one point at infinity for odd degree
	n = 0;

	// t_5 = \Delta^5 f(x) is constant
	t51 = 0;  t50 = (120*f[5]) % p;

	for ( k = 0 ; k < p ; k++ ) {
		// Compute t_j = \Delta^j f(kz) where z = sqrt(s) is our extension element in F_p^2 = F[x]/(z^2-s)
		// We do this the old fashioned way, computing f(kz), f(kz+1), f(kz+2), ..., f(kz+d-1) and taking differences
		// Clumsy but effective - speed is not really an issue here
		x1 = k;  x0 = 0;  y1 = 0;  y0 = f[5];
		for ( i = 4 ; i >= 0 ; i-- ) { t0 = (x0*y0+s*x1*y1+ f[i]) % p;  t1 = (x0*y1+x1*y0) % p;  y0 = t0;  y1 = t1; }
		t00 = (unsigned) t0;  t01 = (unsigned) t1;
		x1 = k;  x0 = 1;  y1 = 0;  y0 = f[5];
		for ( i = 4 ; i >= 0 ; i-- ) { t0 = (x0*y0+s*x1*y1+ f[i]) % p;  t1 = (x0*y1+x1*y0) % p;  y0 = t0;  y1 = t1; }
		t10 = (unsigned) t0;  t11 = (unsigned) t1;
		x1 = k;  x0 = 2;  y1 = 0;  y0 = f[5];
		for ( i = 4 ; i >= 0 ; i-- ) { t0 = (x0*y0+s*x1*y1+ f[i]) % p;  t1 = (x0*y1+x1*y0) % p;  y0 = t0;  y1 = t1; }
		t20 = (unsigned) t0;  t21 = (unsigned) t1;
		x1 = k;  x0 = 3;  y1 = 0;  y0 = f[5];
		for ( i = 4 ; i >= 0 ; i-- ) { t0 = (x0*y0+s*x1*y1+ f[i]) % p;  t1 = (x0*y1+x1*y0) % p;  y0 = t0;  y1 = t1; }
		t30 = (unsigned) t0;  t31 = (unsigned) t1;
		x1 = k;  x0 = 4;  y1 = 0;  y0 = f[5];
		for ( i = 4 ; i >= 0 ; i-- ) { t0 = (x0*y0+s*x1*y1+ f[i]) % p;  t1 = (x0*y1+x1*y0) % p;  y0 = t0;  y1 = t1; }
		t40 = (unsigned) t0;  t41 = (unsigned) t1;
		t40 -= t30;  if ( t40 < 0 ) t40 += p;  t41 -= t31;  if ( t41 < 0 ) t41 += p;
		t30 -= t20;  if ( t30 < 0 ) t30 += p;  t31 -= t21;  if ( t31 < 0 ) t31 += p;
		t20 -= t10;  if ( t20 < 0 ) t20 += p;  t21 -= t11;  if ( t21 < 0 ) t21 += p;
		t10 -= t00;  if ( t10 < 0 ) t10 += p;  t11 -= t01;  if ( t11 < 0 ) t11 += p;
		t40 -= t30;  if ( t40 < 0 ) t40 += p;  t41 -= t31;  if ( t41 < 0 ) t41 += p;
		t30 -= t20;  if ( t30 < 0 ) t30 += p;  t31 -= t21;  if ( t31 < 0 ) t31 += p;
		t20 -= t10;  if ( t20 < 0 ) t20 += p;  t21 -= t11;  if ( t21 < 0 ) t21 += p;
		t40 -= t30;  if ( t40 < 0 ) t40 += p;  t41 -= t31;  if ( t41 < 0 ) t41 += p;
		t30 -= t20;  if ( t30 < 0 ) t30 += p;  t31 -= t21;  if ( t31 < 0 ) t31 += p;
		t40 -= t30;  if ( t40 < 0 ) t40 += p;  t41 -= t31;  if ( t41 < 0 ) t41 += p;
		for ( j = 0 ; j < p ; j++ ) {
			i = p*t01+t00;
			if ( ! i ) m++;
			if ( (map[i>>6] & tab[i&0x3F]) ) n++;
			t00 += t10; if ( t00 >= p ) t00 -= p;  t10 += t20;  if ( t10 >= p ) t10 -= p; t20 += t30; if ( t20 >= p ) t20 -= p; t30 += t40; if ( t30 >= p ) t30 -= p; t40 += t50; if ( t40 >= p ) t40 -= p;
			t01 += t11; if ( t01 >= p ) t01 -= p;  t11 += t21;  if ( t11 >= p ) t11 -= p; t21 += t31; if ( t21 >= p ) t21 -= p; t31 += t41; if ( t31 >= p ) t31 -= p; t41 += t51; if ( t41 >= p ) t41 -= p;
		}
	}
	return m+2*n;
}


// This code is brutally slow and should only be used for very small p
unsigned pointcount_slow_d2 (ff_t f[], int d, unsigned p)
{
	register unsigned m, n;
	register ff_t  x0, x1, y0, y1, t0, t1, t, s;
	register int i;

	
	if ( ! (d&1) ) { err_printf ("degree not odd in pointcount_slow_d2\n");  exit (0); }
	m = pointcount_map_residues_d2 (p);
	_ff_set_ui (s, m);
	m = 1;										// one point at infinity
	n = 0;

	_ff_set_zero(x1);
	do {
		_ff_set_zero(x0);
		do {
			// compute f(x) = f(x1*z+x0) where z = sqrt(s) is our extension element, F_p^2 = F[x]/(z^2-s)
			_ff_set_zero (y1);  _ff_set (y0, f[d]);
			for ( i = d-1 ; i >= 0 ; i-- ) {
				_ff_mult (t0, x0, y0);  _ff_mult (t, x1, y1);  ff_mult (t, t, s);  _ff_addto (t0, t);
				_ff_mult (t1, x0, y1);  _ff_mult (t, x1, y0); _ff_addto (t1, t);
				_ff_set (y1, t1);  _ff_add (y0, t0, f[i]);
			}
			if ( _ff_zero(y1) && _ff_zero(y0) ) m++;
			// pull out of field representation to compute table index
			i = (int) (_ff_get_ui (y1) * p + _ff_get_ui (y0));
			if ( ( map[i>>6] & tab[i&0x3f]) ) n++;
			_ff_inc(x0);
		} while ( _ff_nonzero(x0) );
		_ff_inc(x1);
	} while ( _ff_nonzero(x1) );
	return m+2*n;
}


/*
	This code is not very efficient and should be changed to use routines in ffext.c.
*/
unsigned pointcount_tiny_d2 (unsigned long f[], int d, unsigned p)
{
	register unsigned m, n;
	register unsigned  x0, x1, y0, y1, t0, t1, s;
	register int i;
	
	if ( p == 2 ) { err_printf ("p==2 not supported in pointcount_tiny_d2\n");  exit (0); }
	
	// avoid worrying about overflow
	if ( p > (1<<14) ) { err_printf ("p=%d not tiny enough in pointcount_tiny_d2\n", p);  exit (0); }
	
	s = pointcount_map_residues_d2 (p);
	m = 1;												// one point at infinity
	if ( !(d&1) ) m += ((map[f[d]>>6] & tab[f[d]&0x3F])?1:-1);  	  	// add/sub point at infinity when f has even degree
	n = 0;

	x1 = 0;
	do {
		x0 = 0;
		do {
			// compute f(x) = f(x1*z+x0) where z = sqrt(s) is our extension element, F_p^2 = F[z]/(z^2-s)
			y1 = 0;  y0 = f[d];
			for ( i = d-1 ; i >= 0 ; i-- ) {
				t0 = (x0*y0+s*x1*y1+ f[i]) % p;
				t1 = (x0*y1+x1*y0) % p;
				y0 = t0;  y1 = t1;
			}
			if ( ! y0 && ! y1 ) m++;
			i = p*y1+y0;
			if ( ( map[i>>6] & tab[i&0x3f]) ) n++;
			x0++;
		} while ( x0 < p );
		x1++;
	} while ( x1 < p );
	return m+2*n;
}

/*
	This code is not very efficient and should be changed to use routines in ffext.c.
*/
unsigned pointcount_tiny_d3 (unsigned long f[], int d, unsigned p)
{
	register unsigned m, n;
	register unsigned  x0, x1, x2, y0, y1, y2, t0, t1, t2, r1, r2, s;
	register int i;
	
	if ( p == 2 ) { err_printf ("p==2 not supported in pointcount_tiny_d2\n");  exit (0); }
	
	// avoid worrying about overflow 
	if ( p > (1<<10) ) { err_printf ("p=%d not tiny enough in pointcount_tiny_d2\n", p);  exit (0); }
	
	s = pointcount_map_residues_d3 (p);
	m = 1;
	if ( !(d&1) ) m += ((map[f[d]>>6] & tab[f[d]&0x3F])?1:-1);  	  	// add/sub point at infinity when f has even degree
	n = 0;

	x2 = 0;
	do {
		x1 = 0;
		do {
			x0 = 0;
			do {
				// compute f(x) = f(x2*z^2+x1*z+x0) where z is our extension element in F_p^3 = F[z]/(z^3-z-s)
				y2 = 0;  y1 = 0;  y0 = f[d];
				for ( i = d-1 ; i >= 0 ; i-- ) {
					r1 = x1*y2+x2*y1;  r2 = x2*y2;
					t0 = (x0*y0 + s*r1 + f[i]) % p;
					t1 = (x0*y1 + x1*y0 + r1 + s*r2) % p;
					t2 = (x0*y2 + x1*y1 + x2*y0  + r2) %p;
					y0 = t0;  y1 = t1;  y2 = t2;
				}
				i = p*p*y2 + p*y1+y0;
				if ( ! i ) m++;
				if ( ( map[i>>6] & tab[i&0x3f]) ) n++;
				x0++;
			} while ( x0 < p );
			x1++;
		} while ( x1 < p );
		x2++;
	} while (x2 < p );
	return m+2*n;
}
