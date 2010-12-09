#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include "gmp.h"
#include "mpzutil.h"
#include "jac.h"
#include "hecurve.h"
#include "ff.h"
#include "smalljactab.h"
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

/*
    As used by smalljac, the CONFIDENCE parameter does not affect the correctness of the output - smalljac
    always proves its results unconditionally.  Setting CONFIDENCE too low could cause smalljac to fail, however.
*/

#define CONFIDENCE			40					// must be less than 64 - set a bit high because we don't trust the uniformity of random elements in genus 3

int jac_vector_logarithm (unsigned long e[], jac_t a[], unsigned long ords[], int k, jac_t b[1], curve_t c[1]);

/*
	Given the order (or multiple of the exponent) of the group, N = #J(C/F_q), jac_structure computes the abelian group structure of J(C/F_q)
	as a product of cylclic groups Z_m[0] x ... x Z_m[n-1] with m[0] | m[1] | ... | m[n-1] and returns n.  The flag fExponent, when set, indicates that N is only
	guaranteed to be a multiple of the group exponent, not necessarily its order.
	
	Generators are not computed, but they easily could be (at a slight increase in cost).
*/
int jac_structure (long m[], curve_t c[1], unsigned long N, int fExponent)
{
	unsigned long h[JAC_MAX_FACTORS];
	unsigned long p[JAC_MAX_FACTORS];
	unsigned long t[JAC_MAX_GENERATORS];
	unsigned long ords[JAC_MAX_GENERATORS];
	jac_t x[1], a[JAC_MAX_GENERATORS];
	register int i, j, k, n, r, w, maxj;
	register unsigned long max, q, M;

	w = ui_factor (p, h, N);
	if ( ! w ) { err_printf ("ui_factor_small failed");  exit (0); }
	n = 1;
	for ( i = 0 ; i < c->d-1 ; i++ ) t[i] = 1;
	// could add code here to optimize rank 1 2-Sylow subgroups by counting the factors of f(x)
	for ( i = 0 ; i < w ; i++ ) {
		if ( ! fExponent ) {
			// The optimizations below rely on N being the group order
			if ( h[i] == 1 ) { t[0] *= p[i];  continue; }
			for ( q = p[i]*p[i], j = 2 ; j < h[i] ; j++ ) q *= p[i];
			// Todo: Add genus 1 check for divisibility of _ff_p-1
			// Before performing a p-Sylow subgroup computation, expend some effort trying to
			// prove the p-Sylow subgroup is cyclic - we expect this to often be the case and
			// it only requires O(lg(N)) gops.
			M = N/p[i];		// was N/q which can't be right
			for ( j = 0 ; j < JAC_CYCLIC_TESTS ; j++ ) {
				_jac_random (x[0], c[0]);
				jac_exp_ui (x, x, M, c);
				if ( ! _jac_is_identity (x[0]) ) break;
			}
			if ( j < JAC_CYCLIC_TESTS ) { for ( j = 0 ; j < h[i] ; j++ ) t[0] *= p[i]; continue; }
		} else {
			q = 0;
		}
		r = jac_sylow (a, ords, p[i], N, q, 0, c);
		if ( r <= 0 ) { err_printf ("%lu-Sylow computation failed\n", p[i]);  exit (0); }
		for ( j = 1, M = ords[0] ; j < r ; j++ ) M *= ords[j];
		if ( q && M != q ) {
			err_printf ("jac_sylow computed subgroup of rank %d, size %lu when size %lu was expected (N=%ld)\n", r, M, q, N);
			printf ("p=%ld: ", _ff_p);  _curve_print (c[0]);
			exit (0);
		}
		// pull p-Sylow factors out in decreasing order by size (largest first).  Don't bother sorting, there can only be a few
		if ( n < r ) n = r;
		for ( j = 0 ; j < r ; j++ ) {
			for ( max = 0, k = 0 ; k < r ; k++ ) if ( ords[k] > max ) { max = ords[k];  maxj = k; }
			t[j] *= max;
			ords[maxj] = 0;
		}
	}
// sanity checks - can be removed eventually
if ( n > c->d-1 ) { err_printf ("jac_structure failed - computed rank %d > 2g = %d\n", n, c->d-1);  exit (0); }
if ( ! fExponent ) {
for ( M = 1, i = 0 ; i < n ; i++ ) M *= t[i];
if ( M != N ) { err_printf ("jac_structure failed - computed order %lu not equal to given order %lu\n", M, N);  exit (0); }
}
for ( i = 0 ; i < n-1 ; i++ ) if ( t[i]%t[i+1] ) { err_printf ("jac_structure failed, cyclic order divisibility test failed t[i] = %lu, t[i+1] = %lu\n", t[i], t[i+1]);  exit (0); }
	// copy results into output array, reversing order to put smallest first
	for ( i = 0 ; i < n ; i++ ) m[i] = t[n-i-1];
	return n;
}


/*
    jac_sylow computes the p-Sylow subgroup of the jacobian using a provided exponent E (E can be any multiple of the group exponent)
    M is an upper bound on the maximum size of the p-Sylow subgroup - if this value is too small, jac_sylow will compute a
    subgroup which may not be the full p-Sylow subgroup.    It is assumed that the order of the p-Sylow subgroup fits in an unsigned long.
    
    The optional parameter limit will force jac_sylow to terminate with an error if the p-Sylow subgroup is found to be larger than limit.
    This is used to ensure an O(|G|^{1/4}) running time.  See Proposition 2 of [KedlayaSutherland2007].

    The computed subgroup is in the form of a basis a[] with the order of each basis element in ords[].  The size of the basis
    (the p-rank of the Jacobian) is returned.  If maxgops is exceeded, -1 is returned.

    This algorithm is an impelmentation of Algorithm 9.1 (including Algorithm 9.2) of [SutherlandThesis], specialized to a single Sylow subgroup.
*/  

int jac_sylow (jac_t a[JAC_MAX_GENERATORS], unsigned long ords[JAC_MAX_GENERATORS], unsigned p, unsigned long E, unsigned long M, unsigned long limit, curve_t c[1])
{
	jac_t b;
	unsigned long h;
	unsigned long e[JAC_MAX_GENERATORS];
	unsigned long m, q, x, order;
	int k, n, cnt, retries, min, d, temp;
	int sts;

	for ( h = 0 ; ! (E%p) ; h++ ) E/= p;
	if ( ! h ) { err_printf ("Trivial call to smalljac_sylow E = %lu, p = %u\n", E, p);  return 0; }
	k = 0;
	order = 1;
	// Get initial generator
	for ( cnt = 0 ; cnt < CONFIDENCE ; cnt++ ) {
		_jac_random (a[k], c[0]);
		jac_exp_ui (a+k, a+k, E, c);
		_jac_set (b, a[k]);
		ords[k] = jac_pp_order_ui (&b, p, h, c);
		if ( ords[k] > 1 ) break;
	}
	if ( ords[k] < 2 ) { err_printf ("smalljac_sylow couldn't find any generators for %d-Sylow subgroup - bad exponent %ld?!\n", p, E);  exit (0); }
	dbg_printf ("First generator for %d-sylow subgroup has order %d, E = %ld, M = %ld, limit = %ld\n", p, ords[k], E, M, limit);
	k++;
	// adjust retries based on p to achieve probabilty > 2^confidence
	for ( retries = 1, x = p ; x < (1UL<<CONFIDENCE) ; x *= p, retries++ );
	for ( cnt = 0 ; cnt < retries ; ) {
		for ( n = 0, x = 1 ; n < k ; x *= ords[n], n++ );
		// if we have exceeded the limit, return an error
		if ( limit && x > limit ){ dbg_printf ("Sylow subgroup computation terminated unsuccessfully due to %ld*%ld > %ld\n", p, x, limit);  return -1; }
		// if adding another basis element would make group size > M then we must be done
		if ( M && p*x > M ) { dbg_printf ("Sylow subgroup computation successfully terminated due to bound M = %ld\n", M);  break; }
		_jac_random (a[k], c[0]);
		jac_exp_ui (a+k, a+k, E, c);
repeat:
		if ( _jac_is_identity (a[k]) ) continue;
		_jac_set (b, a[k]);
		_jac_invert (b);
		sts = jac_vector_logarithm (e, a, ords, k, &b, c);
		if ( sts > 0 ) { cnt++;  continue; }
		if ( sts < 0 ) return -1;
		// if an unspanned element forces us past limit, return error
		if ( limit && p*x > limit ){ dbg_printf ("Sylow subgroup computation terminated unsuccessfully due to %ld*%ld > %ld\n", p, x, limit);  return -1; }
		ords[k] = jac_pp_order_ui (&b, p, h, c);
		if ( ords[k] < 2 ) { err_printf ("unexpected order %d returned from jac_pp_order_ui", ords[k]);  exit (0); }
		dbg_printf ("Found unspanned element with order %d\n", ords[k]);
		for ( m = p ; m < ords[k] ; m *= p ) {
			jac_exp_ui (&b, &b, p, c);
			sts = jac_vector_logarithm (e, a, ords, k, &b, c);
			if ( sts > 0 ) break;
			if ( sts < 0 ) return -1;
		}
		if ( m < ords[k] ) {
			for ( n = 0 ; n < k ; n++ ) {
				if ( e[n] ) {
					q = ui_pp_div (e[n], p);
					if ( q != e[n] ) {
						jac_exp_ui (a+n, a+n, e[n]/q, c);
						e[n] = q;
					}
				}
			}	
			min = -1;  d = m;
			for ( n = 0 ; n < k ; n++ )
				if ( e[n] && e[n] < d ) { min = n;  d = e[n]; }
			dbg_printf ("min = %d, d = %d\n", min, d);
			if ( min == -1 ) {
				for ( n = 0 ; n < k ; n++ ) {
					if ( e[n] != 0 ) {
						if ( (e[n]%m) != 0 ) { out_printf ("Problem, %d not divisible by %d\n", e[n], m);  exit (1); }
						jac_exp_ui (&b, a+n, e[n]/m, c);
						_jac_mult (a[k], a[k], b, c[0]);
					}
				}
				ords[k] = m;
				temp = jac_pp_order_ui (a+k, p, h, c);
				if ( temp != ords[k] ) { err_printf ("%lu: Order mismatch 1, %d != %d for element ", _ff_p, temp, ords[k]);  _jac_print(a[k]);  exit (0); }
				dbg_printf ("Added adjusted generator using relation min == -1.\n");
				k++;
			} else {
				if ( d == 1 ) {
					dbg_printf ("Removed generator with exponent 1, recomputing relation...\n", e[n], m);
					n = min;
					while ( n < k ) {
						ords[n] = ords[n+1];
						_jac_set (a[n], a[n+1]);
						temp = jac_pp_order_ui (a+n, p, h, c);
						if ( temp != ords[n] ) { err_printf ("Order mismatch 2, %d != %d for element ", temp, ords[n]);  _jac_print(a[n]);  exit (0); }
						n++;					 	
					}
					k--;
					if ( ! k) {
						dbg_printf ("no more generators left, independent element is new generator\n");
						k++;
					} else {							
						dbg_printf ("repeating with k = %d\n", k);
						goto repeat;
					}
				} else {
					if ( (m%e[min]) != 0 ) { out_printf ("Problem, %d not divisible by %d\n", e[min], m);  exit (1); }
					jac_exp_ui (&b, a+k, m/e[min], c);
					_jac_mult (a[min], a[min], b, c[0]);
					for ( n = 0 ; n < k ; n++ ) {
						if ( n != min && e[n] ) {
							if ( (e[n]%e[min]) != 0 ) { out_printf ("Problem, %d not divisible by %d\n", e[min], m);  exit (1); }
							jac_exp_ui (&b, a+n, e[n]/e[min], c);
							_jac_mult (a[min], a[min], b, c[0]);
						}
					}
					dbg_printf ("Adjusted existing generator using relation, changed order from %d to %d.\n", ords[min], e[min]);
					ords[min] = e[min];
					temp = jac_pp_order_ui (a+min, p, h, c);
					if ( temp != ords[min] ) { err_printf ("Order mismatch 3, %d != %d for element ", temp, ords[min]);  _jac_print(a[min]);  exit (0); }
					goto repeat;
				}
			}
		} else {
			dbg_printf ("Added new independent generator %d with order %d\n", k, ords[k]);
			k++;
		}
//printf ("%d generators: \n", k);
//for ( n = 0 ; n < k ; n++ ) _jac_print(a[n]);  printf ("  order %d\n", ords[n]);
	}
	return k;
}


/*
     jac_vector_logarithm implements Algorithm 9.3 of [SutherlandThesis].  It is assumed the order of the group
     spanned by the provided basis (i.e. the product of the ords[i]) fits in a long.
     
     This algorithm is a straight copy of the generic version and does not include any optimizations for fast inverses
     or parallel group operations.   smalljac uses this very infrequently, so we don't bother optimizing it.
*/

int jac_vector_logarithm (unsigned long e[], jac_t a[], unsigned long ords[], int k, jac_t beta[1], curve_t curve[1])
{
	static int stashsize;
	static jac_t *stash;
	jac_t g, h, ht, betainverse;
	unsigned long E;
	register long i, b, c, ct, d, u, s, t, j, l;
	uint32_t matches[SMALLJAC_MAX_MATCHES];
	long m, n;
	int sts, bits, found;
	
	if ( _jac_is_identity (beta[0]) ) { for ( i = 0 ; i < k ; i++ ) e[i] = 0;  return 1; }
	if ( k == 0 ) { err_printf ("Invalid k value in grp_vector_logarithm");  exit (0); }

	sts = 0;

	// quick sanity check
	for ( i = 0 ; i < k ; i++ ) if ( ords[i] < 2 ) { out_printf ("Invalid order %d specified in grp_vector_logarithm", ords[i]);  exit (0); }
	
//out_printf ("Vector logarithm searching for ");  _jac_print (beta[0]);
//for ( i = 0 ; i < k ; i++ ) { out_printf ("(%d) ",ords[i]); _jac_print (a[i]); }
	
	for ( E = ords[0], i = 1 ; i < k ; i++ ) E *= ords[i];
//out_printf ("E = %lu\n", E);

	c = (long) floor(sqrt(E));
	b = 1;
	for ( d = 0 ; d < k ; d++ ) {
		if ( b*ords[d] > c ) break;			// overflow could happen here if the product of ords[d] doesn't fit in a long
		b *= ords[d];
	}
	c /= b;
	l = b;
	b *= c;

	// Note that we will be giant stepping via alpha[d]^c
	
	// initialize table and init/expand stash if needed
	for ( bits = 10, i = (1<<10) ; i < 2*b ; i <<= 1, bits++ );		// aim for a load factor of 0.5
	smalljac_table_init (bits);
	if ( b >= stashsize ) {
		if ( stash ) {
			for ( stashsize *= 2 ; stashsize < b ; stashsize *= 2 );
			stashsize *= 2;  stash = (jac_t *) realloc (stash, stashsize*sizeof(*stash));
		} else {
			for ( stashsize = 1024 ; stashsize < b ; stashsize *= 2 );
			stash = (jac_t *) malloc (stashsize*sizeof(*stash));
		}
	}

	_jac_set_identity (g);
	// baby steps - we currently don't bother trying to do these in parallel
	for ( i = 1 ; i < b ; i++ ) {
		_jac_mult (g, g, a[0], curve[0]);
		for ( j = 0, m = ords[0] ; (i%m) == 0 ; m *= ords[++j] ) _jac_mult (g, g, a[j+1], curve[0]);
		if ( _jac_cmp (g, beta[0]) == 1 ) break;
		smalljac_table_insert (&g, i);
		_jac_set (stash[i], g);
	}
	if ( _jac_cmp (g, beta[0]) == 1 ) { 
		for ( j = 0 ; j < k ; j++ ) {
			e[j] = i%ords[j];
			i /= ords[j];
		}
		sts = 1;
		goto finish;
	}
	// giant steps - step counter is combination of t for the kth digit and u for the higher digits
	_jac_set (betainverse, beta[0]);
	_jac_invert (betainverse);
	_jac_set (g, betainverse);
	jac_exp_ui (&h, a+d, c, curve);
	// Setup ht to make top step as we cycle on the center digit
	if ( c > 1 ) {
		ct = ords[d]-1 - c*((ords[d]-1)/c);
		jac_exp_ui (&ht, a+d, ct, curve);
	} else {
		_jac_set_identity (ht);
	}
	found = 0;
	for ( s = 1, t = 0, u = 0 ; ! found ; s++ ) {
		// Need to hit top of the center digit every time as we cycle through
		if ( t < ords[d]-1 && t+c >= ords[d] ) {
			_jac_mult (g, g, ht, curve[0]);
			t += ct;
			if ( t != ords[d]-1 ) { err_printf ("ct=%d invalid t = %d, c=%d, ords[d] = %d - didn't hit ords[d]-1 on center digit!", ct, t, c, ords[d]);  exit (0); }
		} else if ( t == ords[d]-1 ) {
			// carry
			_jac_mult (g, g, a[d], curve[0]);
			if ( d == k-1 ) goto done;
			t = 0;			// Note, must set t to 0 here, not reduce mod c
			_jac_mult (g, g, a[d+1], curve[0]);
			u++;
			// propogate carry as required
			for ( j = d+1, m = ords[d+1] ; (u%m) == 0 ; m *= ords[++j] ) {
				if ( j == k-1 ) goto done;
				_jac_mult (g, g, a[j+1], curve[0]);
			}
		} else {
			t += c;
			_jac_mult (g, g, h, curve[0]);
		}
		if ( _jac_is_identity (g) ) {
			m = 0;
			found = 1;
		} else {
			n = smalljac_table_lookup (matches, &g);
			for ( i = 0 ; i < n ; i++ ) if ( _jac_cmp(stash[matches[i]], g) == 1 ) break;
			if ( i < n ) { found = 1;  m = matches[i]; } else found = 0;
		}
	}
done:
	// If not found take one final step to the top
	if ( ! found ) {
		_jac_set (h, a[0]);
		for ( i = 1 ; i < k ; i++ ) {
			_jac_mult (h, h, a[i], curve[0]);
		}
		_jac_invert (h);
		_jac_mult (g, h, betainverse, curve[0]);
		n = smalljac_table_lookup (matches, &g);
		for ( i = 0 ; i < n ; i++ ) if ( _jac_cmp(stash[matches[i]], g) == 1 ) break;
		if ( i < n ) { found = 1;  m = matches[i]; } else found = 0;
		if ( ! found && ! _jac_is_identity (g) ) goto finish;
		for ( i = 0 ; i < k ; i++ ) {
			e[i] = (ords[i] - 1 - (m%ords[i])) % ords[i];
			m /= ords[i];
		} 
	} else {
		j = m%l;
		for ( i = 0 ; i < d ; i++ ) {
			e[i] = (ords[i] - (j%ords[i])) % ords[i];
			j /= ords[i];
		}
		m /= l;
		if ( m <= t ) e[d] = t-m; else e[d] = ords[d]+t-m;
		for ( i = d+1 ; i < k ; i++ ) {
			e[i] = u%ords[i];
			u /= ords[i];
		}
	}
	sts = 1;
finish:

	// verification code - comment this out eventually
	if ( sts ) {
//out_printf ("Vector logarithm found (k=%d) ", k);  _jac_print (beta[0]);
		_jac_set_identity (g);
		for ( i = 0 ; i < k ; i++ ) {
			jac_exp_ui (&h, a+i, e[i], curve);
			_jac_mult (g, g, h, curve[0]);
//out_printf ("%d* ", e[i]);  _jac_print(a[i]);
			if ( e[i] == ords[i] ) { puts ("Unreduced exponent!");  exit (1); }
		}
//puts ("");
		if ( ! _jac_cmp (beta[0], g) == 1 ) { err_printf ("Vector Logarithm invalid - orders suspect!\n"); _jac_print(beta[0]);  _jac_print(g); sts = -1; }
	} else {
//out_printf ("Vector logarithm failed to find ");  _jac_print(beta[0]);
	}

	return sts;
}
