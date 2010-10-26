//
// helper routines for L-function computations
//
// - induced ordering on elements of: F_p[x], F_q, F_q[x]
// - converting elements to unsigned longs: F_p[x], F_q(, F_q[x])
// - enumerating elements: F_p[x], F_q(, F_q[x])
// - evaluating elements of F_q[x] at elements of F_q
// - exponentiation in F_q^*

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <NTL/ZZX.h>
#include <NTL/lzz_p.h>
#include <NTL/lzz_pE.h>
#include <NTL/lzz_pX.h>
#include <NTL/lzz_pEX.h>
#include <NTL/lzz_pXFactoring.h>
#include <NTL/lzz_pEXFactoring.h>
#include <NTL/lzz_pX.h>

#include "helper.h"
#include "lzz_pEExtra.h"

NTL_CLIENT

#define NO_CARRY   0
#define CARRY      1

///////////////////////////////////////////////////////////////////////////

// return a modulus appropriate for F_q=F_{p^d}

void get_modulus(zz_pX& pi, int p, int d)
{
    zz_p::init(p);
    SetX(pi);
    pi <<= d-1;
    while (IterIrredTest(pi) == 0 && inc(pi,d-1) == 0)
        ;
    assert(IterIrredTest(pi) == 1 && deg(pi) == d);
}

// return moduli pi_1 for F_{p^d1}/F_p and pi_2 for F_{p^{d1*d2}}/F_p and
// identify a zero of pi_1 in F_{p^{d1*d2}}

void get_modulus(zz_pX& pi_1, zz_pX& pi_2, zz_pX& a, int p, int d1, int d2)
{
    get_modulus(pi_1, p, d1);
    get_modulus(pi_2, p, d1*d2);

    // find alpha
    zz_pE::init(pi_2);

    zz_pEX   pi;
    conv(pi, pi_1);

    zz_pE    zero;

    FindRoot(zero, pi);

    a = rep(zero);
}

// initialize F_{p^d} with canonical irreducible in F_p[x]

void init_NTL_ff(int p, int d, int precompute_inverses,
    int precompute_square_roots, int precompute_legendre_char,
    int precompute_pth_frobenius_map)
{
    zz_pX   pi;
    get_modulus(pi, p, d);
    zz_pE::init(pi);

    // make sure size of finite field fits in a long
    assert(zz_pE::cardinality().WideSinglePrecision());

    zz_pEExtraContext  c(precompute_inverses,
			 precompute_square_roots,
			 precompute_legendre_char,
                         precompute_pth_frobenius_map);
    c.restore();
}

///////////////////////////////////////////////////////////////////////////

// compare elements in F_p[x] ordered by degree then leading coefficient

long operator<(zz_pX& f, zz_pX& g)
{
    int     i;

    if (deg(f) < deg(g))  return 1;
    if (deg(f) > deg(g))  return 0;
    for (i = deg(f); i >= 0; i--) {
	if (coeff(f,i).LoopHole() < coeff(g,i).LoopHole())  return 1;
	if (coeff(f,i).LoopHole() > coeff(g,i).LoopHole())  return 0;
    }

    return 0;
}
inline long operator<=(zz_pX& f, zz_pX& g)
    { return (f < g) || (f == g); }
inline long operator>( zz_pX& f, zz_pX& g)
    { return (g <= f); }
inline long operator>=(zz_pX& f, zz_pX& g)
    { return (g < f); }

// compare elements of F_q

long operator<(zz_pE& x, zz_pE& y)
{
    return x.LoopHole() < y.LoopHole();
}
inline long operator<=(zz_pE& f, zz_pE& g)
    { return (f < g) || (f == g); }
inline long operator>(zz_pE& f, zz_pE& g)
    { return (g <= f); }
inline long operator>=(zz_pE& f, zz_pE& g)
    { return (g < f); }

// compare elements in F_q[x] ordered by degree then leading coefficient

long operator<(zz_pEX& f, zz_pEX& g)
{
    int     i;

    if (deg(f) < deg(g))  return 1;
    if (deg(f) > deg(g))  return 0;
    for (i = deg(f); i >= 0; i--) {
	zz_pE	cf = coeff(f,i), cg = coeff(g,i);
	if (cf < cg)  return 1;
	if (cf > cg)  return 0;
    }

    return 0;
}
inline long operator<=(zz_pEX& f, zz_pEX& g)
    { return (f < g) || (f == g); }
inline long operator>( zz_pEX& f, zz_pEX& g)
    { return (g <= f); }
inline long operator>=(zz_pEX& f, zz_pEX& g)
    { return (g < f); }

///////////////////////////////////////////////////////////////////////////

// convert elt of F_p[x] to int using coefficients as base-p digits
// used to determine index of elt in table

unsigned long to_ulong(const zz_pX& x)
{
    unsigned long  u;
    int     i;
    static zz_p    c;

    c = LeadCoeff(x);
    static long p;
    p = c.modulus();
    for (u = 0, i = deg(x); i >= 0; i--) {
	GetCoeff(c, x,i);
	u = u * p + rep(c);
    }

    return u;
}

void from_ulong(unsigned long ul, zz_pX& x)
{
    static zz_pX   t_to_i;
    static long    p;
    static zz_p    c;

    // get current field characteristic
    t_to_i = 1;
    c = ConstTerm(t_to_i);
    p = c.modulus();

    // convert base-p number ul to polynomial
    x = 0;
    while (ul != 0) {
	x += (ul%p)*t_to_i;
	t_to_i <<= 1;
	ul /= p;
    }
}

void from_ulong(unsigned long ul, zz_pE& x)
{
    static zz_pX  rep;
    from_ulong(ul, rep);
    x = to_zz_pE(rep);
}

unsigned long to_ulong(zz_pE& x)
{
    static zz_pX  y;
    y = x.LoopHole();
    return to_ulong(y);
}

// increment polynomial as if coefficients were base-p digits of an
// integer.  corresponding sequence of polynomials is a so-called lexical
// ordering.

int inc(zz_pX& x, long max_deg)
{
    zz_pX   tee_to_i(0,1);
    int    i;

    i = 0;
    do {
	x += tee_to_i;
	tee_to_i <<= 1;
    } while (coeff(x,i) == 0 && ++i <= max_deg);

    return (i > max_deg) ? CARRY : NO_CARRY;
}

int inc(zz_pEX& x, long max_deg)
{
    zz_pEX   tee_to_i(0,1);
    int      i, r;
    zz_pE    c;

    i = 0;
    do {
	c = coeff(x, i);
	r = inc(c);
	SetCoeff(x, i, c);
    } while (r == 1 && ++i <= max_deg);

    return (i > max_deg) ? CARRY : NO_CARRY;
}

// replaces x with next element in 'lexical ordering' of F_q

int inc(zz_pE& x)
{
    long    d = deg(x.modulus()), i;
    zz_pE   tee = to_zz_pE(zz_pX(1,1)), tee_to_i = to_zz_pE(zz_pX(0,1));

    i = 0;
    do {
	x += tee_to_i;
	tee_to_i *= tee;
    } while (coeff(x.LoopHole(),i) == 0 && ++i < d);

    return (i >= d) ? 1 : 0;
}

// evaluate polynomial f at x

zz_pE eval(const zz_pX& f, const zz_pE& x)
{
    static zz_pE  y,z;
    static zz_p   t;
    long    i;

    GetCoeff(t, f,deg(f));
    y=t;
    for (i = deg(f)-1; i >= 0; i--) {
	GetCoeff(t, f, i);
	mul(z, y, x);
	add(y, z, t);
    }

    return y;
}

ZZ eval(const ZZX& f, const ZZ& x)
{
    static ZZ  y,z;
    static ZZ   t;
    long    i;

    GetCoeff(t, f, deg(f));
    y=t;
    for (i = deg(f)-1; i >= 0; i--) {
        GetCoeff(t, f, i);
        mul(z, y, x);
        add(y, z, t);
    }
    return y;
}

// convert a long to a string

char *ltoa(long i)
{
    static char   result[20], *str;

    str = result + 19;
    *str-- = 0;
    do {
        *str-- = (i%10) + '0';
        i /= 10;
    } while (i != 0);

    return str+1;
}

// compute x^e

zz_pE operator^(const zz_pE& x, const int e)
{
    zz_pE   result, l;
    long    i;

    result = 1;
    if (e == 0) return result;
    if (e <  0) return 1/(x^(-e));

    // I = 2^i
    for (i = 0, l = x; (1<<i) <= e; i++) {
	if (e & (1<<i))
	    result *= l;
	l = sqr(l);
    }

    return result;
}

// compute x^e

zz_pEX operator^(const zz_pEX& x, const int e)
{
    zz_pEX  result, l;
    long    i;

    assert(e >= 0);

    result = 1;
    if (e == 0) return result;

    // I = 2^i
    for (i = 0, l = x; (1<<i) <= e; i++) {
	if (e & (1<<i))
	    result *= l;
	l = sqr(l);
    }

    return result;
}

#ifdef MAIN
int main(int argc, char **argv)
{
    zz_pX   pi_1, pi_2, a;

    get_modulus(pi_1, pi_2, a, 5, 2, 2);

    long    q = to_long(zz_pE::cardinality());
    cout << "q = " << q << endl;

    cout << "pi_1 = " << pi_1 << endl;
    cout << "pi_2 = " << pi_2 << endl;
    cout << "a = " << a << endl;
}
#endif // MAIN
