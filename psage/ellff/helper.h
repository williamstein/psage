#ifndef HELPER_H
#define HELPER_H

#include <NTL/lzz_p.h>
#include <NTL/lzz_pE.h>
#include <NTL/lzz_pX.h>
#include <NTL/lzz_pEX.h>
#include <NTL/ZZX.h>

NTL_CLIENT

#define NO_CARRY   0
#define CARRY      1

// return a modulus appropriate for F_q=F_{p^d}

void get_modulus(zz_pX& pi, int p, int d);

// return moduli pi_1,pi_2 for extensions F_q/F_p,F_{q^d}/F_q and zero of pi_1
void get_modulus(zz_pX& pi_1, zz_pX& pi_2, zz_pX& a, int p, int d1, int d2);

// setup F_{p^d} with canonical irreducible in F_p[x]

void init_NTL_ff(int p, int d, int precompute_inverses=1,
    int precompute_square_roots=1, int precompute_legendre_char=1,
    int precompute_pth_frobenius_map=1);
void init_NTL_ff(int p, int d1, int d2, int precompute_inverses,
    int precompute_square_roots, int precompute_legendre_char,
    int precompute_pth_frobenius_map);

// compare elements in F_p[x] ordered by degree then leading coefficient   

extern long operator<(zz_pX& f, zz_pX& g);
extern inline long operator<=(zz_pX& f, zz_pX& g);
extern inline long operator>( zz_pX& f, zz_pX& g);
extern inline long operator>=(zz_pX& f, zz_pX& g);

// compare elements of F_q 

extern long operator<(zz_pE& x, zz_pE& y);
extern inline long operator<=(zz_pE& f, zz_pE& g);
extern inline long operator>(zz_pE& f, zz_pE& g);
extern inline long operator>=(zz_pE& f, zz_pE& g);

// convert elt of F_p[x] to int using coefficients as base-p digits
// used to determine index of elt in table

extern unsigned long to_ulong(const zz_pX& x);
extern unsigned long to_ulong(zz_pE& x);

extern void from_ulong(unsigned long ul, zz_pX& x);
extern void from_ulong(unsigned long ul, zz_pE& x);

// increment polynomial as if coefficients were base-p digits of an
// integer.  corresponding sequence of polynomials is a so-called lexical
// ordering.

extern int inc(zz_pX& x, long max_deg);
extern int inc(zz_pEX& x, long max_deg);

// replaces x with next element in 'lexical ordering' of F_q

extern int inc(zz_pE& x);

// evaluate polynomial f at x

extern zz_pE eval(const zz_pX& f, const zz_pE& x);
extern ZZ eval(const ZZX& f, const ZZ& x);

// convert a long to a string

extern char *ltoa(long i);

// compute x^e

extern zz_pE operator^(const zz_pE& x, const int e);
extern zz_pEX operator^(const zz_pEX& x, const int e);

#endif	// HELPER_H
