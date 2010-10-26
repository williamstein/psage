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
#include <string.h>
#include <math.h>
#include <NTL/ZZX.h>
#include <NTL/lzz_p.h>
#include <NTL/lzz_pE.h>
#include <NTL/lzz_pX.h>
#include <NTL/lzz_pEX.h>
#include <NTL/lzz_pXFactoring.h>
#include <NTL/lzz_pX.h>

#include "helper.h"
#include "lzz_pEExtra.h"

NTL_CLIENT

///////////////////////////////////////////////////////////////////////////

zz_pEExtraInfoPtr   zz_pEExtraInfo = NULL;

///////////////////////////////////////////////////////////////////////////

zz_pEExtraInfoT::zz_pEExtraInfoT(int precompute_inverses,
    int precompute_square_roots, int precompute_legendre_char,
    int precompute_pth_frobenius_map)
{
    int p = zz_p::modulus();
    q = to_long(zz_pE::cardinality());

    ref_count = 1;

    inv_table = precompute_inverses ? new invTable(q) : NULL;
    root_table = precompute_square_roots ? new rootTable(q) : NULL;
    legendre_table = precompute_legendre_char ? new legendreChar(q) : NULL;
    frob_map = precompute_pth_frobenius_map ? new frobeniusMap(q) : NULL;

    // precompute a non-square in F_q
    zz_pE   x, y;
    do {
	x = random_zz_pE();
    } while(x == 0 || legendre_char(x) == 1);
    non_square = x;

    // precompute image of basis 1,x^1,...,x^{d-1} under Frobenius
    frob_of_basis = new zz_pE[zz_pE::degree()];
    x = 0;
    SetCoeff(x.LoopHole(), 1, 1);
    frob_of_basis[0] = 1;
    for (int i = 1; i < zz_pE::degree(); i++) {
        power(y, x, i);
        power(frob_of_basis[i], y, p);
    }
}

zz_pEExtraInfoT::~zz_pEExtraInfoT()
{
    if (inv_table != NULL)	delete inv_table;
    if (root_table != NULL)	delete root_table;
    if (legendre_table != NULL)	delete legendre_table;
    if (frob_map != NULL)	delete frob_map;

    delete [] frob_of_basis;
}

void CopyPointer(zz_pEExtraInfoPtr& dst, zz_pEExtraInfoPtr src)
{
   if (src == dst) return;

   if (dst) {
      dst->ref_count--;

      if (dst->ref_count < 0) 
         Error("internal error: negative zz_pEExtraInfoT ref_count");

      if (dst->ref_count == 0) delete dst;
   }

   if (src) {
      if (src->ref_count == NTL_MAX_LONG) 
         Error("internal error: zz_pEExtraInfoT ref_count overflow");

      src->ref_count++;
   }

   dst = src;
}

///////////////////////////////////////////////////////////////////////////
// table of square roots
///////////////////////////////////////////////////////////////////////////

zz_pEExtraInfoT::rootTable::rootTable(long q)
{
    zz_pE  x, x_sqr;
    x = 0;
    table = new zz_pE[q];
    do {
	x_sqr = sqr(x);
	table[to_ulong(x_sqr)] = x;
    } while (inc(x) == NO_CARRY);
}

zz_pEExtraInfoT::rootTable::~rootTable()
{
    delete table;
}

///////////////////////////////////////////////////////////////////////////
// compute square root
///////////////////////////////////////////////////////////////////////////

zz_pE zz_pEExtraInfoT::square_root(zz_pE& x)
{
    if (root_table != NULL)
	return root_table->table[to_ulong(x)];

    if (q % 4 == 3)
	return power(x, (q+1)/4);

    zz_pE    r, s, v, w;

    // compute square root using Tonelli-Shanks
    long    e;
    int     ord_2;
    for (e = (q-1)/2, ord_2 = 1; e % 2 == 0; ord_2++)
	e /= 2;

    v = power(non_square, e);
    r = power(x, (e+1)/2);

    int     i;
    while (1) {
	s = r*r/x;
	for (i = 0; s != 1; i++)
	    s = sqr(s);
	assert(i < ord_2);

	if (i == 0)
	    return r;

	w = v;
	while (i++ < ord_2-1)
	    w = sqr(w);
	r = r*w;
    }

    assert(r*r == x);

    return r;
}

///////////////////////////////////////////////////////////////////////////
// compute Legendre character
///////////////////////////////////////////////////////////////////////////

int zz_pEExtraInfoT::legendre_char(zz_pE& x)
{
    if (legendre_table != NULL)
	return legendre_table->table[to_ulong(x)];

    static zz_pE   y;

    assert(x != 0);
    y = power(x, (q-1)/2);
    assert(y == 1 || y == -1);

    return (y == 1);
}

///////////////////////////////////////////////////////////////////////////
// table of Legendre character values
///////////////////////////////////////////////////////////////////////////

zz_pEExtraInfoT::legendreChar::legendreChar(long q)
{
    zz_pE  x, y;
    table = (char*)malloc(sizeof(char)*q);
    memset(table, 0, sizeof(char)*q);
    x = 0;
    while (inc(x) == NO_CARRY) {  // loop over non-zero values
	y = x*x;
	table[to_ulong(y)] = 1;	// 1 square, 0 non-square (or 0)
    }
}

zz_pEExtraInfoT::legendreChar::~legendreChar()
{
    free(table);
}

///////////////////////////////////////////////////////////////////////////
// table of multiplicative inverses
///////////////////////////////////////////////////////////////////////////

zz_pEExtraInfoT::invTable::invTable(long q)
{
    zz_pE  x;
    table = new zz_pE[q];
    while (inc(x) == NO_CARRY)    // loop over non-zero values
	table[to_ulong(x)] = 1/x;
}

zz_pEExtraInfoT::invTable::~invTable()
{
    delete table;
}

///////////////////////////////////////////////////////////////////////////
// pth Frobenius map
///////////////////////////////////////////////////////////////////////////

// y = x^p

void zz_pEExtraInfoT::frobenius(zz_pE&x, zz_pE&y)
{
    zz_pE  z;
    zz_pX  f = x.LoopHole();
    for (int i = 0; i <= deg(f); i++)
        z += coeff(f, i) * frob_of_basis[i];
    y = z;
}

unsigned long zz_pEExtraInfoT::frobenius(unsigned long x)
{
    if (frob_map != NULL)
        return frob_map->map[x];

    zz_pE  _x, _y;

    from_ulong(x, _x);
    frobenius(_x, _y);
    return to_ulong(_y);
}

zz_pEExtraInfoT::frobeniusMap::frobeniusMap(long q)
{
    zz_pE   x, y;
    int     p = zz_p::modulus();

    map = (unsigned long*)malloc(sizeof(unsigned long)*q);
    do {
        power(y, x, p);
        map[to_ulong(x)] = to_ulong(y);
    } while (inc(x) == NO_CARRY);
}

zz_pEExtraInfoT::frobeniusMap::~frobeniusMap()
{
    free(map);
}

///////////////////////////////////////////////////////////////////////////

zz_pEExtraContext::zz_pEExtraContext(int precompute_inverses,
    int precompute_square_roots, int precompute_legendre_char,
    int precompute_pth_frobenius_map)
{
    ptr = new zz_pEExtraInfoT(precompute_inverses,
                precompute_square_roots,
		precompute_legendre_char,
                precompute_pth_frobenius_map);
}

zz_pEExtraContext::zz_pEExtraContext(const zz_pEExtraContext& a)
{
   ptr = NULL;
   CopyPointer(ptr, a.ptr);
}

zz_pEExtraContext& zz_pEExtraContext::operator=(const zz_pEExtraContext& a)
{
   CopyPointer(ptr, a.ptr);
   return *this;
}

zz_pEExtraContext::~zz_pEExtraContext()
{
   CopyPointer(ptr, NULL);
}

void zz_pEExtraContext::save()
{
   CopyPointer(ptr, zz_pEExtraInfo);
}

void zz_pEExtraContext::restore() const
{
   CopyPointer(zz_pEExtraInfo, ptr);
}

