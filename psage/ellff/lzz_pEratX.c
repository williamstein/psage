//
// rational function field F_q(t)
//

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
#include <NTL/lzz_pX.h>

#include "lzz_pEratX.h"

NTL_CLIENT

///////////////////////////////////////////////////////////////////////////
// basic rational elements of F_q(t)

zz_pEratX::zz_pEratX(const zz_pEX& num, const zz_pEX& den)
{
    init(num, den);
}

zz_pEratX::zz_pEratX(const zz_pEX& num)
{
    zz_pEX   den;
    den = 1;
    init(num, den);
}

zz_pEratX::zz_pEratX()
{
    init();
}

void zz_pEratX::init()
{
    num = 0;
    den = 1;
}

void zz_pEratX::init(const zz_pEX& num)
{
    zz_pEX   one;
    one = 1;
    init(num, one);
}

void zz_pEratX::init(const zz_pEX& num, const zz_pEX& den)
{
    assert(!IsZero(den));
    this->num = num;
    this->den = den;
    zz_pEX	d = GCD(num, den);
    this->num /= d;
    this->den /= d;
    zz_pE	c = LeadCoeff(this->den);
    this->num /= c;
    this->den /= c;
}

// g(t) = t^deg(f)*f(1/t)

void zz_pEratX::flip(zz_pEX& g, const zz_pEX& f)
{
    g = 0;
    for (int i = 0; i <= deg(f); i++)
	g = (g<<1) + coeff(f, i);
}

// perform involution t-->1/t

void zz_pEratX::inv_t()
{
    zz_pEX  new_num, new_den;
    int	deg_num = deg(num), deg_den = deg(den);
    flip(new_num, num);
    flip(new_den, den);
    if (deg_num < deg_den)
	new_num <<= deg_den - deg_num;
    else if (deg_num > deg_den)
	new_den <<= deg_num - deg_den;
    num = new_num;
    den = new_den;
}

int IsZero(const zz_pEratX& f)
{
    return (IsZero(f.num));
}

int IsConstant(const zz_pEratX& f)
{
    if (IsZero(f.num))
	return 1;
    if (deg(f.num) != 0)
	return 0;
    if (deg(f.den) != 0)
	return 0;
    return 1;
}

int deg(const zz_pEratX& f)
{
    return IsZero(f.num) ? -1 : deg(f.num) - deg(f.den);
}

long operator==(const zz_pEratX& a, const zz_pEratX& b)
{
    return IsZero(a.num*b.den - b.num*a.den);
}

// l*x

zz_pEratX operator*(const long l, const zz_pEratX& x)
{
    zz_pEratX	result(x.num, x.den);
    result.num *= l;
    return result;
}

// a*b

zz_pEratX operator*(const zz_pEratX& a, const zz_pEratX& b)
{
    zz_pEratX	result(a.num*b.num, a.den*b.den);
    return result;
}

// a/b

zz_pEratX operator/(const zz_pEratX& a, const zz_pEX& b)
{
    static zz_pEratX	result;
    assert(!IsZero(b));
    result.init(a.num, a.den*b);
    return result;
}

// a/b

zz_pEratX operator/(const zz_pEratX& a, const zz_pEratX& b)
{
    static zz_pEratX	result;
    assert(!IsZero(b));
    result.init(a.num*b.den, a.den*b.num);
    return result;
}

// a+b

zz_pEratX operator+(const zz_pEratX& a, const zz_pE& b)
{
    static zz_pEratX  result;
    result.init(a.num + b*a.den, a.den);
    return result;
}

// a+b

zz_pEratX operator+(const zz_pEratX& a, const zz_pEratX& b)
{
    static zz_pEratX  result;
    result.init(a.num*b.den + b.num*a.den, a.den*b.den);
    return result;
}

// a+b

zz_pEratX operator-(const zz_pEratX& a, const zz_pEratX& b)
{
    static zz_pEratX  result;
    result.init(a.num*b.den - b.num*a.den, a.den*b.den);
    return result;
}

// compute x^e

zz_pEratX operator^(const zz_pEratX& x, const int e)
{
    zz_pEX   one;
    one = 1;

    zz_pEratX  result(one), l;
    long    i;

    assert(e >= 0);

    if (e == 0) return result;

    // I = 2^i
    for (i = 0, l = x; (1<<i) <= e; i++) {
	if (e & (1<<i))
	    result = result * l;
	l = l * l;
    }

    return result;
}

// compose rational functions (f,g)-->f(g)

zz_pEratX eval(const zz_pEratX& f, const zz_pEratX& g)
{
    static zz_pEratX   r;

    r.init();
    if (IsZero(f.num))
	return r;
    
    return eval(f.num, g) / eval(f.den, g);
}

// compose polynomial with rational function (f,g)-->f(g)

zz_pEratX eval(const zz_pEX& f, const zz_pEratX& g)
{
    static zz_pEratX  r, z;
    static zz_pEX t;
    static zz_pE  c;
    long    i;

    c = coeff(f, deg(f));
    t = c;
    r.init(t);
    for (i = deg(f)-1; i >= 0; i--) {
	//r = r*g + coeff(f,i);
	c = coeff(f, i);
	z = r*g;
	r = z+c;
    }

    return r;
}

