#ifndef LZZ_PERATX_H
#define LZZ_PERATX_H

#include <NTL/lzz_p.h>
#include <NTL/lzz_pE.h>
#include <NTL/lzz_pX.h>
#include <NTL/lzz_pEX.h>
#include <NTL/ZZX.h>

NTL_CLIENT

class zz_pEratX {
    private:
	void flip(zz_pEX& g, const zz_pEX& f);

    public:
    	zz_pEratX(const zz_pEX& num, const zz_pEX& den);
    	zz_pEratX(const zz_pEX& num);
    	zz_pEratX();

	void init();
	void init(const zz_pEX& num);
	void init(const zz_pEX& num, const zz_pEX& den);

	void inv_t();

    	zz_pEX	num, den;
};

extern int IsConstant(const zz_pEratX& f);
extern int deg(const zz_pEratX& f);

extern long operator==(const zz_pEratX& a, const zz_pEratX& b);

extern zz_pEratX operator+(const zz_pEratX& a, const zz_pEratX& b);
extern zz_pEratX operator-(const zz_pEratX& a, const zz_pEratX& b);
extern zz_pEratX operator*(const long l, const zz_pEratX& x);
extern zz_pEratX operator*(const zz_pEratX& a, const zz_pEratX& b);
extern zz_pEratX operator/(const zz_pEratX& a, const zz_pEX& b);
extern zz_pEratX operator^(const zz_pEratX& x, const int e);

extern zz_pEratX eval(const zz_pEratX& f, const zz_pEratX& g);
extern zz_pEratX eval(const zz_pEX& f, const zz_pEratX& g);

#endif	// LZZ_PERATX_H
