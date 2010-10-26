/*********************************************************************

 (c) Copyright 2006-2010 Salman Baig and Chris Hall

 This file is part of ELLFF

 ELLFF is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 ELLFF is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

*********************************************************************/

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
