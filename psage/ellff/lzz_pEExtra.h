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

#ifndef LZZ_PEEXTRA_H
#define LZZ_PEEXTRA_H

#include <NTL/lzz_p.h>
#include <NTL/lzz_pE.h>
#include <NTL/lzz_pX.h>
#include <NTL/lzz_pEX.h>
#include <NTL/ZZX.h>

NTL_CLIENT

class zz_pEExtraInfoT {
    private:
    	zz_pEExtraInfoT();	// do not use

	zz_pE	non_square;

    public:
	// table of multiplicative inverses
	class invTable {
	    public:
		invTable(long q);
		~invTable();

		zz_pE    *table;
	};

	// table of square roots
	class rootTable {
	    public:
		rootTable(long q);
		~rootTable();

		zz_pE    *table;
	};

	zz_pE square_root(zz_pE& x);

	// Legendre character
	class legendreChar {
	    public:
		legendreChar(long q);
		~legendreChar();

		char    *table;
	};

	int legendre_char(zz_pE& x);

        // Frobenius map
        class frobeniusMap {
            public:
                frobeniusMap(long q);
                ~frobeniusMap();

                unsigned long    *map;
        };

        void frobenius(zz_pE& x, zz_pE& y);
        unsigned long frobenius(unsigned long x);

        zz_pE   *frob_of_basis;

    public:
	long ref_count;		// for garbage collection

	long q;

	zz_pEExtraInfoT(int precompute_inverses,
	    int precompute_square_roots,
            int precompute_legendre_char,
            int precompute_pth_frobenius_map);

	~zz_pEExtraInfoT();

	invTable	*inv_table;
	rootTable    	*root_table;
	legendreChar 	*legendre_table;
        frobeniusMap    *frob_map;
};

typedef zz_pEExtraInfoT* zz_pEExtraInfoPtr;

class zz_pEExtraContext {
    private:
	zz_pEExtraInfoPtr  ptr;

    public:
	void save();
	void restore() const;

	zz_pEExtraContext() { ptr = NULL; }
	zz_pEExtraContext(int precompute_inverses,
	    int precompute_square_roots,
            int precompute_legendre_char,
            int precompute_pth_frobenius_map);

	zz_pEExtraContext(const zz_pEExtraContext&);

	zz_pEExtraContext& operator=(const zz_pEExtraContext&);

	~zz_pEExtraContext();
};

extern zz_pEExtraInfoPtr   zz_pEExtraInfo;

#endif	// LZZ_PEEXTRA_H
