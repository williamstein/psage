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

#ifndef ELL_H
#define ELL_H

// meant for elliptic curves in char != 2,3

class ell_pEInfoT {
  private:
    ell_pEInfoT();	// do not use
    ell_pEInfoT(const ell_pEInfoT&);	// do not use
    void operator=(const ell_pEInfoT&);	// do not use

  public:
    long    ref_count;	// for garbage collection

    ell_pEInfoT(const zz_pE& new_a4, const zz_pE& new_a6);

    ~ell_pEInfoT();

    zz_pE    a4, a6;
    long     q;
};

typedef ell_pEInfoT  *ell_pEInfoPtr;

class ell_pEContext {
  private:
    ell_pEInfoT   *ptr;

  public:
    void save();
    void restore() const;

    ell_pEContext() { ptr = NULL; }
    ell_pEContext(const zz_pE& a4, const zz_pE& a6);

    ell_pEContext(const ell_pEContext&);

    ell_pEContext& operator=(const ell_pEContext&);

    ~ell_pEContext();
};

extern ell_pEInfoPtr   ell_pEInfo;

class ell_pE {
  public:
    //ell_pE();

    // set curve to E/F_q : y^2 = x^3 + a4*x + a6
    static void init(zz_pE& a4, zz_pE& a6);

    static long order();
};

class ellpoint_pE {
  private:

  public:
    zz_pE   x, y;
    int     identity_f;  

    ellpoint_pE();
    ellpoint_pE(const zz_pE& x, const zz_pE& y);

    inline int on_curve()
	{ return sqr(y) == ((sqr(x)+ell_pEInfo->a4)*x + ell_pEInfo->a6); }

    inline void set() { identity_f = 1; }
    inline void set(const zz_pE& _x, const zz_pE& _y)
	{ identity_f = 0; x = _x; y = _y; }

    long order();
};

void neg(ellpoint_pE& Q, const ellpoint_pE& P);

void add(ellpoint_pE& R, const ellpoint_pE& P, const ellpoint_pE& Q);
void sub(ellpoint_pE& R, const ellpoint_pE& P, const ellpoint_pE& Q);

ellpoint_pE operator-(const ellpoint_pE& P);

ellpoint_pE operator+(const ellpoint_pE& P, const ellpoint_pE& Q);
ellpoint_pE operator-(const ellpoint_pE& P, const ellpoint_pE& Q);

ellpoint_pE& operator+=(ellpoint_pE& P, const ellpoint_pE& Q);
ellpoint_pE& operator-=(ellpoint_pE& P, const ellpoint_pE& Q);

ellpoint_pE operator*(int m, const ellpoint_pE& P);

int is_zero(const ellpoint_pE& P);

long operator==(const ellpoint_pE& P, const ellpoint_pE& Q);
long operator<( ellpoint_pE& P, ellpoint_pE& Q);
long operator>( ellpoint_pE& P, ellpoint_pE& Q);
long operator<=(ellpoint_pE& P, ellpoint_pE& Q);
long operator>=(ellpoint_pE& P, ellpoint_pE& Q);
ostream& operator<<(ostream& s, const ellpoint_pE& P);

unsigned long* to_ulong(ellpoint_pE& P);

// routine for freeing any large chunks of memory allocated
void ell_cleanup();

#endif // ELL_H

