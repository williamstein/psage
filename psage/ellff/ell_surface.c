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

//
// C++ elliptic surface class
//
// this code offers basic routines for elliptic curves over a function
// field F_q(t), where F_q is a finite field.  currently it is limited to
// curves with a model of the form y^2=x^3+a4*x+a6, which excludes
// characteristic two and some curves in characteristic three.  because the
// latter two characteristics present other difficulties for the L-function
// code, we felt it was ok to make this exclusion.
//

#include <assert.h>
#include <string.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/lzz_p.h>
#include <NTL/lzz_pE.h>
#include <NTL/lzz_pX.h>
#include <NTL/lzz_pEX.h>
#include <NTL/lzz_pEXFactoring.h>

NTL_CLIENT

#include "ell_surface.h"
#include "helper.h"

///////////////////////////////////////////////////////////////////////////
// minimal Weierstrass model over affine ring
//  - can be used for F_q[t] or local ring O_pi

ell_surfaceInfoT::affine_model::affine_model()
{
    init();
}

// determine kodaira type and other info of singular fiber over pi
//
//  - should be called exactly *once* because changes epsilon and various
//  divisors describing bad reduction

void ell_surfaceInfoT::affine_model::kodaira(const zz_pEX& pi)
{
    zz_pEX  r;

    // good reduction
    rem(r, disc, pi);
    if (!IsZero(r))
	return;

    // j = 0 (mod pi): I_0*,II,IV,IV*,II*
    rem(r, j.num, pi);
    if (IsZero(r)) {
	static zz_pEX  t;
	int     ord_pi = 0, local_epsilon = 0;

	// determine ord_pi(disc)
	A *= pi;
	t = disc;
	do {
	    t /= pi;
	    ord_pi++;
	    rem(r, t, pi);
	} while (IsZero(r));

	// explicit Kodaira type
	assert(ord_pi > 0 && ord_pi < 12 && ord_pi % 2 == 0);
	switch (ord_pi) {
	    case  2 : II      *= pi; local_epsilon = -1; break;
	    case 10 : II_star *= pi; local_epsilon = -1; break;
	    case  4 : IV      *= pi; local_epsilon = -3; break;
	    case  8 : IV_star *= pi; local_epsilon = -3; break;
	    case  6 : I_star  *= pi; local_epsilon = -1; break;
	    default : assert(0);
	}

	// determine contribution to sign of functional equation
	switch (local_epsilon) {
	    case -1 : epsilon *= (q%4==1 || deg(pi)%2==0) ? +1 : -1; break;
	    case -3 : epsilon *= (q%3==1 || deg(pi)%2==0) ? +1 : -1; break;
	    default : assert(0);
	}

	return;
    }

    // j = 1728 (mod pi): I_0*,III,III*
    rem(r, j.num-1728*j.den, pi);
    if (IsZero(r)) {
	static zz_pEX  t;
	int     ord_pi = 0, local_epsilon = 0;

	// determine ord_pi(disc)
	A *= pi;
	t = disc;
	do {
	    t /= pi;
	    ord_pi++;
	    rem(r, t, pi);
	} while (IsZero(r));

	// explicit Kodaira type
	assert(ord_pi > 0 && ord_pi < 12 && ord_pi % 3 == 0);
	switch (ord_pi) {
	    case  3 : III_star *= pi; local_epsilon = -2; break;
	    case  9 : III      *= pi; local_epsilon = -2; break;
	    case  6 : I_star   *= pi; local_epsilon = -1; break;
	    default : assert(0);
	}

	// contribution to sign of functional equation
	switch (local_epsilon) {
	    case -1 : epsilon *= (q%4==1 || deg(pi)%2==0) ? +1 : -1; break;
	    case -2 : epsilon *= (q%8<=3 || deg(pi)%2==0) ? +1 : -1; break;
	    default : assert(0);
	}

	return;
    }

    // j = infty (mod pi): I_n,I_n^*
    rem(r, j.den, pi);
    if (IsZero(r)) {
	// determine if multiplicate or additive
	rem(r, a6, pi);
	if (IsZero(r)) {	// add've
	    A *= pi;
	    I_star *= pi;
	    if (q % 4 == 3 && deg(pi) % 2 == 1)
		epsilon *= -1;
	} else {		// mult've
	    // determine if split or non-split
	    //  - split <==> 6*a6 (mod pi) is a square

	    ZZ    q_to_d, e;
	    q_to_d = 1;
	    for (int j = deg(pi); j > 0; j--)
		q_to_d *= q;
	    e = (q_to_d - 1) / 2;

	    zz_pEXModulus  mod_pi(pi);
	    zz_pEX  tmp;
	    rem(tmp, 6*a6, mod_pi);
	    PowerMod(r, tmp, e, mod_pi);
	    assert(r == 1 || r == -1);

	    if (r == 1) {
		M_sp *= pi;
		epsilon *= -1;
	    } else
		M_ns *= pi;
	}

	return;
    }

    // I_0^* : ord_pi(j) = ord_pi(j-1728) = 0
    A *= pi;
    I_star *= pi;
    if (q%4==3 && deg(pi)%2==1)
	epsilon *= -1;
}

// initialize

void ell_surfaceInfoT::affine_model::init()
{
    q = to_long(zz_pE::cardinality());

    I_star   = 1;	// I_n^*
    II       = 1;	// II
    II_star  = 1;	// II*
    III      = 1;	// III
    III_star = 1;	// III*
    IV       = 1;	// IV
    IV_star  = 1;	// IV*

    M_sp = 1;		// split multiplicative reduction
    M_ns = 1;		// non-split multiplicative reduction
    A    = 1;		// additive reduction

    epsilon = 1;	// epsilon factor

    disc = 0;
}

void ell_surfaceInfoT::affine_model::init(const zz_pEratX& rat_a4, const zz_pEratX& rat_a6)
{
    init();

    // step 1: clear denominators by map (a4,a6)-->(a4*d^4,a6*d^6) to get a
    //   (not necessarily minimal) Weierstrass over F_q[t]
    static zz_pEX  d;
    d = rat_a4.den * rat_a6.den / GCD(rat_a4.den, rat_a6.den);

    a4 = rat_a4.num * (d^4) / rat_a4.den;
    a6 = rat_a6.num * (d^6) / rat_a6.den;
    assert(!IsZero(a4) && !IsZero(a6));

    // step 2: calculate j-invariant for current Weierstrass model
    disc = 4*a4*a4*a4 + 27*a6*a6;
    j.init(6912*a4*a4*a4, disc);
    assert(!IsConstant(j));

    // step 3: divide a4,a6 by max'l degree polynomial d so d^4|a4 and d^6|a6
    //  - for each prime pi|a4, keep dividing a4,a6 by pi^4,pi^6 until
    //    one or both are not divisible
    //  - result is minimal Weierstrass model over F_q[t]
    vec_pair_zz_pEX_long factors;
    CanZass(factors, a4 / LeadCoeff(a4));
    for (int i = 0; i < factors.length(); i++) {
	static zz_pEX	pi, pi_to_4, pi_to_6;

	pi = factors[i].a;
	pi_to_4 = pi^4;
	pi_to_6 = pi^6;
	while (1) {
	    static zz_pEX   r_a4, r_a6;

	    rem(r_a4, a4, pi_to_4);
	    rem(r_a6, a6, pi_to_6);
	    if (!IsZero(r_a4) || !IsZero(r_a6))
		break;
	    a4 /= pi_to_4;
	    a6 /= pi_to_6;
	}
    }

    // recalculate discriminant for new model
    disc = 4*a4*a4*a4 + 27*a6*a6;
}

void ell_surfaceInfoT::affine_model::init(const zz_pEX& a4, const zz_pEX& a6)
{
    zz_pEratX rat_a4(a4), rat_a6(a6);
    init(rat_a4, rat_a6);
}

///////////////////////////////////////////////////////////////////////////
// elliptic surface parameter class
//   - we keep track of the coeffs for a Weierstrass model over F_q[t],
//   the coeffs for a minimal Weierstrass model over F_q[u] for u=1/t,
//   and finer information about the type of bad reduction, e.g. split vs
//   non-split mult've, and Kodaira type in the case of add've reduction
///////////////////////////////////////////////////////////////////////////

// constructor

ell_surfaceInfoT::ell_surfaceInfoT(const zz_pEratX& a4, const zz_pEratX& a6)
{
    static zz_pEX	pi;

    ref_count = 1;

    q = to_long(zz_pE::cardinality());

    // compute minimal Weierstrass model over F_q[t]
    finite_model.init(a4, a6);

    // determine Kodaira type of bad reduction for minimal Weierstrass model
    vec_pair_zz_pEX_long factors;
    CanZass(factors, finite_model.disc / LeadCoeff(finite_model.disc));
    for (int i = 0; i < factors.length(); i++) {
	pi = factors[i].a;
	finite_model.kodaira(pi);
    }

    // compute minimal Weierstrass model about t=infty
    zz_pEratX	a4_of_u = a4, a6_of_u = a6;
    a4_of_u.inv_t();
    a6_of_u.inv_t();
    infinite_model.init(a4_of_u, a6_of_u);

    // determine Kodaira type of bad reduction
    pi = 1;
    pi <<= 1;
    infinite_model.kodaira(pi);

    // compute invariants for L-function
    sign = finite_model.epsilon * infinite_model.epsilon;
    deg_L  = 2*(deg(finite_model.A)    + deg(infinite_model.A))
	      + deg(finite_model.M_sp) + deg(infinite_model.M_sp)
	      + deg(finite_model.M_ns) + deg(infinite_model.M_ns) - 4;
}

// destructor

ell_surfaceInfoT::~ell_surfaceInfoT()
{
}

zz_pEX ell_surfaceInfoT::finite_A()
{
    return finite_model.A;
}

int ell_surfaceInfoT::infinite_A()
{
    if (deg(infinite_model.A) > 0 && IsZero(coeff(infinite_model.A, 0)))
        return 1;
    return 0;
}

zz_pEX ell_surfaceInfoT::finite_M()
{
    return finite_model.M_sp * finite_model.M_ns;
}

int ell_surfaceInfoT::infinite_M()
{
    if (deg(infinite_model.M_sp) > 0 && IsZero(coeff(infinite_model.M_sp, 0)))
        return 1;
    if (deg(infinite_model.M_ns) > 0 && IsZero(coeff(infinite_model.M_ns, 0)))
        return 1;
    return 0;
}

// copied from NTL

static void CopyPointer(ell_surfaceInfoPtr& dst, ell_surfaceInfoPtr src)
{
   if (src == dst) return;

   if (dst) {
      dst->ref_count--;

      if (dst->ref_count < 0) 
         Error("internal error: negative ell_surfaceContext ref_count");

      if (dst->ref_count == 0) delete dst;
   }

   if (src) {
      if (src->ref_count == NTL_MAX_LONG) 
         Error("internal error: ell_surfaceContext ref_count overflow");

      src->ref_count++;
   }

   dst = src;
}

////////////////////////////////////////////////////////////
// context switching class, shamelessly copied from NTL
////////////////////////////////////////////////////////////

// info for current elliptic surface
ell_surfaceInfoPtr   ell_surfaceInfo = NULL;

ell_surfaceContext::ell_surfaceContext(const zz_pEratX& a4, const zz_pEratX& a6)
{
   ptr = new ell_surfaceInfoT(a4, a6);
}

ell_surfaceContext::ell_surfaceContext(const ell_surfaceContext& a)
{
   ptr = NULL;
   CopyPointer(ptr, a.ptr);
}

ell_surfaceContext& ell_surfaceContext::operator=(const ell_surfaceContext& a)
{
   CopyPointer(ptr, a.ptr);
   return *this;
}

ell_surfaceContext::~ell_surfaceContext()
{
   CopyPointer(ptr, NULL);
}

void ell_surfaceContext::save()
{
   CopyPointer(ptr, ell_surfaceInfo);
}

void ell_surfaceContext::restore() const
{
   CopyPointer(ell_surfaceInfo, ptr);
}

///////////////////////////////////////////////////////////////////////////

void ell_surface::init(const zz_pEratX& a4, const zz_pEratX& a6)
{
    ell_surfaceContext c(a4, a6);

    c.restore();
}

void ell_surface::init(const zz_pEX& a4, const zz_pEX& a6)
{
    zz_pEX	one;
    one = 1;
    zz_pEratX	rat_a4(a4, one), rat_a6(a6, one);
    init(rat_a4, rat_a6);
}

void ell_surface::init(const zz_pEratX& a1, const zz_pEratX& a2,
    const zz_pEratX& a3, const zz_pEratX& a4, const zz_pEratX& a6)
{
    zz_pEratX   new_a4, new_a6;
    zz_pEX	d;
    d = 48;
    new_a4 = (48*a4 - (a1^4) - 8*(a1^2)*a2 + 24*a1*a3 - 16*(a2^2)) / d;
    d = 864;
    new_a6 = (-72*a4*(a1^2) - 288*a4*a2 + 864*a6 + (a1^6)
	    + 12*a2*(a1^4) - 36*a3*(a1^3) + 48*(a2^2)*(a1^2)
	    - 144*a3*a2*a1 + 64*(a2^3) + 216*(a3^2)) / d;
    init(new_a4, new_a6);
}

//////////////////////////////

void to_ZZ_pX(const zz_pX& src, ZZ_pX& dst)
{
    ZZ_p    c_;

    //cout << "to_ZZ_px: src = " << src << endl;
    clear(dst);
    for (int i = 0; i <= deg(src); i++) {
        long    c = rep(coeff(src, i));
        c_ = to_ZZ_p(c);
        //cout << "c = " << c << ", c_ = " << c_ << ", modulus(c_) = " << c_.modulus() << endl;
        SetCoeff(dst, i, c_);
    }
    //cout << "to_ZZ_px: dst = " << dst << endl;
}

void to_ZZ_pEX(const zz_pEX& src, ZZ_pEX& dst)
{
    ZZ   p;
    p = zz_p::modulus();
    ZZ_p::init(p);

    // cout << "src = " << src << endl;
    clear(dst);
    for (int i = 0; i <= deg(src); i++) {
        zz_pE   c = coeff(src, i);

        if (i == 0) {
            zz_pX   m = c.modulus();
            ZZ_pX   M;
            to_ZZ_pX(m, M);
            ZZ_pE::init(M);
        }

        ZZ_pE   C;
        to_ZZ_pX(rep(c), C.LoopHole());
        SetCoeff(dst, i, C);
    }
    // cout << "dst = " << dst << endl;
}

void get_an(ZZ_pEX& a4, ZZ_pEX& a6)
{
    ell_surfaceInfoT    *info = ell_surface::getSurfaceInfo();

    to_ZZ_pEX(info->finite_model.a4, a4);
    to_ZZ_pEX(info->finite_model.a6, a6);
}

void get_reduction(ZZ_pEX **divisors, int finite)
{
    ell_surfaceInfoT    *info = ell_surface::getSurfaceInfo();

    if (finite) {
        to_ZZ_pEX(info->finite_model.M_sp,     *(divisors[0]));
        to_ZZ_pEX(info->finite_model.M_ns,     *(divisors[1]));
        to_ZZ_pEX(info->finite_model.I_star,   *(divisors[2]));
        to_ZZ_pEX(info->finite_model.II,       *(divisors[3]));
        to_ZZ_pEX(info->finite_model.II_star,  *(divisors[4]));
        to_ZZ_pEX(info->finite_model.III,      *(divisors[5]));
        to_ZZ_pEX(info->finite_model.III_star, *(divisors[6]));
        to_ZZ_pEX(info->finite_model.IV,       *(divisors[7]));
        to_ZZ_pEX(info->finite_model.IV_star,  *(divisors[8]));
    } else {
        to_ZZ_pEX(info->infinite_model.M_sp,     *(divisors[0]));
        to_ZZ_pEX(info->infinite_model.M_ns,     *(divisors[1]));
        to_ZZ_pEX(info->infinite_model.I_star,   *(divisors[2]));
        to_ZZ_pEX(info->infinite_model.II,       *(divisors[3]));
        to_ZZ_pEX(info->infinite_model.II_star,  *(divisors[4]));
        to_ZZ_pEX(info->infinite_model.III,      *(divisors[5]));
        to_ZZ_pEX(info->infinite_model.III_star, *(divisors[6]));
        to_ZZ_pEX(info->infinite_model.IV,       *(divisors[7]));
        to_ZZ_pEX(info->infinite_model.IV_star,  *(divisors[8]));
    }
}

void get_disc(ZZ_pEX& disc, int finite)
{
    ell_surfaceInfoT    *info = ell_surface::getSurfaceInfo();

    if (finite)
        to_ZZ_pEX(info->finite_model.disc, disc);
    else
        to_ZZ_pEX(info->infinite_model.disc, disc);
}

void get_j_invariant(ZZ_pEX& j_num, ZZ_pEX& j_den)
{
    ell_surfaceInfoT    *info = ell_surface::getSurfaceInfo();

    to_ZZ_pEX(info->finite_model.j.num, j_num);
    to_ZZ_pEX(info->finite_model.j.den, j_den);
}

#if 0
void get_finite_M(ell_surfaceInfoT *info, ZZ_pEX& M_sp, ZZ_pEX& M_ns)
{
    to_ZZ_pEX(info->finite_model.M_sp, M_sp);
    to_ZZ_pEX(info->finite_model.M_ns, M_ns);
}

void get_finite_A(ell_surfaceInfoT *info, ZZ_pEX& A, ZZ_pEX& I_star, ZZ_pEX& II, ZZ_pEX& II_star, ZZ_pEX& III, ZZ_pEX& III_star, ZZ_pEX& IV, ZZ_pEX& IV_star)
{
    to_ZZ_pEX(info->finite_model.A, A);
    to_ZZ_pEX(info->finite_model.I_star,   I_star);
    to_ZZ_pEX(info->finite_model.II,       II);
    to_ZZ_pEX(info->finite_model.II_star,  II_star);
    to_ZZ_pEX(info->finite_model.III,      III);
    to_ZZ_pEX(info->finite_model.III_star, III_star);
    to_ZZ_pEX(info->finite_model.IV,       IV);
    to_ZZ_pEX(info->finite_model.IV_star,  IV_star);
}
#endif

#if 0
void get_infinite_an(ell_surfaceInfoT *info, ZZ_pEX& a4, ZZ_pEX& a6)
{
    to_ZZ_pEX(info->infinite_model.a4, a4);
    to_ZZ_pEX(info->infinite_model.a6, a6);
}
#endif

#if 0
void get_infinite_M(ell_surfaceInfoT *info, ZZ_pEX& M_sp, ZZ_pEX& M_ns)
{
    to_ZZ_pEX(info->infinite_model.M_sp, M_sp);
    to_ZZ_pEX(info->infinite_model.M_ns, M_ns);
}

void get_infinite_A(ell_surfaceInfoT *info, ZZ_pEX& A, ZZ_pEX& I_star, ZZ_pEX& II, ZZ_pEX& II_star, ZZ_pEX& III, ZZ_pEX& III_star, ZZ_pEX& IV, ZZ_pEX& IV_star)
{
    to_ZZ_pEX(info->infinite_model.A, A);
    to_ZZ_pEX(info->infinite_model.I_star,   I_star);
    to_ZZ_pEX(info->infinite_model.II,       II);
    to_ZZ_pEX(info->infinite_model.II_star,  II_star);
    to_ZZ_pEX(info->infinite_model.III,      III);
    to_ZZ_pEX(info->infinite_model.III_star, III_star);
    to_ZZ_pEX(info->infinite_model.IV,       IV);
    to_ZZ_pEX(info->infinite_model.IV_star,  IV_star);
}
#endif

#if 0
void get_infinite_disc(ell_surfaceInfoT *info, ZZ_pEX& disc)
{
    to_ZZ_pEX(info->infinite_model.disc, disc);
}
#endif

//////////////////////////////

#ifdef MAIN
int main(int argc, char **argv)
{
    long    f = 0, p = 5, d = 1;

    if (argc < 2 || argc > 4) {
	fprintf(stderr, "usage: %s f [p [d]]\n", argv[0]);
	return -1;
    }
    f = atoi(argv[1]);
    if (argc > 2)
	p = atoi(argv[2]);
    if (argc > 3)
	d = atoi(argv[3]);

    if (p < 4) {
	fprintf(stderr, "only characteristic >= 5 supported\n");
	return -1;
    }

    // set current prime field to F_p
    init_NTL_ff(p, d);

    zz_pEX     t, one;
    zz_pEratX  a4, a6;
    t = 1; t <<= 1;
    one = 1;
    switch (f) {
	case 0: a4.init(-27*t*t*((t - 1)*t + 1), one);
		a6.init(27*t*t*t*(((2*t-3)*t-3)*t+2), one);
		break;
	case 1: a4.init(-27*t*t*t*t, one);
		a6.init(54*t*t*t*t*t*(t-2), one);
		break;
	case 2: a4.init(108*t*t*t*(3-4*t), one);
		a6.init(432*t*t*t*t*t*(8*t-9), one);
		break;
	case 3: a4.init(t*t*t*(24-27*t), one);
		a6.init(t*t*t*t*((54*t-72)*t+16), one);
		break;
	case 4: a4.init(3*t*(8-9*t*t*t), one);
		a6.init((54*t*t*t - 72)*t*t*t + 16, one);
		break;
	case 5: a4.init(6*t-27, one);
		a6.init((t-18)*t+54, one);
		break;
	case 6: a4.init(-3*t*t*t*(t+8), one);
		a6.init(-2*t*t*t*t*((t-20)*t-8), one);
		break;
	case 7: a4.init(- 27*t*(t-1728)*(t-1728)*(t-1728), one);
		a6.init(54*t*(t-1728)*(t-1728)*(t-1728)*(t-1728)*(t-1728), one);
		break;
	default : assert(0);
    }

    ell_surface::init(a4, a6);

    int epsilon = 1;
    for (int i = 0; i < 2; i++) {
	ell_surfaceInfoT::affine_model  *model;
	
	model = (i == 0) ? &(ell_surfaceInfo->finite_model)
			 : &(ell_surfaceInfo->infinite_model);
	if (i == 0)
	    cerr << "j = " << model->j.num << "/"
			   << model->j.den << endl;
	cerr << "epsilon = " << model->epsilon << endl;
	cerr << endl;
	cerr << "M_sp = " << model->M_sp << endl;
	cerr << "M_ns = " << model->M_ns << endl;
	cerr << endl;
	cerr << "A = " << model->A << endl;
	cerr << endl;
	cerr << "I_star = " << model->I_star << endl;
	cerr << "II = " << model->II << endl;
	cerr << "II_star = " << model->II_star << endl;
	cerr << "III = " << model->III << endl;
	cerr << "III_star = " << model->III_star << endl;
	cerr << "IV = " << model->IV << endl;
	cerr << "IV_star = " << model->IV_star << endl;
	cerr << endl;
	epsilon *= model->epsilon;
    }
    cerr << "global epsilon = " << epsilon << endl;
    cerr << "deg(L) = " << ell_surfaceInfo->deg_L << endl;

    return 0;
}
#endif // MAIN
