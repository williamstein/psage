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

#include <assert.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <sstream>
using std::ostringstream;
#include <stdio.h>
#include <stdlib.h>
#include <NTL/lzz_p.h>
#include <NTL/lzz_pE.h>
#include <NTL/lzz_pX.h>
#include <NTL/lzz_pEX.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_pEX.h>

#include "helper.h"
#include "ell.h"
#include "ell_surface.h"
#include "euler.h"
#include "lzz_pEExtra.h"

NTL_CLIENT

// given tau in F_q compute the 'trace of Frob_q' for the curve
// y^2=x^3+a4*x+a6 over the residue field F_q(tau).  if D!=0, then
// curve is assumed to be smooth.  otherwise, distinguishes between
// mult've and add've red'n and computes correct trace.  (D=disc)

long trace_tau(zz_pE& tau, ell_surfaceInfoT::affine_model& model)
{
    static zz_pE  a4_tau, a6_tau, b, e;

    // additive reduction
    if (IsZero(eval(model.A, tau)))
	return 0;

    // split multiplicative reduction
    if (IsZero(eval(model.M_sp, tau)))
	return +1;

    // non-split multiplicative reduction
    if (IsZero(eval(model.M_ns, tau)))
	return -1;

    // good red'n
    a4_tau = eval(model.a4, tau);
    a6_tau = eval(model.a6, tau);
    ell_pE::init(a4_tau, a6_tau);

    long   order = ell_pE::order();
    static ZZ  trace;
    trace = zz_pE::cardinality()+1-order;
    assert(trace.WideSinglePrecision());
    return to_long(trace);
}

// compute Euler factors for range of tau in P^1(F_q)
// 
// inputs:
//   E - elliptic surface
//   table - array of length >= max_tau-min_tau+1 to store traces of Frob
//   min_tau..max_tau - subrange of [0..q] correspond to points in P^1(F_q)
//
// TODO:
//   - rewrite to take advantage of (constant) field of definition for
//     (coeffs) of E so that only need to compute one trace in each
//     Galois orbit of taus
//   - allow previous speed-up to be avoided so that we can test that
//     each tau in an orbit gets the same trace

void euler_table(long *table, long min_tau, long max_tau,
        int frob_deg)
{
    static zz_pE   tau, tau_0;
    long  ul_tau, ul_tau_0, q;

    q = ell_surfaceInfo->q;
    assert(min_tau >= 0 && min_tau <= max_tau && max_tau <= q);
    assert(frob_deg == -1 || frob_deg > 0);
    assert(frob_deg == -1 || (min_tau == 0 && max_tau == q));

    // enumerate all elements of P^1(F_q) in given range
    from_ulong(min_tau, tau);
    for (ul_tau = 0; ul_tau <= max_tau-min_tau; ul_tau++) {
	// finite tau use a different model than tau=infty
	if (ul_tau + min_tau < q) {
            if (frob_deg != -1) {
                // test Frobenius map
                if (false) {
                    tau_0 = tau;
                    for (int i = 0; i < zz_pE::degree(); i++)
                        zz_pEExtraInfo->frobenius(tau_0, tau_0);
                    assert(tau_0 == tau);
                }
                if (false && zz_pEExtraInfo->frob_map != NULL) {
                    ul_tau_0 = ul_tau;
                    for (int i = 0; i < zz_pE::degree(); i++)
                        ul_tau_0 = zz_pEExtraInfo->frobenius(ul_tau_0);
                    assert(tau_0 == tau);
                }

                if (zz_pEExtraInfo->frob_map != NULL) {
                    ul_tau_0 = ul_tau;
                    do {
                        for (int i = 0; i < frob_deg; i++)
                            ul_tau_0 =
                                zz_pEExtraInfo->frobenius(ul_tau_0);
                    } while (ul_tau_0 > ul_tau);
                } else {
                    tau_0 = tau;
                    do {
                        for (int i = 0; i < frob_deg; i++)
                            zz_pEExtraInfo->frobenius(tau_0, tau_0);
                        ul_tau_0 = to_ulong(tau_0);
                    } while (ul_tau_0 > ul_tau);
                }
                if (ul_tau_0 < ul_tau)
                    table[ul_tau] = table[ul_tau_0];
                else
                    table[ul_tau]
                        = trace_tau(tau,
                                ell_surfaceInfo->finite_model);
            } else
                table[ul_tau]
                    = trace_tau(tau, ell_surfaceInfo->finite_model);
        } else {
            assert(IsZero(tau));
	    table[ul_tau]
                = trace_tau(tau, ell_surfaceInfo->infinite_model);
	}

	// increment to next tau for next iteration through loop
	if (ul_tau < max_tau-min_tau)
	    inc(tau);
    }
}

// compute Euler factors for current elliptic surface, regarded as a twist

void twist_table(const zz_pEX& f, long *untwisted_table, long *twisted_table,
    long min_tau, long max_tau)
{
    zz_pE   tau, f_of_tau;
    long  ul_tau, q;

    q = ell_surfaceInfo->q;
    assert(min_tau >= 0 && min_tau <= max_tau && max_tau <= q);

    // enumerate all elements of P^1(F_q) in given range
    from_ulong(min_tau, tau);
    for (ul_tau = 0; ul_tau <= max_tau-min_tau; ul_tau++) {
	// finite tau use a different model than tau=infty
	ell_surfaceInfoT::affine_model	*model;
	if (ul_tau + min_tau < q)
	    model = &(ell_surfaceInfo->finite_model);
	else {
	    assert(IsZero(tau));
	    model = &(ell_surfaceInfo->infinite_model);
	}

	if (eval(model->A, tau) == 0)		// additive
	    twisted_table[ul_tau] = 0;
	else if (eval(model->M_sp, tau) == 0)	// split multiplicative
	    twisted_table[ul_tau] = +1;
	else if (eval(model->M_ns, tau) == 0)	// non-split multiplicative
	    twisted_table[ul_tau] = -1;
	else {					// good
	    // determine twisting character at tau
	    // - at finite tau want chi(f(tau))
	    // - at tau=infty, want chi(g(0)) for g(u)=u^(deg(f)+e)*f(1/u)
	    //   where e in {0,1} satisfies deg(f)=e (mod 2)
	    if (ul_tau + min_tau < q)
		f_of_tau = eval(f, tau);
	    //else if (deg(f) == 0)
		//f_of_tau = 0;
	    else
		f_of_tau = LeadCoeff(f);

	    // if untwisted model has bad reduction but twisted has good
	    // reduction, need to calculate Euler factor from scratch
	    if (f_of_tau == 0)
		twisted_table[ul_tau] = trace_tau(tau, *model);
	    else {
		static int  chi, sign;

		chi  = zz_pEExtraInfo->legendre_char(f_of_tau);
		sign = (chi == 1) ? 1 : -1;

		twisted_table[ul_tau] = sign * untwisted_table[ul_tau];
	    }
	}

	// increment to next tau for next iteration through loop
	if (ul_tau < max_tau-min_tau)
	    inc(tau);
    }
}

// compute Euler factors for current elliptic surface, regarded as a pullback

void pullback_table(const zz_pEX& finite_disc, const zz_pEX& infinite_disc,
    const zz_pEratX& f, long *base_table, long *pullback_table,
    long min_tau, long max_tau)
{
    zz_pE   tau, f_of_tau;
    long  ul_tau, q;

    q = ell_surfaceInfo->q;
    assert(min_tau >= 0 && min_tau <= max_tau && max_tau <= q);

    // enumerate all elements of P^1(F_q) in given range
    from_ulong(min_tau, tau);
    for (ul_tau = 0; ul_tau <= max_tau-min_tau; ul_tau++) {
	// finite tau use a different model than tau=infty
	static ell_surfaceInfoT::affine_model	*model;
	static zz_pEX	base_disc;
	if (ul_tau + min_tau < q) {
	    model = &(ell_surfaceInfo->finite_model);
	    base_disc = finite_disc;
	} else {
	    assert(IsZero(tau));
	    model = &(ell_surfaceInfo->infinite_model);
	    base_disc = infinite_disc;
	}

	int	debug = -1;
	if (eval(model->A, tau) == 0) {		// additive
	    pullback_table[ul_tau] = 0;
	    debug = 0;
	}
	else if (eval(model->M_sp, tau) == 0) {	// split multiplicative
	    pullback_table[ul_tau] = +1;
	    debug = 1;
	}
	else if (eval(model->M_ns, tau) == 0) {	// non-split multiplicative
	    pullback_table[ul_tau] = -1;
	    debug = 2;
	}
	else {					// good
	    // if base model has bad reduction but pullback has good
	    // reduction, need to calculate Euler factor from scratch
	    if (eval(base_disc, tau) == 0)
		pullback_table[ul_tau] = trace_tau(tau, *model);
	    else {
		static unsigned long  ul_base;
		static zz_pE  den_of_tau;
		den_of_tau = eval(f.den, tau);
		if (IsZero(den_of_tau))
		    ul_base = q;
		else {
		    f_of_tau = eval(f.num, tau) / den_of_tau;
		    ul_base = to_ulong(f_of_tau);
		}
		pullback_table[ul_tau] = base_table[ul_base];
	    }
	}
#if 0
	long    test_trace = trace_tau(tau, *model);
	assert(pullback_table[ul_tau] == test_trace);
#endif

	// increment to next tau for next iteration through loop
	if (ul_tau < max_tau-min_tau)
	    inc(tau);
    }
}

void compute_c_n(ZZ **b, ZZ **c, int n)
{
    c[n][0] = 0;
    for (int m = 1; m <= n; m++)
	c[n][0] += b[m][0] * c[n-m][0];
    assert(c[n][0] % n == 0);
    c[n][0] /= n;
}

void sum_table(ZZ& sum, long *table, long min, long max)
{
    sum = 0;
    for (long i = min; i <= max; i++)
        sum += table[i];
}

void compute_b_n(ZZ& b_n)
{
    // compute traces over Frobenius
    long    q = to_ulong(zz_pE::cardinality()), i;
    
    // XXX : rewrite the following to allocate less memory and make
    //       extra invocations to euler_table

    long    *trace_table = (long*)malloc(sizeof(long)*(q+1));
    euler_table(trace_table, 0, q);

    b_n = 0;
    for (i = 0; i <= q; i++)
	b_n += trace_table[i];

    free(trace_table);
}

#ifdef MAIN
int main(int argc, char **argv)
{
    long    p, d, s, i;
    ofstream output;

    if (argc != 4) {
       cerr << "Usage: " << argv[0] << " s p d\n";
        return -1;
    }

    s = atol(argv[1]);
    p = atol(argv[2]);
    d = atol(argv[3]);

    if (s >= 8 || s < -1) {
        cerr << "Switch parameter must be in {-1,...,7}. "
	     << "You inputted " << s << ".\n";
        return -1;
    }

    if (p == 3)
        assert(s == 0 || s == 2);

#if VERBOSE
    cout << "--------------------------\n";
    switch(s) {
        case -1: cout << "\tE=custom\n"; break;
        case 0: cout << "\tE=E_Leg\n"; break;
        case 1: cout << "\tE=X_211\n"; break;
        case 2: cout << "\tE=X_321\n"; break;
        case 3: cout << "\tE=X_431\n"; break;
        case 4: cout << "\tE=E_Hes\n"; break;
        case 5: cout << "\tE=E_Vil\n"; break;
        case 6: cout << "\tE=E_Hal\n"; break;
        case 7: cout << "\tE=E_j\n"; break;
    }
#endif	// VERBOSE

#if 1
    int max_d = d;

    // setup finite fields
    zz_pEContext	ff_contexts[d];
    zz_pEExtraContext	ff_extras_contexts[d];
    for (d = 1; d <= max_d; d++) {
	init_NTL_ff(p, d, 1, 1, 1);
	ff_contexts[d-1].save();
	ff_extras_contexts[d-1].save();
    }

    ZZ  untwisted_traces[d+1], twisted_traces[d+1], pullback_traces[d+1];

    untwisted_traces[0] = 0;
    twisted_traces[0] = 0;
    pullback_traces[0] = 0;

    int untwisted_sign = 0, twisted_sign = 0, pullback_sign = 0;
    int untwisted_deg = -1, twisted_deg = -1, pullback_deg = -1;
    for (d = 1; d <= max_d; d++) {
	// set current prime field to F_p
	ff_contexts[d-1].restore();
	ff_extras_contexts[d-1].restore();

	long    q = to_long(zz_pE::cardinality());

	// setup elliptic surface
	zz_pEX     t;
	zz_pEratX  a4, a6;
	t = 1; t <<= 1;
	switch (s) {
	    case 0: a4.init(-27*t*t*((t - 1)*t + 1));
		    a6.init(27*t*t*t*(((2*t-3)*t-3)*t+2));
		    break;
	    case 1: a4.init(-27*t*t*t*t);
		    a6.init(54*t*t*t*t*t*(t-2));
		    break;
	    case 2: a4.init(108*t*t*t*(3-4*t));
		    a6.init(432*t*t*t*t*t*(8*t-9));
		    break;
	    case 3: a4.init(t*t*t*(24-27*t));
		    a6.init(t*t*t*t*((54*t-72)*t+16));
		    break;
	    case 4: a4.init(3*t*(8-9*t*t*t));
		    a6.init((54*t*t*t - 72)*t*t*t + 16);
		    break;
	    case 5: a4.init(6*t-27);
		    a6.init((t-18)*t+54);
		    break;
	    case 6: a4.init(-3*t*t*t*(t+8));
		    a6.init(-2*t*t*t*t*((t-20)*t-8));
		    break;
	    case 7: a4.init(- 27*t*(t-1728)*(t-1728)*(t-1728));
		    a6.init(54*t*(t-1728)*(t-1728)*(t-1728)*(t-1728)*(t-1728));
		    break;
	    default : assert(0);
	}
	ell_surface::init(a4, a6);

	// compute traces over Frobenius
	long    *trace_table = (long*)malloc(sizeof(long)*(q+1));
	euler_table(trace_table, 0, q);

	// compute the sum of the traces
	for (untwisted_traces[d] = 0, i = 0; i <= q; i++)
	    untwisted_traces[d] += trace_table[i];

	// calculate sign of functional equation
	int     sign_fe;
	sign_fe = ell_surfaceInfo->sign;
	if (d == 1) {
	    untwisted_sign = sign_fe;
	    untwisted_deg = ell_surfaceInfo->deg_L;
	}

#if VERBOSE
	cout << "d = " << d << endl;
	cout << "q = " << q << endl;
	cout << "untwisted sum of traces = " << untwisted_traces[d] << " : ";
	cout << "sign_fe = " << sign_fe << " : " << endl;
#endif

	// save info for untwisted curve
	ell_surfaceContext  untwisted_curve;
	untwisted_curve.save();

	// do same for twisted curve
	long    *twist_trace_table = (long*)malloc(sizeof(long)*(q+1));
	zz_pEX  f;
	// f = 1; f <<= 1; f -= 1;	// f = t-1
	// f = 1; f <<= 1; f -= 2;	// f = t-2
	f = 1; f <<= 1; // f = t

	// calculate minimal Weierstrass model, etc. for twisted surface
	zz_pEX  twist_a4 = ell_surfaceInfo->finite_model.a4*f*f;
	zz_pEX  twist_a6 = ell_surfaceInfo->finite_model.a6*f*f*f;
	ell_surface::init(twist_a4, twist_a6);

	twist_table(f, trace_table, twist_trace_table, 0, q);

	// compute the sum of the traces
	twisted_traces[d] = 0;
	for (i = 0; i <= q; i++)
	    twisted_traces[d] += twist_trace_table[i];

	// calculate sign of functional equation
	sign_fe = ell_surfaceInfo->sign;
	if (d == 1) {
	    twisted_sign = sign_fe;
	    twisted_deg = ell_surfaceInfo->deg_L;
	}

#if VERBOSE
	cout << "twisted sum of traces = " << twisted_traces[d] << " : ";
	cout << "sign_fe = " << sign_fe << " : " << endl;
#endif

	// restore info for untwisted curve
	untwisted_curve.restore();

	// do same for pullback curve
	long    *pullback_trace_table = (long*)malloc(sizeof(long)*(q+1));
	zz_pEX  num, den;
	zz_pEratX   g;
	// num = 1; den = num <<= 1; num -= 1;		// (t-1)/t
	// den = num = 1; num <<= 2;			// t^2
	// den = num = 1; num <<= 2; num += 1;		// t^2+1
	// den = num = 1; num <<= 3;			// t^3
	// den = num = 1; num <<= 3; num += 1;		// t^3+1
	// den = num = 1; num <<= 4;			// t^4
	den = num = 1; num <<= 4; num += 1;		// t^4+1
	g.init(num, den);
#if 0
	if (d == 1)
	    cout << "g = " << g.num << "/" << g.den << endl;
#endif

	// keep track of discriminant divisor for untwisted surface
	zz_pEX   finite_disc   = ell_surfaceInfo->finite_model.disc;
	zz_pEX   infinite_disc = ell_surfaceInfo->infinite_model.disc;

	// calculate minimal Weierstrass model, etc. for twisted surface
	zz_pEratX  pullback_a4 = eval(ell_surfaceInfo->finite_model.a4, g);
	zz_pEratX  pullback_a6 = eval(ell_surfaceInfo->finite_model.a6, g);
	ell_surface::init(pullback_a4, pullback_a6);

#if 0
	for (i = 0; d == 1 && i < 2; i++) {
	    ell_surfaceInfoT::affine_model  *model;
	    
	    model = (i == 0) ? &(ell_surfaceInfo->finite_model)
			     : &(ell_surfaceInfo->infinite_model);
	    //if (i == 0)
		//cerr << "j = " << model->j.num << "/"
			       //<< model->j.den << endl;
	    cerr << "M = " << model->M_sp*model->M_ns << endl;
	    cerr << "A = " << model->A << endl;
	    cerr << endl;
	    //cerr << endl;
	    //cerr << "I_star = " << model->I_star << endl;
	    //cerr << "II = " << model->II << endl;
	    //cerr << "II_star = " << model->II_star << endl;
	    //cerr << "III = " << model->III << endl;
	    //cerr << "III_star = " << model->III_star << endl;
	    //cerr << "IV = " << model->IV << endl;
	    //cerr << "IV_star = " << model->IV_star << endl;
	}
#endif

	pullback_table(finite_disc, infinite_disc, g, trace_table,
	    pullback_trace_table, 0, q);

	// compute the sum of the traces
	pullback_traces[d] = 0;
	for (i = 0; i <= q; i++)
	    pullback_traces[d] += pullback_trace_table[i];

	// calculate sign of functional equation
	sign_fe = ell_surfaceInfo->sign;
	if (d == 1) {
	    pullback_sign = sign_fe;
	    pullback_deg = ell_surfaceInfo->deg_L;
	}

	untwisted_curve.restore();

#if VERBOSE
	cout << "pullback sum of traces = " << pullback_traces[d] << " : ";
	cout << "sign_fe = " << sign_fe << " : " << endl;
#endif
    }
    assert(untwisted_sign != 0 && twisted_sign != 0);

    // convert traces to L-functions
    cout << "signs = " << untwisted_sign << ", " << twisted_sign
	<< ", " << pullback_sign;
    cout << "; degrees = " << untwisted_deg << ", " << twisted_deg
	<< ", " << pullback_deg << endl;
    for (i = 0; i < 3; i++) {
	ZZ 	*b;
	int     deg_L, sign;

	switch (i) {
	    case 0 : b = untwisted_traces;
		     deg_L = untwisted_deg;
		     sign = untwisted_sign;
		     break;
	    case 1 : b = twisted_traces;
		     deg_L = twisted_deg;
		     sign = twisted_sign;
		     break;
	    case 2 : b = pullback_traces;
		     deg_L = pullback_deg;
		     sign = pullback_sign;
		     break;
	    default : assert(0);
	}

	ZZ      c[max_d+1];
	c[0] = 1;

	int	n, N = max_d;
	for (n = 1; n <= N; n++)
	    compute_c_n(b, c, n);

	ZZX	L_fcn;
	for (n = 0; n <= N; n++)
	    SetCoeff(L_fcn, n, c[n]);

	// use rest of functional equation to determine remainder of L-fcn
	if (2*N >= deg_L) {
	    ZZ  q;
	    ff_contexts[0].restore();
	    q = zz_pE::cardinality();

	    ZZ  m;
	    if (deg_L % 2 == 1)
		m = q;
	    else
		m = 1;
	    m *= sign;
	    for (d = (deg_L+1)/2; d <= deg_L; d++, m *= q*q) {
		ZZ  actual_c = coeff(L_fcn, d);
		SetCoeff(L_fcn, d, coeff(L_fcn, deg_L-d)*m);
		assert(d > N || actual_c == coeff(L_fcn, d));
	    }
	}

	switch (i) {
	    case 0 : cout << "untwisted "; break;
	    case 1 : cout << "twisted "; break;
	    case 2 : cout << "pullback "; break;
	    default : assert(0);
	}
	cout << "L-function : " << L_fcn << endl;
    }

#else
    init_NTL_ff(p, d);

    for (q = 1, i = 0; i < d; i++)
        q *= p;

    //cout << "pi = " << pi << "\n";

    zz_pE    tau;
    long     e_infty;
    zz_pX    a4_tau, a6_tau, D_tau;

    // precompute a4, a6 of model y^2 = x^3 + a4*x + a6
    static zz_pX   t;
    t = 1; t <<= 1;

    switch (s) {
    case 0: 
        if(p != 3) {
            a4_tau = -27*t*t*((t - 1)*t + 1);
	        a6_tau = 27*t*t*t*(((2*t-3)*t-3)*t+2);
	        e_infty = 1; break;
        }
        else { // not correct model yet
            a4_tau = - 27*t*t*((t-1)*t+1)/81;
            a6_tau = 27*t*t*t*(((2*t-3)*t-3)*t+2)/729;
            e_infty = 1; break;
        }
    case 1: a4_tau = -27*t*t*t*t;
	    a6_tau = 54*t*t*t*t*t*(t-2);
	    e_infty = 1; break;
    case 2: a4_tau = 108*t*t*t*(3-4*t);
	    a6_tau = 432*t*t*t*t*t*(8*t-9);
	    e_infty = 1; break;
    case 3: a4_tau = t*t*t*(24-27*t);
	    a6_tau = t*t*t*t*((54*t-72)*t+16);
	    e_infty = 1; break;
    case 4: a4_tau = 3*t*(8-9*t*t*t);
	    a6_tau = (54*t*t*t - 72)*t*t*t + 16;
	    e_infty = 1; break;
    case 5: a4_tau = 6*t-27;
	    a6_tau = (t-18)*t+54;
	    e_infty = 0; break;
    case 6: a4_tau = -3*t*t*t*(t+8);
	    a6_tau = -2*t*t*t*t*((t-20)*t-8);
	    e_infty = -3; break;
    case 7: a4_tau = - 27*t*(t-1728)*(t-1728)*(t-1728);
        a6_tau = 54*t*(t-1728)*(t-1728)*(t-1728)*(t-1728)*(t-1728);
        e_infty = 1; break;
    case -1: cout << "a4 = "; cin >> a4_tau;
	    cout << "a6 = "; cin >> a6_tau;
	    cout << "e_infty = "; cin >> e_infty; break;
    }

    char filename[35]; // 35 max characters in filename 's-p-d.txt'
    sprintf(filename, "%ld-%ld-%ld.txt", s, p, d);

/*        for(long ctr=0; ctr<=deg(a4_tau); ctr++) {
            if(ctr==0) sprintf(filename, "[");
            strcat(filename, ltoa(coeff(a4_tau, ctr)._zz_p__rep));
            if(ctr < deg(a4_tau)) strcat(filename, "-");
            else if(ctr == deg(a4_tau)) strcat(filename, "]-");
        }
        for(long ctr=0; ctr<=deg(a6_tau); ctr++) {
            if(ctr==0) sprintf(filename, "[");
            strcat(filename, ltoa(coeff(a6_tau, ctr)._zz_p__rep));
            if(ctr < deg(a6_tau)) strcat(filename, "-");
            else if(ctr == deg(a4_tau)) strcat(filename, "]");
        }
        char fend[10];
        sprintf(fend, "-%d-%d.txt", p, d);
        strcat(filename, fend); */

    output.open(filename);

    // compute discriminant.  note, we tacitly assume we have a minimal
    // Weierstrass model over A^1/F_q.
    D_tau = 4*a4_tau*a4_tau*a4_tau + 27*a6_tau*a6_tau;

    long a_q = 0, a_tau_q;

    // compute traces for each finite tau in P^1(F_q)
    tau = 0;
    do {
	a_tau_q = trace_tau(tau, a4_tau, a6_tau, D_tau);
	output << a_tau_q << endl;
	a_q += a_tau_q;
    } while(inc(tau) == NO_CARRY);

    // trace at tau=infty is Legendre character of e_infty
    if (e_infty == 1)
	a_tau_q = 1;
    else if (e_infty != 0) {
	static zz_pE  x;
	x = e_infty;
	a_tau_q = IsOne(x^((q-1)/2)) ? 1 : -1;
    } else
	a_tau_q = 0;
    output << a_tau_q << endl;
    a_q += a_tau_q;

#if VERBOSE
    printf("a_q = %ld\n", a_q);
    if(s!=-1 && s!= 7) assert(a_q == 0);
    else cout << a_q << endl;
#endif	// VERBOSE
#endif

    return 0;
}

#endif
