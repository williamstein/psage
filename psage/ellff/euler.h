#ifndef EULER_H
#define EULER_H

#include <NTL/lzz_p.h>
#include <NTL/lzz_pE.h>
#include <NTL/lzz_pX.h>
#include <NTL/lzz_pEX.h>
#include <NTL/ZZX.h>

#include "ell_surface.h"
#include "lzz_pEratX.h"

NTL_CLIENT

long trace_tau(zz_pE& tau, ell_surfaceInfoT::affine_model& model);

void euler_table(long *table, long min_tau, long max_tau, int euler_deg=-1);
void euler_table(long *table, long min_tau, long max_tau, int euler_deg);

void twist_table(const zz_pEX& f, long *untwisted_table, long *twisted_table,
    long min_tau, long max_tau);

void pullback_table(const zz_pEX& finite_disc, const zz_pEX& infinite_disc,
    const zz_pEratX& f, long *base_table, long *pullback_table,
    long min_tau, long max_tau);

void sum_table(ZZ& sum, long *table, long min, long max);

void compute_c_n(ZZ **b, ZZ **c, int n);
void compute_b_n(ZZ& b_n);

#endif	// EULER_H
