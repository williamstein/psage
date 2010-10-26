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
#include <fstream>
#include <iostream>
#include <sstream>
using std::ostringstream;
#include <NTL/lzz_p.h>
#include <NTL/lzz_pE.h>
#include <NTL/lzz_pX.h>
#include <NTL/ZZ.h>

#include "helper.h"
#include "lzz_pEExtra.h"

NTL_CLIENT

static long gcd(long a, long b)
{
    if (b < 0)  return gcd(a, -b);
    if (a < b)  return gcd(b,  a);
    if (b == 0) return a;
    return gcd(b, a%b);

}

// suppose d divides q-1 and let chi : F_q^* --> Z/d be the
// composition of x|-->x^{(q-1)/d} and a chosen isomorphism
// F_q^*/F_q^{*d} --> Z/d;
// let z1,z2 be independent generic elements
// satisfying z1^d=z2^d=1.  we calculate the 'abstract'
// Jacobi sum
//
//   J(d,q) = sum_{x!=0,1} z1^{chi(x)}*z2^{chi(1-x)}
//
// we return the result as a (d x d) array of ZZ's
// (i.e. ZZ[d][d]).
//
// - we canonically choose the isomorphism used to construct
//   chi
// - a different choice of chi leads to J(d,q) with z1,z2
//   replaced by z1^e,z2^e, for some e in (Z/d)^*.
// - by choosing dth roots of unity z1,z2, we get an actual
//   Jacobi sum

void jacobi_sum(ZZ ***sum, long d, long *chi_of_16)
{
    long    q, e;
    q = to_long(zz_pE::cardinality());

    assert((q-1) % d == 0);
    e = (q-1)/d;

    // find an element x so x^{(q-1)/d} generates mu_d
    zz_pE   x, y, z;
    x = 1;
    while (true) {
        int r;
        r = inc(x);
        assert(r == NO_CARRY);

        power(y, x, e);
        power(z, y, d);
        assert(z == 1);

        int i;
        z = y;
        for (i = 1; z != 1; i++)
            z *= y;
        assert(d % i == 0);

        if (i < d)
            continue;

        break;
    }

    // enumerate elements of mu_d
    unsigned long  ul_mu[d];

    y = 1;
    power(z, x, e);
    for (int i = 0; i < d; i++) {
        assert((i == 0 && y == 1) || (i != 0 && y != 1));
        ul_mu[i] = to_ulong(y);
        y *= z;
    }

    // calculate log of Artin character y|-->y^{(q-1)/d_}
    unsigned long  ul_y;
    int     log[q];

    x = 1;
    do {
        power(y, x, e);
        power(z, y, d);
        assert(z == 1);

        ul_y = to_ulong(y);

        int     i;
        for (i = 0; i < d; i++)
            if (ul_mu[i] == ul_y)
                break;
        assert(i < d);

        log[to_ulong(x)] = i;
    } while (inc(x) == NO_CARRY);

    for (int i = 0; i < d; i++)
        for (int j = 0; j < d; j++)
            sum[i][j][0] = 0;

    // compute Jacobi sum as element in Z[z]/(z^d-1)
    unsigned long  u1, u2;
    long    l1, l2;
    x = 0;
    do {
        // if (true) {
            // add(y, x, 1);
            // mul(z, y, -1);
        // } else {
            sub(y, x, 1);
            mul(z, y, -1);
        // }

        if (x == 0 || z == 0)
            continue;

        u1 = to_ulong(x);
        u2 = to_ulong(z);

        l1 = log[u1];
        l2 = log[u2];
        sum[l1][l2][0]++;
    } while (inc(x) == NO_CARRY);

    // calculate chi(16)
    x = 16;
    *chi_of_16 = log[to_ulong(x)];
}

#ifdef MAIN
int main(int argc, char **argv)
{
    long    p, m, e1, e2, d, s, i;
    ofstream output;

    if (argc != 6) {
       cerr << "Usage: " << argv[0] << " p m e1 e2 d\n";
        return -1;
    }

    p  = atol(argv[1]);
    m  = atol(argv[2]);
    e1 = atol(argv[3]);
    e2 = atol(argv[4]);
    d  = atol(argv[5]);

    assert(d % p != 0);

    init_NTL_ff(p, m, 1, 1, 1);

    ZZ  sum[d];
    
    jacobi_sum(sum, e1, e2, d);

    for (int i = 0; i < d; i++) {
        if (i)
            cout << ", ";
        cout << sum[i];
    }
    cout << endl;

    return 0;
}

#endif
