#################################################################################
#
# (c) Copyright 2011 William Stein
#
#  This file is part of PSAGE.
#
#  PSAGE is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
# 
#  PSAGE is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################

"""
We study low zeros.
"""
from past.builtins import cmp
from builtins import range
from builtins import object
from sage.all import is_fundamental_discriminant, ZZ, parallel, var, sgn, fundamental_discriminant, kronecker_character,\
                       log,find_root
import os
from sage.libs.lcalc.lcalc_Lfunction import (
    Lfunction_from_character,
    Lfunction_from_elliptic_curve)



class LowZeros(object):
    def __init__(self, num_zeros, params, ncpus=None):
        self.num_zeros = num_zeros
        self.params = params
        self.ncpus = ncpus
        self.zeros = self.compute_all_low_zeros(self.ncpus)

    def compute_all_low_zeros(self, ncpus=None):
        if ncpus is not None:
            if ncpus == 1:
                return [self.compute_low_zeros(x) for x in self.params]
            else:
                @parallel(ncpus)
                def f(z):
                    return self.compute_low_zeros(z)
        else:
            @parallel
            def f(z):
                return self.compute_low_zeros(z)

        # Get the answers back, and sort them in the same order
        # as the input params.  This is *very* important to do, i.e.,
        # to get the order right!
        Z = dict([(x[0][0][0], x[1]) for x in f(self.params)])
        return [Z[x] for x in self.params]

    def ith_zero_means(self, i=0):
        z = self.zeros
        s = 0.0
        m = []
        for j in range(len(self.params)):
            s += z[j][i]   # i-th zero for j-th parameter (e.g., j-th discriminant)
            # The mean is s/(j+1)
            m.append((self.params[j], s/(j+1)))
        return m

    def ith_zeros(self, i=0):
        return [(self.params[j], self.zeros[j][i]) for j in range(len(self.params))]

def fundamental_discriminants(A, B):
    """Return the fundamental discriminants between A and B (inclusive), as Sage integers,
    ordered by absolute value, with negatives first when abs same."""
    v = [ZZ(D) for D in range(A, B+1) if is_fundamental_discriminant(D)]
    v.sort(lambda x,y: cmp((abs(x),sgn(x)),(abs(y),sgn(y))))
    return v

class RealQuadratic(LowZeros):
    def __init__(self, num_zeros, max_D, **kwds):
        self.max_D  = max_D
        params = fundamental_discriminants(2,max_D)
        super(RealQuadratic, self).__init__(num_zeros, params, **kwds)

    def __repr__(self):
        return "Family of real quadratic zeta functions with discriminant <= %s"%self.max_D

    def compute_low_zeros(self, D):
        return quadratic_twist_zeros(D, self.num_zeros)

class QuadraticImaginary(LowZeros):
    def __init__(self, num_zeros, min_D, **kwds):
        self.min_D  = min_D
        params = fundamental_discriminant(min_D, -1)
        super(QuadraticImaginary, self).__init__(num_zeros, params, **kwds)

    def __repr__(self):
        return "Family of quadratic imaginary zeta functions with discriminant <= %s"%self.max_D

    def compute_low_zeros(self, D):
        return quadratic_twist_zeros(D, self.num_zeros)


class EllCurveZeros(LowZeros):
    def __init__(self, num_zeros, curves, **kwds):
        # sort the curves into batches by conductor
        d = {}
        for E in curves:
            N = E.conductor()
            if N in d:
                d[N].append(E)
            else:
                d[N] = [E]
        params = list(sorted(d.keys()))
        self.d = d
        super(EllCurveZeros, self).__init__(num_zeros, params, **kwds)

    def __repr__(self):
        return "Family of %s elliptic curve L functions"%len(self.params)

    def compute_low_zeros(self, N):
        a = []
        for E in self.d[N]:
            L = Lfunction_from_elliptic_curve(E)
            a.append(L.find_zeros_via_N(self.num_zeros))
        num_curves = len(a)
        return [sum(a[i][j] for i in range(num_curves))/num_curves
                for j in range(self.num_zeros)]

class EllQuadraticTwists(LowZeros):
    def __init__(self, num_zeros, curve, discs, number_of_coeffs=10000, **kwds):
        self.curve = curve
        self.number_of_coeffs = number_of_coeffs
        params = discs
        super(EllQuadraticTwists, self).__init__(num_zeros, params, **kwds)

    def __repr__(self):
        return "Family of %s elliptic curve L functions"%len(self.params)

    def compute_low_zeros(self, D):
        L = Lfunction_from_elliptic_curve(self.curve.quadratic_twist(D),
                                          number_of_coeffs=self.number_of_coeffs)
        return L.find_zeros_via_N(self.num_zeros)


def quadratic_twist_zeros(D, n, algorithm='clib'):
    """
    Return imaginary parts of the first n zeros of all the Dirichlet
    character corresponding to the quadratic field of discriminant D.

    INPUT:
        - D -- fundamental discriminant
        - n -- number of zeros to find
        - algorithm -- 'clib' (use C library) or 'subprocess' (open another process)
    """
    if algorithm == 'clib':
        L = Lfunction_from_character(kronecker_character(D), type="int")
        return L.find_zeros_via_N(n)
    elif algorithm == 'subprocess':
        assert is_fundamental_discriminant(D)
        cmd = "lcalc -z {0} --twist-quadratic --start {1} --finish {2}".format(n, D, D)
        out = os.popen(cmd).read().split()
        return [float(out[i]) for i in range(len(out)) if i%2!=0]
    else:
        raise ValueError("unknown algorithm '{0}'".format(algorithm))
    
        
        
        
        
def fit_to_power_of_log(v):
    """
    INPUT:
        - v -- a list of (x,y) values, with x increasing.
    OUTPUT:
        - number a such that data is "approximated" by b*log(x)^a.
    """
    # ALGORITHM: transform data to (log(log(x)), log(y)) and find the slope
    # of the best line that fits this transformed data. 
    # This is the right thing to do, since if y = log(x)^a, then log(y) = a*log(log(x)).
    from math import log, exp
    w = [(log(log(x)), log(y)) for x,y in v]
    a, b = least_squares_fit(w)
    return float(a), exp(float(b))


def least_squares_fit(v):
    """
    INPUT:
        - v -- a list of (x,y) values that are floats
    OUTPUT:
        - a and b such that the line y=a*x + b is the least squares fit for the data.

    All computations are done using floats. 
    """
    import numpy
    x = numpy.array([a[0] for a in v])
    y = numpy.array([a[1] for a in v])
    A = numpy.vstack([x, numpy.ones(len(x))]).T
    a, b = numpy.linalg.lstsq(A,y)[0]    
    return a, b
    
    
class XLogInv(object):
    def __init__(self, a):
        self.a = a
        x = var('x')
        self.f = (x * (log(x) - a))._fast_float_(x)

    def __repr__(self):
        return "The inverse of the function x*(log(x) - %s), for sufficiently large positive x"%self.a

    def __call__(self, y):
        return find_root(self.f - float(y), 1, 10e9)
