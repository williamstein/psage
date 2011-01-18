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

from sage.all import is_fundamental_discriminant
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
        return [x[1] for x in f(self.params)]

    def ith_zero_means(self, i=0):
        z = self.zeros
        s = 0.0
        m = []
        for j in range(len(self.params)):
            s += z[j][i]   # i-th zero for j-th parameter (e.g., j-th discriminant)
            # The mean is s/(j+1)
            m.append( (self.params[j], s/(j+1)) )
        return m

    def ith_zeros(self, i=0):
        return [(self.params[j], self.zeros[j][i]) for j in range(len(self.params))]

class RealQuadratic(LowZeros):
    def __init__(self, num_zeros, max_D, **kwds):
        self.max_D  = max_D
        params = [D for D in range(2,max_D+1) if is_fundamental_discriminant(D)]
        super(RealQuadratic, self).__init__(num_zeros, params, **kwds)

    def __repr__(self):
        return "Family of real quadratic zeta functions with discriminant <= %s"%self.max_D

    def compute_low_zeros(self, D):
        return quadratic_twist_zeros(D, self.num_zeros)

class QuadraticImaginary(LowZeros):
    def __init__(self, num_zeros, min_D, **kwds):
        self.min_D  = min_D
        params = list(reversed([D for D in range(min_D, 0) if is_fundamental_discriminant(D)]))
        super(QuadraticImaginary, self).__init__(num_zeros, params, **kwds)

    def __repr__(self):
        return "Family of quadratic imaginary zeta functions with discriminant <= %s"%self.max_D

    def compute_low_zeros(self, D):
        return quadratic_twist_zeros(D, self.num_zeros)


class EllCurveZeros(LowZeros):
    def __init__(self, num_zeros, curves, **kwds):
        params = list(curves)
        super(EllCurveZeros, self).__init__(num_zeros, params, **kwds)

    def __repr__(self):
        return "Family of %s elliptic curve L functions"%len(self.params)

    def compute_low_zeros(self, E):
        L = Lfunction_from_elliptic_curve(E)
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
        cmd = "lcalc -z %s --twist-quadratic --start %s --finish %s"%(n, D, D)
        out = os.popen(cmd).read().split()
        return [float(out[i]) for i in range(len(out)) if i%2!=0]
    else:
        raise ValueError, "unknown algorithm '%s'"%algorithm
    
        
        
        
        
