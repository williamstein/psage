#################################################################################
#
# (c) Copyright 2010 William Stein
#
#  This file is part of PSAGE
#
#  PSAGE is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
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
Highly optimized computation of the modular symbols map associated to
a simple new space of modular symbols.
"""

include 'stdsage.pxi'

cdef enum:
    MAX_CONTFRAC = 100
    MAX_DEG = 10000

def test_contfrac_q(a, b):
    cdef long qi[MAX_CONTFRAC]
    cdef int n = contfrac_q(qi, a, b)
    return [qi[i] for i in range(n)]

cdef long contfrac_q(long qi[MAX_CONTFRAC], long a, long b) except -1:
    """
    Given coprime integers a, b, compute the sequence q_i of
    denominators in the continued fraction expansion of a/b, storing
    the results in the array q.
    
    Returns the number of q_i.  The rest of the q array is garbage.
    """
    cdef long c

    qi[0] = 1
    if a == 0:
        return 1
    if b == 0:
        raise ZeroDivisionError

    # one iteration isn't needed, since we are only computing the q_i.
    a,b = b, a%b
    cdef int i=1
    while b:
        if i >= 2:
            qi[i] = qi[i-1]*(a//b) + qi[i-2]
        else:
            qi[1] = a//b
        a,b = b, a%b
        i += 1

    return i

from sage.modular.modsym.p1list cimport P1List

cdef class ModularSymbolMap:
    cdef long d, denom, N
    cdef long* X  # coefficients of linear map from P^1 to Q^d.
    cdef P1List P1
    
    def __cinit__(self):
        self.X = NULL

    def __repr__(self):
        return "Modular symbols map for modular symbols factor of dimension %s and level %s"%(self.d, self.N)
        
    def __init__(self, A):
        # Very slow generic setup code.  This is O(1) since we care
        # about evaluation time being fast.
        from sage.matrix.all import matrix
        from sage.rings.all import ZZ
        self.N = A.level()
        self.P1 = P1List(self.N)
        if A.sign() == 0:
            raise ValueError, "A must have sign +1 or -1"
        M = A.ambient_module()        
        W = matrix([M.manin_symbol(x).element() for x in self.P1])
        B = A.dual_free_module().basis_matrix().transpose()
        X, self.denom = (W * B)._clear_denom()
        self.d = A.dimension()
        # Now store in a C data structure the entries of X (as long's), and d
        self.X = <long*>sage_malloc(sizeof(long*)*X.nrows()*X.ncols())
        cdef Py_ssize_t i, j, n
        n = 0
        for a in X.list():
            self.X[n] = a
            n += 1

    def __dealloc__(self):
        if self.X:
            sage_free(self.X)

    cdef int evaluate(self, long v[MAX_DEG], long a, long b) except -1:
        cdef long q[MAX_CONTFRAC]
        cdef int i, j, k, n, sign=1
        cdef long* x
        
        # initialize answer vector to 0
        for i in range(self.d):
            v[i] = 0
            
        # compute continued fraction
        n = contfrac_q(q, a, b)

        # compute corresponding modular symbols, mapping over...
        for i in range(1,n):
            j = self.P1.index(sign*q[i], q[i-1])
            # map over, adding a row of the matrix self.X
            # to the answer vector v.
            x = self.X + j*self.d
            for k in range(self.d):
                v[k] += x[k]
            # change sign, so q[i] is multiplied by (-1)^(i-1)
            sign *= -1

    def denominator(self):
        return self.d
        
    def _eval0(self, a, b):
        cdef long v[MAX_DEG]
        self.evaluate(v, a, b)

    def _eval1(self, a, b):
        """
        EXAMPLE::
        
            sage: A = ModularSymbols(188,sign=1).cuspidal_subspace().new_subspace().decomposition()[-1]
            sage: f = ModularSymbolMap(A)
            sage: f._eval1(-3,7)
            [-3, 0]
        """
        cdef long v[MAX_DEG]
        self.evaluate(v, a, b)
        cdef int i
        return [v[i] for i in range(self.d)]
        
            
        

    
