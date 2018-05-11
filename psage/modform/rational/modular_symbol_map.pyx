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

AUTHOR:
    - William Stein
"""


from sage.matrix.all import matrix
from sage.rings.all import QQ, ZZ, Integer

from cysignals.memory cimport sig_free,sig_malloc
from cysignals.signals cimport sig_on,sig_off

def test_contfrac_q(a, b):
    """
    EXAMPLES::
    
        sage: import psage.modform.rational.modular_symbol_map
        sage: psage.modform.rational.modular_symbol_map.test_contfrac_q(7, 97)
        [1, 13, 14, 97]
        sage: continued_fraction(7/97).convergents()
        [0, 1/13, 1/14, 7/97]
        sage: psage.modform.rational.modular_symbol_map.test_contfrac_q(137, 93997)
        [1, 686, 6175, 43911, 93997]
        sage: continued_fraction(137/93997).convergents()
        [0, 1/686, 9/6175, 64/43911, 137/93997]
    """
    cdef long qi[MAX_CONTFRAC]
    cdef int n = contfrac_q(qi, a, b)
    return [qi[i] for i in range(n)]

cdef long contfrac_q(long qi[MAX_CONTFRAC], long a, long b) except -1:
    """
    Given coprime integers `a`, `b`, compute the sequence `q_i` of
    denominators in the continued fraction expansion of `a/b`, storing
    the results in the array `q`.
    
    Returns the number of `q_i`.  The rest of the `q` array is
    garbage.
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

cdef class ModularSymbolMap:
    def __cinit__(self):
        self.X = NULL

    def __repr__(self):
        return "Modular symbols map for modular symbols factor of dimension %s and level %s"%(self.d, self.N)
        
    def __init__(self, A):
        """
        EXAMPLES::

        We illustrate a few ways to construct the modular symbol map for 389a.
        
            sage: from psage.modform.rational.modular_symbol_map import ModularSymbolMap
            sage: A = ModularSymbols(389,sign=1).cuspidal_subspace().new_subspace().decomposition()[0]
            sage: f = ModularSymbolMap(A)
            sage: f._eval1(-3,7)
            [-2]
            sage: f.denom
            2

            sage: E = EllipticCurve('389a'); g = E.modular_symbol()
            sage: h = ModularSymbolMap(g)
            sage: h._eval1(-3,7)
            [-2]
            sage: h.denom
            1
            sage: g(-3/7)
            -2

            sage: f.denom*f.C == h.denom*h.C
            True
        """
        import sage.schemes.elliptic_curves.ell_modular_symbols as ell_modular_symbols
        from sage.modular.modsym.space import ModularSymbolsSpace

        if isinstance(A, ell_modular_symbols.ModularSymbol):
            C = self._init_using_ell_modular_symbol(A)

        elif isinstance(A, ModularSymbolsSpace):
            C = self._init_using_modsym_space(A)

        else:

            raise NotImplementedError, "Creating modular symbols from object of type %s not implemented"%type(A)

        # useful for debugging only -- otherwise is a waste of memory/space
        self.C = C
        
        X, self.denom = C._clear_denom()
        # Now store in a C data structure the entries of X (as long's)
        self.X = <long*>sig_malloc(sizeof(long*)*X.nrows()*X.ncols())
        cdef Py_ssize_t i, j, n
        n = 0
        for a in X.list():
            self.X[n] = a
            n += 1

    def _init_using_ell_modular_symbol(self, f):
        # Input f is an elliptic curve modular symbol map
        assert f.sign() != 0
        self.d = 1  # the dimension

        E = f.elliptic_curve()
        self.N = E.conductor()
        self.P1 = P1List(self.N)

        # Make a matrix whose rows are the images of the Manin symbols
        # corresponding to the elements of P^1 under f.
        n = len(self.P1)
        C = matrix(QQ, n, 1)
        for i in range(n):
            # If the ith element of P1 is (u,v) for some u,v modulo N.
            # We need to turn this into something we can evaluate f on.
            # 1. Lift to a 2x2 SL2Z matrix M=(a,b;c,d) with c=u, d=v (mod N).
            # 2. The modular symbol is M{0,oo} = {b/d,a/c}.
            # 3. So {b/d, a/c} = {b/d,oo} + {oo,a/c} = -{oo,b/d} + {oo,a/c} = f(a/c)-f(b/d).
            # 4. Thus x |--> f(a/c)-f(b/d).
            a,b,c,d = self.P1.lift_to_sl2z(i)  # output are Python ints, so careful!
            C[i,0] = (f(Integer(a)/c) if c else 0) - (f(Integer(b)/d) if d else 0)
        return C

    def _init_using_modsym_space(self, A):
        # Very slow generic setup code.  This is "O(1)" in that we
        # care mainly about evaluation time being fast, at least in
        # this code.
        if A.sign() == 0:
            raise ValueError, "A must have sign +1 or -1"
        self.d = A.dimension()
        
        self.N = A.level()
        self.P1 = P1List(self.N)

        # The key data we need from the modular symbols space is the
        # map that assigns to an element of P1 the corresponding
        # element of ZZ^n.  That's it.  We forget everything else.
        M = A.ambient_module()        
        W = matrix([M.manin_symbol(x).element() for x in self.P1])
        B = A.dual_free_module().basis_matrix().transpose()

        # Return matrix whose rows are the images of the Manin symbols
        # corresponding to the elements of P^1 under the modular
        # symbol map.
        return W*B


    def __dealloc__(self):
        if self.X:
            sig_free(self.X)

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
            j = self.P1.index((sign*q[i])%self.N, q[i-1]%self.N)
            # map over, adding a row of the matrix self.X
            # to the answer vector v.
            x = self.X + j*self.d
            for k in range(self.d):
                v[k] += x[k]
            # change sign, so q[i] is multiplied by (-1)^(i-1)
            sign *= -1

    def dimension(self):
        return self.d
        
    def _eval0(self, a, b):
        cdef long v[MAX_DEG]
        self.evaluate(v, a, b)

    def _eval1(self, a, b):
        """
        EXAMPLE::

            sage: from psage.modform.rational.modular_symbol_map import ModularSymbolMap
            sage: A = ModularSymbols(188,sign=1).cuspidal_subspace().new_subspace().decomposition()[-1]
            sage: f = ModularSymbolMap(A)
            sage: f._eval1(-3,7)
            [-3, 0]

        """
        cdef long v[MAX_DEG]
        self.evaluate(v, a, b)
        cdef int i
        return [v[i] for i in range(self.d)]
        
            
        

    
