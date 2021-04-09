r"""
Classes describing the Fourier expansion of Siegel modular forms of genus n.

AUTHORS:

- Martin Raum (2009 - 05 - 10) Initial version
"""
from __future__ import division

#===============================================================================
# 
# Copyright (C) 2010 Martin Raum
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================

from past.builtins import cmp
from builtins import map
from builtins import range
from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids \
              import TrivialCharacterMonoid, TrivialRepresentation
from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring \
              import EquivariantMonoidPowerSeriesRing
from operator import xor
from sage.matrix.constructor import diagonal_matrix, matrix, zero_matrix, identity_matrix
from sage.matrix.matrix import is_Matrix
from sage.misc.flatten import flatten
from sage.misc.functional import isqrt
from sage.misc.latex import latex
from sage.arith.all import gcd
from sage.rings.infinity import infinity
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.structure.sage_object import SageObject
import itertools
from functools import reduce

#===============================================================================
# SiegelModularFormGnIndices_diagonal_lll
#===============================================================================

class SiegelModularFormGnIndices_diagonal_lll ( SageObject ) :
    r"""
    All positive definite quadratic forms filtered by their discriminant.
    """
    def __init__(self, n, reduced = True) :
        self.__n = n
        self.__reduced = reduced

    def ngens(self) :
        ##FIXME: This is most probably not correct
        return (self.__n * (self.__n + 1)) // 2 \
               if self.__reduced else \
               self.__n**2
    
    def gen(self, i = 0) :
        if i < self.__n :
            t = diagonal_matrix(ZZ, i * [0] + [2] + (self.__n - i - 1) * [0])
            t.set_immutable()
            
            return t
        elif i >= self.__n and i < (self.__n * (self.__n + 1)) // 2 :
            i = i - self.__n
            
            for r in range(self.__n) :
                if i >=  self.__n - r - 1 :
                    i = i - (self.__n - r - 1)
                    continue
                
                c = i + r + 1
                break
            
            t = zero_matrix(ZZ, self.__n)
            t[r,c] = 1
            t[c,r] = 1
            
            t.set_immutable()
            return t
        elif not self.__reduced and i >= (self.__n * (self.__n + 1)) // 2 \
             and i < self.__n**2 :
            i = i - (self.__n * (self.__n + 1)) // 2
            
            for r in range(self.__n) :
                if i >=  self.__n - r - 1 :
                    i = i - (self.__n - r - 1)
                    continue
                
                c = i + r + 1
                break
            
            t = zero_matrix(ZZ, self.__n)
            t[r,c] = -1
            t[c,r] = -1
            
            t.set_immutable()
            return t
            
        raise ValueError("Generator not defined")
    
    def gens(self) :    
        return [self.gen(i) for i in range(self.ngens())]

    def genus(self) :
        return self.__n
    
    def is_commutative(self) :
        return True
    
    def is_reduced(self) :
        return self.__reduced
    
    def monoid(self) :
        return SiegelModularFormGnIndices_diagonal_lll(self.__n, False) 

    def group(self) :
        return "GL(%s,ZZ)" % (self.__n,)
        
    def is_monoid_action(self) :
        """
        ``True`` if the representation respects the monoid structure.
        """
        return True
    
    def filter(self, bound) :
        return SiegelModularFormGnFilter_diagonal_lll(self.__n, bound, self.__reduced)
        
    def filter_all(self) :
        return SiegelModularFormGnFilter_diagonal_lll(self.__n, infinity, self.__reduced)
    
    def minimal_composition_filter(self, ls, rs) :
        if len(ls) == 0 or len(rs) == 0 :
            return SiegelModularFormGnFilter_diagonal_lll(self.__n, 0, self.__reduced)

        maxd = flatten( [[ml_mr[0][i,i] + ml_mr[1][i,i] for i in range(self.__n)] for ml_mr in itertools.product(ls, rs)] )

        return SiegelModularFormGnFilter_diagonal_lll(self.__n, maxd + 1, self.__reduced)
  
    def _gln_lift(self, v, position = 0) :
        """
        Create a unimodular matrix which contains v as a row or column.
        
        INPUT:
        
        - `v`           -- A primitive vector over `\ZZ`.
        - ``position``  -- (Integer, default: 0) Determines the position of `v` in the result.
                           - `0` -- left column
                           - `1` -- right column
                           - `2` -- upper row
                           - `3` -- bottom row
        """
        
        n = len(v)
        u = identity_matrix(n)
        if n == 1 :
            return u
        
        for i in range(n) :
            if v[i] < 0 :
                v[i] = -v[i]
                u[i,i] = -1
        
        while True :
            (i, e) = min(enumerate(v), key = lambda e: e[1])
            
            if e == 0 :
                cur_last_entry = n - 1
                for j in range(n - 1, -1, -1) :
                    if v[j] == 0 :
                        if cur_last_entry != j :
                            us = identity_matrix(n)
                            us[cur_last_entry,cur_last_entry] = 0
                            us[j,j] = 0
                            us[cur_last_entry,j] = 1
                            us[j,cur_last_entry] = 1
                            u = us * u
                            
                            v[j] = v[cur_last_entry]
                            v[cur_last_entry] = 0

                        cur_last_entry = cur_last_entry - 1
            
                
                us = identity_matrix(n,n) 
                us.set_block(0, 0, self._gln_lift(v[:cur_last_entry + 1]))
                
                
                u = u.inverse() * us
                u = matrix(ZZ, u)
                
                if position == 1 or position == 3 :
                    us = identity_matrix(n)
                    us[0,0] = 0
                    us[n-1,n-1] = 0
                    us[0,n-1] = 1
                    us[n-1,0] = 1
                    
                    u = u * us
                    
                if position == 0 or position == 1 :
                    return u
                elif position == 2 or position == 3 :
                    return u.transpose()
                else :
                    raise ValueError("Unknown position")
            
            if i != 0 :
                us = identity_matrix(n)
                us[0,0] = 0
                us[i,i] = 0
                us[0,i] = 1
                us[i,0] = 1
                u = us * u
                
                v[i] = v[0]
                v[0] = e

            for j in range(1,n) :
                h = - (v[j] // e)
                us = identity_matrix(n)
                us[j,0] = h
                u = us * u
                
                v[j] = v[j] + h*v[0]
  
    def _check_definiteness(self, t) :
        """
        OUTPUT:
        
            An integer. `1` if `t` is positive definite, `0` if it is semi definite and `-1` in all other cases.
            
        NOTE:
        
            We have to use this implementations since the quadratic forms' implementation has a bug. 
        """
        ## We compute the rational diagonal form and whenever there is an indefinite upper
        ## left submatrix we abort.
        t = matrix(QQ, t)
        n = t.nrows()
        
        indefinite = False
        
        for i in range(n) :
            if t[i,i] < 0 :
                return -1
            elif t[i,i] == 0 :
                indefinite = True
                for j in range(i + 1, n) :
                    if t[i,j] != 0 :
                        return -1
            else :
                for j in range(i + 1, n) :
                    t.add_multiple_of_row(j, i, -t[j,i]/t[i,i])
                    t.add_multiple_of_column(j, i, -t[i,j]/t[i,i])
                
        return 0 if indefinite else 1
        
    def _reduction_function(self) :
        return self.reduce
    
    def reduce(self, t) :
        ## We compute the rational diagonal form of t. Whenever a zero entry occures we
        ## find a primitive isotropic vector and apply a base change, such that t finally
        ## has the form diag(0,0,...,P) where P is positive definite. P will then be a
        ## LLL reduced.
        
        ot = t
        t = matrix(QQ, t)
        n = t.nrows()
        u = identity_matrix(QQ, n)
        
        for i in range(n) :
            if t[i,i] < 0 :
                return None
            elif t[i,i] == 0 :
                ## get the isotropic vector
                v = u.column(i)
                v = v.denominator() * v
                v = v / gcd(list(v))
                u = self._gln_lift(v, 0)
                
                t = u.transpose() * ot * u
                ts = self.reduce(t.submatrix(1,1,n-1,n-1))
                if ts is None :
                    return None
                t.set_block(1, 1, ts[0])
                
                t.set_immutable()
                
                return (t,1)
            else :
                for j in range(i + 1, n) :
                    us = identity_matrix(QQ, n, n)
                    us[i,j] = -t[i,j]/t[i,i]
                    u = u * us
                     
                    t.add_multiple_of_row(j, i, -t[j,i]/t[i,i])
                    t.add_multiple_of_column(j, i, -t[i,j]/t[i,i])
        
        u = ot.LLL_gram()
        t = u.transpose() * ot * u
        t.set_immutable()
        
        return (t, 1)

    def decompositions(self, t) :
        for diag in itertools.product(*[range(t[i,i] + 1) for i in range(self.__n)]) :
            for subents in  itertools.product( *[ range(-2 * isqrt(diag[i] * diag[j]), 2 * isqrt(diag[i] * diag[j]) + 1)
                                                  for i in range(self.__n - 1) 
                                                  for j in range(i + 1, self.__n) ] ) :
                t1 = matrix(ZZ, [[ 2 * diag[i]
                                   if i == j else
                                      (subents[self.__n * i - (i * (i + 1)) // 2 + j - i - 1]
                                       if i < j else
                                       subents[self.__n * j - (j * (j + 1)) // 2 + i - j - 1])
                                  for i in range(self.__n) ]
                                 for j in range(self.__n) ] )
 
                if self._check_definiteness(t1) == -1 :
                    continue
 
                t2 = t - t1
                if self._check_definiteness(t2) == -1 :
                    continue
                
                t1.set_immutable()
                t2.set_immutable()
                
                yield (t1, t2)

        raise StopIteration
    
    def zero_element(self) :
        t = zero_matrix(ZZ, self.__n)
        t.set_immutable()
        return t

    def __contains__(self, x) :
        return is_Matrix(x) and x.base_ring() is ZZ and x.nrows() == self.__n and x.ncols() == self.__n
        
    def __cmp__(self, other) :
        c = cmp(type(self), type(other))
        if c == 0 :
            c = cmp(self.__reduced, other.__reduced)
        if c == 0 :
            c = cmp(self.__n, other.__n)
            
        return c
    
    def __hash__(self) :
        return xor(hash(self.__reduced), hash(self.__n))
    
    def _repr_(self) :
        if self.__reduced :
            return "Reduced quadratic forms of rank %s over ZZ" % (self.__n,)
        else :
            return "Quadratic forms over of rank %s ZZ" % (self.__n,)
            
    def _latex_(self) :
        if self.__reduced :
            return "Reduced quadratic forms of rank %s over ZZ" % (latex(self.__n),)
        else :
            return "Quadratic forms over of rank %s ZZ" % (latex(self.__n),)

#===============================================================================
# SiegelModularFormGnFilter_diagonal_lll
#===============================================================================

class SiegelModularFormGnFilter_diagonal_lll ( SageObject ) :
    def __init__(self, n, bound, reduced = True) :
        if isinstance(bound, SiegelModularFormGnFilter_diagonal_lll) :
            bound = bound.index()

        self.__n = n
        self.__bound = bound
        self.__reduced = reduced
        
        self.__ambient = SiegelModularFormGnIndices_diagonal_lll(n, reduced)

    def filter_all(self) :
        return SiegelModularFormGnFilter_diagonal_lll(self.__n, infinity, self.__reduced)
    
    def zero_filter(self) :
        return SiegelModularFormGnFilter_diagonal_lll(self.__n, 0, self.__reduced) 

    def is_infinite(self) :
        return self.__bound is infinity
    
    def is_all(self) :
        return self.is_infinite()

    def genus(self) :
        return self.__n

    def index(self) :
        return self.__bound

    def is_reduced(self) :
        return self.__reduced
    
    def __contains__(self, t) :
        if self.__bound is infinity :
            return True

        for i in range(self.__n) :
            if t[i,i] >= 2 * self.__bound :
                return False
            
        return True
    
    def __iter__(self) :
        if self.__bound is infinity :
            raise ValueError("infinity is not a true filter index")

        if self.__reduced :
            ##TODO: This is really primitive. We barely reduce arbitrary forms.
            for diag in itertools.product(*[range(self.__bound) for _ in range(self.__n)]) :
                for subents in  itertools.product( *[ range(-2 * isqrt(diag[i] * diag[j]), 2 * isqrt(diag[i] * diag[j]) + 1)
                                                      for i in range(self.__n - 1) 
                                                      for j in range(i + 1, self.__n) ] ) :
                    t = matrix(ZZ, [[ 2 * diag[i]
                                      if i == j else
                                      (subents[self.__n * i - (i * (i + 1)) // 2 + j - i - 1]
                                       if i < j else
                                       subents[self.__n * j - (j * (j + 1)) // 2 + i - j - 1])
                                     for i in range(self.__n) ]
                                    for j in range(self.__n) ] )

                    t = self.__ambient.reduce(t)
                    if t is None : continue
                    t = t[0]
                    
                    t.set_immutable()
                    
                    yield t
        else :
            for diag in itertools.product(*[range(self.__bound) for _ in range(self.__n)]) :
                for subents in  itertools.product( *[ range(-2 * isqrt(diag[i] * diag[j]), 2 * isqrt(diag[i] * diag[j]) + 1)
                                                      for i in range(self.__n - 1) 
                                                      for j in range(i + 1, self.__n) ] ) :
                    t = matrix(ZZ, [[ 2 * diag[i]
                                      if i == j else
                                      (subents[self.__n * i - (i * (i + 1)) // 2 + j - i - 1]
                                       if i < j else
                                       subents[self.__n * j - (j * (j + 1)) // 2 + i - j - 1])
                                     for i in range(self.__n) ]
                                    for j in range(self.__n) ] )
                    if self._check_definiteness(t) == -1 :
                        continue
                    
                    t.set_immutable()
                    
                    yield t

        raise StopIteration
    
    def iter_positive_forms(self) :
        if self.__reduced :
            ##TODO: This is really primitive. We barely reduce arbitrary forms.
            for diag in itertools.product(*[range(self.__bound) for _ in range(self.__n)]) :
                for subents in  itertools.product( *[ range(-2 * isqrt(diag[i] * diag[j]), 2 * isqrt(diag[i] * diag[j]) + 1)
                                                      for i in range(self.__n - 1) 
                                                      for j in range(i + 1, self.__n) ] ) :
                    t = matrix(ZZ, [[ 2 * diag[i]
                                      if i == j else
                                      (subents[self.__n * i - (i * (i + 1)) // 2 + j - i - 1]
                                       if i < j else
                                       subents[self.__n * j - (j * (j + 1)) // 2 + i - j - 1])
                                     for i in range(self.__n) ]
                                    for j in range(self.__n) ] )
    
                    t = self.__ambient.reduce(t)
                    if t is None : continue
                    t = t[0]
                    if t[0,0] == 0 : continue
                    
                    t.set_immutable()
                    
                    yield t
        else :
            for diag in itertools.product(*[range(self.__bound) for _ in range(self.__n)]) :
                for subents in  itertools.product( *[ range(-2 * isqrt(diag[i] * diag[j]), 2 * isqrt(diag[i] * diag[j]) + 1)
                                                      for i in range(self.__n - 1) 
                                                      for j in range(i + 1, self.__n) ] ) :
                    t = matrix(ZZ, [[ 2 * diag[i]
                                      if i == j else
                                      (subents[self.__n * i - (i * (i + 1)) // 2 + j - i - 1]
                                       if i < j else
                                       subents[self.__n * j - (j * (j + 1)) // 2 + i - j - 1])
                                     for i in range(self.__n) ]
                                    for j in range(self.__n) ] )
    
                    if self.__ambient._check_definiteness(t) != 1 :
                        continue
                    
                    t.set_immutable()
                    
                    yield t
    
    def __cmp__(self, other) :
        c = cmp(type(self), type(other))
        if c == 0 :
            c = cmp(self.__reduced, other.__reduced)
        if c == 0 :
            c = cmp(self.__n, other.__n)
        if c == 0 :
            c = cmp(self.__bound, other.__bound)
            
        return c

    def __hash__(self) :
        return reduce(xor, list(map(hash, [type(self), self.__reduced, self.__bound])))
                   
    def _repr_(self) :
        return "Diagonal filter (%s)" % self.__bound
    
    def _latex_(self) :
        return "Diagonal filter (%s)" % latex(self.__bound)

def SiegelModularFormGnFourierExpansionRing(K, genus) :

    R = EquivariantMonoidPowerSeriesRing(
         SiegelModularFormGnIndices_diagonal_lll(genus),
         TrivialCharacterMonoid("GL(%s,ZZ)" % (genus,), ZZ),
         TrivialRepresentation("GL(%s,ZZ)" % (genus,), K) )

    return R
