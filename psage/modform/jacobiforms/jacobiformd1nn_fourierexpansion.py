r"""
Classes describing the Fourier expansion of Jacobi forms of degree `1`
with indices in `\mathbf{N}`.

AUTHOR :

    - Martin Raum (2010 - 04 - 04) Initial version
"""

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

## even weight
## c(n, r) = c(n', r') <=> r' \equiv \pm r (2m) and r'**2 - 4 n' m = r**2 - 4 n m

from past.builtins import cmp
from builtins import map
from builtins import range
from operator import xor
from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import TrivialCharacterMonoid,\
    TrivialRepresentation
from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import EquivariantMonoidPowerSeriesModule
from psage.modform.jacobiforms.jacobiformd1nn_fourierexpansion_cython import creduce, \
                          mult_coeff_int, mult_coeff_int_weak, \
                          mult_coeff_generic, mult_coeff_generic_weak
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_method, cached_function
from sage.misc.functional import isqrt
from sage.misc.latex import latex
from sage.rings.infinity import infinity
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.structure.sage_object import SageObject
import itertools
import operator
from functools import reduce


#===============================================================================
# JacobiFormD1NNIndices
#===============================================================================

class JacobiFormD1NNIndices ( SageObject ) :
    def __init__(self, m, reduced = True, weak_forms = False) :
        r"""
        INPUT:

            - `m`                -- The index of the associated Jacobi forms.
            - ``reduced``        -- If True the reduction of Fourier indices
                                    with respect to the full Jacobi group
                                    will be considered. Otherwise, only the
                                    restriction `r**2 =< 4 n m`  or `r**2 =< 4 m n + m**2`
                                    will be considered.
            - ``weak_forms``     -- If True the weak condition
                                    `r**2 =< 4 m n` will be imposed on the
                                    indices.
        NOTE:

            The Fourier expansion of a form is assumed to be indexed
            `\sum c(n,r) z^n \zeta^r` . The indices are pairs `(n, r)`.
        """
        self.__m = m
        self.__reduced = reduced
        self.__weak_forms = weak_forms
        
    def ngens(self) :
        return len(self.gens())
    
    def gen(self, i = 0) :
        if i < self.ngens() :
            return self.gens()[i]
        
        raise ValueError("There is no generator %s" % (i,))
        
    @cached_method
    def gens(self) :
        # FIXME: This is incorrect for almost all indices m 
        return [(1,0), (1,1)]
    
    def jacobi_index(self) :
        return self.__m
    
    def is_commutative(self) :
        return True
    
    def monoid(self) :
        return JacobiFormD1NNIndices(self.__m, False, self.__weak_forms)
    
    def group(self) :
        r"""
        These are the invertible, integral lower triogonal matrices
        with bottom right entry `1`.
        """
        return "L^1_2(ZZ)"
    
    def is_monoid_action(self) :
        r"""
        True if the representation respects the monoid structure.
        """
        return False
    
    def filter(self, bound) :
        return JacobiFormD1NNFilter(bound, self.__m, self.__reduced, self.__weak_forms)
    
    def filter_all(self) :
        return JacobiFormD1NNFilter(infinity, self.__m, self.__reduced, self.__weak_forms)
    
    def minimal_composition_filter(self, ls, rs) :
        return JacobiFormD1NNFilter( min([k[0] for k in ls])
                               + min([k[0] for k in rs]),
                               self.__reduced, self.__weak_forms ) 
    
    def _reduction_function(self) :
        return lambda k: creduce(k, self.__m)
    
    def reduce(self, s) :
        return creduce(s, self.__m)
    
    def decompositions(self, s) :
        (n, r) = s
        
        fm = 4 * self.__m
        if self.__weak_forms :
            yield ((0,0), (n,r)) 
            yield ((n,r), (0,0))
            
            msq = self.__m**2
            for n1 in range(1, n) :
                n2 = n - n1
                for r1 in range( max(r - isqrt(fm * n2 + msq),
                                      isqrt(fm * n1 + msq - 1) + 1),
                                  min( r + isqrt(fm * n2 + msq) + 1,
                                       isqrt(fm * n1 + msq) + 1 ) ) :
                    yield ((n1, r1), (n2, r - r1))
        else :
            yield ((0,0), (n,r)) 
            yield ((n,r), (0,0))
            
            for n1 in range(1, n) :
                n2 = n - n1
                ##r = r1 + r2
                ##r1**2 <= 4 n1 m
                ## (r - r1)**2 <= 4 n2 m
                ## r1**2 - 2*r1*r + r**2 - 4 m n2 <= 0
                ## r1 <-> r \pm \sqrt{r**2 - r**2 + 4 m n2}
                for r1 in range( max(r - isqrt(fm * n2),
                                      isqrt(fm * n1 - 1) + 1),
                                  min( r + isqrt(fm * n2) + 1,
                                       isqrt(fm * n1) + 1 ) ) :
                    yield ((n1, r1), (n2, r - r1))
                
        raise StopIteration
            
    def zero_element(self) :
        return (0,0)
    
    def __contains__(self, k) :
        try :
            (n, r) = k
        except TypeError:
            return False
        
        return isinstance(n, (int, Integer)) and isinstance(r, (int,Integer))
    
    def __cmp__(self, other) :
        c = cmp(type(self), type(other))
        
        if c == 0 :
            c = cmp(self.__reduced, other.__reduced)
        if c == 0 :
            c = cmp(self.__weak_forms, other.__weak_forms)
        if c == 0 :
            c = cmp(self.__m, other.__m)
        
        return c

    def __hash__(self) :
        return reduce(xor, [hash(self.__m), hash(self.__reduced),
                            hash(self.__weak_forms)])
    
    def _repr_(self) :
        return "Jacobi Fourier indices for index %s forms" % (self.__m,)
    
    def _latex_(self) :
        return r"\text{Jacobi Fourier indices for index $%s$ forms}" % (latex(self.__m),)
        
#===============================================================================
# JacobiFormD1NNFilter
#===============================================================================

class JacobiFormD1NNFilter ( SageObject ) :
    r"""
    The filter which will consider the index `n` in the normal
    notation `\sum c(n,r) z^n \zeta^r`.
    """
    
    def __init__(self, bound, m, reduced = True, weak_forms = False) :
        r"""
        INPUT:

            - ``bound``          -- A natural number or exceptionally
                                    infinity reflection the bound for n.
            - `m`                -- The index of the associated Jacobi forms.
            - ``reduced``        -- If True the reduction of Fourier indices
                                    with respect to the full Jacobi group
                                    will be considered. Otherwise, only the
                                    restriction `r**2 \le 4 n m`  or `r**2 \le 4 m n + m^2`
                                    will be considered.
            - ``weak_forms``     -- If True the weak condition
                                    `r**2 \le 4 m n` will be imposed on the
                                    indices.
        NOTE:

            The Fourier expansion of a form is assumed to be indexed
            `\sum c(n,r) z^n \zeta^r` . The indices are pairs `(n, r)`.
        """
        self.__m = m
        if isinstance(bound, JacobiFormD1NNFilter) :
            bound = bound.index()
        self.__bound = bound
        self.__reduced = reduced
        self.__weak_forms = weak_forms
        
    def jacobi_index(self) :
        return self.__m
    
    def is_reduced(self) :
        return self.__reduced
    
    def is_weak_filter(self) :
        """
        Return whether this is a filter for weak Jacobi forms or not.
        """
        return self.__weak_jacobi_forms
    
    def filter_all(self) :
        return JacobiFormD1NNFilter(infinity, self.__m, self.__reduced, self.__weak_forms)

    def zero_filter(self) :
        return JacobiFormD1NNFilter(0, self.__m, self.__reduced, self.__weak_forms)
    
    def is_infinite(self) :
        return self.__bound is infinity
    
    def is_all(self) :
        return self.is_infinite()
    
    def index(self) :
        return self.__bound
    
    def __contains__(self, k) :
        m = self.__m
        
        if ( k[1]**2 > 4 * m * k[0] + m**2
             if self.__weak_forms
             else k[1]**2 > 4 * m * k[0] ) :
            return False
        
        if k[0] < self.__bound :
            return True
        elif self.__reduced :
            return creduce(k, m)[0][0] < self.__bound
        
        return False
        
    def __iter__(self) :
        return itertools.chain(self.iter_indefinite_forms(),
                               self.iter_positive_forms())
    
    def iter_positive_forms(self) :
        fm = 4 * self.__m
        if self.__reduced :
            if self.__weak_forms :
                msq = self.__m**2
                for n in range(1, self.__bound) :
                    for r in range(min(self.__m + 1, isqrt(fm * n + msq - 1) + 1)) :
                        yield (n, r)
            else :
                for n in range(1, self.__bound) :
                    for r in range(min(self.__m + 1, isqrt(fm * n - 1) + 1)) :
                        yield (n, r)
        else :
            if self.__weak_forms :
                msq = self.__m**2
                for n in range(1, self.__bound) :
                    for r in range(isqrt(fm * n + msq - 1) + 1) :
                        yield (n, r)
            else :
                for n in range(1, self.__bound) :
                    for r in range(isqrt(fm * n - 1) + 1) :
                        yield (n, r)
                        
        raise StopIteration
    
    def iter_indefinite_forms(self) :
        fm = Integer(4 * self.__m)
        
        if self.__reduced :
            if self.__weak_forms :
                msq = self.__m**2
                for n in range(0, min(self.__m // 4 + 1, self.__bound)) :
                    for r in range( isqrt(fm * n - 1) + 1 if n != 0 else 0,
                                     isqrt(fm * n + msq + 1) ) :
                        yield (n, r)
            else :

                for r in range(0, min(self.__m + 1,
                                       isqrt((self.__bound - 1) * fm) + 1) ) :
                    if fm.divides(r**2) :
                        yield (r**2 // fm, r)
        else :
            if self.__weak_forms :
                msq = self.__m**2
                for n in range(0, self.__bound) :
                    for r in range( isqrt(fm * n - 1) + 1 if n != 0 else 0,
                                     isqrt(fm * n + msq + 1) ) :
                        yield (n, r)
            else :
                for n in range(0, self.__bound) :
                    if (fm * n).is_square() :
                        yield(n, isqrt(fm * n))

        raise StopIteration
    
    def __cmp__(self, other) :
        c = cmp(type(self), type(other))
        
        if c == 0 :
            c = cmp(self.__reduced, other.__reduced)
        if c == 0 :
            c = cmp(self.__weak_forms, other.__weak_forms)
        if c == 0 :
            c = cmp(self.__m, other.__m)
        if c == 0 :
            c = cmp(self.__bound, other.__bound)
        
        return c
    
    def __hash__(self) :
        return reduce( xor, list(map(hash, [ self.__reduced, self.__weak_forms,
                                        self.__m, self.__bound ])) )

    def _repr_(self) :
        return "Jacobi precision %s" % (self.__bound,)
    
    def _latex_(self) :
        return r"\text{Jacobi precision $%s$}" % (latex(self.__bound),)

#===============================================================================
# JacobiD1NNFourierExpansionModule
#===============================================================================

def JacobiD1NNFourierExpansionModule(K, m, weak_forms = False) :
        r"""
        INPUT:

            - `m`                -- The index of the associated Jacobi forms.
            - `weak_forms`       -- If True the weak condition
                                    `r^2 \le 4 m n`n will be imposed on the
                                    indices.
        """
        
        R = EquivariantMonoidPowerSeriesModule(
             JacobiFormD1NNIndices(m, weak_forms = weak_forms),
             TrivialCharacterMonoid("L^1_2(ZZ)", ZZ),
             TrivialRepresentation("L^1_2(ZZ)", K) )
    
        if K is ZZ :
            R._set_multiply_function( lambda k, d1,d2, ch1,ch2, res : mult_coeff_int(k, d1, d2, ch1, ch2, res, m)
                                      if not weak_forms
                                      else lambda k, d1,d2, ch1,ch2, res : mult_coeff_int_weak(k, d1, d2, ch1, ch2, res, m) )
        else :
            R._set_multiply_function( lambda k, d1,d2, ch1,ch2, res : mult_coeff_generic(k, d1, d2, ch1, ch2, res, m)
                                      if not weak_forms
                                      else lambda k, d1,d2, ch1,ch2, res : mult_coeff_generic_weak(k, d1, d2, ch1, ch2, res, m) )
            
        return R
