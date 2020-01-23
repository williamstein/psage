r"""
Classes describing the fourier expansion of Siegel modular forms genus 2.

AUTHORS:

- Martin Raum (2009 - 07 - 28) Initial version
"""

#===============================================================================
# 
# Copyright (C) 2009 Martin Raum
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

from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import TrivialCharacterMonoid, \
                                           TrivialRepresentation
from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import EquivariantMonoidPowerSeriesRing
from operator import xor
from sage.functions.other import floor
from sage.functions.other import sqrt
from sage.misc.functional import isqrt
from sage.misc.latex import latex
from sage.arith.all import gcd
from sage.rings.infinity import infinity
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.structure.sage_object import SageObject
from psage.modform.paramodularforms.siegelmodularformg2_fourierexpansion_cython import \
                     mult_coeff_int_without_character, \
                     mult_coeff_generic_without_character, \
                     reduce_GL, xreduce_GL
from sage.rings.integer import Integer
from functools import reduce

#===============================================================================
# SiegelModularFormG2Filter_discriminant
#===============================================================================

class SiegelModularFormG2Filter_discriminant ( SageObject ) :
    def __init__(self, disc, reduced = True) :
        if isinstance(disc, SiegelModularFormG2Filter_discriminant) :
            disc = disc.index()

        if disc is infinity :
            self.__disc = disc
        else :
            Dmod = disc % 4
            if Dmod >= 2 :
                self.__disc = disc - Dmod + 1
            else :
                self.__disc = disc

        self.__reduced = reduced

    def filter_all(self) :
        return SiegelModularFormG2Filter_discriminant(infinity, self.__reduced)
    
    def zero_filter(self) :
        return SiegelModularFormG2Filter_discriminant(0, self.__reduced) 

    def is_infinite(self) :
        return self.__disc is infinity

    def is_all(self) :
        return self.is_infinite()

    def index(self) :
        return self.__disc

    def discriminant(self) :
        return self.index()
    
    def is_reduced(self) :
        return self.__reduced
    
    def __contains__(self, f) :
        if self.__disc is infinity :
            return True
        
        (a, b, c) = f
        disc = 4*a*c - b**2
        if disc == 0 :
            return gcd([a, b, c]) < self._indefinite_content_bound()
        else :
            return disc < self.__disc
    
    def _indefinite_content_bound(self) :
        r"""
        Return the maximal trace for semi definite forms, which are considered
        to be below this precision.
        """
        return max(1, 2 * self.index() // 3) if self.index() != 0 else 0
    
    def __iter__(self) :
        if self.__disc is infinity :
            raise ValueError("infinity is not a true filter index")
        
        if self.__reduced :
            for c in xrange(0, self._indefinite_content_bound()) :
                yield (0,0,c)
                
            for a in xrange(1,isqrt(self.__disc // 3) + 1) :
                for b in xrange(a+1) :
                    for c in xrange(a, (b**2 + (self.__disc - 1))//(4*a) + 1) :
                        yield (a,b,c)
        else :
            ##FIXME: These are not all matrices
            for a in xrange(0, self._indefinite_content_bound()) :
                yield (a,0,0)
            for c in xrange(1, self._indefinite_content_bound()) :
                yield (0,0,c)
            
            maxtrace = floor(5*self.__disc / 15 + sqrt(self.__disc)/2)
            for a in xrange(1, maxtrace + 1) :
                for c in xrange(1, maxtrace - a + 1) :
                    Bu = isqrt(4*a*c - 1)

                    di = 4*a*c - self.__disc
                    if di >= 0 :
                        Bl = isqrt(di) + 1 
                    else :
                        Bl = 0
                    
                    for b in xrange(-Bu, -Bl + 1) :
                        yield (a,b,c)
                    for b in xrange(Bl, Bu + 1) :
                        yield (a,b,c)
        #! if self.__reduced
        
        raise StopIteration

    def iter_positive_forms_with_content(self) :
        if self.__disc is infinity :
            raise ValueError("infinity is not a true filter index")
        
        
        if self.__reduced :        
            for a in xrange(1,isqrt(self.__disc // 3) + 1) :
                for b in xrange(a+1) :
                    g = gcd(a, b)
                    for c in xrange(a, (b**2 + (self.__disc - 1))//(4*a) + 1) :
                        yield (a,b,c), gcd(g,c)
        else :
            maxtrace = floor(5*self.__disc / 15 + sqrt(self.__disc)/2)
            for a in xrange(1, maxtrace + 1) :
                for c in xrange(1, maxtrace - a + 1) :
                    g = gcd(a,c)
                    
                    Bu = isqrt(4*a*c - 1)

                    di = 4*a*c - self.__disc
                    if di >= 0 :
                        Bl = isqrt(di) + 1 
                    else :
                        Bl = 0
                    
                    for b in xrange(-Bu, -Bl + 1) :
                        yield (a,b,c), gcd(g,b)
                    for b in xrange(Bl, Bu + 1) :
                        yield (a,b,c), gcd(g,b)
        #! if self.__reduced

        raise StopIteration
    
    def iter_indefinite_forms(self) :
        if self.__disc is infinity :
            raise ValueError("infinity is not a true filter index")
        
        
        if self.__reduced :
            for c in xrange(self._indefinite_content_bound()) :
                yield (0, 0, c)
        else :
            raise NotImplementedError
        
        raise StopIteration
    
    def _hecke_operator(self, n) :
        return SiegelModularFormG2Filter_discriminant(self.__disc // n**2, self.__reduced) 
    
    def __cmp__(self, other) :
        c = cmp(type(self), type(other))
        if c == 0 :
            c = cmp(self.__reduced, other.__reduced)
        if c == 0 :
            c = cmp(self.__disc, other.__disc)
            
        return c

    def __hash__(self) :
        return reduce(xor, map(hash, [type(self), self.__reduced, self.__disc]))
                                         
    def _repr_(self) :
        return "Discriminant filter (%s)" % self.__disc
    
    def _latex_(self) :
        return "Discriminant filter (%s)" % latex(self.__disc)
    
#===============================================================================
# SiegelModularFormG2Indices_discriminant_xreduce
#===============================================================================

class SiegelModularFormG2Indices_discriminant_xreduce( SageObject ) :
    r"""
    All positive definite quadratic forms filtered by their discriminant.
    """
    def __init__(self, reduced = True) :
        self.__reduced = reduced

    def ngens(self) :
        return 4 if not self.__reduced else 3
    
    def gen(self, i = 0) :
        if i == 0 :
            return (1, 0, 0)
        elif i == 1 :
            return (0, 0, 1)
        elif i == 2 :
            return (1, 1, 1)
        elif not self.__reduced and i == 3 :
            return (1, -1, 1)
        
        raise ValueError("Generator not defined")
    
    def gens(self) :    
        return [self.gen(i) for i in xrange(self.ngens())]

    def is_commutative(self) :
        return True
    
    def monoid(self) :
        return SiegelModularFormG2Indices_discriminant_xreduce(False) 

    def group(self) :
        return "GL(2,ZZ)"
        
    def is_monoid_action(self) :
        r"""
        True if the representation respects the monoid structure.
        """
        return True
    
    def filter(self, disc) :
        return SiegelModularFormG2Filter_discriminant(disc, self.__reduced)
        
    def filter_all(self) :
        return SiegelModularFormG2Filter_discriminant(infinity, self.__reduced)
    
    def minimal_composition_filter(self, ls, rs) :
        if len(ls) == 0 or len(rs) == 0 :
            return SiegelModularFormG2Filter_discriminant(0, self.__reduced)
        
        if len(ls) == 1 and ls[0] == (0,0,0) :
            return SiegelModularFormG2Filter_discriminant(min(4*a*c - b**2 for (a,b,c) in rs) + 1,
                                      self.__reduced)
        if len(rs) == 1 and rs[0] == (0,0,0) :
            return SiegelModularFormG2Filter_discriminant(min(4*a*c - b**2 for (a,b,c) in ls) + 1,
                                      self.__reduced)
        
        raise ArithmeticError("Discriminant filter does not " + \
                               "admit minimal composition filters")
        
    def _reduction_function(self) :
        return xreduce_GL
    
    def reduce(self, s) :
        return xreduce_GL(s)
  
    def decompositions(self, s) :
        (a, b, c) = s
        
        for a1 in xrange(a + 1) :
            a2 = a - a1
            for c1 in xrange(c + 1) :
                c2 = c - c1
                
                B1 = isqrt(4*a1*c1)
                B2 = isqrt(4*a2*c2)
                for b1 in xrange(max(-B1, b - B2), min(B1 + 1, b + B2 + 1)) :
                    yield ((a1, b1, c1), (a2,b - b1, c2))
        
        raise StopIteration
    
    def zero_element(self) :
        return (0,0,0)

    def __contains__(self, x) :
        return isinstance(x, tuple) and len(x) == 3 and \
               all(isinstance(e, (int,Integer)) for e in x)
        
    def __cmp__(self, other) :
        c = cmp(type(self), type(other))
        if c == 0 :
            c = cmp(self.__reduced, other.__reduced)
            
        return c
    
    def __hash__(self) :
        return hash(self.__reduced)
    
    def _repr_(self) :
        if self.__reduced :
            return "Reduced quadratic forms over ZZ"
        else :
            return "Quadratic forms over ZZ"
            
    def _latex_(self) :
        return self._repr_() 

#===============================================================================
# SiegelModularFormG2VVRepresentation
#===============================================================================

class SiegelModularFormG2VVRepresentation ( SageObject ) :
    def __init__(self, K) :
        self.__K = K
        self.__R = PolynomialRing(K, ['x', 'y'])
        self.__x = self.__R.gen(0)
        self.__y = self.__R.gen(1)
    
    def from_module(self, R) :
        assert R == PolynomialRing(R.base_ring(), ['x', 'y'])
        
        return SiegelModularFormG2VVRepresentation(R.base_ring())
        
    def base_ring(self) :
        return self.__K
    
    def codomain(self) :
        return self.__R
    
    def base_extend(self, L) :
        if L.has_coerce_map_from(self.__K) :
            return SiegelModularFormG2VVRepresentation( L )
        
        raise ValueError("Base extension of representation is not defined")
        
    def extends(self, other) :
        if isinstance(other, TrivialRepresentation) :
            return self.__K.has_coerce_map_from(other.codomain())
        elif type(self) != type(other) :
            return False
        
        return self.__K.has_coerce_map_from(other.__K)
        
    def group(self) :
        return "GL(2,ZZ)"
    
    def _apply_function(self) :
        return self.apply
    
    def apply(self, g, a) :
        return a(g[0]*self.__x + g[1]*self.__y, g[2]*self.__x + g[3]*self.__y)
    
    def __cmp__(self, other) :
        c = cmp(type(self), type(other))
        
        if c == 0 :
            c = cmp(self.__K, other.__K)
            
        return c
    
    def __hash__(self) :
        return hash(self.__K)

    def _repr_(self) :
        return "Siegel modular form degree 2 vector valued representation on %s" % (self.__K)
    
    def _latex_(self) :
        return "Siegel modular form degree 2 vector valued representation on %s" % (latex(self.__K))

#===============================================================================
# SiegelModularFormG2FourierExpansionRing
#===============================================================================

def SiegelModularFormG2FourierExpansionRing(K, with_character = False) :
    
    if with_character :
        raise NotImplementedError
    
        #R = EquivariantMonoidPowerSeriesRing(
        #     SiegelModularFormG2Indices_discriminant_xreduce(),
        #     GL2CharacterMonoid(K),
        #     TrivialRepresentation("GL(2,ZZ)", K) )
        # Characters in GL2CharacterMonoid should accept (det, sgn)
        #R._set_reduction_function(sreduce_GL)

        #if K is ZZ :
        #    R._set_multiply_function(mult_coeff_int_character)
        #else :
        #    R._set_multiply_function(mult_coeff_generic_character)
    else : 
        R = EquivariantMonoidPowerSeriesRing(
             SiegelModularFormG2Indices_discriminant_xreduce(),
             TrivialCharacterMonoid("GL(2,ZZ)", ZZ),
             TrivialRepresentation("GL(2,ZZ)", K) )

        R._set_reduction_function(reduce_GL)
    
        if K is ZZ :
            R._set_multiply_function(mult_coeff_int_without_character)
        else :
            R._set_multiply_function(mult_coeff_generic_without_character)
        
    return R

#===============================================================================
# SiegelModularFormG2VVFourierExpansionRing
#===============================================================================

## TODO: This is far from optimal. For specific weights we can use
##       free modules.
def SiegelModularFormG2VVFourierExpansionRing(K) :
    R = EquivariantMonoidPowerSeriesRing(
         SiegelModularFormG2Indices_discriminant_xreduce(),
         TrivialCharacterMonoid("GL(2,ZZ)", ZZ),
         SiegelModularFormG2VVRepresentation(K) )
        
    return R
