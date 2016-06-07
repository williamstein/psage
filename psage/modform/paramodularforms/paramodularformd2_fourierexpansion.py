r"""
Classes describing the Fourier expansion of Paramodular modular forms.

AUTHORS:

- Martin Raum (2010 - 04 - 09) Initial version
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

from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import TrivialCharacterMonoid, \
                                           TrivialRepresentation
from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import EquivariantMonoidPowerSeriesRing
from psage.modform.jacobiforms.jacobiformd1nn_types import JacobiFormD1NN_Gamma, JacobiFormsD1NN
from operator import xor
from psage.modform.paramodularforms.paramodularformd2_fourierexpansion_cython import apply_GL_to_form, reduce_GL#, \
from sage.functions.other import ceil, floor
from sage.functions.other import sqrt
from sage.misc.functional import isqrt
from sage.misc.latex import latex
from sage.modular.modsym.p1list import P1List
from sage.rings.all import Mod
from sage.arith.all import gcd, kronecker_symbol
from sage.arith.all import legendre_symbol
from sage.rings.infinity import infinity
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.structure.sage_object import SageObject
import itertools

#===============================================================================
# ParamodularFormD2Indices_discriminant
#===============================================================================

class ParamodularFormD2Indices_discriminant( SageObject ) :
    """
    Associated with a form of level `N`.
    Indices are pairs of positive binary quadratic forms `(a, b, Nc)` and an
    element of the projective line over `\ZZ/ N \ZZ` . The last is modeled by an
    index of the element with in ``P1List(N)``.
    
    A pair `(t, v)` of a quadratic form and a left coset representative represents
    the equivalence class `v^{-1} t v^{-tr} . Every `v` represented by its first column
    of its inverse which is well defined up to units of `\ZZ / N \ZZ`.
    """
    def __init__(self, level, reduced = True) :
        if level == 1 :
            ## P1List(1) does not accept 1 as an index which is what we consider the
            ## unit element in GL(2, ZZ)
            raise NotImplementedError( "Level must not be 1")
        self.__level = level
        self.__reduced = reduced
        
        self.__p1list = P1List(level)

    def ngens(self) :
        ## FIXME: This is not correct if N != 1
        return 4 if not self.__reduced else 3
    
    def gen(self, i = 0) :
        ## FIXME: This is not correct if N != 1
        if i == 0 :
            return ((1, 0, 0), 0)
        elif i == 1 :
            return ((0, 0, 1), 0)
        elif i == 2 :
            return ((1, 1, 1), 0)
        elif not self.__reduced and i == 3 :
            return ((1, -1, 1), 0)
        
        raise ValueError, "Generator not defined"
    
    def gens(self) :    
        return [self.gen(i) for i in xrange(self.ngens())]

    def is_commutative(self) :
        return True
    
    def monoid(self) :
        return ParamodularFormD2Indices_discriminant(self.__level, False) 

    def group(self) :
        return "GL(2,ZZ)_0 (%s)" % (self.__level,)
    
    def level(self) :
        return self.__level
    
    def _p1list(self) :
        return self.__p1list
        
    def is_monoid_action(self) :
        """
        True if the representation respects the monoid structure.
        """
        return True
    
    def filter(self, disc) :
        return ParamodularFormD2Filter_discriminant(disc, self.__level, self.__reduced)
        
    def filter_all(self) :
        return ParamodularFormD2Filter_discriminant(infinity, self.__level, self.__reduced)
    
    def zero_filter(self) :
        return ParamodularFormD2Filter_discriminant(0, self.__level, self.__reduced) 
    
    def minimal_composition_filter(self, ls, rs) :
        if len(ls) == 0 or len(rs) == 0 :
            return ParamodularFormD2Filter_discriminant(0, self.__reduced)
        
        if len(ls) == 1 and ls[0][0] == (0,0,0) :
            return ParamodularFormD2Filter_discriminant(
                    min(4*self.__level*a*c - b**2 for (a,b,c) in rs),
                    self.__reduced)
        if len(rs) == 1 and rs[0][0] == (0,0,0) :
            return ParamodularFormD2Filter_discriminant(
                    min(4*self.__level*a*c - b**2 for (a,b,c) in ls),
                    self.__reduced)
        
        raise ArithmeticError, "Discriminant filter does not " + \
                               "admit minimal composition filters"
        
    def _reduction_function(self) :
        return lambda s: reduce_GL(s, self.__p1list)
    
    def reduce(self, s) :
        return reduce_GL(s, self.__p1list)
  
    def decompositions(self, s) :
        ((a, b, c), l) = s
        
        #                  r t r^tr = t_1 + t_2
        # \Rightleftarrow  t = r t_1 r^tr + r t_2 r^tr
        for a1 in xrange(a + 1) :
            a2 = a - a1
            for c1 in xrange(c + 1) :
                c2 = c - c1
                
                B1 = isqrt(4*a1*c1)
                B2 = isqrt(4*a2*c2)
                for b1 in xrange(max(-B1, b - B2), min(B1 + 1, b + B2 + 1)) :
                    h1 = apply_GL_to_form(self.__p1list[l], (a1, b1, c1))
                    h2 = apply_GL_to_form(self.__p1list[l], (a2, b - b1, c2))
                    if h1[2] % self.__level == 0 and h2[2] % self.__level == 0:
                        yield ((h1, 1), (h2,1))
        
        raise StopIteration
    
    def zero_element(self) :
        return ((0,0,0), 1)

    def __contains__(self, x) :
        return isinstance(x, tuple) and len(x) == 2 and \
               isinstance(x[0], tuple) and len(x[0]) == 3 and \
               all(isinstance(e, (int,Integer)) for e in x[0]) and \
               isinstance(x[1], (int,Integer))
        
    def __cmp__(self, other) :
        c = cmp(type(self), type(other))
        
        if c == 0 :
            c = cmp(self.__level, other.__level)
        if c == 0 :
            c = cmp(self.__reduced, other.__reduced)
            
        return c
    
    def __hash__(self) :
        return xor(hash(self.__level), hash(self.__reduced))
    
    def _repr_(self) :
        if self.__reduced :
            return "Paramodular indices for level %s" % (self.__level,)
        else :
            return "Quadratic forms over ZZ with P^1(ZZ/%s ZZ) structure" % (self.__level,)
            
    def _latex_(self) :
        if self.__reduced :
            return "Paramodular indices for level %s" % (latex(self.__level),)
        else :
            return "Quadratic forms over $\ZZ$ with $\mathbb{P}^1(\ZZ/%s \ZZ)$ structure" \
                    % (latex(self.__level),)

#===============================================================================
# ParamodularFormD2Filter_discriminant
#===============================================================================

class ParamodularFormD2Filter_discriminant ( SageObject ) :
    def __init__(self, disc, level, reduced = True) :
        self.__level = level
        
        if isinstance(disc, ParamodularFormD2Filter_discriminant) :
            disc = disc.index()

        if disc is infinity :
            self.__disc = disc
        else :
            oDmod = (-disc + 1) % (4 * level)
            Dmod = oDmod
            while Dmod > 0 :
                if not Mod(Dmod, 4 * level) :
                    Dmod = Dmod - 1
                else :
                    break
            
            self.__disc = disc - (oDmod - Dmod)
             
        self.__reduced = reduced
        
        self.__p1list = P1List(level)

    def filter_all(self) :
        return ParamodularFormD2Filter_discriminant(infinity, self.__level, self.__reduced)
    
    def zero_filter(self) :
        return ParamodularFormD2Filter_discriminant(0, self.__level, self.__reduced) 

    def is_infinite(self) :
        return self.__disc is infinity
    
    def is_all(self) :
        return self.is_infinite()

    def level(self) :
        return self.__level
        
    def index(self) :
        return self.__disc

    def _p1list(self) :
        return self.__p1list

    def _contained_trace_bound(self) :
        if self.index() == 3 :
            return 2
        else :
            return 2 * self.index() // 3

    def _contained_discriminant_bound(self) :
        return self.index()
    
    def _enveloping_discriminant_bound(self) :
        return self.index()
    
    def is_reduced(self) :
        return self.__reduced
    
    def __contains__(self, f) :
        """
        Check whether an index has discriminant less than ``self.__index`` and 
        whether its bottom right entry is divisible by ``self.__level``.
        """
        if self.__disc is infinity :
            return True
        
        (s, l) = f

        (a, b, c) = apply_GL_to_form(self.__p1list[l], s)
        if not c % self.__level == 0 :
            return False
        
        disc = 4*a*c - b**2
        if disc == 0 :
            return gcd([a,b,c]) < self._indefinite_content_bound()
        else :
            return disc < self.__disc
    
    def _indefinite_content_bound(self) :
        r"""
        Return the maximal trace for semi definite forms, which are considered
        to be below this precision.
        
        NOTE:
        
            The optimal value would be `\sqrt{D / 3}`, but this leads to very small precisions
            on converting to trace bounds.
        """
        return 2 * self.index() // 3
    
    def __iter__(self) :
        return itertools.chain( self.iter_positive_forms(),
                                self.iter_indefinite_forms() )
        
    def iter_positive_forms(self) :
        if self.__disc is infinity :
            raise ValueError, "infinity is not a true filter index"
        
        if self.__reduced :
            for (l, (u,x)) in enumerate(self.__p1list) :
                if u == 0 :
                    for a in xrange(self.__level, isqrt(self.__disc // 4) + 1, self.__level) :
                        for b in xrange(a+1) :
                            for c in xrange(a, (b**2 + (self.__disc - 1))//(4*a) + 1) :
                                yield ((a,b,c), l)
                else :
                    for a in xrange(1, isqrt(self.__disc // 3) + 1) :
                        for b in xrange(a+1) :
                            ## We need x**2 * a + x * b + c % N == 0
                            h = (-((x**2 + 1) * a + x * b)) % self.__level
                            for c in xrange( a + h,
                                             (b**2 + (self.__disc - 1))//(4*a) + 1, self.__level ) :
                                yield ((a,b,c), l)
        #! if self.__reduced
        else :
            maxtrace = floor(self.__disc / Integer(3) + sqrt(self.__disc / Integer(3)))
            for (l, (u,x)) in enumerate(self.__p1list) :
                if u == 0 :
                    for a in xrange(self.__level, maxtrace + 1, self.__level) :
                        for c in xrange(1, maxtrace - a + 1) :
                            Bu = isqrt(4*a*c - 1)
        
                            di = 4*a*c - self.__disc
                            if di >= 0 :
                                Bl = isqrt(di) + 1 
                            else :
                                Bl = 0
                            
                            for b in xrange(-Bu, -Bl + 1) :
                                yield ((a,b,c), l)
                            for b in xrange(Bl, Bu + 1) :
                                yield ((a,b,c), l)
                else :
                    for a in xrange(1, maxtrace + 1) :
                        for c in xrange(1, maxtrace - a + 1) :
                            Bu = isqrt(4*a*c - 1)
        
                            di = 4*a*c - self.__disc
                            if di >= 0 :
                                Bl = isqrt(di) + 1 
                            else :
                                Bl = 0
                            
                            h = (-x * a - int(Mod(x, self.__level)**-1) * c - Bu) % self.__level \
                                if x != 0 else \
                                (-Bu) % self.__level
                            for b in xrange(-Bu + self.__level - h, -Bl + 1, self.__level) :
                                yield ((a,b,c), l)
                            h = (-x * a - int(Mod(x, self.__level)**-1) * c + Bl) % self.__level \
                                if x !=0 else \
                                Bl % self.__level
                            for b in xrange(Bl + self.__level - h, Bu + 1, self.__level) :
                                yield ((a,b,c), l)
        #! else self.__reduced
        
        raise StopIteration
    
    def iter_indefinite_forms(self) :
        if self.__disc is infinity :
            raise ValueError, "infinity is not a true filter index"
        
        
        if self.__reduced :
            for (l, (u,_)) in enumerate(self.__p1list) :
                if u == 0 :
                    for c in xrange(self._indefinite_content_bound()) :
                        yield ((0,0,c), l)
                else :
                    for c in xrange(0, self._indefinite_content_bound(), self.__level) :
                        yield ((0,0,c), l)
        else :
            raise NotImplementedError
        
        raise StopIteration

    def _iter_positive_forms_with_content_and_discriminant(self) :
        if self.__disc is infinity :
            raise ValueError, "infinity is not a true filter index"
        
        if self.__reduced :
            for (l, (u,x)) in enumerate(self.__p1list) :
                if u == 0 :
                    for a in xrange(self.__level, isqrt(self.__disc // 4) + 1, self.__level) :
                        frpa = 4 * a
                        for b in xrange(a+1) :
                            g = gcd(a // self.__level,b)
                            bsq = b**2

                            for c in xrange(a, (b**2 + (self.__disc - 1))//(4*a) + 1) :
                                yield (((a,b,c), l), gcd(g,c), frpa*c - bsq)
                else :
                    for a in xrange(1, isqrt(self.__disc // 3) + 1) :
                        frpa = 4 * a
                        for b in xrange(a+1) :
                            g = gcd(a, b)
                            bsq = b**2
                                                        
                            ## We need x**2 * a + x * b + c % N == 0
                            h = (-((x**2 + 1) * a + x * b)) % self.__level
                            for c in xrange( a + h,
                                             (b**2 + (self.__disc - 1))//(4*a) + 1, self.__level ) :
                                yield (((a,b,c), l), gcd(g,(x**2 * a + x * b + c) // self.__level), frpa*c - bsq)
        #! if self.__reduced
        else :
            raise NotImplementedError
        
        raise StopIteration
    
    def _hecke_operator(self, n) :
        if gcd(n, self.__level) != 1 :
            raise NotImplementedError
        else :
            return ParamodularFormD2Filter_discriminant(self.__disc // n**2, self.__level, self.__reduced) 
    
    def __cmp__(self, other) :
        c = cmp(type(self), type(other))
        
        if c == 0 :
            c = cmp(self.__level, other.__level)
        if c == 0 :
            c = cmp(self.__reduced, other.__reduced)
        if c == 0 :
            c = cmp(self.__disc, other.__disc)
            
        return c

    def __hash__(self) :
        return reduce(xor, map(hash, [type(self), self.__level, self.__disc]))
                   
    def _repr_(self) :
        return "Discriminant filter (%s) for paramodular forms of level %s " \
                 % (self.__disc, self.__level)
    
    def _latex_(self) :
        return "Discriminant filter (%s) for paramodular forms of level %s " \
                 % (latex(self.__disc), latex(self.__level))


class ParamodularFormD2Filter_trace (SageObject) :
    def __init__(self, precision, level, reduced = True) :
        self.__level = level
        
        if isinstance(precision, ParamodularFormD2Filter_trace) :
            precision = precision.index()
        elif isinstance(precision, ParamodularFormD2Filter_discriminant) :
            if precision.index() is infinity :
                precision = infinity
            else :
                precision = isqrt(precision.index() - 1) + 1 

        self.__trace = precision
        self.__reduced = reduced
        self.__p1list = P1List(level)

    def filter_all(self) :
        return ParamodularFormD2Filter_trace(infinity, self.__level, self.__reduced)
    
    def zero_filter(self) :
        return ParamodularFormD2Filter_trace(0, self.__level, self.__reduced) 

    def is_infinite(self) :
        return self.__trace is infinity
    
    def is_all(self) :
        return self.is_infinite()

    def is_all(self) :
        return self.is_infinite()

    def level(self) :
        return self.__level
    
    def index(self) :
        return self.__trace

    def _contained_trace_bound(self) :
        return self.index()
        
    def _contained_discriminant_bound(self) :
        return floor( (self.index() / ( Integer(1) + self.__level + self.__level**2) )**2 )
    
    def _enveloping_discriminant_bound(self) :
        return ceil(3 * self.index() / Integer(2))
    
    def is_reduced(self) :
        return self.__reduced
    
    def __contains__(self, f) :
        r"""
        Check whether an index has discriminant less than ``self.__index`` and 
        whether its bottom right entry is divisible by ``self.__level``.
        """
        if self.__disc is infinity :
            return True
        
        (s, l) = f

        (a, _, c) = apply_GL_to_form(self.__p1list[l], s)
        if not c % self.__level == 0 :
            return False

        return a + c < self.index()
    
    def __iter__(self) :
        return itertools.chain( self.iter_positive_forms(),
                                self.iter_indefinite_forms() )
        
    def iter_positive_forms(self) :
        if self.__disc is infinity :
            raise ValueError, "infinity is not a true filter index"
        
        if self.__reduced :
            for (l, (u,x)) in enumerate(self.__p1list) :
                if u == 0 :
                    for a in xrange(self.__level, self.__trace, self.__level) :
                        for c in xrange(a, self.__trace - a) :
                            for b in xrange( 2 * isqrt(a * c - 1) + 2 if a*c != 1 else 1,
                                             2 * isqrt(a * c) ) :
                                yield ((a,b,c), l)
                else :
                    for a in xrange(1, self.__trace) :
                        for b in xrange(a+1) :
                            ## We need x**2 * a + x * b + c % N == 0
                            h = ((x**2 + 1) * a + x * b) % self.__level
                            if x == 0 and h == 0 : h = 1
                            for c in xrange( h, self.__trace - x * b - (x**2  + 1) * a, self.__level ) :
                                yield ((a,b,c), l)
        #! if self.__reduced
        else :
            for (l, (u,x)) in enumerate(self.__p1list) :
                if u == 0 :
                    for a in xrange(self.__level, self.__trace, self.__level) :
                        for c in xrange(1, self.__trace - a) :
                            for b in xrange( 2 * isqrt(a * c - 1) + 2 if a*c != 1 else 1,
                                             2 * isqrt(a * c) ) :
                                yield ((a,b,c), l)
                else :
                    for a in xrange(1, self.__trace) :
                        for c in xrange(self.__level, self.__trace - a, self.__level) :
                            for b in xrange( 2 * isqrt(a * c - 1) + 2 if a*c != 1 else 1,
                                             2 * isqrt(a * c) ) :
                                yield ((a, b - 2 * x * a, c - x * b - x**2 * a), l)
        #! else self.__reduced
        
        raise StopIteration
    
    def iter_indefinite_forms(self) :
        if self.__disc is infinity :
            raise ValueError, "infinity is not a true filter index"
        
        if self.__reduced :
            for (l, (u,_)) in enumerate(self.__p1list) :
                if u == 0 :
                    for c in xrange(self.__trace) :
                        yield ((0,0,c), l)
                else :
                    for c in xrange(0, self.__trace, self.__level) :
                        yield ((0,0,c), l)
        else :
            raise NotImplementedError
        
        raise StopIteration
    
    def _hecke_operator(self, n) :
        raise ValueError( "Non-GL(2, ZZ) invariant filter does not admit Hecke action.\n" +
                          "Use conversion to discriminant filters first.")
        
    def __cmp__(self, other) :
        c = cmp(type(self), type(other))
        
        if c == 0 :
            c = cmp(self.__level, other.__level)
        if c == 0 :
            c = cmp(self.__reduced, other.__reduced)
        if c == 0 :
            c = cmp(self.__trace, other.__trace)
            
        return c

    def __hash__(self) :
        return reduce(xor, map(hash, [type(self), self.__level, self.__trace]))
                   
    def _repr_(self) :
        return "Trace filter (%s) for paramodular forms of level %s " \
                 % (self.__disc, self.__level)
    
    def _latex_(self) :
        return "Trace filter (%s) for paramodular forms of level %s " \
                 % (latex(self.__disc), latex(self.__level))

#===============================================================================
# SiegelModularFormG2FourierExpansionRing
#===============================================================================

def ParamodularFormD2FourierExpansionRing(K, level) :
    
    R = EquivariantMonoidPowerSeriesRing(
         ParamodularFormD2Indices_discriminant(level),
         TrivialCharacterMonoid("GL(2,ZZ)_0 (%s)" % (level,), ZZ),
         TrivialRepresentation("GL(2,ZZ)_0 (%s)" % (level,), K) )

#    R._set_reduction_function(reduce_GL)

#    if K is ZZ :
#        R._set_multiply_function(mult_coeff_int)
#    else :
#        R._set_multiply_function(mult_coeff_generic)
        
    return R
