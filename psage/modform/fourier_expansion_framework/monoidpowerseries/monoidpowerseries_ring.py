r"""
Rings of monoid power series and rings of equivariant monoid power series.

AUTHOR :
    -- Martin Raum (2009 - 07 - 25) Initial version
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

from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ambient import MonoidPowerSeriesAmbient_abstract, \
                                      EquivariantMonoidPowerSeriesAmbient_abstract
from sage.algebras.algebra import Algebra
from sage.rings.all import Integer
from sage.structure.element import Element
from sage.structure.parent import Parent

#===============================================================================
# MonoidPowerSeriesRing
#===============================================================================

_monoidpowerseries_ring_cache = dict()

def MonoidPowerSeriesRing(A, S) :
    r"""
    Return the globally unique monoid power series ring with indices
    over the filtered monoid `S` and coefficients in `A`.
    
    INPUT:
        - `A` -- A ring.
        - `S` -- A monoid as implemented in :class:~`fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids.NNMonoid`.

    OUTPUT:
        An instance of :class:~`.MonoidPowerSeriesRing_generic`.

    EXAMPLES::
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import MonoidPowerSeriesRing
        sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
        sage: mps is MonoidPowerSeriesRing(QQ, NNMonoid(False))
        True
    """
    global _monoidpowerseries_ring_cache
    key = (A, S)
    
    try :
        return _monoidpowerseries_ring_cache[key]
    except KeyError :
        P = MonoidPowerSeriesRing_generic(A, S)
        _monoidpowerseries_ring_cache[key] = P
        
        return P
    
#===============================================================================
# MonoidPowerSeriesRing_generic
#===============================================================================

class MonoidPowerSeriesRing_generic ( MonoidPowerSeriesAmbient_abstract, Algebra ) :
    r"""
    Given some `K` algebra `A` and a monoid `S` filtered over
    a net `\Lambda` construct a ring of monoid power series.
    
    Set `R = B[S]`. Then the projective limit of `R / R_\lambda` for
    `\lambda \in \Lambda \rightarrow \infty` considered as a
    `K` algebra is implemented by this class.
    """
    
    def __init__(self, A, S) :
        r"""
        INPUT:
            - `A` -- A ring.
            - `S` -- A monoid as implemented in :class:~`fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids.NNMonoid`.

        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import MonoidPowerSeriesRing_generic
            sage: mps = MonoidPowerSeriesRing_generic(QQ, NNMonoid(False))
            sage: mps = MonoidPowerSeriesRing_generic(ZZ, NNMonoid(False))
            sage: mps.base_ring()
            Integer Ring
            sage: (1 / 2) * mps.0
            Monoid power series in Ring of monoid power series over NN
        """
        Algebra.__init__(self, A)        
        MonoidPowerSeriesAmbient_abstract.__init__(self, A, S)

        self.__monoid_gens = \
          [ self._element_class(self, dict([(s, A.one_element())]),
                                    self.monoid().filter_all() )
            for s in S.gens() ]
    
        from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import MonoidPowerSeriesBaseRingInjection

        self._populate_coercion_lists_(
          coerce_list = [MonoidPowerSeriesBaseRingInjection(self.base_ring(), self)] + \
                        ([S] if isinstance(S, Parent) else []) )
    
    def ngens(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import MonoidPowerSeriesRing_generic
            sage: mps = MonoidPowerSeriesRing_generic(QQ, NNMonoid(False))
            sage: mps.ngens()
            1
        """
        return len(self.__monoid_gens)

    def gen(self, i = 0) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import MonoidPowerSeriesRing_generic
            sage: mps = MonoidPowerSeriesRing_generic(QQ, NNMonoid(False))
            sage: mps.gen().coefficients()
            {1: 1}
        """
        return self.gens()[i]
    
    def gens(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import MonoidPowerSeriesRing_generic
            sage: mps = MonoidPowerSeriesRing_generic(QQ, NNMonoid(False))
            sage: mps.gens()
            [Monoid power series in Ring of monoid power series over NN]
        """
        return self.__monoid_gens
    
    def construction(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import MonoidPowerSeriesRing_generic
            sage: mps = MonoidPowerSeriesRing_generic(QQ, NNMonoid(False))
            sage: (f, a) = mps.construction()
            sage: (f, a)
            (MonoidPowerSeriesRingFunctor, Rational Field)
            sage: f(a) == mps
            True
        """
        from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import MonoidPowerSeriesRingFunctor

        return MonoidPowerSeriesRingFunctor(self.monoid()), self.coefficient_domain()
    
    def _element_constructor_(self, x) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import MonoidPowerSeriesRing_generic
            sage: mps = MonoidPowerSeriesRing_generic(QQ, NNMonoid(False))
            sage: h = mps(1) # indirect doctest
            sage: h = mps(mps.monoid().zero_element())
            sage: h = mps.zero_element()
            sage: K.<rho> = CyclotomicField(6)
            sage: mps = MonoidPowerSeriesRing_generic(K, NNMonoid(False))
            sage: h = mps(rho)
            sage: h = mps(1)
        """
        if isinstance(x, int) :
            x = Integer(x)
            
        if isinstance(x, Element) :
            P = x.parent()

            if P is self.coefficient_domain() :
                return self._element_class( self, {self.monoid().zero_element(): x},
                                                self.monoid().filter_all() )
            elif self.coefficient_domain().has_coerce_map_from(P) :
                return self._element_class( self, {self.monoid().zero_element(): self.coefficient_domain()(x)},
                                                self.monoid().filter_all() )
            elif P is self.monoid() :
                return self._element_class( self, {x: self.base_ring().one_element},
                                                self.monoid().filter_all() )
                
        return MonoidPowerSeriesAmbient_abstract._element_constructor_(self, x)
    
    def _repr_(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import MonoidPowerSeriesRing_generic
            sage: MonoidPowerSeriesRing_generic(QQ, NNMonoid(False))
            Ring of monoid power series over NN
        """
        return "Ring of monoid power series over " + self.monoid()._repr_()

    def _latex_(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import MonoidPowerSeriesRing_generic
            sage: latex(MonoidPowerSeriesRing_generic(QQ, NNMonoid(False)))
            \text{Ring of monoid power series over }\Bold{N}
        """
        return r"\text{Ring of monoid power series over }" + self.monoid()._latex_()
    
###############################################################################
###############################################################################
###############################################################################

#===============================================================================
# EquivariantMonoidPowerSeriesRing
#===============================================================================

_equivariantmonoidpowerseries_ring_cache = dict()

def EquivariantMonoidPowerSeriesRing(O, C, R) :
    r"""
    Return the globally unique ring of equivariant monoid power
    over the monoid with action `O` with coefficients in the codomain `R`
    with a representation and a set of virtual characters `C`.
    
    INPUT:
        - `O` -- A monoid with an action of a group; As implemented in
                 :class:~`fourier_expansion_framework.monoidpowerseries.NNMonoid`.
        - `C` -- A monoid of characters; As implemented in ::class:~`fourier_expansion_framework.monoidpowerseries.CharacterMonoid_class`.
        - `R` -- A representation on an algebra; As implemented
                 in :class:~`fourier_expansion_framework.monoidpowerseries.TrivialRepresentation`.
        
    EXAMPLES::
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import EquivariantMonoidPowerSeriesRing
        sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
        sage: emps is EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
        True
    """

    ## TODO: Implement optional checking of the relations of the characters
    if O.group() != C.group() :
        raise ValueError("The action on S and the characters must have the same group")
    if R.base_ring() != C.codomain() :
        if C.codomain().has_coerce_map_from(R.base_ring()) :
            K = C.codomain()
            R = R.base_extend(K)
        elif R.base_ring().has_coerce_map_from(C.codomain()) :
            K = R.base_ring()
        else :
            from sage.categories.pushout import pushout
            
            try :
                K = pushout(C.codomain(), R.base_ring())
                R = R.base_extend(K)
            except :
                raise ValueError("character codomain and representation base ring have no common extension")
    
    global _equivariantmonoidpowerseries_ring_cache
    key = (O, C, R)
    
    try :
        return _equivariantmonoidpowerseries_ring_cache[key]
    except KeyError :
        P = EquivariantMonoidPowerSeriesRing_generic(O, C, R)
        _equivariantmonoidpowerseries_ring_cache[key] = P
        
        return P

#===============================================================================
# EquivariantMonoidPowerSeriesRing_generic
#===============================================================================

class EquivariantMonoidPowerSeriesRing_generic ( EquivariantMonoidPowerSeriesAmbient_abstract, Algebra ) :
    r"""
    Given some ring `A`, a monoid `S` filtered over some originated
    net `\Lambda` such that all induced submonoids are finite, a group `G`, a
    semigroup `C` with a map `c \rightarrow \mathrm{Hom}(G, Aut_K(A))`, a
    homomorphism `\phi : G -> Aut(S)` and a homomorphism `\eta : G -> C`, where
    `K` is the base ring of `A`.
    
    Suppose for every `c, c'` in `C`, and `g` in `G`, and `a, a'` in `A` we have
      `(c c') (g) (a a') = c(g)(a) c'(g)(a')`. 
    Set `R = B[C][S]`. Then the projective limit of
      `R / R_\lambda` for `\lambda \in \Lambda \rightarrow \infinity` is a
      `K`-algebra.
      
    The set of generators is the set of generators of the underlying
    monoidal power series ring and does not take into account the
    group action
    """
        
    def __init__(self, O, C, R) :
        r"""
        INPUT:
            - `O` -- A monoid with an action of a group; As implemented in
                     :class:~`fourier_expansion_framework.monoidpowerseries.NNMonoid`.
            - `C` -- A monoid of characters; As implemented in ::class:~`fourier_expansion_framework.monoidpowerseries.CharacterMonoid_class`.
            - `R` -- A representation on an algebra; As implemented
                     in :class:~`fourier_expansion_framework.monoidpowerseries.TrivialRepresentation`.

        EXAMPLES::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import EquivariantMonoidPowerSeriesRing_generic
            sage: emps = EquivariantMonoidPowerSeriesRing_generic(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ)) # indirect doctest
            sage: emps = EquivariantMonoidPowerSeriesRing_generic(NNMonoid(True), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ)) # indirect doctest
            sage: emps.base_ring()
            Integer Ring
            sage: (1 / 2) * emps.0
            Equivariant monoid power series in Ring of equivariant monoid power series over NN
        """
        
        Algebra.__init__(self, R.base_ring())
        EquivariantMonoidPowerSeriesAmbient_abstract.__init__(self, O, C, R)
    
        self.__monoid_gens = \
          [self._element_class(self,
            dict([( C.one_element(), dict([(s, self.coefficient_domain().one_element())]) )]),
            self.monoid().filter_all() )
           for s in self.action().gens()]
        self.__character_gens = \
          [self._element_class(self,
            dict([( c, dict([(self.monoid().zero_element(), self.coefficient_domain().one_element())]) )]),
            self.monoid().filter_all() )
           for c in C.gens()]
        self.__coefficient_gens = \
          [self._element_class(self,
            dict([( C.one_element(), dict([(self.monoid().zero_element(), g)]))]),
            self.monoid().filter_all() )
           for g in self.coefficient_domain().gens()]
        
        from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import MonoidPowerSeriesBaseRingInjection

        self._populate_coercion_lists_(
          coerce_list = [MonoidPowerSeriesBaseRingInjection(R.codomain(), self)] ,
          convert_list = ([O.monoid()] if isinstance(O.monoid(), Parent) else []) )

    def ngens(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import EquivariantMonoidPowerSeriesRing_generic
            sage: emps = EquivariantMonoidPowerSeriesRing_generic(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
            sage: emps.ngens()
            3
        """
        return len(self.__monoid_gens) + len(self.__character_gens) + len(self.__coefficient_gens)
    
    def gen(self, i = 0) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import EquivariantMonoidPowerSeriesRing_generic
            sage: emps = EquivariantMonoidPowerSeriesRing_generic(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
            sage: emps.gen()
            Equivariant monoid power series in Ring of equivariant monoid power series over NN
        """
        return self.gens()[i]
  
    def gens(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import EquivariantMonoidPowerSeriesRing_generic
            sage: emps = EquivariantMonoidPowerSeriesRing_generic(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
            sage: emps.gens()
            [Equivariant monoid power series in Ring of equivariant monoid power series over NN, Equivariant monoid power series in Ring of equivariant monoid power series over NN, Equivariant monoid power series in Ring of equivariant monoid power series over NN]
        """
        return self.__monoid_gens + self.__character_gens + self.__coefficient_gens

    def construction(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import EquivariantMonoidPowerSeriesRing_generic
            sage: emps = EquivariantMonoidPowerSeriesRing_generic(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
            sage: (f, a) = emps.construction()
            sage: (f, a)
            (EquivariantMonoidPowerSeriesRingFunctor, Rational Field)
            sage: f(a) == emps
            True
        """
        from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import EquivariantMonoidPowerSeriesRingFunctor

        return EquivariantMonoidPowerSeriesRingFunctor(self.action(), self.characters(), self.representation()), \
               self.coefficient_domain()

    def _element_constructor_(self, x) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import EquivariantMonoidPowerSeriesRing_generic
            sage: emps = EquivariantMonoidPowerSeriesRing_generic(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
            sage: h = emps(1)
            sage: h = emps(emps.monoid().zero_element())
            sage: h = emps.zero_element()
            sage: K.<rho> = CyclotomicField(6)
            sage: emps = EquivariantMonoidPowerSeriesRing_generic(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", K))
            sage: h = emps(rho)
            sage: h = emps(1)
        """
        if isinstance(x, int) :
            x = Integer(x)
        
        if isinstance(x, Element) :
            P = x.parent()

            if P is self.coefficient_domain() :
                return self._element_class( self,
                        {self.characters().one_element():
                                {self.monoid().zero_element(): x}},
                        self.action().filter_all() )
            elif self.coefficient_domain().has_coerce_map_from(P) :
                return self._element_class( self,
                        {self.characters().one_element():
                                {self.monoid().zero_element(): self.coefficient_domain()(x)}},
                        self.action().filter_all() )                
            elif P is self.monoid() :
                return self._element_class( self,
                        {self.characters().one_element():
                                {x: self.base_ring().one_element()}},
                        self.action().filter_all(),
                        symmetrise = True  )
            
        return EquivariantMonoidPowerSeriesAmbient_abstract._element_constructor_(self, x)

    def _cmp_(self, other) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import MonoidPowerSeriesRing
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: mps2 = MonoidPowerSeriesRing(ZZ, NNMonoid(False))
            sage: mps == MonoidPowerSeriesRing(QQ, NNMonoid(False))
            True
            sage: mps == mps2
            False
        """
        c = cmp(type(self), type(other))
        
        if c == 0 :
            c = cmp(self.base_ring(), other.base_ring())
        if c == 0 :
            c = cmp(self.monoid(), other.monoid())
            
        return c

    def _repr_(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import EquivariantMonoidPowerSeriesRing_generic
            sage: EquivariantMonoidPowerSeriesRing_generic(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
            Ring of equivariant monoid power series over NN
        """
        return "Ring of equivariant monoid power series over " + self.monoid()._repr_()

    def _latex_(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import EquivariantMonoidPowerSeriesRing_generic
            sage: latex(EquivariantMonoidPowerSeriesRing_generic(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ)))
            \text{Ring of equivariant monoid power series over }\Bold{N}
        """
        return r"\text{Ring of equivariant monoid power series over }" + self.monoid()._latex_()
