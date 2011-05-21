"""
Modules of monoid power series and modules of equivariant monoid power series.

AUTHOR :
    -- Martin Raum (2010 - 02 - 10) Initial version
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

from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ambient import MonoidPowerSeriesAmbient_abstract,\
                                      EquivariantMonoidPowerSeriesAmbient_abstract
from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import TrivialRepresentation
from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import EquivariantMonoidPowerSeriesFunctor,\
    MonoidPowerSeriesFunctor
from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import MonoidPowerSeriesRing, \
                                   EquivariantMonoidPowerSeriesRing
from sage.modules.module import Module
from sage.rings.integer import Integer
from sage.structure.element import Element

_monoidpowerseries_module_cache = dict()
_equivariantmonoidpowerseries_module_cache = dict()

#===============================================================================
# MonoidPowerSeriesModule
#===============================================================================

def MonoidPowerSeriesModule(A, S) :
    """
    Return the globally unique monoid power series ring with indices
    in the filtered monoid `S` and coefficients in `A`.
    
    INPUT:
        - `A` -- A module.
        - `S` -- A monoid as implemented in :class:~`fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids.NNMonoid`.

    EXAMPLES::
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import MonoidPowerSeriesModule
        sage: mps = MonoidPowerSeriesModule(QQ, NNMonoid(False))
        sage: mps is MonoidPowerSeriesModule(QQ, NNMonoid(False))
        True
    """
    global _monoidpowerseries_module_cache
    key = (A, S)
    
    try :
        return _monoidpowerseries_module_cache[key]
    except KeyError :
        P = MonoidPowerSeriesModule_generic(A, S)
        _monoidpowerseries_module_cache[key] = P
        
        return P

#===============================================================================
# MonoidPowerSeriesModule_generic
#===============================================================================

class MonoidPowerSeriesModule_generic ( MonoidPowerSeriesAmbient_abstract, Module ) :
    """
    Given some `K` module `A` and a monoid `S` filtered over
    a net `\Lambda` construct a module of monoid power series.
    
    Set `R = B[S]`. Then the projective limit of `R / R_\lambda` for
    `\lambda \in \Lambda \rightarrow \infty` considered as a
    `K` module is implemented by this class.

    NOTE:
        The implementation respects left and right modules.
    """

    def __init__(self, A, S) :
        """
        INPUT:
            - `A` -- A module.
            - `S` -- A monoid as implemented in :class:~`fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids.NNMonoid`.

        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import MonoidPowerSeriesModule_generic
            sage: mps = MonoidPowerSeriesModule_generic(FreeModule(QQ,2), NNMonoid(False))
        """       
        Module.__init__(self, MonoidPowerSeriesRing(A.base_ring(), S))
        MonoidPowerSeriesAmbient_abstract.__init__(self, A, S)
        
        self.__coeff_gens = \
          [ self._element_class(self, dict([(S.zero_element(), a)]),
                                    self.monoid().filter_all() )
            for a in A.gens() ]

    def ngens(self) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import MonoidPowerSeriesModule_generic
            sage: mps = MonoidPowerSeriesModule_generic(FreeModule(QQ,2), NNMonoid(False))
            sage: mps.ngens()
            2
        """
        return len(self.__coeff_gens)
    
    def gen(self, i = 0) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import MonoidPowerSeriesModule_generic
            sage: mps = MonoidPowerSeriesModule_generic(FreeModule(QQ,2), NNMonoid(False))
            sage: mps.gen()
            Monoid power series in Module of monoid power series over NN
        """
        return self.gens()[i]
        
    def gens(self) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import MonoidPowerSeriesModule_generic
            sage: mps = MonoidPowerSeriesModule_generic(FreeModule(QQ,2), NNMonoid(False))
            sage: mps.gens()
            [Monoid power series in Module of monoid power series over NN, Monoid power series in Module of monoid power series over NN]
        """
        return self.__coeff_gens
    
    def construction(self) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import MonoidPowerSeriesModule_generic
            sage: mps = MonoidPowerSeriesModule_generic(FreeModule(QQ,2), NNMonoid(False))
            sage: (f, a) = mps.construction()
            sage: (f, a)
            (MonoidPowerSeriesFunctor, Vector space of dimension 2 over Rational Field)
            sage: f(a) == mps
            True
        """
        return MonoidPowerSeriesFunctor(self.monoid()), self.coefficient_domain()

    def zero_element(self) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import MonoidPowerSeriesModule_generic
            sage: mps = MonoidPowerSeriesModule_generic(FreeModule(QQ,2), NNMonoid(False))
            sage: h = mps.zero_element()
        """
        return self(0)
    
    def _element_constructor_(self, x) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import MonoidPowerSeriesModule_generic
            sage: mps = MonoidPowerSeriesModule_generic(FreeModule(QQ,2), NNMonoid(False))
            sage: h = mps(0) # indirect doctest
            sage: h = mps(int(0)) # indirect doctest
        """
        if isinstance(x, int) and x == 0 :
            return self._element_class( self, dict(),
                        self.monoid().filter_all() )
        if isinstance(x, Element) and x.is_zero() :
            P = x.parent()

            if self.base_ring().base_ring().has_coerce_map_from(P) :
                return self._element_class( self, dict(),
                    self.monoid().filter_all() )
    
        return MonoidPowerSeriesAmbient_abstract._element_constructor_(self, x)
 
    def _repr_(self) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import MonoidPowerSeriesModule_generic
            sage: MonoidPowerSeriesModule_generic(FreeModule(QQ,2), NNMonoid(False))
            Module of monoid power series over NN
        """
        return "Module of monoid power series over " + self.monoid()._repr_()

    def _latex_(self) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import MonoidPowerSeriesModule_generic
            sage: latex(MonoidPowerSeriesModule_generic(FreeModule(QQ,2), NNMonoid(False)))
            Module of monoid power series over $\NN$
        """
        return "Module of monoid power series over " + self.monoid()._latex_()
    
###############################################################################
###############################################################################
###############################################################################

#===============================================================================
# EquivariantMonoidPowerSeriesModule
#===============================================================================

def EquivariantMonoidPowerSeriesModule(O, C, R) :
    """
    Return the globally unique module of equivariant monoid power
    over the monoid with action `O` with coefficients in the codomain `R`
    with a representation and a set of virtual characters `C`.

    INPUT:
        - `O` -- A monoid with an action of a group; As implemented in
                 :class:~`fourier_expansion_framework.monoidpowerseries.NNMonoid`.
        - `C` -- A monoid of characters; As implemented in ::class:~`fourier_expansion_framework.monoidpowerseries.CharacterMonoid_class`.
        - `R` -- A representation on a module; As implemented
                 in :class:~`fourier_expansion_framework.monoidpowerseries.TrivialRepresentation`.
        
    EXAMPLES::
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import EquivariantMonoidPowerSeriesModule
        sage: emps = EquivariantMonoidPowerSeriesModule(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", FreeModule(QQ, 2)))
        sage: emps is EquivariantMonoidPowerSeriesModule(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", FreeModule(QQ, 2)))
        True
    """

    ## TODO: Implement optional checking of the relation (*)
    if O.group() != C.group() :
        raise ValueError, "The action on S and the characters must have the same group"
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
                raise ValueError, "character codomain and representation base ring have no common extension"
    
    global _equivariantmonoidpowerseries_module_cache
    key = (O, C, R)
    
    try :
        return _equivariantmonoidpowerseries_module_cache[key]
    except KeyError :
        P = EquivariantMonoidPowerSeriesModule_generic(O, C, R)
        _equivariantmonoidpowerseries_module_cache[key] = P
        
        return P

#===============================================================================
# EquivariantMonoidPowerSeriesModule_generic
#===============================================================================

class EquivariantMonoidPowerSeriesModule_generic ( EquivariantMonoidPowerSeriesAmbient_abstract, Module ) :
    """
    Given some module `A`, a monoid `S` filtered over some originated
    net `\Lambda` such that all induced submonoids are finite, a group `G`, a
    semigroup `C` with a map `c \rightarrow \mathrm{Hom}(G, Aut_K(A))`, a
    homomorphism `\phi : G -> Aut(S)` and a homomorphism `\eta : G -> C`, where
    `K` is the base ring of `A`.
    
    Suppose for every `c, c'` in `C`, and `g` in `G`, and `a, a'` in `A` we have
      `(c c') (g) (a a') = c(g)(a) c'(g)(a')`. 
    Set `R = B[C][S]`. Then the projective limit of
      `R / R_\lambda` for `\lambda \in \Lambda \rightarrow \infinity` is a
      `K`-module.
      
    The set of generators is the set of generators of the underlying
    monoidal power series module and does not take into account the
    group action
    """
        
    def __init__(self, O, C, R) :
        """
        INPUT:
            - `O` -- A monoid with an action of a group; As implemented in
                     :class:~`fourier_expansion_framework.monoidpowerseries.NNMonoid`.
            - `C` -- A monoid of characters; As implemented in ::class:~`fourier_expansion_framework.monoidpowerseries.CharacterMonoid_class`.
            - `R` -- A representation on a module; As implemented
                     in :class:~`fourier_expansion_framework.monoidpowerseries.TrivialRepresentation`.

        EXAMPLES::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import EquivariantMonoidPowerSeriesModule
            sage: emps = EquivariantMonoidPowerSeriesModule_generic(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", FreeModule(QQ, 2))) # indirect doctest
        """
        
        # If the representation O respects the monoid structure of S
        # the base ring should be the associated power series ring.
        if O.is_monoid_action() :
            Module.__init__(self, EquivariantMonoidPowerSeriesRing(O,C,TrivialRepresentation(R.group(), R.base_ring())))
        else :
            Module.__init__(self, R.codomain())
        EquivariantMonoidPowerSeriesAmbient_abstract.__init__(self, O, C, R)        
        
        self.__coeff_gens = \
          [self._element_class( self,
            dict([( C.one_element(), dict([(self.monoid().zero_element(), a)]) )]),
            self.monoid().filter_all() )
           for a in self.coefficient_domain().gens()]
          
    def ngens(self) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import EquivariantMonoidPowerSeriesModule
            sage: emps = EquivariantMonoidPowerSeriesModule_generic(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", FreeModule(QQ, 2)))
            sage: emps.ngens()
            2
        """
        return len(self.__coeff_gens)
    
    def gen(self, i = 0) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import EquivariantMonoidPowerSeriesModule
            sage: emps = EquivariantMonoidPowerSeriesModule_generic(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", FreeModule(QQ, 2)))
            sage: emps.gen()
            Equivariant monoid power series in Module of equivariant monoid power series over NN
        """
        return self.gens()[i]
    
    def gens(self) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import EquivariantMonoidPowerSeriesModule
            sage: emps = EquivariantMonoidPowerSeriesModule_generic(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", FreeModule(QQ, 2)))
            sage: emps.gens()
            [Equivariant monoid power series in Module of equivariant monoid power series over NN, Equivariant monoid power series in Module of equivariant monoid power series over NN]
        """
        return self.__coeff_gens
    
    def construction(self) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import EquivariantMonoidPowerSeriesModule_generic
            sage: emps = EquivariantMonoidPowerSeriesModule_generic(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", FreeModule(QQ, 2)))
            sage: (f, a) = emps.construction()
            sage: (f, a)
            (EquivariantMonoidPowerSeriesFunctor, Rational Field)
            sage: f(a) == emps
            True
        """
        return EquivariantMonoidPowerSeriesFunctor(self.action(), self.characters(), self.representation()), \
               self.coefficient_domain().base_ring()
    
    def zero_element(self) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import EquivariantMonoidPowerSeriesModule
            sage: emps = EquivariantMonoidPowerSeriesModule_generic(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", FreeModule(QQ, 2)))
            sage: h = emps.zero_element()
        """
        return self(0)
    
    def _element_constructor_(self, x) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import EquivariantMonoidPowerSeriesModule
            sage: emps = EquivariantMonoidPowerSeriesModule_generic(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", FreeModule(QQ, 2)))
            sage: h = emps(0) # indirect doctest
            sage: h = emps(int(0)) # indirect doctest
        """
        if isinstance(x, int) and x == 0 :
            return self._element_class( self,
                dict( [(self.characters().one_element(), dict())] ),
                self.action().filter_all() )
        elif isinstance(x, Element) and x.is_zero() :
            P = x.parent()
            
            if self.action().is_monoid_action() and \
               self.base_ring().base_ring().has_coerce_map_from(P) or \
               not self.action().is_monoid_action() and \
               self.base_ring().has_coerce_map_from(P) :
                return self._element_class( self,
                    dict( [(self.characters().one_element(), dict())] ),
                    self.action().filter_all() )
    
        return EquivariantMonoidPowerSeriesAmbient_abstract._element_constructor_(self, x)

    def _repr_(self) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import EquivariantMonoidPowerSeriesModule
            sage: EquivariantMonoidPowerSeriesModule_generic(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", FreeModule(QQ, 2)))
            Module of equivariant monoid power series over NN
        """
        return "Module of equivariant monoid power series over " + self.monoid()._repr_()

    def _latex_(self) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import EquivariantMonoidPowerSeriesModule
            sage: latex( EquivariantMonoidPowerSeriesModule_generic(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", FreeModule(QQ, 2))) )
            Module of equivariant monoid power series over $\NN$
        """
        return "Module of equivariant monoid power series over " + self.monoid()._latex_()
