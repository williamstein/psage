r"""
Functor creating rings or modules of (equivariant) monoid power series. 

AUTHOR :
    -- Martin Raum (2009 - 07 - 25) Initial version
"""
from __future__ import absolute_import

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

from sage.categories.commutative_additive_groups import CommutativeAdditiveGroups
from sage.categories.rings import Rings
from sage.categories.modules import Modules
from sage.categories.morphism import Morphism
from sage.categories.pushout import pushout, ConstructionFunctor
from sage.rings.ring import Ring
from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import MonoidPowerSeriesRing

#===============================================================================
# MonoidPowerSeriesRingFunctor
#===============================================================================

class MonoidPowerSeriesRingFunctor ( ConstructionFunctor ) :
    r"""
    Functor mapping a coefficient ring to a monoid power series ring
    over a given monoid.
    """
    
    rank = 9
    
    def __init__(self, S) :
        r"""
        INPUT:
            - `S` -- A monoid as in :class:`~from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids.NNMonoid`
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import MonoidPowerSeriesRingFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: F = MonoidPowerSeriesRingFunctor(NNMonoid(False))
        """
        self.__S = S
        
        ConstructionFunctor.__init__(self, Rings(), Rings())
        
    def monoid(self) :
        r"""
        Return the monoid associated to this functor.
        
        OUTPUT:
            A monoid as in :class:`~from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids.NNMonoid`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import MonoidPowerSeriesRingFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: F = MonoidPowerSeriesRingFunctor(NNMonoid(False))
            sage: F.monoid()
            NN
        """
        return self.__S

    def __call__(self, A) :
        r"""
        Apply a ``self`` to a coefficient domain `A`.
        
        INPUT:
            - `A` -- A ring or module. The domain of coefficients for
                     a ring or module of monoid power series.
        
        OUTPUT:
            An instance of :class:`~from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ambient.MonoidPowerSeriesAmbient_abstract`.
            A ring if `A` is a ring, a module if `A` is a module.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import MonoidPowerSeriesRingFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: F = MonoidPowerSeriesRingFunctor(NNMonoid(False))
            sage: mps = F(ZZ)
            sage: mps.monoid() == NNMonoid(False)
            True
            sage: mps.coefficient_domain()
            Integer Ring
        """
        from .monoidpowerseries_ring import MonoidPowerSeriesRing

        return MonoidPowerSeriesRing(A, self.__S)
        
    def __cmp__(self, other) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import MonoidPowerSeriesRingFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: F = MonoidPowerSeriesRingFunctor(NNMonoid(False))
            sage: F == MonoidPowerSeriesRingFunctor(NNMonoid(False))
            True
        """
        c = cmp(type(self), type(other))
    
        if c == 0 :
            c = cmp(self.__S, other.__S)
            
        return c

#===============================================================================
# MonoidPowerSeriesModuleFunctor
#===============================================================================

class MonoidPowerSeriesModuleFunctor ( ConstructionFunctor ) :
    r"""
    Functor mapping a coefficient module to a monoid power series module
    over a given monoid.
    """
    
    rank = 9
    
    def __init__(self, B, S) :
        r"""
        INPUT:
            - `B` -- A ring.
            - `S` -- A monoid as in :class:`~from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids.NNMonoid`
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import MonoidPowerSeriesModuleFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: F = MonoidPowerSeriesModuleFunctor(ZZ, NNMonoid(False))
        """
        self.__S = S
        
        ConstructionFunctor.__init__(self, Modules(B), CommutativeAdditiveGroups())
        
    def monoid(self) :
        r"""
        Return the monoid associated to this functor.
        
        OUTPUT:
            A monoid as in :class:`~from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids.NNMonoid`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import MonoidPowerSeriesModuleFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: F = MonoidPowerSeriesModuleFunctor(ZZ, NNMonoid(False))
            sage: F.monoid()
            NN
        """
        return self.__S

    def __call__(self, A) :
        r"""
        Apply a ``self`` to a coefficient domain `A`.
        
        INPUT:
            - `A` -- A ring or module. The domain of coefficients for
                     a ring or module of monoid power series.
        
        OUTPUT:
            An instance of :class:`~from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ambient.MonoidPowerSeriesAmbient_abstract`.
            A ring if `A` is a ring, a module if `A` is a module.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import MonoidPowerSeriesModuleFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: F = MonoidPowerSeriesModuleFunctor(ZZ, NNMonoid(False))
            sage: mps = F(FreeModule(ZZ, 1))
            sage: mps.monoid() == NNMonoid(False)
            True
            sage: mps.coefficient_domain()
            Ambient free module of rank 1 over the principal ideal domain Integer Ring
            sage: F(FreeModule(QQ, 3)).coefficient_domain()
            Vector space of dimension 3 over Rational Field
        """
        from .monoidpowerseries_module import MonoidPowerSeriesModule
            
        return MonoidPowerSeriesModule(A, self.__S)
        
    def __cmp__(self, other) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import MonoidPowerSeriesModuleFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: F = MonoidPowerSeriesModuleFunctor(ZZ, NNMonoid(False))
            sage: F == MonoidPowerSeriesModuleFunctor(ZZ, NNMonoid(False))
            True
            sage: F == MonoidPowerSeriesModuleFunctor(QQ, NNMonoid(False))
            False
        """
        c = cmp(type(self), type(other))
        
        if c == 0 :
            c = cmp(self.domain().base_ring(), other.domain().base_ring())
        if c == 0 :
            c = cmp(self.__S, other.__S)
            
        return c
    
#===============================================================================
# EquivariantMonoidPowerSeriesRingFunctor
#===============================================================================

class EquivariantMonoidPowerSeriesRingFunctor ( ConstructionFunctor ) :
    r"""
    Functor mapping a coefficient domain to a equivariant monoid power series
    ring or module over a given monoid action with character and representation.
    The representation will be extended by scalars if the coefficient domain's
    base ring is to big.
    """
    
    rank = 9
    
    def __init__(self, O, C, R) :
        r"""
        INPUT:
            - `O` -- An action of a group `G` on a monoid as implemented in :class:`~fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basismonoids.NNMonoid`.
            - `C` -- A monoid of charcters `G -> K` for a ring `K`. As implemented in :class:`~fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basismonoids.CharacterMonoid_class`.
            - `R` -- A representation of `G` on some `K`-algebra or module `A`.
                     As implemented in :class:`~fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basismonoids.TrivialRepresentation`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import EquivariantMonoidPowerSeriesRingFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid, TrivialRepresentation, TrivialCharacterMonoid
            sage: F = EquivariantMonoidPowerSeriesRingFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
        """
        if O.group() != C.group() :
            raise ValueError( "The action on S and the characters must have the same group" )
        if R.base_ring() != C.codomain() :
#            if C.codomain().has_coerce_map_from(R.base_ring()) :
#                pass
#            el
            if R.base_ring().has_coerce_map_from(C.codomain()) :
                pass
            else :
                raise ValueError( "character codomain and representation base ring must be coercible" )
        if not O.is_monoid_action() :
            raise ValueError( "monoid structure must be compatible with group action" )
        
        self.__O = O
        self.__C = C
        self.__R = R
        
        ConstructionFunctor.__init__(self, Rings(), Rings())

    def __call__(self, K) :
        r"""
        Apply a ``self`` to a coefficient domain `A`.
        
        INPUT:
            - `A` -- A ring or module. The domain of coefficients for
                     a ring or module of equivariant monoid power series.
        
        OUTPUT:
            An instance of :class:`~from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ambient.EquivariantMonoidPowerSeriesAmbient_abstract`.
            A ring if the representation's extension by `A` is a ring, a module if this extension is a module.
            
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import EquivariantMonoidPowerSeriesRingFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid, TrivialRepresentation, TrivialCharacterMonoid
            sage: F = EquivariantMonoidPowerSeriesRingFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: emps = F(QQ)
            sage: emps.action() == NNMonoid()
            True
            sage: emps.characters()
            Character monoid over Trivial monoid
            """
        if not self.__R.base_ring().has_coerce_map_from(K) : 
            R = self.__R.base_extend( pushout(self.__R.base_ring(), K) )
        else :
            R = self.__R
        
        from .monoidpowerseries_ring import EquivariantMonoidPowerSeriesRing

        return EquivariantMonoidPowerSeriesRing( self.__O, self.__C, R )

    def __cmp__(self, other) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import EquivariantMonoidPowerSeriesRingFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid, TrivialRepresentation, TrivialCharacterMonoid
            sage: F = EquivariantMonoidPowerSeriesRingFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: F == EquivariantMonoidPowerSeriesRingFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            True
            sage: F == EquivariantMonoidPowerSeriesRingFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", FreeModule(QQ, 3)))
            False
        """
        c = cmp(type(self), type(other))
        if c == 0 :
            c = cmp(self.__R, other.__R)
        if c == 0 :
            c = cmp(self.__O, other.__O)
        if c == 0 :
            c = cmp(self.__C, other.__C)
        
        return c
    
    def expand(self) :
        r"""
        An equivariant monoid power series can be constructed by first
        constructing a monoid power series and then symmetrising it with
        respect to the group action.
        
        OUTPUT:
            A list of functors.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import EquivariantMonoidPowerSeriesRingFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid, TrivialRepresentation, TrivialCharacterMonoid
            sage: F = EquivariantMonoidPowerSeriesRingFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: F.expand()
            [MonoidPowerSeriesSymmetrisationRingFunctor, MonoidPowerSeriesRingFunctor]
        """
        return [ MonoidPowerSeriesSymmetrisationRingFunctor(self.__O, self.__C, self.__R),
                 MonoidPowerSeriesRingFunctor(self.__O.monoid()) ]
        
    def merge(self, other) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import EquivariantMonoidPowerSeriesRingFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid, TrivialRepresentation, TrivialCharacterMonoid
            sage: F = EquivariantMonoidPowerSeriesRingFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: G = EquivariantMonoidPowerSeriesRingFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", QQ))
            sage: F.merge(G) == G
            True
            sage: G = EquivariantMonoidPowerSeriesRingFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", FreeModule(ZZ, 3)))
            sage: F.merge(G) is None
            True
        """
        if self == other :
            return self
        elif type(self) == type(other) and \
           self.__O == self.__O and \
           self.__C == self.__C :
            try :
                if self.__R.extends(other.__R) :
                    return self
                elif other.__R.extends(self.__R) :
                    return other
            except AttributeError :
                return None
        
        return None

#===============================================================================
# EquivariantMonoidPowerSeriesFunctor
#===============================================================================

class EquivariantMonoidPowerSeriesModuleFunctor ( ConstructionFunctor ) :
    r"""
    Functor mapping a coefficient domain to a equivariant monoid power series
    ring or module over a given monoid action with character and representation.
    The representation will be extended by scalars if the coefficient domain's
    base ring is to big.
    """
    
    rank = 9
    
    def __init__(self, O, C, R) :
        r"""
        INPUT:
            - `O` -- An action of a group `G` on a monoid as implemented in :class:`~fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basismonoids.NNMonoid`.
            - `C` -- A monoid of charcters `G -> K` for a ring `K`. As implemented in :class:`~fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basismonoids.CharacterMonoid_class`.
            - `R` -- A representation of `G` on some `K`-algebra or module `A`.
                     As implemented in :class:`~fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basismonoids.TrivialRepresentation`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import EquivariantMonoidPowerSeriesModuleFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid, TrivialRepresentation, TrivialCharacterMonoid
            sage: F = EquivariantMonoidPowerSeriesModuleFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
        """
        if O.group() != C.group() :
            raise ValueError( "The action on S and the characters must have the same group" )
        if R.base_ring() != C.codomain() :
#            if C.codomain().has_coerce_map_from(R.base_ring()) :
#                pass
#            el
            if R.base_ring().has_coerce_map_from(C.codomain()) :
                pass
            else :
                raise ValueError( "character codomain and representation base ring must be coercible" )
        
        self.__O = O
        self.__C = C
        self.__R = R
        
        ConstructionFunctor.__init__(self, Modules(R.base_ring()), CommutativeAdditiveGroups())

    def __call__(self, K) :
        r"""
        Apply a ``self`` to a coefficient domain `A`.
        
        INPUT:
            - `A` -- A ring or module. The domain of coefficients for
                     a ring or module of equivariant monoid power series.
        
        OUTPUT:
            An instance of :class:`~from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ambient.EquivariantMonoidPowerSeriesAmbient_abstract`.
            A ring if the representation's extension by `A` is a ring, a module if this extension is a module.
            
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import EquivariantMonoidPowerSeriesModuleFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid, TrivialRepresentation, TrivialCharacterMonoid
            sage: F = EquivariantMonoidPowerSeriesModuleFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", FreeModule(ZZ, 1)))
            sage: emps = F(QQ)
            sage: emps.action() == NNMonoid()
            True
            sage: emps.characters()
            Character monoid over Trivial monoid
            sage: F = EquivariantMonoidPowerSeriesModuleFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", FreeModule(ZZ, 3)))
            sage: emps = F(QQ)
            sage: emps.coefficient_domain()
            Vector space of dimension 3 over Rational Field
            """
        if not self.__R.base_ring().has_coerce_map_from(K) : 
            R = self.__R.base_extend( pushout(self.__R.base_ring(), K) )
        else :
            R = self.__R
        
        from .monoidpowerseries_module import EquivariantMonoidPowerSeriesModule

        return EquivariantMonoidPowerSeriesModule( self.__O, self.__C, R )

    def __cmp__(self, other) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import EquivariantMonoidPowerSeriesModuleFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid, TrivialRepresentation, TrivialCharacterMonoid
            sage: F = EquivariantMonoidPowerSeriesModuleFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", FreeModule(ZZ, 1)))
            sage: F == EquivariantMonoidPowerSeriesModuleFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", FreeModule(ZZ, 1)))
            True
            sage: F == EquivariantMonoidPowerSeriesModuleFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", FreeModule(QQ, 3)))
            False
        """
        c = cmp(type(self), type(other))
        if c == 0 :
            c = cmp(self.__R, other.__R)
        if c == 0 :
            c = cmp(self.__O, other.__O)
        if c == 0 :
            c = cmp(self.__C, other.__C)
        
        return c
    
    def expand(self) :
        r"""
        An equivariant monoid power series can be constructed by first
        constructing a monoid power series and then symmetrising it with
        respect to the group action.
        
        OUTPUT:
            A list of functors.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import EquivariantMonoidPowerSeriesModuleFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid, TrivialRepresentation, TrivialCharacterMonoid
            sage: F = EquivariantMonoidPowerSeriesModuleFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", FreeModule(ZZ, 1)))
            sage: F.expand()
            [MonoidPowerSeriesSymmetrisationModuleFunctor, MonoidPowerSeriesModuleFunctor]
        """
        return [ MonoidPowerSeriesSymmetrisationModuleFunctor(self.__O, self.__C, self.__R),
                 MonoidPowerSeriesModuleFunctor(self.domain().base_ring(), self.__O.monoid()) ]
        
    def merge(self, other) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import EquivariantMonoidPowerSeriesModuleFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid, TrivialRepresentation, TrivialCharacterMonoid
            sage: F = EquivariantMonoidPowerSeriesModuleFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", FreeModule(ZZ, 1)))
            sage: G = EquivariantMonoidPowerSeriesModuleFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", FreeModule(QQ, 1)))
            sage: F.merge(G) == G
            True
            sage: G = EquivariantMonoidPowerSeriesModuleFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", FreeModule(ZZ, 3)))
            sage: F.merge(G) is None
            True
        """
        if self == other :
            return self
        elif type(self) == type(other) and \
           self.__O == self.__O and \
           self.__C == self.__C :
            try :
                if self.__R.extends(other.__R) :
                    return self
                elif other.__R.extends(self.__R) :
                    return other
            except AttributeError :
                return None
        
        return None

#===============================================================================
# MonoidPowerSeriesSymmetrisationRingFunctor
#===============================================================================

class MonoidPowerSeriesSymmetrisationRingFunctor ( ConstructionFunctor) :
    r"""
    A functor mapping rings or modules of monoid power series
    to a the an equivariant power series via symmetrisation.
    """
    
    rank = 9
    
    def __init__(self, O, C, R) :
        r"""
        INPUT:
            - `O` -- An action of a group `G` on a monoid as implemented in :class:`~fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basismonoids.NNMonoid`.
            - `C` -- A monoid of charcters `G -> K` for a ring `K`. As implemented in :class:`~fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basismonoids.CharacterMonoid_class`.
            - `R` -- A representation of `G` on some `K`-algebra or module `A`.
                     As implemented in :class:`~fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basismonoids.TrivialRepresentation`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import MonoidPowerSeriesSymmetrisationRingFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid, TrivialRepresentation, TrivialCharacterMonoid
            sage: F = MonoidPowerSeriesSymmetrisationRingFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
        """
        if O.group() != C.group() :
            raise ValueError( "The action on S and the characters must have the same group" )
        if R.base_ring() != C.codomain() :
#            if C.codomain().has_coerce_map_from(R.base_ring()) :
#                pass
#            el
            if R.base_ring().has_coerce_map_from(C.codomain()) :
                pass
            else :
                raise ValueError( "character codomain and representation base ring must be coercible" )
        if not O.is_monoid_action() :
            raise ValueError( "monoid structure must be compatible with group action" )
        
        self.__O = O
        self.__C = C
        self.__R = R
        
        ConstructionFunctor.__init__(self, Rings(), Rings())
        
    def __call__(self, P) :
        r"""
        Map a monoid power series to a symmetrisation extending the associated representation.

        INPUT:
            - `P` -- An instance of :class:`~from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ambient.MonoidPowerSeriesAmbient_abstract`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import MonoidPowerSeriesSymmetrisationRingFunctor, MonoidPowerSeriesRingFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid, TrivialRepresentation, TrivialCharacterMonoid
            sage: mps = MonoidPowerSeriesRingFunctor(NNMonoid(False))(QQ)
            sage: F = MonoidPowerSeriesSymmetrisationRingFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: emps = F(mps)
            sage: emps.representation() == TrivialRepresentation("1", QQ)
            True
            sage: mps = MonoidPowerSeriesRingFunctor(NNMonoid(False))(ZZ)
            sage: F = MonoidPowerSeriesSymmetrisationRingFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", QQ))
            sage: emps = F(mps)
            sage: emps.representation() == TrivialRepresentation("1", QQ)
            True
        """
        if self.__O.monoid() != P.monoid() :
            raise ValueError( "Action has to be defined on the monoid associated to P." )
        
        PR = self.__R.from_module(P.coefficient_domain())
        
        if not self.__R.base_ring().has_coerce_map_from(PR.base_ring()) : 
            R = self.__R.base_extend( pushout(self.__R.base_ring(), PR.base_ring()) )
        else :
            R = self.__R
        
        from .monoidpowerseries_ring import EquivariantMonoidPowerSeriesRing

        return EquivariantMonoidPowerSeriesRing( self.__O, self.__C, R )
        
    def __cmp__(self, other) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import MonoidPowerSeriesSymmetrisationRingFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid, TrivialRepresentation, TrivialCharacterMonoid
            sage: F = MonoidPowerSeriesSymmetrisationRingFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: F == MonoidPowerSeriesSymmetrisationRingFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            True
            sage: F == MonoidPowerSeriesSymmetrisationRingFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", QQ))
            False
        """
        c = cmp(type(self), type(other))
        if c == 0 :
            c = cmp(self.__R, other.__R)
        if c == 0 :
            c = cmp(self.__O, other.__O)
        if c == 0 :
            c = cmp(self.__C, other.__C)
           
        return c

    def merge(self, other) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import MonoidPowerSeriesSymmetrisationRingFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid, TrivialRepresentation, TrivialCharacterMonoid
            sage: F = MonoidPowerSeriesSymmetrisationRingFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: G = MonoidPowerSeriesSymmetrisationRingFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", QQ))
            sage: F.merge(G) == G
            True
            sage: G = MonoidPowerSeriesSymmetrisationRingFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", FreeModule(QQ, 3)))
            sage: F.merge(G) is None
            True
        """
        if self == other :
            return self
        elif type(self) == type(other) and \
           self.__O == self.__O :
            try :
                if self.__C.extends(other.__C) and self.__R.extends(other.__R) :
                    return self
                elif other.__C.extends(self.__C) and other.__R.extends(self.__R) :
                    return other
            except AttributeError :
                return None
        
        return None

#===============================================================================
# MonoidPowerSeriesSymmetrisationModuleFunctor
#===============================================================================

class MonoidPowerSeriesSymmetrisationModuleFunctor ( ConstructionFunctor) :
    r"""
    A functor mapping rings or modules of monoid power series
    to a the an equivariant power series via symmetrisation.
    """
    
    rank = 9
    
    def __init__(self, O, C, R) :
        r"""
        INPUT:
            - `O` -- An action of a group `G` on a monoid as implemented in :class:`~fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basismonoids.NNMonoid`.
            - `C` -- A monoid of charcters `G -> K` for a ring `K`. As implemented in :class:`~fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basismonoids.CharacterMonoid_class`.
            - `R` -- A representation of `G` on some `K`-algebra or module `A`.
                     As implemented in :class:`~fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basismonoids.TrivialRepresentation`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import MonoidPowerSeriesSymmetrisationModuleFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid, TrivialRepresentation, TrivialCharacterMonoid
            sage: F = MonoidPowerSeriesSymmetrisationModuleFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", FreeModule(ZZ, 1)))
        """
        if O.group() != C.group() :
            raise ValueError("The action on S and the characters must have the same group")
        if R.base_ring() != C.codomain() :
            if C.codomain().has_coerce_map_from(R.base_ring()) :
                pass
            elif R.base_ring().has_coerce_map_from(C.codomain()) :
                pass
            else :
                raise ValueError("character codomain and representation base ring must be coercible")        
        
        self.__O = O
        self.__C = C
        self.__R = R
        
        ConstructionFunctor.__init__(self, CommutativeAdditiveGroups(), CommutativeAdditiveGroups())
        
    def __call__(self, P) :
        r"""
        Map a monoid power series to a symmetrisation extending the associated representation.

        INPUT:
            - `P` -- An instance of :class:`~from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ambient.MonoidPowerSeriesAmbient_abstract`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import MonoidPowerSeriesSymmetrisationModuleFunctor, MonoidPowerSeriesModuleFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid, TrivialRepresentation, TrivialCharacterMonoid
            sage: mps = MonoidPowerSeriesModuleFunctor(ZZ, NNMonoid(False))(FreeModule(QQ, 1))
            sage: F = MonoidPowerSeriesSymmetrisationModuleFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", FreeModule(ZZ, 1)))
            sage: emps = F(mps)
            sage: emps.representation() == TrivialRepresentation("1", FreeModule(QQ, 1))
            True
            sage: mps = MonoidPowerSeriesModuleFunctor(ZZ, NNMonoid(False))(FreeModule(ZZ, 1))
            sage: F = MonoidPowerSeriesSymmetrisationModuleFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", FreeModule(QQ, 1)))
            sage: emps = F(mps)
            sage: emps.representation() == TrivialRepresentation("1", FreeModule(QQ, 1))
            True
        """
        if self.__O.monoid() != P.monoid() :
            raise ValueError( "Action has to be defined on the monoid associated to P." )
        
        PR = self.__R.from_module(P.coefficient_domain())
        
        if not self.__R.base_ring().has_coerce_map_from(PR.base_ring()) : 
            R = self.__R.base_extend( pushout(self.__R.base_ring(), PR.base_ring()) )
        else :
            R = self.__R
        
        from .monoidpowerseries_module import EquivariantMonoidPowerSeriesModule

        return EquivariantMonoidPowerSeriesModule( self.__O, self.__C, R )

    def __cmp__(self, other) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import MonoidPowerSeriesSymmetrisationModuleFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid, TrivialRepresentation, TrivialCharacterMonoid
            sage: F = MonoidPowerSeriesSymmetrisationModuleFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", FreeModule(ZZ, 1)))
            sage: F == MonoidPowerSeriesSymmetrisationModuleFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", FreeModule(ZZ, 1)))
            True
            sage: F == MonoidPowerSeriesSymmetrisationModuleFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", FreeModule(QQ, 1)))
            False
        """
        c = cmp(type(self), type(other))
        if c == 0 :
            c = cmp(self.__R, other.__R)
        if c == 0 :
            c = cmp(self.__O, other.__O)
        if c == 0 :
            c = cmp(self.__C, other.__C)
           
        return c

    def merge(self, other) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import MonoidPowerSeriesSymmetrisationModuleFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid, TrivialRepresentation, TrivialCharacterMonoid
            sage: F = MonoidPowerSeriesSymmetrisationModuleFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", FreeModule(ZZ, 1)))
            sage: G = MonoidPowerSeriesSymmetrisationModuleFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", FreeModule(QQ, 1)))
            sage: F.merge(G) == G
            True
            sage: G = MonoidPowerSeriesSymmetrisationModuleFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", FreeModule(QQ, 3)))
            sage: F.merge(G) is None
            True
        """
        if self == other :
            return self
        elif type(self) == type(other) and \
           self.__O == self.__O :
            try :
                if self.__C.extends(other.__C) and self.__R.extends(other.__R) :
                    return self
                elif other.__C.extends(self.__C) and other.__R.extends(self.__R) :
                    return other
            except AttributeError :
                return None
        
        return None
        
#===============================================================================
# MonoidPowerSeriesBaseRingInjection
#===============================================================================

class MonoidPowerSeriesBaseRingInjection ( Morphism ) :
    r"""
    The injection of the base ring into a (equivariant) monoid power
    series ring.
    """
    
    def __init__(self, domain, codomain) :
        r"""
        INPUT:
            - ``domain``   -- A ring; The base ring.
            - ``codomain`` -- A ring; The ring of monoid power series.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import MonoidPowerSeriesRingFunctor, MonoidPowerSeriesModuleFunctor, MonoidPowerSeriesBaseRingInjection
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: mps = MonoidPowerSeriesRingFunctor(NNMonoid(False))(ZZ)
            sage: binj = MonoidPowerSeriesBaseRingInjection(ZZ, mps)
            sage: mps = MonoidPowerSeriesModuleFunctor(ZZ,NNMonoid(False))(FreeModule(ZZ, 1))
            sage: binj = MonoidPowerSeriesBaseRingInjection(ZZ, mps)

        """
        Morphism.__init__(self, domain, codomain)
        
        self._repr_type_str = "MonoidPowerSeries base injection"

    def _call_(self, x) :
        r"""
        Coerce an element into the ring of monoid power series.
        
        INPUT:
            - `x` -- An element of a ring; An element of the base ring.
        
        OUTPUT:
            An element of a ring. 
            
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import MonoidPowerSeriesRingFunctor, MonoidPowerSeriesBaseRingInjection
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: mps = MonoidPowerSeriesRingFunctor(NNMonoid(False))(ZZ)
            sage: binj = MonoidPowerSeriesBaseRingInjection(ZZ, mps)
            sage: binj(1)
            Monoid power series in Ring of monoid power series over NN
        """
        return self.codomain()._element_constructor_(x)
            
    def _call_with_args(self, x, *args, **kwds):
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import MonoidPowerSeriesRingFunctor, MonoidPowerSeriesBaseRingInjection
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: mps = MonoidPowerSeriesRingFunctor(NNMonoid(False))(ZZ)
            sage: binj = MonoidPowerSeriesBaseRingInjection(ZZ, mps)
            sage: binj._call_with_args(1)
            Monoid power series in Ring of monoid power series over NN
        """
        return self.codomain()._element_constructor_(x, *args, **kwds)
