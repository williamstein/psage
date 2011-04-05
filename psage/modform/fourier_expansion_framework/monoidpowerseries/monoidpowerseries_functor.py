"""
Functor creating rings or modules of (equivariant) monoid power series. 

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

from sage.categories.rings import Rings
from sage.categories.morphism import Morphism
from sage.categories.pushout import ConstructionFunctor
from sage.rings.ring import Ring

#===============================================================================
# MonoidPowerSeriesFunctor
#===============================================================================

class MonoidPowerSeriesFunctor ( ConstructionFunctor ) :
    """
    Functor mapping a coefficient ring to a monoid power series ring
    over a given monoid.
    """
    
    rank = 9
    
    def __init__(self, S) :
        """
        INPUT:
            - `S` -- A monoid as in :class:`~from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids.NNMonoid`
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import MonoidPowerSeriesFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: F = MonoidPowerSeriesFunctor(NNMonoid(False))
        """
        self.__S = S
        
        ConstructionFunctor.__init__(self, Rings(), Rings())
        
    def monoid(self) :
        """
        Return the monoid associated to this functor.
        
        OUTPUT:
            A monoid as in :class:`~from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids.NNMonoid`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import MonoidPowerSeriesFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: F = MonoidPowerSeriesFunctor(NNMonoid(False))
            sage: F.monoid()
            NN
        """
        return self.__S

    def __call__(self, A) :
        """
        Apply a ``self`` to a coefficient domain `A`.
        
        INPUT:
            - `A` -- A ring or module. The domain of coefficients for
                     a ring or module of monoid power series.
        
        OUTPUT:
            An instance of :class:`~from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ambient.MonoidPowerSeriesAmbient_abstract`.
            A ring if `A` is a ring, a module if `A` is a module.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import MonoidPowerSeriesFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: F = MonoidPowerSeriesFunctor(NNMonoid(False))
            sage: mps = F(ZZ)
            sage: mps.monoid() == NNMonoid(False)
            True
            sage: mps.coefficient_domain()
            Integer Ring
            sage: F(FreeModule(QQ, 3)).coefficient_domain()
            Vector space of dimension 3 over Rational Field
        """
        if isinstance(A, Ring) :
            from monoidpowerseries_ring import MonoidPowerSeriesRing

            return MonoidPowerSeriesRing(A, self.__S)
        else :
            from monoidpowerseries_module import MonoidPowerSeriesModule
            
            return MonoidPowerSeriesModule(A, self.__S)
        
    def __cmp__(self, other) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import MonoidPowerSeriesFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: F = MonoidPowerSeriesFunctor(NNMonoid(False))
            sage: F == MonoidPowerSeriesFunctor(NNMonoid(False))
            True
        """
        c = cmp(type(self), type(other))
    
        if c == 0 :
            c = cmp(self.__S, other.__S)
            
        return c
    
#===============================================================================
# EquivariantMonoidPowerSeriesFunctor
#===============================================================================

class EquivariantMonoidPowerSeriesFunctor ( ConstructionFunctor ) :
    """
    Functor mapping a coefficient domain to a equivariant monoid power series
    ring or module over a given monoid action with character and representation.
    The representation will be extended by scalars if the coefficient domain's
    base ring is to big.
    """
    
    rank = 9
    
    def __init__(self, O, C, R) :
        """
        INPUT:
            - `O` -- An action of a group `G` on a monoid as implemented in :class:`~fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basismonoids.NNMonoid`.
            - `C` -- A monoid of charcters `G -> K` for a ring `K`. As implemented in :class:`~fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basismonoids.CharacterMonoid_class`.
            - `R` -- A representation of `G` on some `K`-algebra or module `A`.
                     As implemented in :class:`~fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basismonoids.TrivialRepresentation`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import EquivariantMonoidPowerSeriesFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid, TrivialRepresentation, TrivialCharacterMonoid
            sage: F = EquivariantMonoidPowerSeriesFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
        """
        if O.group() != C.group() :
            raise ValueError( "The action on S and the characters must have the same group" )
        if R.base_ring() != C.codomain() :
            if C.codomain().has_coerce_map_from(R.base_ring()) :
                pass
            elif R.base_ring().has_coerce_map_from(C.codomain()) :
                pass
            else :
                raise ValueError( "character codomain and representation base ring must be coercible" )
        
        self.__O = O
        self.__C = C
        self.__R = R
        
        ## TODO: replace Rings by an suitable category
        ConstructionFunctor.__init__(self, Rings(), Rings())

    def __call__(self, K) :
        """
        Apply a ``self`` to a coefficient domain `A`.
        
        INPUT:
            - `A` -- A ring or module. The domain of coefficients for
                     a ring or module of equivariant monoid power series.
        
        OUTPUT:
            An instance of :class:`~from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ambient.EquivariantMonoidPowerSeriesAmbient_abstract`.
            A ring if the representation's extension by `A` is a ring, a module if this extension is a module.
            
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import EquivariantMonoidPowerSeriesFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid, TrivialRepresentation, TrivialCharacterMonoid
            sage: F = EquivariantMonoidPowerSeriesFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: emps = F(QQ)
            sage: emps.action() == NNMonoid()
            True
            sage: emps.characters()
            Character monoid over Trivial monoid
            sage: F = EquivariantMonoidPowerSeriesFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", FreeModule(ZZ, 3)))
            sage: emps = F(QQ)
            sage: emps.coefficient_domain()
            Vector space of dimension 3 over Rational Field
            """
        R = self.__R.base_extend(K)
        if self.__O.is_monoid_action() and isinstance(R.codomain(), Ring) :
            from monoidpowerseries_ring import EquivariantMonoidPowerSeriesRing

            return EquivariantMonoidPowerSeriesRing( self.__O, self.__C, R )
        else :
            from monoidpowerseries_module import EquivariantMonoidPowerSeriesModule

            return EquivariantMonoidPowerSeriesModule( self.__O, self.__C, R )

    def __cmp__(self, other) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import EquivariantMonoidPowerSeriesFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid, TrivialRepresentation, TrivialCharacterMonoid
            sage: F = EquivariantMonoidPowerSeriesFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: F == EquivariantMonoidPowerSeriesFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            True
            sage: F == EquivariantMonoidPowerSeriesFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", FreeModule(QQ, 3)))
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
        """
        An equivariant monoid power series can be constructed by first
        constructing a monoid power series and then symmetrising it with
        respect to the group action.
        
        OUTPUT:
            A list of functors.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import EquivariantMonoidPowerSeriesFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid, TrivialRepresentation, TrivialCharacterMonoid
            sage: F = EquivariantMonoidPowerSeriesFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: F.expand()
            [MonoidPowerSeriesSymmetrisationFunctor, MonoidPowerSeriesFunctor]
        """
        return [ MonoidPowerSeriesSymmetrisationFunctor(self.__O, self.__C, self.__R),
                 MonoidPowerSeriesFunctor(self.__O.monoid()) ]
        
    def merge(self, other) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import EquivariantMonoidPowerSeriesFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid, TrivialRepresentation, TrivialCharacterMonoid
            sage: F = EquivariantMonoidPowerSeriesFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: G = EquivariantMonoidPowerSeriesFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", QQ))
            sage: F.merge(G) == G
            True
            sage: G = EquivariantMonoidPowerSeriesFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", FreeModule(ZZ, 3)))
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
# MonoidPowerSeriesSymmetrisationFunctor
#===============================================================================

class MonoidPowerSeriesSymmetrisationFunctor ( ConstructionFunctor) :
    """
    A functor mapping rings or modules of monoid power series
    to a the an equivariant power series via symmetrisation.
    """
    
    rank = 9
    
    def __init__(self, O, C, R) :
        """
        INPUT:
            - `O` -- An action of a group `G` on a monoid as implemented in :class:`~fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basismonoids.NNMonoid`.
            - `C` -- A monoid of charcters `G -> K` for a ring `K`. As implemented in :class:`~fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basismonoids.CharacterMonoid_class`.
            - `R` -- A representation of `G` on some `K`-algebra or module `A`.
                     As implemented in :class:`~fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basismonoids.TrivialRepresentation`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import MonoidPowerSeriesSymmetrisationFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid, TrivialRepresentation, TrivialCharacterMonoid
            sage: F = MonoidPowerSeriesSymmetrisationFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
        """
        if O.group() != C.group() :
            raise ValueError, "The action on S and the characters must have the same group"
        if R.base_ring() != C.codomain() :
            if C.codomain().has_coerce_map_from(R.base_ring()) :
                pass
            elif R.base_ring().has_coerce_map_from(C.codomain()) :
                pass
            else :
                raise ValueError, "character codomain and representation base ring must be coercible"        
        
        self.__O = O
        self.__C = C
        self.__R = R
        
        ## TODO: replace Rings by an suitable category
        ConstructionFunctor.__init__(self, Rings(), Rings())
        
    def __call__(self, P) :
        """
        Map a monoid power series to a symmetrisation extending the associated representation.

        INPUT:
            - `P` -- An instance of :class:`~from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ambient.MonoidPowerSeriesAmbient_abstract`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import MonoidPowerSeriesSymmetrisationFunctor, MonoidPowerSeriesFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid, TrivialRepresentation, TrivialCharacterMonoid
            sage: mps = MonoidPowerSeriesFunctor(NNMonoid(False))(QQ)
            sage: F = MonoidPowerSeriesSymmetrisationFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: emps = F(mps)
            sage: emps.representation() == TrivialRepresentation("1", QQ)
            True
            sage: mps = MonoidPowerSeriesFunctor(NNMonoid(False))(ZZ)
            sage: F = MonoidPowerSeriesSymmetrisationFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", QQ))
            sage: emps = F(mps)
            sage: emps.representation() == TrivialRepresentation("1", QQ)
            True
            sage: mps = MonoidPowerSeriesFunctor(NNMonoid(False))(FreeModule(ZZ, 3))
            sage: F(mps)
            Traceback (most recent call last):
            ...
            TypeError: no base extension defined
        """
        if self.__O.monoid() != P.monoid() :
            raise ValueError( "Action has to be defined on the monoid associated to P." )
        if not self.__R.base_ring().has_coerce_map_from(P.base_ring()) : 
            R = self.__R.base_extend(P.base_ring())
        else :
            R = self.__R
        
        if self.__O.is_monoid_action() and isinstance(R.codomain(), Ring) :
            from monoidpowerseries_ring import EquivariantMonoidPowerSeriesRing

            return EquivariantMonoidPowerSeriesRing( self.__O, self.__C, R )
        else :
            from monoidpowerseries_module import EquivariantMonoidPowerSeriesModule

            return EquivariantMonoidPowerSeriesModule( self.__O, self.__C, R )

    def __cmp__(self, other) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import MonoidPowerSeriesSymmetrisationFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid, TrivialRepresentation, TrivialCharacterMonoid
            sage: F = MonoidPowerSeriesSymmetrisationFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: F == MonoidPowerSeriesSymmetrisationFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            True
            sage: F == MonoidPowerSeriesSymmetrisationFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", QQ))
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
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import MonoidPowerSeriesSymmetrisationFunctor
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid, TrivialRepresentation, TrivialCharacterMonoid
            sage: F = MonoidPowerSeriesSymmetrisationFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: G = MonoidPowerSeriesSymmetrisationFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", QQ))
            sage: F.merge(G) == G
            True
            sage: G = MonoidPowerSeriesSymmetrisationFunctor(NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", FreeModule(QQ, 3)))
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
    """
    The injection of the base ring into a (equivariant) monoid power
    series ring.
    """
    
    def __init__(self, domain, codomain) :
        """
        INPUT:
            - ``domain``   -- A ring; The base ring.
            - ``codomain`` -- A ring; The ring of monoid power series.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import MonoidPowerSeriesFunctor, MonoidPowerSeriesBaseRingInjection
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: mps = MonoidPowerSeriesFunctor(NNMonoid(False))(ZZ)
            sage: binj = MonoidPowerSeriesBaseRingInjection(ZZ, mps)
        """
        Morphism.__init__(self, domain, codomain)
        
        self._repr_type_str = "MonoidPowerSeries base injection"

    def _call_(self, x) :
        """
        Coerce an element into the ring of monoid power series.
        
        INPUT:
            - `x` -- An element of a ring; An element of the base ring.
        
        OUTPUT:
            An element of a ring. 
            
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import MonoidPowerSeriesFunctor, MonoidPowerSeriesBaseRingInjection
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: mps = MonoidPowerSeriesFunctor(NNMonoid(False))(ZZ)
            sage: binj = MonoidPowerSeriesBaseRingInjection(ZZ, mps)
            sage: binj(1)
            Monoid power series in Ring of monoid power series over NN
        """
        return self.codomain()._element_constructor_(x)
            
    def _call_with_args(self, x, *args, **kwds):
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_functor import MonoidPowerSeriesFunctor, MonoidPowerSeriesBaseRingInjection
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: mps = MonoidPowerSeriesFunctor(NNMonoid(False))(ZZ)
            sage: binj = MonoidPowerSeriesBaseRingInjection(ZZ, mps)
            sage: binj._call_with_args(1)
            Monoid power series in Ring of monoid power series over NN
        """
        return self.codomain()._element_constructor_(x, *args, **kwds)
