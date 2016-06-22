r"""
Module abstractly spanned by Fourier expansions. 

AUTHOR:
    - Martin Raum (2010 - 05 - 15) Initial version
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

from psage.modform.fourier_expansion_framework.gradedexpansions.expansion_lazy_evaluation import LazyFourierExpansionEvaluation
from psage.modform.fourier_expansion_framework.gradedexpansions.fourierexpansionwrapper import FourierExpansionWrapper
from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import EquivariantMonoidPowerSeriesAmbient_abstract, MonoidPowerSeriesAmbient_abstract
from sage.categories.pushout import pushout
from sage.interfaces.magma import magma
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_method
from sage.misc.flatten import flatten
from sage.misc.latex import latex
from sage.misc.misc import union
from sage.modules.free_module import FreeModule, FreeModule_generic, \
      FreeModule_ambient_pid, FreeModule_submodule_pid 
from sage.modules.free_module import is_FreeModule
from sage.modules.free_module_element import FreeModuleElement_generic_dense
from sage.modules.free_module_element import vector
from sage.arith.all import random_prime
from sage.rings.integer_ring import ZZ
from sage.rings.number_field.order import Order
from sage.rings.padics.factory import Qp
from sage.rings.ring import PrincipalIdealDomain
from sage.rings.rational_field import QQ
from sage.rings.ring import Ring
from sage.structure.element import Element
from sage.structure.sequence import Sequence, Sequence_generic
import itertools

#===============================================================================
# ExpansionModule
#===============================================================================

def ExpansionModule(forms) :
    r"""
    Construct a module of over the forms' base ring of rank ``len(forms)``
    with underlying expansions associated to.
    
    INPUT:
        - ``forms`` -- A sequence or nonempty list of (equivariant) monoid power series.
    
    OUTPUT:
        An instance of :class:~`.ExpansionModule_abstract`.

    EXAMPLES::
        sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
        sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
        sage: m = FreeModule(QQ, 3)
        sage: empsm = EquivariantMonoidPowerSeriesModule(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", m))
        sage: em = ExpansionModule([emps.one_element(), emps.one_element()])
        sage: em = ExpansionModule([empsm.zero_element(), empsm.zero_element()])
        sage: mps = MonoidPowerSeriesRing(ZZ, NNMonoid())
        sage: em = ExpansionModule([mps.one_element()]) 
        sage: mpsm = MonoidPowerSeriesModule(m, NNMonoid())
        sage: em = ExpansionModule([mpsm.zero_element()]) 
    
    TESTS::
        sage: em = ExpansionModule(Sequence([], universe = emps))
        sage: em = ExpansionModule(Sequence([], universe = empsm))
        sage: h = EquivariantMonoidPowerSeries(emps, {emps.characters().one_element(): {1: 1, 3: 2}}, emps.action().filter_all())
        sage: hv = EquivariantMonoidPowerSeries(empsm, {empsm.characters().one_element(): {1: m([1,0,0]), 2: m([0,0,1])}}, empsm.action().filter_all())
        sage: em = ExpansionModule([emps.one_element(), h, h^2])
        sage: em = ExpansionModule([empsm.zero_element(), hv, hv * h])
        sage: em = ExpansionModule([])
        Traceback (most recent call last):
        ...
        ValueError: Empty modules must be constructed with a universe.
        sage: em = ExpansionModule(Sequence([], universe = m))
        Traceback (most recent call last):
        ...
        ValueError: Common parent of all forms must be a monoid power series ring or module.
        sage: qa = QuaternionAlgebra(QQ, -1, -1)
        sage: em = ExpansionModule(Sequence([], universe = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", qa)) ))
        Traceback (most recent call last):
        ...
        TypeError: The forms' base ring must be a commutative ring.
    """
    if not isinstance(forms, Sequence_generic) :
        forms = Sequence(forms)
        if len(forms) == 0 :
            raise ValueError( "Empty modules must be constructed with a universe." )
                
    if not isinstance(forms.universe(), (MonoidPowerSeriesAmbient_abstract,
                                         EquivariantMonoidPowerSeriesAmbient_abstract)) :
        raise ValueError( "Common parent of all forms must be a monoid power series ring or module." )

    if isinstance(forms.universe(), EquivariantMonoidPowerSeriesAmbient_abstract) \
       and forms.universe().representation().base_ring() == forms.universe().representation().codomain() : 
        base_ring = forms.universe().base_ring()
    elif isinstance(forms.universe(), MonoidPowerSeriesAmbient_abstract) \
       and isinstance(forms.universe().coefficient_domain(), Ring) :
        base_ring = forms.universe().base_ring()
    else :
        base_ring = forms.universe().base_ring().base_ring()
        
    if not base_ring.is_commutative():
        raise TypeError( "The forms' base ring must be a commutative ring." )

    if base_ring.is_field() \
       or isinstance(base_ring, PrincipalIdealDomain) \
       or ( isinstance(base_ring, Order) \
            and base_ring.is_maximal() and base_ring.class_number() == 1 ) :
        return ExpansionModule_ambient_pid(forms)
    else :
        return ExpansionModule_generic(forms)

#===============================================================================
# ExpansionModule_abstract
#===============================================================================

class ExpansionModule_abstract :
    r"""
    An abstract implementation of a module with expansions associated to its basis elements.
    """
    
    def __init__(self, basis, **kwds) :
        r"""
        INPUT:
            - ``basis`` -- A sequence of (equivariant) monoid power series.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.expansion_module import ExpansionModule_abstract
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
            sage: em = ExpansionModule_abstract(Sequence([emps.one_element(), emps.one_element()]))
        """
        self.__abstract_basis = basis

    def _abstract_basis(self) :
        r"""
        Return a basis in terms of elements of the graded ambient.
        
        OUTPUT:
            A sequence of (equivariant) monoid power series.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
            sage: basis = Sequence([emps.one_element(), emps.one_element()])
            sage: em = ExpansionModule(basis)
            sage: em._abstract_basis() == basis
            True
        """
        return self.__abstract_basis

    @cached_method
    def precision(self) :
        r"""
        A common precision of the expansions associated to the basis elements.
        
        OUTPUT:
            A filter for the expansions' parent's action or monoid.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
            sage: em = ExpansionModule(Sequence([emps.one_element(), emps.one_element()]))
            sage: em.precision()
            Filtered NN with action up to +Infinity
            sage: em = ExpansionModule(Sequence([emps.one_element().truncate(3), emps.one_element()]))
            sage: em.precision()
            Filtered NN with action up to 3
            sage: em = ExpansionModule(Sequence([], universe = emps))
            sage: em.precision()
            Filtered NN with action up to +Infinity
        """
        if len(self.__abstract_basis) != 0 :
            return min([b.precision() for b in self.__abstract_basis])
        else :
            if isinstance(self.__abstract_basis.universe(), EquivariantMonoidPowerSeriesAmbient_abstract) :
                return self.__abstract_basis.universe().action().filter_all()
            else :
                return self.__abstract_basis.universe().monoid().filter_all()
    
    def _bounding_precision(self) :
        r"""
        A common precision of the expansions associtated to the basis elements
        if it is finite or in case it is not a filter that comprises all
        nonzero coefficient of the expansion.
        
        OUTPUT:
            A filter for the expansions' parent's action or monoid.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
            sage: em = ExpansionModule(Sequence([emps.one_element(), emps.one_element()]))
            sage: em._bounding_precision()
            Filtered NN with action up to 1
            sage: em = ExpansionModule(Sequence([emps.one_element().truncate(4), emps.one_element()]))
            sage: em._bounding_precision()
            Filtered NN with action up to 4
            sage: em = ExpansionModule(Sequence([], universe = emps))
            sage: em._bounding_precision()
            Filtered NN with action up to 0
        """
        if len(self.__abstract_basis) != 0 :
            return max([b._bounding_precision() for b in self.__abstract_basis])
        else :
            if isinstance(self.__abstract_basis.universe(), EquivariantMonoidPowerSeriesAmbient_abstract) :
                return self.__abstract_basis.universe().action().zero_filter()
            else :
                return self.__abstract_basis.universe().monoid().zero_filter()

    @cached_method
    def _check_precision(self, precision = None, lazy_rank_check = False) :
        r"""
        Check whether the elements of this module are uniquely determined
        by their Fourier expansions up to ``precision``. If ``precision`` is ``None``
        the precision of this module will be used.
        
        INPUT:
            - ``precision``        -- A filter for the expansions' parent monoid or action or ``None`` (default: ``None``).
            - ``lazy_rank_check``  -- A boolean (default: ``False``); If ``True`` the involved rank checks will
                                      be done over `Q_p` (with finite precision).
        
        OUTPUT:
            A boolean.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
            sage: em = ExpansionModule(Sequence([emps.one_element(), emps.one_element()]))
            sage: em._check_precision()
            False
            sage: em = ExpansionModule(Sequence([], universe = emps))
            sage: em._check_precision()
            True
            sage: em._check_precision(lazy_rank_check = True)
            True
            sage: h = EquivariantMonoidPowerSeries(emps, {emps.characters().one_element(): {1: 1, 3: 2}}, emps.action().filter_all())
            sage: em = ExpansionModule(Sequence([emps.one_element(), h]))
            sage: em._check_precision()
            True
            sage: em._check_precision(emps.action().filter(1))
            False
            sage: m = FreeModule(QQ, 3)
            sage: empsm = EquivariantMonoidPowerSeriesModule(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", m))
            sage: hv = EquivariantMonoidPowerSeries(empsm, {empsm.characters().one_element(): {1: m([1,0,0])}}, empsm.action().filter_all())
            sage: hv2 = EquivariantMonoidPowerSeries(empsm, {empsm.characters().one_element(): {1: m([1,1,0])}}, empsm.action().filter_all())
            sage: em = ExpansionModule(Sequence([hv, hv2]))
            sage: em._check_precision()
            True
        """
        if lazy_rank_check and not( self.base_ring() is ZZ or self.base_ring() is QQ) :
            raise NotImplemented, "lazy rank checks only implemented for ZZ and QQ" 

        if precision is None :
            precision = self.precision()
        elif not precision <= self.precision() :
            raise ValueError, "precison must be less equal self.__graded_ambient.precision()"
        
        if precision.is_infinite() :
            precision = self._bounding_precision()
        if self.precision().is_infinite() :
            total_precision = self._bounding_precision()
        else :
            total_precision = self.precision()

        basis_fe_expansion = self.fourier_expansion_homomorphism().matrix()
        
        if precision != self.precision() :
            if isinstance(self._abstract_basis().universe().coefficient_domain(), Ring) :
                indices = [i for (i,k) in enumerate(total_precision) if k in precision]
            else :
                indices = [i for (i,(k,_)) in enumerate(
                  itertools.product(total_precision, range(self._abstract_basis().universe().coefficient_domain().rank())) )
                             if k in precision]
            basis_fe_expansion = basis_fe_expansion.matrix_from_columns(indices)
        
        if lazy_rank_check :
            basis_fe_expansion = matrix(Qp(random_prime(10**9), 10), basis_fe_expansion)
        
        return basis_fe_expansion.rank() >= self.rank()
    
    def _fourier_expansion_of_element(self, e) :
        r"""
        The Fourier expansion of an element optained via linear combinations
        of the basis' Fourier expansions.
        
        INPUT:
            - `e` -- An element of ``self``.
        
        OUTPUT:
            A (equivariant) monoid power series.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
            sage: em = ExpansionModule(Sequence([emps.one_element(), emps.one_element()]))
            sage: fe = em._fourier_expansion_of_element(em([1,0]))
            sage: m = FreeModule(QQ, 3)
            sage: empsm = EquivariantMonoidPowerSeriesModule(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", m))
            sage: hv = EquivariantMonoidPowerSeries(empsm, {empsm.characters().one_element(): {1: m([1,0,0])}}, empsm.action().filter_all())
            sage: em = ExpansionModule(Sequence([hv]))
            sage: fe = em._fourier_expansion_of_element(em([1]))            
        """
        return sum( (e[k]*b for (k,b) in enumerate(self.__abstract_basis) if e[k] != 0),
                    self.__abstract_basis.universe()(0) )

    @cached_method
    def _non_zero_characters(self) :
        r"""
        Return those characters which cannot be guaranteed to have vanishing Fourier
        expansion associated with for all basis elements of ``self``.
        
        OUTPUT:
            A list of characters of the expansions' parent if it is equivariant
            or otherwise None.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
            sage: em = ExpansionModule(Sequence([emps.one_element(), emps.one_element()]))
            sage: em._non_zero_characters()
            [1]
            sage: em = ExpansionModule(Sequence([EquivariantMonoidPowerSeries(emps, {}, emps.action().filter_all())]))
            sage: em._non_zero_characters()
            []
        """
        if isinstance(self.__abstract_basis.universe(), EquivariantMonoidPowerSeriesAmbient_abstract) :
            return list( reduce(union, [ set(b.non_zero_components())
                                         for b in flatten(self.__abstract_basis, tuple) ], set()) )
        else :
            return None 

    @cached_method
    def _fourier_expansion_indices(self) :
        r"""
        A list of Fourier indices which are considered by the Fourier expansion morphism.
        
        OUTPUT:
            A list of monoid elements or components indices and monoid elements
            in case the expansions' parent is not equivariant. Otherwise, it
            is pairs of characters and monoid elements possibly with a component
            index in case the coefficient domain of the expansions parent is only
            a module.  
        
        SEE:
            :meth:~`.fourier_expansion_homomorphism`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
            sage: em = ExpansionModule(Sequence([emps.one_element(), emps.one_element()]))
            sage: em._fourier_expansion_indices()
            [(1, 0)]
            sage: m = FreeModule(QQ, 3)
            sage: empsm = EquivariantMonoidPowerSeriesModule(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", m))
            sage: hv = EquivariantMonoidPowerSeries(empsm, {empsm.characters().one_element(): {1: m([1,0,0]), 2: m([0,0,1])}}, empsm.action().filter_all())
            sage: em = ExpansionModule(Sequence([hv]))
            sage: em._fourier_expansion_indices()
            [(0, (1, 0)), (1, (1, 0)), (2, (1, 0)), (0, (1, 1)), (1, (1, 1)), (2, (1, 1)), (0, (1, 2)), (1, (1, 2)), (2, (1, 2))]
            sage: mps = MonoidPowerSeriesRing(ZZ, NNMonoid())
            sage: em = ExpansionModule([mps.one_element()])
            sage: em._fourier_expansion_indices()
            [0]
            sage: mpsm = MonoidPowerSeriesModule(m, NNMonoid())
            sage: hv = MonoidPowerSeries(mpsm, {1: m([1,0,0]), 2: m([0,0,1])}, mpsm.monoid().filter_all())
            sage: em = ExpansionModule([hv])
            sage: em._fourier_expansion_indices()
            [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1), (0, 2), (1, 2), (2, 2)]
        """
        characters = self._non_zero_characters()
        
        if isinstance(self.__abstract_basis.universe().coefficient_domain(), Ring) :
            if characters is None :
                return list(self._bounding_precision())
            else :
                return [ (ch, k) for k in self._bounding_precision()
                                 for ch in characters ]
        else :
            if characters is None :
                return [ (i, k) for k in self._bounding_precision()
                            for i in range(self.__abstract_basis.universe().coefficient_domain().rank()) ]
            else :
                return [ (i,(ch, k)) for k in self._bounding_precision()
                            for ch in characters
                            for i in range(self.__abstract_basis.universe().coefficient_domain().rank()) ]

    @cached_method
    def fourier_expansion_homomorphism(self, precision = None) :
        r"""
        A morphism mapping elements of the underlying module to
        a module such that each component of the image corresponds to an
        Fourier coefficient of the Fourier expansion associated with this
        element.
        
        INPUT:
            - ``precision``  -- A fitler for the expansions monoid or action or ``None``
                                (default: ``None``); If not ``None`` the Fourier expansions
                                of the basis will be truncated to this precision.
        
        OUTPUT:
            A morphism from ``self`` to a free module over the expansions' base ring.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
            sage: em = ExpansionModule(Sequence([emps.one_element(), emps.one_element()]))
            sage: em.fourier_expansion_homomorphism()
            Free module morphism defined by the matrix
            [1]
            [1]
            Domain: Module of Fourier expansions in Ring of equivariant monoid power ...
            Codomain: Vector space of dimension 1 over Rational Field
            sage: m = FreeModule(QQ, 3)
            sage: empsm = EquivariantMonoidPowerSeriesModule(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", m))
            sage: hv = EquivariantMonoidPowerSeries(empsm, {empsm.characters().one_element(): {1: m([1,0,0]), 2: m([0,0,1])}}, empsm.action().filter_all())
            sage: em = ExpansionModule(Sequence([hv]))
            sage: em.fourier_expansion_homomorphism()
            Free module morphism defined by the matrix
            [0 0 0 1 0 0 0 0 1]
            Domain: Module of Fourier expansions in Module of equivariant monoid ...
            Codomain: Vector space of dimension 9 over Rational Field
        """
        keys = self._fourier_expansion_indices()
        
        codomain = FreeModule(self.base_ring(), len(keys))
        if isinstance(self.__abstract_basis.universe().coefficient_domain(), Ring) :
            basis_images = [ codomain([b[k] for k in keys])
                                            for b in self.__abstract_basis ]
        else :
            basis_images = [ codomain([b[k][i] for (i,k) in keys])
                                               for b in self.__abstract_basis ]
                            
        return self.Hom(codomain)(basis_images)
    
    @cached_method
    def _fourier_expansion_matrix_over_fraction_field(self) :
        r"""
        The matrix associated with the Fourier expansion homomorphism
        such that its base ring is a field.
        
        OUTPUT:
            A matrix over the fraction field of the expansions' base ring.
        
        SEE:
            :meth:~`.fourier_expansion_homomorphism`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: em = ExpansionModule(Sequence([emps.one_element(), emps.one_element()]))
            sage: em._fourier_expansion_matrix_over_fraction_field().base_ring()
            Rational Field
        """
        if self.base_ring().is_field() :
            return self.fourier_expansion_homomorphism().matrix()
        else :
            return self.fourier_expansion_homomorphism().matrix(). \
                    base_extend(self.base_ring().fraction_field())
                  
    @cached_method
    def pivot_elements(self) :
        r"""
        Determine a set of generators, which minimally span the image of the Fourier expansion
        homomorphism.
        
        SEE:
            :meth:~`.fourier_expansion_homomorphism`.
        
        OUTPUT:
            A list of indices.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: em = ExpansionModule(Sequence([emps.one_element(), emps.one_element()]))
            sage: em.pivot_elements()
            [0]
            sage: em = ExpansionModule(Sequence([emps.one_element().truncate(0), emps.one_element().truncate(0)]))
            sage: em.pivot_elements()
            []
            sage: em = ExpansionModule(Sequence([], universe = emps))
            sage: em.pivot_elements()
            []
            sage: m = FreeModule(QQ, 3)
            sage: empsm = EquivariantMonoidPowerSeriesModule(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", m))
            sage: hv = EquivariantMonoidPowerSeries(empsm, {empsm.characters().one_element(): {1: m([1,0,0]), 2: m([0,0,1])}}, empsm.action().filter_all())
            sage: em = ExpansionModule(Sequence([hv,hv]))
            sage: em.pivot_elements()
            [0]
        """
        expansion_matrix = self.fourier_expansion_homomorphism().matrix().transpose()

        if expansion_matrix.rank() == self.rank() :
            return range(self.rank())
        else :
            return list(expansion_matrix.pivots())
                  
    def _element_to_fourier_expansion_generator(self, e) :
        r"""
        Given a monoid power series `e` return a generator iterating over the components of
        the image of `e` under the Fourier expansion morphism.

        INTPUT:
            - `e` -- An element of the ``self``.
            
        OUTPUT:
            A generator over elements of the expansions' base ring.

        SEE::
            :meth:~`.fourier_expansion_homomorphism`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: em = ExpansionModule(Sequence([emps.one_element(), emps.one_element()]))
            sage: list(em._element_to_fourier_expansion_generator(em([3,0]).fourier_expansion()))
            [3]
            sage: m = FreeModule(QQ, 3)
            sage: empsm = EquivariantMonoidPowerSeriesModule(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", m))
            sage: hv = EquivariantMonoidPowerSeries(empsm, {empsm.characters().one_element(): {1: m([1,0,0]), 2: m([0,0,1])}}, empsm.action().filter_all())
            sage: em = ExpansionModule(Sequence([hv]))
            sage: list(em._element_to_fourier_expansion_generator(em([1]).fourier_expansion()))
            [0, 0, 0, 1, 0, 0, 0, 0, 1]
            sage: mps = MonoidPowerSeriesRing(ZZ, NNMonoid())
            sage: em = ExpansionModule([mps.one_element()])
            sage: list(em._element_to_fourier_expansion_generator(em([2]).fourier_expansion()))
            [2]
            sage: mpsm = MonoidPowerSeriesModule(m, NNMonoid())
            sage: hv = MonoidPowerSeries(mpsm, {1: m([1,0,0]), 2: m([0,0,1])}, mpsm.monoid().filter_all())
            sage: em = ExpansionModule([hv])
            sage: list(em._element_to_fourier_expansion_generator(em([1]).fourier_expansion()))
            [0, 0, 0, 1, 0, 0, 0, 0, 1]
        """
        keys = self._fourier_expansion_indices()
        
        if isinstance(self.__abstract_basis.universe().coefficient_domain(), Ring) :
            return (e[k] if k in e else 0 for k in keys) 
        else :
            return (e[k][i] if k in e else 0 for (i,k) in keys)

    def coordinates(self, x, in_base_ring = True, force_ambigous = False) :
        r"""
        The coordinates in ``self`` of an element either of the following:
            - An element of a submodule.
            - An expansion in the parent of the basis' expansions.
        
        INPUT:
            - `x` -- Either of the types listed above.
            - ``in_base_ring`` -- A boolean (default: ``True``); If ``True``
                                  enforce the result to be definied over the
                                  base ring of ``self``.
            - ``force_ambigous`` -- A boolean (default: ``False``); If ``True``
                                    also return the solutions that are not unique.
        
        OUTPUT:
            A list of elements in the base ring or and extension of it.
        
        NOTE:
            If the Fourier expansion of `x` lakes sufficient precision the
            expansions associated to ``self`` will be truncated.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: em = ExpansionModule(Sequence([emps.one_element()]))
            sage: em.coordinates(2 * emps.one_element())
            [2]
            sage: m = FreeModule(QQ, 3)
            sage: empsm = EquivariantMonoidPowerSeriesModule(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", m))
            sage: hv = EquivariantMonoidPowerSeries(empsm, {empsm.characters().one_element(): {1: m([1,0,0]), 2: m([0,0,1])}}, empsm.action().filter_all())
            sage: em = ExpansionModule(Sequence([hv]))
            sage: em.coordinates(2 * hv)
            [2]
            sage: mps = MonoidPowerSeriesRing(ZZ, NNMonoid())
            sage: em = ExpansionModule([mps.one_element()])
            sage: em.coordinates(3 * mps.one_element())
            [3]
            sage: mpsm = MonoidPowerSeriesModule(m, NNMonoid())
            sage: hv = MonoidPowerSeries(mpsm, {1: m([1,0,0]), 2: m([0,0,1])}, mpsm.monoid().filter_all())
            sage: em = ExpansionModule([hv])
            sage: em.coordinates(2 * hv)
            [2]
            sage: em = ExpansionModule(Sequence([emps.one_element()]))
            sage: h = EquivariantMonoidPowerSeries(emps, {}, emps.action().filter_all())
            sage: em.coordinates(h)
            [0]
            sage: em.coordinates(1 / 2 * emps.one_element())
            Traceback (most recent call last):
            ...
            ArithmeticError: Equivariant monoid power series in Ring of equivariant monoid power series over NN is not contained in this space.
            sage: em.coordinates(1 / 2 * emps.one_element(), in_base_ring = False )
            [1/2]
            sage: K.<rho> = CyclotomicField(6)
            sage: em.coordinates(rho * emps.one_element(), in_base_ring = False)
            [zeta6]
            sage: h = EquivariantMonoidPowerSeries(emps, {emps.characters().one_element(): {1: 2, 2: 3}}, emps.action().filter_all())
            sage: em.coordinates(h)
            Traceback (most recent call last):
            ...
            ArithmeticError: Equivariant monoid power series in Ring of equivariant monoid power series over NN is not contained in this space.
            sage: em.coordinates(hv)
            Traceback (most recent call last):
            ...
            ArithmeticError: No coordinates for Monoid power series in Module of monoid power series over NN with action.
            sage: em = ExpansionModule(Sequence([emps.one_element(), emps.one_element()]))
            sage: em.coordinates(2 * emps.one_element())
            Traceback (most recent call last):
            ...
            ValueError: Not unambigous coordinates available in this submodule.
            sage: h = EquivariantMonoidPowerSeries(emps, {}, emps.action().filter_all())
            sage: em = ExpansionModule(Sequence([h]))
            sage: em.coordinates(emps.one_element())
            Traceback (most recent call last):
            ...
            ValueError: Not unambigous coordinates available in this submodule.
        """
        if isinstance(x, Element) :
            P = x.parent()
            if P is self :
                if self.ambient_module() is self :
                    return x.list()
                else :
                    return list(self.coordinate_vector(x))

                        
            if isinstance(P, (EquivariantMonoidPowerSeriesAmbient_abstract,
                                MonoidPowerSeriesAmbient_abstract)) :
                             
                if not force_ambigous and \
                   not self._check_precision() :
                    raise ValueError( "Not unambigous coordinates available in this submodule." )

                if P != self.__abstract_basis.universe() :                    
                    try :
                        Pnew = pushout(self.__abstract_basis.universe(), P)
                    except TypeError :
                        raise ArithmeticError( "No coordinates for %s." % (x,) )
                    
                    if isinstance(Pnew, EquivariantMonoidPowerSeriesAmbient_abstract) \
                       and Pnew.representation().base_ring() == Pnew.representation().codomain() : 
                        A = Pnew.base_ring()
                    elif isinstance(Pnew, MonoidPowerSeriesAmbient_abstract) \
                       and isinstance(Pnew.coefficient_domain(), Ring) :
                        A = Pnew.base_ring()
                    else :
                        A = Pnew.base_ring().base_ring()

                    try :
                        x = Pnew(x)
                    except TypeError :
                        raise ArithmeticError( "No coordinates for %s." % (x,) )
                else :
                    A = self.base_ring()

                # check those components of x, which have to be zero
                if not self._non_zero_characters() is None :
                    x._cleanup_coefficients()
                    if set( x.non_zero_components() ) - set( self._non_zero_characters() ) != set() :
                        raise ArithmeticError( "%s is not contained in this space." % (x,) )
                                                
                if in_base_ring :
                    fe_matrix = self.fourier_expansion_homomorphism().matrix().transpose()
                    x_fe_vector = matrix( A, fe_matrix.nrows(),
                                   list(self._element_to_fourier_expansion_generator(x)) )
                else :
                    fe_matrix = self._fourier_expansion_matrix_over_fraction_field().transpose()
                    x_fe_vector = matrix( A, fe_matrix.nrows(),
                                   list(self._element_to_fourier_expansion_generator(x)) )
                
                if self._non_zero_characters() is None :
                    if isinstance(self.__abstract_basis.universe().coefficient_domain(), Ring) :
                        valid_indices = [ i for (i,k) in enumerate(self._fourier_expansion_indices())
                                          if k in x.precision() ]
                    else :
                        valid_indices = [ i for (i,(_,k)) in enumerate(self._fourier_expansion_indices())
                                          if k in x.precision() ]
                else :
                    if isinstance(self.__abstract_basis.universe().coefficient_domain(), Ring) :
                        valid_indices = [ i for (i,(_,k)) in enumerate(self._fourier_expansion_indices())
                                          if k in x.precision() ]
                    else :
                        valid_indices = [ i for (i,(_,(_,k))) in enumerate(self._fourier_expansion_indices())
                                          if k in x.precision() ]
                    
                fe_matrix = fe_matrix.matrix_from_rows(valid_indices)
                x_fe_vector = x_fe_vector.matrix_from_rows(valid_indices)
                fe_matrix = matrix(A, fe_matrix)
                        
                ## TODO: use linbox
                try :
                    ## TODO: We deactivate the magma interface as it is almost never used and
                    ##       it is in bad shape
                    if True or (A is not ZZ and A is not QQ) :
                        raise TypeError
                    
                    magma_fe_matrix = magma(fe_matrix.transpose())
                    
                    if not force_ambigous and \
                       len(valid_indices) != self.fourier_expansion_homomorphism().matrix().ncols() and \
                       magma_fe_matrix.Rank() != self.rank() :
                        raise ValueError( "No unambigous coordinates available." )
                
                    coords = magma(fe_matrix.transpose()).Solution(magma(x_fe_vector.transpose())).sage()
                    coords = coords.row(0)
                except TypeError, msg :
                    if "Runtime error in 'Solution': No solution exists" in msg :
                        raise ArithmeticError( "%s is not contained in this space." % (x,) )

                    if not force_ambigous and \
                       len(valid_indices) != self.fourier_expansion_homomorphism().matrix().ncols() and \
                       fe_matrix.rank() != self.rank() :
                        raise ValueError( "No unambigous coordinates available." )
                    
                    try :
                        coords = fe_matrix.solve_right(x_fe_vector)
                    except ValueError, msg :
                        raise ArithmeticError( "%s is not contained in this space, %s" % (x, msg) )
                    coords = coords.column(0)
                    
                    if self.precision() != self._bounding_precision() and \
                       not x.precision() < self._bounding_precision() :
                        if not self.change_ring(A)(list(coords)).fourier_expansion() == x :
                            raise ArithmeticError( "%s is not contained in this space." % (x,) )   
                
                if in_base_ring :
                    try :
                        return [self.base_ring()(c) for c in coords]
                    except TypeError :
                        raise ArithmeticError( "%s is not contained in this space." % (x,) )
                else :
                    return list(coords)
            
            #! elif isinstance(P, (EquivariantMonoidPowerSeriesAmbient_abstract,
            #                      MonoidPowerSeriesAmbient_abstract)) :
        #! if isinstance(x, (tuple, Element)) :

        raise ArithmeticError( "No coordinates for %s." % (x,) )
        
    def _sparse_module(self):
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: em = ExpansionModule(Sequence([emps.one_element()]))
            sage: em._sparse_module()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def _dense_module(self):
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: em = ExpansionModule(Sequence([emps.one_element()]))
            sage: em._dense_module() is em
            True
        """
        return self

    def _repr_(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: ExpansionModule(Sequence([emps.one_element()]))
            Module of Fourier expansions in Ring of equivariant monoid power series over NN
        """
        return "Module of Fourier expansions in %s" % (self.__abstract_basis.universe(),)            

    def _latex_(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: latex( ExpansionModule(Sequence([emps.one_element()])) )
            \text{Module of Fourier expansions in }\text{Ring of equivariant monoid power series over }\Bold{N}
        """
        return r"\text{Module of Fourier expansions in }%s" % (latex(self.__abstract_basis.universe()),)


class ExpansionModule_generic ( ExpansionModule_abstract, FreeModule_generic ) :
    r"""
    A generic module of abstract elements with Fourier expansions attached to.
    The base ring has to be an integral domain.
    """
    
    def __init__(self, basis, degree, **kwds) :
        r"""
        INPUT:
            - ``basis``  -- A sequence of (equivariant) monoid power series.
            - ``degree`` -- An integer; The ambient's module dimension.
            - ``kwds``   -- A keyword dictionary that will be forwarded to
                            initialization of the underlying free module.

        NOTE:
            The base ring of the expansions' parent must be an integral
            domain.

        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.expansion_module import ExpansionModule_generic, ExpansionModuleVector_generic
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: em = ExpansionModule_generic(Sequence([emps.one_element()]), 2)
            sage: em = ExpansionModule_generic(Sequence([emps.one_element()]), 2, _element_class = ExpansionModuleVector_generic )
        """
        if not hasattr(self, '_element_class') :
            if '_element_class' in kwds :
                self._element_class = kwds['_element_class']
            else :
                self._element_class = ExpansionModuleVector_generic

        ExpansionModule_abstract.__init__(self, basis)
        FreeModule_generic.__init__(self, basis.universe().base_ring(), len(basis), degree, sparse = False)

    def gen(self, i) :
        r"""
        The `i`-th generator of the module.
        
        INPUT:
            - `i` -- An integer.
            
        OUTPUT:
            An element of ``self``.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.expansion_module import ExpansionModule_generic
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: em = ExpansionModule_generic(Sequence([emps.one_element()]), 1)
            sage: em.gen(0)
            Traceback (most recent call last):
            ...
            NotImplementedError
            """
        raise NotImplementedError

    def basis(self) :
        r"""
        A basis of ``self``.
        
        OUTPUT:
            A list of elements of ``self``.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.expansion_module import ExpansionModule_generic
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: em = ExpansionModule_generic(Sequence([emps.one_element()]), 1)
            sage: em.basis()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError
    
    def change_ring(self, R):
        r"""
        Return the ambient expansion module over `R` with the same basis as ``self``.
        
        INPUT:
            - `R` -- A ring.
        
        OUTPUT:
            An instance of :class:~`.ExpansionModule_generic`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.expansion_module import ExpansionModule_generic
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: em = ExpansionModule_generic(Sequence([emps.one_element()]), 2)
            sage: em.change_ring(ZZ) is em
            True
            sage: em.change_ring(QQ).base_ring()
            Rational Field
        """
        if self.base_ring() == R :
            return self
        
        R = pushout(self._abstract_basis().universe(), R)

        return ExpansionModule_generic(Sequence(self._abstract_basis(), universe = R), self.degree(), _element_class = self._element_class)
        
#===============================================================================
# _fourier_expansion_kernel
#
# This will be used in ExpansionModule_ambient_pid and
# ExpansionModule_submodule_pid 
#===============================================================================

def _fourier_expansion_kernel(self) :
    r"""
    The kernel of the Fourier expansion morphism.
    
    OUTPUT:
        A module over the basis ring.
    
    SEE:
        :meth:~`.ExpansionModule_abstract.fourier_expansion_homomorphism`.
    
    TESTS::
        sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
        sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
        sage: em = ExpansionModule(Sequence([emps.one_element(), emps.one_element(), emps.one_element()]))
        sage: em.fourier_expansion_kernel()
        Free module of degree 3 and rank 2 over Integer Ring
        Echelon basis matrix:
        [ 1  0 -1]
        [ 0  1 -1]
        sage: em.span([em([1,0,0]), em([1,2,0])]).fourier_expansion_kernel()
        Free module of degree 2 and rank 1 over Integer Ring
        Echelon basis matrix:
        [ 3 -1]
    """
    return self.fourier_expansion_homomorphism().matrix().left_kernel()

#===============================================================================
# _span
#
# This will be used in ExpansionModule_ambient_pid and
# ExpansionModule_submodule_pid 
#===============================================================================

def _span( self, gens, base_ring = None, check = True, already_echelonized = False ) :
    r"""
    The expansion submodule spanned by ``gens``.
    
    INPUT:
        - ``gens``                -- A list, tuple or sequence of module elements.
        - ``base_ring``           -- A ring or ``None`` (default: ``None``); If ``None``
                                     the base ring of ``self`` will be used.
        - ``check``               -- A boolean (default: ``True``); If ``True`` check
                                     whether the generators are appropriately coerced.
        - ``already_echelonized`` -- A boolean (default: ``False``); If ``True``
                                     the generators are already in echelon form
                                     with respect to the ambient's basis.
    
    OUTPUT:
        An instance of :class:~`.ExpansionModule_submodule_pid`.
    
    TESTS::
        sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
        sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
        sage: em = ExpansionModule(Sequence([emps.one_element(), emps.one_element(), emps.one_element()]))
        sage: em.span([em([1,0,0]), em([1,2,0])])
        Module of Fourier expansions in Ring of equivariant monoid power series over NN
        sage: em.span([em([1,0,0]), em([1,2,0])]).span([em([1,0,0])])
        Module of Fourier expansions in Ring of equivariant monoid power series over NN
        """
    if is_FreeModule(gens):
        gens = gens.gens()
    if not isinstance(gens, (list, tuple, Sequence)):
        raise TypeError, "Argument gens (= %s) must be a list, tuple, or sequence." % (gens,)

    if base_ring is None or base_ring == self.base_ring() :
        gens = Sequence(gens, check = check, universe = self.ambient_module())
        
        return ExpansionModule_submodule_pid(self.ambient_module(), gens)
    else:
        try:
            M = self.ambient_module().change_ring(base_ring)
        except TypeError:
            raise ValueError, "Argument base_ring (= %s) is not compatible " % (base_ring,) + \
                "with the base field (= %s)." % (self.base_field(),)
        try: 
            return M.span(gens)
        except TypeError:
            raise ValueError, "Argument gens (= %s) is not compatible " % (gens,) + \
                "with base_ring (= %s)." % (base_ring,)

#===============================================================================
# ExpansionModule_ambient_pid
#===============================================================================

class ExpansionModule_ambient_pid ( ExpansionModule_abstract, FreeModule_ambient_pid ) :
    r"""
    An ambient module of expansions over a principal ideal domain.
    """
    
    def __init__(self, basis, _element_class = None, **kwds) :
        r"""
        INPUT:
            - ``basis``         -- A list or sequence of (equivariant) monoid power series.
            - ``element_class`` -- A type or ``None`` (default: ``None``); The element class
                                   attached to ``self``. If ``None`` :class:~`.ExpansionModuleVector_generic`
                                   will be used.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.expansion_module import ExpansionModule_ambient_pid, ExpansionModuleVector_generic
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: em = ExpansionModule_ambient_pid(Sequence([emps.one_element(), emps.one_element(), emps.one_element()]))
            sage: em = ExpansionModule_ambient_pid(Sequence([emps.one_element(), emps.one_element(), emps.one_element()]), _element_class = ExpansionModuleVector_generic)
        """
        if not hasattr(self, '_element_class') :
            if not _element_class is None :
                self._element_class = _element_class
            else :
                self._element_class = ExpansionModuleVector_generic

        if isinstance(basis.universe(), EquivariantMonoidPowerSeriesAmbient_abstract) \
           and basis.universe().representation().base_ring() == basis.universe().representation().codomain() : 
            base_ring = basis.universe().base_ring()
        elif isinstance(basis.universe(), MonoidPowerSeriesAmbient_abstract) \
           and isinstance(basis.universe().coefficient_domain(), Ring) :
            base_ring = basis.universe().base_ring()
        else :
            base_ring = basis.universe().base_ring().base_ring()
        
        ExpansionModule_abstract.__init__(self, basis)
        FreeModule_ambient_pid.__init__(self, base_ring, len(basis))
        
    global _fourier_expansion_kernel, _span
    fourier_expansion_kernel = _fourier_expansion_kernel
    span = _span
    
    def gen(self, i) :
        r"""
        The `i`-th generator of the module.
        
        INPUT:
            - `i` -- An integer.
            
        OUTPUT:
            An element of ``self``.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.expansion_module import ExpansionModule_ambient_pid
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: em = ExpansionModule_ambient_pid(Sequence([emps.one_element(), emps.one_element(), emps.one_element()]))
            sage: em.gen(0)
            (1, 0, 0) 
        """
        return self._element_class(self, FreeModule_ambient_pid.gen(self, i))

    def basis(self) :
        r"""
        A basis of ``self``.
        
        OUTPUT:
            A list of elements of ``self``.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.expansion_module import ExpansionModule_ambient_pid
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: em = ExpansionModule_ambient_pid(Sequence([emps.one_element(), emps.one_element(), emps.one_element()]))
            sage: em.basis()
            [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
        """
        return map(lambda b: self._element_class(self, b), FreeModule_ambient_pid.basis(self))    
        
    def change_ring(self, R):
        r"""
        Return the ambient expansion module over `R` with the same basis as ``self``.
        
        INPUT:
            - `R` -- A ring.
        
        OUTPUT:
            An instance of :class:~`.ExpansionModule_ambient_pid`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.expansion_module import ExpansionModule_ambient_pid
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: em = ExpansionModule_ambient_pid(Sequence([emps.one_element(), emps.one_element(), emps.one_element()]))
            sage: em.change_ring(ZZ) is em
            True
            sage: emq = em.change_ring(QQ)
            sage: emq.base_ring()
            Rational Field
            sage: emq.ambient_module() is emq
            True
        """
        if self.base_ring() == R :
            return self
        
        R = pushout(self._abstract_basis().universe(), R)

        return ExpansionModule_ambient_pid(Sequence(self._abstract_basis(), universe = R), _element_class = self._element_class)

    def ambient_module(self) :
        r"""
        Return the ambient module of ``self``.
        
        OUTPUT:
            An instance of :class:~`.ExpansionModule_ambient_pid`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.expansion_module import ExpansionModule_ambient_pid
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: em = ExpansionModule_ambient_pid(Sequence([emps.one_element(), emps.one_element(), emps.one_element()]))
            sage: em.ambient_module() is em
            True
        """
        return self

#===============================================================================
# ExpansionModule_submodule_pid
#===============================================================================

class ExpansionModule_submodule_pid ( ExpansionModule_abstract, FreeModule_submodule_pid ) :
    r"""
    A submodule of another module of expansions over a principal ideal domain.
    """

    def __init__(self, ambient, gens, _element_class = None, **kwds) :
        r"""
        INPUT:
            - ``ambient``       -- An instance of :class:~`.ExpansionModule_ambient_pid` or :class:~`.ExpansionModule_submodule_pid`.
            - ``gens``          -- A list, tuple or sequence of elements of ``ambient``.
            - ``element_class`` -- A type or ``None`` (default: ``None``); The element class
                                   attached to ``self``. If ``None`` :class:~`.ExpansionModuleVector_generic`
                                   will be used.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.expansion_module import ExpansionModule_ambient_pid, ExpansionModule_submodule_pid, ExpansionModuleVector_generic
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: em = ExpansionModule_ambient_pid(Sequence([emps.one_element(), emps.one_element(), emps.one_element()]))
            sage: ems = ExpansionModule_submodule_pid(em, [em([1,0,0]), em([1,2,0])])
            sage: ems = ExpansionModule_submodule_pid(em, [em([1,0,0]), em([1,2,0])], _element_class = ExpansionModuleVector_generic)
        """
        if not hasattr(self, '_element_class') :
            if not _element_class is None :
                self._element_class = _element_class
            else :
                self._element_class = ExpansionModuleVector_generic

        if ambient.precision().is_infinite() :
            ExpansionModule_abstract.__init__(self, Sequence( map( lambda g: g.fourier_expansion(), gens ),
                                                              universe = ambient._abstract_basis().universe() ))
        else :
            ExpansionModule_abstract.__init__(self, Sequence( map( lambda g: LazyFourierExpansionEvaluation( ambient._abstract_basis().universe(), g,
                                                                                                             ambient.precision() ),
                                                                   gens ),
                                                              universe = ambient._abstract_basis().universe() ))
        FreeModule_submodule_pid.__init__(self, ambient, gens)

    global _fourier_expansion_kernel, _span
    fourier_expansion_kernel = _fourier_expansion_kernel
    span = _span

    def gen(self, i) :
        r"""
        The `i`-th generator of the module.
        
        INPUT:
            - `i` -- An integer.
            
        OUTPUT:
            An element of ``self``.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.expansion_module import ExpansionModule_ambient_pid, ExpansionModule_submodule_pid
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: em = ExpansionModule_ambient_pid(Sequence([emps.one_element(), emps.one_element(), emps.one_element()]))
            sage: ems = ExpansionModule_submodule_pid(em, [em([1,0,0]), em([1,2,0])])
            sage: ems.gen(0)
            (1, 0, 0)
        """
        return self._element_class(self, super(ExpansionModule_submodule_pid, self).gen(i).list())

    # TODO: The module should be hidden so that we can adopt the category framework
    # def basis(self) :
    #     r"""
    #     A basis of ``self``.
        
    #     OUTPUT:
    #         A list of elements of ``self``.
        
    #     TESTS::
    #         sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
    #         sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
    #         sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
    #         sage: from psage.modform.fourier_expansion_framework.gradedexpansions.expansion_module import ExpansionModule_ambient_pid, ExpansionModule_submodule_pid
    #         sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
    #         sage: em = ExpansionModule_ambient_pid(Sequence([emps.one_element(), emps.one_element(), emps.one_element()]))
    #         sage: ems = ExpansionModule_submodule_pid(em, [em([1,0,0]), em([1,2,0])])
    #         sage: ems.basis()
    #         [(1, 0, 0), (0, 2, 0)]
    #     """
    #     return [self._element_class(self, b.list()) for b in super().basis()]
    
    def change_ring(self, R):
        r"""
        Return the ambient expansion module over `R` with the same basis as ``self``.
        
        INPUT:
            - `R` -- A ring.
        
        OUTPUT:
            An instance of :class:~`.ExpansionModule_ambient_pid`.
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.expansion_module import ExpansionModule_ambient_pid, ExpansionModule_submodule_pid
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: em = ExpansionModule_ambient_pid(Sequence([emps.one_element(), emps.one_element(), emps.one_element()]))
            sage: ems = ExpansionModule_submodule_pid(em, [em([1,0,0]), em([1,2,0])])
            sage: emc = ems.change_ring(ZZ)
            sage: emc.ambient_module() is emc
            True
            sage: emc = ems.change_ring(QQ)
            sage: emc.ambient_module() is emc
            True
            sage: emc.base_ring()
            Rational Field
        """
        if self.base_ring() == R :
            return ExpansionModule_ambient_pid(self._abstract_basis(), _element_class = self._element_class)
        
        R = pushout(self._abstract_basis().universe(), R)

        return ExpansionModule_ambient_pid(Sequence(self._abstract_basis(), universe = R), _element_class = self._element_class)

#===============================================================================
# ExpansionModuleVector_generic
#===============================================================================

class ExpansionModuleVector_generic ( FreeModuleElement_generic_dense, FourierExpansionWrapper ) :
    r"""
    A generic implementation of an element in a module of expansions.
    """
    
    def __init__(self, parent, x, coerce = True, copy = True) :
        r"""
        INPUT:
            - ``parent``  -- An instance of :class:~`.ExpansionModule_abstract`.
            - `x`         -- A list or tuple of integers or an element that admits coordinates
                             in ``parent``.
            - ``coerce``  -- A boolean (default: ``True``); If ``True`` coerce coordinates into the base ring.
            - ``copy``    -- A boolean (default: ``True``); If ``True`` store a copy of the coordinates.

        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.expansion_module import ExpansionModuleVector_generic
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: em = ExpansionModule(Sequence([emps.one_element(), emps.one_element(), emps.one_element()]))
            sage: ExpansionModuleVector_generic(em, [1,0,2])
            (1, 0, 2)
        """
        if not isinstance(x, (list, tuple)) and x != 0 :
            x = parent.ambient_module().coordinates(x)
            coerce = False
            copy = False
            
        FreeModuleElement_generic_dense.__init__(self, parent, x, coerce, copy)
    
    def _add_(left, right) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: em = ExpansionModule(Sequence([emps.one_element(), emps.one_element(), emps.one_element()]))
            sage: em.0 + em.1
            (1, 1, 0)
        """
        return left.parent()._element_class(left.parent(), FreeModuleElement_generic_dense._add_(left, right))
    
    def __copy__(self) :
        r"""
        Return a copy of ``self``.
        
        OUTPUT:
            An instance of :class:~`.ExpansionModuleVector_generic`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.expansion_module import ExpansionModuleVector_generic
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: em = ExpansionModule(Sequence([emps.one_element(), emps.one_element(), emps.one_element()]))
            sage: copy(ExpansionModuleVector_generic(em, [1,0,2])) == ExpansionModuleVector_generic(em, [1,0,2])
            True
        """
        return ExpansionModuleVector_generic( self.parent(), self.list(), coerce = False, copy = True)
