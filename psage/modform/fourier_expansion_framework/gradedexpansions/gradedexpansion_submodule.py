"""
Finite dimensional submodules of a ring or module of graded expansions. 

AUTHOR :
    -- Martin Raum (2009 - 07 - 27) Initial version
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

from psage.modform.fourier_expansion_framework.gradedexpansions.expansion_lazy_evaluation import LazyFourierExpansionEvaluation
from psage.modform.fourier_expansion_framework.gradedexpansions.expansion_module import ExpansionModule_abstract
from psage.modform.fourier_expansion_framework.gradedexpansions.expansion_module import ExpansionModule_generic, ExpansionModule_ambient_pid, \
      ExpansionModule_submodule_pid, ExpansionModuleVector_generic
from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ambient import GradedExpansionAmbient_abstract
from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ambient import MonoidPowerSeriesAmbient_abstract
from sage.categories.pushout import pushout
from sage.interfaces.magma import magma
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_method
from sage.misc.flatten import flatten
from sage.misc.latex import latex
from sage.misc.misc import mul
from sage.misc.misc import union
from sage.modules.free_module import FreeModule, FreeModule_generic, \
      FreeModule_ambient_pid, FreeModule_submodule_pid 
from sage.modules.free_module import is_FreeModule
from sage.modules.free_module_element import FreeModuleElement_generic_dense
from sage.modules.free_module_element import vector
from sage.rings.arith import random_prime
from sage.rings.integer_ring import ZZ
from sage.rings.number_field.order import Order as NumberFieldOrder
from sage.rings.padics.factory import Qp
from sage.rings.principal_ideal_domain import PrincipalIdealDomain
from sage.rings.rational_field import QQ
from sage.rings.ring import Ring
from sage.structure.element import Element
from sage.structure.sequence import Sequence

#===============================================================================
# GradedExpansionSubmodule
#===============================================================================

def GradedExpansionSubmodule(arg1, arg2) :
    """
    A submodule of either a graded ring, module or submodule.
    
    INPUT (first possibility):
        - ``arg1`` -- An instance of :class:~`fourier_expansion_framework.gradedexpansions.gradedexpansion_ambient.GradedExpansionAmbient_abstract`.
        - ``arg2`` -- A tuple, list or sequence of elements of ``arg1``.
    INPUT (second possibility):
        - ``arg1`` -- An instance of :class:~`.GradedExpansionSubmodule_ambient_pid` or :class:~`.GradedExpansionSubmodule_submodule_pid`.
        - ``arg2`` -- A tuple, list or sequence of elements of ``arg1``.

    NOTE:
        The base ring of the graded expansion ambient must be an integral domain.

    OUTPUT:
        An instance of :class:~`.GradedExpansionSubmodule_abstract`.
    
    TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_module import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_submodule import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: ger = GradedExpansionRing_class(None, Sequence([MonoidPowerSeries(mps, {1 : 4, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {1 : 1, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {2 : 1, 3: 6, 4 : 9}, mps.monoid().filter_all())]), PolynomialRing(QQ, ['a', 'b', 'c']).ideal(0), DegreeGrading((1,2,3)))
            sage: sm = GradedExpansionSubmodule(ger, [ger.0, ger.1])
            sage: sm = GradedExpansionSubmodule(ger, (ger.0, ger.1))
            sage: sm = GradedExpansionSubmodule(ger, Sequence([ger.0, ger.1]))
            sage: sm2 = GradedExpansionSubmodule(sm, [sm.0])
            sage: m = FreeModule(QQ, 3)
            sage: mpsm = MonoidPowerSeriesModule(m, NNMonoid(False))
            sage: mps = mpsm.base_ring()
            sage: ger = GradedExpansionModule_class(None, Sequence([MonoidPowerSeries(mpsm, {1 : m([1,2,3]), 2 : m([3,-3,2])}, mpsm.monoid().filter(4)), MonoidPowerSeries(mpsm, {1 : m([2,-1,-1]), 2 : m([1,0,0])}, mpsm.monoid().filter(4))]), PolynomialRing(ZZ, ['a', 'b']).ideal(0), DegreeGrading((1,2)))
            sage: sm = GradedExpansionSubmodule(ger, [ger.0])
    """
    if isinstance(arg1, GradedExpansionAmbient_abstract) :
        base_ring = arg1.relations().base_ring()
        
        if base_ring.is_field() or \
           isinstance(base_ring, PrincipalIdealDomain) or \
           isinstance(base_ring, NumberFieldOrder) \
            and base_ring.is_maximal() and base_ring.class_number() == 1 :
                return GradedExpansionSubmodule_ambient_pid(arg1, arg2)
        else :
                return GradedExpansionSubmodule_generic(arg1, arg2, len(arg2))
    elif isinstance(arg1, GradedExpansionSubmodule_ambient_pid) \
       or isinstance(arg1, GradedExpansionSubmodule_submodule_pid) :
        return GradedExpansionSubmodule_submodule_pid(arg1, arg2)

    raise ValueError( "Cannot construct a new subspace from %s, %s" % (arg1, arg2) )

#===============================================================================
# GradedExpansionSubmodule_abstract
#===============================================================================

class GradedExpansionSubmodule_abstract ( ExpansionModule_abstract ) :
    """
    Abstract implementation of a finite dimensional module of graded expansions
    within an ambient.
    
    SEE:
        :class:~`fourier_expansion_framework.gradedexpansions.gradedexpansion_ambient.GradedExpansion_ambient`.
    """
    
    def __init__(self, graded_ambient, basis, degree, **kwds) :
        """
        INPUT:
            - ``graded_ambient`` -- An instance of :class:~`fourier_expansion_framework.gradedexpansions.gradedexpansion_ambient.GradedExpansionAmbient_abstract`.
            - ``basis``          -- A tuple, list or sequence of elements of the graded ambient.
            - ``degree``         -- An integer; The degree of the module within its ambient module.
        
        NOTE:
            The base ring of the graded expansion ambient must be an integral domain.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_module import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_submodule import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: ger = GradedExpansionRing_class(None, Sequence([MonoidPowerSeries(mps, {1 : 4, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {1 : 1, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {2 : 1, 3: 6, 4 : 9}, mps.monoid().filter_all())]), PolynomialRing(QQ, ['a', 'b', 'c']).ideal(0), DegreeGrading((1,2,3)))
            sage: sm = GradedExpansionSubmodule_abstract(ger, Sequence([ger.0, ger.1]), 3)
            sage: sm = GradedExpansionSubmodule_abstract(ger, Sequence([ger.0, ger.1]), 2, no_expansion_init = True)
            sage: m = FreeModule(QQ, 3)
            sage: mpsm = MonoidPowerSeriesModule(m, NNMonoid(False))
            sage: mps = mpsm.base_ring()
            sage: ger = GradedExpansionModule_class(None, Sequence([MonoidPowerSeries(mpsm, {1 : m([1,2,3]), 2 : m([3,-3,2])}, mpsm.monoid().filter(4)), MonoidPowerSeries(mpsm, {1 : m([2,-1,-1]), 2 : m([1,0,0])}, mpsm.monoid().filter(4))]), PolynomialRing(ZZ, ['a', 'b']).ideal(0), DegreeGrading((1,2)))
            sage: sm = GradedExpansionSubmodule_abstract(ger, [ger.0], 1)
        """
        self.__graded_ambient = graded_ambient
        self.__basis_in_graded_ambient = basis
        
        if not "no_expansion_init" in kwds or not kwds["no_expansion_init"] :
            if graded_ambient.fourier_expansion_precision().is_infinite() or isinstance(self.__graded_ambient.fourier_ring(), MonoidPowerSeriesAmbient_abstract) :
                ExpansionModule_abstract.__init__(self, Sequence( map( lambda b: b.fourier_expansion(), basis ),
                                                                  universe = graded_ambient.fourier_ring() ))
            else :
                ExpansionModule_abstract.__init__(self, Sequence( map( lambda b: LazyFourierExpansionEvaluation( graded_ambient.fourier_ring(), b,
                                                                                                       graded_ambient.fourier_expansion_precision() ),
                                                                       basis ),
                                                                  universe = graded_ambient.fourier_ring(), check = False ))
        
    def graded_ambient(self) :
        """
        The graded ambientm, namely, the graded ring or module
        ``self`` is a submodule of.
        
        OUTPUT:
            An instance of :class:~`fourier_expansion_framework.gradedexpansions.gradedexpansion_ambient.GradedExpansionAmbient_abstract`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_module import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_submodule import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: ger = GradedExpansionRing_class(None, Sequence([MonoidPowerSeries(mps, {1 : 4, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {1 : 1, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {2 : 1, 3: 6, 4 : 9}, mps.monoid().filter_all())]), PolynomialRing(QQ, ['a', 'b', 'c']).ideal(0), DegreeGrading((1,2,3)))
            sage: sm = GradedExpansionSubmodule_abstract(ger, Sequence([ger.0, ger.1]), 3)
            sage: sm.graded_ambient() is ger
            True
        """
        return self.__graded_ambient

    def _basis_in_graded_ambient(self) :
        """
        A basis for ``self`` in terms of elements of the graded ambient.
        
        OUTPUT:
            A sequence of elements of the graded ambient.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_module import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_submodule import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: ger = GradedExpansionRing_class(None, Sequence([MonoidPowerSeries(mps, {1 : 4, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {1 : 1, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {2 : 1, 3: 6, 4 : 9}, mps.monoid().filter_all())]), PolynomialRing(QQ, ['a', 'b', 'c']).ideal(0), DegreeGrading((1,2,3)))
            sage: sm = GradedExpansionSubmodule_abstract(ger, Sequence([ger.0, ger.1]), 3)
            sage: sm._basis_in_graded_ambient() == Sequence([ger.0, ger.1])
            True
            """
        return self.__basis_in_graded_ambient

    @cached_method
    def _reduced_basis_polynomials(self, coerce_to = None) :
        """
        A list of reduced polynomials associated to the basis of ``self``
        within the graded ambient. If coerce_to is not None, these elements
        will first be coerced into ``coerce_to`` and the resulting polynomials
        will be returned.
        
        INPUT:
            - ``coerce_to`` -- A graded ambient or ``None`` (default: ``None``);
                               See the discription above.
        
        OUTPUT:
            A Sequence of polynomials.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_module import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_submodule import *
            sage: K.<rho> = CyclotomicField(6)
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: P.<a, b, c> = PolynomialRing(K)
            sage: ger = GradedExpansionRing_class(None, Sequence([MonoidPowerSeries(mps, {1 : 4, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {1 : 1, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {2 : 1, 3: 6, 4 : 9}, mps.monoid().filter_all())]), PolynomialRing(QQ, ['a', 'b', 'c']).ideal(0), DegreeGrading((1,2,3)))
            sage: ger2 = GradedExpansionRing_class(None, Sequence([MonoidPowerSeries(mps, {1 : 4, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {1 : 1, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {2 : 1, 3: 6, 4 : 9}, mps.monoid().filter_all())]), P.ideal(a - b), DegreeGrading((1,2,3)))
            sage: sm = GradedExpansionSubmodule_abstract(ger, Sequence([ger.0, ger.1]), 3)
            sage: sm._reduced_basis_polynomials()
            [a, b]
            
        We need to introduce a workaround for coercions as long as graded expansion ambients do not
        support coercion from one to another. Since the doctests to not allow for class definitions,
        we hack an instance of an arbitrary Python class.
        ::
            sage: from htmllib import HTMLParser
            sage: coerce_workaround = HTMLParser('')
            sage: setattr(coerce_workaround, '__call__', lambda e : ger2(P(e.polynomial())))
            sage: setattr(coerce_workaround, 'relations', lambda : ger2.relations()) 
            sage: sm._reduced_basis_polynomials(coerce_to = coerce_workaround)[0].parent().base_ring()
            Cyclotomic Field of order 6 and degree 2
        """
        if coerce_to is None :
            coerced_basis = self.__basis_in_graded_ambient
            relations = self.graded_ambient().relations()
        else :
            coerced_basis = map(coerce_to, self.__basis_in_graded_ambient)
            relations = coerce_to.relations()
        
        return Sequence( [ relations.reduce(b.polynomial())
                           for b in coerced_basis ],
                         universe = relations.ring() )
        
    @cached_method
    def _non_zero_monomials(self, coerce_to = None) :
        """
        A list of monomials which occur in the reduced polynomials
        associated with the basis of ``self``.
        
        INPUT:
            - ``coerce_to`` -- A graded ambient or ``None`` (default: ``None``);
                               Will be forwarded to :meth:~`._reduced_basis_polynomials`.
        
        OUTPUT:
            A sequence of monomials. 
        
        SEE::
            :meth:~`._reduced_basis_polynomials`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_module import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_submodule import *
            sage: K.<rho> = CyclotomicField(6)
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: P.<a, b, c> = PolynomialRing(K)
            sage: ger = GradedExpansionRing_class(None, Sequence([MonoidPowerSeries(mps, {1 : 4, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {1 : 1, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {2 : 1, 3: 6, 4 : 9}, mps.monoid().filter_all())]), PolynomialRing(QQ, ['a', 'b', 'c']).ideal(0), DegreeGrading((1,2,3)))
            sage: ger2 = GradedExpansionRing_class(None, Sequence([MonoidPowerSeries(mps, {1 : 4, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {1 : 1, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {2 : 1, 3: 6, 4 : 9}, mps.monoid().filter_all())]), P.ideal(a - b), DegreeGrading((1,2,3)))
            sage: sm = GradedExpansionSubmodule_abstract(ger, Sequence([ger.0, ger.1]), 3)
            sage: sm._non_zero_monomials()
            [a, b]
            sage: from htmllib import HTMLParser
            sage: coerce_workaround = HTMLParser('')
            sage: setattr(coerce_workaround, '__call__', lambda e : ger2(P(e.polynomial())))
            sage: setattr(coerce_workaround, 'relations', lambda : ger2.relations()) 
            sage: sm._non_zero_monomials(coerce_to = coerce_workaround)
            [b]
        """
        red_basis = self._reduced_basis_polynomials(coerce_to = coerce_to)
        
        return Sequence( reduce(union, [set(b.monomials()) for b in red_basis], set()),
                         universe = red_basis.universe() )

    @cached_method
    def _monomial_homomorphism(self, coerce_to = None) :
        """
        A homomorphism that maps the underlying module to a vector space
        where each component corresponds to a coefficient of the polynomial
        associated to elements of ``self`` within the graded ambient or
        ``coerce_to``. 

        INPUT:
            - ``coerce_to`` -- A graded ambient or ``None`` (default: ``None``);
                               Will be forwarded to :meth:~`._reduced_basis_polynomials`.
        
        OUTPUT:
            A morphism from ``self`` to a vector space.
        
        NOTE:
            The homomorphism is defined over a fraction field of the base ring
            of the monomials.
            
        SEE:
            :meth:~`._non_zero_monomials` and :meth:~`._reduced_basis_polynomials`.

        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_module import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_submodule import *
            sage: K.<rho> = CyclotomicField(6)
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: P.<a, b, c> = PolynomialRing(K)
            sage: ger = GradedExpansionRing_class(None, Sequence([MonoidPowerSeries(mps, {1 : 4, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {1 : 1, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {2 : 1, 3: 6, 4 : 9}, mps.monoid().filter_all())]), PolynomialRing(QQ, ['a', 'b', 'c']).ideal(0), DegreeGrading((1,2,3)))
            sage: ger2 = GradedExpansionRing_class(None, Sequence([MonoidPowerSeries(mps, {1 : 4, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {1 : 1, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {2 : 1, 3: 6, 4 : 9}, mps.monoid().filter_all())]), P.ideal(a - b), DegreeGrading((1,2,3)))
            sage: sm = GradedExpansionSubmodule_ambient_pid(ger, Sequence([ger.0, ger.1]), 3)
            sage: sm._monomial_homomorphism()
            Free module morphism defined by the matrix
            [1 0]
            [0 1]
            Domain: Submodule of Graded expansion ring with generators a, b, c
            Codomain: Vector space of dimension 2 over Rational Field
            sage: from htmllib import HTMLParser
            sage: coerce_workaround = HTMLParser('')
            sage: setattr(coerce_workaround, '__call__', lambda e : ger2(P(e.polynomial())))
            sage: setattr(coerce_workaround, 'relations', lambda : ger2.relations()) 
            sage: sm._monomial_homomorphism(coerce_to = coerce_workaround)
            Free module morphism defined by the matrix
            [1]
            [1]
            Domain: Submodule of Graded expansion ring with generators a, b, c
            Codomain: Vector space of dimension 1 over Cyclotomic Field of order 6 ...
        """
        reduced_basis = self._reduced_basis_polynomials(coerce_to = coerce_to)
        all_mons = self._non_zero_monomials(coerce_to = coerce_to)
        
        codomain = FreeModule(all_mons.universe().base_ring().fraction_field(), len(all_mons))
        basis_images = [ codomain([b.monomial_coefficient(m) for m in all_mons])
                         for b in reduced_basis ]

        return self.hom(basis_images)
          
          
    def _ambient_element_to_monomial_coefficients_generator(self, x, reduce = False, coerce_basis_to = None) :
        """
        Given an element `x` of the graded ambient ring or space return a generator
        corresponding to the image of `x` with respect to the morphism returned
        by :meth:~`._monomial_homomorphism`.
        
        INPUT:
            - `x`                 -- An element of the graded ambient.
            - ``reduce``          -- A boolean (default: ``False``); If ``True`` the polynomial
                                     attached to `x` will be reduced.
            - ``coerce_basis_to`` -- A graded ambient or ``None`` (default: ``None``);
                                     If not ``None`` the basis of ``self`` and `x` will be
                                     coerced before determining the monomials.
        
        OUTPUT:
            A generator over elements in the monomials' base ring.
          
        SEE:
            :meth:~`._monomial_homomorphism`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_module import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_submodule import *
            sage: K.<rho> = CyclotomicField(6)
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: P.<a, b, c> = PolynomialRing(K)
            sage: ger = GradedExpansionRing_class(None, Sequence([MonoidPowerSeries(mps, {1 : 4, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {1 : 1, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {2 : 1, 3: 6, 4 : 9}, mps.monoid().filter_all())]), PolynomialRing(QQ, ['a', 'b', 'c']).ideal(0), DegreeGrading((1,2,3)))
            sage: ger2 = GradedExpansionRing_class(None, Sequence([MonoidPowerSeries(mps, {1 : 4, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {1 : 1, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {2 : 1, 3: 6, 4 : 9}, mps.monoid().filter_all())]), P.ideal(a - b), DegreeGrading((1,2,3)))
            sage: sm = GradedExpansionSubmodule_ambient_pid(ger, Sequence([ger.0, ger.1]), 3)
            sage: list( sm._ambient_element_to_monomial_coefficients_generator(ger.2) )
            [0, 0]
            sage: from htmllib import HTMLParser
            sage: coerce_workaround = HTMLParser('')
            sage: setattr(coerce_workaround, '__call__', lambda e : ger2(P(e.polynomial())))
            sage: setattr(coerce_workaround, 'relations', lambda : ger2.relations()) 
            sage: list( sm._ambient_element_to_monomial_coefficients_generator(ger.2, coerce_basis_to = coerce_workaround) )
            [0]
            sage: list( sm._ambient_element_to_monomial_coefficients_generator(ger.1, reduce = True, coerce_basis_to = coerce_workaround) )
            [1]
        """ 
        if not coerce_basis_to is None :
            x = coerce_basis_to(x)
        if reduce :
            x = x._reduce_polynomial()
        else :
            x = x.polynomial()
        
        return ( x.monomial_coefficient(m)
                 for m in self._non_zero_monomials(coerce_to = coerce_basis_to) ) 

    def coordinates(self, x, in_base_ring = True, force_ambigous = False) :
        """
        The coordinates in ``self`` of an element either of the following:
            - The graded ambient,
            - An element of a submodule,
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
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_module import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_submodule import *
            sage: K.<rho> = CyclotomicField(6)
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: ger = GradedExpansionRing_class(None, Sequence([MonoidPowerSeries(mps, {1 : 4, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {1 : 1, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {2 : 1, 3: 6, 4 : 9}, mps.monoid().filter_all())]), PolynomialRing(QQ, ['a', 'b', 'c']).ideal(0), DegreeGrading((1,2,3)))
            sage: ger2 = GradedExpansionRing_class(None, Sequence([MonoidPowerSeries(mps, {1 : 4, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {1 : 1, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {2 : 1, 3: 6, 4 : 9}, mps.monoid().filter_all())]), PolynomialRing(K, ['a', 'b', 'c']).ideal(0), DegreeGrading((1,2,3)))
            sage: sm = GradedExpansionSubmodule_ambient_pid(ger, Sequence([ger.0, ger.1]))
            sage: sm2 = GradedExpansionSubmodule_ambient_pid(ger, Sequence([ger.0, ger.2]))
            sage: sm3 = GradedExpansionSubmodule_ambient_pid(ger2, Sequence([ger2.0, ger2.2]))
            sage: sm.coordinates(ger.2)
            Traceback (most recent call last):
            ...
            ArithmeticError: Graded expansion c is not contained in this space.
            sage: sm.coordinates(sm2.0)
            [1, 0]
            sage: sm.coordinates(sm2.1)
            Traceback (most recent call last):
            ...
            ArithmeticError: Graded expansion c is not contained in this space.
            sage: sm.coordinates(sm3.0)
            Traceback (most recent call last):
            ...
            ArithmeticError: No coordinates for (1, 0).
        """
        if isinstance(x, Element) :
            P = x.parent()
            if P is self :
                if self.ambient_module() is self :
                    return x.list()
                else :
                    return list(self.coordinate_vector(x))
            
            if self.graded_ambient().has_coerce_map_from(P) :
                if not P is self.graded_ambient() :
                    x = self.graded_ambient()(x)
                #x = x.polynomial()
                
                if not set(x._reduce_polynomial().monomials()).issubset(set(self._non_zero_monomials())) :
                    raise ArithmeticError( "%s is not contained in this space." % (x,) )
                
                mon_matrix = self._monomial_homomorphism().matrix().transpose()
                x_mon_vector = vector( self.base_ring().fraction_field(),
                                self._ambient_element_to_monomial_coefficients_generator(x, True) )
                
                    
                ## TODO: use linbox
                try :
                    coords = magma(mon_matrix.transpose()).Solution(magma(matrix(x_mon_vector))).sage()
                    coords = coords.row(0)
                except TypeError, msg :
                    if "Runtime error in 'Solution': No solution exists" in msg :
                        raise ArithmeticError( "%s is not contained in this space." % (x,) )
                    
                    try :
                        coords = mon_matrix.solve_right(x_mon_vector)
                    except ValueError :
                        raise ArithmeticError( "%s is not contained in this space." % (x,) )

                try :
                    # we used the base graded_ambient's fraction field, so convert it
                    return [self.base_ring()(c) for c in coords]
                except TypeError :
                    if in_base_ring :
                        raise ArithmeticError( "%s is not contained in this space." % (x,) )
                    else :
                        return coords.list()
                     
            #! elif self.graded_ambient().has_coerce_map_from(P) :
            elif P.has_coerce_map_from(self.graded_ambient()) :
                #x = x.polynomial()
                    
                mon_matrix = self._monomial_homomorphism().matrix().transpose()
                mon_matrix = matrix(P.base_ring().fraction_field(), mon_matrix)
                
                x_mon_vector = vector( P.base_ring().fraction_field(),
                                self._ambient_element_to_monomial_coefficients_generator(x, True, P) )
                
                try :
                    coords = mon_matrix.solve_right(x_mon_vector)
                except ValueError :
                    raise ArithmeticError( "%s is not contained in the image of this space." % (x,) )

                try :
                    # we used the base graded_ambient's fraction field, so convert it
                    return [self.base_ring()(c) for c in coords]
                except TypeError :
                    if in_base_ring :
                        raise ArithmeticError( "%s is not contained in the image of this space." % (x,) )
                    else :
                        return coords.list()
                     
            #! elif P.has_coerce_map_from(self.graded_ambient()) :
        #! if isinstance(x, (tuple, Element)) :
            
        return ExpansionModule_abstract.coordinates(self, x, in_base_ring, force_ambigous)
        
    def _graded_expansion_submodule_to_graded_ambient_(self, x) :
        """
        The element `x` of ``self`` as an element of the graded ambient ring or space.
        
        INPUT:
            - `x` -- An element of self.
        
        OUTPUT:
            An element of the graded ambient.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_module import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_submodule import *
            sage: K.<rho> = CyclotomicField(6)
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: P.<a, b, c> = PolynomialRing(K)
            sage: ger = GradedExpansionRing_class(None, Sequence([MonoidPowerSeries(mps, {1 : 4, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {1 : 1, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {2 : 1, 3: 6, 4 : 9}, mps.monoid().filter_all())]), PolynomialRing(QQ, ['a', 'b', 'c']).ideal(0), DegreeGrading((1,2,3)))
            sage: ger2 = GradedExpansionRing_class(None, Sequence([MonoidPowerSeries(mps, {1 : 4, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {1 : 1, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {2 : 1, 3: 6, 4 : 9}, mps.monoid().filter_all())]), P.ideal(a - b), DegreeGrading((1,2,3)))
            sage: sm = GradedExpansionSubmodule_ambient_pid(ger, Sequence([ger.0, ger.1]))
            sage: sm._graded_expansion_submodule_to_graded_ambient_(2*sm.0 + 3*sm.1)
            Graded expansion 2*a + 3*b
        """
        return sum( map(mul, self.coordinates(x), self._basis_in_graded_ambient()) )

    def _repr_(self) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_module import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_submodule import *
            sage: K.<rho> = CyclotomicField(6)
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: P.<a, b, c> = PolynomialRing(K)
            sage: ger = GradedExpansionRing_class(None, Sequence([MonoidPowerSeries(mps, {1 : 4, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {1 : 1, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {2 : 1, 3: 6, 4 : 9}, mps.monoid().filter_all())]), PolynomialRing(QQ, ['a', 'b', 'c']).ideal(0), DegreeGrading((1,2,3)))
            sage: ger2 = GradedExpansionRing_class(None, Sequence([MonoidPowerSeries(mps, {1 : 4, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {1 : 1, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {2 : 1, 3: 6, 4 : 9}, mps.monoid().filter_all())]), P.ideal(a - b), DegreeGrading((1,2,3)))
            sage: GradedExpansionSubmodule_ambient_pid(ger, Sequence([ger.0, ger.1]))
            Submodule of Graded expansion ring with generators a, b, c
        """
        return "Submodule of %s" % (self.__graded_ambient,)
    
    def _latex_(self) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_module import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_submodule import *
            sage: K.<rho> = CyclotomicField(6)
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: P.<a, b, c> = PolynomialRing(K)
            sage: ger = GradedExpansionRing_class(None, Sequence([MonoidPowerSeries(mps, {1 : 4, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {1 : 1, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {2 : 1, 3: 6, 4 : 9}, mps.monoid().filter_all())]), PolynomialRing(QQ, ['a', 'b', 'c']).ideal(0), DegreeGrading((1,2,3)))
            sage: ger2 = GradedExpansionRing_class(None, Sequence([MonoidPowerSeries(mps, {1 : 4, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {1 : 1, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {2 : 1, 3: 6, 4 : 9}, mps.monoid().filter_all())]), P.ideal(a - b), DegreeGrading((1,2,3)))
            sage: latex( GradedExpansionSubmodule_ambient_pid(ger, Sequence([ger.0, ger.1])) )
            Submodule of Graded expansion ring with generators a, b, c
        """
        return "Submodule of %s" % (latex(self.__graded_ambient),)

#===============================================================================
# GradedExpansionSubmodule_generic
#===============================================================================

class GradedExpansionSubmodule_generic ( GradedExpansionSubmodule_abstract, ExpansionModule_generic ) :
    """
    A finite dimensional module of graded expansions within an ambient.
    
    SEE:
        :class:~`fourier_expansion_framework.gradedexpansions.gradedexpansion_ambient.GradedExpansion_ambient`.
    """
    
    def __init__(self, graded_ambient, basis, degree, **kwds) :
        """
        INPUT:
            - ``graded_ambient`` -- An instance of :class:~`fourier_expansion_framework.gradedexpansions.gradedexpansion_ambient.GradedExpansionAmbient_abstract`.
            - ``basis``          -- A tuple, list or sequence of elements of the graded ambient.
            - ``degree``         -- An integer; The degree of the module within its ambient module.
        
        NOTE:
            The base ring of the graded expansion ambient must be an integral domain.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_module import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_submodule import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: ger = GradedExpansionRing_class(None, Sequence([MonoidPowerSeries(mps, {1 : 4, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {1 : 1, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {2 : 1, 3: 6, 4 : 9}, mps.monoid().filter_all())]), PolynomialRing(QQ, ['a', 'b', 'c']).ideal(0), DegreeGrading((1,2,3)))
            sage: sm = GradedExpansionSubmodule_generic(ger, Sequence([ger.0, ger.1]), 3)
            sage: sm = GradedExpansionSubmodule_generic(ger, Sequence([ger.0, ger.1]), 2, _element_class = GradedExpansionSubmoduleVector_generic, no_expansion_init = True)
        """
        if not hasattr(self, '_element_class') :
            if '_element_class' in kwds :
                self._element_class = kwds['_element_class']
            else :
                self._element_class = GradedExpansionSubmoduleVector_generic

        GradedExpansionSubmodule_abstract.__init__(self, graded_ambient, basis, degree)
        ExpansionModule_generic.__init__(self, self._abstract_basis(), degree, **kwds)

    def change_ring(self, R):
        """
        Return the ambient module of graded expansions over `R` with the same basis as ``self``.
        
        INPUT:
            - `R` -- A ring.
        
        OUTPUT:
            An instance of :class:~`.GradedExpansionSubmodule_generic`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_module import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_submodule import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: ger = GradedExpansionRing_class(None, Sequence([MonoidPowerSeries(mps, {1 : 4, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {1 : 1, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {2 : 1, 3: 6, 4 : 9}, mps.monoid().filter_all())]), PolynomialRing(QQ, ['a', 'b', 'c']).ideal(0), DegreeGrading((1,2,3)))
            sage: sm = GradedExpansionSubmodule_generic(ger, Sequence([ger.0, ger.1]), 3)
            sage: sm.change_ring(QQ) is sm
            True
            sage: sm.change_ring(PolynomialRing(QQ, 'x'))
            Traceback (most recent call last):
            ...
            ValueError: Associated expansion of graded ambient not defined over Univariate Polynomial Ring in x over Rational Field.
        """
        if self.base_ring() == R:
            return self

        raise ValueError( "Associated expansion of graded ambient not defined over %s." % R )

    def basis(self) :
        """
        A basis of ``self``.
        
        OUTPUT:
            A list of elements of ``self``.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_module import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_submodule import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: ger = GradedExpansionRing_class(None, Sequence([MonoidPowerSeries(mps, {1 : 4, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {1 : 1, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {2 : 1, 3: 6, 4 : 9}, mps.monoid().filter_all())]), PolynomialRing(QQ, ['a', 'b', 'c']).ideal(0), DegreeGrading((1,2,3)))
            sage: sm = GradedExpansionSubmodule_generic(ger, Sequence([ger.0, ger.1]), 3)
            sage: sm.basis()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError    

    def hom(self, other) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_module import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_submodule import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: ger = GradedExpansionRing_class(None, Sequence([MonoidPowerSeries(mps, {1 : 4, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {1 : 1, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {2 : 1, 3: 6, 4 : 9}, mps.monoid().filter_all())]), PolynomialRing(QQ, ['a', 'b', 'c']).ideal(0), DegreeGrading((1,2,3)))
            sage: sm = GradedExpansionSubmodule_generic(ger, Sequence([ger.0, ger.1]), 3)
            sage: m = FreeModule(QQ, 2)
            sage: sm.hom([m([1,0]),m([1,1])])
            Free module morphism defined by the matrix
            [1 0]
            [1 1]
            Domain: Submodule of Graded expansion ring with generators a, b, c
            Codomain: Vector space of dimension 2 over Rational Field
            sage: sm.hom(sm)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if isinstance(other, GradedExpansionSubmodule_abstract) and \
           self.graded_ambient() == other.graded_ambient() :
            raise NotImplementedError
        else :
            images = other
        
        return FreeModule_generic.hom(self, images)

#===============================================================================
# _span
# This will be used in GradedExpansionSubmodule_ambient_pid and
# GradedExpansionSubmodule_submodule_pid 
#===============================================================================

def _span( self, gens, base_ring = None, check = True, already_echelonized = False ) :
    """
    The submodule of graded expansions spanned by ``gens``.
    
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
        An instance of :class:~`.GradedExpansionSubmodule_submodule_pid`.
    
    TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_module import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_submodule import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: ger = GradedExpansionRing_class(None, Sequence([MonoidPowerSeries(mps, {1 : 4, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {1 : 1, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {2 : 1, 3: 6, 4 : 9}, mps.monoid().filter_all())]), PolynomialRing(QQ, ['a', 'b', 'c']).ideal(0), DegreeGrading((1,2,3)))
            sage: sm = GradedExpansionSubmodule_ambient_pid(ger, Sequence([ger.0, ger.1]))
            sage: sms = sm.span([sm.0])
            sage: sms = sm.span([sm.0], base_ring = PolynomialRing(QQ, 'x'))
            Traceback (most recent call last):
            ...
            ValueError: Associated expansion of graded ambient not defined over Univariate Polynomial Ring in x over Rational Field.
    """
    if is_FreeModule(gens):
        gens = gens.gens()
    if not isinstance(gens, (list, tuple, Sequence)):
        raise TypeError, "Argument gens (= %s) must be a list, tuple, or sequence."%gens
    
    if base_ring is None or base_ring == self.base_ring():
        gens = Sequence(gens, check = check, universe = self.ambient_module())
        
        return GradedExpansionSubmodule_submodule_pid( self.ambient_module(), gens )
    else:
        try:
            M = self.change_ring(base_ring)
        except TypeError:
            raise ValueError, "Argument base_ring (= %s) is not compatible " % (base_ring,) + \
                "with the base field (= %s)." % (self.base_field(),)
        try: 
            return M.span(gens)
        except TypeError:
            raise ValueError, "Argument gens (= %s) is not compatible " % (gens,) + \
                "with base_ring (= %s)." % (base_ring,)

#===============================================================================
# GradedExpansionSubmodule_ambient_pid
#===============================================================================

class GradedExpansionSubmodule_ambient_pid ( GradedExpansionSubmodule_abstract, ExpansionModule_ambient_pid ) :
    """
    An ambient module of graded expansions expansions over a principal ideal domain.
    """
    
    def __init__(self, graded_ambient, basis, _element_class = None, **kwds) :
        """
        INPUT:
            - ``graded_ambient`` -- An instance of :class:~`fourier_expansion_framework.gradedexpansions.gradedexpansion_ambient.GradedExpansionAmbient_abstract`.
            - ``basis``          -- A list or sequence of (equivariant) monoid power series.
            - ``element_class``  -- A type or ``None`` (default: ``None``); The element class
                                    attached to ``self``. If ``None`` :class:~`.GradedExpansionSubmoduleVector_generic`
                                    will be used.
            - ``kwds``           -- Will be forwarded to :class:~`fourier_expansion_module.gradedexpansions.expansion_module.ExpansionModule_ambient_pid`
                                    and :class:~`.GradedExpansionSubmodule_abstract`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_module import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_submodule import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: ger = GradedExpansionRing_class(None, Sequence([MonoidPowerSeries(mps, {1 : 4, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {1 : 1, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {2 : 1, 3: 6, 4 : 9}, mps.monoid().filter_all())]), PolynomialRing(QQ, ['a', 'b', 'c']).ideal(0), DegreeGrading((1,2,3)))
            sage: sm = GradedExpansionSubmodule_ambient_pid(ger, Sequence([ger.0, ger.1]))
            sage: sm = GradedExpansionSubmodule_ambient_pid(ger, Sequence([ger.0, ger.1]), _element_class = GradedExpansionSubmoduleVector_generic, no_expansion_init = True)
            Traceback (most recent call last):
            ...
            AttributeError: 'GradedExpansionSubmodule_ambient_pid' object has no attribute '_ExpansionModule_abstract__abstract_basis'
        """
        if not hasattr(self, '_element_class') :
            if not _element_class is None :
                self._element_class = _element_class
            else :
                self._element_class = GradedExpansionSubmoduleVector_generic

        GradedExpansionSubmodule_abstract.__init__(self, graded_ambient, basis, len(basis), **kwds)
        ExpansionModule_ambient_pid.__init__(self, self._abstract_basis(), **kwds)
        
    global _span
    span = _span
    
    def change_ring(self, R):
        """
        Return the ambient module of graded expansions over `R` with the
        same basis as ``self``.
        
        INPUT:
            - `R` -- A ring.
        
        OUTPUT:
            An instance of :class:~`.GradedExpansionSubmodule_ambient_pid`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_module import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_submodule import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: ger = GradedExpansionRing_class(None, Sequence([MonoidPowerSeries(mps, {1 : 4, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {1 : 1, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {2 : 1, 3: 6, 4 : 9}, mps.monoid().filter_all())]), PolynomialRing(QQ, ['a', 'b', 'c']).ideal(0), DegreeGrading((1,2,3)))
            sage: sm = GradedExpansionSubmodule_generic(ger, Sequence([ger.0, ger.1]), 3)
            sage: sm.change_ring(QQ) is sm
            True
            sage: sm.change_ring(PolynomialRing(QQ, 'x'))
            Traceback (most recent call last):
            ...
            ValueError: Associated expansion of graded ambient not defined over Univariate Polynomial Ring in x over Rational Field.
        """
        if self.base_ring() == R:
            return self
            
        raise ValueError( "Associated expansion of graded ambient not defined over %s." % R )

    def hom(self, other) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_module import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_submodule import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: ger = GradedExpansionRing_class(None, Sequence([MonoidPowerSeries(mps, {1 : 4, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {1 : 1, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {2 : 1, 3: 6, 4 : 9}, mps.monoid().filter_all())]), PolynomialRing(QQ, ['a', 'b', 'c']).ideal(0), DegreeGrading((1,2,3)))
            sage: sm = GradedExpansionSubmodule_ambient_pid(ger, Sequence([ger.0, ger.1]))
            sage: sm2 = GradedExpansionSubmodule_ambient_pid(ger, Sequence([ger.0, ger.2, ger.1]))
            sage: m = FreeModule(QQ, 2)
            sage: sm.hom(sm2)
            Free module morphism defined by the matrix
            [1 0 0]
            [0 0 1]
            Domain: Submodule of Graded expansion ring with generators a, b, c
            Codomain: Submodule of Graded expansion ring with generators a, b, c
            sage: sm.hom([m([1,0]),m([1,1])])
            Free module morphism defined by the matrix
            [1 0]
            [1 1]
            Domain: Submodule of Graded expansion ring with generators a, b, c
            Codomain: Vector space of dimension 2 over Rational Field
        """
        if isinstance(other, GradedExpansionSubmodule_abstract) and \
           self.graded_ambient() == other.graded_ambient() :
            images = [other(b) for b in self._basis_in_graded_ambient()]
        else :
            images = other
        
        return FreeModule_ambient_pid.hom(self, images)            

#===============================================================================
# GradedExpansionSubmodule_submodule_pid
#===============================================================================

class GradedExpansionSubmodule_submodule_pid ( GradedExpansionSubmodule_abstract, ExpansionModule_submodule_pid ) :
    """
    A submodule of another module of graded expansions over a principal ideal domain.
    """

    def __init__(self, ambient, gens, element_class = None, **kwds) :
        """
        INPUT:
            - ``ambient``       -- An instance of :class:~`.GradedExpansionSubmodule_submodule_pid` :class:~`.GradedExpansionSubmodule_submodule_pid` or .
            - ``gens``          -- A tuple, list or sequence of elements of ``ambient``.
            - ``element_class`` -- A type or ``None`` (default: ``None``); The element class
                                   attached to ``self``. If ``None`` :class:~`.GradedExpansionSubmoduleVector_generic`
                                   will be used.
            - ``kwds``          -- Will be forwarded to :class:~`fourier_expansion_module.gradedexpansions.expansion_module.ExpansionModule_submodule_pid`
                                   and :class:~`.GradedExpansionSubmodule_abstract`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_module import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_submodule import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: ger = GradedExpansionRing_class(None, Sequence([MonoidPowerSeries(mps, {1 : 4, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {1 : 1, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {2 : 1, 3: 6, 4 : 9}, mps.monoid().filter_all())]), PolynomialRing(QQ, ['a', 'b', 'c']).ideal(0), DegreeGrading((1,2,3)))
            sage: sm = GradedExpansionSubmodule_ambient_pid(ger, Sequence([ger.0, ger.1]))
            sage: sms = GradedExpansionSubmodule_submodule_pid(sm, [sm.0, sm.1])
            sage: sms = GradedExpansionSubmodule_submodule_pid(sm, [sm.0, sm.1], _element_class = GradedExpansionSubmoduleVector_generic, no_expansion_init = True)
        """
        if not hasattr(self, '_element_class') :
            if not element_class is None :
                self._element_class = element_class
            else :
                self._element_class = GradedExpansionSubmoduleVector_generic

        GradedExpansionSubmodule_abstract.__init__( self, ambient.graded_ambient(),
                                                    [ambient.graded_ambient()(g) for g in gens], ambient.dimension(),
                                                    **kwds )
        ExpansionModule_submodule_pid.__init__(self, ambient, gens, **kwds)
        
    global _span
    span = _span
    
    def change_ring(self, R):
        """
        Return the ambient module of graded expansions over `R` with the
        same basis as ``self``.
        
        INPUT:
            - `R` -- A ring.
        
        OUTPUT:
            An instance of :class:~`.GradedExpansionSubmodule_ambient_pid`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_module import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_submodule import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: ger = GradedExpansionRing_class(None, Sequence([MonoidPowerSeries(mps, {1 : 4, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {1 : 1, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {2 : 1, 3: 6, 4 : 9}, mps.monoid().filter_all())]), PolynomialRing(QQ, ['a', 'b', 'c']).ideal(0), DegreeGrading((1,2,3)))
            sage: sm = GradedExpansionSubmodule_ambient_pid(ger, Sequence([ger.0, ger.1]))
            sage: sms = GradedExpansionSubmodule_submodule_pid(sm, [sm.0, sm.1])
            sage: sms.change_ring(QQ) is sms
            True
            sage: sms.change_ring(PolynomialRing(QQ, 'x'))
            Traceback (most recent call last):
            ...
            ValueError: Associated expansion of graded ambient not defined over Univariate Polynomial Ring in x over Rational Field.
        """
        if self.base_ring() == R:
            return self
                    
        raise ValueError( "Associated expansion of graded ambient not defined over %s." % R )

    def hom(self, other) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_module import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_submodule import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: ger = GradedExpansionRing_class(None, Sequence([MonoidPowerSeries(mps, {1 : 4, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {1 : 1, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {2 : 1, 3: 6, 4 : 9}, mps.monoid().filter_all())]), PolynomialRing(QQ, ['a', 'b', 'c']).ideal(0), DegreeGrading((1,2,3)))
            sage: sm = GradedExpansionSubmodule_ambient_pid(ger, Sequence([ger.0, ger.1]))
            sage: sms = GradedExpansionSubmodule_submodule_pid(sm, [sm.0, sm.1])
            sage: sms2 = GradedExpansionSubmodule_submodule_pid(sm, [sm.0, sm.0 + sm.1])
            sage: m = FreeModule(QQ, 2)
            sage: sms.hom(sms2)
            Free module morphism defined by the matrix
            [1 0]
            [0 1]
            Domain: Submodule of Graded expansion ring with generators a, b, c
            Codomain: Submodule of Graded expansion ring with generators a, b, c
            sage: sms.hom([m([1,0]),m([1,1])])
            Free module morphism defined by the matrix
            [1 0]
            [1 1]
            Domain: Submodule of Graded expansion ring with generators a, b, c
            Codomain: Vector space of dimension 2 over Rational Field
        """
        if isinstance(other, GradedExpansionSubmodule_abstract) and \
           self.graded_ambient() == other.graded_ambient() :
            images = [other(b) for b in self._basis_in_graded_ambient()]
        else :
            images = other
        
        return FreeModule_submodule_pid.hom(self, images)
    
#===============================================================================
# GradedExpansionSubmoduleVector_generic
#===============================================================================

class GradedExpansionSubmoduleVector_generic ( ExpansionModuleVector_generic ) :
    """
    A generic implementation of an element in a module of graded expansions.
    """
    
    def __copy__(self) :
        """
        Return a copy of ``self``.
        
        OUTPUT:
            An instance of :class:~`.GradedExpansionSubmoduleVector_generic`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_module import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_submodule import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: ger = GradedExpansionRing_class(None, Sequence([MonoidPowerSeries(mps, {1 : 4, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {1 : 1, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {2 : 1, 3: 6, 4 : 9}, mps.monoid().filter_all())]), PolynomialRing(QQ, ['a', 'b', 'c']).ideal(0), DegreeGrading((1,2,3)))
            sage: sm = GradedExpansionSubmodule_ambient_pid(ger, Sequence([ger.0, ger.1]))
            sage: h = GradedExpansionSubmoduleVector_generic(sm, [1,0])
            sage: copy(h)
            (1, 0)
        """
        return GradedExpansionSubmoduleVector_generic( self.parent(), self.list(), coerce = False, copy = True )


    def polynomial(self) :
        """
        Return the polynomial associated to ``self`` in the graded ambient
        of its parent.
        
        OUTPUT:
            A polynomial in the polynomial ring attached to the graded ambient.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_module import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_submodule import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: ger = GradedExpansionRing_class(None, Sequence([MonoidPowerSeries(mps, {1 : 4, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {1 : 1, 2 : 3}, mps.monoid().filter_all()), MonoidPowerSeries(mps, {2 : 1, 3: 6, 4 : 9}, mps.monoid().filter_all())]), PolynomialRing(QQ, ['a', 'b', 'c']).ideal(0), DegreeGrading((1,2,3)))
            sage: sm = GradedExpansionSubmodule_ambient_pid(ger, Sequence([ger.0, ger.1]))
            sage: h = GradedExpansionSubmoduleVector_generic(sm, [1,0])
            sage: h.polynomial()
            a
        """
        return sum( self[k]*self.parent()._basis_in_graded_ambient()[k].polynomial()
                    for k in xrange(self.parent().rank()) )
