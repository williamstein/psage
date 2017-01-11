r"""
Elements wrapping a Fourier expansion which partially known relations.

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

from psage.modform.fourier_expansion_framework.gradedexpansions.fourierexpansionwrapper import FourierExpansionWrapper
from itertools import groupby
from sage.structure.element import AlgebraElement,ModuleElement
from sage.misc.all import prod
from sage.misc.latex import latex
from sage.rings.infinity import infinity
import operator

#===============================================================================
# GradedExpansion_abstract
#===============================================================================

class GradedExpansion_abstract ( FourierExpansionWrapper ) :
    r"""
    An abstract graded expansion.
    """
    
    def __init__(self, parent, polynomial) :
        r"""
        INPUT:
            - ``parent``     -- An instance of :class:~`fourier_expansion_framework.gradedexpansions.gradedexpansion_ambient.GradedExpansionAmbient_abstract`.
            - ``polynomial`` -- A polynomial in the polynomial ring underlying the parent.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_element import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: P.<a,b> = QQ[]
            sage: ger = GradedExpansionRing_class(Sequence([MonoidPowerSeries(mps, {1: 1}, mps.monoid().filter(4))]), Sequence([MonoidPowerSeries(mps, {1: 1, 2: 3}, mps.monoid().filter(4))]), P.ideal(0), DegreeGrading((1,2)))
            sage: h = GradedExpansion_abstract(ger, P(0))
            sage: h = GradedExpansion_abstract(ger, a*b)
        """
        self.__polynomial = polynomial
        
    def polynomial(self) :
        r"""
        The underlying polynomial.
        
        OUTPUT:
            A polynomial in the polynomial ring underlying the parent.
        
        NOTE:
            This polynomial might change by elements of the parent's
            relation ideal.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_element import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: P.<a,b> = QQ[]
            sage: ger = GradedExpansionRing_class(Sequence([MonoidPowerSeries(mps, {1: 1}, mps.monoid().filter(4))]), Sequence([MonoidPowerSeries(mps, {1: 1, 2: 3}, mps.monoid().filter(4))]), P.ideal(0), DegreeGrading((1,2)))
            sage: h = GradedExpansion_abstract(ger, a * b)
            sage: h.polynomial()
            a*b
        """
        return self.__polynomial
    
    def _reduce_polynomial(self, set_polynomial = True) :
        r"""
        The Groebner reduced underlying polynomial.

        INPUT:
            - ``set_polynomial`` -- A boolean (default: ``True``); If ``True``
                                    the underlying polynomial will be set to the
                                    reduced polynomial.
        OUTPUT:
            A polynomial in the polynomial ring underlying the parent.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_element import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: P.<a,b> = QQ[]
            sage: ger = GradedExpansionRing_class(Sequence([MonoidPowerSeries(mps, {1: 1}, mps.monoid().filter(4))]), Sequence([MonoidPowerSeries(mps, {1: 1, 2: 3}, mps.monoid().filter(4))]), P.ideal(a - b), DegreeGrading((1,2)))
            sage: h = GradedExpansion_class(ger, a - b)
            sage: h._reduce_polynomial(set_polynomial = False)
            0
            sage: h.polynomial()
            a - b
            sage: h._reduce_polynomial()
            0
            sage: h.polynomial()
            0
        """
        if set_polynomial :
            self.__polynomial = self.parent().relations().reduce(self.__polynomial)
        
            return self.__polynomial
        else :
            return self.parent().relations().reduce(self.__polynomial)
    
    def homogeneous_components(self) :
        r"""
        Split the underlying polynomial with respect to grading imposed on
        the parent.
        
        OUTPUT:
           A dictionary with key the grading values and values the polynomials.
        
        EXAMPLES::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_element import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: P.<a,b> = QQ[]
            sage: ger = GradedExpansionRing_class(Sequence([MonoidPowerSeries(mps, {1: 1}, mps.monoid().filter(4))]), Sequence([MonoidPowerSeries(mps, {1: 1, 2: 3}, mps.monoid().filter(4))]), P.ideal(0), DegreeGrading((1,2)))
            sage: h = GradedExpansion_class(ger, a - b)
            sage: h.homogeneous_components()
            {1: a, 2: -b}

        TESTS::
            sage: h = GradedExpansion_class(ger, P(0))
            sage: h.homogeneous_components()
            {}
        """
        if self.__polynomial.is_zero() :
            return dict()
        
        components = dict()
        grading = self.parent().grading().index
        if self.__polynomial.parent().ngens() == 1 :
            normalize = lambda e: (e,)
            var = self.__polynomial.parent().gen(0)
            monomial = lambda e: var**e
        else :
            normalize = lambda e: e
            vars = self.__polynomial.parent().gens()
            monomial = lambda e: prod( map(operator.pow, vars, e) )
             
        for e in self.__polynomial.exponents() :
            ne = normalize(e)
            m = monomial(e)
            try :
                components[grading(ne)] = components[grading(ne)] + m * self.__polynomial[e]
            except KeyError :
                components[grading(ne)] = m * self.__polynomial[e]
        
        return components

    def exponents(self) :
        r"""
        The exponents of the underlying polynomial.
        
        OUTPUT:
            A list of tuples.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_element import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: P.<a,b> = QQ[]
            sage: ger = GradedExpansionRing_class(Sequence([MonoidPowerSeries(mps, {1: 1}, mps.monoid().filter(4))]), Sequence([MonoidPowerSeries(mps, {1: 1, 2: 3}, mps.monoid().filter(4))]), P.ideal(a - b), DegreeGrading((1,2)))
            sage: h = GradedExpansion_abstract(ger, a - b)
            sage: h.exponents()
            [(1, 0), (0, 1)]
            sage: P.<a> = QQ[]
            sage: ger = GradedExpansionRing_class(None, Sequence([MonoidPowerSeries(mps, {1: 1, 2: 3}, mps.monoid().filter(4))]), P.ideal(0), DegreeGrading((1,)))
            sage: h = GradedExpansion_abstract(ger, a)
            sage: h.exponents()
            [(1,)]
        """
        if self.polynomial().parent().ngens() == 1 :
            return map(lambda e: (e,), self.polynomial().exponents())
        else :
            return map(tuple, self.polynomial().exponents())

    def grading_index(self) :
        r"""
        If ``self`` is homogeneous, return the index of ``self`` with respect to
        the grading imposed on the parent. Otherwise, a ``ValueError`` is raised.
        
        OUTPUT:
            A grading value.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_element import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: P.<a,b> = QQ[]
            sage: ger = GradedExpansionRing_class(Sequence([MonoidPowerSeries(mps, {1: 1}, mps.monoid().filter(4))]), Sequence([MonoidPowerSeries(mps, {1: 1, 2: 3}, mps.monoid().filter(4))]), P.ideal(a - b), DegreeGrading((1,2)))
            sage: h = GradedExpansion_class(ger, a^2)
            sage: h.grading_index()
            2
            sage: h = GradedExpansion_class(ger, P(0))
            sage: h.grading_index()
            +Infinity
            sage: h = GradedExpansion_class(ger, a - b)
            sage: h.grading_index()
            Traceback (most recent call last):
            ...
            ValueError: No homogeneous weight.
        """
        cps = self.homogeneous_components()
        if len(cps) == 0 :
            return infinity
        elif len(cps) != 1 :
            raise ValueError( "No homogeneous weight." )
        
        return cps.keys()[0]
    
    def _add_(left, right) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_element import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: P.<a,b> = QQ[]
            sage: ger = GradedExpansionRing_class(Sequence([MonoidPowerSeries(mps, {1: 1}, mps.monoid().filter(4))]), Sequence([MonoidPowerSeries(mps, {1: 1, 2: 3}, mps.monoid().filter(4))]), P.ideal(a - b), DegreeGrading((1,2)))
            sage: h = GradedExpansion_class(ger, a^2)
            sage: l = GradedExpansion_class(ger, b)
            sage: (h + l).polynomial()
            a^2 + b
        """
        return left.__class__( left.parent(),
                left.polynomial() + right.polynomial() )
        
    def _lmul_(self, c) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_element import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: P.<a,b> = QQ[]
            sage: ger = GradedExpansionRing_class(Sequence([MonoidPowerSeries(mps, {1: 1}, mps.monoid().filter(4))]), Sequence([MonoidPowerSeries(mps, {1: 1, 2: 3}, mps.monoid().filter(4))]), P.ideal(a - b), DegreeGrading((1,2)))
            sage: h = GradedExpansion_class(ger, a - b)
            sage: (h*ger.base_ring()(4)).polynomial()
            4*a - 4*b
        """
        if self.parent().nbasegens() != 0 :
            return self.__class__( self.parent(),
                    self.polynomial() * c.polynomial() )
        else :
            return self.__class__( self.parent(),
                    self.polynomial() * c )
        
    def _rmul_(self, c) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_element import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: P.<a,b> = QQ[]
            sage: ger = GradedExpansionRing_class(Sequence([MonoidPowerSeries(mps, {1: 1}, mps.monoid().filter(4))]), Sequence([MonoidPowerSeries(mps, {1: 1, 2: 3}, mps.monoid().filter(4))]), P.ideal(a - b), DegreeGrading((1,2)))
            sage: h = GradedExpansion_class(ger, a^2)
            sage: (ger.base_ring()(4) * h).polynomial()
            4*a^2
        """
        if self.parent().nbasegens() != 0 :
            return self.__class__( self.parent(),
                    c.polynomial() * self.polynomial() )
        else :
            return self.__class__( self.parent(),
                    c * self.polynomial() )
            
    def __cmp__(self, other) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_element import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: P.<a,b> = QQ[]
            sage: ger = GradedExpansionRing_class(Sequence([MonoidPowerSeries(mps, {1: 1}, mps.monoid().filter(4))]), Sequence([MonoidPowerSeries(mps, {1: 1, 2: 3}, mps.monoid().filter(4))]), P.ideal(a - b), DegreeGrading((1,2)))
            sage: h = GradedExpansion_class(ger, a^2)
            sage: h == GradedExpansion_class(ger, a^2)
            True
            sage: h == GradedExpansion_class(ger, a - b)
            False
            sage: h = GradedExpansion_class(ger, P(0))
            sage: h == GradedExpansion_class(ger, a - b)
            True
            sage: ger = GradedExpansionRing_class(Sequence([MonoidPowerSeries(mps, {1: 1}, mps.monoid().filter(4))]), Sequence([MonoidPowerSeries(mps, {1: 1, 2: 3}, mps.monoid().filter(4))]), P.ideal(0), DegreeGrading((1,2)))
            sage: h = GradedExpansion_class(ger, a^2)
            sage: h == GradedExpansion_class(ger, a^2)
            True
            sage: h = GradedExpansion_class(ger, P(0))
            sage: h == GradedExpansion_class(ger, a - b)
            False
            sage: ger = GradedExpansionRing_class(Sequence([MonoidPowerSeries(mps, {1: 1}, mps.monoid().filter(4))]), Sequence([MonoidPowerSeries(mps, {1: 1, 2: 3}, mps.monoid().filter(4))]), P.ideal(a - b), DegreeGrading((1,2)), all_relations = False)
            sage: h = GradedExpansion_class(ger, a^2)
            sage: h == GradedExpansion_class(ger, a^2)
            False
        """
        c = cmp(type(self), type(other))
        
        if c == 0 :
            if self.parent().all_relations_known() :
                if self.parent().has_relation_free_generators() :
                    c = cmp(self.polynomial(), other.polynomial())
                else :
                    c = cmp(self.parent().relations().reduce(self.polynomial() - other.polynomial()), 0)
            else :
                c = 1

        return c

#===============================================================================
# GradedExpansion_class
#===============================================================================

class GradedExpansion_class ( GradedExpansion_abstract, AlgebraElement ) :
    r"""
    A graded expansion which is element of an algebra.
    
    SEE:
        :class:`fourier_expansion_framework.gradedexpansions.gradedexpansion_ring.GradedExpansionRing_class`.
    """
    def __init__(self, parent, polynomial) :
        r"""
        INPUT:
            - ``parent``     -- An instance of :class:~`fourier_expansion_framework.gradedexpansions.gradedexpansion_ring.GradedExpansionRing_class`.
            - ``polynomial`` -- A polynomial in the polynomial ring underlying the parent.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_element import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: P.<a,b> = QQ[]
            sage: ger = GradedExpansionRing_class(Sequence([MonoidPowerSeries(mps, {1: 1}, mps.monoid().filter(4))]), Sequence([MonoidPowerSeries(mps, {1: 1, 2: 3}, mps.monoid().filter(4))]), P.ideal(a - b), DegreeGrading((1,2)))
            sage: h = GradedExpansion_class(ger, a^2)
        """
        AlgebraElement.__init__(self, parent)
        GradedExpansion_abstract.__init__(self, parent, polynomial)        
        
    def _mul_(left, right) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_element import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: P.<a,b> = QQ[]
            sage: ger = GradedExpansionRing_class(Sequence([MonoidPowerSeries(mps, {1: 1}, mps.monoid().filter(4))]), Sequence([MonoidPowerSeries(mps, {1: 1, 2: 3}, mps.monoid().filter(4))]), P.ideal(a - b), DegreeGrading((1,2)))
            sage: h = GradedExpansion_class(ger, a^2)
            sage: (h * h).polynomial()
            a^4
        """
        return left.__class__( left.parent(),
                left.polynomial() * right.polynomial() )

    def __hash__(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_element import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: P.<a,b> = QQ[]
            sage: ger = GradedExpansionRing_class(Sequence([MonoidPowerSeries(mps, {1: 1}, mps.monoid().filter(4))]), Sequence([MonoidPowerSeries(mps, {1: 1, 2: 3}, mps.monoid().filter(4))]), P.ideal(a - b), DegreeGrading((1,2)))
            sage: h = GradedExpansion_class(ger, a^2)
            sage: hash(h)
            -1239234133 # 32-bit
            12415370528999851 # 64-bit
        """
        return hash(self.polynomial())
    
    def _repr_(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_element import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: P.<a,b> = QQ[]
            sage: ger = GradedExpansionRing_class(Sequence([MonoidPowerSeries(mps, {1: 1}, mps.monoid().filter(4))]), Sequence([MonoidPowerSeries(mps, {1: 1, 2: 3}, mps.monoid().filter(4))]), P.ideal(a - b), DegreeGrading((1,2)))
            sage: GradedExpansion_class(ger, a^2)
            Graded expansion a^2
        """
        return "Graded expansion %s" % (self.polynomial(),)
    
    def _latex_(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_element import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: P.<a,b> = QQ[]
            sage: ger = GradedExpansionRing_class(Sequence([MonoidPowerSeries(mps, {1: 1}, mps.monoid().filter(4))]), Sequence([MonoidPowerSeries(mps, {1: 1, 2: 3}, mps.monoid().filter(4))]), P.ideal(a - b), DegreeGrading((1,2)))
            sage: latex( GradedExpansion_class(ger, a^2) )
            \text{Graded expansion }a^{2}
        """
        return r"\text{Graded expansion }%s" % latex(self.polynomial())

class GradedExpansionVector_class ( GradedExpansion_abstract, ModuleElement ) :
    r"""
    A graded expansion which is element of a module.
        
    SEE:
        :class:`fourier_expansion_framework.gradedexpansions.gradedexpansion_module.GradedExpansionModule_class`.
    """

    def __init__(self, parent, polynomial) :
        r"""
        INPUT:
            - ``parent``     -- An instance of :class:~`fourier_expansion_framework.gradedexpansions.gradedexpansion_module.GradedExpansionModule_class`.
            - ``polynomial`` -- A polynomial in the polynomial ring underlying the parent.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_module import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_element import *
            sage: m = FreeModule(QQ, 3)
            sage: mpsm = MonoidPowerSeriesModule(m, NNMonoid(False))
            sage: mps = mpsm.base_ring()
            sage: P.<a,b> = QQ[]
            sage: ger = GradedExpansionModule_class(Sequence([MonoidPowerSeries(mps, {1: 1}, mps.monoid().filter(4))]), Sequence([MonoidPowerSeries(mpsm, {1: m([1,1,1]), 2: m([1,3,-3])}, mpsm.monoid().filter(4))]), P.ideal(0), DegreeGrading((1,2)))
            sage: h = GradedExpansionVector_class(ger, a)
        """
        ModuleElement.__init__(self, parent)
        GradedExpansion_abstract.__init__(self, parent, polynomial)

    def __hash__(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_module import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_element import *
            sage: m = FreeModule(QQ, 3)
            sage: mpsm = MonoidPowerSeriesModule(m, NNMonoid(False))
            sage: mps = mpsm.base_ring()
            sage: P.<a,b> = QQ[]
            sage: ger = GradedExpansionModule_class(Sequence([MonoidPowerSeries(mps, {1: 1}, mps.monoid().filter(4))]), Sequence([MonoidPowerSeries(mpsm, {1: m([1,1,1]), 2: m([1,3,-3])}, mpsm.monoid().filter(4))]), P.ideal(0), DegreeGrading((1,2)))
            sage: h = GradedExpansionVector_class(ger, a)
            sage: hash(h)
            -1239234136 # 32-bit
            12415370528999848 # 64-bit
        """
        return hash(self.polynomial())
    
    def _repr_(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_module import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_element import *
            sage: m = FreeModule(QQ, 3)
            sage: mpsm = MonoidPowerSeriesModule(m, NNMonoid(False))
            sage: mps = mpsm.base_ring()
            sage: P.<a,b> = QQ[]
            sage: ger = GradedExpansionModule_class(Sequence([MonoidPowerSeries(mps, {1: 1}, mps.monoid().filter(4))]), Sequence([MonoidPowerSeries(mpsm, {1: m([1,1,1]), 2: m([1,3,-3])}, mpsm.monoid().filter(4))]), P.ideal(0), DegreeGrading((1,2)))
            sage: GradedExpansionVector_class(ger, a*b)
            Graded expansion vector (a,)
        """
        poly = self.polynomial()
        
        exps = poly.exponents()
        if poly.parent().ngens() == 1 :
            exps = map(lambda e: (e,), exps)
        exps = sorted(map(tuple, exps), reverse = True)
        
        coefficient_monomials = poly.parent().gens()[:self.parent().nbasegens()]
        exps = groupby(exps, lambda e: e[self.parent().nbasegens():])
        vector_repr = [poly.parent()(0) for _ in range(self.parent().ngens())]
        for v, es in exps :
            if v.count(1) != 1 :
                continue
            
            ind = v.index(1)
            c = 0
            for e in es :
                c += poly[e] * prod(map(operator.pow, coefficient_monomials, e[:self.parent().nbasegens()]))
            vector_repr[ind] += c
        
        return "Graded expansion vector %s" % (tuple(vector_repr),)
    
    def _latex_(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_module import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_element import *
            sage: m = FreeModule(QQ, 3)
            sage: mpsm = MonoidPowerSeriesModule(m, NNMonoid(False))
            sage: mps = mpsm.base_ring()
            sage: P.<a,b> = QQ[]
            sage: ger = GradedExpansionModule_class(Sequence([MonoidPowerSeries(mps, {1: 1}, mps.monoid().filter(4))]), Sequence([MonoidPowerSeries(mpsm, {1: m([1,1,1]), 2: m([1,3,-3])}, mpsm.monoid().filter(4))]), P.ideal(0), DegreeGrading((1,2)))
            sage: latex( GradedExpansionVector_class(ger, a*b) )
            \text{Graded expansion vector }\left(a\right)
        """
        poly = self.polynomial()
        
        exps = poly.exponents()
        if poly.parent().ngens() == 1 :
            exps = map(lambda e: (e,), exps)
        exps = sorted(map(tuple, exps), reverse = True)
        
        coefficient_monomials = poly.parent().gens()[:self.parent().nbasegens()]
        exps = groupby(exps, lambda e: e[self.parent().nbasegens():])
        vector_repr = [poly.parent()(0) for _ in range(self.parent().ngens())]
        for v, es in exps :
            if v.count(1) != 1 :
                continue
            
            ind = v.index(1)
            c = 0
            for e in es :
                c += poly[e] * prod(map(operator.pow, coefficient_monomials, e[:self.parent().nbasegens()]))
            vector_repr[ind] += c
                
        return r"\text{Graded expansion vector }%s" % (latex(tuple(vector_repr)),)
    
