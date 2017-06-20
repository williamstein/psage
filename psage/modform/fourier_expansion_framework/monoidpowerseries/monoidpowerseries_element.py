r"""
Monoid power series and equivariant monoid power series.

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

from copy import copy
from sage.structure.element import AlgebraElement
from sage.misc.latex import latex
from sage.misc.misc import union
from sage.modules.module import Module 
from sage.structure.element import ModuleElement
from sage.rings.ring import Ring

#===============================================================================
# MonoidPowerSeries
#===============================================================================

def MonoidPowerSeries( parent, coefficients, precision, cleanup_coefficients = False) :
    r"""
    Create a monoid power series within a given parent.
    
    INPUT:
        - ``parent``       -- A ring or module of monoid power series.
        - ``coefficients`` -- A dictionary with keys in the parent's monoid and values
                              in the parent coefficient domain.
        - ``precision``    -- A filter for the parent's monoid.
        - ``cleanup_coefficients`` -- A boolean (default: ``False``); If ``True`` zero
                                      coefficients will be erased from the dictionary.
    
    OUTPUT:
        An instance of :class:~`.MonoidPowerSeries_abstract`.
    
    EXAMPLES::
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
        sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
        sage: h = MonoidPowerSeries(mps, {1 : 1}, mps.monoid().filter(3))
        
    TESTS::
        sage: h = MonoidPowerSeries(mps, {1 : 1}, mps.monoid().zero_filter(), True)
        sage: h = MonoidPowerSeries(mps, {1 : 4, 0 : 3}, mps.monoid().filter_all())
        sage: mps = MonoidPowerSeriesModule(FreeModule(QQ, 2), NNMonoid(False))
        sage: h = MonoidPowerSeries(mps, {1 : 1}, mps.monoid().filter(3))
    """
    if isinstance(parent, Module) : 
        return MonoidPowerSeries_moduleelement(parent, coefficients, precision, cleanup_coefficients)
    if isinstance(parent, Ring) :
        return MonoidPowerSeries_algebraelement(parent, coefficients, precision, cleanup_coefficients)

    raise TypeError, "Unexpected type of parent"

#===============================================================================
# MonoidPowerSeries_abstract
#===============================================================================

class MonoidPowerSeries_abstract :
    r"""
    An element of the monoid power series ring or module up to
    given precision.
    """
    
    def __init__(self, parent, precision) :
        r"""
        INPUT:
            - ``parent``       -- An instance of :class:~`fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ambient.MonoidPowerSeriesAmbient_abstract`.
            - ``precision``    -- A filter associated to the parent's monoid.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: h = MonoidPowerSeries_abstract(mps, mps.monoid().zero_filter())
        """
        self.__precision = parent.monoid().filter(precision)
        
    def precision(self) :
        r"""
        The series' precision.
        
        OUTPUT:
            A filter for the parent's monoid.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: MonoidPowerSeries_abstract(mps, mps.monoid().zero_filter()).precision() == mps.monoid().zero_filter()
            True
        """
        return self.__precision

    def _set_precision(self, precision) :
        r"""
        Set the series' precision.
        
        INPUT:
            - ``precision`` -- A filter for the parent's monoid or a an object that can be converted
                               to a filter.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: h = copy(mps.one_element())
            sage: h._set_precision(mps.monoid().filter(2))
            sage: h.precision() == mps.monoid().filter(2)
            True
            sage: h._set_precision(3)
            sage: h.precision() == mps.monoid().filter(3)
            True
        """
        self.__precision = self.parent().monoid().filter(precision)

    def _bounding_precision(self) :
        r"""
        If ``self.precision()`` is an infinite filter,  return a filter
        which contains all non zero coefficients of this series. Otherwise,
        return ``self.precision()``
        
        OUTPUT:
            A filter for the parent's monoid.
            
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: h = MonoidPowerSeries(mps, dict(), mps.monoid().filter(2))
            sage: h._bounding_precision() == mps.monoid().filter(2)
            True
            sage: h = MonoidPowerSeries(mps, dict(), mps.monoid().filter_all())
            sage: h._bounding_precision() == mps.monoid().zero_filter()
            True
        """
        if not self.precision().is_infinite() :
            return self.precision()
         
        return self.parent().monoid().minimal_composition_filter( self.coefficients().keys(),
                                                                 [self.parent().monoid().zero_element()] )

    def coefficients(self) :
        r"""
        The coefficients of ``self``.
        
        OUTPUT:
            A dictionary with keys the elements of the parent's monoid and values in the
            parent's coefficient domain.

        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: MonoidPowerSeries_abstract(mps, mps.monoid().filter_all()).coefficients()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError
    
    def _truncate_in_place(self, precision) :
        r"""
        Truncate ``self`` modifying the coefficient dictionary directly.
        
        INPUT:
            - ``precision``    -- A filter for the parent's monoid or a an object that can be converted
                                  to a filter.

        OUTPUT:
            ``None``
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: MonoidPowerSeries_abstract(mps, mps.monoid().filter_all())._truncate_in_place(mps.monoid().zero_filter())
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError
    
    def truncate(self, precision) :
        r"""
        Truncate a copy of ``self``.

        INPUT:
            - ``precision``    -- A filter for the parent's monoid or a an object that can be converted
                                  to a filter.

        OUTPUT:
            An instance of :class:~`.MonoidPowerSeries_abstract_nonlazy`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: MonoidPowerSeries_abstract(mps, mps.monoid().filter_all()).truncate(mps.monoid().zero_filter())
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError
    
    def _add_(left, right) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: mps.one_element() == mps.one_element() + MonoidPowerSeries(mps, dict(), mps.monoid().filter_all())
            True
        """
        prec = min(left.__precision, right.__precision)
        lcoeffs = left.coefficients()
        rcoeffs = right.coefficients()
        
        lkeys = set(lcoeffs)
        rkeys = set(rcoeffs)

        d = dict()
        for k in lkeys - rkeys :
            d[k] = lcoeffs[k]
        for k in rkeys - lkeys :
            d[k] = rcoeffs[k]
        for k in lkeys.intersection(rkeys) :
            d[k] = lcoeffs[k] + rcoeffs[k]
                
        return MonoidPowerSeries(left.parent(), d, prec)
    
    def _mul_(left, right, switch_factors = False) :
        r"""
        NOTE:
            This function has to accept algebra and module elements and will
            also compute the action of the base ring (which might be a monoid power series)
            on a module.

        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: mps.one_element() * mps.one_element() == mps.one_element()
            True
        """
        mul_fc = left.parent()._multiply_function()
        
        if not switch_factors :
            lcoeffs = left.coefficients()
            rcoeffs = right.coefficients()
        else :
            lcoeffs = right.coefficients()
            rcoeffs = left.coefficients()
        
        prec = min(left.precision(), right.precision())
        if prec.is_infinite() :
            if len(lcoeffs) == 0 or len(rcoeffs) == 0:
                return MonoidPowerSeries(left.parent(), dict(), prec)
            iter_prec = left.parent().monoid(). \
                         minimal_composition_filter(set(lcoeffs), set(rcoeffs))
        else :
            iter_prec = prec
        
        d = dict()
        for k in iter_prec :
            v = mul_fc( k, lcoeffs, rcoeffs,
                           left.parent().coefficient_domain().zero_element() )
            if not v.is_zero() :
                d[k] = v
            
        return MonoidPowerSeries(left.parent(), d, prec)
    
    def _lmul_(self, c) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: (mps.one_element() * 2) * 3 == mps.one_element() * 6
            True
            sage: m = FreeModule(QQ, 3)
            sage: mpsm = MonoidPowerSeriesModule(m, NNMonoid(False))
            sage: h = MonoidPowerSeries(mps, {1 : 1, 3 : 2}, mps.monoid().filter(4))
            sage: hv = MonoidPowerSeries(mpsm, {1 : m([1,0,0]), 2 : m([0,1,0])}, mps.monoid().filter(4))
            sage: hh = hv * h
        """
        if c.is_zero() :
            return MonoidPowerSeries(self.parent(), dict(), None)
        
        if isinstance(c, MonoidPowerSeries_abstract) and not isinstance(self, AlgebraElement) :
            coeffs = c.coefficients()
            if len(coeffs) == 1 and self.parent().monoid().zero_element() in coeffs :
                c = coeffs[self.parent().monoid().zero_element()]
                d = dict((k, c*v) for (k,v) in self.coefficients().iteritems())
            
                return MonoidPowerSeries(self.parent(), d, self.precision())
        
            return self._mul_(c, False)
        else :
            d = dict((k, c*v) for (k,v) in self.coefficients().iteritems())
            
            return MonoidPowerSeries(self.parent(), d, self.precision())
    
    def _rmul_(self, c) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: 3 * (2 * mps.one_element())  == 6 * mps.one_element()
            True
            sage: m = FreeModule(QQ, 3)
            sage: mpsm = MonoidPowerSeriesModule(m, NNMonoid(False))
            sage: h = MonoidPowerSeries(mps, {1 : 1, 3 : 2}, mps.monoid().filter(4))
            sage: hv = MonoidPowerSeries(mpsm, {1 : m([1,0,0]), 2 : m([0,1,0])}, mps.monoid().filter(4))
            sage: hh = h * hv
        """

        if c.is_zero() :
            return MonoidPowerSeries(self.parent(), dict(), None)

        if isinstance(c, MonoidPowerSeries_abstract) and not isinstance(self, AlgebraElement) :
            coeffs = c.coefficients()
            if len(coeffs) == 1 and self.parent().monoid().zero_element() in coeffs :
                c = coeffs[self.parent().monoid().zero_element()]
                d = dict((k, v*c) for (k,v) in self.coefficients().iteritems())
            
                return MonoidPowerSeries(self.parent(), d, self.precision())
            
            return self._mul_(c, True)
        else :
            d = dict((k, v*c) for (k,v) in self.coefficients().iteritems())
            
            return MonoidPowerSeries(self.parent(), d, self.precision())

    def __contains__(self, k) :
        r"""
        Check whether `k` is below the series' precision.
        
        INPUT:
            - `k` -- An element of the parent's monoid.
        
        OUTPUT:
            A boolean.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: 2 in mps.one_element()
            True
            sage: 1 in mps.one_element().truncate(mps.monoid().zero_filter())
            False
        """
        return k in self.precision()

    def __getitem__(self, s) :
        r"""
        Return the `k`-th coefficient if it below the series' precision.
        
        INPUT:
            - `k` -- An element of the parent's monoid.
        
        OUTPUT:
            An element of the parent's coefficient domain.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: MonoidPowerSeries_abstract(mps, mps.monoid().filter_all())[0]
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError
        
    def __cmp__(self, other) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: mps.one_element() == mps.one_element()
            True
            sage: mps.one_element() == mps.zero_element()
            False
        """

        c = cmp(self.precision(), other.precision())
        if c != 0 : return c
        
        self_coeffs = self.coefficients()
        other_coeffs = other.coefficients() 
        self_keys = set(self_coeffs)
        other_keys = set(other_coeffs)

        for k in self_keys - other_keys :
            if not self_coeffs[k].is_zero() :
                return -1
        for k in other_keys - self_keys :
            if not other_coeffs[k].is_zero() :
                return -1
        for k in self_keys.intersection(other_keys) :
            if self_coeffs[k] != other_coeffs[k] :
                return -1

        return 0

    def _repr_(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: mps.one_element() # indirect doctest
            Monoid power series in Ring of monoid power series over NN
        """
        return "Monoid power series in %s" % self.parent()
    
    def _latex_(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: latex(mps.one_element()) # indirect doctest
            \text{Monoid power series in } \text{Ring of monoid power series over }\Bold{N}
        """
        return r"\text{Monoid power series in }" + latex(self.parent())

#===============================================================================
# MonoidPowerSeries_abstract_nonlazy
#===============================================================================

class MonoidPowerSeries_abstract_nonlazy (MonoidPowerSeries_abstract) :
    r"""
    A abstract implementation of monoid power series that store their coefficients.
    """
        
    def __init__(self, parent, coefficients, precision, cleanup_coefficients) :
        r"""
        INPUT:
            - ``parent``       -- An instance of :class:~`fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ambient.MonoidPowerSeriesAmbient_abstract`.
            - ``coefficients`` -- A dictionary with keys in the parent's monoid and values
                                  in the parent coefficient domain.
            - ``precision``    -- A filter associated to the parent's monoid.
            - ``cleanup_coefficients`` -- A boolean; If ``True`` zero coefficients will be
                                          erased from the dictionary.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: h = MonoidPowerSeries_abstract_nonlazy(mps, dict(), mps.monoid().zero_filter(), False)
        """
        MonoidPowerSeries_abstract.__init__(self, parent, precision)

        if cleanup_coefficients and len(coefficients) != 0 :
            coefficients = self._cleanup_coefficients( coefficients, in_place = True )
            
        self.__coefficients = coefficients
    
    def coefficients(self) :
        r"""
        The coefficients of ``self``.
        
        OUTPUT:
            A dictionary with keys the elements of the parent's monoid and values in the
            parent's coefficient domain. 

        NOTE:
            Some keys may be invalid. To get an exact result
            call ``_cleanup_coefficients(in_place = True)`` before.

        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: mps.one_element().coefficients()
            {0: 1}
        """
        return self.__coefficients
    
    def _cleanup_coefficients(self, coefficients = None, in_place = True) :
        r"""
        Remove zero entries and entries not below ``self.precision()`` from a coefficient dictionary.
        
        INPUT:
            - ``coefficients`` -- ``None`` or a dictionary (default: ``None``); If ``None`` the
                                  coefficient dictionary assigned to ``self`` will be cleaned.
            - ``in_place``     -- A boolean (default: ``True``); If ``False`` a copy of the coefficient
                                  dictionary will me cleaned and returned.
        
        OUTPUT:
            A dictionary with keys the elements of the parent's monoid and values in the
            parent's coefficient domain. 

        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: mps.one_element()._cleanup_coefficients()
            {0: 1}
            sage: d = {1 : 1}
            sage: tmp = MonoidPowerSeries(mps, dict(), mps.monoid().zero_filter())._cleanup_coefficients(d)
            sage: d
            {}
            sage: mps.zero_element()._cleanup_coefficients({1 : 0}, False)
            {}
            sage: h = copy(mps.one_element())
            sage: h._set_precision(mps.monoid().zero_filter())
            sage: h._cleanup_coefficients(in_place = False)
            {}
        """
        if coefficients is None :
            coefficients = self.__coefficients
        
        if in_place :
            for s in coefficients.keys() :
                if not s in self.precision() or coefficients[s].is_zero() :
                    del coefficients[s]
        else :
            ncoefficients = dict()
            
            for s in coefficients :
                if not s in self.precision() : continue

                v = coefficients[s]
                if v.is_zero() : continue

                ncoefficients[s] = v

        if in_place :
            return coefficients
        else :
            return ncoefficients

    def _truncate_in_place(self, precision) :
        r"""
        Truncate ``self`` modifying the coefficient dictionary directly.
        
        INPUT:
            - ``precision``    -- A filter for the parent's monoid or a an object that can be converted
                                  to a filter.

        OUTPUT:
            ``None``
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: h = copy(mps.one_element())
            sage: h._truncate_in_place(mps.monoid().zero_filter())
            sage: h.coefficients()
            {}
            sage: h = copy(mps.one_element())
            sage: h._truncate_in_place(0)
            sage: h.coefficients()
            {}
        """
        precision = self.parent().monoid().filter(precision)
        nprec = min(self.precision(), precision)

        if nprec != self.precision() :
            coefficients = self.__coefficients
            for k in coefficients.keys() :
                if not k in nprec :
                    del coefficients[k]
            
        self._set_precision(nprec)
    
    def truncate(self, precision) :
        r"""
        Truncate a copy of ``self``.
        
        INPUT:
            - ``precision``    -- A filter for the parent's monoid or a an object that can be converted
                                  to a filter.

        OUTPUT:
            An instance of :class:~`.MonoidPowerSeries_abstract_nonlazy`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: mps.one_element().truncate(mps.monoid().zero_filter()).coefficients()
            {}
            sage: mps.one_element().truncate(0).coefficients()
            {}
        """
        precision = self.parent().monoid().filter(precision)
        nprec = min(self.precision(), precision)

        ncoefficients = copy(self.__coefficients)
        return MonoidPowerSeries( self.parent(), ncoefficients, nprec, cleanup_coefficients = True )
    
    def __getitem__(self, s) :
        r"""
        Return the `k`-th coefficient if it below the series' precision.
        
        INPUT:
            - `k` -- An element of the parent's monoid.
        
        OUTPUT:
            An element of the parent's coefficient domain.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: MonoidPowerSeries_abstract_nonlazy(mps, { 0 : 10, 1 : 1 }, mps.monoid().filter_all(), False)[0]
            10
        """
        try :
            return self.coefficients()[s]
        except KeyError :
            return self.parent().coefficient_domain().zero_element()

#===============================================================================
# MonoidPowerSeries_moduleelement
#===============================================================================

class MonoidPowerSeries_moduleelement ( MonoidPowerSeries_abstract_nonlazy, ModuleElement ) :
    r"""
    An element of a module of monoid power series.
    """

    def __init__(self, parent, coefficients, precision, cleanup_coefficients) :
        r"""
        INPUT:
            - ``parent``       -- An instance of :class:~`fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ambient.MonoidPowerSeriesAmbient_abstract`.
            - ``coefficients`` -- A dictionary with keys in the parent's monoid and values
                                  in the parent coefficient domain.
            - ``precision``    -- A filter associated to the parent's monoid.
            - ``cleanup_coefficients`` -- A boolean; If ``True`` zero coefficients will be
                                          erased from the dictionary.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: mps = MonoidPowerSeriesModule(FreeModule(QQ, 2), NNMonoid(False))
            sage: h = MonoidPowerSeries_moduleelement(mps, dict(), mps.monoid().zero_filter(), False)
        """
        ModuleElement.__init__(self, parent)
        MonoidPowerSeries_abstract_nonlazy.__init__(self, parent, coefficients, precision, cleanup_coefficients)        

#===============================================================================
# MonoidPowerSeries_algebraelement
#===============================================================================

class MonoidPowerSeries_algebraelement ( MonoidPowerSeries_abstract_nonlazy, AlgebraElement ) :
    r"""
    An element of a algebra of monoid power series.
    """
    
    def __init__(self, parent, coefficients, precision, cleanup_coefficients) :
        r"""
        INPUT:
            - ``parent``       -- An instance of :class:~`fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ambient.MonoidPowerSeriesAmbient_abstract`.
            - ``coefficients`` -- A dictionary with keys in the parent's monoid and values
                                  in the parent coefficient domain.
            - ``precision``    -- A filter associated to the parent's monoid.
            - ``cleanup_coefficients`` -- A boolean; If ``True`` zero coefficients will be
                                          erased from the dictionary.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: h = MonoidPowerSeries_algebraelement(mps, dict(), mps.monoid().zero_filter(), False)
        """
        AlgebraElement.__init__(self, parent)
        MonoidPowerSeries_abstract_nonlazy.__init__(self, parent, coefficients, precision, cleanup_coefficients)

#===============================================================================
# EquivariantMonoidPowerSeries
#===============================================================================

def EquivariantMonoidPowerSeries( parent, coefficients, precision, symmetrise = False,
                                  cleanup_coefficients = False) :
    r"""
    Create an equivariant monoid power series within a given parent.
    
    INPUT:
        - ``parent``       -- A ring or module of equivariant monoid power series.
        - ``coefficients`` -- A dictionary with keys in the parent's monoid and values
                              in the parent coefficient domain.
        - ``precision``    -- A filter for the parent's action.
        - ``symmetrise``   -- A boolean (default: ``False``); If ``True`` every enty in
                              ``coefficients`` will contribute to its whole orbit.
        - ``cleanup_coefficients`` -- A boolean (default: ``False``); If ``True`` zero
                                      coefficients will be erased from the dictionary.
    
    OUTPUT:
        An instance of :class:~`.EquivariantMonoidPowerSeries_abstract`.
    
    EXAMPLES::
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
        sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
        sage: h = EquivariantMonoidPowerSeries(emps, {emps.characters().one_element() : {1 : 1}}, emps.action().filter(3))
        
    TESTS::
        sage: h = EquivariantMonoidPowerSeries(emps, {emps.characters().one_element() : {1 : 1}}, emps.action().zero_filter(), True)
        sage: h = EquivariantMonoidPowerSeries(emps, {emps.characters().one_element() : {1 : 4, 0 : 3}}, emps.action().filter_all())
        sage: h = EquivariantMonoidPowerSeries(emps, {emps.characters().one_element() : {1 : 4, 0 : 3}}, emps.action().filter_all(), symmetrise = True)
        sage: emps = EquivariantMonoidPowerSeriesModule( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", FreeModule(QQ, 2)) )
        sage: h = EquivariantMonoidPowerSeries(emps, {emps.characters().one_element() : {1 : 1}}, emps.action().filter(3))
    """
    if isinstance(parent, Module) : 
        return EquivariantMonoidPowerSeries_moduleelement( parent, coefficients, precision, symmetrise,
                                                           cleanup_coefficients )
    if isinstance(parent, Ring) :
        return EquivariantMonoidPowerSeries_algebraelement( parent, coefficients, precision, symmetrise,
                                                            cleanup_coefficients )

    raise TypeError, "Unexpected type of parent"
    
#===============================================================================
# EquivariantMonoidPowerSeries_abstract
#===============================================================================

class EquivariantMonoidPowerSeries_abstract :
    r"""
    An abstract element of an equivariant monoid power series ring up to
    given precision.
    """
    
    def __init__( self, parent, precision ) :
        r"""
        INPUT:
            - ``parent``       -- A ring or module of equivariant monoid power series.
            - ``precision``    -- A filter for the parent's action.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
            sage: h = EquivariantMonoidPowerSeries_abstract(emps, emps.action().zero_filter()) # indirect doctest
        """
        self.__precision = parent.action().filter(precision)

    def precision(self) :
        r"""
        The series' precision.
        
        OUTPUT:
            A filter for the parent's action.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
            sage: EquivariantMonoidPowerSeries_abstract(emps, emps.action().filter(3)).precision() == emps.action().filter(3)
            True
        """
        return self.__precision
    
    def _set_precision(self, precision) :
        r"""
        Set the series' precision.
    
        INPUT:
            - ``precision`` -- A filter for the parent's action.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
            sage: e = copy(emps.one_element())
            sage: e._set_precision(emps.action().filter(3))
            sage: e.precision() == emps.action().filter(3)
            True
            sage: e._set_precision(2)
            sage: e.precision() == emps.action().filter(2)
            True
        """
        self.__precision = self.parent().action().filter(precision)

    def non_zero_components(self) :
        r"""
        Return all those characters  which are not guaranteed to have only
        vanishing coefficients associated to.
        
        OUTPUT:
            A list of elements of the character monoid.
        
        NOTE:
            The components associated to characters this function returns can vanish. 
            For exact results use ``_cleanup_coefficients(in_place = True)`` before.
            
        EXAMPLES::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
            sage: e = EquivariantMonoidPowerSeries(emps, dict(), emps.monoid().zero_filter())
            sage: e.non_zero_components()
            []
            sage: emps.one_element().non_zero_components()
            [1] 
        """
        return list(self.parent().characters())
    
    def _bounding_precision(self) :
        r"""
        If ``self.precision()`` is an infinite filter,  return a filter
        which contains all non zero coefficients of this series. Otherwise,
        return ``self.precision()``
        
        OUTPUT:
            A filter for the parent's action.
            
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
            sage: e = EquivariantMonoidPowerSeries(emps, dict(), emps.action().filter(2))
            sage: e._bounding_precision() == emps.action().filter(2)
            True
            sage: e = EquivariantMonoidPowerSeries(emps, dict(), emps.action().filter_all())
            sage: e._bounding_precision() == emps.action().zero_filter()
            True
        """
        if not self.precision().is_infinite() :
            return self.precision()
         
        coeffs = self.coefficients(True)
        m = self.parent().action().zero_filter()
        for c in self.non_zero_components() :
            m = max(m, self.parent().action().minimal_composition_filter( coeffs[c].keys(),
                                                                          [self.parent().action().zero_element()] ))
        return m
    
    def coefficients(self, force_characters = False) :
        r"""
        The coefficients of ``self``. 

        INPUT:
            - ``force_characters`` -- A boolean (default: ``False``); If ``True`` the
                                      the dictionary returned will have characters as keys
                                      in any cases.
        
        OUTPUT:
            Either of the following two:
                - A dictionary with keys the elements of the parent's monoid and values in the
                  parent's coefficient domain.
                - A dictionary with keys the parent's characters and values the a dictionary as follows. This
                  dictionary has keys the elements of the parent's monoid and values in the parent's
                  coefficient domain. 
        
        EXAMPLES::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
            sage: EquivariantMonoidPowerSeries_abstract(emps, emps.action().zero_filter()).coefficients()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError
        
    def _truncate_in_place(self, precision) :
        r"""
        Truncate ``self`` modifying the coefficient dictionary directly.
        
        INPUT:
            - ``precision``    -- A filter for the parent's monoid or a an object that can be converted
                                  to a filter.

        OUTPUT:
            ``None``
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
            sage: e = emps.one_element()
            sage: EquivariantMonoidPowerSeries_abstract(emps, emps.action().filter_all())._truncate_in_place(emps.monoid().zero_filter())
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError
        
    def truncate(self, precision) :
        r"""
        Truncate a copy of ``self``.

        INPUT:
            - ``precision``    -- A filter for the parent's monoid or a an object that can be converted
                                  to a filter.

        OUTPUT:
            An instance of :class:~`.EquivariantMonoidPowerSeries_abstract`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
            sage: EquivariantMonoidPowerSeries_abstract(emps, emps.action().filter_all()).truncate(emps.action().zero_filter()) 
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError
    
    def _add_(left, right) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
            sage: e = EquivariantMonoidPowerSeries(emps, {emps.characters().one_element() : { 1 : 1 }}, emps.action().filter_all())
            sage: (e + e)[1]
            2
        """
        prec = min(left.precision(), right.precision())
        
        left_coefficients = left.coefficients(True)
        right_coefficients = right.coefficients(True)
        
        left_characters = set(left_coefficients)
        right_characters = set(right_coefficients)
        coefficients = dict()
        
        for ch in left_characters - right_characters :
            coefficients[ch] = copy(left_coefficients[ch])
        for ch in right_characters - left_characters :
            coefficients[ch] = copy(right_coefficients[ch])
            
        for ch in left_characters.intersection(right_characters) :
            lcoeffs = left_coefficients[ch]
            rcoeffs = right_coefficients[ch]

            lcoeff_keys = set(lcoeffs)
            rcoeff_keys = set(rcoeffs)
            
            nd = dict()
            for k in lcoeff_keys - rcoeff_keys :
                nd[k] = lcoeffs[k]
            for k in rcoeff_keys - lcoeff_keys :
                nd[k] = rcoeffs[k]
            for k in lcoeff_keys.intersection(rcoeff_keys) :
                nd[k] = lcoeffs[k] + rcoeffs[k]
                
            coefficients[ch] = nd
                                    
        return EquivariantMonoidPowerSeries(left.parent(),
                coefficients, prec)
    
    def _mul_(left, right, switch_factors = False) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
            sage: e = EquivariantMonoidPowerSeries(emps, {emps.characters().one_element() : { 1 : 1 }}, emps.action().filter_all())
            sage: (e * e)[2]
            1
        """            
        mul_fc = left.parent()._multiply_function()
        coefficient_domain = left.parent().coefficient_domain()

        prec = min(left.precision(), right.precision())
        if not switch_factors :
            left_coefficients = left.coefficients(True)
            right_coefficients = right.coefficients(True)
        else :
            right_coefficients = left.coefficients(True)
            left_coefficients = right.coefficients(True)

        if prec.is_infinite() :
            left_keys  = reduce(union, (set(left_coefficients[c]) for c in left_coefficients), set())
            right_keys = reduce(union, (set(right_coefficients[c]) for c in right_coefficients), set())
            
            if len(left_keys) == 0 or len(right_keys) == 0:
                return EquivariantMonoidPowerSeries(left.parent(), dict(), prec)
            iter_prec = left.parent().action(). \
                         minimal_composition_filter(left_keys, right_keys)
        else :
            iter_prec = prec
                
        coefficients = dict()
        for c1 in left_coefficients :
            lcoeffs = left_coefficients[c1]
            if len(lcoeffs) == 0 : continue
            
            for c2 in right_coefficients :
                rcoeffs = right_coefficients[c2]
                if len(rcoeffs) == 0 : continue

                try :
                    d = coefficients[c1 * c2]
                except KeyError :
                    d = dict()
                    coefficients[c1 * c2] = d
                    
                for k in iter_prec :
                    v = mul_fc( k, lcoeffs, rcoeffs, c1, c2, coefficient_domain(0) )
                        
                    if not v.is_zero() :
                        try :
                            d[k] += v
                            if d[k].is_zero() :
                                del d[k]
                        except KeyError :
                            d[k] = v

        return EquivariantMonoidPowerSeries( left.parent(), coefficients, prec )

    def _lmul_(self, c) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
            sage: (emps.one_element() * 2) * 3 == emps.one_element() * 6
            True
            sage: m = FreeModule(QQ, 3)
            sage: empsm = EquivariantMonoidPowerSeriesModule(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", m))
            sage: h = EquivariantMonoidPowerSeries(emps, {emps.characters().one_element(): {1: 1, 3: 2}}, emps.action().filter_all())
            sage: hv = EquivariantMonoidPowerSeries(empsm, {empsm.characters().one_element(): {1: m([1,0,0]), 2: m([0,0,1])}}, empsm.action().filter_all())
            sage: hh = hv * h
        """
        if c.is_zero() :
            return EquivariantMonoidPowerSeries(self.parent(), dict(), self.parent().monoid().filter_all())

        if isinstance(c, EquivariantMonoidPowerSeries_abstract) and not isinstance(self, AlgebraElement) :
            nzc = c.non_zero_components()
            if len(nzc) == 1 :
                coeffs = c.coefficients(True)[nzc[0]]
    
                if len(coeffs) == 1 and self.parent().action().zero_element() in coeffs :
                    c = coeffs[self.parent().action().zero_element()]
                    
                    self_coefficients = self.coefficients(True)
                    coefficients = dict()
                    for ch in self_coefficients :
                        coefficients[ch] = dict((k, c*v) for (k,v) in self_coefficients[ch].iteritems())
                    
                    return EquivariantMonoidPowerSeries(self.parent(),
                            coefficients, self.precision())
    
            return self._mul_(c, False)
        else :
            self_coefficients = self.coefficients(True)
            coefficients = dict()
            for ch in self_coefficients :
                coefficients[ch] = dict((k, c*v) for (k,v) in self_coefficients[ch].iteritems())
            
            return EquivariantMonoidPowerSeries(self.parent(),
                                                coefficients, self.precision())
            
    def _rmul_(self, c) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
            sage: 3 * (2 * emps.one_element())  == 6 * emps.one_element()
            True
            sage: m = FreeModule(QQ, 3)
            sage: empsm = EquivariantMonoidPowerSeriesModule(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", m))
            sage: h = EquivariantMonoidPowerSeries(emps, {emps.characters().one_element(): {1: 1, 3: 2}}, emps.action().filter_all())
            sage: hv = EquivariantMonoidPowerSeries(empsm, {empsm.characters().one_element(): {1: m([1,0,0]), 2: m([0,0,1])}}, empsm.action().filter_all())
            sage: hh = h * hv
        """
        if c.is_zero() :
            return EquivariantMonoidPowerSeries(self.parent(), dict(), self.parent().monoid().filter_all())

        if isinstance(c, EquivariantMonoidPowerSeries_abstract) and not isinstance(self, AlgebraElement) :
            nzc = c.non_zero_components()
            if len(nzc) == 1 :
                coeffs = c.coefficients(True)[nzc[0]]
    
                if len(coeffs) == 1 and self.parent().action().zero_element() in coeffs :
                    c = coeffs[self.parent().action().zero_element()]
                    
                    self_coefficients = self.coefficients(True)
                    coefficients = dict()
                    for ch in self_coefficients :
                        coefficients[ch] = dict((k, v*c) for (k,v) in self_coefficients[ch].iteritems())
                    
                    return EquivariantMonoidPowerSeries(self.parent(),
                            coefficients, self.precision())
    
            return self._mul_(c, True)
        else :
            self_coefficients = self.coefficients(True)
            coefficients = dict()
            for ch in self_coefficients :
                coefficients[ch] = dict((k, v*c) for (k,v) in self_coefficients[ch].iteritems())
                    
            return EquivariantMonoidPowerSeries(self.parent(),
                                                coefficients, self.precision())

    def __contains__(self, k) :
        r"""
        Check whether an index or a pair of character and index
        is containted in the precision.
        
        EXAMPLES:
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
            sage: e = EquivariantMonoidPowerSeries(emps, {emps.characters().one_element() : {1 : 1, 2 : 1}}, emps.action().filter(3))
            sage: 3 in e
            False
            sage: 2 in e
            True
            sage: (emps.characters().one_element(),2) in e
            True 
        """
        try :
            (ch, k) = k
            if k not in self.parent().monoid() :
                k = (ch, k)
        except TypeError :
            pass
            
        return k in self.precision()

    def __getitem__(self, k) :
        r"""
        Return the `k`-th coefficient if it below the series' precision. If no character is contained
        in the key ``self`` must have only one nonvanishing component.
        
        INPUT:
            - `k` -- A pair of an element of the parent's character monoid and
                     and element of the parent's monoid or an element of the parent's monoid.
        
        OUTPUT:
            An element of the parent's coefficient domain.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
            sage: EquivariantMonoidPowerSeries_abstract(emps, emps.action().filter_all())[(emps.characters().one_element(), 0)]
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def __cmp__(self, other) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
            sage: e = EquivariantMonoidPowerSeries(emps, {emps.characters().one_element() : { 1 : 1 }}, emps.action().zero_filter())
            sage: e == EquivariantMonoidPowerSeries(emps, {emps.characters().one_element() : { 1 : 1 }}, emps.action().zero_filter())
            True
            sage: e == 2 * e
            False
            sage: e == EquivariantMonoidPowerSeries(emps, {}, emps.action().zero_filter())
            False
        """
        c = cmp(self.precision(), other.precision())
        if c != 0 : return c

        self_coefficients = self.coefficients(True)
        other_coefficients = other.coefficients(True)
        for ch in set(self_coefficients) - set(other_coefficients) :
            d = self_coefficients[ch]
            for k in d:
                if not d[k] == 0 :
                    return -1
        for ch in set(other_coefficients) - set(self_coefficients) :
            d = other_coefficients[ch]
            for k in d:
                if not d[k] == 0 :
                    return -1
                
        for ch in set(self_coefficients).intersection(set(other_coefficients)) :
            s = self_coefficients[ch]
            o = other_coefficients[ch]
            self_keys = set(s)
            other_keys = set(o)

            for k in self_keys - other_keys :
                if not s[k] == 0 :
                    return -1
            for k in other_keys - self_keys :
                if not o[k] == 0 :
                    return -1
            
            for k in self_keys.intersection(other_keys) :
                if s[k] != o[k] :
                    return -1
    
        return 0

    def _repr_(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
            sage: EquivariantMonoidPowerSeries(emps, {emps.characters().one_element() : { 1 : 1 }}, emps.action().zero_filter())
            Equivariant monoid power series in Ring of equivariant monoid power series over NN
        """
        return "Equivariant monoid power series in %s" % (self.parent(),)
    
    def _latex_(self) :
        r"""
        EXAMPLES:
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
            sage: latex( EquivariantMonoidPowerSeries(emps, {emps.characters().one_element() : { 1 : 1 }}, emps.action().zero_filter()) )
            \text{Equivariant monoid power series in }\text{Ring of equivariant monoid power series over }\Bold{N}
        """
        return r"\text{Equivariant monoid power series in }%s" % latex(self.parent())

#===============================================================================
# EquivariantMonoidPowerSeries_abstract_nonlazy
#===============================================================================

class EquivariantMonoidPowerSeries_abstract_nonlazy ( EquivariantMonoidPowerSeries_abstract ) :
    r"""
    A abstract implementation of equiavariant monoid power series that store their coefficients.
    """
    
    def __init__( self, parent, coefficients, precision, symmetrise,
                  cleanup_coefficients ) :
        r"""
        INPUT:
            - ``parent``       -- A ring or module of equivariant monoid power series.
            - ``coefficients`` -- A dictionary with keys in the parent's monoid and values
                                  in the parent coefficient domain.
            - ``precision``    -- A filter for the parent's action.
            - ``symmetrise``   -- A boolean (default: ``False``); If ``True`` every enty in
                                  ``coefficients`` will contribute to its whole orbit.
            - ``cleanup_coefficients`` -- A boolean (default: ``False``); If ``True`` zero
                                          coefficients will be erased from the dictionary.
                
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
            sage: h = EquivariantMonoidPowerSeries_abstract_nonlazy(emps, {emps.characters().one_element() : {1 : 1}}, emps.action().filter(3), False, False)
            sage: h = EquivariantMonoidPowerSeries_abstract_nonlazy(emps, {emps.characters().one_element() : {1 : 1}}, emps.action().zero_filter(), False, True)
            sage: h = EquivariantMonoidPowerSeries_abstract_nonlazy(emps, {emps.characters().one_element() : {1 : 4, 0 : 3}}, emps.action().filter_all(), False, False,)
            sage: h = EquivariantMonoidPowerSeries_abstract_nonlazy(emps, {emps.characters().one_element() : {1 : 4, 0 : 3}}, emps.action().filter_all(), True, False)
            sage: emps = EquivariantMonoidPowerSeriesModule( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", FreeModule(QQ, 2)) )
            sage: h = EquivariantMonoidPowerSeries_abstract_nonlazy(emps, {emps.characters().one_element() : {1 : 1}}, emps.action().filter(3), False, False)
        """
        EquivariantMonoidPowerSeries_abstract.__init__(self, parent, precision)
        
        if cleanup_coefficients and len(coefficients) != 0 :
            ncoefficients = self._cleanup_coefficients( coefficients,
                                  in_place = True )
        else :
            ncoefficients = coefficients

        if symmetrise :
            ## for symmetrisation we are guaranteed that
            ## the representation acts trivially on all coefficients
            reduction = parent._reduction_function()
            character_eval = parent._character_eval_function()
            
            self.__coefficients = dict()
            
            for ch in ncoefficients :
                d = ncoefficients[ch]
                nd = dict()
                self.__coefficients[ch] = nd
                
                for s in d :
                    rs, g = reduction(s)
                    try :
                        nd[rs] += character_eval(g, ch) * d[s]
                    except KeyError :
                        nd[rs] = character_eval(g, ch) * d[s]
        else :
            self.__coefficients = ncoefficients

    def non_zero_components(self) :
        r"""
        Return all those characters  which are not guaranteed to have only
        vanishing coefficients associated to.
        
        OUTPUT:
            A list of elements of the character monoid.
        
        NOTE:
            The components associated to characters this function returns can vanish. 
            For exact results use ``_cleanup_coefficients(in_place = True)`` before.
            
        EXAMPLES::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
            sage: e = EquivariantMonoidPowerSeries_abstract_nonlazy(emps, {emps.characters().one_element() : {1 : 0}}, emps.action().zero_filter(), False, False)
            sage: e.non_zero_components()
            [1]
            sage: e._cleanup_coefficients(in_place = True)
            {}
            sage: e.non_zero_components()
            []
        """
        return self.__coefficients.keys()
    
    def coefficients(self, force_characters = False) :
        r"""
        The coefficients of ``self``. 

        INPUT:
            - ``force_characters`` -- A boolean (default: ``False``); If ``True`` the
                                      the dictionary returned will have characters as keys
                                      in any cases.
        
        OUTPUT:
            Either of the following two:
                - A dictionary with keys the elements of the parent's monoid and values in the
                  parent's coefficient domain.
                - A dictionary with keys the parent's characters and values the a dictionary as follows. This
                  dictionary has keys the elements of the parent's monoid and values in the parent's
                  coefficient domain. 

        NOTE:
            Some keys may be invalid. To get an exact result
            call ``_cleanup_coefficients(in_place = True)`` before.
        
        EXAMPLES::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
            sage: EquivariantMonoidPowerSeries_abstract_nonlazy(emps, {emps.characters().one_element() : {1 : 1}}, emps.action().zero_filter(), False, False).coefficients()
            {1: 1}
            sage: EquivariantMonoidPowerSeries_abstract_nonlazy(emps, {emps.characters().one_element() : {1 : 1}}, emps.action().zero_filter(), False, False).coefficients(True)
            {1: {1: 1}}
        """
        if len(self.__coefficients) == 0 :
            return dict()
        elif not force_characters and len(self.__coefficients) == 1 :
            return self.__coefficients.values()[0] 
        else :
            return self.__coefficients

    def _cleanup_coefficients(self, coefficients = None, in_place = True) :
        r"""
        Remove zero entries and entries not below ``self.precision()`` from a coefficient dictionary.
        
        INPUT:
            - ``coefficients`` -- ``None`` or a dictionary (default: ``None``); If ``None`` the
                                  coefficient dictionary assigned to ``self`` will be cleaned.
            - ``in_place``     -- A boolean (default: ``True``); If ``False`` a copy of the coefficient
                                  dictionary will me cleaned and returned.
        
        OUTPUT:
            A dictionary with keys the parent's characters and values the a dictionary as follows. This
            dictionary has keys the elements of the parent's monoid and values in the parent's
            coefficient domain. 

        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
            sage: emps.one_element()._cleanup_coefficients()
            {1: {0: 1}}
            sage: d = {emps.characters().one_element() : {1 : 1}}
            sage: tmp = EquivariantMonoidPowerSeries(emps, dict(), emps.action().zero_filter())._cleanup_coefficients(d)
            sage: d
            {}
            sage: emps.zero_element()._cleanup_coefficients({emps.characters().one_element() : {1 : 0}}, False)
            {}
            sage: h = copy(emps.one_element())
            sage: h._set_precision(emps.action().zero_filter())
            sage: h._cleanup_coefficients(in_place = False)
            {}
        """
        if coefficients is None :
            coefficients = self.__coefficients
        
        if not in_place :
            ncoefficients = dict()
            
        for ch in coefficients.keys() :
            d = coefficients[ch]
            
            if in_place :
                for s in d.keys() :
                    if not s in self.precision() or d[s].is_zero() :
                        del d[s]
                        
                if len(d) == 0 :
                    del coefficients[ch]
            else :
                nd = dict()
                
                for s in d :
                    if not s in self.precision() : continue

                    v = d[s]
                    if v.is_zero() : continue

                    nd[s] = v
                
                if len(nd) != 0 :
                    ncoefficients[ch] = nd

        if in_place :
            return coefficients
        else :
            return ncoefficients

    def _truncate_in_place(self, precision) :
        r"""
        Truncate ``self`` modifying the coefficient dictionary directly.

        INPUT:
            - ``precision``    -- A filter for the parent's monoid or a an object that can be converted
                                  to a filter.

        OUTPUT:
            ``None``
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
            sage: e = copy(emps.one_element())
            sage: e._truncate_in_place(emps.monoid().zero_filter())
            sage: e.coefficients()
            {}
            sage: e = EquivariantMonoidPowerSeries(emps, {emps.characters().one_element() : {1 : 1, 2 : 0}}, emps.action().filter(3))
            sage: e._truncate_in_place(2)
            sage: 2 in e.coefficients()
            False
        """
        precision = self.parent().action().filter(precision)
        nprec = min(self.precision(), precision)

        if nprec != self.precision() :
            for c in self.__coefficients :
                d = self.__coefficients[c]
                for k in d.keys() :
                    if not k in nprec :
                        del d[k]
            
            self._set_precision(nprec)
    
    def truncate(self, precision) :
        r"""
        Truncate a copy of ``self``.

        INPUT:
            - ``precision``    -- A filter for the parent's monoid or a an object that can be converted
                                  to a filter.

        OUTPUT:
            An instance of :class:~`.EquivariantMonoidPowerSeries_abstract`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
            sage: emps.one_element().truncate(emps.monoid().zero_filter()).coefficients()
            {}
            sage: e = EquivariantMonoidPowerSeries(emps, {emps.characters().one_element() : {1 : 1, 2 : 0}}, emps.action().filter(3))
            sage: 2 in e.truncate(2)
            False
        """
        precision = self.parent().action().filter(precision)
        nprec = min(self.precision(), precision)

        ncoefficients = dict( (ch, copy(self.__coefficients[ch]))
                               for ch in self.__coefficients )
        return EquivariantMonoidPowerSeries( self.parent(),
                 ncoefficients, nprec, cleanup_coefficients = True )

    def __getitem__(self, k) :
        r"""
        Return the `k`-th coefficient if it below the series' precision. If no character is contained
        in the key ``self`` must have only one nonvanishing component.
        
        INPUT:
            - `k` -- A pair of an element of the parent's character monoid and
                     and element of the parent's monoid or an element of the parent's monoid.
        
        OUTPUT:
            An element of the parent's coefficient domain.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
            sage: EquivariantMonoidPowerSeries(emps, {emps.characters().one_element() : { 0 : 10, 1 : 1 }}, emps.action().filter_all())[(emps.characters().one_element(), 0)]
            10
            sage: EquivariantMonoidPowerSeries(emps, {emps.characters().one_element() : { 0 : 10, 1 : 1 }}, emps.action().filter_all())[1]
            1
        """
        try :
            if not isinstance(k, tuple) :
                raise ValueError
            
            (ch, k) = k
            if k not in self.parent().monoid() :
                s = (ch, k)
                ch = None
            else :
                s = k
        except ValueError :
            s = k
            ch = None
        
        try :
            if not ch.parent() == self.parent().characters() :
                ch = None
        except AttributeError :
            ch = None
            
        if ch is None :
            ns = self.non_zero_components()
            if len(ns) == 0 :
                return 0
            elif len(ns) == 1 :
                ch = ns[0]
            else :
                raise ValueError, "you must specify a character"
        
        if not s in self.precision() :
            raise ValueError, "%s out of bound" % (s,)

        try :
            return self.__coefficients[ch][s]
        except KeyError :
            (rs, g) = self.parent()._reduction_function()(s)
            
            try :
                return self.parent()._character_eval_function()(g, ch) \
                 * self.parent()._apply_function()(g, self.__coefficients[ch][rs])
            except KeyError :
                return self.parent().coefficient_domain().zero_element()

#===============================================================================
# EquivariantMonoidPowerSeries_moduleelement
#===============================================================================

class EquivariantMonoidPowerSeries_moduleelement ( EquivariantMonoidPowerSeries_abstract_nonlazy, ModuleElement ) :
    r"""
    An element of a module of equivariant monoid power series.
    """

    def __init__(self, parent, coefficients, precision, symmetrise,
                 cleanup_coefficients ) :
        r"""
        INPUT:
            - ``parent``       -- A ring or module of equivariant monoid power series.
            - ``coefficients`` -- A dictionary with keys in the parent's monoid and values
                                  in the parent coefficient domain.
            - ``precision``    -- A filter for the parent's action.
            - ``symmetrise``   -- A boolean (default: ``False``); If ``True`` every enty in
                                  ``coefficients`` will contribute to its whole orbit.
            - ``cleanup_coefficients`` -- A boolean (default: ``False``); If ``True`` zero
                                          coefficients will be erased from the dictionary.

        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesModule( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", FreeModule(QQ, 2)) )
            sage: h = EquivariantMonoidPowerSeries_moduleelement(emps, dict(), emps.action().zero_filter(), False, False)
        """
        ModuleElement.__init__(self, parent)
        EquivariantMonoidPowerSeries_abstract_nonlazy.__init__(self, parent, coefficients, precision, symmetrise,
                                            cleanup_coefficients)

#===============================================================================
# EquivariantMonoidPowerSeries_algebraelement
#===============================================================================

class EquivariantMonoidPowerSeries_algebraelement ( EquivariantMonoidPowerSeries_abstract_nonlazy, AlgebraElement ) :
    r"""
    An element of an algebra of equivariant monoid power series.
    """

    def __init__(self, parent, coefficients, precision, symmetrise,
                 cleanup_coefficients ) :
        r"""
        INPUT:
            - ``parent``       -- A ring or module of equivariant monoid power series.
            - ``coefficients`` -- A dictionary with keys in the parent's monoid and values
                                  in the parent coefficient domain.
            - ``precision``    -- A filter for the parent's action.
            - ``symmetrise``   -- A boolean (default: ``False``); If ``True`` every enty in
                                  ``coefficients`` will contribute to its whole orbit.
            - ``cleanup_coefficients`` -- A boolean (default: ``False``); If ``True`` zero
                                          coefficients will be erased from the dictionary.

        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
            sage: h = EquivariantMonoidPowerSeries_algebraelement(emps, dict(), emps.action().zero_filter(), False, False)
        """
        AlgebraElement.__init__(self, parent)
        EquivariantMonoidPowerSeries_abstract_nonlazy.__init__(self, parent, coefficients, precision, symmetrise,
                                            cleanup_coefficients)
