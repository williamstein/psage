r"""
An equivariant monoid power series with lazy evaluation.

AUTHOR :
    -- Martin Raum (2009 - 08 - 05) Initial version
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

from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ambient import EquivariantMonoidPowerSeriesAmbient_abstract
from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import EquivariantMonoidPowerSeries_abstract
from sage.structure.element import AlgebraElement,ModuleElement
from sage.misc.misc import union
from sage.modules.module import Module 
from sage.rings.ring import Ring

#===============================================================================
# EquivariantMonoidPowerSeries_lazy
#===============================================================================

def EquivariantMonoidPowerSeries_lazy(parent, precision, coefficient_function, components = None, bounding_precision = None) :
    r"""
    Construct a equivariant monoid power series, which calculates its coefficients
    on demand.
    
    INPUT:
        - ``parent``       -- An instance of :class:~`fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ambient.MonoidPowerSeriesAmbient_abstract`.
        - ``precision``    -- A filter for the parent's action.
        - ``coefficient_function`` -- A function returning for each pair of characters and
                                      monoid element a Fourier coefficients.
        - ``components``   -- ``None`` or a list of characters (default: ``None``). A list of components
                              that do not vanish. If ``None`` no component is assumed to be zero.
        - ``bounding_precision``   -- ``None`` or a filter for the parent's action. If not ``None``
                                      coefficients will be assumed to vanish outside the filter.
    
    OUTPUT:
        An instance of :class:~`EquivariantMonoidPowerSeries_abstract_lazy`.
    
    EXAMPLES::
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_lazyelement import *
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
        sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
        sage: h = EquivariantMonoidPowerSeries_lazy(emps, emps.action().filter(3), lambda (ch, k) : 1)
        sage: m = FreeModule(QQ, 3)
        sage: emps = EquivariantMonoidPowerSeriesModule( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", m) )
        sage: h = EquivariantMonoidPowerSeries_lazy(emps, emps.action().filter(3), lambda (ch, k) : m([1,2,3]))
        sage: h = EquivariantMonoidPowerSeries_lazy(emps, emps.action().filter_all(), lambda (ch, k) : m([1,2,3]), bounding_precision = emps.action().filter(2))
        sage: h = EquivariantMonoidPowerSeries_lazy(emps, emps.action().filter_all(), lambda (ch, k) : m([1,2,3]))
        Traceback (most recent call last):
        ...
        ValueError: Lazy equivariant monoid power series cannot have infinite precision.
    """
    if (bounding_precision is None or bounding_precision.is_infinite()) and \
       precision.is_infinite() and not (isinstance(components, list) and len(components) == 0) :
        raise ValueError( "Lazy equivariant monoid power series cannot have infinite precision." )
    if isinstance(parent, Module) :
        return EquivariantMonoidPowerSeries_moduleelement_lazy(parent, precision, coefficient_function, components, bounding_precision)
    if isinstance(parent, Ring) :
        return EquivariantMonoidPowerSeries_algebraelement_lazy(parent, precision, coefficient_function, components, bounding_precision)

#===============================================================================
# EquivariantMonoidPowerSeries_abstract_lazy
#===============================================================================

class EquivariantMonoidPowerSeries_abstract_lazy (EquivariantMonoidPowerSeries_abstract) :
    r"""
    This class implements an equivariant monoid power series which calculates the
    coefficients on demand.
    """

    def __init__(self, parent, precision, coefficient_function, components, bounding_precision = None ) :
        r"""
        INPUT:
            - ``parent``               -- An instance of :class:~`fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ambient.MonoidPowerSeriesAmbient_abstract`.
            - ``precision``            -- A filter for the parent's action.
            - ``coefficient_function`` -- A function returning for each pair of characters and
                                          monoid element a Fourier coefficients.
            - ``components``           -- ``None`` or a list of characters. A list of components that do not
                                          vanish. If ``None`` no component is assumed to be zero.
            - ``bounding_precision``   -- ``None`` or a filter for the parent's action. If not ``None``
                                          coefficients will be assumed to vanish outside the filter.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_lazyelement import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
            sage: h = EquivariantMonoidPowerSeries_abstract_lazy(emps, emps.action().filter(3), lambda (ch, k) : 1, None)
            sage: h = EquivariantMonoidPowerSeries_abstract_lazy(emps, emps.action().filter(3), lambda (ch, k) : 1, [])
            sage: h = EquivariantMonoidPowerSeries_abstract_lazy(emps, emps.action().filter_all(), lambda (ch, k) : 1, [], emps.action().filter(4))
        """
        EquivariantMonoidPowerSeries_abstract.__init__(self, parent, precision)

        self.__bounding_precision = bounding_precision
        self.__coefficient_function = coefficient_function
        self.__coefficients = dict()
        self.__coefficients_complete = False
        self.__components = components
        
    def non_zero_components(self) :
        r"""
        Return all those characters which are not guaranteed to have only
        vanishing coefficients associated with.

        OUTPUT:
            A list of elements of the character monoid.
        
        NOTE:
            This will only return the list passed during the construction
            of self or a list of all characters.
             
        EXAMPLES::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_lazyelement import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
            sage: e = EquivariantMonoidPowerSeries(emps, {emps.characters().one_element() : { 1 : 1}}, emps.action().filter_all())
            sage: EquivariantMonoidPowerSeries_LazyMultiplication(e, e).non_zero_components()
            [1]
            sage: EquivariantMonoidPowerSeries_abstract_lazy(emps, emps.action().filter(3), lambda (ch, k) : 1, []).non_zero_components()
            []
        """
        if self.__components is None :
            self.__components = [c for c in self.parent().characters()] 

        return self.__components
    
    def _bounding_precision(self) :
        r"""
        If a filter for the vanishing of coefficients is given return this. Otherwise,
        return the precision, which is then guaranteed to be finite.

        OUTPUT:
            A filter for the parent's action.
            
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_lazyelement import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
            sage: e = EquivariantMonoidPowerSeries(emps, {}, emps.action().filter_all())
            sage: EquivariantMonoidPowerSeries_LazyMultiplication(e, e)._bounding_precision()
            Filtered NN with action up to 0
            sage: EquivariantMonoidPowerSeries_abstract_lazy(emps, emps.action().filter(3), lambda (ch, k) : 1, [emps.characters().one_element()])._bounding_precision()
            Filtered NN with action up to 3
            sage: h = EquivariantMonoidPowerSeries_abstract_lazy(emps, emps.action().filter_all(), lambda (ch, k) : 1, [])._bounding_precision() ## This would call the parent
            Traceback (most recent call last):
            ...
            AttributeError: EquivariantMonoidPowerSeries_abstract_lazy instance has no attribute 'parent'
            sage: h = EquivariantMonoidPowerSeries_abstract_lazy(emps, emps.action().filter_all(), lambda (ch, k) : 1, [], emps.action().filter(4))._bounding_precision() ## This would call the parent
            Traceback (most recent call last):
            ...
            AttributeError: EquivariantMonoidPowerSeries_abstract_lazy instance has no attribute 'parent'
            sage: EquivariantMonoidPowerSeries_abstract_lazy(emps, emps.action().filter_all(), lambda (ch, k) : 1, [emps.characters().one_element()], emps.action().filter(4))._bounding_precision()
            Filtered NN with action up to 4
            sage: EquivariantMonoidPowerSeries_abstract_lazy(emps, emps.action().filter_all(), lambda (ch, k) : 1, [emps.characters().one_element()])._bounding_precision()
            Traceback (most recent call last):
            ...
            ValueError: No bounding precision for ...
        """
        if len(self.non_zero_components()) == 0 :
            return self.parent().action().zero_filter() 
        elif not self.__bounding_precision is None :
            return min(self.__bounding_precision, self.precision())
        elif self.precision().is_infinite() :
            raise ValueError( "No bounding precision for %s." % (self,) )
        
        return self.precision()
    
    def coefficients(self, force_characters = False) :
        r"""
        Evaluate all coefficients within the precision bounds and return a
        dictionary which saves all coefficients of this element.
    
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
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_lazyelement import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
            sage: e = EquivariantMonoidPowerSeries_abstract_lazy(emps, emps.action().filter(3), lambda (ch, k) : 1, [emps.characters().one_element()])
            sage: e.coefficients()
            {0: 1, 1: 1, 2: 1}
            sage: e.coefficients(True)
            {1: {0: 1, 1: 1, 2: 1}}
        """
        if not self.__coefficients_complete :
            self.__coefficients_complete = True
            
            coefficient_function = self.__coefficient_function
            
            for ch in self.non_zero_components() :
                if not ch in self.__coefficients :
                    coeffs = dict()
                    self.__coefficients[ch] = coeffs
                    for k in self._bounding_precision() :
                        coeffs[k] = coefficient_function((ch, k))
                else :
                    coeffs = self.__coefficients[ch]
                    for k in self._bounding_precision() :
                        if not k in coeffs :
                            coeffs[k] = coefficient_function((ch, k))

        if len(self.__coefficients) == 0 and not force_characters :
            return dict()
        elif len(self.__coefficients) == 1 and not force_characters :
            return self.__coefficients.values()[0]
        else :
            return self.__coefficients

    def _truncate_in_place(self, precision) :
        r"""
        Truncate ``self`` modifying also the coefficient cache.
        
        INPUT:
            - ``precision``    -- A filter for the parent's action or a an object that can be converted
                               to a filter.

        OUTPUT:
            ``None``
        
        EXAMPLE:
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_lazyelement import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
            sage: h = EquivariantMonoidPowerSeries_abstract_lazy(emps, emps.action().filter(3), lambda (ch, k) : 1, None)
            sage: e = EquivariantMonoidPowerSeries(emps, {emps.characters().one_element() : { 1 : 1}}, emps.action().filter_all())
            sage: e = EquivariantMonoidPowerSeries_LazyMultiplication(e, e)
            sage: e._truncate_in_place(emps.monoid().filter(2))
            sage: e.coefficients()
            {0: 0, 1: 0}
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

    def __getitem__(self, k) :
        r"""
        Evaluate and return the `k`-th coefficient if it below the series' precision. If no character is contained
        in the key ``self`` must have only one nonvanishing component.
        
        INPUT:
            - `k` -- A pair of an element of the parent's character monoid and
                     and element of the parent's monoid or an element of the parent's monoid.
        
        OUTPUT:
            An element of the parent's coefficient domain.

        EXAMPLES::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_lazyelement import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
            sage: e = EquivariantMonoidPowerSeries(emps, {emps.characters().one_element() : { 1 : 1, 2 : 4}}, emps.action().filter(3))
            sage: e = EquivariantMonoidPowerSeries_LazyMultiplication(e, e)
            sage: e[(emps.characters().one_element(),2)]
            1
            sage: e[1], e[2]
            (0, 1)
        """
        try :
            (ch, s) = k
        except TypeError :
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
            
            s = k
            
        if not s in self.precision() :
            raise ValueError, "%s out of bound" % s

        try :
            return self.__coefficients[ch][s]
        except KeyError :
            (rs, g) = self.parent()._reduction_function()(s)
            try :
                return self.parent()._character_eval_function()(g, ch) \
                 * self.parent()._apply_function()(g, self.__coefficients[ch][rs])
            except KeyError :
                e = self.__coefficient_function((ch, rs))

                try :
                    self.__coefficients[ch][rs] = e
                except KeyError :
                    self.__coefficients[ch] = dict()
                    self.__coefficients[ch][rs] = e
                
                return e

#===============================================================================
# EquivariantMonoidPowerSeries_modulelement_lazy
#===============================================================================

class EquivariantMonoidPowerSeries_moduleelement_lazy (EquivariantMonoidPowerSeries_abstract_lazy, ModuleElement) :
    r"""
    A lazy element of a module of equivariant monoid power series.
    """

    def __init__(self, parent, precision, coefficient_function, components, bounding_precision ) :
        r"""
        INPUT:
            - ``parent``       -- An instance of :class:~`fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ambient.MonoidPowerSeriesAmbient_abstract`.
            - ``precision``    -- A filter for the parent's action.
            - ``coefficient_function`` -- A function returning for each pair of characters and
                                          monoid element a Fourier coefficients.
            - ``components``   -- ``None`` or a list of characters. A list of components that do not
                                  vanish. If ``None`` no component is assumed to be zero.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_lazyelement import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: m = FreeModule(QQ, 3)
            sage: emps = EquivariantMonoidPowerSeriesModule( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", m) )
            sage: h = EquivariantMonoidPowerSeries_moduleelement_lazy(emps, emps.action().filter(3), lambda (ch, k) : m([1,2,3]), None, None)
        """
        ModuleElement.__init__(self, parent)
        EquivariantMonoidPowerSeries_abstract_lazy.__init__(self, parent, precision,
                 coefficient_function, components, bounding_precision)

#===============================================================================
# EquivariantMonoidPowerSeries_algebraelement_lazy
#===============================================================================

class EquivariantMonoidPowerSeries_algebraelement_lazy (EquivariantMonoidPowerSeries_abstract_lazy, AlgebraElement) :
    r"""
    A lazy element of an algebra of equivariant monoid power series.
    """

    def __init__(self, parent, precision, coefficient_function, components, bounding_precision ) :
        r"""
        INPUT:
            - ``parent``       -- An instance of :class:~`fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ambient.MonoidPowerSeriesAmbient_abstract`.
            - ``precision``    -- A filter for the parent's action.
            - ``coefficient_function`` -- A function returning for each pair of characters and
                                          monoid element a Fourier coefficients.
            - ``components``   -- ``None`` or a list of characters. A list of components that do not
                                  vanish. If ``None`` no component is assumed to be zero.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_lazyelement import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
            sage: h = EquivariantMonoidPowerSeries_algebraelement_lazy(emps, emps.action().filter(3), lambda (ch, k) : 1, None, None)
        """
        AlgebraElement.__init__(self, parent)
        EquivariantMonoidPowerSeries_abstract_lazy.__init__(self, parent, precision,
                 coefficient_function, components, bounding_precision)

#===============================================================================
# EquivariantMonoidPowerseries_MultiplicationDelayedFactory
#===============================================================================

class EquivariantMonoidPowerseries_MultiplicationDelayedFactory :
    r"""
    A helper class for lazy multplication of equivariant monoid power series.
    """
    
    def __init__( self, left, right ) :
        r"""
        INPUT:
            - ``left``  -- An instance of :class:~`fourier_expansion_framework.monoidpowerseries.EquivariantMonoidPowerSeries_abstract`.
            - ``right`` -- An instance of :class:~`fourier_expansion_framework.monoidpowerseries.EquivariantMonoidPowerSeries_abstract`.

        NOTE:
            ``left`` and ``right`` must have the same parents.

        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_lazyelement import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
            sage: e = EquivariantMonoidPowerSeries(emps, {emps.characters().one_element() : { 1 : 1}}, emps.action().filter_all())
            sage: EquivariantMonoidPowerSeries_LazyMultiplication(e, e).coefficients() # indirect doctest
            {0: 0, 1: 0, 2: 1}
        """
        assert left.parent() == right.parent()
        
        self.__left = left
        self.__right = right
        
        self.__mul_fc = left.parent()._multiply_function()
        self.__coefficient_ring = left.parent().coefficient_domain()
                
        self.__left_coefficients = None
        self.__right_coefficients = None

    def getcoeff(self, (ch, k)) :
        r"""
        Return the `k`-th coefficient of the component ``ch`` of the product.
        
        INPUT:
            - ``(ch, k)`` -- A pair of character and monoid element.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_lazyelement import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
            sage: e = EquivariantMonoidPowerSeries(emps, {emps.characters().one_element() : { 1 : 1}}, emps.action().filter_all())
            sage: EquivariantMonoidPowerSeries_LazyMultiplication(e, e).coefficients() # indirect doctest
            {0: 0, 1: 0, 2: 1}
        """        
        res = self.__coefficient_ring(0)
        
        if self.__left_coefficients is None :
            self.__left_coefficients = self.__left.coefficients(True)
        if self.__right_coefficients is None :
            self.__right_coefficients = self.__right.coefficients(True)
        
        for c1 in self.__left_coefficients :
            lcoeffs = self.__left_coefficients[c1]
            if len(lcoeffs) == 0 : continue
            
            for c2 in self.__right_coefficients :
                rcoeffs = self.__right_coefficients[c2]
                if len(rcoeffs) == 0 : continue
                
                if c1 * c2 != ch :
                    continue
                
                res += self.__mul_fc( k, lcoeffs, rcoeffs, c1, c2, self.__coefficient_ring(0) )

        return res

#===============================================================================
# EquivariantMonoidPowerSeries_LazyMultiplication
#===============================================================================

def EquivariantMonoidPowerSeries_LazyMultiplication(left, right) :
    r"""
    The product of two equivariant monoid power series which is only evaluated
    if a coefficient is demanded.
    
    INPUT:
        - ``left``  -- An instance of :class:~`fourier_expansion_framework.monoidpowerseries.EquivariantMonoidPowerSeries_abstract`.
        - ``right`` -- An instance of :class:~`fourier_expansion_framework.monoidpowerseries.EquivariantMonoidPowerSeries_abstract`.
    
    EXAMPLES:
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_lazyelement import *
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import *
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import *
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
        sage: emps = EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ) )
        sage: e = EquivariantMonoidPowerSeries(emps, {emps.characters().one_element() : { 1 : 1}}, emps.action().filter_all())
        sage: EquivariantMonoidPowerSeries_LazyMultiplication(e, e).coefficients() # indirect doctest
        {0: 0, 1: 0, 2: 1}
    """
    # TODO: Insert coercing
    if not isinstance(left.parent(), EquivariantMonoidPowerSeriesAmbient_abstract) \
       or not isinstance(right.parent(), EquivariantMonoidPowerSeriesAmbient_abstract) :
        raise TypeError, "both factors must be power series"
       
    if left.parent() != right.parent() :
        if left.parent() is right.parent().base_ring() :
            parent = right.parent()
        elif left.parent().base_ring() is right.parent() :
            parent = left.parent()
        
        raise ValueError, "incorrect parents of the factors"
    else :
        parent = left.parent()
    
    coefficients_factory = \
     EquivariantMonoidPowerseries_MultiplicationDelayedFactory( left, right )

    precision = min(left.precision(), right.precision())
    bounding_precision = None
    if precision.is_infinite() :
        left_coefficients = left.coefficients(True)
        right_coefficients = right.coefficients(True)
        
        left_keys  = reduce(union, (set(c) for c in left_coefficients.itervalues()), set())
        right_keys = reduce(union, (set(c) for c in right_coefficients.itervalues()), set())
        
        bounding_precision = left.parent().action(). \
                        minimal_composition_filter(left_keys, right_keys)
    
    return EquivariantMonoidPowerSeries_lazy( parent,
                                              precision,
                                              coefficients_factory.getcoeff,
                                              bounding_precision = bounding_precision )
