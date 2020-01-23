"""
Lazy monoid power series which wrap graded expansions.

AUTHOR :
    -- Martin Raum (2010 - 05 - 23) Initial version
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

from builtins import object
from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_lazyelement import EquivariantMonoidPowerSeries_lazy

def LazyFourierExpansionEvaluation(parent, element, precision) :
    """
    Create an lazy equivaraint monoid power series which evaluates a
    graded expansion.
    
    INPUT:
        - ``parent``    -- A ring or module of equivariant monoid power series. 
        - ``element``   -- A graded expansion element.
        - ``precision`` -- A filter for the parent's action.
    
    OUTPUT:
        An instance of :class:~`fourier_expansion_framework.monoidpowerseries.monoidpowerseries_lazyelement.EquivariantMonoidPowerSeries_abstract_lazy`.
    
    TESTS::
        sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
        sage: from psage.modform.fourier_expansion_framework.gradedexpansions.expansion_lazy_evaluation import LazyFourierExpansionEvaluation
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
        sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
        sage: em = ExpansionModule(Sequence([emps.one_element(), emps.one_element(), emps.one_element()]))
        sage: h = LazyFourierExpansionEvaluation(emps, em([1,2,3]), emps.action().filter(3))
        sage: h.coefficients()
        {0: 6, 1: 0, 2: 0}
    """
    if precision.is_infinite() :
        raise ValueError( "Lazy evaluation of infinite expansions is not possible")
    
    delayed_coeffs = DelayedEvaluation_fourier_expansion(parent, element)
    
    return EquivariantMonoidPowerSeries_lazy(parent, precision, delayed_coeffs.getcoeff)


class DelayedEvaluation_fourier_expansion(object) :
    """
    Helper class which evaluates the a graded expansion on and returns its
    coefficients.
    """
    
    def __init__(self, parent, element) :
        """
        INPUT:
            - ``parent``    -- A ring or module of equivariant monoid power series. 
            - ``element``   -- A graded expansion element.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.expansion_lazy_evaluation import DelayedEvaluation_fourier_expansion
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: em = ExpansionModule(Sequence([emps.one_element(), emps.one_element(), emps.one_element()]))
            sage: de = DelayedEvaluation_fourier_expansion(emps, em([1,2,5]))
        """
        self.__parent = parent
        self.__element = element
        
    def getcoeff(self, key) :
        """
        Return a coefficient of the Fourier expansion of a graded expansion.
        
        INPUT:
            - ``key`` -- A pair ``(ch, k)`` of a character and a monoid element.
        
        OUTPUT:
            An element of the parents coefficient domain.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.expansion_lazy_evaluation import DelayedEvaluation_fourier_expansion
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: em = ExpansionModule(Sequence([emps.one_element(), emps.one_element(), emps.one_element()]))
            sage: de = DelayedEvaluation_fourier_expansion(emps, em([1,2,-2]))
            sage: (de.getcoeff(0), de.getcoeff(1))
            (1, 0)
        """
        try :
            return self.__fourier_expansion[key]
        except AttributeError :
            self.__fourier_expansion = self.__element.fourier_expansion()
            if self.__fourier_expansion.parent().coefficient_domain() != self.__parent.coefficient_domain() :
                self.__fourier_expansion = self.__parent(self.__fourier_expansion)
            
            return self.__fourier_expansion[key]
