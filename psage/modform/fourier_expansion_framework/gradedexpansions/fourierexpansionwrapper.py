r"""
Interface for elements which provide a Fourier expansion.

AUTHOR :
    -- Martin Raum (2009 - 08 - 03) Initial version
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

class FourierExpansionWrapper :
    r"""
    Abstract class for elements, which do not represent Fourier
    expansions on their own, but encapsulate one.
    
    SEE:
        :class:~`fourier_expansion_framework.gradedexpansions.GradedExpansion_class`.
    """
    
    def fourier_expansion(self, cache = True) :
        r"""
        The Fourier expansion which is associated with ``self``.
        
        INPUT:
            - ``cache`` -- A boolean (default: ``True``); If ``True`` the return value
                           will be cached.
        
        OUTPUT:
            A (equivariant) monoid power series.
        
        NOTE:
            The parent has to implement the function ``_fourier_expansion_of_element``.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: em = ExpansionModule(Sequence([emps.one_element().truncate(3), emps.one_element(), emps.one_element()]))
            sage: h = em([1,2,-3])
            sage: fe = h.fourier_expansion()
        """
        try :
            return self.__fourier_expansion
        except AttributeError :
            if cache :
                self.__fourier_expansion = \
                  self.parent()._fourier_expansion_of_element(self)
                return self.__fourier_expansion
            else :
                return self.parent()._fourier_expansion_of_element(self)
            
    def _set_fourier_expansion(self, expansion ) :
        r"""
        Set the cache for the Fourier expansion to an given monoid power series.
        
        INPUT:
            - ``expansion`` -- A (equivariant) monoid power series.
        
        OUTPUT:
            ``None``
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: em = ExpansionModule(Sequence([emps.one_element().truncate(3), emps.one_element(), emps.one_element()]))
            sage: h = em([1,2,5])
            sage: fe = h.fourier_expansion()
            sage: h = em([1,2,0])
            sage: h._set_fourier_expansion(fe)
            sage: h.fourier_expansion() == fe
            True
        """
        self.__fourier_expansion = expansion
