r"""
A functor creating rings of orthogonal modular forms.

AUTHOR:
    -- Martin Raum (2009 - 07 - 30) Initial version
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

from sage.categories.rings import Rings
from sage.categories.pushout import ConstructionFunctor

class ModularFormsFunctor (ConstructionFunctor):
    
    rank = 10
    
    def __init__(self, type, precision):
        """
        A functor constructing a ring or module of modular forms.
        
        INPUT:
            - ``type``      -- A type of modular forms.
            - ``precision`` -- A precision.
        
        NOTE:
            This does not respect keyword and has to be extended as soon as subclasses of
            ModularFormsAmbient_abstract demand for it.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_functor import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: F = ModularFormsFunctor( ModularFormTestType_scalar(), NNFilter(5) )
        """
        self.__type = type
        self.__precision = precision
        
        ConstructionFunctor.__init__(self, Rings(), Rings())

    def __call__(self, A):
        """
        INPUT:
            - `A` -- A ring.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_functor import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: F = ModularFormsFunctor( ModularFormTestType_scalar(), NNFilter(5) )
            sage: F(QQ)
            Graded expansion ring with generators g1, g2, g3, g4, g5
        """
        from .modularform_ambient import ModularFormsAmbient

        return ModularFormsAmbient(A, self.__type, self.__precision)
        
    def merge(self, other):
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_functor import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: F = ModularFormsFunctor( ModularFormTestType_scalar(), NNFilter(5) )
            sage: G = ModularFormsFunctor( ModularFormTestType_vectorvalued(), NNFilter(5) )
            sage: F.merge(F) is F
            True
            sage: F.merge(G) is None
            True
            sage: G.merge(F) is G 
            True
        """
        if type(other) != type(self):
            return None
        
        if self.__type == other.__type and \
           self.__precision == other.__precision:
            return self
        else:
            try:
                if other.__type.vector_valued() == self.__type and \
                   self.__precision == other.__precision:
                    return self
            except (AttributeError,NotImplementedError):
                return None 
    
        return None
