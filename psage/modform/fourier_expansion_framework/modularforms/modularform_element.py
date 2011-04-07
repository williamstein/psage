from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ambient import GradedExpansionAmbient_abstract
r"""
A orthogonal modular form, namely a graded expansion providing additional features.

AUTHOR :
    -- Martin Raum (2009 - 07 - 30) Initial version
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

from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_element import GradedExpansion_class, \
                            GradedExpansionVector_class, GradedExpansion_abstract

class ModularForm_abstract (object) :
    """
    NOTE:
        We assume that the deriving classes also derive (indirectly)
        from GradedExpansion_abstract.
    """
    
    def is_cusp_form(self) :
        """
        Whether ``self`` is a cusp form or not.
        
        OUTPUT:
            A boolean.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: ma = ModularFormsAmbient( QQ, ModularFormTestType_scalar(), NNFilter(5) )
            sage: ma.0.is_cusp_form()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def is_eisenstein_series(self) :
        """
        Whether ``self`` is an Eisenstein series or not.
        
        OUTPUT:
            A boolean.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: ma = ModularFormsAmbient( QQ, ModularFormTestType_scalar(), NNFilter(5) )
            sage: ma.0.is_eisenstein_series()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def _lmul_(self, c) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: ma = ModularFormsAmbient( QQ, ModularFormTestType_scalar(), NNFilter(5) )
            sage: mavv = ModularFormsAmbient( QQ, ModularFormTestType_vectorvalued(), NNFilter(5) )
            sage: (mavv.0 * ma.0).polynomial()
            g1*v1
            sage: (ma.0 * 5).polynomial()
            5*g1
        """
        ## For vector valued forms we use an extended base ring, which
        ## needs conversion before we can multiply
        if self.parent().type().is_vector_valued() :
            return self.parent()._element_class( self.parent(),
               c.polynomial().subs(c.parent().type()._hom_to_vector_valued(self.parent().relations().base_ring())) \
             * self.polynomial() )
        else :
            return super(ModularForm_abstract, self)._lmul_(c)

    def _rmul_(self, c) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: ma = ModularFormsAmbient( QQ, ModularFormTestType_scalar(), NNFilter(5) )
            sage: mavv = ModularFormsAmbient( QQ, ModularFormTestType_vectorvalued(), NNFilter(5) )
            sage: (ma.0 * mavv.0).polynomial()
            g1*v1
            sage: (5 * ma.0).polynomial()
            5*g1
        """
        ## For vector valued forms we use an extended base ring, which
        ## needs conversion before we can multiply
        if self.parent().type().is_vector_valued() :
            return self.parent()._element_class( self.parent(),
               c.polynomial().subs(c.parent().type()._hom_to_vector_valued(self.parent().relations().base_ring())) \
             * self.polynomial() )
        else :
            return super(ModularForm_abstract, self)._rmul_(c)

class ModularForm_generic ( ModularForm_abstract, GradedExpansion_class ) :
    weight = GradedExpansion_class.grading_index
    

class ModularFormVector_generic ( ModularForm_abstract, GradedExpansionVector_class ) :
    weight = GradedExpansionVector_class.grading_index
