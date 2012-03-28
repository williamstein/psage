r"""
Abstract subspaces of Siegel modular forms.

AUTHORS:

- Martin Raum (2009 - 08 - ??) Initial version
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

from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_submodule import GradedExpansionSubmoduleVector_generic
from psage.modform.fourier_expansion_framework.modularforms.modularform_submodule import ModularFormsSubmodule_singleweight_ambient_pid, \
                                               ModularFormsSubmodule_heckeinvariant_submodule
from sage.misc.cachefunc import cached_method

#===============================================================================
# DiscrimiantPrecisionSubmodule_abstract
#===============================================================================

class DiscrimiantPrecisionSubmodule_abstract :
    def _minimal_discriminant_precision(self, minimal_precision = 0, lazy_rank_check = False) :
        last_precision = None
        for p in xrange(minimal_precision, self.graded_ambient().precision().discriminant() + 1) :
            precision = self.graded_ambient().fourier_ring().filter(p)
            if not precision < last_precision : continue
             
            if self._check_precision(self, precision, lazy_rank_check) :
                return precision
        else :
            raise ValueError( "This basis is not determined completely by its Fourier expansions." )

class SiegelModularFormG2WeightSubmodule_class ( ModularFormsSubmodule_singleweight_ambient_pid,
                                                 DiscrimiantPrecisionSubmodule_abstract ) :
    @cached_method
    def maass_space(self) :
        return SiegelModularFormG2Submodule_maassspace(
            self, map(self, self.graded_ambient().type(). \
                            _maass_generators( self.weight(),
                                               self.graded_ambient().fourier_expansion_precision())) )    
    
class SiegelModularFormG2Submodule_maassspace (
        ModularFormsSubmodule_heckeinvariant_submodule, DiscrimiantPrecisionSubmodule_abstract ) :
 
    def __init__(self, ambient, basis, **kwds) :
        self.__provided_basis = basis
        ModularFormsSubmodule_heckeinvariant_submodule.__init__(self, ambient, basis, **kwds)
    
    def weight(self) :
        return self.ambient_module().weight()

    def _provided_basis(self) :
        return self.__provided_basis

class SiegelModularFormG2SubmoduleVector_generic ( GradedExpansionSubmoduleVector_generic ) :

    def is_cusp_form(self) :
        r"""
        Check wheater self is a cusp form or not.
        """
        if self.parent().ring().type().group() != "Sp(2,ZZ)" :
            raise NotImplementedError
        
        try :
            weight = self.parent().weight() 
        except AttributeError :
            raise TypeError, "Can check cusps only for spaces of homogeneous weight"
        
        neccessary_precision = weight // 12 + 1 if weight % 12 != 2 \
                                                else weight // 12
        if self.parent().ring().fourier_expansion_precision()._indefinite_content_bound() < neccessary_precision :
            raise ValueError, "the parents precision doesn't suffice"
        
        evc = self.fourier_expansion()
        return all([evc[(0,0,l)] == 0 for l in xrange(neccessary_precision)])
