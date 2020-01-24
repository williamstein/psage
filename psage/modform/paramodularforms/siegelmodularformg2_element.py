r"""
Elements of rings of Siegel modular forms of degree 2.

AUTHOR:

- Martin Raum (2009 - 08 - ??) Initial version.
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

from builtins import range
from builtins import object
from psage.modform.fourier_expansion_framework.modularforms.modularform_element import ModularForm_generic, \
                                             ModularFormVector_generic

class SiegelModularFormG2_generic(object) :

    def is_cusp_form(self) :
        r"""
        Check wheater self is a cusp form or not.
        """
        
        if self.parent().type().group() != "Sp(2,ZZ)" :
            raise NotImplementedError
        
        try :
            weight = self.weight() 
        except ValueError :
            return all([ f.is_cusp_form()
                         for f in self.homogeneous_components().values() ])
        
        neccessary_precision = weight // 12 + 1 if weight % 12 != 2 \
                                                else weight // 12
        if self.parent().fourier_expansion_precision()._indefinite_content_bound() < neccessary_precision :
            raise ValueError("the parents precision doesn't suffice")
        
        evc = self.fourier_expansion()
        return all([evc[(0,0,l)] == 0 for l in range(neccessary_precision)])
        
class SiegelModularFormG2_classical ( SiegelModularFormG2_generic, ModularForm_generic ) :

    def is_maass_form(self) :    
        r"""
        Check whether self is in the Maass spezialschar.
        """
        hc = self.homogeneous_components()
        if len(hc) == 0 : return True
        if len(hc) > 1 : return False
        
        
        grading = list(hc.keys())[0]
        
        ss = self.parent().graded_submodule(grading)

        return ss(self) in ss.maass_space()

class SiegelModularFormG2_vectorvalued ( SiegelModularFormG2_generic,
                                         ModularFormVector_generic ) :
    pass
