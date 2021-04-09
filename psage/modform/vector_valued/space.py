# -*- coding: utf-8 -*-
#*****************************************************************************
#  Copyright (C) 2013 
#
#  Distributed under the terms of the GNU General Public License (GPLv2)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
r"""
 Implements a class for a (vector) space of vector valued modular forms.

 AUTHORS::

  - Stephan Ehlen (Initial version)

 EXAMPLES::
"""

from sage.all import SageObject
from sage.modules.free_module import *
from sage.modular.arithgroup.congroup_sl2z import *

class VectorValuedModularFormsSpace_generic_class(FreeModule_generic):
 
    def __init__(self, weight, module, ambient_module=None):
        self._k = weight
        self._M = module
        self._ambient_module = ambient_module
        rank = self.__calculate_rank()
        super(JacobiFormsModule_class,self).__init__(ModularForms(SL2Z), rank, rank)

    def ambient_space(self):
        return self

    def ambient_module(self):
        r"""
          Returns the ambient module of self.
        """
        return NotImplementedError("This method is currently not implemented. It should be overridden by the specific subclasses.")

    def weight(self):
        return self._k

    def module(self):
        return self._M

    def dimension(self):
        return NotImplementedError("Currently not implemented.")

    def dimension_cusp_forms(self):
        return NotImplementedError("Currently not implemented.")

    def dimension_eisenstein_forms(self):
        return NotImplementedError("Currently not implemented.")

    def __repr__(self):
        return "Space of vector valud modular forms of weight {0} and representation given by the module {1}, of dimension {2}"\
               .format(self._k, self._M, self.dimension())

class VectorValuedModularFormsSubspace_generic_class(VectorValuedModularFormsSpace_generic_class):

    def __init__(ambient_space):
        self._ambient_space = ambient_space

    def ambient_space(self):
        return self._ambient_space

class VectorValuedCuspFormsSpace_generic_class(VectorValuedModularFormsSubspace_generic_class):

    def __init__(self, weight, lattice, character, ambient_space=None):
        super(VectorValuedCuspFormsSpace_generic_class, self).__init__(ambient_space)
        super(VectorValuedCuspFormsSpace_generic_class, self).__init__(weight, module)

    def dimension(self):
        return self.dimension_cusp_forms()

    def __repr__(self):
        return "Vector valued cusp forms of weight {0} and representation given by the module {1}, of dimension {2}"\
               .format(self._k, self._M, self.dimension())


class VectorValuedEisensteinFormsSpace_generic_class(VectorValuedModularFormsSubspace_generic_class):

    def __init__(self, weight, module, ambient_space=None):
        super(VectorValuedEisensteinFormsSpace_generic_class, self).__init__(ambient_space)
        super(VectorValuedEisensteinFormsSpace_generic_class, self).__init__(weight, module)

    def dimension(self):
        return self.dimension_eisenstein_forms()

    def __repr__(self):
        return "Vector valued Eisenstein forms of weight {0} and representation given by the module {1}, of dimension {2}"\
               .format(self._k, self._M, self.dimension())

def VectorValuedModularForms(weight, module):
    return VectorValuedModularFormsSpace_class(weight, module)

def VectorValuedEisensteinForms(weight, lattice, character):
    return VectorValuedEisensteinFormsSpace_class(weight, module)

def VectorValuedCuspForms(weight, module):
    return VectorValuedCuspFormsSpace_class(weight, module)
