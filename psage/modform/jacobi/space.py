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
 Implements a class for a (vector) space of Jacobi forms.

 AUTHORS::

  - Stephan Ehlen
  - Nils Skoruppa
  - Fredrik Strömberg

 EXAMPLES::
"""

from sage.modules.free_module import *
from sage.all import QQ, Parent, SageObject
from psage.modform.jacobi.jacobiform import *

class JacobiFormsSpace_class(Parent):

    #Element = JacobiForm_space_element
 
    def __init__(self, weight, lattice, character, ambient_module=None):
        self._k = weight
        self._L = lattice
        self._h = character
        self._ambient_module = ambient_module

    def ambient_space(self):
        return self

    def ambient_module(self):
        r"""
          Returns the ambient module of self.
        """
        if self._ambient_module == None:
            self._ambient_module = JacobiFormsModule(self._L, self._h)
        return self._ambient_module

    def weight(self):
        return self._k

    def lattice(self):
        return self._L

    def index(self):
        return self._L

    def character(self):
        return self._h

    def dimension(self):
        raise NotImplementedError("Currently not implemented.")

    def dimension_cusp_forms(self):
        raise NotImplementedError("Currently not implemented.")

    def dimension_eisenstein_forms(self):
        raise NotImplementedError("Currently not implemented.")

    def vector_valued_forms(self):
        raise NotImplementedError("Currently not implemented.")

    def __repr__(self):
        return "Space of Jacobi forms of weight {0}, index {1}, character epsilon^{2} of dimension {3}"\
               .format(self._k, self._L, self._h, self.dimension())

class JacobiFormsSubspace_class(JacobiFormsSpace_class):

    def __init__(self,ambient_space):
        self._ambient_space = ambient_space

    def ambient_space(self):
        return self._ambient_space

class JacobiCuspFormsSpace_class(JacobiFormsSubspace_class):

    def __init__(self, weight, lattice, character, ambient_space=None):
        super(JacobiCuspFormsSpace_class, self).__init__(ambient_space)
        super(JacobiCuspFormsSpace_class, self).__init__(weight, lattice, character)

    def dimension(self):
        return self.dimension_cusp_forms()
    
    def dimension_eisenstein_forms(self):
        raise ValueError("This space is cuspidal.")

    def __repr__(self):
        return "Space of Jacobi cusp forms of weight {0}, index {1}, character epsilon^{2} of dimension {3}"\
               .format(self._k, self._L, self._h, self.dimension())


class JacobiEisensteinFormsSpace_class(JacobiFormsSpace_class):

    def __init__(self, weight, lattice, character, ambient_space=None):
        super(JacobiEisensteinFormsSpace_class, self).__init__(ambient_space)
        super(JacobiEisensteinFormsSpace_class, self).__init__(weight, lattice, character)

    def dimension(self):
        return self.dimension_eisenstein_forms()

    def dimension_cusp_forms(self):
        raise ValueError("This space is an Eisenstein space.")

    def __repr__(self):
        return "Space of Jacobi Eisenstein forms of weight {0}, index {1}, character epsilon^{2} of dimension {3}"\
               .format(self._k, self._L, self._h, self.dimension())

def JacobiForms(weight, lattice, character):
    return JacobiFormsSpace_class(weight, lattice, character)

def JacobiEisensteinForms(weight, lattice, character):
    return JacobiEisensteinFormsSpace_class(weight, lattice, character)

def JacobiCuspForms(weight, lattice, character):
    return JacobiCuspFormsSpace_class(weight, lattice, character)
