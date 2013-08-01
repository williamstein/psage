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
 Implements a class for a module of Jacobi forms (direct sum of the spaces of all weights).

 AUTHORS::

  - Stephan Ehlen
  - Nils Skoruppa
  - Fredrik Strömberg

 EXAMPLES::
"""

# python imports

from sage.modular.arithgroup.congroup_sl2z import *
from sage.all import Integer, ModularFormsRing, CommutativeRing, Infinity, SageObject, Parent
from psage.modform.jacobi.space import *
from psage.modform.jacobi.jacobiform import *

class ModularFormsForSL2Z(ModularFormsRing, CommutativeRing):

    def __init__(self):
        super(ModularFormsForSL2Z, self).__init__(SL2Z)

class JacobiFormsModule_class(Module):

    Element = JacobiForm_module_element
 
    def __init__(self, lattice, character):
        self._L = lattice
        self._h = character
        self._rank = self.__calculate_rank()
        super(JacobiFormsModule_class,self).__init__(ModularFormsForSL2Z())

    def lattice(self):
        return self._L

    def index(self):
        return self._L

    def character(self):
        return self._h

    def graded_component(self, k):
        return JacobiFormsSpace_class(k, self._L, self._h, self)

    def subspace (self, k):
        return self.graded_component(k)

    def generators(self):
        raise NotImplementedError("Currently not implemented")

    def rank(self):
        return self._rank

    def __calculate_rank(self):
        # this will work once the lattice class is used and Nils's
        # implementation of the module class
        if self._L != None:
            L = self._L
            o_inv = L.o_invariant()
            V2 = L.bullet_vectors_of_order_2()
            return (Integer(L.det()) + Integer(len(V2)*(-1)**(self.__par + o_inv)).divide_knowing_divisible_by(2))
        else:
        # stupid... just for testing purposes
            return Integer(1)

    def dimension(self, k):
        return self.graded_component(k).dimension()

    def zero(self):
        return JacobiForm_module_element_zero(self)

    def __repr__(self):
        return "Module of Jacobi forms of index {0}, character epsilon^{1}.".format(self._L, self._h)


