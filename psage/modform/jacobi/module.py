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

from builtins import map
from sage.modular.arithgroup.congroup_sl2z import *
from sage.all import Integer, ModularFormsRing, CommutativeRing, Infinity, SageObject, Parent, IntegerRing, EisensteinForms
from psage.modform.jacobi.space import *
from psage.modform.jacobi.jacobiform import *
from sage.structure.unique_representation import *
from sage.structure.element_wrapper import *
from sage.categories.graded_modules_with_basis import *
from sage.categories.commutative_rings import *
from sage.modules.module import *
from sage.structure.element import *
from sage.modular.modform.element import ModularFormElement

class ModularFormsForSL2Z(UniqueRepresentation, Parent):

    def __init__(self):
        Parent.__init__(self, category=CommutativeRings(IntegerRing()))

    def product(self, x, y):
        x = x.value
        y = y.value
        return self(x*y)

    def an_element(self):
        E4 = self(EisensteinForms(1,4).basis()[0])
        return E4

    def __repr__(self):
        return "Ring of modular forms for the full modular group."

    class Element(ElementWrapper):
        wrapped_class = ModularFormElement
            
        
class JacobiFormsModule_generic_class(UniqueRepresentation, Module):
 
    def __init__(self, lattice, character):
        self._L = lattice
        self._h = character
        self._rank = self.__calculate_rank()
        Module.__init__(self, ModularFormsForSL2Z())
        self._base = ModularFormsForSL2Z()
        #super(JacobiFormsModule_class,self).__init__(ModularFormsForSL2Z())

    def an_element(self):
        A = ModularFormsForSL2Z()
        f = A.an_element()
        return self(4,definition=[(f,1)])
        

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

    def gens(self):
        return self.generators()

    def rank(self):
        return self._rank

    def __calculate_rank(self):
        # this will work once the lattice class is used and Nils's
        # implementation of the module class
        if self._L != None:
            L = self._L
            o_inv = L.o_invariant()
            V2 = L.bullet_vectors_of_order_2()
            return L.det() + Integer(len(V2)*(-1)**(self.__par + o_inv)).divide_knowing_divisible_by(2)
        else:
        # stupid... just for testing purposes
            return Integer(1)

    def dimension(self, k):
        return self.graded_component(k).dimension()

    def zero(self):
        return JacobiForm_module_element_zero(self)

    def __repr__(self):
        return "Module of Jacobi forms of index {0}, character epsilon^{1}.".format(self._L, self._h)

    class Element(JacobiForm_class, ModuleElement):

        def __init__(self, weight, parent, prec=10, definition=list()):
            if not isinstance(definition, list):
                raise ValueError("The argument `definition` has to be instance of type::`list`.")
            self._definition = definition
            ModuleElement.__init__(self, parent = parent)
            JacobiForm_class.__init__(self, weight, parent.lattice(), parent.character(), prec)

        def definition(self):
            return self._definition

        def __add__(self, other):
            f = other
            if not self._k == f.weight():
                raise ValueError("The weights have to agree.")
            return self.parent()(self._k, definition = self._definition + f.definition())

        def __mul__(self, f):
            kk = self._k + f.value.weight()
            newdef_map = lambda t: (f*t[0], t[1])
            newdef = list(map(newdef_map, self._definition))
            return self.parent()(kk, definition = newdef)

        def __lmul__(self, f):
            return self.__mul__(f)

        def __rmul__(self, f):
            return self.__mul__(f)

        def _acted_upon_(self, f, self_on_left=False):
            return self.__mul__(f)

        def _neg_(self):
            return self.parent()(self._k, definition = [(t[0],-t[1]) for t in self._definition])
