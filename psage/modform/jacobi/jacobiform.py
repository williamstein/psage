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
 Implements classes for (individual) Jacobi forms.

 AUTHORS::

  - Stephan Ehlen
  - Nils Skoruppa
  - Fredrik Strömberg

 EXAMPLES::
"""

from sage.all import SageObject, PowerSeriesRing, QQbar, Infinity
#from sage.modules.free_module_element import FreeModuleElement
from sage.modular.modform.element import ModularFormElement

class JacobiForm_class(SageObject):
 
    def __init__(self, weight, lattice, character, prec=10):
        print "Init of JacobiForm_class"
        self._k = weight
        self._L = lattice
        self._h = character
        self._prec = prec
        # _feparent: this is the ambient space in which the Fourier expansion of self lives
        # we have to think about this...
        self._feparent = PowerSeriesRing(QQbar, 'q,z')

    def fourier_expansion(self, prec=None):
        r"""
          Returns the Fourier expansion of self.
        """
        if prec == None:
            prec = self._prec
        raise NotImplementedError("This method is currently not implemented. It should be overriden by the specific subclasses.")

    def weight(self):
        return self._k

    def lattice(self):
        return self._L

    def index(self):
        return self._L

    def character(self):
        return self._h

    def __repr__(self):
        return "Jacobi form of weight {0}, index {1} and character epsilon^{2}".format(self._k, self._L, self._h)

class JacobiForm_space_element(JacobiForm_class, ModuleElement):

    def __init__(self, ambient_space, prec=10):
        print "Init of JacobiForm_space_element"
        self._ambient_space = J = ambient_space
        if ambient_space != None:
            try:
                self._weight = weight = J.weight()
                self._lattice = lattice = J.lattice()
                self._character = character = J.character()
            except:
                raise ValueError("ambient_space has to be a JacobiFormsSpace")
            super(JacobiForm_space_element, self).__init__(weight, lattice, character, prec)
        else:
            super(JacobiForm_space_element, self).__init__(-Infinity, None, None, prec)

    def ambient_space(self):
        r"""
          Returns the ambient space of Jacobi forms.
        """
        if self._ambient_space == None:
            self._ambient_space = JacobiFormsSpace(weight, lattice, character)
        return self._ambient_space

class JacobiForm_module_element(JacobiForm_class, ModuleElement):

    def __init__(self, weight, ambient_module, construction=None, prec=10):
        try:
            if not weight == -Infinity:
                ambient_space = ambient_module.graded_component(weight)
                super(JacobiForm_module_element, self).__init__(ambient_space, prec)
            else:
                ambient_space = None
                print "here"
                super(JacobiForm_module_element, self).__init__(ambient_space, prec)
                print "here2"
                L = ambient_module.lattice()
                c = ambient_module.character()
                print "here3"
                #super(JacobiForm_module_element, self).__init__(weight, L, c, prec)
        except:
            raise ValueError("ambient_module has to be a JacobiFormsModule_class.")
        self._ambient_module = ambient_module
        self._construction = construction

    def ambient_module(self):
        r"""
          Returns the ambient module of Jacobi forms.
        """
        if self._ambient_module == None:
            self._ambient_module = JacobiFormsSpace(weight, lattice, character)
        return self._ambient_module

    def space_element(self):
        r"""
          Returns self as an element of the ambient space (forget the module).
          Does this make sense?
        """
        raise NotImplementedError()

    def construction(self):
        return self._construction

    def __add__(self, phi):
        return JacobiForm_module_element(self._k + phi.weight(), self._ambient_module, self._construction + phi.construction())

    def __rmul__(self, other):
        if isinstance(other, ModularFormElement):
            f = other
            if f.level() != 1:
                raise ValueError("Can only multiply by modular forms of level 1.")
            construction = map(lambda t: (f*t[0], t[1]), self._construction)
            return JacobiForm_module_element(self._k + f.weight(), self._ambient_module, construction)
        elif isinstance(other, JacobiForm_module_element):
            raise NotImplementedError("Currently, only multiplication by modular forms for SL_2(Z) is implemented")

class JacobiForm_module_element_zero(JacobiForm_module_element):

    def __init__(self, ambient_module):
        super(JacobiForm_module_element_zero, self).__init__(-Infinity, ambient_module, construction = [(0,0)])

    def fourier_expansion(self, prec=None):
        return self._feparent.zero()
        
class JacobiForm_space_element_zero(JacobiForm_space_element):

    def fourier_expansion(self, prec=None):
        return self._feparent.zero()

    
