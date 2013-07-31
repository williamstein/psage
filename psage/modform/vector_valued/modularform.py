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
 Implements a base class for vector valued modular forms.

 AUTHORS::

  - Stephan Ehlen (initial version)

 EXAMPLES::
"""

class VectorValuedModularForm_generic_class(SageObject):
 
    def __init__(self, weight, module, ambient_space=None, ambient_module=None):
        self._k = weight
        self._M = module
        self._ambient_space = ambient_space
        self._ambient_module = ambient_module
        
    def fourier_expansion(self):
        r"""
          Returns the Fourier expansion of self.
        """
        return NotImplementedError("This method is currently not implemented. It should be overriden by the specific subclasses.")

    def ambient_space(self):
        r"""
          Returns the ambient space of vector valued modular forms.
        """
        return NotImplementedError("This method is currently not implemented. It should be overriden by the specific subclasses.")

    def ambient_module(self):
        r"""
          Returns the ambient module of vector valued modular forms.
        """
        return NotImplementedError("This method is currently not implemented. It should be overriden by the specific subclasses.")

    def weight(self):
        return self._k

    def module(self):
        return self._M

    def __repr__(self):
        return "Modular form of weight %s for the representation given by %s".format(self._k, self._M)
