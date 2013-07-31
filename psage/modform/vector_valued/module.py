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
 Implements a class for a module of vector valued forms (the direct sum of the spaces of all weights).

 AUTHORS::

  - Stephan Ehlen (Initial version)

 EXAMPLES::
"""

class VectorValuedModularFormsModule_generic_class(SageObject):
 
    def __init__(self, module):
        self._M = module

    def graded_component(self, k):
        return VectorValuedModularFormsSpace(k, self._M, self)

    def subspace (self, k):
        return self.graded_component(k)

    def __repr__(self):
        return "Module of vector valued modular forms for the representation given by the module %s.".format(self._M)
    
