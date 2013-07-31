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

class JacobiFormsModule_class(SageObject):
 
    def __init__(self, lattice, character):
        self._L = lattice
        self._h = character

    def lattice(self):
        return self._L

    def index(self):
        return self._L

    def character(self):
        return self._h

    def graded_component(self, k):
        return JacobiFormsSpace(k, self._L, self._h, self)

    def subspace (self, k):
        return self.graded_component(k)

    def __repr__(self):
        return "Module of Jacobi forms of index %s, character epsilon^%s.".format(self._L, self._h)
