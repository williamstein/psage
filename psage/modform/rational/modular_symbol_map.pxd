#################################################################################
#
# (c) Copyright 2010 William Stein
#
#  This file is part of PSAGE
#
#  PSAGE is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  PSAGE is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################

cdef enum:
    MAX_CONTFRAC = 100
    MAX_DEG = 10000

from sage.modular.modsym.p1list cimport P1List

cdef class ModularSymbolMap:
    cdef long d, N
    cdef public long denom
    cdef long* X  # coefficients of linear map from P^1 to Q^d.
    cdef public object C
    cdef P1List P1
    cdef int evaluate(self, long v[MAX_DEG], long a, long b) except -1
