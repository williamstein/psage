#################################################################################
#
# (c) Copyright 2011 William Stein
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

cdef class Prime:
    cdef public long p, r
    cdef bint first
    cpdef long norm(self)
    cpdef bint is_split(self)
    cpdef bint is_inert(self)
    cpdef bint is_ramified(self)
