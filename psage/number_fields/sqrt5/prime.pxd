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

cdef class PrimesOfBoundedNorm:
    cdef public long bound

    # size of the array of primes -- it has this many entries
    cdef long table_size

    # fast *UNSAFE* access to the i-th prime and root through these C arrays
    cdef long* prime
    cdef long* root

    # fast SAFE access to the i-th prime and the i-th root
    cdef long get_prime(self, Py_ssize_t i) except -1
    cdef long get_root(self, Py_ssize_t i) except 1099511627776

    # only used internally
    cdef long _allocate_memory(self) except -1
    # only used internally
    cdef int _reallocate_memory(self, long max_size) except -1

    
    
