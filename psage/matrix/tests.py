from __future__ import absolute_import
#################################################################################
#
# (c) Copyright 2010 Fredrik Stroemberg
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

from sage.rings.complex_mpc import MPComplexField
from sage.all import ZZ,MatrixSpace,Matrix
from .linalg_complex_dense import *
from .matrix_complex_dense import *


def test_eigenvalues(prec=100,nmax=10,dimmax=10):
    r"""
    Test the eigenvalue computations for some random matrices.
    """
    F = MPComplexField(prec)
    dim = ZZ.random_element(2,dimmax)
    M = MatrixSpace(F,dim)
    for n in range(nmax):
	A,U,l=random_matrix_eigenvalues(F,dim)
	ev = A.eigenvalues()
	ev.sort(); l.sort()
	test = max([abs(ev[j]-l[j]) for j in range(len(ev))])
	assert test < A.eps()*100

##
## Helper functions
##

def random_matrix_eigenvalues(F,n):
    r"""
    Give a random matrix together with its eigenvalues.
    """
    l=list()
    M = MatrixSpace(F,n)
    U = Matrix(F,n)
    D = Matrix(F,n)
    for i in xrange(n):
        x = F.random_element()
        l.append(x)
        D[i,i]=x
    # now we need a unitary matrix:
    # make a random choice of vectors and use Gram-Schmidt to orthogonolize
    U = random_unitary_matrix(F,n)
    UT = U.transpose().conjugate()
    A = U*D*UT
    l.sort(cmp=my_abscmp)
    return   A,U,l
