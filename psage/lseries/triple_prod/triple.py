"""
Triple product L-series.
"""

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


from sage.all import ZZ

class TripleProductLseries(object):
    """
    A triple product `L`-series attached to three newforms of level `N`.
    """
    def __init__(self, N, f, g, h):
        """
        INPUT:
            - N -- squarefree integer
            - f -- an object such that for n>=0, we have
              f[n] = a_n(f) = n-th Fourier coefficient of
              a newform on Gamma_0(N).
            - g -- like f
            - h -- like f
        """
        self._N = ZZ(N)
        if not (self._N.is_squarefree() and self._N > 0):
            raise ValueError, "N (=%s) must be a squarefree positive integer"%self._N
        self._f = f
        self._g = g
        self._h = h

    def __repr__(self):
        """
        Return text representation of this triple product `L`-function.
        """
        return "Triple product L-function L(s,f,g,h) of three newforms on Gamma_0(%s)"%self._N

    def level(self):
        """
        Return the common level `N` of the newforms in the triple product.

        OUTPUT:
            - Integer
        
        EXAMPLES::
        """
        return self._N

    def newforms(self):
        """
        Return 3-tuple (f,g,h) of the data that defines the newforms
        in the triple product.

        OUTPUT:
            - 3-tuple
        
        EXAMPLES::
        """
        return self._f, self._g, self._h

    def charpoly(self, p):
        """
        Return the denominator as a polynomial in `X` (=`p^{-s}`) of the local
        factor at `p` of this `L`-series.

        The degree of the polynomial is ???? [[todo]].

        INPUT:
            - `p` -- prime

        OUTPUT:
            - polynomial in `X`

        EXAMPLES::

        """
        if self._N % p == 0:
            return self._charpoly_bad(p)
        else:
            return self._charpoly_good(p)

    def _charpoly_good(self, p):
        """
        Internal function that returns the local charpoly at a good prime.
        
        INPUT:
            - `p` -- prime

        OUTPUT:
            - polynomial in `X`

        EXAMPLES::

        """

    def _charpoly_bad(self, p):
        """
        Internal function that returns the local charpoly at a bad prime.
        
        INPUT:
            - `p` -- prime

        OUTPUT:
            - polynomial in `X`

        EXAMPLES::

        """
        X = self._gen
        a_p, b_p, c_p = self._f[p], self._g[p], self._h[p]
        return (1 - a_p*b_p*c_p * X) * (1 - a_p*b_p*c_p*p*X)**2
    
        
    
