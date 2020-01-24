r"""
We provide methods to create Fourier expansions of (weak) Jacobi forms `\mathrm{mod} p`.
"""

#===============================================================================
# 
# Copyright (C) 2010 Martin Raum
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================

from builtins import map
from builtins import range
from sage.misc.all import prod
from sage.rings.all import PowerSeriesRing, GF
from sage.rings.all import binomial, factorial
from sage.structure.sage_object import SageObject
import operator

#===============================================================================
# JacobiFormD1NNModularFactory_class
#===============================================================================

class JacobiFormD1NNModularFactory_class (SageObject) :
    
    def __init__(self, precision, p) :
        self.__precision = precision
        self.__power_series_ring = PowerSeriesRing(GF(p), 'q')
        self.__p = int(p)

    def index(self) :
        return self.__precision.jacobi_index()

    def power_series_ring(self) :
        return self.__power_series_ring
    
    def _qexp_precision(self) :
        return self.__precision.index()

    def _set_theta_factors(self, theta_factors) :
        self.__theta_factors = theta_factors
    
    def _theta_factors(self) :
        try :
            return self.__theta_factors
        except AttributeError :
            raise RuntimeError( "Theta factors have to be set first" )
    
    def _set_eta_factor(self, eta_factor) :
        self.__eta_factor = eta_factor
    
    def _eta_factor(self) :
        try :
            return self.__eta_factor
        except AttributeError :
            raise RuntimeError( "Eta factor have to be set first" )
 
    def by_taylor_expansion(self, fs, k) :
        """
        We combine the theta decomposition and the heat operator as in [Sko].
        This yields a bijections of Jacobi forms of weight `k` and
        `M_k \times S_{k+2} \times .. \times S_{k+2m}`.
        
        NOTE:

            To make ``phi_divs`` integral we introduce an extra factor
            `2^{\mathrm{index}} * \mathrm{factorial}(k + 2 * \mathrm{index} - 1)`.
        """
        ## we introduce an abbreviations
        p = self.__p
        PS = self.power_series_ring()
            
        if not len(fs) == self.__precision.jacobi_index() + 1 :
            raise ValueError("fs must be a list of m + 1 elliptic modular forms or their fourier expansion")
        
        qexp_prec = self._qexp_precision()
        if qexp_prec is None : # there are no forms below the precision
            return dict()
        
        f_divs = dict()
        for (i, f) in enumerate(fs) :
            f_divs[(i, 0)] = PS(f(qexp_prec), qexp_prec)
                
        for i in range(self.__precision.jacobi_index() + 1) :
            for j in range(1, self.__precision.jacobi_index() - i + 1) :
                f_divs[(i,j)] = f_divs[(i, j - 1)].derivative().shift(1)
            
        phi_divs = list()
        for i in range(self.__precision.jacobi_index() + 1) :
            ## This is the formula in Skoruppas thesis. He uses d/ d tau instead of d / dz which yields
            ## a factor 4 m
            phi_divs.append( sum( f_divs[(j, i - j)] * (4 * self.__precision.jacobi_index())**i
                                  * binomial(i,j) * ( 2**self.index() // 2**i)
                                  * prod(2*(i - l) + 1 for l in range(1, i))
                                  * (factorial(k + 2*self.index() - 1) // factorial(i + k + j - 1))
                                  * factorial(2*self.__precision.jacobi_index() + k - 1)
                                  for j in range(i + 1) ) )
            
        phi_coeffs = dict()
        for r in range(self.index() + 1) :
            series = sum( map(operator.mul, self._theta_factors()[r], phi_divs) )
            series = self._eta_factor() * series

            for n in range(qexp_prec) :
                phi_coeffs[(n, r)] = int(series[n].lift()) % p

        return phi_coeffs
