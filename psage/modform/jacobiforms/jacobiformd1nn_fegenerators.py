"""
We provide methods to create Fourier expansions of (weak) Jacobi forms.
"""
from __future__ import division

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
from builtins import object
from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_lazyelement import \
                                        EquivariantMonoidPowerSeries_lazy
from psage.modform.jacobiforms.jacobiformd1nn_fourierexpansion import JacobiD1NNFourierExpansionModule
from psage.modform.jacobiforms.jacobiformd1nn_fourierexpansion import JacobiFormD1NNFilter
from sage.combinat.partition import number_of_partitions
from sage.libs.flint.fmpz_poly import Fmpz_poly  
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_function
from sage.misc.functional import isqrt
from sage.misc.misc import prod
from sage.modular.modform.constructor import ModularForms
from sage.modular.modform.element import ModularFormElement
from sage.modules.free_module_element import vector
from sage.rings.all import GF
from sage.arith.all import binomial, factorial
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.rational_field import QQ
from sage.structure.sage_object import SageObject
import operator

#===============================================================================
# _jacobi_forms_by_taylor_expansion_coords
#===============================================================================

_jacobi_forms_by_taylor_expansion_coords_cache = dict()

def _jacobi_forms_by_taylor_expansion_coords(index, weight, precision) :
    global _jacobi_forms_by_taylor_expansion_coords_cache
    
    key = (index, weight)
    try :
        return _jacobi_forms_by_taylor_expansion_coords_cache[key]
    except KeyError :
        if precision < (index - 1) // 4 + 1 :
            precision = (index - 1) // 4 + 1
        
        weak_forms = _all_weak_jacobi_forms_by_taylor_expansion(index, weight, precision)
        weak_index_matrix = matrix(ZZ, [ [ f[(n,r)] for n in range((index - 1) // 4 + 1)
                                                    for r in range(isqrt(4 * n * index) + 1,  index + 1) ]
                                         for f in weak_forms] )
        _jacobi_forms_by_taylor_expansion_coords_cache[key] = \
          weak_index_matrix.left_kernel().echelonized_basis()
          
        return _jacobi_forms_by_taylor_expansion_coords_cache[key]

#===============================================================================
# jacobi_forms_by_taylor_expansion
#===============================================================================

def jacobi_form_by_taylor_expansion(i, index, weight, precision) :
    r"""
    We first lift an echelon basis of elliptic modular forms to weak Jacobi forms.
    Then we return an echelon basis with respect enumeration of this first echelon
    basis of the '\ZZ'-submodule of Jacobi forms.
    """
    expansion_ring = JacobiD1NNFourierExpansionModule(ZZ, index)
        
    coefficients_factory = DelayedFactory_JacobiFormD1NN_taylor_expansion( i, index, weight, precision )
    
    return EquivariantMonoidPowerSeries_lazy(expansion_ring, precision, coefficients_factory.getcoeff)

#===============================================================================
# DelayedFactory_JacobiFormD1NN_taylor_expansion
#===============================================================================

class DelayedFactory_JacobiFormD1NN_taylor_expansion(object) :
    def __init__(self, i, index, weight, precision) :
        self.__i = i
        self.__index = index
        self.__precision = precision
        self.__weight = weight
        
        self.__series = None
    
    def getcoeff(self, key, **kwds) :
        (_, k) = key
        # for speed we ignore the character 
        if self.__series is None :
            self.__series = \
              sum( map( operator.mul,
                       _jacobi_forms_by_taylor_expansion_coords(self.__index, self.__weight, self.__precision)[self.__i],
                       _all_weak_jacobi_forms_by_taylor_expansion(self.__index, self.__weight, self.__precision) ) )

        try :
            return self.__series[k]
        except KeyError :
            return 0

#===============================================================================
# _theta_decomposition_indices
#===============================================================================

@cached_function
def _theta_decomposition_indices(index, weight) :
    r"""
    A list of possible indices of the echelon bases `M_k, S_{k+2}` etc. or `-1`
    if a component is zero." 
    """
    dims = [ ModularForms(1, weight).dimension()] + \
           [ ModularForms(1, weight + 2*i).dimension() - 1
             for i in range(1, index + 1) ]
           
    return [ (i,j) for (i,d) in enumerate(dims) for j in range(d)]

#===============================================================================
# _all_weak_jacobi_forms_by_taylor_expansion
#===============================================================================

@cached_function
def _all_weak_jacobi_forms_by_taylor_expansion(index, weight, precision) :
    """
    TESTS:
    
    We compute the Fourier expansion of a Jacobi form of weight `4` and index `2`.  This
    is denoted by ``d``.  Moreover, we span the space of all Jacobi forms of weight `8` and
    index `2`.  Multiplying the weight `4` by the Eisenstein series of weight `4` must
    yield an element of the weight `8` space.  Note that the multiplication is done using
    a polynomial ring, since no native multiplication for Jacobi forms is implemented.
    
    ::
    
        sage: from psage.modform.jacobiforms.jacobiformd1nn_fourierexpansion import *
        sage: from psage.modform.jacobiforms.jacobiformd1nn_fegenerators import _all_weak_jacobi_forms_by_taylor_expansion
        sage: from psage.modform.fourier_expansion_framework import *
        sage: prec = 20
        sage: d = JacobiFormD1NNFilter(prec, 2)
        sage: f = _all_weak_jacobi_forms_by_taylor_expansion(2, 4, prec)[0]
        sage: (g1,g2) = tuple(_all_weak_jacobi_forms_by_taylor_expansion(2, 8, prec))
        sage: em = ExpansionModule([g1, g2])
        sage: P.<q, zeta> = PolynomialRing(QQ, 2)
        sage: fp = sum(f[k] * q**k[0] * zeta**k[1] for k in JacobiFormD1NNFilter(prec, 2, reduced = False))
        sage: mf4 = ModularForms(1, 4).0.qexp(prec).polynomial()
        sage: h = mf4 * fp
        sage: eh = EquivariantMonoidPowerSeries(g1.parent(), {g1.parent().characters().gen(0) : dict((k, h[k[0],k[1]]) for k in d)}, d)
        sage: em.coordinates(eh.truncate(5), in_base_ring = False)
        [7/66, 4480]
        
    According to this we express ``eh`` in terms of the basis of the weight `8` space.
    
    ::
    
        sage: neh = eh - em.0.fourier_expansion() * 7 / 66 - em.1.fourier_expansion() * 4480
        sage: neh.coefficients()
        {(18, 0): 0, (12, 1): 0, (9, 1): 0, (3, 0): 0, (11, 2): 0, (8, 0): 0, (16, 2): 0, (2, 1): 0, (15, 1): 0, (6, 2): 0, (14, 0): 0, (19, 0): 0, (5, 1): 0, (7, 2): 0, (4, 0): 0, (1, 2): 0, (12, 2): 0, (9, 0): 0, (8, 1): 0, (18, 2): 0, (15, 0): 0, (17, 2): 0, (14, 1): 0, (11, 1): 0, (18, 1): 0, (5, 0): 0, (2, 2): 0, (10, 0): 0, (4, 1): 0, (1, 1): 0, (3, 2): 0, (0, 0): 0, (13, 2): 0, (8, 2): 0, (7, 1): 0, (6, 0): 0, (17, 1): 0, (11, 0): 0, (19, 2): 0, (16, 0): 0, (10, 1): 0, (4, 2): 0, (1, 0): 0, (14, 2): 0, (0, 1): 0, (13, 1): 0, (7, 0): 0, (15, 2): 0, (12, 0): 0, (9, 2): 0, (6, 1): 0, (3, 1): 0, (16, 1): 0, (2, 0): 0, (19, 1): 0, (5, 2): 0, (17, 0): 0, (13, 0): 0, (10, 2): 0}
    """
    modular_form_bases = \
      [ ModularForms(1, weight + 2*i) \
          .echelon_basis()[(0 if i == 0 else 1):]
        for i in range(index + 1) ]

    factory = JacobiFormD1NNFactory(precision, index)

    return [ weak_jacbi_form_by_taylor_expansion(
                i*[0] + [modular_form_bases[i][j]] + (index - i)*[0],
                precision, True, weight = weight,
                factory = factory )
             for (i,j) in _theta_decomposition_indices(index, weight) ]

#===============================================================================
# weak_jacbi_form_by_taylor_expansion
#===============================================================================

def weak_jacbi_form_by_taylor_expansion(fs, precision, is_integral = False, weight = None, factory = None) :
    if factory is None :
        factory = JacobiFormD1NNFactory(precision, len(fs) - 1)
            
    if is_integral :
        expansion_ring = JacobiD1NNFourierExpansionModule(ZZ, len(fs) - 1, True)
    else :
        expansion_ring = JacobiD1NNFourierExpansionModule(QQ, len(fs) - 1, True)

    f_exps = list()
    for (i,f) in enumerate(fs) :
        if f == 0 :
            f_exps.append(lambda p : 0)
        elif isinstance(f, ModularFormElement) :
            f_exps.append(f.qexp)
            
            if weight is None :
                weight = f.weight() - 2 * i
            else :
                if not weight == f.weight() - 2 * i :
                    ValueError("Weight of the i-th form must be k + 2*i.")
            if i != 0 and not f.is_cuspidal() :
                ValueError("All but the first form must be cusp forms.")    
        else :
            f_exps.append(f)

    if weight is None :
        raise ValueError("Either one element of fs must be a modular form or " + \
                         "the weight must be passed.")
    
    coefficients_factory = DelayedFactory_JacobiFormD1NN_taylor_expansion_weak( factory, f_exps, weight )
    return EquivariantMonoidPowerSeries_lazy(expansion_ring, expansion_ring.monoid().filter(precision), coefficients_factory.getcoeff)

#===============================================================================
# DelayedFactory_JacobiFormD1NN_taylor_expansion
#===============================================================================

class DelayedFactory_JacobiFormD1NN_taylor_expansion_weak(object) :
    def __init__(self, factory, fs, weight) :
        self.__factory = factory
        self.__fs = fs
        self.__weight = weight
        
        self.__series = None
    
    def getcoeff(self, key, **kwds) :
        (_, k) = key
        # for speed we ignore the character 
        if self.__series is None :
            self.__series = \
             self.__factory.by_taylor_expansion( self.__fs, self.__weight,
                                                    is_integral = True )
        
        try :
            return self.__series[k]
        except KeyError :
            return 0

#===============================================================================
# JacobiFormD1NNFactory
#===============================================================================

_jacobi_form_d1nn_factory_cache = dict()

def JacobiFormD1NNFactory(precision, m=None) :
    if not isinstance(precision, JacobiFormD1NNFilter) :
        if m is None :
            raise ValueError("if precision is not filter the index m must be passed.")
        precision = JacobiFormD1NNFilter(precision, m)
    
    global _jacobi_form_d1nn_factory_cache
        
    try :
        return _jacobi_form_d1nn_factory_cache[precision]
    except KeyError :
        tmp = JacobiFormD1NNFactory_class(precision)
        _jacobi_form_d1nn_factory_cache[precision] = tmp
        
        return tmp

#===============================================================================
# JacobiFormD1NNFactory_class
#===============================================================================

class JacobiFormD1NNFactory_class (SageObject) :
    
    def __init__(self, precision) :
        self.__precision = precision
        
        self.__power_series_ring_ZZ = PowerSeriesRing(ZZ, 'q')
        self.__power_series_ring = PowerSeriesRing(QQ, 'q')

    def index(self) :
        return self.__precision.jacobi_index()

    def power_series_ring(self) :
        return self.__power_series_ring
    
    def integral_power_series_ring(self) :
        return self.__power_series_ring_ZZ

    def _qexp_precision(self) :
        return self.__precision.index()

    def _set_theta_factors(self, theta_factors) :
        self.__theta_factors = theta_factors
    
    def _theta_factors(self, p = None) :
        r"""
        Return the factor `W^\# (\theta_0, .., \theta_{2m - 1})^{\mathrm{T}}` as a list.
        The `q`-expansion is shifted by `-(m + 1)(2*m + 1) / 24` which will be compensated
        for by the eta factor.
        """
        try :
            if p is None :
                return self.__theta_factors
            else :
                P = PowerSeriesRing(GF(p), 'q')
                
                return [list(map(P, facs)) for facs in self.__theta_factors] 

        except AttributeError :
            qexp_prec = self._qexp_precision()
            if p is None :
                PS = self.integral_power_series_ring()
            else :
                PS = PowerSeriesRing(GF(p), 'q')
            m = self.__precision.jacobi_index()
            
            twom = 2 * m
            frmsq = twom ** 2
            
            thetas = dict( ((i, j), dict())
                           for i in range(m + 1) for j in range(m + 1) )
            
            ## We want to calculate \hat \theta_{j,l} = sum_r (2 m r + j)**2l q**(m r**2 + j r).

            for r in range(0, isqrt((qexp_prec - 1 + m)//m) + 2) :
                for j in [0,m] :
                    fact = (twom*r + j)**2
                    coeff = 2
                    for l in range(0, m + 1) :
                        thetas[(j,l)][m*r**2 + r*j] = coeff
                        coeff = coeff * fact
            thetas[(0,0)][0] = 1
                    
            for r in range(0, isqrt((qexp_prec - 1 + m)//m) + 2) :
                for j in range(1, m) :
                    fact_p = (twom*r + j)**2
                    fact_m = (twom*r - j)**2
                    coeff_p = 2
                    coeff_m = 2
                    
                    for l in range(0, m + 1) :
                        thetas[(j,l)][m*r**2 + r*j] = coeff_p 
                        thetas[(j,l)][m*r**2 - r*j] = coeff_m
                        coeff_p = coeff_p * fact_p
                        coeff_m = coeff_m * fact_m 
                                    
            thetas = dict( ( k, PS(th).add_bigoh(qexp_prec) )
                           for (k,th) in thetas.items() )
            
            W = matrix(PS, m + 1, [ thetas[(j, l)]
                                    for j in range(m + 1) for l in range(m + 1) ])
            
            
            ## Since the adjoint of matrices with entries in a general ring
            ## is extremely slow for matrices of small size, we hard code the
            ## the cases `m = 2` and `m = 3`.  The expressions are obtained by
            ## computing the adjoint of a matrix with entries `w_{i,j}` in a
            ## polynomial algebra.
            if m == 2 and qexp_prec > 10**5 :
                adj00 =   W[1,1] * W[2,2] - W[2,1] * W[1,2]
                adj01 = - W[1,0] * W[2,2] + W[2,0] * W[1,2]
                adj02 =   W[1,0] * W[2,1] - W[2,0] * W[1,1]
                adj10 = - W[0,1] * W[2,2] + W[2,1] * W[0,2]
                adj11 =   W[0,0] * W[2,2] - W[2,0] * W[0,2]
                adj12 = - W[0,0] * W[2,1] + W[2,0] * W[0,1]
                adj20 =   W[0,1] * W[1,2] - W[1,1] * W[0,2]
                adj21 = - W[0,0] * W[1,2] + W[1,0] * W[0,2]
                adj22 =   W[0,0] * W[1,1] - W[1,0] * W[0,1]

                Wadj = matrix(PS, [ [adj00, adj01, adj02],
                                    [adj10, adj11, adj12],
                                    [adj20, adj21, adj22] ])
                  
            elif m == 3 and qexp_prec > 10**5 :
                adj00 = -W[0,2]*W[1,1]*W[2,0] + W[0,1]*W[1,2]*W[2,0] + W[0,2]*W[1,0]*W[2,1] - W[0,0]*W[1,2]*W[2,1] - W[0,1]*W[1,0]*W[2,2] + W[0,0]*W[1,1]*W[2,2]
                adj01 = -W[0,3]*W[1,1]*W[2,0] + W[0,1]*W[1,3]*W[2,0] + W[0,3]*W[1,0]*W[2,1] - W[0,0]*W[1,3]*W[2,1] - W[0,1]*W[1,0]*W[2,3] + W[0,0]*W[1,1]*W[2,3]
                adj02 = -W[0,3]*W[1,2]*W[2,0] + W[0,2]*W[1,3]*W[2,0] + W[0,3]*W[1,0]*W[2,2] - W[0,0]*W[1,3]*W[2,2] - W[0,2]*W[1,0]*W[2,3] + W[0,0]*W[1,2]*W[2,3]
                adj03 = -W[0,3]*W[1,2]*W[2,1] + W[0,2]*W[1,3]*W[2,1] + W[0,3]*W[1,1]*W[2,2] - W[0,1]*W[1,3]*W[2,2] - W[0,2]*W[1,1]*W[2,3] + W[0,1]*W[1,2]*W[2,3]

                adj10 = -W[0,2]*W[1,1]*W[3,0] + W[0,1]*W[1,2]*W[3,0] + W[0,2]*W[1,0]*W[3,1] - W[0,0]*W[1,2]*W[3,1] - W[0,1]*W[1,0]*W[3,2] + W[0,0]*W[1,1]*W[3,2]
                adj11 = -W[0,3]*W[1,1]*W[3,0] + W[0,1]*W[1,3]*W[3,0] + W[0,3]*W[1,0]*W[3,1] - W[0,0]*W[1,3]*W[3,1] - W[0,1]*W[1,0]*W[3,3] + W[0,0]*W[1,1]*W[3,3]
                adj12 = -W[0,3]*W[1,2]*W[3,0] + W[0,2]*W[1,3]*W[3,0] + W[0,3]*W[1,0]*W[3,2] - W[0,0]*W[1,3]*W[3,2] - W[0,2]*W[1,0]*W[3,3] + W[0,0]*W[1,2]*W[3,3]
                adj13 = -W[0,3]*W[1,2]*W[3,1] + W[0,2]*W[1,3]*W[3,1] + W[0,3]*W[1,1]*W[3,2] - W[0,1]*W[1,3]*W[3,2] - W[0,2]*W[1,1]*W[3,3] + W[0,1]*W[1,2]*W[3,3]

                adj20 = -W[0,2]*W[2,1]*W[3,0] + W[0,1]*W[2,2]*W[3,0] + W[0,2]*W[2,0]*W[3,1] - W[0,0]*W[2,2]*W[3,1] - W[0,1]*W[2,0]*W[3,2] + W[0,0]*W[2,1]*W[3,2]
                adj21 = -W[0,3]*W[2,1]*W[3,0] + W[0,1]*W[2,3]*W[3,0] + W[0,3]*W[2,0]*W[3,1] - W[0,0]*W[2,3]*W[3,1] - W[0,1]*W[2,0]*W[3,3] + W[0,0]*W[2,1]*W[3,3]
                adj22 = -W[0,3]*W[2,2]*W[3,0] + W[0,2]*W[2,3]*W[3,0] + W[0,3]*W[2,0]*W[3,2] - W[0,0]*W[2,3]*W[3,2] - W[0,2]*W[2,0]*W[3,3] + W[0,0]*W[2,2]*W[3,3]
                adj23 = -W[0,3]*W[2,2]*W[3,1] + W[0,2]*W[2,3]*W[3,1] + W[0,3]*W[2,1]*W[3,2] - W[0,1]*W[2,3]*W[3,2] - W[0,2]*W[2,1]*W[3,3] + W[0,1]*W[2,2]*W[3,3]

                adj30 = -W[1,2]*W[2,1]*W[3,0] + W[1,1]*W[2,2]*W[3,0] + W[1,2]*W[2,0]*W[3,1] - W[1,0]*W[2,2]*W[3,1] - W[1,1]*W[2,0]*W[3,2] + W[1,0]*W[2,1]*W[3,2]
                adj31 = -W[1,3]*W[2,1]*W[3,0] + W[1,1]*W[2,3]*W[3,0] + W[1,3]*W[2,0]*W[3,1] - W[1,0]*W[2,3]*W[3,1] - W[1,1]*W[2,0]*W[3,3] + W[1,0]*W[2,1]*W[3,3]
                adj32 = -W[1,3]*W[2,2]*W[3,0] + W[1,2]*W[2,3]*W[3,0] + W[1,3]*W[2,0]*W[3,2] - W[1,0]*W[2,3]*W[3,2] - W[1,2]*W[2,0]*W[3,3] + W[1,0]*W[2,2]*W[3,3]
                adj33 = -W[1,3]*W[2,2]*W[3,1] + W[1,2]*W[2,3]*W[3,1] + W[1,3]*W[2,1]*W[3,2] - W[1,1]*W[2,3]*W[3,2] - W[1,2]*W[2,1]*W[3,3] + W[1,1]*W[2,2]*W[3,3]

                Wadj = matrix(PS, [ [adj00, adj01, adj02, adj03],
                                    [adj10, adj11, adj12, adj13],
                                    [adj20, adj21, adj22, adj23],
                                    [adj30, adj31, adj32, adj33] ])
            else :
                Wadj = W.adjoint()
            
            theta_factors = [ [ Wadj[i,r] for i in range(m + 1) ]
                              for r in range(m + 1) ]
            
            if p is None :
                self.__theta_factors = theta_factors
                
            return theta_factors
    
    def _set_eta_factor(self, eta_factor) :
        self.__eta_factor = eta_factor
    
    def _eta_factor(self) :
        r"""
        The inverse determinant of `W`, which in these cases is always a negative
        power of the eta function. 
        """
        try :
            return self.__eta_factor
        except AttributeError :
            m = self.__precision.jacobi_index()
            pw = (m + 1) * (2 * m + 1)
            qexp_prec = self._qexp_precision()
            
            self.__eta_factor = self.integral_power_series_ring() \
                 ( [ number_of_partitions(n) for n in range(qexp_prec) ] ) \
                 .add_bigoh(qexp_prec) ** pw
                 
            return self.__eta_factor
 
    def by_taylor_expansion(self, fs, k, is_integral=False) :
        r"""
        We combine the theta decomposition and the heat operator as in [Sko].
        This yields a bijections of Jacobi forms of weight `k` and
        `M_k \times S_{k+2} \times .. \times S_{k+2m}`.
        """
        ## we introduce an abbreviations
        if is_integral :
            PS = self.integral_power_series_ring()
        else :
            PS = self.power_series_ring()
            
        if not len(fs) == self.__precision.jacobi_index() + 1 :
            raise ValueError("fs must be a list of m + 1 elliptic modular forms or their fourier expansion")
        
        qexp_prec = self._qexp_precision()
        if qexp_prec is None : # there are no forms below the precision
            return dict()
        
        f_divs = dict()
        for (i, f) in enumerate(fs) :
            f_divs[(i, 0)] = PS(f(qexp_prec), qexp_prec)
        
        if self.__precision.jacobi_index() == 1 :
            return self._by_taylor_expansion_m1(f_divs, k, is_integral)
        
        for i in range(self.__precision.jacobi_index() + 1) :
            for j in range(1, self.__precision.jacobi_index() - i + 1) :
                f_divs[(i,j)] = f_divs[(i, j - 1)].derivative().shift(1)
            
        phi_divs = list()
        for i in range(self.__precision.jacobi_index() + 1) :
            ## This is the formula in Skoruppas thesis. He uses d/ d tau instead of d / dz which yields
            ## a factor 4 m
            phi_divs.append( sum( f_divs[(j, i - j)] * (4 * self.__precision.jacobi_index())**i
                                  * binomial(i,j) / 2**i#2**(self.__precision.jacobi_index() - i + 1)
                                  * prod(2*(i - l) + 1 for l in range(1, i))
                                  / factorial(i + k + j - 1)
                                  * factorial(2*self.__precision.jacobi_index() + k - 1) 
                                  for j in range(i + 1) ) )
            
        phi_coeffs = dict()
        for r in range(self.__precision.jacobi_index() + 1) :
            series = sum( map(operator.mul, self._theta_factors()[r], phi_divs) )
            series = self._eta_factor() * series

            for n in range(qexp_prec) :
                phi_coeffs[(n, r)] = series[n]

        return phi_coeffs

    def _by_taylor_expansion_m1(self, f_divs, k, is_integral=False) :
        r"""
        This provides special, faster code in the Jacobi index `1` case.
        """
        if is_integral :
            PS = self.integral_power_series_ring()
        else :
            PS = self.power_series_ring()
            
        qexp_prec = self._qexp_precision()
        
        
        fderiv = f_divs[(0,0)].derivative().shift(1)
        f = f_divs[(0,0)] * Integer(k/2)
        gfderiv = f_divs[(1,0)] - fderiv
        
        ab_prec = isqrt(qexp_prec + 1)
        a1dict = dict(); a0dict = dict()
        b1dict = dict(); b0dict = dict()
    
        for t in range(1, ab_prec + 1) :
            tmp = t**2
            a1dict[tmp] = -8*tmp
            b1dict[tmp] = -2
        
            tmp += t
            a0dict[tmp] = 8*tmp + 2
            b0dict[tmp] = 2
        b1dict[0] = -1
        a0dict[0] = 2; b0dict[0] = 2 
        
        a1 = PS(a1dict); b1 = PS(b1dict)
        a0 = PS(a0dict); b0 = PS(b0dict)

        Ifg0 = (self._eta_factor() * (f*a0 + gfderiv*b0)).list()
        Ifg1 = (self._eta_factor() * (f*a1 + gfderiv*b1)).list()

        if len(Ifg0) < qexp_prec :
            Ifg0 += [0]*(qexp_prec - len(Ifg0))
        if len(Ifg1) < qexp_prec :
            Ifg1 += [0]*(qexp_prec - len(Ifg1))
                    
        Cphi = dict([(0,0)])
        for i in range(qexp_prec) :
            Cphi[-4*i] = Ifg0[i]
            Cphi[1-4*i] = Ifg1[i]

        del Ifg0[:], Ifg1[:]

        phi_coeffs = dict()
        m = self.__precision.jacobi_index()
        for r in range(2 * self.__precision.jacobi_index()) :
            for n in range(qexp_prec) :
                k = 4 * m * n - r**2
                if k >= 0 :
                    phi_coeffs[(n, r)] = Cphi[-k]
                               
        return phi_coeffs
