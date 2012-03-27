r"""
Construction of Fourier expansions of Siegel modular forms of arbitrary genus.

AUTHORS:

- Martin Raum (2009 - 05 - 11) Initial version
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

from sage.matrix.constructor import diagonal_matrix
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_method
from sage.misc.functional import isqrt
from sage.misc.misc import prod
from sage.modules.all import vector
from sage.modules.free_module import FreeModule
from sage.quadratic_forms.quadratic_form import QuadraticForm
from sage.rings.arith import gcd, fundamental_discriminant
from sage.rings.all import ZZ, Integer, GF
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.structure.sage_object import SageObject
from psage.modform.paramodularforms.siegelmodularformgn_fourierexpansion import SiegelModularFormGnFourierExpansionRing,\
    SiegelModularFormGnFilter_diagonal_lll
import itertools
from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import EquivariantMonoidPowerSeries_abstract

#===============================================================================
# DukeImamogulLift
#===============================================================================

def DukeImamogluLift(f, f_weight, precision) :
    r"""
    INPUT:
    
    - `f` -- Fourier expansion of a Jacobi form of index `1` and weight ``k + 1 - precision.genus / 2``.
        
    NOTE:
    
        We follow Kohnen's paper "Lifting modular forms of half-integral weight to Siegel
        modular forms of even genus", Section 5, Theorem 1. 
    """
    if not isinstance(precision, SiegelModularFormGnFilter_diagonal_lll) :
        raise TypeError, "precision must be an instance of SiegelModularFormGnFilter_diagonal_lll"

    ## TODO: Neatly check the input type.
    half_integral_weight = not isinstance(f, EquivariantMonoidPowerSeries_abstract)
    if half_integral_weight :
        f_k = Integer(f_weight - Integer(1)/2)
    else :
        f_k = f_weight
    
    fe_ring = SiegelModularFormGnFourierExpansionRing( 
                            f.parent().base_ring() if half_integral_weight else f.parent().coefficient_domain(),
                            precision.genus() )
    
    form = fe_ring(DukeImamogluLiftFactory().duke_imamoglu_lift( f if half_integral_weight else f.coefficients(),
                                                                 f_k, precision, half_integral_weight))
    form._set_precision(precision)
    
    return form


#===============================================================================
# DukeImamogluLiftFactory
#===============================================================================

_duke_imamoglu_lift_factory_cache = None

def DukeImamogluLiftFactory() :
    global _duke_imamoglu_lift_factory_cache
    
    if _duke_imamoglu_lift_factory_cache is None :
        _duke_imamoglu_lift_factory_cache = DukeImamogluLiftFactory_class()
        
    return _duke_imamoglu_lift_factory_cache

#===============================================================================
# DukeImamogluLiftFactory_class
#===============================================================================

class DukeImamogluLiftFactory_class (SageObject) :
    
    @cached_method
    def _kohnen_phi(self, a, t) :
        ## We use a modified power of a, namely for each prime factor p of a
        ## we use p**floor(v_p(a)/2). This is compensated for in the routine
        ## of \rho
        a_modif = 1
        for (p,e) in a.factor() :
            a_modif = a_modif * p**(e // 2)
        
        res = 0
        for dsq in filter(lambda d: d.is_square(), a.divisors()) :
            d = isqrt(dsq)
            
            for g_diag in itertools.ifilter( lambda diag: prod(diag) == d,
                                  itertools.product(*[d.divisors() for _ in xrange(t.nrows())]) ) :
                for subents in itertools.product(*[xrange(r) for (j,r) in enumerate(g_diag) for _ in xrange(j)]) :
                    columns = [subents[(j * (j - 1)) // 2:(j * (j + 1)) // 2] for j in xrange(t.nrows())]
                    g = diagonal_matrix(list(g_diag))
                    for j in xrange(t.nrows()) :
                        for i in xrange(j) :
                            g[i,j] = columns[j][i]
                         
                    ginv = g.inverse()   
                    tg = ginv.transpose() * t * ginv
                    try :
                        tg= matrix(ZZ, tg)
                    except :
                        continue
                    if any(tg[i,i] % 2 == 1 for i in xrange(tg.nrows())) :
                        continue
                    
                    tg.set_immutable()
                        
                    res = res + self._kohnen_rho(tg, a // dsq)
        
        return a_modif * res

    @cached_method
    def _kohnen_rho(self, t, a) :
        dt = (-1)**(t.nrows() // 2) * t.det()
        d = fundamental_discriminant(dt)
        eps = isqrt(dt // d)
        
        res = 1
        for (p,e) in gcd(a, eps).factor() :
            pol = self._kohnen_rho_polynomial(t, p).coefficients()
            if e < len(pol) :
                res = res * pol[e]
            
        return res
    
    @cached_method
    def _kohnen_rho_polynomial(self, t, p) :
        P = PolynomialRing(ZZ, 'x')
        x = P.gen(0)
        
        #tc = self._minimal_isotropic_subspace_complement_mod_p(t, p)
        
        if t[0,0] != 0 :
            return (1 - x**2)
        
        for i in xrange(t.nrows()) :
            if t[i,i] != 0 :
                break
        
        if i % 2 == 0 :
            ## Since we have semi definite indices Kohnen's lambda simplifies
            lambda_p = 1 if t == 0 else -1
            ## For we multiply x in the second factor by p**(1/2) to make all coefficients rational.
            return (1 - x**2) * (1 + lambda_p * p**((i + 1) // 2) * x) * \
                   prod(1 - p**(2*j - 1) * x**2 for j in xrange(1, i // 2))
        else :
            return (1 - x**2) * \
                   prod(1 - p**(2*j - 1) * x**2 for j in xrange(1, i // 2))
    
    #===========================================================================
    # ## TODO: Speadup by implementation in Cython and using enumeration mod p**n and fast
    # ##       evaluation via shifts
    # def _minimal_isotropic_subspace_complement_mod_p(self, t, p, cur_form_subspace = None, cur_isotropic_space = None) :
    #    if not t.base_ring() is GF(p) :
    #        t = matrix(GF(p), t)
    #    #q = QuadraticForm(t)
    #    
    #    if cur_form_subspace is None :
    #        cur_form_subspace = FreeModule(GF(p), t.nrows())
    #    if cur_isotropic_space is None :
    #        cur_isotropic_space = cur_form_subspace.ambient_module().submodule([])
    #    
    #    ## We find a zero length vector. In the orthogonal complement we can then proceed.
    #    
    #    for v in cur_form_subspace :
    #        if v not in cur_isotropic_space \
    #          and (v * t * v).is_zero() :
    #            cur_form_subspace = cur_form_subspace.intersection(matrix(GF(p), v * t).right_kernel())
    #            cur_isotropic_space = cur_isotropic_space + cur_isotropic_space.ambient_module().submodule([v])
    #            
    #            return self._minimal_isotropic_subspace_complement_mod_p(t, p, cur_form_subspace, cur_isotropic_space)
    #    # return complement
    #    complement_coords =   set(xrange(cur_form_subspace.ambient_module().rank())) \
    #                        - set(cur_isotropic_space.echelonized_basis_matrix().pivots())
    #    complement_basis = [ vector(GF(p), [0]*i + [1] + [0]*(cur_form_subspace.ambient_module().rank() - i - 1))
    #                         for i in complement_coords ]
    #    
    #    return matrix(GF(p), [[u * t * v for u in complement_basis] for v in complement_basis])
    #===========================================================================
    
    def duke_imamoglu_lift(self, f, f_k, precision, half_integral_weight = False) :
        """
        INPUT:
        
        - ``half_integral_weight``   -- If ``False`` we assume that `f` is the Fourier expansion of a
                                        Jacobi form. Otherwise we assume it is the Fourier expansion
                                        of a half integral weight elliptic modular form.
        """
        
        if half_integral_weight :
            coeff_index = lambda d : d
        else :
            coeff_index = lambda d : ((d + (-d % 4)) // 4, (-d) % 4)
        
        coeffs = dict()
        
        for t in precision.iter_positive_forms() :
            dt = (-1)**(precision.genus() // 2) * t.det()
            d = fundamental_discriminant(dt)
            eps = Integer(isqrt(dt / d))
    
            coeffs[t] = 0 
            for a in eps.divisors() :
                d_a = abs(d * (eps // a)**2)
                 
                coeffs[t] = coeffs[t] + a**(f_k - 1) * self._kohnen_phi(a, t) \
                                        * f[coeff_index(d_a)]

        return coeffs
    