r"""
Generator functions for the fourier expansion of paramodular modular forms.
Therefore apply Gritsenko's lift to a Jacobi form of index N. 

AUTHORS:

- Martin Raum (2010 - 04 - 09) Initial version.
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

from psage.modform.jacobiforms.jacobiformd1nn_fegenerators import JacobiFormD1NNFactory
from psage.modform.jacobiforms.jacobiformd1nn_types import JacobiFormsD1NN,\
    JacobiFormD1NN_Gamma
from psage.modform.paramodularforms.paramodularformd2_fourierexpansion import ParamodularFormD2Filter_discriminant
from psage.modform.paramodularforms.paramodularformd2_fourierexpansion import ParamodularFormD2FourierExpansionRing
from psage.modform.paramodularforms.paramodularformd2_fourierexpansion_cython import apply_GL_to_form
from psage.modform.paramodularforms.siegelmodularformg2_types import SiegelModularFormsG2,\
    SiegelModularFormG2_Classical_Gamma
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_function
from sage.misc.cachefunc import cached_method
from sage.arith.all import bernoulli, sigma
from sage.arith.all import random_prime
from sage.rings.all import GF
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.structure.all import Sequence
from sage.structure.sage_object import SageObject
from psage.modform.paramodularforms.siegelmodularformg2_misc_cython import divisor_dict

#===============================================================================
# _minimal_paramodular_precision
#===============================================================================

def _minimal_paramodular_precision(level, weight) :
    ## This is the bound, that Poor and Yuen provide in Paramodular
    ## cusp forms, Corollary 5.7

    if Integer(level).is_prime() :
        return ParamodularFormD2Filter_discriminant( 2 * weight**2 / Integer(225)
                                                     * (1 + level**2)**2 / Integer(1 + level)**2,
                                                     2 )
    else :
        raise NotImplementedError( "Minimal bound can only be provided for prime level" )

def symmetrise_siegel_modular_form(f, level, weight) :
    """
    ALGORITHM:
        We use the Hecke operators \delta_p and [1,level,1,1/level], where the last one
        can be chosen, since we assume that f is invariant with respect to the full
        Siegel modular group. 
    """
    if not Integer(level).is_prime() :
        raise NotImplementedError( "Symmetrisation is only implemented for prime levels" )
    fr = ParamodularFormD2FourierExpansionRing(f.parent().coefficient_domain(), level)
    
    precision = ParamodularFormD2Filter_discriminant(f.precision().discriminant() // level**2, level)

    coeffs = dict()
    for (k,l) in precision :
        kp = apply_GL_to_form(precision._p1list()[l], k)
        coeffs[(k,l)] = f[kp]

    g1 = fr( {fr.characters().one_element() : coeffs} )
    g1._set_precision(precision)
    g1._cleanup_coefficients()

    coeffs = dict()
    for (k,l) in precision :
        kp = apply_GL_to_form(precision._p1list()[l], k)
        if kp[1] % level == 0 and kp[2] % level**2 == 0 :
            coeffs[(k,l)] = f[(kp[0], kp[1] // level, kp[2] // level**2)]

    g2 = fr( {fr.characters().one_element() : coeffs} )
    g2._set_precision(precision)
    g2._cleanup_coefficients()

    return g1 + level**(weight-1) * g2

#===============================================================================
# symmetrised_siegel_modular_forms
#===============================================================================

def symmetrised_siegel_modular_forms(level, weight, precision) :
    min_precision = _minimal_paramodular_precision(level, weight)
    if precision is None :
        precision = min_precision
    else :
        precision = max(precision, min_precision)

    sr = SiegelModularFormsG2( QQ, SiegelModularFormG2_Classical_Gamma(),
                               precision._enveloping_discriminant_bound() * level**2 )
    return map( lambda b: symmetrise_siegel_modular_form(b.fourier_expansion(), level, weight),
                sr.graded_submodule(weight).basis() )
    
#===============================================================================
# gritsenko_lift_fourier_expansion
#===============================================================================

_gritsenko_lift_factory_cache = dict()

def gritsenko_lift_fourier_expansion(f, precision, is_integral = False) :
    """
    Given initial data return the Fourier expansion of a Gritsenko lift
    for given level and weight.
    
    INPUT:
    
        - `f` -- A Jacobi form.
        
    OUTPUT:
    
        A dictionary of coefficients.
        
    TODO:
    
        Make this return lazy power series.
    """
    global _gritsenko_lift_factory_cache
    N = f.parent().type().index()
    
    try :
        fact = _gritsenko_lift_factory_cache[(precision, N)]
    except KeyError :
        fact = ParamodularFormD2Factory(precision, N)
        _gritsenko_lift_factory_cache[(precision, N)] = fact
    
    
    fe = ParamodularFormD2FourierExpansionRing(ZZ, N)(
           fact.gritsenko_lift(f.fourier_expansion(), f.parent().type().weight(), is_integral) )
    fe._set_precision(precision)
    return fe

#===============================================================================
# gritsenko_lift_subspace
#===============================================================================

@cached_function
def gritsenko_lift_subspace(N, weight, precision) :
    """
    Return a list of data, which can be used by the Fourier expansion
    generator, for the space of Gritsenko lifts.
    
    INPUT:
        - `N`           -- the level of the paramodular group 
        - ``weight``    -- the weight of the space, that we consider
        - ``precision`` -- A precision class 
        
    NOTE:
        The lifts returned have to have integral Fourier coefficients.
    """
    jf = JacobiFormsD1NN( QQ, JacobiFormD1NN_Gamma(N, weight),
                          (4 * N**2 + precision._enveloping_discriminant_bound() - 1)//(4 * N) + 1)

    return Sequence( [ gritsenko_lift_fourier_expansion( g if i != 0 else (bernoulli(weight) / (2 * weight)).denominator() * g,
                                                         precision, True )
                       for (i,g) in enumerate(jf.gens()) ],
                       universe = ParamodularFormD2FourierExpansionRing(ZZ, N), immutable = True,
                       check = False )
    
#===============================================================================
# gritsenko_products
#===============================================================================

def gritsenko_products(N, weight, dimension, precision = None) :
    """
    INPUT:
        - `N`           -- the level of the paramodular group.
        - ``weight``    -- the weight of the space, that we consider.
        - ``precision`` -- the precision to be used for the underlying
                           Fourier expansions.
    """
    min_precision = _minimal_paramodular_precision(N, weight)
    if precision is None :
        precision = min_precision
    else :
        precision = max(precision, min_precision)
    
    # - Lifts betrachten
    # - Zerlegung des Gewichts benutzen, dort Lifts nach und nach
    #   multiplizieren, um jeweils den Rang der unterliegenden Matrix
    #   zu bestimmen. Dabei mod p mit zufaelligem p arbeiten
        
    fourier_expansions = list(gritsenko_lift_subspace(N, weight, precision))
    products = [ [(weight, i)] for i in xrange(len(fourier_expansions)) ]

    fourier_indices = list(precision)
    fourier_matrix = matrix(ZZ, [ [e[k] for k in fourier_indices]
                                  for e in fourier_expansions] )
    
    for k in xrange(min((weight // 2) - ((weight // 2) % 2), weight - 2), 1, -2) :
        space_k = list(enumerate(gritsenko_lift_subspace(N, k, precision)))
        space_wk = list(enumerate(gritsenko_lift_subspace(N, weight - k, precision)))

        if len(space_wk) == 0 :
            continue
        
        for (i,f) in space_k :
            
            for (j,g) in space_wk :
                if fourier_matrix.nrows() == dimension :
                    return products, fourier_expansions
                
                cur_fe = f * g
                # We are using lazy rank checks
                for p in [random_prime(10**10) for _ in range(2)] :
                    fourier_matrix_cur = matrix(GF(p), fourier_matrix.nrows() + 1,
                                                       fourier_matrix.ncols())
                    fourier_matrix_cur.set_block(0, 0, fourier_matrix)
                    fourier_matrix_cur.set_block( fourier_matrix.nrows(), 0,
                                                  matrix(GF(p), 1, [cur_fe[fi] for fi in fourier_indices]) )
                    
                    if fourier_matrix_cur.rank() == fourier_matrix_cur.nrows() :
                        break
                else :
                    continue
                
                products.append( [(k, i), (weight - k, j)] )
                fourier_expansions.append(cur_fe)
                
                fourier_matrix_cur = matrix( ZZ, fourier_matrix.nrows() + 1,
                                                 fourier_matrix.ncols() )
                fourier_matrix_cur.set_block(0, 0, fourier_matrix)
                fourier_matrix_cur.set_block( fourier_matrix.nrows(), 0,
                                              matrix(ZZ, 1, [cur_fe[fi] for fi in fourier_indices]) )
                fourier_matrix = fourier_matrix_cur
    
    return products, fourier_expansions

#===============================================================================
# ParamodularFormD2Factory
#===============================================================================

_paramodular_form_d2_factory_cache = dict()

def ParamodularFormD2Factory(precision, level) :
    global _paramodular_form_d2_factory_cache

    if not isinstance(precision, ParamodularFormD2Filter_discriminant) :
        precision = ParamodularFormD2Filter_discriminant(precision, level)
        
    try :
        return _paramodular_form_d2_factory_cache[precision]
    except KeyError :
        tmp = ParamodularFormD2Factory_class(precision)
        _paramodular_form_d2_factory_cache[precision] = tmp
        
        return tmp
    
class ParamodularFormD2Factory_class ( SageObject ) :
    
    def __init__(self, precision) :
        self.__precision = precision
        self.__jacobi_factory = None

    def precision(self) :
        return self.__precision

    def level(self) :
        return self.__precision.level()

    def _discriminant_bound(self) :
        return self.__precision._contained_discriminant_bound()
    
    @cached_method
    def _divisor_dict(self) :
        return divisor_dict(self._discriminant_bound())
    
    def gritsenko_lift(self, f, k, is_integral = False) :
        """
        INPUT:
            - `f` -- the Fourier expansion of a Jacobi form as a dictionary
        """
        N = self.level()
        frN = 4 * N
        p1list = self.precision()._p1list()

        divisor_dict = self._divisor_dict()
        
        ##TODO: Precalculate disc_coeffs or at least use a recursive definition
        disc_coeffs = dict()
        coeffs = dict()
        f_index = lambda d,b : ((d + b**2)//frN, b)
        
        for (((a,b,c),l), eps, disc) in self.__precision._iter_positive_forms_with_content_and_discriminant() :
            (_,bp,_) = apply_GL_to_form(p1list[l], (a,b,c))
            try :
                coeffs[((a,b,c),l)] = disc_coeffs[(disc, bp, eps)]
            except KeyError :
                disc_coeffs[(disc, bp, eps)] = \
                    sum(   t**(k-1) * f[ f_index(disc//t**2, (bp // t) % (2 * N)) ]
                           for t in divisor_dict[eps] )
 
                if disc_coeffs[(disc, bp, eps)]  != 0 :
                    coeffs[((a,b,c),l)] = disc_coeffs[(disc, bp, eps)]

        for ((a,b,c), l) in self.__precision.iter_indefinite_forms() :
            if l == 0 :
                coeffs[((a,b,c),l)] = ( sigma(c, k-1) * f[(0,0)]
                                        if c != 0 else 
                                         ( Integer(-bernoulli(k) / Integer(2 * k) * f[(0,0)])
                                           if is_integral else
                                           -bernoulli(k) / Integer(2 * k) * f[(0,0)] ) )
            else :
                coeffs[((a,b,c),l)] = ( sigma(c//self.level(), k-1) * f[(0,0)]
                                        if c != 0 else 
                                         ( Integer(-bernoulli(k) / Integer(2 * k) * f[(0,0)])
                                           if is_integral else
                                           -bernoulli(k) / Integer(2 * k) * f[(0,0)] ) )
        
        return coeffs
