r"""
Types for paramodular forms of fixed level and weight.

AUTHOR :
    -- Martin Raum (2010 - 04 - 15) Initial version.
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

from past.builtins import cmp
from builtins import map
from builtins import range
from psage.modform.fourier_expansion_framework.gradedexpansions.expansion_module import ExpansionModule
from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import TrivialGrading
from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import ModularFormsModule_generic
from psage.modform.fourier_expansion_framework.modularforms.modularform_types import ModularFormType_abstract
from psage.modform.jacobiforms.jacobiformd1nn_types import JacobiFormD1NN_Gamma
from operator import xor
from psage.modform.paramodularforms.paramodularformd2_element import ParamodularFormD2_classical
from psage.modform.paramodularforms.paramodularformd2_fegenerators import gritsenko_products
from psage.modform.paramodularforms.paramodularformd2_fourierexpansion import ParamodularFormD2Filter_discriminant,\
    ParamodularFormD2FourierExpansionRing
from psage.modform.paramodularforms.paramodularformd2_heckeaction import ParamodularFormD2FourierExpansionHeckeAction
from psage.modform.paramodularforms.paramodularformd2_fegenerators import symmetrised_siegel_modular_forms
from sage.categories.number_fields import NumberFields
from sage.misc.cachefunc import cached_method
from sage.arith.all import kronecker_symbol
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.rational_field import QQ
from sage.structure.sequence import Sequence
from sage.modular.all import ModularForms
from functools import reduce

#===============================================================================
# ParamodularFormsD2
#===============================================================================

_paramodularforms_cache = dict()

def ParamodularFormsD2(A, type, precision, *args, **kwds) :
    global _paramodularforms_cache
    if isinstance(precision, (int, Integer)) :
        precision = ParamodularFormD2Filter_discriminant(precision, type.level()) 
    k = (A, type, precision)
    
    try :
        return _paramodularforms_cache[k]
    except KeyError :
        if isinstance(type,ParamodularFormD2_Gamma) :
            M = ModularFormsModule_generic(A, type, precision)
        else :
            raise TypeError("{0} must be a paramodular form type".format(type))
        
        _paramodularforms_cache[k] = M
        return M

#===============================================================================
# ParamodularFormD2_Gamma
#===============================================================================

class ParamodularFormD2_Gamma ( ModularFormType_abstract ) :
    """
    Type of paramodular forms of degree 2. 
    """
    def __init__(self, level, weight) :
        if weight % 2 != 0 :
            raise NotImplementedError("Only even weight forms are implemented.")
  
        if level != 1 and not Integer(level).is_prime() :
            raise NotImplementedError("Only prime level or level 1 is implemented.")
  
        self.__level = Integer(level)
        self.__weight = weight        
    
    def level(self) :
        return self.__level
    
    def _ambient_construction_function(self) :
        return ParamodularFormsD2
    
    def _ambient_element_class(self) :
        return ParamodularFormD2_classical
    
    def group(self) :
        return "\Gamma_2 ^para"
    
    def _rank(self, K) :
        """
        Return the dimension of the level N space of given weight. 
        """
        if not K is QQ :
            raise NotImplementedError
        
        if not self.__level.is_prime() :
            raise NotImplementedError
    
        if self.__weight % 2 != 0 :
            raise NotImplementedError( "Only even weights available")
    
        N = self.__level
        
        if N == 1 :
            ## By Igusa's theorem on the generators of the graded ring
            P = PowerSeriesRing(ZZ, 't')
            t = P.gen(0)
                                
            dims = 1 / ((1 - t**4) * (1 - t**6) * (1 - t**10) * (1 - t**12)).add_bigoh(self.__weight + 1)
            
            return dims[self.__weight]
        if N == 2 :
            ## As in Ibukiyama, Onodera - On the graded ring of modular forms of the Siegel
            ## paramodular group of level 2, Proposition 2
            P = PowerSeriesRing(ZZ, 't')
            t = P.gen(0)
                                
            dims = ((1 + t**10) * (1 + t ** 12) * (1 + t**11))
            dims = dims / ((1 - t**4) * (1 - t**6) * (1 - t**8) * (1 - t**12)).add_bigoh(self.__weight + 1)
            
            return dims[self.__weight]
        if N == 3 :
            ## By Dern, Paramodular forms of degree 2 and level 3, Corollary 5.6
            P = PowerSeriesRing(ZZ, 't')
            t = P.gen(0)
            
            dims = ((1 + t**12) * (1 + t**8 + t**9 + t**10 + t**11 + t**19))
            dims = dims / ((1 - t**4) * (1 - t**6)**2 * (1 - t**12)).add_bigoh(self.__weight + 1)
            
            return dims[self.__weight]
        if N == 5 :
            ## By Marschner, Paramodular forms of degree 2 with particular emphasis on level t = 5,
            ## Corollary 7.3.4. PhD thesis electronically available via the library of
            ## RWTH University, Aachen, Germany  
            P = PowerSeriesRing(ZZ, 't')
            t = P.gen(0)
            
            dims =   t**30 + t**24 + t**23 + 2*t**22 + t**21 + 2*t**20 + t**19 + 2*t**18 \
                   + 2*t**16 + 2*t**14 + 2*t**12 + t**11 + 2*t**10 + t**9 + 2*t**8 + t**7 + t**6 + 1
            dims = dims / ((1 - t**4) * (1 - t**5) * (1 - t**6) * (1 - t**12)).add_bigoh(self.__weight + 1)
            
            return dims[self.__weight]
        
        if self.__weight == 2 :
            ## There are only cuspforms, since there is no elliptic modular form
            ## of weight 2.
            if N < 277 :
                ## Poor, Yuen - Paramodular cusp forms tells us that all forms are
                ## Gritsenko lifts
                return JacobiFormD1NN_Gamma(self.__level, 2)._rank(QQ)

            raise NotImplementedError
        elif self.__weight == 4 :
            ## This is the formula cited by Poor and Yuen in Paramodular cusp forms
            cuspidal_dim =  Integer(   (N**2 - 143) / Integer(576) + N / Integer(8)
                                     + kronecker_symbol(-1, N) * (N - 12) / Integer(96)
                                     + kronecker_symbol(2, N) / Integer(8)
                                     + kronecker_symbol(3, N) / Integer(12)
                                     + kronecker_symbol(-3, N) * N / Integer(36) )
        else :
            ## This is the formula given by Ibukiyama in
            ## Relations of dimension of automorphic forms of Sp(2,R) and its compact twist Sp(2),
            ## Theorem 4
            p = N
            k = self.__weight
            
            ## This is the reversed Ibukiyama symbol [.., .., ..; ..]
            def ibukiyama_symbol(modulus, *args) :
                return args[k % modulus]

            ## if p == 2 this formula is wrong. If the weight is even it differs by
            ## -3/16 from the true dimension and if the weight is odd it differs by
            ## -1/16 from the true dimension. 
            H1 = (p**2 + 1) * (2 * k - 2) * (2 * k - 3) * (2 * k - 4) / Integer(2**9 * 3**3 * 5)
            H2 = (-1)**k * (2 * k - 2) * (2 * k - 4) / Integer(2**8 * 3**2) \
                 + ( (-1)**k * (2 * k - 2) * (2 * k - 4) / Integer(2**7 * 3)
                     if p !=2 else
                     (-1)**k * (2 * k - 2) * (2 * k - 4) / Integer(2**9) )
            H3 = ( ibukiyama_symbol(4, k - 2, -k + 1, -k + 2, k - 1) / Integer(2**4 * 3)
                   if p != 3 else
                   5 * ibukiyama_symbol(4, k - 2, -k + 1, -k + 2, k - 1) / Integer(2**5 * 3) )
            H4 = ( ibukiyama_symbol(3, 2 * k - 3, -k + 1, -k + 2) / Integer(2**2 * 3**3)
                   if p != 2 else
                   5 * ibukiyama_symbol(3, 2 * k - 3, -k + 1, -k + 2) / Integer(2**2 * 3**3) )
            H5 = ibukiyama_symbol(6, -1, -k + 1, -k + 2, 1, k - 1, k - 2) / Integer(2**2 * 3**2)
            if p % 4 == 1 :
                H6 = 5 * (2 * k - 3) * (p + 1) / Integer(2**7 * 3) + (-1)**k * (p + 1) / Integer(2**7)
            elif p % 4 == 3 :
                H6 = (2 * k - 3) * (p - 1) / Integer(2**7) + 5 * (-1)**k * (p - 1) / Integer(2**7 * 3)
            else :
                H6 = 3 * (2 * k - 3) / Integer(2**7) + 7 * (-1)**k / Integer(2**7 * 3)
            if p % 3 == 1 :
                H7 =   (2 * k - 3) * (p + 1) / Integer(2 * 3**3) \
                     + (p + 1) * ibukiyama_symbol(3, 0, -1, 1) / Integer(2**2 * 3**3)
            elif p % 3 == 2 :
                H7 =   (2 * k - 3) * (p - 1) / Integer(2**2 * 3**3) \
                     + (p - 1) * ibukiyama_symbol(3, 0, -1, 1) / Integer(2 * 3**3)
            else :
                H7 =   5 * (2 * k - 3) / Integer(2**2 * 3**3) \
                     + ibukiyama_symbol(3, 0, -1, 1) / Integer(3**3)
            H8 = ibukiyama_symbol(12, 1, 0, 0, -1, -1, -1, -1, 0, 0, 1, 1, 1) / Integer(2 * 3)
            H9 = ( 2 * ibukiyama_symbol(6, 1, 0, 0, -1, 0, 0) / Integer(3**2)
                   if p != 2 else
                   ibukiyama_symbol(6, 1, 0, 0, -1, 0, 0) / Integer(2 * 3**2) )
            H10 = (1 + kronecker_symbol(5, p)) * ibukiyama_symbol(5, 1, 0, 0, -1, 0) / Integer(5)
            H11 = (1 + kronecker_symbol(2, p)) * ibukiyama_symbol(4, 1, 0, 0, -1) / Integer(2**3)
            if p % 12 == 1 :
                H12 = ibukiyama_symbol(3, 0, 1, -1) / Integer(2 * 3)
            elif p % 12 == 11 :
                H12 = (-1)**k / Integer(2 * 3)
            elif p == 2 or p == 3 :
                H12 = (-1)**k / Integer(2**2 * 3)
            else :
                H12 = 0
                
            I1 = ibukiyama_symbol(6, 0, 1, 1, 0, -1, -1) / Integer(6)
            I2 = ibukiyama_symbol(3, -2, 1, 1) / Integer(2 * 3**2)
            if p == 3 :
                I3 = ibukiyama_symbol(3, -2, 1, 1)/ Integer(3**2)
            elif p % 3 == 1 :
                I3 = 2 * ibukiyama_symbol(3, -1, 1, 0) / Integer(3**2)
            else :
                I3 = 2 * ibukiyama_symbol(3, -1, 0, 1) / Integer(3**2)
            I4 = ibukiyama_symbol(4, -1, 1, 1, -1) / Integer(2**2)
            I5 = (-1)**k / Integer(2**3)
            I6 = (-1)**k * (2 - kronecker_symbol(-1, p)) / Integer(2**4)
            I7 = -(-1)**k * (2 * k - 3) / Integer(2**3 * 3)
            I8 = -p * (2 * k - 3) / Integer(2**4 * 3**2)
            I9 = -1 / Integer(2**3 * 3)
            I10 = (p + 1) / Integer(2**3 * 3)
            I11 = -(1 + kronecker_symbol(-1, p)) / Integer(8)
            I12 = -(1 + kronecker_symbol(-3, p)) / Integer(6)
            
            cuspidal_dim =   H1 + H2 + H3 + H4 + H5 + H6 + H7 + H8 + H9 + H10 + H11 + H12 \
                           + I1 + I2 + I3 + I4 + I5 + I6 + I7 + I8 + I9 + I10 + I11 + I12
                       

        mfs = ModularForms(1, self.__weight)
        
        return cuspidal_dim + mfs.dimension() + mfs.cuspidal_subspace().dimension()

    def _gritsenko_products(self, precision = None) :
        try :
            if precision is None :
                return self.__gritsenko_products
            elif precision == self.__gritsenko_products_precision :
                return self.__gritsenko_products
            elif precision < self.__gritsenko_products_precision :
                return ( self.__gritsenko_products[0],
                         [p.truncate(precision) for p in self.__gritsenko_products[1]] )
        except AttributeError :
            pass
        
        self.__gritsenko_products = \
          gritsenko_products(self.__level, self.__weight, self._rank(QQ), precision)
        self.__gritsenko_products_precision = precision
    
        return self.__gritsenko_products
        
    def _symmetrised_siegel_modular_forms(self, precision = None) :
        try :
            if precision is None :
                return self.__symmetrised_siegel_modular_forms
            elif precision == self.__symmetriesed_siegel_modular_forms_precision :
                return self.__symmetrised_siegel_modular_forms
            elif precision < self.__symmetriesed_siegel_modular_forms_precision :
                return [s.truncate(precision) for s in self.__symmetrised_siegel_modular_forms]
        except AttributeError :
            pass
        
        self.__symmetrised_siegel_modular_forms = \
          symmetrised_siegel_modular_forms(self.__level, self.__weight, precision)
        self.__symmetrised_siegel_modular_forms_precision = precision
        
        return self.__symmetrised_siegel_modular_forms 
        
    @cached_method
    def generators(self, K, precision) :
        if K is QQ or K in NumberFields() :
            gps = self._gritsenko_products(precision)[1]
            if len(gps) != self._rank(QQ) :
                syms = self._symmetrised_siegel_modular_forms(precision)
                em = ExpansionModule(Sequence(gps + syms, universe = ParamodularFormD2FourierExpansionRing(QQ, self.__level) ) )
                gens = [e.fourier_expansion() for e in em.pivot_elements()]
            else :
                gens = gps
                            
            if len(gens) == self._rank(QQ) :
                return Sequence( gens, universe = ParamodularFormD2FourierExpansionRing(QQ, self.__level) )
            
            raise ArithmeticError("Gritsenko products do not span this space.")
            
        raise NotImplementedError
    
    def grading(self, K) :
        if K is QQ or K in NumberFields() :
            return TrivialGrading( self._rank(K), self.__weight )
        
        raise NotImplementedError

    def _generator_names(self, K) :
        if K is QQ or K in NumberFields() :
            ## We assume that the space is spanned by Gritsenko products
            ## Introduce new names, as soon as new cases are implemented
            nmb_gps = len(self._gritsenko_products(None)[0])
            return [ "GP_%s" % (i,) for i in range(nmb_gps)] + \
                   [ "SymS_%s" % (i,) for i in range(self._rank(K) - nmb_gps) ]

        raise NotImplementedError
    
    def _generator_by_name(self, K, name) :
        if K is QQ or K in NumberFields() :
            R = self.generator_relations(K).ring()
            try :
                return R.gen(self._generator_names(K).index(name))
            except ValueError :
                raise ValueError("name {0} doesn't exist for {1}".format(name, K))
        
        raise NotImplementedError
    
    @cached_method
    def generator_relations(self, K) :
        r"""
        An ideal I in a polynomial ring R, such that the associated module
        is (R / I)_1. 
        """
        if K is QQ or K in NumberFields() :
            R = PolynomialRing(K, self._generator_names(K))
            return R.ideal(0)
            
        raise NotImplementedError

    def reduce_before_evaluating(self, K) :
        return False

    def weights(self, K) :
        """
        A list of integers corresponding to the weights.
        """
        if K is QQ or K in NumberFields() :
            return self._rank(K) * [self.__weight]
            
        raise NotImplementedError
    
    def graded_submodules_are_free(self) :
        return True
    
    def _hecke_operator_class(self) :
        return ParamodularFormD2FourierExpansionHeckeAction

    def __cmp__(self, other) :
        c = cmp(type(self), type(other))
        
        if c == 0 :
            c = cmp(self.__level, other.__level)
        if c == 0 :
            c = cmp(self.__weight, other.__weight)
            
        return c

    def __hash__(self) :
        return reduce(xor, list(map(hash, [self.__level, self.__weight])))
