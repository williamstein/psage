"""
Types of Jacobi forms of fixed index and weight.

AUTHOR :
    - Martin Raum (2010 - 04 - 07) Initial version.
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
from operator import xor
from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import TrivialGrading
from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import ModularFormsModule_generic
from psage.modform.fourier_expansion_framework.modularforms.modularform_types import ModularFormType_abstract
from psage.modform.jacobiforms.jacobiformd1nn_fegenerators import jacobi_form_by_taylor_expansion,\
    _jacobi_forms_by_taylor_expansion_coords
from psage.modform.jacobiforms.jacobiformd1nn_fourierexpansion import JacobiD1NNFourierExpansionModule, \
                                                                      JacobiFormD1NNFilter
from sage.categories.number_fields import NumberFields
from sage.matrix.constructor import diagonal_matrix, matrix, zero_matrix,\
    identity_matrix
from sage.misc.cachefunc import cached_method
from sage.misc.mrange import mrange
from sage.modular.modform.constructor import ModularForms
from sage.rings.all import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.number_field.number_field import CyclotomicField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.structure.sequence import Sequence
from functools import reduce


#===============================================================================
# JacobiFormsD1NN
#===============================================================================

_jacobiforms_cache = dict()

def JacobiFormsD1NN(A, type, precision, *args, **kwds) :
    global _jacobiforms_cache
    
    if isinstance(precision, (int, Integer)) :
        precision = JacobiFormD1NNFilter(precision, type.index())
    
    k = (A, type, precision)
    
    try :
        return _jacobiforms_cache[k]
    except KeyError :
        if isinstance(type, JacobiFormD1NN_Gamma) :
            M = ModularFormsModule_generic(A, type, precision)
        else :
            raise TypeError("{0} must be a Jacobi form type".format(type))
        
        _jacobiforms_cache[k] = M
        return M
    
#===============================================================================
# JacobiFormD1NN_Gamma
#===============================================================================


class JacobiFormD1NN_Gamma(ModularFormType_abstract):
    r"""
    Type of Jacobi forms of degree `1` and index in `\mathbb{N}` associated with
    the full modular group.
    
    TESTS::
    
        sage: from psage.modform.jacobiforms import *
        sage: from psage.modform.jacobiforms.jacobiformd1nn_fourierexpansion import JacobiFormD1NNFilter 
        sage: JR = JacobiFormsD1NN(QQ, JacobiFormD1NN_Gamma(3, 6), JacobiFormD1NNFilter(10, 3))
        sage: JR.gens()
        (Graded expansion TDE_0, Graded expansion TDE_1)
        sage: JR.0 + 2 * JR.1
        Graded expansion TDE_0 + 2*TDE_1
    """
    def __init__(self, index, weight) :
        if weight % 2 != 0 :
            raise NotImplementedError("Only even weight forms are implemented.")
  
        self.__index = index
        self.__weight = weight        
    
    def index(self) :
        return self.__index
    
    def weight(self) :
        return self.__weight
    
    def _ambient_construction_function(self) :
        return JacobiFormsD1NN
    
    def group(self) :
        return "Sp(2, ZZ)_\infty"
    
    @cached_method
    def _rank(self, K) :
        
        if K is QQ or K in NumberFields():
            return len(_jacobi_forms_by_taylor_expansion_coords(self.__index, self.__weight, 0))

            ## This is the formula used by Poor and Yuen in Paramodular cusp forms
            if self.__weight == 2 :
                delta = len(self.__index.divisors()) // 2 - 1
            else :
                delta = 0

            return sum(ModularForms(1, self.__weight + 2 * j).dimension() + j**2 // (4 * self.__index)
                        for j in range(self.__index + 1)) \
                   + delta


            ## This is the formula given by Skoruppa in
            ## Jacobi forms of critical weight and Weil representations
            ##FIXME: There is some mistake here
            if self.__weight % 2 != 0 :
                ## Otherwise the space X(i**(n - 2 k)) is different
                ## See: Skoruppa, Jacobi forms of critical weight and Weil representations
                raise NotImplementedError

            m = self.__index
            K = CyclotomicField(24 * m, 'zeta')
            zeta = K.gen(0)

            quadform = lambda x : 6 * x**2
            bilinform = lambda x,y : quadform(x + y) - quadform(x) - quadform(y)

            T = diagonal_matrix([zeta**quadform(i) for i in range(2*m)])
            S =   sum(zeta**(-quadform(x)) for x in range(2 * m)) / (2 * m) \
                * matrix([[zeta**(-bilinform(j,i)) for j in range(2*m)] for i in range(2*m)])
            subspace_matrix_1 = matrix([[1 if j == i or j == 2*m - i else 0 for j in range(m + 1)]
                                        for i in range(2*m)])
            subspace_matrix_2 = zero_matrix(ZZ, m + 1, 2*m)
            subspace_matrix_2.set_block(0,0,identity_matrix(m+1))

            T = subspace_matrix_2 * T * subspace_matrix_1
            S = subspace_matrix_2 * S * subspace_matrix_1

            sqrt3 = (zeta**(4*m) - zeta**(-4*m)) * zeta**(-6*m)
            rank =   (self.__weight - QQ(1)/QQ(2) - 1) / QQ(2) * (m + 1) \
                   + 1/8 * (zeta**(3*m * (2*self.__weight - 1)) * S.trace()
                             + zeta**(3*m * (1 - 2*self.__weight)) * S.trace().conjugate()) \
                   + 2/(3*sqrt3) * (zeta**(4 * m * self.__weight) * (S*T).trace()
                                     + zeta**(-4 * m * self.__weight) * (S*T).trace().conjugate()) \
                   - sum((j**2 % (m+1))/QQ(m+1) -QQ(1)/QQ(2) for j in range(0,m+1))

            if self.__weight > 5 / 2:
                return rank
            else:
                raise NotImplementedError

        raise NotImplementedError
    
    @cached_method
    def generators(self, K, precision):
        if K is QQ or K in NumberFields():
            return Sequence([jacobi_form_by_taylor_expansion(i, self.__index, self.__weight, precision)
                               for i in range(self._rank(K))],
                             universe = JacobiD1NNFourierExpansionModule(QQ, self.__index))
        
        raise NotImplementedError
    
    def grading(self, K):
        if K is QQ or K in NumberFields():
            return TrivialGrading(self._rank(K),
                                  (self.__index, self.__weight))
        
        raise NotImplementedError

    def _generator_names(self, K):
        if K is QQ or K in NumberFields():
            return ["TDE_{0}".format(i) for i in range(self._rank(K))]
        
        raise NotImplementedError

    def _generator_by_name(self, K, name):
        if K is QQ or K in NumberFields():
            R = self.generator_relations(K).ring()
            try:
                return R.gen(self._generator_names(K).index(name))
            except ValueError:
                raise ValueError("name {0} doesn't exist for {1}".format(name, K))
        
        raise NotImplementedError
    
    @cached_method
    def generator_relations(self, K):
        r"""
        An ideal I in a polynomial ring R, such that the associated module
        is (R / I)_1. 
        """
        if K is QQ or K in NumberFields():
            R = PolynomialRing(K, self._generator_names(K))
            return R.ideal(0)
            
        raise NotImplementedError

    def weights(self, K):
        r"""
        A list of integers corresponding to the weights.
        """
        if K is QQ or K in NumberFields():
            return len(self._theta_decomposition_indices()) \
                    * [(self.__index, self.__weight)]
            
        raise NotImplementedError
    
    def graded_submodules_are_free(self):
        return True

    def __cmp__(self, other):
        c = cmp(type(self), type(other))
        
        if c == 0:
            c = cmp(self.__index, other.__index)
        if c == 0:
            c = cmp(self.__weight, other.__weight)
            
        return c

    def __hash__(self):
        return reduce(xor, list(map(hash, [type(self), self.__index, self.__weight])))
