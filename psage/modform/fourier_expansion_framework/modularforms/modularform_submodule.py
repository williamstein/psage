r"""
Submodules of rings of orthogonal modular forms.

AUTHOR :
    -- Martin Raum (2009 - 07 - 30) Initial version
"""

#===============================================================================
# 
# Copyright (C) 2009 Martin Raum
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
from builtins import object
from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ambient import GradedExpansionAmbient_abstract
from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_submodule import GradedExpansionSubmodule_ambient_pid, \
                                                       GradedExpansionSubmodule_submodule_pid, \
                                                       GradedExpansionSubmodule_abstract, \
                                                       GradedExpansionSubmoduleVector_generic
from sage.misc.cachefunc import cached_method
from sage.structure.sequence import Sequence
import operator

#===============================================================================
# HeckeInvariantSubmodule_abstract
#===============================================================================

class HeckeInvariantSubmodule_abstract(object) :       
    @cached_method
    def _hecke_action(self, n) :
        """
        Calculate the action of `T(n)` on the basis.
        
        INPUT:
            - `n`      -- A Hecke modulus. Probably an integer.
        
        OUTPUT:
            A sequence of elements of self.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_submodule import *
            sage: ma = ModularFormsModule_withheckeaction( QQ, ModularFormTestType_vectorvalued(), NNFilter(5) )
            sage: sm = ma.graded_submodule(3)
            sage: sm._hecke_action(7)
            [(1)]
        """
        T = self.graded_ambient().type()._hecke_operator_class()(n)
        
        return Sequence( [ self( T.eval(b.fourier_expansion(), self.graded_ambient()(b).weight()) )
                           for b in self.basis() ],
                         universe = self )
        
    @cached_method
    def hecke_homomorphism(self, n) :
        """
        Calculate the matrix corresponding to the action of `T(n)` on the basis.

        INPUT:
            - `n`      -- A Hecke modulus. Probably an integer.

        OUTPUT:
            A homomorphism from ``self`` to ``self``.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_submodule import *
            sage: ma = ModularFormsModule_withheckeaction( QQ, ModularFormTestType_vectorvalued(), NNFilter(5) )
            sage: sm = ma.graded_submodule(3)
            sage: sm.hecke_homomorphism(7)
            Free module morphism defined by the matrix
            [1]
            Domain: Submodule of Graded expansion module with generators v1, v2, ...
            Codomain: Submodule of Graded expansion module with generators v1, v2, ...
        """
        return self.hom(self._hecke_action(n))
    
    @cached_method
    def hecke_eigenforms(self, n ) :
        """
        Return a basis of eigenforms with respect to the Hecke operator `T(n)`.

        INPUT:
            - `n`      -- A Hecke modulus. Probably an integer.

        OUTPUT:
            A list of elements in the graded ambient.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_submodule import *
            sage: ma = ModularFormsModule_withheckeaction( QQ, ModularFormTestType_vectorvalued(), NNFilter(5) )
            sage: sm = ma.graded_submodule(3)
            sage: sm.hecke_eigenforms(7)
            [Graded expansion v1]
        """
        hm = self.hecke_homomorphism(n).matrix().transpose()
        hes = hm.eigenspaces_right()
        
        efs = []
        for e in hes :
            ring_elements = Sequence( [sum(map(operator.mul, self._basis_in_graded_ambient(), b.list()))
                                       for b in e[1].basis() ] )

            efs += list(ring_elements)
                 
        return efs 

#===============================================================================
# ModularFormsSubmoduleHeckeInvariant
#===============================================================================

def ModularFormsSubmoduleHeckeInvariant(arg1, arg2) :
    """
    INPUT:
        - ``arg1``     -- A graded ambient or an ambient module.
        - ``arg2``     -- A list of elements in ``arg1``. The basis.
    
    OUTPUT:
        An submodule of modular forms.
    
    TESTS::
        sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import *
        sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
        sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_submodule import *
        sage: ma = ModularFormsModule_withheckeaction( QQ, ModularFormTestType_vectorvalued(), NNFilter(5) )
        sage: sm = ma.graded_submodule(6)
        sage: ssm = ModularFormsSubmoduleHeckeInvariant(sm, [sm.0])
    """
    if isinstance(arg1, GradedExpansionAmbient_abstract) :
        return ModularFormsSubmodule_heckeinvariant_ambient(arg1, arg2)
    elif isinstance(arg1, GradedExpansionSubmodule_abstract) :
        return ModularFormsSubmodule_heckeinvariant_submodule(arg1, arg2)
    else :
        raise ValueError( "Cannot construct subspace in %s with basis %s." % (arg1, arg2) )

    return ModularFormsSubmodule_heckeinvariant_ambient

#===============================================================================
# ModularFormsSubmodule_heckeinvariant_ambient
#===============================================================================

class ModularFormsSubmodule_heckeinvariant_ambient ( 
        GradedExpansionSubmodule_ambient_pid,
        HeckeInvariantSubmodule_abstract ) :
    
    def __init__(self, graded_ambient, basis, **kwds) :
        """
        INPUT:
            - ``graded_ambient`` -- A graded ambient.
            - ``basis``          -- A list of elements in the graded ambient.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_submodule import *
            sage: ma = ModularFormsModule_withheckeaction( QQ, ModularFormTestType_vectorvalued(), NNFilter(5) )
            sage: sm = ModularFormsSubmodule_heckeinvariant_ambient(ma, [ma.0])
        """
        if not hasattr(self, '_element_class') :
            try :
                _element_class = graded_ambient.type()._space_element_class()
            except NotImplementedError :
                _element_class = GradedExpansionSubmoduleVector_generic

        GradedExpansionSubmodule_ambient_pid.__init__(self, graded_ambient, basis,
                                                      **kwds)

#===============================================================================
# ModularFormsSubmodule_heckeinvariant_submodule
#===============================================================================

class ModularFormsSubmodule_heckeinvariant_submodule (
        GradedExpansionSubmodule_submodule_pid,
        HeckeInvariantSubmodule_abstract ) :
    
    def __init__(self, ambient, basis, **kwds) :
        """
        INPUT:
            - ``ambient`` -- A submodule of modular forms.
            - ``basis``   -- A list of elements in the ambient.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_submodule import *
            sage: ma = ModularFormsModule_withheckeaction( QQ, ModularFormTestType_vectorvalued(), NNFilter(5) )
            sage: sm = ModularFormsSubmodule_heckeinvariant_ambient(ma, [ma.0, ma.1])
            sage: ssm = ModularFormsSubmodule_heckeinvariant_submodule(sm, [sm.1])
        """
        if not hasattr(self, '_element_class') :
            try :
                _element_class = ambient.graded_ambient().type()._space_element_class()
            except NotImplementedError :
                _element_class = GradedExpansionSubmoduleVector_generic
        
        self.__integral_basis = None

        GradedExpansionSubmodule_submodule_pid.__init__(self, ambient, basis,
                                                        **kwds)
        
    def _set_integral_basis(self, basis) :
        """
        An additional basis, that will not be echelonized.
        
        INPUT:
            - ``basis`` -- A list of elements of ``self``. A basis for ``self``.
        
        OUTPUT:
            ``None``.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_submodule import *
            sage: ma = ModularFormsModule_withheckeaction( QQ, ModularFormTestType_vectorvalued(), NNFilter(5) )
            sage: sm = ModularFormsSubmodule_heckeinvariant_ambient(ma, [ma.0, ma.1])
            sage: ssm = ModularFormsSubmodule_heckeinvariant_submodule(sm, [sm.1])
            sage: ssm._set_integral_basis([2 * sm.1])
            sage: ssm._integral_basis()
            [(0, 2)]
        """
        self.__integral_basis = basis
        
    def _integral_basis(self) :
        """
        An additional basis, that will not be echelonized.

        OUTPUT:
            A list of elements of ``self``.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_submodule import *
            sage: ma = ModularFormsModule_withheckeaction( QQ, ModularFormTestType_vectorvalued(), NNFilter(5) )
            sage: sm = ModularFormsSubmodule_heckeinvariant_ambient(ma, [ma.0, ma.1])
            sage: ssm = ModularFormsSubmodule_heckeinvariant_submodule(sm, [sm.1])
            sage: ssm._integral_basis()
            Traceback (most recent call last):
            ...
            RuntimeError: Integral basis is not set.
        """
        if self.__integral_basis is None :
            raise RuntimeError( "Integral basis is not set." )
        
        return self.__integral_basis

#===============================================================================
# ModularFormsWeightSubmodule
#===============================================================================

def ModularFormsWeightSubmodule(graded_ambient, basis, weight) :
    """
    A submodule of modular forms, that have the same weights.
    
    INPUT:
        - ``graded_ambient`` -- A graded ambient.
        - ``basis``          -- A list of elements in the graded ambient.
        - ``weight``         -- A grading index. The commen weight of all basis elements.
    
    OUTPUT:
        An instance of :class:~`.ModularFormsSubmodule_singleweight_ambient_pid`.
    
    TESTS::
        sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import *
        sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
        sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_submodule import *
        sage: ma = ModularFormsModule_withheckeaction( QQ, ModularFormTestType_vectorvalued(), NNFilter(5) )
        sage: sm = ModularFormsWeightSubmodule(ma, [ma.0], 3)
    """
    return ModularFormsSubmodule_singleweight_ambient_pid(graded_ambient, basis, weight)

#===============================================================================
# ModularFormSubmodule_singleweight
#===============================================================================

class ModularFormsSubmodule_singleweight_ambient_pid ( 
        GradedExpansionSubmodule_ambient_pid,
        HeckeInvariantSubmodule_abstract ) :
    
    def __init__(self, graded_ambient, basis, weight, **kwds) :
        """
        INPUT:
            - ``graded_ambient``   -- A graded ambient.
            - ``basis``            -- A list of elements in the graded ambient.
            - ``weight``           -- A grading index. The commen weight of all basis elements.
            - ``kwds``             -- Will be forwarded to :class:~`fourier_expansion_framework.gradedexpansions.gradedexpansion_submodule.GradedExpansionSubmodule_ambient_pid`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_submodule import *
            sage: ma = ModularFormsModule_withheckeaction( QQ, ModularFormTestType_vectorvalued(), NNFilter(5) )
            sage: sm = ModularFormsSubmodule_singleweight_ambient_pid(ma, [ma.0], 3)
        """
        self.__weight = weight
        if not hasattr(self, '_element_class') :
            try :
                _element_class = graded_ambient.type()._space_element_class()
            except NotImplementedError :
                _element_class = GradedExpansionSubmoduleVector_generic
             
        GradedExpansionSubmodule_ambient_pid.__init__(self, graded_ambient, basis,
                                                      **kwds)

    def weight(self) :
        return self.__weight

## TODO: Provide standard implementation of is_cusp_form for vectors
## TODO: Provide implementation of homogeneous weight vectors and
##       make Hecke operators work for this
