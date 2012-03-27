r"""
Interfaces for modular forms which admit Hecke actions or ring which have
Maass lifts.

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

from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_element import GradedExpansion_abstract
from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import MonoidPowerSeries_abstract, \
                                                        EquivariantMonoidPowerSeries_abstract
from psage.modform.fourier_expansion_framework.modularforms.modularform_submodule import ModularFormsSubmoduleHeckeInvariant
from psage.modform.fourier_expansion_framework.gradedexpansions.fourierexpansionwrapper import FourierExpansionWrapper

#===============================================================================
# ModularFormsAmbientWithHeckeAction_abstract
#===============================================================================

class ModularFormsAmbientWithHeckeAction_abstract :
    """
    The standard implementation assumes that the action only depends on the
    modulus and the weight.
    The deriving class must override self._hecke_operator_class or it will
    be derived from the type. 
    """
    
    def __init__(self, type) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: ma = ModularFormsModule_withheckeaction( QQ, ModularFormTestType_vectorvalued(), NNFilter(5) )
            sage: sm = ma.graded_submodule(3)
            sage: sm._hecke_action(7)
            [(1)]
            sage: ma = ModularFormsRing_withheckeaction( QQ, ModularFormTestType_scalar(), NNFilter(5) )
            Traceback (most recent call last):
            ...
            ValueError: Type for modular forms ambients with Hecke action must support Hecke operators.
        """
        if not hasattr(self, '_hecke_operator_class') :
            if not type.has_hecke_action() :
                raise ValueError( "Type for modular forms ambients with Hecke action must support Hecke operators." )
            self._hecke_operator_class = type._hecke_operator_class()
        
        hecke_invariant_pred = lambda basis, **kwds : "is_hecke_invariant" in kwds and kwds["is_hecke_invariant"]
        def hecke_invariant_fcn(basis, **kwds) :
            try :
                return self.type()._submodule_heckeinvariant_class(self, basis, **kwds)
            except NotImplementedError :
                return ModularFormsSubmoduleHeckeInvariant(self, basis, **kwds)
        
        self._submodule_classes.insert(-2, (hecke_invariant_pred, hecke_invariant_fcn))

    def _hecke_action(self, n, form) :
        """
        The image of ``form`` under the `n`-th Hecke operator.
        
        INPUT:
            - `n`      -- A Hecke modulus. Probably an integer.
            - ``form`` -- An element whose parent supports conversion from
                          Fourier expansions.
        
        OUTPUT:
            An element in the parent of ``form``.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: ma = ModularFormsModule_withheckeaction( QQ, ModularFormTestType_vectorvalued(), NNFilter(5) )
            sage: sm = ma.graded_submodule(3)
            sage: ma._hecke_action(2, sm.0)
            Equivariant monoid power series in Module of equivariant monoid power series over NN
            sage: ma._hecke_action(7, sm.0)
            (1)
            sage: ma._hecke_action(2, 2)
            Traceback (most recent call last):
            ...
            TypeError: Form must be a Fourier expansion or wrap one.
        """
        T = self._hecke_operator_class(n)
        
        if isinstance(form, FourierExpansionWrapper) :
            expansion = form.fourier_expansion()
            try :
                weight = form.weight()
            except AttributeError :
                weight = None 
        elif isinstance(form, (MonoidPowerSeries_abstract, EquivariantMonoidPowerSeries_abstract)) :
            expansion = form
            weight = None
        else :
            raise TypeError( "Form must be a Fourier expansion or wrap one." ) 
        
        hecke_expansion = T.eval(expansion, weight)
        
        if isinstance(form, FourierExpansionWrapper) :
            try :
                return form.parent()(hecke_expansion)
            except (ValueError, ArithmeticError) :
                pass
        
        return hecke_expansion
