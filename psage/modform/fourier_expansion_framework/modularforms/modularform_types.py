r"""
An abstract class for modular form types.

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

from sage.structure.sage_object import SageObject

#===============================================================================
# ModularFormType_abstract
#===============================================================================

class ModularFormType_abstract ( SageObject ) :
    r"""
    Types should be globally unique.
    """
    
    def _ambient_construction_function(self) :
        """
        Return a function that will construct the ambient ring or module
        of modular forms.
        
        OUTPUT:
            A function with INPUT:
                - `A`      -- A ring or module; The Fourier coefficients' domain.
                - ``type`` -- A type of modular forms.
                - ``precision`` -- A precision class; The Fourier expansion's precision.
            and OUTPUT:
                A ring or module of modular forms.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_types import ModularFormType_abstract
            sage: ModularFormType_abstract()._ambient_construction_function()
            Traceback (most recent call last):
            ...
            NotImplementedError: Subclass has to implement this function.
        """
        raise NotImplementedError( "Subclass has to implement this function." )
    
    def _ambient_element_class(self) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_types import ModularFormType_abstract
            sage: ModularFormType_abstract()._ambient_element_class()
            Traceback (most recent call last):
            ...
            NotImplementedError: Subclass has to implement this function.
        """
        raise NotImplementedError( "Subclass has to implement this function." )
    
    def _space_element_class(self) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_types import ModularFormType_abstract
            sage: ModularFormType_abstract()._space_element_class()
            Traceback (most recent call last):
            ...
            NotImplementedError: Subclass has to implement this function.
        """
        raise NotImplementedError( "Subclass has to implement this function." )
    
    def _weight_submodule_class(self) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_types import ModularFormType_abstract
            sage: ModularFormType_abstract()._weight_submodule_class()
            Traceback (most recent call last):
            ...
            NotImplementedError: Subclass has to implement this function.
        """
        raise NotImplementedError( "Subclass has to implement this function." )
    
    def _submodule_heckeinvariant_class(self) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_types import ModularFormType_abstract
            sage: ModularFormType_abstract()._submodule_heckeinvariant_class()
            Traceback (most recent call last):
            ...
            NotImplementedError: Subclass has to implement this function.
        """
        raise NotImplementedError( "Subclass has to implement this function." )
    
    def group(self) :
        """
        The modular group which ``self`` is associated with.
        
        OUTPUT:
            An arbitrary type.
        
        NOTE:
            The framwork might change later such that this function has
            to return a group.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_types import ModularFormType_abstract
            sage: ModularFormType_abstract().group()
            Traceback (most recent call last):
            ...
            NotImplementedError: Subclass has to implement this function.
        """
        raise NotImplementedError( "Subclass has to implement this function." )
    
    def base_ring_generators(self, K, precision) :
        """
        If the ring of modular forms can be interpreted as an algebra
        over a ring of modular forms with much simpler Fourier coefficient
        domains, it is a good idea to implement this here.
        Return the Fourier expansions.
        The parent of these expansions is expected to admit coercion
        into the ring of the generators. 
        
        INPUT:
            - `K`           -- A ring.
            - ``precision`` -- A precision instance.
        
        OUTPUT:
            ``None`` or a sequence of equivariant monoid power series.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_types import ModularFormType_abstract
            sage: ModularFormType_abstract().base_ring_generators(QQ, None) is None
            True
        """
        return None
    
    def generators(self, K, precision) :
        """
        A list of Fourier expansions of forms that generate the ring
        or module of modular forms.
        
        INPUT:
            - `K`           -- A ring or module; The ring of Fourier coefficients.
            - ``precision`` -- A precision class; The precision of the Fourier
                               expansions.
        
        OUTPUT:
            A sequence of equivariant monoid power series.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_types import ModularFormType_abstract
            sage: ModularFormType_abstract().generators(QQ, None)
            Traceback (most recent call last):
            ...
            NotImplementedError: Subclass has to implement this function.
        """
        raise NotImplementedError( "Subclass has to implement this function." )

    def grading(self, K) :
        """
        A grading for the ring or module of modular forms.
        
        INPUT:
            - `K` -- A ring or module; The domain of Fourier coefficients.
        
        OUTPUT:
            A grading class.
        
        NOTE:
            This will usually coincide with the weight grading.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_types import ModularFormType_abstract
            sage: ModularFormType_abstract().grading(QQ)
            Traceback (most recent call last):
            ...
            NotImplementedError: Subclass has to implement this function.
        """
        raise NotImplementedError( "Subclass has to implement this function." )
    
    def _generator_names(self, K) :
        """
        Names of the generators returned by :meth:~`.generators` within the
        attached polynomial ring.
        
        INPUT:
            - `K` -- A ring.
        
        OUTPUT:
            A list of strings.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_types import ModularFormType_abstract
            sage: ModularFormType_abstract()._generator_names(QQ)
            Traceback (most recent call last):
            ...
            NotImplementedError: Subclass has to implement this function.
        """
        raise NotImplementedError( "Subclass has to implement this function." )
    
    def _generator_by_name(self, K, name) :
        """
        Return the generator ``name`` as an element of the attached
        polynomial ring.
        
        INPUT:
            - `K`      -- A ring or module; The ring of Fourier coefficients.
            - ``name`` -- A string; The generator's name.
        
        OUTPUT:
            An element of a polynomial ring.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_types import ModularFormType_abstract
            sage: ModularFormType_abstract()._generator_names(QQ)
            Traceback (most recent call last):
            ...
            NotImplementedError: Subclass has to implement this function.
        """
        raise NotImplementedError( "Subclass has to implement this function." )
    
    def generator_relations(self, K) :
        """
        An ideal `I` in the attach polynomial ring `R`, such that the ring or module of
        modular forms is a subquotient of `R / I`. This ideal must be unique for `K`.
        
        INPUT:
            - `K`      -- A ring or module; The ring of Fourier coefficients.
        
        OUTPUT:
            An ideal in a polynomial ring.
            
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_types import ModularFormType_abstract
            sage: ModularFormType_abstract().generator_relations(QQ)
            Traceback (most recent call last):
            ...
            NotImplementedError: Subclass has to implement this function.
        """
        raise NotImplementedError( "Subclass has to implement this function." )
    
    def reduce_before_evaluating(self, K) :
        """
        Determine whether polynomials in the generators should first be
        Groebner reduced, before they are evaluated.
        
        INPUT:
            - `K` -- A ring.
        
        OUTPUT:
            A boolean.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_types import ModularFormType_abstract
            sage: ModularFormType_abstract().reduce_before_evaluating(QQ)
            True
        """
        return True
    
    def weights(self, K) :
        """
        The weights of the generators returned by :meth:~`.generators`.
        
        INPUT:
            - `K`      -- A ring or module; The ring of Fourier coefficients.
        
        OUTPUT:
            A list of integers.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_types import ModularFormType_abstract
            sage: ModularFormType_abstract().weights(QQ)
            Traceback (most recent call last):
            ...
            NotImplementedError: Subclass has to implement this function.
        """
        raise NotImplementedError( "Subclass has to implement this function." )

    def is_vector_valued(self) :
        """
        ``True`` if this is the vector valued version of a scalar valued type of modular forms.
        
        OUTPUT:
            A boolean.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_types import ModularFormType_abstract
            sage: ModularFormType_abstract().is_vector_valued()
            False
        """
        return False
    
    def _hom_base_extension(self, K, L) :
        """
        Images of generators over `K` in terms of those over `L`.
        
        INPUT:
            - `K` -- A ring.
            - `L- -- A ring.
        
        OUTPUT:
            A dictionary with keys the generators over `K` and values in the set of generators
            over `L`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import ModularFormTestType_scalar
            sage: ModularFormTestType_scalar()._hom_base_extension(ZZ, QQ)
            {g4: g4, g5: g5, g2: g2, g3: g3, g1: g1}
        """
        return dict( (self._generator_by_name(K, g), self._generator_by_name(L, g))
                     for g in self._generator_names(K) )
        
    def non_vector_valued(self) :
        """
        Return the non vector values version of this type. 
        
        OUTPUT:
            A type of modular forms.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_types import ModularFormType_abstract
            sage: ModularFormType_abstract().weights(QQ)
            Traceback (most recent call last):
            ...
            NotImplementedError: Subclass has to implement this function.
        """
        raise NotImplementedError
    
    def _hom_to_vector_valued(self, K) :
        """
        This should be a homomorphism of the underlying polynomial rings.
        
        INPUT:
            - `K` -- A ring.
        
        OUTPUT:
            A dictionary with keys the non vector valued generators over `K` and values in the
            set of generators of the vector valued type.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import ModularFormTestType_scalar
            sage: ModularFormTestType_scalar()._hom_to_vector_valued(QQ)
            {g4: g4, g5: g5, g2: g2, g3: g3, g1: g1}
        """
        if self.is_vector_valued() :
            raise ValueError( "This type is already vector valued." )
        
        nvvtype = self.non_vector_valued()
        
        return dict( (nvvtype._generator_by_name(K, g), self._generator_by_name(K, g))
                     for g in self._generator_names(K) )
    
    def graded_submodules_are_free(self, K = None) :
        """
        Whether the modules of elements of fixed grading are free over
        their base ring `K' or over all base rings, respectively.
        
        INPUT:
            - `K` -- A ring or module or None (default: None) 
        
        OUTPUT:
            A boolean.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_types import ModularFormType_abstract
            sage: ModularFormType_abstract().graded_submodules_are_free()
            Traceback (most recent call last):
            ...
            NotImplementedError: Subclass has to implement this function.
        """
        raise NotImplementedError( "Subclass has to implement this function." )

    def has_hecke_action(self) :
        """
        Whether the associated modular forms are equipped with a Hecke action.
        
        OUTPUT:
            A boolean.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_types import ModularFormType_abstract
            sage: ModularFormType_abstract().has_hecke_action()
            False
        """
        return False
    
    def _hecke_operator_class(self) :
        """
        A class that implements the Hecke operation.
        
        OUTPUT:
            An class or function that constructs an instance.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_types import ModularFormType_abstract
            sage: ModularFormType_abstract()._hecke_operator_class()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError
