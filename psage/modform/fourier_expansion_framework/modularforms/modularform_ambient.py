r"""
Rings of orthogonal modular forms.

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

from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_module import GradedExpansionModule_class
from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_ring import GradedExpansionRing_class
from psage.modform.fourier_expansion_framework.modularforms.modularform_element import ModularForm_generic
from psage.modform.fourier_expansion_framework.modularforms.modularform_functor import ModularFormsFunctor
from psage.modform.fourier_expansion_framework.modularforms.modularform_interfaces import ModularFormsAmbientWithHeckeAction_abstract
from psage.modform.fourier_expansion_framework.modularforms.modularform_submodule import ModularFormsWeightSubmodule
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.structure.element import Element

#===============================================================================
# ModularFormsAmbient
#===============================================================================

def ModularFormsAmbient( A, type, precision, *args, **kwds) :
    """
    Create a ring or module of modular forms of given type.  The underlying Fourier
    expansions are calculated up to ``precision``.
    
    INPUT:
        - `A`             -- A ring; The base ring for the modular forms.
        - ``type``        -- An inystance of :class:~`fourier_expansion_framework.modularforms.modularform_types.ModularFormType_abstract`.
        - ``precision``   -- A precision class.
        - ``*arg``        -- Will be forwarded to the type's construction function.
        - ``**kwds``      -- Will be forwarded to the type's construction function.
    
    OUTPUT:
        An instance of :class:~`ModularFormsAmbient_abstract`.
    
    TESTS::
        sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import *
        sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
        sage: ma = ModularFormsAmbient( QQ, ModularFormTestType_scalar(), NNFilter(5) )
        sage: ma = ModularFormsAmbient( QQ, ModularFormTestType_vectorvalued(), NNFilter(5), reduce_before_evaluating = False )
    """
    if not 'reduce_before_evaluating' in kwds :
        kwds['reduce_before_evaluating'] = type.reduce_before_evaluating(A)
        
    return type._ambient_construction_function()(A, type, precision, *args, **kwds)    

#===============================================================================
# ModularFormsAmbient_abstract
#===============================================================================

class ModularFormsAmbient_abstract :
    """
    An abstract implementation of a graded expansion ambient, that deduced its structure from
    data stored by a type of modular forms.
    """
    
    def __init__(self, type, precision) :
        """
        INPUT:
            - ``type``        -- An inystance of :class:~`fourier_expansion_framework.modularforms.modularform_types.ModularFormType_abstract`.
            - ``precision``   -- A precision class.
        
        NOTE:
            - The attribute ``_extended_base_ring`` must be set before calling the constructor or 
              will be ignored.
            - The attribute ``_submodule_classes`` will be overwritten and has to be populated after calling the
              constructor. See :meth:~`._submodule` for its description.
            - The attribute ``_element_class`` may not be set if it has to be adopted from the type. 
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: ma = ModularFormsAmbient_abstract( ModularFormTestType_scalar(), NNFilter(5) )        
        """
        self.__type = type
        self.__precision = precision
        
        if not hasattr(self, '_element_class') :
            try :
                self._element_class = type._ambient_element_class()
            except NotImplementedError :
                self._element_class = ModularForm_generic
                
        single_weight_pred = lambda basis, **kwds: "grading_indices" in kwds and len(kwds["grading_indices"]) == 1
        def single_weight_function(basis, **kwds) :
            try :
                return self.__type._weight_submodule_class()(self, basis, kwds["grading_indices"][0], **kwds)
            except NotImplementedError :
                return ModularFormsWeightSubmodule(self, basis, kwds["grading_indices"][0]) 

        self._submodule_classes = [( single_weight_pred, single_weight_function ),
                                   ( lambda _, **kwds : True,
                                     lambda basis, **kwds : self._graded_ambient_class._submodule(self, basis, **kwds) ) ]
 
#   This couldn't be refactored completely and will most likely cause problems in the 
#   the concrete implementations. A complete replacement has taken place and this is
#   probably no issue.
#   def precision(self) :
    def fourier_expansion_precision(self) :
        """
        The common precision of the underlying Fourier expansions.
        
        OUTPUT:
            A filter for the Fourier expansion ambient's monoid or action.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: ma = ModularFormsAmbient_abstract( ModularFormTestType_scalar(), NNFilter(5) )        
            sage: ma.fourier_expansion_precision()
            Filtered NN with action up to 5
        """
        return self.__precision

    def type(self) :
        """
        The type of modular forms this ambient contains.
        
        OUTPUT:
            An instance of :class:~`fourier_expansion_framework.modularforms.modularform_types.ModularFormType_abstract`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: ma = ModularFormsAmbient_abstract( ModularFormTestType_scalar(), NNFilter(5) )            
            sage: ma.type()
            Test type of modular forms with 5 generators
        """
        return self.__type

    def group(self) :
        """
        The modular group the modular forms in this ambient are attached to.
        
        OUTPUT:
            An arbitrary type.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: ma = ModularFormsAmbient_abstract( ModularFormTestType_scalar(), NNFilter(5) )            
            sage: ma.group()
            '1'
        """
        return self.__type.group()
    
    def weights(self) :
        """
        The generators' weights.
        
        OUTPUT:
            A tuple of (generalized) weights.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: ma = ModularFormsAmbient( QQ, ModularFormTestType_scalar(), NNFilter(5) )
            sage: ma.weights()
            (1, 2, 3, 4, 5)
        """
        return self.__type.weights(self.relations().base_ring())

    def graded_submodules_are_free(self) :
        """
        Whether the modules of elements of fixed grading are free.
        
        OUTPUT:
            A boolean.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: ma = ModularFormsAmbient( QQ, ModularFormTestType_scalar(), NNFilter(5) )
            sage: ma.graded_submodules_are_free()
            True
        """
        return self.__type.graded_submodules_are_free(self.relations().base_ring())

    def _submodule(self, basis, **kwds) :
        """
        A submodule with given basis.
        
        INPUT:
            - ``basis``  -- A list of elements of ``self``.
            - ``**kwds`` -- Will be forwarded to the submodule construction function.

        OUTPUT:
            A submodule of graded expansions.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: ma = ModularFormsAmbient( QQ, ModularFormTestType_scalar(), NNFilter(5) )
            sage: sm = ma._submodule([ma.0], grading_indices = [1])
        """
        for pred, fcn in self._submodule_classes :
            if pred(basis, **kwds) : return fcn(basis, **kwds)
        
        raise RuntimeError("submodule classes do not match {0}, {1}".format(basis, kwds))

    def construction(self):
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: ma = ModularFormsAmbient( QQ, ModularFormTestType_scalar(), NNFilter(5) )
            sage: (F, A) = ma.construction()
            sage: F(A) == ma
            True
        """
        return ModularFormsFunctor(self.__type, self.__precision), self.relations().base_ring()

    def _coerce_map_from_(self, other) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: ma = ModularFormsAmbient( QQ, ModularFormTestType_scalar(), NNFilter(5) )
            sage: ma2 = ModularFormsAmbient( ZZ, ModularFormTestType_scalar(), NNFilter(5) )
            sage: ma._coerce_map_from_(ma2)
            Conversion via _element_constructor_ map:
              From: Graded expansion ring with generators g1, g2, g3, g4, g5
              To:   Graded expansion ring with generators g1, g2, g3, g4, g5
        """
        from sage.structure.coerce_maps import CallableConvertMap
        
        if isinstance(other, ModularFormsAmbient_abstract) and \
           self.relations().base_ring().has_coerce_map_from(other.relations().base_ring()) and \
           self.type() == other.type() :
            return CallableConvertMap(other, self, self._element_constructor_)

        return self._graded_ambient_class._coerce_map_from_(self, other)

    def _element_constructor_(self, x) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: ma = ModularFormsAmbient( QQ, ModularFormTestType_scalar(), NNFilter(5) )
            sage: ma2 = ModularFormsAmbient( ZZ, ModularFormTestType_scalar(), NNFilter(5) )
            sage: ma(ma2.0)
            Graded expansion g1
        """
        if isinstance(x, (int, Integer)) and x == 0 :
            return self._element_class(self, self.relations().ring().zero())
        
        if isinstance(x, Element) :
            P = x.parent()
            if isinstance(P, ModularFormsAmbient_abstract) :
                if P.type() == self.type() and \
                   self.relations().base_ring().has_coerce_map_from(P.relations().base_ring()) :
                    return self._element_class( self,
                            self.relations().ring()( x.polynomial(). \
                              subs(self.type()._hom_base_extension(P.relations().base_ring(), self.relations().base_ring())) )
                           )

        return self._graded_ambient_class._element_constructor_(self, x)
    
    def __cmp__(self, other) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: ma = ModularFormsAmbient( QQ, ModularFormTestType_scalar(), NNFilter(5) )
            sage: ma == ModularFormsAmbient( QQ, ModularFormTestType_scalar(), NNFilter(5) )
            True
            sage: ma == ModularFormsAmbient( ZZ, ModularFormTestType_scalar(), NNFilter(5) )
            False
            sage: ma == ModularFormsAmbient( QQ, ModularFormTestType_vectorvalued(), NNFilter(5) )
            False
        """
        c = cmp(type(self), type(other)) 
        
        if c == 0 :
            c = cmp(self.type(), other.type())
        if c == 0 :
            c = cmp(self.relations().base_ring(), other.relations().base_ring())
        
        return c
    
#===============================================================================
# ModularFormsRing_generic
#===============================================================================

class ModularFormsRing_generic ( ModularFormsAmbient_abstract, GradedExpansionRing_class ) :
    def __init__(self, K, type, precision, **kwds) :
        """
        INPUT:
            - `K`             -- A ring.
            - ``type``        -- An inystance of :class:~`fourier_expansion_framework.modularforms.modularform_types.ModularFormType_abstract`.
            - ``precision``   -- A precision class.
            - ``**kwds``      -- Will be forwardd to :class:~`fourier_expansion_framework.gradedexpansions.gradedexpansion_ring.GradedExpansionRing_class`.
        
        NOTE:
            - The attribute ``_extended_base_ring`` must be set before calling the constructor or 
              will be ignored.
            - The attribute ``_graded_ambient_class`` may not be set.
            - The attribute ``_submodule_classes`` will be overwritten and has to be populated after calling the
              constructor. See :meth:~`fourier_expansion_framework.modularforms.modularform_ambient._submodule` for its description.
            - The attribute ``_element_class`` may not be set if it has to be adopted from the type. 
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: ma = ModularFormsRing_generic( QQ, ModularFormTestType_scalar(), NNFilter(5) )
            sage: ma.has_coerce_map_from(QQ)
            True    
        """
        if not hasattr(self, '_extended_base_ring') :
            try :
                if type.is_vector_valued() :
                    nvv_type = type.non_vector_valued()
                    self._extended_base_ring = ModularFormsAmbient(K, nvv_type, precision)
            except NotImplementedError :
                del self._extended_base_ring
        
        if not hasattr(self, '_graded_ambient_class') :
            self._graded_ambient_class = GradedExpansionRing_class
        
        if not 'all_relations' in kwds :
            kwds['all_relations'] = True
            
        ModularFormsAmbient_abstract.__init__(self, type, precision)
        GradedExpansionRing_class.__init__(self, type.base_ring_generators(K, precision),
         type.generators(K, precision), type.generator_relations(K), type.grading(K), **kwds)

        #=======================================================================
        # self._populate_coercion_lists_(
        #  coerce_list = [GradedExpansionBaseringInjection(self.base_ring(), self)], 
        #  convert_list = [self.relations().ring()],
        #  convert_method_name = "_graded_expansion_submodule_to_graded_ambient_" )
        #=======================================================================
        
    def _coerce_map_from_(self, other):
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: ma = ModularFormsAmbient( QQ, ModularFormTestType_scalar(), NNFilter(5) )
            sage: ma._coerce_map_from_(ma)
            Conversion via _element_constructor_ map:
              From: Graded expansion ring with generators g1, g2, g3, g4, g5
              To:   Graded expansion ring with generators g1, g2, g3, g4, g5
        """
        from sage.structure.coerce_maps import CallableConvertMap
        
        if isinstance(other, ModularFormsRing_generic) and \
           self.base_ring().has_coerce_map_from(other.base_ring()) and \
           self.type().is_vector_valued() and \
           self.type().non_vector_valued() == other.type() :
            return CallableConvertMap(other, self, self._element_constructor_)
        
        return ModularFormsAmbient_abstract._coerce_map_from_(self, other)

    def _element_constructor_(self, x) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: ma = ModularFormsAmbient( QQ, ModularFormTestType_scalar(), NNFilter(5) )
            sage: ma2 = ModularFormsAmbient( ZZ, ModularFormTestType_scalar(), NNFilter(5) )
            sage: ma(ma2.0)
            Graded expansion g1
        """
        if isinstance(x, Element) :
            P = x.parent()
            if isinstance(P, ModularFormsRing_generic) :
                try :
                    if self.type().is_vector_valued() and \
                       self.type().non_vector_values() == P.type() and \
                       self.base_ring().has_coerce_map_from(P.base_ring()) :
                        if self.base_ring() != P.base_ring() :
                            from sage.categories.pushout import pushout
                            x = pushout(P, self.base_ring())(x)
                        
                        return self._element_class( self,
                                self.relations().ring()( x.polynomial(). \
                                 subs(self.type()._hom_to_vector_valued(self.base_ring())) )
                               )
                except NotImplementedError :
                    pass
                
        return ModularFormsAmbient_abstract._element_constructor_(self, x)

#===============================================================================
# ModularFormsModule_generic
#===============================================================================

class ModularFormsModule_generic ( ModularFormsAmbient_abstract, GradedExpansionModule_class ) :
    def __init__(self, K, type, precision, **kwds) :
        """
        INPUT:
            - `K`             -- A ring.
            - ``type``        -- An inystance of :class:~`fourier_expansion_framework.modularforms.modularform_types.ModularFormType_abstract`.
            - ``precision``   -- A precision class.
            - ``**kwds``      -- Will be forwardd to :class:~`fourier_expansion_framework.gradedexpansions.gradedexpansion_ring.GradedExpansionRing_class`.
        
        NOTE:
            - The attribute ``_extended_base_ring`` must be set before calling the constructor or 
              will be ignored.
            - The attribute ``_graded_ambient_class`` may not be set.
            - The attribute ``_submodule_classes`` will be overwritten and has to be populated after calling the
              constructor. See :meth:~`fourier_expansion_framework.modularforms.modularform_ambient._submodule` for its description.
            - The attribute ``_element_class`` may not be set if it has to be adopted from the type. 
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: ma = ModularFormsModule_generic( QQ, ModularFormTestType_vectorvalued(), NNFilter(5) )
        """
        if not hasattr(self, '_extended_base_ring') :
            try :
                if type.is_vector_valued() :
                    nvv_type = type.non_vector_valued()
                    self._extended_base_ring = ModularFormsAmbient(K, nvv_type, precision)
            except NotImplementedError :
                del self._extended_base_ring

        if not hasattr(self, '_graded_ambient_class') :
            self._graded_ambient_class = GradedExpansionModule_class
        
        if not 'all_relations' in kwds :
            kwds['all_relations'] = True
            
        ModularFormsAmbient_abstract.__init__(self, type, precision)
        GradedExpansionModule_class.__init__(self, type.base_ring_generators(K, precision),
         type.generators(K, precision), type.generator_relations(K), type.grading(K), **kwds)

        #=======================================================================
        # self._populate_coercion_lists_(
        #  convert_list = [self.relations().ring()],
        #  convert_method_name = "_graded_expansion_submodule_to_graded_ambient_" )
        #=======================================================================

#===============================================================================
# ModularFormsRing_withheckeaction
#===============================================================================

class ModularFormsRing_withheckeaction(ModularFormsAmbientWithHeckeAction_abstract, ModularFormsRing_generic ) :
    def __init__(self, K, type, precision, **kwds) :
        """
        INPUT:
            - `K`             -- A ring.
            - ``type``        -- An inystance of :class:~`fourier_expansion_framework.modularforms.modularform_types.ModularFormType_abstract`.
            - ``precision``   -- A precision class.
            - ``**kwds``      -- Will be forwardd to :class:~`fourier_expansion_framework.gradedexpansions.gradedexpansion_ring.GradedExpansionRing_class`.
        
        NOTE:
            - The attribute ``_extended_base_ring`` must be set before calling the constructor or 
              will be ignored.
            - The attribute ``_graded_ambient_class`` may not be set.
            - The attribute ``_submodule_classes`` will be overwritten and has to be populated after calling the
              constructor. See :meth:~`fourier_expansion_framework.modularforms.modularform_ambient._submodule` for its description.
            - The attribute ``_element_class`` may not be set if it has to be adopted from the type. 
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: ma = ModularFormsRing_generic( QQ, ModularFormTestType_scalar(), NNFilter(5) )    
        """
        ModularFormsRing_generic.__init__(self, K, type, precision, **kwds)
        ModularFormsAmbientWithHeckeAction_abstract.__init__(self, type)

#===============================================================================
# ModularFormsModule_withheckeaction
#===============================================================================

class ModularFormsModule_withheckeaction(ModularFormsAmbientWithHeckeAction_abstract, ModularFormsModule_generic ) :
    def __init__(self, K, type, precision, **kwds) :
        """
        INPUT:
            - `K`             -- A ring.
            - ``type``        -- An inystance of :class:~`fourier_expansion_framework.modularforms.modularform_types.ModularFormType_abstract`.
            - ``precision``   -- A precision class.
            - ``**kwds``      -- Will be forwardd to :class:~`fourier_expansion_framework.gradedexpansions.gradedexpansion_ring.GradedExpansionRing_class`.
        
        NOTE:
            - The attribute ``_extended_base_ring`` must be set before calling the constructor or 
              will be ignored.
            - The attribute ``_graded_ambient_class`` may not be set.
            - The attribute ``_submodule_classes`` will be overwritten and has to be populated after calling the
              constructor. See :meth:~`fourier_expansion_framework.modularforms.modularform_ambient._submodule` for its description.
            - The attribute ``_element_class`` may not be set if it has to be adopted from the type. 
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: ma = ModularFormsModule_withheckeaction( QQ, ModularFormTestType_vectorvalued(), NNFilter(5) )
        """
        ModularFormsModule_generic.__init__(self, K, type, precision, **kwds)
        ModularFormsAmbientWithHeckeAction_abstract.__init__(self, type)
        