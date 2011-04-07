"""
The implementation of a type of modular forms, that implements all potential features.
It is meant to be used in doctests.

AUTHOR :
    - Martin Raum (2010 - 09 - 26) Initial version.
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

from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import ModularFormsRing_generic,\
    ModularFormsModule_generic
from psage.modform.fourier_expansion_framework.modularforms.modularform_types import ModularFormType_abstract
from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid, \
                                              NNFilter, TrivialCharacterMonoid, TrivialRepresentation
from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import EquivariantMonoidPowerSeries
from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_module import EquivariantMonoidPowerSeriesModule
from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import EquivariantMonoidPowerSeriesRing
from operator import xor
from sage.misc.cachefunc import cached_method
from sage.modules.all import FreeModule
from sage.rings.all import Integer
from sage.rings.all import ZZ, PolynomialRing
from sage.structure.all import Sequence, SageObject

#===============================================================================
# ModularFormTestType_scalar
#===============================================================================

class ModularFormTestType_scalar ( ModularFormType_abstract ) :

    nmb_gens = 5
    
    def __init__(self) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: t = ModularFormTestType_scalar()
        """
        ModularFormType_abstract.__init__(self)
        
    def _ambient_construction_function(self) :
        """
        Return a function which will can be called by :function:~`fourier_expansion_framework.modularforms.modularform_ambient.ModularFormsAmbient`.
        
        OUTPUT:
            A function.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: t = ModularFormTestType_scalar()
            sage: h = t._ambient_construction_function()(QQ, t, NNFilter(4))
        """
        return lambda A, type, precision, **kwds : \
                 ModularFormsRing_generic(A, type,
                    NNFilter(precision) if isinstance(precision, (int, Integer)) else precision, **kwds)

    def group(self) :
        """
        Return the modular group this type. Here it will act trivially
        
        OUTPUT:
            A string.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: t = ModularFormTestType_scalar()
            sage: t.group()
            '1'
        """
        return "1"
    
    @cached_method
    def _g(self, i, precision) :
        """
        Return the Fourier expansion of the i-th generator.
        
        INPUT:
            - ``precision`` -- An instance of :class:~`fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids.NNFilter`.
            
        OUTPUT:
            An element of :class:`~fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring.EquivariantMonoidPowerSeriesRing`.

        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: t = ModularFormTestType_scalar()
            sage: t._g(2, NNFilter(3)).coefficients()
            {2: 1}
        """
        ea = TestExpansionAmbient(ZZ)
        return EquivariantMonoidPowerSeries(ea , {ea.characters().one_element(): {i: 1}}, NNFilter(precision) )
    
    def generators(self, K, precision) :
        """
        The generators for modular forms with coefficients in `K`.
        
        INPUT:
            - `K`      -- A ring.
            - ``precision`` -- An instance of :class:~`fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids.NNFilter`.

        OUTPUT:
            A sequence of elements of :class:`~fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring.EquivariantMonoidPowerSeriesRing`. 

        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: t = ModularFormTestType_scalar()
            sage: t.generators(QQ, NNFilter(4))
            [Equivariant monoid power series in Ring of equivariant monoid power series over NN, Equivariant monoid power series in Ring of equivariant monoid power series over NN, Equivariant monoid power series in Ring of equivariant monoid power series over NN, Equivariant monoid power series in Ring of equivariant monoid power series over NN, Equivariant monoid power series in Ring of equivariant monoid power series over NN]
        """
        if K.has_coerce_map_from(ZZ) :
            return Sequence( [ self._g(i, precision) for i in range(self.nmb_gens) ],
                             universe = TestExpansionAmbient(K) )
                             
        raise NotImplementedError

    def grading(self, K) :
        """
        The weight grading of the underlying ring of modular form with
        coefficients in `K`.
        
        INPUT:
            - `K` -- A ring.
        
        OUTPUT:
            An instance of :class:`~fourier_expansion_framework.gradedexpansions.gradedexpansion_grading.Grading_abstract`.

        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: t = ModularFormTestType_scalar()
            sage: t.grading(QQ)
            Degree grading (1, 2, 3, 4, 5)
        """
        if K.has_coerce_map_from(ZZ) :
            return DegreeGrading(range(1, self.nmb_gens + 1))

        raise NotImplementedError

    def _generator_names(self, K) :
        """
        The generators' names for modular forms with coefficients in `K`.
        
        INPUT:
            - `K` -- A ring.
        
        OUTPUT:
            A list of strings.

        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: t = ModularFormTestType_scalar()
            sage: t._generator_names(QQ)
            ['g1', 'g2', 'g3', 'g4', 'g5']
        """
        return ["g%s" % (i,) for i in range(1, self.nmb_gens + 1)]
    
    def _generator_by_name(self, K, name) :
        """
        Given a name return the associated generator of modular forms with
        coefficients in `K` in the underlying polynomial algebra.
        
        INPUT:
            - `K`      -- A ring.
            - ``name`` -- A string. The generator's name.
            
        OUTPUT:
            An element of a polynomial ring.
            
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: t = ModularFormTestType_scalar()
            sage: t._generator_by_name(QQ, 'g3')
            g3
        """
        if K.has_coerce_map_from(ZZ) :
            R = self.generator_relations(K).ring()
            try :
                return R.gens()[self._generator_names(K).index(name)]
            except ValueError :
                raise ValueError( "Generator name %s doesn't exist for %s" % (name, K))
            
        raise NotImplementedError

    @cached_method
    def generator_relations(self, K) :
        """
        An ideal `I` in a polynomial ring `R` over `K`, such that the
        associated ring is `R / I` surjects onto the ring of modular forms
        with coefficients in `K`.
        
        INPUT:
            - `K` -- A ring.
            
        OUTPUT:
            An ideal in a polynomial ring.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: t = ModularFormTestType_scalar()
            sage: t.generator_relations(QQ)
            Ideal (g1^2 - g2, g1^3 - g3, g1^4 - g4, g1^5 - g5) of Multivariate Polynomial Ring in g1, g2, g3, g4, g5 over Rational Field
        """
        if K.has_coerce_map_from(ZZ) :
            R = PolynomialRing(K, self._generator_names(K))
            g1 = R.gen(0)
            return R.ideal([g1**i - g for (i,g) in list(enumerate([None] + list(R.gens())))[2:]])
            
        raise NotImplementedError

    def weights(self, K) :
        """
        The weights of the generators of the ring of modular forms
        with coefficients in `K`.
        
        INPUT:
            - `K` -- A ring.
            
        OUTPUT:
            A tuple of weights.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: t = ModularFormTestType_scalar()
            sage: t.weights(QQ)
            (1, 2, 3, 4, 5)
        """
        return self.grading(K).gens()

    def non_vector_valued(self) :
        """
        Return the non vector values version of this type.
        
        OUTPUT:
            An instance of :class:`~fourier_expansion_framework.modularforms.modularform_types.ModularFormType_abstract`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: t = ModularFormTestType_scalar()
            sage: t.non_vector_valued() is t
            True
        """
        return self

    def vector_valued(self) :
        """
        Return the vector values version of this type.
        
        OUTPUT:
            An instance of :class:`~fourier_expansion_framework.modularforms.modularform_types.ModularFormType_abstract`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: t = ModularFormTestType_scalar()
            sage: t.vector_valued()
            Test type of vector valued modular forms of rank 3
        """
        return ModularFormTestType_vectorvalued()
    
    def graded_submodules_are_free(self, K = None) :
        """
        Return True if all submodules of forms of fixed grading according to
        the grading given by :meth:~`.grading` are free.
        
        OUTPUT:
            A boolean.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: t = ModularFormTestType_scalar()
            sage: t.graded_submodules_are_free()
            True
        """
        return True

    def __cmp__(self, other) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: t = ModularFormTestType_scalar()
            sage: t == ModularFormTestType_scalar()
            True
            sage: t == ModularFormTestType_vectorvalued()
            False
        """
        return cmp(type(self), type(other))

    def __hash__(self) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: hash( ModularFormTestType_scalar() )
            ?? # 32-bit
            -4492425810583750348 # 64-bit
        """
        return reduce(xor, map(hash, [self.nmb_gens, self._repr_()] ) )
    
    def _repr_(self) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: ModularFormTestType_scalar()
            Test type of modular forms with 5 generators
        """
        return "Test type of modular forms with %s generators" % (self.nmb_gens,)


class ModularFormTestType_vectorvalued ( ModularFormType_abstract ) :

    rank = 3
    
    def __init__(self) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: t = ModularFormTestType_vectorvalued()
        """
        ModularFormType_abstract.__init__(self)
        
    def _ambient_construction_function(self) :
        """
        Return a function which will can be called by :function:~`fourier_expansion_framework.modularforms.modularform_ambient.ModularFormsAmbient`.
        
        OUTPUT:
            A function.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: t = ModularFormTestType_vectorvalued()
            sage: h = t._ambient_construction_function()(QQ, t, NNFilter(4, False))
        """
        return lambda A, type, precision, **kwds : \
                 ModularFormsModule_generic(A, type,
                    NNFilter(precision) if isinstance(precision, (int, Integer)) else precision, **kwds)

    def group(self) :
        """
        Return the modular group this type. Here it will act trivially
        
        OUTPUT:
            A string.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: t = ModularFormTestType_vectorvalued()
            sage: t.group()
            '1'
        """
        return "1"
    
    @cached_method
    def _g(self, i, precision) :
        """
        Return the Fourier expansion of the i-th generator.
        
        INPUT:
            - ``precision`` -- An instance of :class:~`fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids.NNFilter`.
            
        OUTPUT:
            An element of :class:`~fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring.EquivariantMonoidPowerSeriesModule`.

        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: t = ModularFormTestType_vectorvalued()
            sage: t._g(2, NNFilter(3)).coefficients()
            {1: (0, 0, 1)}
        """
        if i < self.rank :
            ea = TestExpansionAmbient_vv(ZZ)
            cd = ea.coefficient_domain()
            v = [0 for _ in range(self.rank)]
            v[i] = 1
            return EquivariantMonoidPowerSeries(ea , {ea.characters().one_element(): {1: cd(v)}}, NNFilter(precision) )
        else :
            raise ValueError( "%s-th generator is not definied" % (i,) )
    
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
        return self.non_vector_valued().generators(K, precision)
    
    def generators(self, K, precision) :
        """
        The generators for modular forms with coefficients in `K`.
        
        INPUT:
            - `K`           -- A ring.
            - ``precision`` -- An instance of :class:~`fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids.NNFilter`.
        
        OUTPUT:
            A sequence of elements of :class:`~fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring.EquivariantMonoidPowerSeriesModule`. 

        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: t = ModularFormTestType_vectorvalued()
            sage: t.generators(QQ, NNFilter(4))
             [Equivariant monoid power series in Module of equivariant monoid power series over NN, Equivariant monoid power series in Module of equivariant monoid power series over NN, Equivariant monoid power series in Module of equivariant monoid power series over NN]
        """
        if K.has_coerce_map_from(ZZ) :
            return Sequence( [ self._g(i, precision) for i in range(self.rank) ],
                             universe = TestExpansionAmbient_vv(K) )
                             
        raise NotImplementedError

    def grading(self, K) :
        """
        The weight grading of the underlying ring of modular form with
        coefficients in `K`.
        
        INPUT:
            K - A ring.
        
        OUTPUT:
            An instance of :class:`~fourier_expansion_framework.gradedexpansions.gradedexpansion_grading.Grading_abstract`.

        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: t = ModularFormTestType_vectorvalued()
            sage: t.grading(QQ)
            Degree grading (1, 2, 3, 4, 5, 3, 6, 9)
        """
        if K.has_coerce_map_from(ZZ) :
            return DegreeGrading(list(self.non_vector_valued().weights(K)) + range(3, 3 * self.rank + 1, 3))

        raise NotImplementedError

    def _generator_names(self, K) :
        """
        The generators' names for modular forms with coefficients in `K`.
        
        INPUT:
            - `K` -- A ring.
        
        OUTPUT:
            A list of strings.

        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: t = ModularFormTestType_vectorvalued()
            sage: t._generator_names(QQ)
            ['v1', 'v2', 'v3']
        """
        return ["v%s" % (i,) for i in range(1, self.rank + 1)]
    
    def _generator_by_name(self, K, name) :
        """
        Given a name return the associated generator of modular forms with
        coefficients in `K` in the underlying polynomial algebra.
        
        INPUT:
            - `K`      -- A ring.
            - ``name`` -- A string. The generator's name.
            
        OUTPUT:
            An element of a polynomial ring.
            
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: t = ModularFormTestType_vectorvalued()
            sage: t._generator_by_name(QQ, 'v3')
            v3
        """
        if K.has_coerce_map_from(ZZ) :
            R = self.generator_relations(K).ring()
            try :
                return R.gens()[len(self.non_vector_valued()._generator_names(K)) + self._generator_names(K).index(name)]
            except ValueError :
                raise ValueError( "Generator name %s doesn't exist for %s" % (name, K))
            
        raise NotImplementedError

    @cached_method
    def generator_relations(self, K) :
        """
        An ideal `I` in a polynomial ring `R` over `K`, such that the
        associated ring is `R / I` surjects onto the ring of modular forms
        with coefficients in `K`.
        
        INPUT:
            - `K` -- A ring.
            
        OUTPUT:
            An ideal in a polynomial ring.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: t = ModularFormTestType_vectorvalued()
            sage: t.generator_relations(QQ)
            Ideal (g1^2 - g2, g1^3 - g3, g1^4 - g4, g1^5 - g5) of Multivariate Polynomial Ring in g1, g2, g3, g4, g5, v1, v2, v3 over Rational Field
        """
        if K.has_coerce_map_from(ZZ) :
            R = PolynomialRing(K, self.non_vector_valued()._generator_names(K) + self._generator_names(K))
            return R.ideal().parent()(self.non_vector_valued().generator_relations(K))
            
        raise NotImplementedError

    def weights(self, K) :
        """
        The weights of the generators of the ring of modular forms
        with coefficients in `K`.
        
        INPUT:
            - `K` -- A ring.
            
        OUTPUT:
            A tuple of weights.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: t = ModularFormTestType_vectorvalued()
            sage: t.weights(QQ)
            (1, 2, 3, 4, 5, 3, 6, 9)
        """
        return self.grading(K).gens()

    def is_vector_valued(self) :
        """
        ``True`` if this is the vector valued version of a scalar valued type of modular forms.
        
        OUTPUT:
            A boolean.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: t = ModularFormTestType_vectorvalued()
            sage: t.is_vector_valued()
            True
        """
        return True
    
    def non_vector_valued(self) :
        """
        Return the non vector values version of this type.
        
        OUTPUT:
            An instance of :class:`~fourier_expansion_framework.modularforms.modularform_types.ModularFormType_abstract`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: t = ModularFormTestType_scalar()
            sage: t.non_vector_valued()
            Test type of modular forms with 5 generators
        """
        return ModularFormTestType_scalar()

    def vector_valued(self) :
        """
        Return the vector values version of this type.
        
        OUTPUT:
            An instance of :class:`~fourier_expansion_framework.modularforms.modularform_types.ModularFormType_abstract`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: t = ModularFormTestType_vectorvalued()
            sage: t.vector_valued() is t
            True
        """
        return self
    
    def graded_submodules_are_free(self, K = None) :
        """
        Return True if all submodules of forms of fixed grading according to
        the grading given by :meth:~`.grading` are free.
        
        OUTPUT:
            A boolean.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: t = ModularFormTestType_vectorvalued()
            sage: t.graded_submodules_are_free()
            True
        """
        return True
    
    def has_hecke_action(self) :
        """
        Whether the associated modular forms are equipped with a Hecke action.
        
        OUTPUT:
            A boolean.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: t = ModularFormTestType_vectorvalued()
            sage: t.has_hecke_action()
            True
        """
        return True
    
    def _hecke_operator_class(self) :
        """
        A class that implements the Hecke operation.
        
        OUTPUT:
            An class or function that constructs an instance.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: t = ModularFormTestType_vectorvalued()
            sage: t._hecke_operator_class()(3)
            Test Hecke operator with modulus 3
        """
        return HeckeOperatorNN_test
    
    def __cmp__(self, other) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: t = ModularFormTestType_vectorvalued()
            sage: t == ModularFormTestType_vectorvalued()
            True
            sage: t == ModularFormTestType_scalar()
            False
        """
        return cmp(type(self), type(other))
    
    def __hash__(self) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: hash( ModularFormTestType_vectorvalued() )
            ?? # 32-bit
            -3841460515652797985 # 64-bit
        """
        return reduce(xor, map(hash, [self.rank, self._repr_()] ) )
    
    def _repr_(self) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: ModularFormTestType_vectorvalued()
            Test type of vector valued modular forms of rank 3
        """
        return "Test type of vector valued modular forms of rank %s" % (self.rank,)
    

def TestExpansionAmbient(A) :
    """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: ea = TestExpansionAmbient(QQ)    
    """
    return EquivariantMonoidPowerSeriesRing( NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", A) )

def TestExpansionAmbient_vv(A) :
    """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: ea = TestExpansionAmbient_vv(QQ)    
    """
    return EquivariantMonoidPowerSeriesModule( NNMonoid(), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", FreeModule(A, ModularFormTestType_vectorvalued.rank)) )

class HeckeOperatorNN_test ( SageObject ) :
    
    def __init__(self, l) :
        """
        INPUT:
            - `l` -- An integer.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: h = HeckeOperatorNN_test(2)
        """
        self.__l = l
    
    def eval(self, expansion, weight = None) :
        """
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import *
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: h = HeckeOperatorNN_test(2)
            sage: ma = ModularFormsModule_withheckeaction( QQ, ModularFormTestType_vectorvalued(), NNFilter(5) )
            sage: h.eval(ma.0.fourier_expansion()).coefficients()
            {1: {1: (1, 0, 0), 3: (1, 0, 0)}}
        """
        precision = expansion.precision()
        if precision.is_infinite() :
            precision = expansion._bounding_precision()
        characters = expansion.non_zero_components()


        hecke_expansion = dict()
        for ch in characters :
            res = dict()
            for (n, v) in expansion.coefficients(True)[ch].iteritems() :
                for m in range(n, precision.index(), self.__l) :
                    try :
                        res[m] += v
                    except KeyError :
                        res[m] = v
            
            hecke_expansion[ch] = res

        result = expansion.parent()._element_constructor_(hecke_expansion)
        result._set_precision(expansion.precision())

        return result

    def _repr_(self) :
        """
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.modularforms.modularform_testtype import *
            sage: HeckeOperatorNN_test(2)
            Test Hecke operator with modulus 2
        """
        return "Test Hecke operator with modulus %s" % (self.__l,)
    