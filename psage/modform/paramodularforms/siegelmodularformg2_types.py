r"""
Types of Siegel modular forms of genus 2.

AUTHORS:

- Martin Raum (2009 - 08 - 03) Initial version.
- Martin Raum (2010 - 03 - 16) Added types for vector valued forms.
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

from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import DegreeGrading
from psage.modform.fourier_expansion_framework.modularforms.modularform_ambient import ModularFormsRing_withheckeaction,\
                                             ModularFormsModule_generic
from psage.modform.fourier_expansion_framework.modularforms.modularform_types import ModularFormType_abstract
from operator import xor
from sage.categories.number_fields import NumberFields
from sage.misc.cachefunc import cached_method
from sage.modular.modform.constructor import ModularForms
from sage.arith.all import bernoulli
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.structure.sequence import Sequence
from psage.modform.paramodularforms.siegelmodularformg2_element import SiegelModularFormG2_classical, SiegelModularFormG2_vectorvalued
from psage.modform.paramodularforms.siegelmodularformg2_fegenerators import SiegelModularFormG2MaassLift
from psage.modform.paramodularforms.siegelmodularformg2vv_fegenerators import SiegelModularFormG2SatohBracket
from psage.modform.paramodularforms.siegelmodularformg2_fourierexpansion import SiegelModularFormG2Filter_discriminant
from psage.modform.paramodularforms.siegelmodularformg2_heckeaction import SiegelModularFormG2FourierExpansionHeckeAction
from psage.modform.paramodularforms.siegelmodularformg2_submodule import SiegelModularFormG2WeightSubmodule_class, \
                                          SiegelModularFormG2SubmoduleVector_generic

#===============================================================================
# SiegelModularFormsG2
#===============================================================================

_siegelmodularforms_cache = dict()

def SiegelModularFormsG2(A, type, precision, *args, **kwds) :
    r"""
    TESTS:
    
    Some brief tests to check that the classical modular forms and the Hecke action work.
    
    ::
    
        sage: from psage.modform.paramodularforms import SiegelModularFormsG2, SiegelModularFormG2_Classical_Gamma
        sage: sr = SiegelModularFormsG2(QQ, SiegelModularFormG2_Classical_Gamma(), 200)
        sage: sr.0         
        Graded expansion I4
        sage: sr.gens()   
        (Graded expansion I4, Graded expansion I6, Graded expansion I10, Graded expansion I12)
        sage: sr.2.fourier_expansion().coefficients()[(1,1,32)]
        319360
        sage: sr.weights()
        [4, 6, 10, 12]
        sage: map(sr, sr.graded_submodule(20).basis()) 
        [Graded expansion I4^5, Graded expansion I4^2*I6^2, Graded expansion I4^2*I12, Graded expansion I4*I6*I10, Graded expansion I10^2]
        sage: sr.graded_submodule(26).hecke_eigenforms(2)
        [Graded expansion I4^5*I6 + 265000/392931*I4^2*I6^3 - 1437170847265656463317060480000/19802288209643185928499101*I4^4*I10 - 71078691383556603157606297600000/1524776192142525316494430777*I4*I6^2*I10 - 1984664637742925286551976960000/19802288209643185928499101*I4^2*I6*I12 + 175115938885406662127531537203200000/217825170306075045213490111*I6*I10^2 + 75976549463485417333446829670400000/19802288209643185928499101*I4*I10*I12, Graded expansion I4^4*I10 - 1255/973*I4*I6^2*I10 + 390/973*I4^2*I6*I12 - 4438886400/973*I6*I10^2 + 561116160/139*I4*I10*I12, Graded expansion I4^4*I10 - 31/22*I4*I6^2*I10 + 3/22*I4^2*I6*I12 + 6903360/11*I6*I10^2 + 48304512/11*I4*I10*I12, Graded expansion I4^5*I6 - I4^2*I6^3 - 60187630099968/579984481*I4^4*I10 + 660226634035200/6379829291*I4*I6^2*I10 + 310518955994112/6379829291*I4^2*I6*I12 - 15571964273098752000/6379829291*I6*I10^2 - 14567211596670566400/6379829291*I4*I10*I12, Graded expansion I4^4*I10 + (3551/1520783457816268800*a4^2 - 7540557907/31682988704505600*a4 + 10297152853879/2750259436155)*I4*I6^2*I10 + (-13/3448488566476800*a4^2 + 27240601/71843511801600*a4 - 45812312422/6236415955)*I4^2*I6*I12 + (-98507/2095435760880*a4^2 + 204121897999/43654911685*a4 - 766406631173075712/8730982337)*I6*I10^2 + (6301/997826552800*a4^2 - 42899467371/62364159550*a4 + 89460377464535424/6236415955)*I4*I10*I12]
    
    The vector valued Siegel modular forms are implemented quite roughly.
    
    ::
    
        sage: from psage.modform.paramodularforms import SiegelModularFormsG2, SiegelModularFormG2_VectorValuedW2_Gamma
        sage: sf = SiegelModularFormsG2(QQ, SiegelModularFormG2_VectorValuedW2_Gamma(), 20)
        sage: sf.0
        Graded expansion vector (1, 0, 0, 0, 0, 0)
        sage: sf.0.fourier_expansion().coefficients()
        {(1, 0, 3): 56198016*x^2 - 370758528*y^2, (0, 0, 9): -16364592*x^2, (0, 0, 8): 12165120*x^2, (0, 0, 11): 76984128*x^2, (1, 0, 1): -30240*x^2 - 30240*y^2, (0, 0, 12): -53415936*x^2, (0, 0, 7): -2411136*x^2, (1, 0, 4): 910465920*x^2 - 4281167520*y^2, (0, 0, 6): -870912*x^2, (0, 0, 5): 695520*x^2, (0, 0, 4): -211968*x^2, (0, 0, 3): 36288*x^2, (0, 0, 2): -3456*x^2, (2, 0, 2): 202023360*x^2 + 202003200*y^2, (0, 0, 1): 144*x^2, (0, 0, 0): 0, (2, 1, 2): 104686848*x^2 + 122863104*x*y + 108476928*y^2, (1, 1, 1): -4032*x^2 - 4032*x*y - 4032*y^2, (2, 2, 2): 1886976*x^2 - 169344*x*y + 35332416*y^2, (1, 1, 3): -22680000*x^2 - 176904000*x*y - 176904000*y^2, (1, 1, 4): -104799744*x^2 - 2477454336*x*y - 2477454336*y^2, (1, 1, 2): -317952*x^2 - 3787776*x*y - 3787776*y^2, (1, 1, 5): -76507200*x^2 - 18406059840*x*y - 18406059840*y^2, (1, 0, 2): 1669248*x^2 - 11866176*y^2, (0, 0, 10): -16692480*x^2}
    """
    global _siegelmodularforms_cache
    k = (A, type, precision)
    try :
        return _siegelmodularforms_cache[k]
    except KeyError :
        if isinstance(type, SiegelModularFormG2_Classical_Gamma) :
            precision = SiegelModularFormG2Filter_discriminant(precision)
            R = ModularFormsRing_withheckeaction(A, type, precision)
        elif isinstance(type, SiegelModularFormG2_VectorValuedW2_Gamma) :
            precision = SiegelModularFormG2Filter_discriminant(precision)
            R = ModularFormsModule_generic(A, type, precision)
        else :
            raise TypeError("{0} must be an Siegel modular form type".format(type))
                
        _siegelmodularforms_cache[k] = R
        return R


#===============================================================================
# SiegelModularFormG2_Classical_Gamma
#===============================================================================

class SiegelModularFormG2_Classical_Gamma ( ModularFormType_abstract ) :
    def __init__(self) :
        ModularFormType_abstract.__init__(self)
    
    def _ambient_construction_function(self) :
        return SiegelModularFormsG2
    
    def _ambient_element_class(self) :
        return SiegelModularFormG2_classical
    
    def _space_element_class(self) :
        return SiegelModularFormG2SubmoduleVector_generic
    
    def _weight_submodule_class(self) :
        return SiegelModularFormG2WeightSubmodule_class
        
    def group(self) :
        return "Sp(2,ZZ)"

    def _I4(self, precision) :
        E4 = ModularForms(1,4).gen(0)
            
        return SiegelModularFormG2MaassLift(lambda p: 60*(E4.qexp(p)), 0, precision, True, weight = 4)

    def _I6(self, precision) :
        E6 = ModularForms(1,6).gen(0)

        return SiegelModularFormG2MaassLift(lambda p: -84*(E6.qexp(p)), 0, precision, True, weight = 6)

    def _I10(self, precision) :
        # we use a standard generator, since its evaluation is much faster
        Delta = ModularForms(1,12).gen(0)
        assert Delta == ModularForms(1,12).cuspidal_subspace().gen(0)
        
        return SiegelModularFormG2MaassLift(0, lambda p: -(Delta.qexp(p)), precision, True, weight = 10)
        
    def _I12(self, precision) :
        Delta = ModularForms(1,12).gen(0)
        assert Delta == ModularForms(1,12).cuspidal_subspace().gen(0)

        return SiegelModularFormG2MaassLift(lambda p: Delta.qexp(p), 0, precision, True, weight = 12)

    def generators(self, K, precision) :
        if K is QQ or K in NumberFields() :            
            return Sequence([ self._I4(precision), self._I6(precision),
                              self._I10(precision), self._I12(precision) ])
            
        raise NotImplementedError

    def grading(self, K) :
        if K is QQ or K in NumberFields() :
            return DegreeGrading([4,6,10,12])
        
        raise NotImplementedError    

    def __maass_lifts(self, k, precision, return_value) :
        r"""
        Return the Fourier expansion of all Maass forms of weight `k`.
        """
        result = []
        
        if k < 4 or k % 2 != 0 :
            return []
                
        mf = ModularForms(1,k).echelon_basis()
        cf = ModularForms(1,k + 2).echelon_basis()[1:]
        integrality_factor = 2*k * bernoulli(k).denominator()

        for c in [(integrality_factor * mf[0],0)] \
                     + [ (f,0) for f in mf[1:] ] + [ (0,g) for g in cf ] :
            if return_value == "lifts" :
                result.append(SiegelModularFormG2MaassLift(c[0],c[1], precision, True))
            else :
                result.append(c)
        
        return result
        
    def _maass_generators(self, k, precision) :
        return self.__maass_lifts(k, precision, "lifts")
    
    def _maass_generator_preimages(self, k) :
        return self.__maass_lifts(k, 0, "preimages")        
    
    def _generator_names(self, K) :
        return ["I4", "I6", "I10", "I12"]
    
    def _generator_by_name(self, K, name) :
        if K is QQ or K in NumberFields() :
            R = self.generator_relations(K).ring()
            if name == "I4" : return R.gen(0)
            elif name == "I6" : return R.gen(1)
            elif name == "I10" : return R.gen(2)
            elif name == "I12" : return R.gen(3)
            
            raise ValueError("name {0} doesn't exist for {1}".format(name, K))
            
        raise NotImplementedError
        
    @cached_method
    def generator_relations(self, K) :
        r"""
        An ideal `I` in a polynomial ring `R`, such that the associated ring
        is `R / I`. This ideal must be unique for `K`. 
        """
        if K is QQ or K in NumberFields() :
            R = PolynomialRing(K, self._generator_names(K))
            return R.ideal(0)
            
        raise NotImplementedError
    
    def weights(self, K) :
        r"""
        A list of integers corresponding to the weights.
        """
        if K is QQ or K in NumberFields() :
            return [4,6,10,12]
        
        raise NotImplementedError

    def non_vector_valued(self) :
        r"""
        Return the non vector values version of this type. 
        """
        return self
        
    def vector_valued(self) :
        r"""
        Return the vector values version of this type.
        """
        raise NotImplementedError    

    def has_hecke_action(self) :
        return True
    
    def graded_submodules_are_free(self) :
        return True
    
    def _hecke_operator_class(self) :
        return SiegelModularFormG2FourierExpansionHeckeAction
    
    def __cmp__(self, other) :
        return cmp(type(self), type(other))

    def __hash__(self) :
        return xor(hash(type(self)), hash(self.grading(QQ)))

#===============================================================================
# SiegelModularFormG2_VectorValuedW2_Gamma
#===============================================================================

class SiegelModularFormG2_VectorValuedW2_Gamma ( ModularFormType_abstract ) :
    def __init__(self) :
        ModularFormType_abstract.__init__(self)
    
    def _ambient_construction_function(self) :
        return SiegelModularFormsG2 
 
    def _ambient_element_class(self) :
        return SiegelModularFormG2_vectorvalued
    
    def _space_element_class(self) :
        return SiegelModularFormG2SubmoduleVector_generic
    
    def group(self) :
        return "Sp(2,ZZ)"

    def _I4(self, precision) :
        E4 = ModularForms(1,4).gen(0)
            
        return SiegelModularFormG2MaassLift(lambda p: 60*(E4.qexp(p)), 0, precision, True, weight = 4)
 
    def _I6(self, precision) :
        E6 = ModularForms(1,6).gen(0)
 
        return SiegelModularFormG2MaassLift(lambda p: -84*(E6.qexp(p)), 0, precision, True, weight = 6)
 
    def _I10(self, precision) :
        # we use a standard generator, since its evaluation is much faster
        Delta = ModularForms(1,12).gen(0)
        assert Delta == ModularForms(1,12).cuspidal_subspace().gen(0)
        
        return SiegelModularFormG2MaassLift(0, lambda p: -(Delta.qexp(p)), precision, True, weight = 10)
        
    def _I12(self, precision) :
        Delta = ModularForms(1,12).gen(0)
        assert Delta == ModularForms(1,12).cuspidal_subspace().gen(0)
 
        return SiegelModularFormG2MaassLift(lambda p: Delta.qexp(p), 0, precision, True, weight = 12)

    def _satoh_I4_I6(self, precision) :
        return SiegelModularFormG2SatohBracket(self._I4(precision), self._I6(precision), 4, 6)

    def _satoh_I4_I10(self, precision) :
        return SiegelModularFormG2SatohBracket(self._I4(precision), self._I10(precision), 4, 10)
    
    def _satoh_I4_I12(self, precision) :
        return SiegelModularFormG2SatohBracket(self._I4(precision), self._I12(precision), 4, 12)

    def _satoh_I6_I10(self, precision) :
        return SiegelModularFormG2SatohBracket(self._I6(precision), self._I10(precision), 6, 10)

    def _satoh_I6_I12(self, precision) :
        return SiegelModularFormG2SatohBracket(self._I6(precision), self._I12(precision), 6, 12)

    def _satoh_I10_I12(self, precision) :
        return SiegelModularFormG2SatohBracket(self._I10(precision), self._I12(precision), 10, 12)

    def base_ring_generators(self, K, precision) :
        if K is QQ or K in NumberFields() :            
            return Sequence([ self._I4(precision), self._I6(precision),
                              self._I10(precision), self._I12(precision) ])
        raise NotImplementedError

    def generators(self, K, precision) :
        if K is QQ or K in NumberFields() :
            return Sequence([ self._satoh_I4_I6(precision), self._satoh_I4_I10(precision),
                              self._satoh_I4_I12(precision), self._satoh_I6_I10(precision),
                              self._satoh_I6_I12(precision), self._satoh_I10_I12(precision) ])
        
        raise NotImplementedError            

    def grading(self, K) :
        if K is QQ or K in NumberFields() :
            return DegreeGrading([4,6,10,12, 10,14,16,16,18,22])
        
        raise NotImplementedError    
        
    def _generator_names(self, K) :
        if K is QQ or K in NumberFields() :
            return [ "I4", "I6", "I10", "I12", "SB_I4_I6", "SB_I4_I10",
                     "SB_I4_I12", "SB_I6_I10", "SB_I6_I12", "SB_I10_I12" ]
    
    def _generator_by_name(self, K, name) :
        if K is QQ or K in NumberFields() :
            R = self.generator_relations(K).ring()
            try :
                i = self._generator_names(K).index(name)
                return R.gen(i)
            except ValueError:
                raise ValueError("name {0} doesn't exist for {1}".format(name, K))

            #===================================================================
            # if name == "I4" : return R.gen(0)
            # elif name == "I6" : return R.gen(1)
            # elif name == "I10" : return R.gen(2)
            # elif name == "I12" : return R.gen(3)
            # elif name == "SB_I4_I6" : return R.gen(4)
            # elif name == "SB_I4_I10" : return R.gen(5)
            # elif name == "SB_I4_I12" : return R.gen(6)
            # elif name == "SB_I6_I10" : return R.gen(7)
            # elif name == "SB_I6_I12" : return R.gen(8)
            # elif name == "SB_I10_I12" : return R.gen(9)
            #===================================================================
            
        raise NotImplementedError
        
    @cached_method
    def generator_relations(self, K) :
        r"""
        An ideal `I` in a polynomial ring `R`, such that the associated ring
        is `R / I`. This ideal must be unique for `K`. 
        """
        if K is QQ or K in NumberFields() :
            R = PolynomialRing(K, self._generator_names(K))
            ##FIXME: There are relations. Find and implement them.
            return R.ideal(0)
            
        raise NotImplementedError
    
    def weights(self, K) :
        r"""
        A list of integers corresponding to the weights.
        """
        if K is QQ or K in NumberFields() :
            return [10,14,16,16,18,22]
        
        raise NotImplementedError

    def is_vector_valued(self) :
        return True

    def non_vector_valued(self) :
        r"""
        Return the non vector values version of this type. 
        """
        return SiegelModularFormG2_Classical_Gamma()
        
    def vector_valued(self) :
        r"""
        Return the vector values version of this type.
        """
        raise self

    def graded_submodules_are_free(self) :
        return True

    def __cmp__(self, other) :
        return cmp(type(self), type(other))
    
    def __hash__(self) :
        return xor(hash(type(self)), hash(self.grading(QQ)))
