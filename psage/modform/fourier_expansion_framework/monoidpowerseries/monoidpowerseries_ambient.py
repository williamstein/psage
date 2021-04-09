"""
Ambients of monoid power series and ambients of equivariant monoid power series.

AUTHOR :
    - Martin Raum (2010 - 02 - 10) Initial version
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

from past.builtins import cmp
from builtins import object
from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import MonoidPowerSeries, EquivariantMonoidPowerSeries
from sage.rings.all import Integer
from sage.structure.element import Element

#===============================================================================
# MonoidPowerSeriesAmbient_abstract
#===============================================================================

class MonoidPowerSeriesAmbient_abstract(object) :
    r"""
    Given some `K` module or algebra `A` and a monoid `S` filtered over
    a net `\Lambda` construct a module or ring of monoid power series.
    
    Set `R = B[S]`. Then the projective limit of `R / R_\lambda` for
    `\lambda \in \Lambda \rightarrow \infty` considered as a
    `K` module or algebra is implemented by this class.
    """
    
    def __init__(self, A, S) :
        r"""
        INPUT:
            - `A` -- A ring or module.
            - `S` -- A monoid as implemented in :class:~`from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids.NNMonoid`.

        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ambient import MonoidPowerSeriesAmbient_abstract
            sage: mps = MonoidPowerSeriesAmbient_abstract(QQ, NNMonoid(False))
        """
        self.__monoid = S
        
        self._set_multiply_function()
        self.__coefficient_domain = A

        if not hasattr(self, "_element_class") :
            self._element_class = MonoidPowerSeries

    def is_exact(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import MonoidPowerSeriesRing
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: mps.is_exact()
            False
        """
        return False
        
    def monoid(self) :
        r"""
        Return the index monoid of ``self``.
        
        OUTPUT:
            A monoid.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ambient import MonoidPowerSeriesAmbient_abstract
            sage: mps = MonoidPowerSeriesAmbient_abstract(QQ, NNMonoid(False))
            sage: mps.monoid()
            NN
        """
        return self.__monoid
    
    def coefficient_domain(self) :
        r"""
        The coefficient domain of ``self``.
        
        OUTPUT:
            A ring or module.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ambient import MonoidPowerSeriesAmbient_abstract
            sage: mps = MonoidPowerSeriesAmbient_abstract(QQ, NNMonoid(False))
            sage: mps.coefficient_domain()
            Rational Field
            sage: mps = MonoidPowerSeriesAmbient_abstract(FreeModule(ZZ, 3), NNMonoid(False))
            sage: mps.coefficient_domain()
            Ambient free module of rank 3 over the principal ideal domain Integer Ring 
        """
        return self.__coefficient_domain
    
    def _multiply_function(self) :
        r"""
        Return the currect multiply function.
        
        The standard implementation of elements asks its parent to
        provide a multiplication function, which has the following signature:
        
        multiply function INPUT:
          - `k`          -- Element of the index monoid; An index.
          - ``lcoeffs``  -- A dictionary; Coefficient dictionary of the left factor.
          - ``rcoeffs``  -- A dictionary; Coefficient dictionary of the right factor.
          - ``null``     -- An element of a ring or module; An initialized zero object
                            of the coefficient domain.
        multiply function OUTPUT:
          A ring element. The `k`-th coefficent of the product ``left * right``.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import MonoidPowerSeriesRing
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: h = mps._multiply_function()
        """
        return self.__multiply_function
        
    def _set_multiply_function(self, f = None) :
        r"""
        Set the multiply function that is decribed in :meth:~`._multiply_function`.
        If `f` is ``None`` an iteration over the decompositions in the
        monoid is used.
        
        INPUT:
            - `f` -- A function or ``None`` (default: ``None``).
        
        OUTPUT:
            ``None``.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import MonoidPowerSeriesRing
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_element import MonoidPowerSeries
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: e = MonoidPowerSeries( mps, { 4 : 1, 5 : 2}, mps.monoid().filter_all() )
            sage: h = e * e # indirect doctest
            sage: h = lambda : None
            sage: mps._set_multiply_function(h)
            sage: h == mps._multiply_function()
            True
        """
        if not f is None :
            self.__multiply_function = f
            return
        
        def mul(s, lcoeffs, rcoeffs, res) :
            for s1, s2 in self.__monoid.decompositions(s) :
                try :
                    v1 = lcoeffs[s1]
                    v2 = rcoeffs[s2]
                except KeyError :
                    continue
                
                res += v1 * v2
                    
            return res
        #! def mul
        
        self.__multiply_function = mul
        
    def _coerce_map_from_(self, other) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import MonoidPowerSeriesRing
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: mps.has_coerce_map_from( MonoidPowerSeriesRing(ZZ, NNMonoid(False)) ) # indirect doctest
            True
        """
        if isinstance(other, MonoidPowerSeriesAmbient_abstract) :
            if self.monoid() == other.monoid() and \
               self.coefficient_domain().has_coerce_map_from(other.coefficient_domain()) :
                from sage.structure.coerce_maps import CallableConvertMap
                return CallableConvertMap(other, self, self._element_constructor_)
    
    def _element_constructor_(self, x) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import MonoidPowerSeriesRing
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: h = mps(dict()) # indirect doctest
            sage: h = mps(h) # indirect doctest
        """
        if isinstance(x, dict) :
            return self._element_class(self, x, self.monoid().filter_all())
        if isinstance(x, Element) :
            P = x.parent()
            if P is self :
                return x
            elif isinstance(P, MonoidPowerSeriesAmbient_abstract) :
                if self.coefficient_domain() is P.coefficient_domain() :
                    return self._element_class( self,
                            x.coefficients(), x.precision() )
                else :
                    coefficient_domain = self.coefficient_domain()
                    return self._element_class( self,
                            dict( (k,coefficient_domain(c)) for (k,c) in x.coefficients().items() ),
                            x.precision() )
        
        raise TypeError("Cannot construct an element of {0}".format(x))
    
    def __cmp__(self, other) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ambient import MonoidPowerSeriesAmbient_abstract
            sage: mps = MonoidPowerSeriesAmbient_abstract(QQ, NNMonoid(False))
            sage: mps2 = MonoidPowerSeriesAmbient_abstract(ZZ, NNMonoid(False))
            sage: mps == MonoidPowerSeriesAmbient_abstract(QQ, NNMonoid(False))
            True
            sage: mps == mps2
            False
        """
        c = cmp(type(self), type(other))
        
        if c == 0 :
            c = cmp(self.coefficient_domain(), other.coefficient_domain())
        if c == 0 :
            c = cmp(self.monoid(), other.monoid())
            
        return c

###############################################################################
###############################################################################
###############################################################################

#===============================================================================
# EquivariantMonoidPowerSeriesAmbient_abstract
#===============================================================================

class EquivariantMonoidPowerSeriesAmbient_abstract(object) :
    r"""
    Given some ring or module `A`, a monoid `S` filtered over some originated
    net `\Lambda` such that all induced submonoids are finite, a group `G`, a
    semigroup `C` with a map `c \rightarrow \mathrm{Hom}(G, Aut_K(A))`, a
    homomorphism `\phi : G -> Aut(S)` and a homomorphism `\eta : G -> C`, where
    `K` is the base ring of `A`.
    
    Suppose for every `c, c'` in `C`, and `g` in `G`, and `a, a'` in `A` we have
      `(c c') (g) (a a') = c(g)(a) c'(g)(a')`. 
    Set `R = B[C][S]`. Then the projective limit of
      `R / R_\lambda` for `\lambda \in \Lambda \rightarrow \infinity` is a
      `K`-algebra or -module.
    """
        
    def __init__(self, O, C, R) :
        r"""
        INPUT:
            - `O` -- A monoid with an action of a group; As implemented in
                     :class:~`fourier_expansion_framework.monoidpowerseries.NNMonoid`.
            - `C` -- A monoid of characters; As implemented in ::class:~`fourier_expansion_framework.monoidpowerseries.CharacterMonoid_class`.
            - `R` -- A representation on an algebra or module; As implemented
                     in :class:~`fourier_expansion_framework.monoidpowerseries.TrivialRepresentation`.
        
        EXAMPLES::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import EquivariantMonoidPowerSeriesRing
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ)) # indirect doctest
        """
        self.__action = O
        self.__characters = C
        self.__representation = R

        self.__coefficient_domain = R.codomain()
        self.__monoid = O.monoid()
        
        self._set_reduction_function()
        self._set_character_eval_function()
        self._set_apply_function()
        self._set_multiply_function()
        
        if not hasattr(self, "_element_class") :
            self._element_class = EquivariantMonoidPowerSeries
    
    def is_exact(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import EquivariantMonoidPowerSeriesRing
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
            sage: emps.is_exact()
            False
        """
        return False
    
    def coefficient_domain(self) :
        r"""
        The domain of coefficients.
        
        OUTPUT:
            Either a ring or a module.
        
        EXAMPLES::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import EquivariantMonoidPowerSeriesRing
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
            sage: emps.coefficient_domain() == QQ
            True
        """
        return self.__coefficient_domain
    
    def group(self) :
        r"""
        The group acting on the index monoid.
        
        OUTPUT:
            Of arbitrary type.
        
        NOTE:
            The framework might change at a later time such that this
            function returns a group.
        
        EXAMPLES::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import EquivariantMonoidPowerSeriesRing
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
            sage: emps.group()
            '1'
        """
        return self.__action.group()
        
    def monoid(self) :
        r"""
        The index monoid.
        
        OUTPUT:
            A monoid as implemented in :class:~`fourier_expansion_framework.monoidpowerseries.NNMonoid`.
        
        EXAMPLES::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import EquivariantMonoidPowerSeriesRing
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
            sage: emps.monoid()
            NN
        """
        return self.__action.monoid()
    
    def action(self) :
        r"""
        The index monoid with the action of a group.
        
        OUTPUT:
            A monoid with action as implemented in :class:~`fourier_expansion_framework.monoidpowerseries.NNMonoid`.
        
        EXAMPLES::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import EquivariantMonoidPowerSeriesRing
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
            sage: emps.action()
            NN with action
        """
        return self.__action
    
    def characters(self) :
        r"""
        The monoid of characters associated to the monoid index' group.
        
        OUTPUT:
            A monoid as implemented in :class:~`fourier_expansion_framework.monoidpowerseries.CharacterMonoid_class`.
        
        EXAMPLES::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import EquivariantMonoidPowerSeriesRing
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
            sage: emps.characters()
            Character monoid over Trivial monoid
        """
        return self.__characters
        
    def representation(self) :
        r"""
        The representation on the coefficient domain.
        
        OUTPUT:
            A representation as implemented in :class:~`fourier_expansion_framework.monoidpowerseries.TrivialRepresentation`.

        EXAMPLES::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import EquivariantMonoidPowerSeriesRing
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
            sage: emps.representation()
            Trivial representation of 1 on Rational Field
        """
        return self.__representation
    
    def _reduction_function(self) :
        r"""
        The reduction function accepts an index `s`. It returns the pair
        reduction `(rs = g^-1 s, g)` of `s` with a group element `g`.
        
        OUTPUT:
            A function.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import EquivariantMonoidPowerSeriesRing
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
            sage: red_function = emps._reduction_function()
            sage: red_function(2)
            (2, 1)
        """
        return self.__reduction_function
    
    def _set_reduction_function(self, f = None) :
        r"""
        Set the reduction function, explained in :class:~`._reduction_function`.
        If `f` is ``None`` the reduction function of the action is used.
        
        INPUT:
            - `f` -- A function or ``None`` (default: ``None``).
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import EquivariantMonoidPowerSeriesRing
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
            sage: h = lambda : None
            sage: emps._set_reduction_function(h)
            sage: h == emps._reduction_function()
            True
            sage: emps._set_reduction_function() ## This is important since emps is globally unique
        """
        if not f is None :
            self.__reduction_function = f
            return
        
        self.__reduction_function = self.__action._reduction_function()
    
    def _character_eval_function(self) :
        r"""
        The character evaluation function. It accepts a character `c` and 
        a group element `g` and returns `c(g)`.
        
        OUTPUT:
            A function.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import EquivariantMonoidPowerSeriesRing
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
            sage: eval_func = emps._character_eval_function()
            sage: eval_func(emps.characters().one_element(), 1)
            1 
        """
        return self.__character_eval_function 
    
    def _set_character_eval_function(self, f = None) :
        r"""
        Set the character evaluation function. If `f` is ``None``, the 
        implementation of the character monoid is used.
        
        INPUT:
            - `f` -- A function or ``None`` (default: ``None``).
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import EquivariantMonoidPowerSeriesRing
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
            sage: h = lambda c,g: 1
            sage: emps._set_character_eval_function(h)
            sage: h == emps._character_eval_function()
            True
            sage: emps._set_character_eval_function() ## This is important since emps is globally unique
        """
        if not f is None :
            self.__character_eval_function = f
            return
        
        self.__character_eval_function = self.__characters._eval_function() 
    
    def _apply_function(self) :
        r"""
        The apply function. It applies a group element `g` to an element `v`
        of the coefficient domain, the base space of the representation.
        
        OUTPUT:
            A function.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import EquivariantMonoidPowerSeriesRing
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
            sage: app_func = emps._apply_function()
            sage: app_func(1, 1/2)
            1/2
        """
        return self.__apply_function
        
    def _set_apply_function(self, f = None) :
        r"""
        Set the apply function. If `f` is ``None``, the implementation of
        the representation is used.
        
        INPUT:
            - `f` -- A function or ``None`` (default: ``None``).
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import EquivariantMonoidPowerSeriesRing
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
            sage: h = lambda : None
            sage: emps._set_apply_function(h)
            sage: h == emps._apply_function()
            True
            sage: emps._set_apply_function() ## This is important since emps is globally unique
        """
        if not f is None :
            self.__apply_function = f
            return
        
        self.__apply_function = self.__representation._apply_function()
        
    def _multiply_function(self) :
        r"""
        The standard implementation of elements of this ring will ask its
        parent to provide multplication function, which has the
        following signature:
        
        multiply function INPUT :
          - `k`          -- An index.
          - ``lcoeffs``  -- A dictionary; The coefficient dictionary of the left factor.
          - ``rcoeffs``  -- A dictionary; The coefficient dictionary of the right factor.
          - ``lch``      -- A character; The character of the left factor.
          - ``rch``      -- A character; The character of the right factor.
          - ``null``     -- A ring or module element; TA initialized zero object of
                            the coefficient domain.
        multiply function OUTPUT :
          A ring or module element. The `k`-th coefficent of the product
          ``left * right``, where the power series transform with character
          ``lch`` and ``rch`` respectively under the action of the power series'
          symmetry group.
          
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import EquivariantMonoidPowerSeriesRing
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
            sage: h = emps._multiply_function()
            sage: e = emps( {emps.characters().one_element() : {4 : 3, 5 : 3} } )
            sage: (e * e).coefficients() # indirect doctest
            {8: 9, 9: 18, 10: 9}
        """
        return self.__multiply_function
    
    def _set_multiply_function(self, f = None) :
        r"""
        Set the multiply function. If `f` is ``None`` an iteration over the
        decompositions in the monoid is used.
        
        INPUT:
            - `f` -- A function or ``None`` (default: ``None``).

        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import EquivariantMonoidPowerSeriesRing
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
            sage: h = lambda : None
            sage: emps._set_multiply_function(h)
            sage: h == emps._multiply_function()
            True
        """
        if not f is None :
            self.__multiply_function = f
            return
        
        def mul(s, lcoeffs, rcoeffs, cl, cr, res) :
            reduction = self._reduction_function()
            apply = self._apply_function()
            character_ev = self._character_eval_function()
            
            rs, g = reduction(s)
            
            for s1, s2 in self.__monoid.decompositions(rs) :
                rs1, g1 = reduction(s1)
                rs2, g2 = reduction(s2)
                
                try :
                    v1 = lcoeffs[rs1]
                    v2 = rcoeffs[rs2]
                except KeyError :
                    continue
                v1 = g1(*v1)
                v2 = g2(*v2)
    
                res += (character_ev(g1, cl) * character_ev(g2, cr)) * v1 * v2
                
            return character_ev(g, cl*cr) * g(*res)
        #! def mul
        
        self.__multiply_function = mul

    def _coerce_map_from_(self, other) :
        r"""
        TODO:
          This is a stub. The dream is that representations know about
          compatible coercions and so would actions and characters. Then
          every equivariant monoid power series ring would be a functorial
          construction in all three parameters (The functor would then be
          applied to a representation within a universe of representations).

        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import EquivariantMonoidPowerSeriesRing
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
            sage: emps.has_coerce_map_from( EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", ZZ)) ) # indirect doctest
            True
        """
        if isinstance(other, EquivariantMonoidPowerSeriesAmbient_abstract) :
            if self.action() == other.action() and \
               self.characters() == other.characters() :
                if self.representation().extends(other.representation()) :
                    from sage.structure.coerce_maps import CallableConvertMap
                    return CallableConvertMap(other, self, self._element_constructor_)

    def _element_constructor_(self, x) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring import MonoidPowerSeriesRing, EquivariantMonoidPowerSeriesRing
            sage: emps = EquivariantMonoidPowerSeriesRing(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
            sage: mps = MonoidPowerSeriesRing(QQ, NNMonoid(False))
            sage: h = emps(dict()) # indirect doctest
            sage: h = emps(1)
            sage: h = emps(h)
            sage: h = emps(mps(1))
        """
        if isinstance(x, Element) :
            P = x.parent()
            if P is self :
                return x
            elif isinstance(P, EquivariantMonoidPowerSeriesAmbient_abstract) :
                if self.coefficient_domain() is P.coefficient_domain() :
                    return self._element_class( self,
                                                x.coefficients(True), x.precision() )
                else :
                    coefficient_domain = self.coefficient_domain()
                    
                    return self._element_class( self,
                                                dict( (ch, dict( (k,coefficient_domain(c)) for (k,c) in coeffs.items()) )
                                                       for (ch, coeffs) in x.coefficients(True).items() ),
                                                x.precision() )
                    
            elif isinstance(P, MonoidPowerSeriesAmbient_abstract) :
                if self.coefficient_domain() is P.coefficient_domain() :
                    return self._element_class(
                            self, dict([(self.__characters.one_element(),
                                         x.coefficients())]),
                            self.action().filter(x.precision()),
                            symmetrise = True  )
                else :
                    return self._element_class(
                            self, dict([(self.__characters.one_element(),
                                         dict( (k,self.coefficient_domain(c)) for (k,c) in x.coefficients().items()) )]),
                            self.action().filter(x.precision()),
                            symmetrise = True  )
        elif isinstance(x, dict) :
            if len(x) != 0 :
                try :
                    if list(x.keys())[0].parent() is self.__characters :
                        return self._element_class( self,
                                x, self.action().filter_all() )
                except AttributeError :
                    pass
            
            return self._element_class( self,
                    dict([(self.__characters.one_element(), x)]),
                    self.action().filter_all() )
        
        raise TypeError("can't convert {0} into {1}".format(x, self))
        
    def __cmp__(self, other) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import *
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ambient import EquivariantMonoidPowerSeriesAmbient_abstract
            sage: emps = EquivariantMonoidPowerSeriesAmbient_abstract(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
            sage: emps == EquivariantMonoidPowerSeriesAmbient_abstract(NNMonoid(True), TrivialCharacterMonoid("1", QQ), TrivialRepresentation("1", QQ))
            True
            sage: emps2 = EquivariantMonoidPowerSeriesAmbient_abstract(NNMonoid(True), TrivialCharacterMonoid("1", ZZ), TrivialRepresentation("1", ZZ))
            sage: emps == emps2
            False
        """
        c = cmp(type(self), type(other))
        
        if c == 0 :
            c = cmp(self.action(), other.action())
        if c == 0 :
            c = cmp(self.representation(), other.representation())
        if c == 0 :
            c = cmp(self.characters(), other.characters())
            
        return c
