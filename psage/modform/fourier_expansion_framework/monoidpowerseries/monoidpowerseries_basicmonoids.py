r"""
Implementation of base classes used as parameters for a ring of monoid 
power series.

AUTHOR :
    - Martin Raum (2009 - 07 - 25) Initial version
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

from operator import xor
from sage.misc.latex import latex
from sage.categories.all import Rings
from sage.rings.infinity import infinity
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.structure.sage_object import SageObject

#===============================================================================
# CharacterMonoid_class
#===============================================================================

class CharacterMonoid_class ( SageObject ) :
    r"""
    The characters for an equivariant monoid power series must
    form a monoid. A basic implementation is given here.
    """
    
    def __init__(self, G, C, codomain, eval, C_multiplicative = True) :
        r"""
        INPUT:
        
            - `G`          -- Every type accepted; Representative of a group which
                              the characters can be evaluated on.
            - `C`          -- A class implementing all functions of :class:`~.TrivialMonoid_class`;
                              A monoid whose elements represent one character each.
            - ``codomain`` -- A ring; A common codomain for all characters.
            - ``eval``     -- A function; It accepts a group element and a character,
                              returning the evaluation of this character at the 
                              group element.
            - ``C_multiplicative`` -- A boolean; If true the monoid `C` is assumed
                              to be multiplicative. Otherwise it is assumed to be additive.
        
        NOTE:
            The interface may change in the future, enforcing `G` to be a group.
        
        EXAMPLES::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import CharacterMonoid_class, TrivialMonoid
            sage: cm = CharacterMonoid_class("ZZ", TrivialMonoid(), QQ, lambda g, c: 1)
            
        TESTS::
            sage: cm = CharacterMonoid_class(None, TrivialMonoid(), AbelianGroup([3,3]), lambda g, c: 1)
            sage: cm = CharacterMonoid_class("ZZ", ZZ, AbelianGroup([3,3]), lambda g, c: 1, True)
        """
        self.__G = G
        self.__C = C
        self.__codomain = codomain
        self.__eval = eval
        self.__C_multiplicative = C_multiplicative

    def ngens(self) :
        r"""
        OUTPUT:
            An integer.
            
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import CharacterMonoid_class, TrivialMonoid
            sage: cm = CharacterMonoid_class("ZZ", TrivialMonoid(), QQ, lambda g, c: 1)
            sage: cm.ngens()
            1
            sage: cm = CharacterMonoid_class("ZZ", AbelianGroup([3, 3]), ZZ, lambda g, c: 1)
            sage: cm.ngens()
            2
        """
        return self.__C.ngens()
    
    def gen(self, i = 0) :
        r"""
        OUTPUT:
            An instance of :class:`~.CharacterMonoidElement_class`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import CharacterMonoid_class, TrivialMonoid
            sage: cm = CharacterMonoid_class("ZZ", TrivialMonoid(), QQ, lambda g, c: 1)
            sage: cm.gen()
            1
        """
        return CharacterMonoidElement_class(self, self.__C.gen(i))
    
    def gens(self) :
        r"""
        OUTPUT:
            A tuple of instances of :class:`~.CharacterMonoidElement_class`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import CharacterMonoid_class, TrivialMonoid
            sage: cm = CharacterMonoid_class("ZZ", TrivialMonoid(), QQ, lambda g, c: 1)
            sage: cm.gens()
            (1,)
        """
        return tuple(map(lambda c: CharacterMonoidElement_class(self, c), self.__C.gens()))
    
    def group(self) :
        r"""
        Return the group, which is the common domain of all characters
        of this monoid of characters.
        
        OUTPUT:
            Of arbitrary type.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import CharacterMonoid_class, TrivialMonoid
            sage: cm = CharacterMonoid_class("ZZ", TrivialMonoid(), QQ, lambda g, c: 1)
            sage: cm.group()
            'ZZ'
        """
        return self.__G

    def codomain(self) :
        r"""
        Return the common codomain of all characters of this monoid of characters.

        OUTPUT:
            A ring.
        
        EXAMPLES::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import CharacterMonoid_class, TrivialMonoid
            sage: cm = CharacterMonoid_class("ZZ", TrivialMonoid(), QQ, lambda g, c: 1)
            sage: cm.codomain()
            Rational Field
        """
        return self.__codomain
  
    def monoid(self) :
        r"""
        Return the abstract monoid underlying this monoid of characters.
        
        OUTPUT:
            A monoid.
        
        EXAMPLES::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import CharacterMonoid_class, TrivialMonoid
            sage: cm = CharacterMonoid_class("ZZ", TrivialMonoid(), QQ, lambda g, c: 1)
            sage: cm.monoid()
            Trivial monoid
            sage: cm = CharacterMonoid_class("ZZ", AbelianGroup([3, 3]), ZZ, lambda g, c: 1)
            sage: cm.monoid()
            Multiplicative Abelian Group isomorphic to C3 x C3
        """
        return self.__C
    
    def extends(self, other) :
        r"""
        Decide whether ``self`` extends ``other``. Namely, whether there is a is an embedding of ``other``
        into ``self`` that is compatible with a common codomain. A negative answer does not mean
        that there is no such embedding.
        
        INPUT:
            - other    -- An instance of :class:.`~CharacterMonoid_class`.
        
        OUTPUT:
            A boolean.
        
        EXAMPLES::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import CharacterMonoid_class, TrivialMonoid
            sage: cm1 = CharacterMonoid_class("ZZ", TrivialMonoid(), QQ, lambda g, c: 1)
            sage: cm2 = CharacterMonoid_class("ZZ", TrivialMonoid(), QQ, lambda g, c: 1)
            sage: cm1 == cm2
            False 
        """
        return self == other
    
    def _is_C_multiplicative(self) :
        r"""
        Return whether the underlying monoid is multiplicative.
        
        OUTPUT:
            A boolean.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import CharacterMonoid_class, TrivialMonoid
            sage: cm = CharacterMonoid_class("ZZ", TrivialMonoid(), QQ, lambda g, c: 1)
            sage: cm._is_C_multiplicative()
            True
        """
        return self.__C_multiplicative
    
    def _eval_function(self) :
        r"""
        Return the evaluation function, mapping an element of the associated group and an element of ``self``
        to an element of the codomain.
        
        OUTPUT:
            A function with signature `g, c \mapsto e`. 
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import CharacterMonoid_class, TrivialMonoid
            sage: e = lambda g, c: 1
            sage: cm = CharacterMonoid_class("ZZ", TrivialMonoid(), QQ, e)
            sage: e == cm._eval_function()
            True
        """
        return self.__eval
        
    def _apply(self, g, c, b) :
        r"""
        Apply `c(g)` to some element `b` by multiplication.
        
        INPUT:
            - g    -- Arbitrary type; A representative of the group.
            - c    -- An element of self; A character.
            - b    -- An element of a ring with embedding from the codomain. An element to
                      apply `c(g)` to.
        
        OUTPUT:
           An element of a ring. 
        
        EXAMPLES::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import CharacterMonoid_class, TrivialMonoid
            sage: cm = CharacterMonoid_class("ZZ", TrivialMonoid(), QQ, lambda g, c: 1)
            sage: cm._apply(1, 1, 2)
            2
        """
        return self.__eval(g, c) * b

    def __cmp__(self, other) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import CharacterMonoid_class, TrivialMonoid
            sage: cm1 = CharacterMonoid_class("ZZ", TrivialMonoid(), QQ, lambda g, c: 1)
            sage: cm2 = CharacterMonoid_class("ZZ", TrivialMonoid(), QQ, lambda g, c: 1)
            sage: cm1 == cm2
            False
        """
        c = cmp(type(self), type(other))
        if c == 0 :
            c = cmp(self.__G, other.__G)
        if c == 0 :
            c = cmp(self.__codomain, other.__codomain)
        if c == 0 :
            c = cmp(self.__C, other.__C)
        if c == 0 and not self.__eval is other.__eval :
            return -1
        
        return c
    
    def __iter__(self) :
        r"""
        TEST::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import CharacterMonoid_class, TrivialMonoid
            sage: cm = CharacterMonoid_class("ZZ", TrivialMonoid(), QQ, lambda g, c: 1)
            sage: list(cm)
            [1]
        """
        for c in self.__C :
            yield CharacterMonoidElement_class(self, c)
            
        raise StopIteration
    
    def one_element(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import CharacterMonoid_class, TrivialMonoid
            sage: cm = CharacterMonoid_class("ZZ", TrivialMonoid(), QQ, lambda g, c: 1)
            sage: cm.one_element()
            1
        """
        return CharacterMonoidElement_class(self, self.__C.one_element())

    def __call__(self, x) :
        r"""
        Convert an element to ``self``.
        
        INPUT:
            - `x` -- An element that is convertible to :meth:~`.monoid()`.
        
        OUTPUT:
            An element of ``self``.
        
        TESTS::
           sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import CharacterMonoid_class, TrivialMonoid
            sage: cm = CharacterMonoid_class("ZZ", TrivialMonoid(), QQ, lambda g, c: 1)
            sage: cm(1)
            1 
        """
        return CharacterMonoidElement_class(self, self.__C(x))
    
    def __hash__(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import CharacterMonoid_class, TrivialMonoid
            sage: cm = CharacterMonoid_class("ZZ", TrivialMonoid(), QQ, lambda g, c: 1)
            sage: d = dict([(cm.one_element, 0)]) # indirect doctest
        """
        return xor(hash(self.__C), hash(self.__G))
    
    def _repr_(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import CharacterMonoid_class, TrivialMonoid
            sage: CharacterMonoid_class("ZZ", TrivialMonoid(), QQ, lambda g, c: 1)
            Character monoid over Trivial monoid
        """
        return "Character monoid over %s" % self.__C
    
    def _latex_(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import CharacterMonoid_class, TrivialMonoid
            sage: latex(CharacterMonoid_class("ZZ", TrivialMonoid(), QQ, lambda g, c: 1))
            \text{Character monoid over } \text{Trivial monoid}
        """
        return r"\text{Character monoid over }" + latex(self.__C)

#===============================================================================
# CharacterMonoidElement_class
#===============================================================================

class CharacterMonoidElement_class ( SageObject ) :
    """
    The element class of :class:`~.CharacterMonoid_class`.
    """
    
    def __init__(self, parent, c) :
        r"""
        INPUT:
            - ``parent``  -- An instance of :class:`~.CharacterMonoid_class`; The parent of ``self``.
            - `c`         -- An element of a monoid; The element in the underlying monoid of ``parent``
                             associated to ``self``.
                
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import CharacterMonoid_class, CharacterMonoidElement_class, TrivialMonoid
            sage: cm = CharacterMonoid_class("ZZ", TrivialMonoid(), QQ, lambda g, c: 1)
            sage: c = CharacterMonoidElement_class(cm, cm.monoid().one_element())
        """
        self.__parent = parent
        self.__c = c
        
    def parent(self) :
        r"""
        OUTPUT:
            An instance of :class:`~.CharacterMonoid_class`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import CharacterMonoid_class, CharacterMonoidElement_class, TrivialMonoid
            sage: cm = CharacterMonoid_class("ZZ", TrivialMonoid(), QQ, lambda g, c: 1)
            sage: c = CharacterMonoidElement_class(cm, cm.monoid().one_element())
            sage: c.parent() is cm
            True
        """
        return self.__parent
    
    def _monoid_element(self) :
        r"""
        The underlying element of the monoid.
        
        OUTPUT:
            An element of a character.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import CharacterMonoid_class, CharacterMonoidElement_class, TrivialMonoid
            sage: cm = CharacterMonoid_class("ZZ", TrivialMonoid(), QQ, lambda g, c: 1)
            sage: c = CharacterMonoidElement_class(cm, cm.monoid().one_element())
            sage: c._monoid_element() == cm.monoid().one_element()
            True
        """
        return self.__c
            
    def __call__(self, g) :
        r"""
        OUTPUT:
            An element of the codomain of the parent.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import CharacterMonoid_class, CharacterMonoidElement_class, TrivialMonoid
            sage: cm = CharacterMonoid_class("ZZ", TrivialMonoid(), QQ, lambda g, c: 1)
            sage: c = CharacterMonoidElement_class(cm, cm.monoid().one_element())
            sage: c(2)
            1
            
        Test old bug, where elements of the underlying monoid were passed to the eval function::
            sage: cm = CharacterMonoid_class("ZZ", TrivialMonoid(), QQ, lambda g, c: 1 if isinstance(c, CharacterMonoidElement_class) else -1)
            sage: c = CharacterMonoidElement_class(cm, cm.monoid().one_element())
            sage: c(2)
            1
        """
        return self.parent()._eval_function()(g, self)

    def __mul__(left, right) :
        r"""
        NOTE:
            If the underlying monoid is additive the character are nevertheless multiplied.
        
        OUTPUT:
            An instance of :class:`~.CharacterMonoidElement_class`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import CharacterMonoid_class, CharacterMonoidElement_class, TrivialMonoid
            sage: cm = CharacterMonoid_class("ZZ", TrivialMonoid(), QQ, lambda g, c: 1)
            sage: c = CharacterMonoidElement_class(cm, cm.monoid().one_element())
            sage: c * c     # indirect doctest
            1
            sage: cm = CharacterMonoid_class("ZZ", ZZ, QQ, lambda g, c: 1, False)
            sage: c = CharacterMonoidElement_class(cm, cm.monoid().one_element())
            sage: c * c # indirect doctest
            2
        """
        if not isinstance(right, CharacterMonoidElement_class) or \
           not left.parent() == right.parent() :
            raise TypeError( "Right factor should be an element of the monoid." )
        
        if left.parent()._is_C_multiplicative() :
            return CharacterMonoidElement_class(left.parent(), left.__c * right.__c)
        else :
            return CharacterMonoidElement_class(left.parent(), left.__c + right.__c)
    
    def __cmp__(self, other) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import CharacterMonoid_class, CharacterMonoidElement_class, TrivialMonoid
            sage: cm = CharacterMonoid_class("ZZ", ZZ, QQ, lambda g, c: 1, False)
            sage: c1 = CharacterMonoidElement_class(cm, cm.monoid().one_element())
            sage: c2 = CharacterMonoidElement_class(cm, cm.monoid().one_element())
            sage: c1 == c2    # indirect doctest
            True
            sage: c1 * c1 == c1    # indirect doctest
            False
        """
        c = cmp(type(self), type(other))
        if c == 0 :
            c = cmp(self.__parent, other.__parent)
        if c == 0 :
            c = cmp(self.__c, other.__c)
            
        return c
    
    def __hash__(self) :
        r"""
        OUTPUT:
            An integer.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import CharacterMonoid_class, CharacterMonoidElement_class, TrivialMonoid
            sage: cm = CharacterMonoid_class("ZZ", TrivialMonoid(), QQ, lambda g, c: 1)
            sage: c1 = CharacterMonoidElement_class(cm, cm.monoid().one_element())
            sage: hash(c1)    # indirect doctest
            -582796950           # 32-bit
            -3589969894075844246 # 64-bit
        """
        return xor(hash(self.__parent), hash(self.__c))
    
    def _repr_(self) :
        r"""
        OUTPUT:
            A string.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import CharacterMonoid_class, CharacterMonoidElement_class, TrivialMonoid
            sage: cm = CharacterMonoid_class("ZZ", TrivialMonoid(), QQ, lambda g, c: 1)
            sage: c1 = CharacterMonoidElement_class(cm, cm.monoid().one_element())
            sage: repr(c1)    # indirect doctest
            '1'
        """
        return repr(self.__c)
    
    def _latex_(self) :
        r"""
        OUTPUT:
            A string.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import CharacterMonoid_class, CharacterMonoidElement_class, TrivialMonoid
            sage: cm = CharacterMonoid_class("ZZ", TrivialMonoid(), QQ, lambda g, c: 1)
            sage: c1 = CharacterMonoidElement_class(cm, cm.monoid().one_element())
            sage: latex(c1)    # indirect doctest
            1
        """
        return latex(self.__c)

#===============================================================================
# TrivialMonoid
#===============================================================================

class TrivialMonoid ( SageObject ) :
    r"""
    A monoid with one element, which implements only functions that are necessary for
    being an arguemnt of :meth:`~.CharacterMonoid_class.__init__`.
    """
    
    def __init__(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import TrivialMonoid
            sage: m = TrivialMonoid()
        """
        SageObject.__init__(self)
    
    def ngens(self) :
        r"""
        OUTPUT:
            An integer.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import TrivialMonoid
            sage: m = TrivialMonoid()
            sage: m.ngens()
            1
        """
        return 1
    
    def gen(self, i = 0) :
        r"""
        OUTPUT:
            An integer.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import TrivialMonoid
            sage: m = TrivialMonoid()
            sage: m.gen()
            1
            sage: m.gen(1)
            Traceback (most recent call last):
            ...
            ValueError: Generator not defined.
        """
        if i == 0 :
            return Integer(1)
        
        raise ValueError("Generator not defined.")
    
    def gens(self) :
        r"""
        OUTPUT:
            A tuple of integers.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import TrivialMonoid
            sage: m = TrivialMonoid()
            sage: m.gens()
            (1,)
        """
        return (Integer(1), )
    
    def one_element(self) :
        r"""
        OUTPUT:
            An integer.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import TrivialMonoid
            sage: m = TrivialMonoid()
            sage: m.one_element()
            1
        """
        return Integer(1)
    
    def __call__(self, x) :
        r"""
        INPUT:
            - `x` -- An integer.
        
        OUTPUT:
            An integer.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import TrivialMonoid
            sage: m = TrivialMonoid()
            sage: m(1)
            1
            sage: m(2)
            Traceback (most recent call last):
            ...
            TypeError: Cannot convert 2 into Trivial monoid.
        """
        if x == 1 :
            return 1
        
        raise TypeError( "Cannot convert %s into Trivial monoid." % (x,) )
    
    def __cmp__(self, other) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import TrivialMonoid
            sage: TrivialMonoid() == TrivialMonoid()    # indirect doctest
            True
        """
        return cmp(type(self), type(other))
    
    def __iter__(self) :
        r"""
        OUTPUT:
            A generator over integers.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import TrivialMonoid
            sage: m = TrivialMonoid()
            sage: list(m)    # indirect doctest
            [1]
        """
        yield Integer(1)
        
        raise StopIteration
    
    def _repr_(self) :
        r"""
        OUTPUT:
            A string.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import TrivialMonoid
            sage: m = TrivialMonoid()
            sage: repr(m)    # indirect doctest
            'Trivial monoid'
        """
        return "Trivial monoid"
    
    def _latex_(self) :
        r"""
        OUTPUT:
            A string.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import TrivialMonoid
            sage: m = TrivialMonoid()
            sage: latex(m)    # indirect doctest
            \text{Trivial monoid}
        """
        return r"\text{Trivial monoid}"

#===============================================================================
# TrivialCharacterMonoid
#===============================================================================

_trivial_evaluations = dict()

def TrivialCharacterMonoid(G, K) :
    r"""
    Return the monoid of characters with codomain `K` with one element.
    
    INPUT:
        - `K`    -- A ring; The codomain for the characters.
        
    OUTPUT:
        An instance of :class:`~.CharacterMonoid_class`.
    
    TESTS::
        sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import TrivialCharacterMonoid
        sage: h = TrivialCharacterMonoid("ZZ", QQ)
        sage: h
        Character monoid over Trivial monoid
        sage: h == TrivialCharacterMonoid("ZZ", QQ)
        True
    """ 
    global _trivial_evaluations
    
    try :
        eval = _trivial_evaluations[K]
    except KeyError :
        eval = lambda g, c : K.one_element()
        _trivial_evaluations[K] = eval
    
    return CharacterMonoid_class(G, TrivialMonoid(), K, eval)

#===============================================================================
# TrivialRepresentation
#===============================================================================

class TrivialRepresentation ( SageObject ) :
    r"""
    A trivial representation over a group `G` with codomain `K` occurring as a representation for
    :class:`~fourier_expansion_framework.monoidpowerseries.EquivariantMonoidPowerSeriesAmbient_abstract`.
    """
    
    def __init__(self, G, K) :
        r"""
        INPUT:
            - `G`    -- Arbitrary type; A group or representative of a group.
            - `K`    -- A module or ring; The codomain of the representation.
        
        OUTPUT:
            An instance of :class:`~.TrivialRepresentation`.
            
        NOTE:
            The interface may change late, enforcing `G` to be an actual group.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import TrivialRepresentation
            sage: rep = TrivialRepresentation("ZZ", QQ)    # indirect doctest
        """
        self.__G = G
        self.__K = K
        
    def base_ring(self) :
        r"""
        The base ring of the representation, commuting with the action.
        
        OUTPUT:
            A ring.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import TrivialRepresentation
            sage: rep = TrivialRepresentation("ZZ", QQ)
            sage: rep.base_ring()
            Rational Field
            sage: rep = TrivialRepresentation("ZZ", FreeModule(ZZ, 3))
            sage: rep.base_ring()
            Integer Ring
        """
        if self.__K in Rings() :
            return self.__K
        else :
            return self.__K.base_ring()
    
    def codomain(self) :
        r"""
        The codomain of the representation.
        
        OUTPUT:
            A ring or a module.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import TrivialRepresentation
            sage: rep = TrivialRepresentation("ZZ", FreeModule(ZZ, 3))
            sage: rep.codomain()
            Ambient free module of rank 3 over the principal ideal domain Integer Ring
        """
        return self.__K
    
    def from_module(self, R) :
        r"""
        Construct a representation of the same type with codomain `R`.
        
        INPUT:
        
        - `R`        -- A module that the representation can be restricted or
                        extended to.
        
        OUTPUT:
        
        - An instance of :class:`~.TrivialRepresentation`.
        
        EXAMPLES::
        
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import TrivialRepresentation
            sage: rep = TrivialRepresentation("ZZ", FreeModule(ZZ, 3))
            sage: rep2 = rep.from_module(FreeModule(QQ, 2))
            sage: rep2.codomain()
            Vector space of dimension 2 over Rational Field
        """
        return TrivialRepresentation( self.__G, R )
    
    def base_extend(self, L) :
        r"""
        Extend the representation's codomain by `L`.
        
        INPUT:
            - `L`    -- A ring; A new base ring.
        
        OUTPUT:
            An instance of :class:`~.TrivialRepresentation`.
        
        EXAMPLES::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import TrivialRepresentation
            sage: rep = TrivialRepresentation("ZZ", FreeModule(ZZ, 3))
            sage: repQQ = rep.base_extend(QQ)
            sage: repQQ.codomain()
            Vector space of dimension 3 over Rational Field
            sage: repQQ.base_extend(GF(3))
            Traceback (most recent call last):
            ...
            TypeError: Base extension of self (over 'Rational Field') to ring 'Finite Field of size 3' not defined.
            sage: rep = TrivialRepresentation("ZZ", ZZ)
            sage: repQQ = rep.base_extend(QQ)
            sage: repQQ.codomain()
            Rational Field
            sage: repQQ.base_extend(GF(3))
            Traceback (most recent call last):
            ...
            AssertionError
        """
        if self.__K in Rings() :
            assert L.has_coerce_map_from(self.__K)
            return TrivialRepresentation( self.__G, L )
        else :
            return TrivialRepresentation( self.__G, self.__K.base_extend(L) )
        
    def extends(self, other) :
        r"""
        Wheter ``self`` is an extension of ``other``.
        
        INPUT:
            - ``other``    -- A representation.
        
        OUTPUT:
            A boolean.
        
        NOTE:
            This does not necessarily mean that ``self`` ocurres as a base extension of
            ``other``.
        
        EXAMPLES::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import TrivialRepresentation
            sage: rep1 = TrivialRepresentation("ZZ", ZZ)
            sage: rep2 = TrivialRepresentation("ZZ", QQ)
            sage: rep1.extends(rep2)
            False
            sage: rep2.extends(rep1)
            True
        """
        if type(self) != type(other) :
            return False
        
        return self.__K.has_coerce_map_from(other.__K)

    def group(self) :
        r"""
        The group that acts.
        
        OUTPUT:
            Arbitrary type.
        
        EXAMPLES::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import TrivialRepresentation
            sage: rep = TrivialRepresentation("ZZ", QQ)
            sage: rep.group()
            'ZZ'
        """
        return self.__G
    
    def _apply_function(self) :
        r"""
        A function that maps a group element and an element of the codomain to its image.
        
        OUTPUT:
            A function with signature `(g,a) \mapsto g.a`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import TrivialRepresentation
            sage: rep = TrivialRepresentation("ZZ", QQ)
            sage: rep._apply_function() == rep.apply
            True
        """
        return self.apply
    
    def apply(self, g, a) :
        r"""
        Return the image of `a` under the transformation `g`.
        
        INPUT:
            - `g`    -- Arbitrary type; A group element.
            - `a`    -- An element of a ring of module; An elemento the codomain.
        
        OUTPUT:
            A element of the codomain.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import TrivialRepresentation
            sage: rep = TrivialRepresentation("ZZ", QQ)
            sage: rep.apply(2, 1/2)
            1/2
        """
        return a
    
    def __cmp__(self, other) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import TrivialRepresentation
            sage: rep1 = TrivialRepresentation("ZZ", QQ)
            sage: rep2 = TrivialRepresentation("ZZ", FreeModule(QQ, 1))
            sage: rep1 == rep1    # indirect doctest
            True
            sage: rep1 == rep2    # indirect doctest
            False
        """
        c = cmp(type(self), type(other))
        
        if c == 0 :
            c = cmp(self.__G, self.__G)
        if c == 0 :
            c = cmp(self.__K, other.__K)
            
        return c
    
    def __hash__(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import TrivialRepresentation
            sage: rep = TrivialRepresentation("ZZ", QQ)
            sage: hash(rep)    # indirect doctest
            -564657610            # 32-bit
            -11520069220171210    # 64-bit
        """
        return xor(hash(self.__G), hash(self.__K))

    def _repr_(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import TrivialRepresentation
            sage: rep = TrivialRepresentation("ZZ", QQ)
            sage: repr(rep)
            'Trivial representation of ZZ on Rational Field'
        """
        return "Trivial representation of %s on %s" % (self.__G, self.__K)
    
    def _latex_(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import TrivialRepresentation
            sage: rep = TrivialRepresentation("ZZ", QQ)
            sage: latex(rep)
            \text{Trivial representation of $\verb|ZZ|$ on $\Bold{Q}$}
        """
        return r"\text{Trivial representation of $%s$ on $%s$}" % (latex(self.__G), latex(self.__K))

#===============================================================================
# NNMonoid
#===============================================================================

class NNMonoid( SageObject ) :
    r"""
    Monoid of all natural numbers with zero as is can occure as an arguement of
    :class:`~fourier_expansion_framework.monoidpowerseries.MonoidPowerSeriesAmbient_abstract`.
    """
    
    def __init__(self, reduced = True) :
        r"""
        INPUT:
            - ``reduce``    -- A boolean (default: True); Wheter ``self`` is equipped with an action
                               or not.
        
        OUTPUT:
            An instance of :class:`~.NNmonoid`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: m = NNMonoid() 
        """
        self.__reduced = reduced

    def ngens(self) :
        r"""
        OUTPUT:
            An integer.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: m = NNMonoid()
            sage: m.ngens()
            1
        """
        return 1
    
    def gen(self, i = 0) :
        r"""
        OUTPUT:
            An integer.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: m = NNMonoid()
            sage: m.gen()
            1
            sage: m.gen(1)
            Traceback (most recent call last):
            ...
            ValueError: Generator not defined.
        """

        if i == 0 :
            return 1
        
        raise ValueError( "Generator not defined." )
    
    def gens(self) :
        r"""
        OUTPUT:
            A tuple of integers.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: m = NNMonoid()
            sage: m.gens()
            (1,)
        """    
        return tuple([self.gen(i) for i in xrange(self.ngens())])

    def is_commutative(self) :
        r"""
        Whether the monoid is commutative or not.
        
        OUTPUT:
            A boolean.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: m = NNMonoid()
            sage: m.is_commutative()
            True
        """
        return True
    
    def monoid(self) :
        r"""
        If ``self`` respects the action of the trivial group return the underlying
        monoid without this action. Otherwise return a copy of ``self``.

        OUTPUT:
            An instance of :class:`~.NNMonoid`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: m_a = NNMonoid()
            sage: m_woa = NNMonoid(False)
            sage: m_a.monoid() == m_woa
            True
            sage: m_woa.monoid() == m_woa
            True
            sage: m_woa.monoid() is m_woa
            False
        """
        return NNMonoid(False) 

    def group(self) :
        r"""
        If ``self`` respects the action of a group return representative.
        
        OUTPUT:
            Arbitrary type.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: m_a = NNMonoid()
            sage: m_woa = NNMonoid(False)
            sage: m_a.group()
            '1'
            sage: m_woa.group()
            Traceback (most recent call last):
            ...
            ArithmeticError: Monoid is not equipped with a group action.
        """
        if self.__reduced :
            return "1"
        else :
            raise ArithmeticError( "Monoid is not equipped with a group action.")

    def is_monoid_action(self) :
        r"""
        In case ``self`` respects the action of a group, decide whether this action is a monoid action
        on the underlying monoid.

        OUTPUT:
            A boolean.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: m_a = NNMonoid()
            sage: m_woa = NNMonoid(False)
            sage: m_a.is_monoid_action()
            1
            sage: m_woa.is_monoid_action()
            Traceback (most recent call last):
            ...
            ArithmeticError: Monoid is not equipped with a group action.
        """
        if self.__reduced :
            return True
        else :
            raise ArithmeticError( "Monoid is not equipped with a group action.")
    
    def filter(self, bound) :
        r"""
        Return a filter with given bound associated to this monoid.
        
        INPUT:
            - ``bound`` -- An integer; An upper bound for natural numbers.
        
        OUTPUT:
            An instance of :class:`~.NNFilter`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: m_a = NNMonoid()
            sage: m_woa = NNMonoid(False)
            sage: f = m_a.filter(5)
            sage: f.index()
            5
            sage: f.is_reduced()
            True
            sage: f = m_woa.filter(7)
            sage: f.is_reduced()
            False
        """
        return NNFilter(bound, self.__reduced)

    def filter_all(self) :
        r"""
        Return the filter associated to this monoid which contains all elements.
        
        OUTPUT:
            An instance of :class:`~.NNFilter`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: m_a = NNMonoid()
            sage: m_woa = NNMonoid(False)
            sage: f = m_a.filter_all()
            sage: f.index()
            +Infinity
            sage: f.is_reduced()
            True
            sage: m_woa.filter_all().is_reduced()
            False
        """
        return NNFilter(infinity, self.__reduced)
    
    def zero_filter(self) :
        r"""
        Return the filter associated to this monoid which contains no elements.
        
        OUTPUT:
            An instance of :class:`~.NNFilter`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: m_a = NNMonoid()
            sage: m_woa = NNMonoid(False)
            sage: f = m_a.zero_filter()
            sage: f.index()
            0
            sage: f.is_reduced()
            True
            sage: m_woa.zero_filter().is_reduced()
            False
        """
        return NNFilter(0, self.__reduced)         

    def minimal_composition_filter(self, ls, rs) :
        r"""
        Given two lists `ls` and `rs` of natural numbers return a filter that contains
        all the sums `l + r` of elements `l \in ls,\, r \in rs`.
        
        INPUT:
            - `ls` -- A list of integers.
            - `rs` -- A list of integers.
        
        OUTPUT:
            An instance of :class:`~.NNFilter`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: m = NNMonoid()
            sage: m.minimal_composition_filter([], []).is_reduced()
            True
            sage: m = NNMonoid(False)
            sage: m.minimal_composition_filter([], []).is_reduced()
            False
            sage: m = NNMonoid()
            sage: m.minimal_composition_filter([], []).index()
            0
            sage: m.minimal_composition_filter([], [1]).index()
            0
            sage: m.minimal_composition_filter([1], []).index()
            0
            sage: m.minimal_composition_filter([1], [1]).index()
            3
            sage: m.minimal_composition_filter([1,2,4], [1,5]).index()
            10
        """
        if len(ls) == 0 or len(rs) == 0 :
            return NNFilter(0, self.__reduced)

        return NNFilter(max(0, max(ls) + max(rs) + 1), self.__reduced)

    def _reduction_function(self) :
        r"""
        In case ``self`` respects the action of a group, return the
        reduction funtion for elements of this monoid.
        
        SEE::
            :meth:`~.reduce`
        
        OUTPUT:
            A function accepting one argument.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: m = NNMonoid()
            sage: m._reduction_function() == m.reduce
            True
            sage: m = NNMonoid(False)
            sage: m._reduction_function()
            Traceback (most recent call last):
            ...
            ArithmeticError: Monoid is not equipped with a group action.
        """
        if self.__reduced :
            return self.reduce
        else :
            raise ArithmeticError( "Monoid is not equipped with a group action." )

    def reduce(self, s) :
        r"""
        Reduce a natural number with respect to the trivial groups.
        
        INPUT:
            - `s` -- An integer.
        
        OUTPUT:
            The pair `(s, 1)`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: m = NNMonoid()
            sage: m.reduce(7)
            (7, 1)
        """
        return (s, 1)
  
    def decompositions(self, s) :
        r"""
        Decompose a natural number `s` in to a sum of two.
        
        INPUT:
            - `s` -- An integer.
        
        OUTPUT:
            A generator of pairs of integers.
        
        EXAMPLES::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: m = NNMonoid()
            sage: list(m.decompositions(4))
            [(0, 4), (1, 3), (2, 2), (3, 1), (4, 0)]
        """
        for n in xrange(s+1) :
            yield (n, s-n)
        
        raise StopIteration
    
    def zero_element(self) :
        r"""
        The zero element of this monoid.
        
        OUTPUT:
            An integer.
        
        TESTS:
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: m = NNMonoid()
            sage: m.zero_element()
            0
        """
        return 0

    def __contains__(self, x) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: m = NNMonoid()
            sage: 2 in m
            True
            sage: 'a' in m
            False
            sage: -1 in m
            False
        """
        return isinstance(x, (int, Integer)) and x >= 0

    def __cmp__(self, other) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: m = NNMonoid()
            sage: m == NNMonoid()
            True
            sage: m == NNMonoid(False)
            False
        """
        c = cmp(type(self), type(other))
        if c == 0 :
            c = cmp(self.__reduced, other.__reduced)
            
        return c
    
    def __hash__(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: hash(NNMonoid())
            1
            sage: hash(NNMonoid(False))
            0
        """
        return hash(self.__reduced)
    
    def _repr_(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: repr(NNMonoid())
            'NN with action'
            sage: repr(NNMonoid(False))
            'NN'
        """
        if not self.__reduced :
            return "NN"
        else :
            return "NN with action"
            
    def _latex_(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNMonoid
            sage: latex(NNMonoid())
            \Bold{N}\text{ with action}
            sage: latex(NNMonoid(False))
            \Bold{N}
        """
        if not self.__reduced :
            return r"\Bold{N}"
        else :
            return r"\Bold{N}\text{ with action}"

#===============================================================================
# NNFilter
#===============================================================================

class NNFilter ( SageObject ) :
    r"""
    A filter for the monoid of natrual numbers bounding elements by their value.
    """
    
    def __init__(self, bound, reduced = True) :
        r"""
        INPUT:
            - ``bound``   -- An integer; A bound for the natural numbers.
            - ``reduced`` -- A boolean (default: True); If True the action of the trivial
                             group is respected by the filter.
        
        OUTPUT:
            An instance of :class:`~.NNFilter`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNFilter
            sage: f = NNFilter(3) # indirect doctest
            sage: f = NNFilter(7, False) # indirect doctest
        """
        if isinstance(bound, NNFilter) :
            bound = bound.index()

        self.__bound = bound
        self.__reduced = reduced

    def is_infinite(self) :
        r"""
        Whether the filter contains infinitely many elements.
        
        OUTPUT:
            A boolean.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNFilter
            sage: NNFilter(3).is_infinite()
            False
            sage: NNFilter(infinity).is_infinite()
            True
        """
        return self.__bound is infinity

    def is_all(self) :
        r"""
        Whether the filter contains all elements of the associated monoid.
        
        OUTPUT:
            A boolean.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNFilter
            sage: NNFilter(3).is_all()
            False
            sage: NNFilter(infinity).is_all()
            True
        """
        return self.is_infinite()

    def index(self) :
        r"""
        Return the bound for this filter.
        
        OUTPUT:
            An integer.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNFilter
            sage: NNFilter(3).index()
            3
        """
        return self.__bound
    
    def is_reduced(self) :
        r"""
        Whether the filter respects the action of a group.
        
        OUTPUT:
            A boolean.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNFilter
            sage: NNFilter(3).is_reduced()
            True
            sage: NNFilter(7, False).is_reduced()
            False
        """
        return self.__reduced
    
    def __contains__(self, n) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNFilter
            sage: f = NNFilter(10)
            sage: 3 in f # indirect doctest
            True
            sage: 10 in f # indirect doctest
            False
            sage: 78 in NNFilter(infinity) # indirect doctest
            True
        """
        if self.__bound is infinity :
            return True
        
        return n < self.__bound
    
    def __iter__(self) :
        r"""
        OUTPUT:
            A generator over integers.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNFilter
            sage: list(NNFilter(4)) == range(4) # indirect doctest
            True
            sage: list(NNFilter(infinity))
            Traceback (most recent call last):
            ...
            ArithmeticError: Cannot iterate over infinite filters.
        """
        if self.__bound is infinity :
            raise ArithmeticError( "Cannot iterate over infinite filters." )
        
        for n in xrange(self.__bound) :
            yield n
        
        raise StopIteration
    
    def __cmp__(self, other) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNFilter
            sage: NNFilter(4) == NNFilter(4) # indirect doctest
            True
            sage: NNFilter(4) == NNFilter(5) # indirect doctest
            False
            sage: NNFilter(4) == NNFilter(4, False) # indirect doctest
            False
        """
        c = cmp(type(self), type(other))
        if c == 0 :
            c = cmp(self.__reduced, other.__reduced)
        if c == 0 :
            c = cmp(self.__bound, other.__bound)
                    
        return c

    def __hash__(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNFilter
            sage: hash(NNFilter(4)) # indirect doctest
            21
            sage: hash(NNFilter(7, False)) # indirect doctest
            7
        """
        return hash(self.__bound) + 17 * hash(self.__reduced)
                   
    def _repr_(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNFilter
            sage: repr(NNFilter(3))
            'Filtered NN with action up to 3'
            sage: repr(NNFilter(4, False))
            'Filtered NN up to 4'
        """
        if self.__reduced :
            return "Filtered NN with action up to %s" % self.__bound
        else :
            return "Filtered NN up to %s" % self.__bound
    
    def _latex_(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids import NNFilter
            sage: latex(NNFilter(3))
            \text{Filtered $\Bold{N}$ with action up to $3$}
            sage: latex(NNFilter(4, False))
            \text{Filtered $\Bold{N}$ up to $4$}
        """
        if self.__reduced :
            return r"\text{Filtered $\Bold{N}$ with action up to $%s$}" % latex(self.__bound)
        else :
            return r"\text{Filtered $\Bold{N}$ up to $%s$}" % latex(self.__bound)
