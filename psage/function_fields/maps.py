#################################################################################
#
# (c) Copyright 2010 William Stein
#
#  This file is part of PSAGE
#
#  PSAGE is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  PSAGE is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################


r"""

"""

from sage.categories.map import Map
from sage.categories.homset import Hom

from function_field_order import FunctionFieldOrder

class FunctionFieldIsomorphism(Map):
    r"""
    A base class for various isomorphisms between function fields and
    vector spaces.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<Y> = K.extension(y^2 - x*y + 4*x^3)
        sage: V, f, t = L.vector_space()
        sage: isinstance(f, sage.rings.function_field.maps.FunctionFieldIsomorphism)
        True
    """
    def _repr_type(self):
        """
        Return the type of this map (an isomorphism), for the purposes of printing out self.  
        
        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<Y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: f._repr_type()
            'Isomorphism'
        """
        return "Isomorphism"

    def is_injective(self):
        """
        Return True, since this isomorphism is injective.
        
        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<Y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()            
            sage: f.is_injective()
            True
        """
        return True

    def is_surjective(self):
        """
        Return True, since this isomorphism is surjective.
        
        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<Y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()            
            sage: f.is_surjective()
            True
        """
        return True

class MapVectorSpaceToFunctionField(FunctionFieldIsomorphism):
    r"""
    An isomorphism from a vector space and a function field.
    
    EXAMPLES:

        sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<Y> = K.extension(y^2 - x*y + 4*x^3)
        sage: V, f, t = L.vector_space(); f
        Isomorphism map:
          From: Vector space of dimension 2 over Rational function field in x over Rational Field
          To:   Function field in Y defined by y^2 - x*y + 4*x^3
    """
    def __init__(self, V, K):
        """
        EXAMPLES::
        
            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<Y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space(); type(f)
            <class 'sage.rings.function_field.maps.MapVectorSpaceToFunctionField'>
        """
        self._V = V
        self._K = K
        self._R = K.polynomial_ring()
        FunctionFieldIsomorphism.__init__(self, Hom(V, K))

    def _call_(self, v):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<Y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()                        
            sage: f(x*V.0 + (1/x^3)*V.1)         # indirect doctest
            1/x^3*Y + x
        """
        f = self._R(self._V(v).list())
        return self._K(f)

    def domain(self):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<Y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()                                    
            sage: f.domain()
            Vector space of dimension 2 over Rational function field in x over Rational Field
        """
        return self._V

    def codomain(self):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<Y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()                                    
            sage: f.codomain()
            Function field in Y defined by y^2 - x*y + 4*x^3
        """
        return self._K


class MapFunctionFieldToVectorSpace(FunctionFieldIsomorphism):
    """
    An isomorphism from a function field to a vector space.
    
    EXAMPLES::

        sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<Y> = K.extension(y^2 - x*y + 4*x^3)
        sage: V, f, t = L.vector_space(); t
        Isomorphism map:
          From: Function field in Y defined by y^2 - x*y + 4*x^3
          To:   Vector space of dimension 2 over Rational function field in x over Rational Field
    """
    def __init__(self, K, V):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<Y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space(); type(t)
            <class 'sage.rings.function_field.maps.MapFunctionFieldToVectorSpace'>
        """
        self._V = V
        self._K = K
        self._zero = K.base_ring()(0)
        self._n = K.degree()
        FunctionFieldIsomorphism.__init__(self, Hom(K, V))

    def domain(self):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<Y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()                                    
            sage: t.domain()
            Function field in Y defined by y^2 - x*y + 4*x^3
        """
        return self._K

    def codomain(self):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<Y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()                                    
            sage: t.codomain()
            Vector space of dimension 2 over Rational function field in x over Rational Field        
        """
        return self._V

    def _repr_type(self):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<Y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: t._repr_type()
            'Isomorphism'
        """
        return "Isomorphism"
    
    def _call_(self, x):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<Y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()                                    
            sage: t(x + (1/x^3)*Y)                       # indirect doctest
            (x, 1/x^3)
        """
        y = self._K(x)
        v = y.list()
        w = v + [self._zero]*(self._n - len(v))
        return self._V(w)


##########################################################################
# Morphisms between function fields

class FunctionFieldMorphism(Map):
    """
    EXAMPLES::

        sage: R.<x> = FunctionField(QQ); S.<y> = R[]
        sage: L.<y> = R.extension(y^2 - x^2)
        sage: f = L.hom(-y); f
        Morphism of function fields defined by y |--> -y    
    """
    def __init__(self, parent, im_gen, base_morphism):
        """
        EXAMPLES::

            sage: R.<x> = FunctionField(GF(7)); S.<y> = R[]
            sage: L.<y> = R.extension(y^3 + 6*x^3 + x); f = L.hom(y*2); f
            Morphism of function fields defined by y |--> 2*y
            sage: type(f)
            <class 'sage.rings.function_field.maps.FunctionFieldMorphism'>
            sage: factor(L.polynomial())
            y^3 + 6*x^3 + x
            sage: f(y).charpoly('y')
            y^3 + 6*x^3 + x
        """
        self._im_gen = im_gen
        self._base_morphism = base_morphism
        Map.__init__(self, parent)
        # Verify that the morphism is valid:
        R = self.codomain()['X']
        v = parent.domain().polynomial().list()
        if base_morphism is not None:
            v = [base_morphism(a) for a in v]
        f = R(v)
        if f(im_gen):
            raise ValueError, "invalid morphism"
    
    def is_injective(self):
        """
        EXAMPLES::

            sage: R.<x> = FunctionField(GF(7)); S.<y> = R[]
            sage: L.<y> = R.extension(y^3 + 6*x^3 + x); f = L.hom(y*2)
            sage: f.is_injective()
            True
        """
        return True

    def __repr__(self):
        """
        EXAMPLES::

            sage: R.<x> = FunctionField(GF(7)); S.<y> = R[]
            sage: L.<y> = R.extension(y^3 + 6*x^3 + x); f = L.hom(y*2)
            sage: f.__repr__()      
            'Morphism of function fields defined by y |--> 2*y'
        """
        return "Morphism of function fields defined by %s"%self._short_repr()

    def __nonzero__(self):
        """
        EXAMPLES::

            sage: R.<x> = FunctionField(GF(7)); S.<y> = R[]
            sage: L.<y> = R.extension(y^3 + 6*x^3 + x); f = L.hom(y*2)
            sage: f.__nonzero__()
            True
            sage: bool(f)
            True
        """
        return True

    def _short_repr(self):
        """
        EXAMPLES::

            sage: R.<x> = FunctionField(GF(7)); S.<y> = R[]
            sage: L.<y> = R.extension(y^3 + 6*x^3 + x); f = L.hom(y*2)
            sage: f._short_repr()
            'y |--> 2*y'
        """
        a = '%s |--> %s'%(self.domain().gen(), self._im_gen)
        if self._base_morphism is not None:
            a += ',  ' + self._base_morphism._short_repr()
        return a
    
    def _call_(self, x):
        """
        EXAMPLES::

            sage: R.<x> = FunctionField(GF(7)); S.<y> = R[]
            sage: L.<y> = R.extension(y^3 + 6*x^3 + x); f = L.hom(y*2)
            sage: f(y/x + x^2/(x+1))            # indirect doctest
            2/x*y + x^2/(x + 1)
            sage: f(y)
            2*y
        """
        v = x.list()
        if self._base_morphism is not None:
            v = [self._base_morphism(a) for a in v]
        f = v[0].parent()['X'](v)
        return f(self._im_gen)

class FunctionFieldMorphism_rational(FunctionFieldMorphism):
    def __init__(self, parent, im_gen):
        """
        EXAMPLES::

            sage: R.<x> = FunctionField(GF(7)); f = R.hom(1/x); f
            Morphism of function fields defined by x |--> 1/x
            sage: type(f)
            <class 'sage.rings.function_field.maps.FunctionFieldMorphism_rational'>
        """
        Map.__init__(self, parent)
        self._im_gen = im_gen
        self._base_morphism = None

    def _call_(self, x):
        """
        EXAMPLES::

            sage: R.<x> = FunctionField(GF(7)); f = R.hom(1/x); f
            Morphism of function fields defined by x |--> 1/x
            sage: f(x+1)                          # indirect doctest
            (x + 1)/x
            sage: 1/x + 1
            (x + 1)/x        
        """
        a = x.element()
        return a.subs({a.parent().gen():self._im_gen})
