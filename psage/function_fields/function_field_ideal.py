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


from sage.rings.ideal import Ideal_generic

class FunctionFieldIdeal(Ideal_generic):
    """
    A fractional ideal of a function field.

    EXAMPLES::

        sage: R.<x> = FunctionField(GF(7))
        sage: S = R.maximal_order(); I = S.ideal(x^3+1)
        sage: isinstance(I, sage.rings.function_field.function_field_ideal.FunctionFieldIdeal)
        True    
    """
    pass

class FunctionFieldIdeal_module(FunctionFieldIdeal):
    """
    A fractional ideal specified by a finitely generated module over
    the integers of the base field.

    EXAMPLES::

    An ideal in a rational function field::

        sage: R.<x> = FunctionField(QQ); S.<y> = R[]
        sage: L.<y> = R.extension(y^2 - x^3 - 1); M = L.equation_order()   
        sage: I = M.ideal(y)
        sage: I
        Ideal (x^3 + 1, -y) of Order in Function field in y defined by y^2 - x^3 - 1
        sage: I^2
        Ideal (x^3 + 1, (-x^3 - 1)*y) of Order in Function field in y defined by y^2 - x^3 - 1    
    """
    def __init__(self, ring, module):
        """
        INPUT:

            - ``ring`` -- an order in a function field
            - ``module`` -- a module

        EXAMPLES::

            sage: R.<x> = FunctionField(QQ); S.<y> = R[]
            sage: L.<y> = R.extension(y^2 - x^3 - 1); M = L.equation_order()   
            sage: I = M.ideal(y)
            sage: type(I)
            <class 'sage.rings.function_field.function_field_ideal.FunctionFieldIdeal_module'>
        """
        self._ring = ring
        self._module = module
        self._structure = ring.fraction_field().vector_space()
        V, from_V, to_V = self._structure
        gens = tuple([from_V(a) for a in module.basis()])
        Ideal_generic.__init__(self, ring, gens, coerce=False)

    def __contains__(self, x):
        """
        Return True if x is in this ideal.

        EXAMPLES::

            sage: R.<x> = FunctionField(GF(7)); S.<y> = R[]
            sage: L.<y> = R.extension(y^2 - x^3 - 1); M = L.equation_order()
            sage: I = M.ideal_with_gens_over_base([1, y]);  I
            Ideal (1, y) of Order in Function field in y defined by y^2 + 6*x^3 + 6
            sage: y in I
            True
            sage: y/x in I
            False
            sage: y^2 - 2 in I
            True
        """
        return self._structure[2](x) in self._module

    def module(self):
        """
        Return module over the maximal order of the base field that
        underlies self.

        The formation of this module is compatible with the vector
        space corresponding to the function field.

        OUTPUT:

            - a module over the maximal order of the base field of self

        EXAMPLES::

            sage: R.<x> = FunctionField(GF(7))
            sage: S = R.maximal_order(); S
            Maximal order in Rational function field in x over Finite Field of size 7
            sage: R.polynomial_ring()
            Univariate Polynomial Ring in x over Rational function field in x over Finite Field of size 7
            sage: I = S.ideal_with_gens_over_base([x^2 + 1, x*(x^2+1)])
            sage: I.gens()
            (x^2 + 1,)
            sage: I.module()
            Free module of degree 1 and rank 1 over Maximal order in Rational function field in x over Finite Field of size 7
            User basis matrix:
            [x^2 + 1]
            sage: V, from_V, to_V = R.vector_space(); V
            Vector space of dimension 1 over Rational function field in x over Finite Field of size 7
            sage: I.module().is_submodule(V)
            True        
        """
        return self._module

    def __add__(self, other):
        """
        Add together two ideals.
        
        EXAMPLES::
        
            sage: R.<x> = FunctionField(GF(7)); S.<y> = R[]
            sage: L.<y> = R.extension(y^2 - x^3 - 1); M = L.equation_order()   
            sage: I = M.ideal(y); J = M.ideal(y+1)
            sage: Z = I + J; Z
            Ideal (y + 1, 6*y) of Order in Function field in y defined by y^2 + 6*x^3 + 6
            sage: 1 in Z
            True
            sage: M.ideal(y^2) + M.ideal(y^3) == M.ideal(y^2,y^3)
            True
        """
        if not isinstance(other, FunctionFieldIdeal_module):
            other = self.ring().ideal(other)
        return FunctionFieldIdeal_module(self.ring(), self.module() + other.module())

    def intersection(self, other):
        """
        Return the intersection of the ideals self and other.
        
        EXAMPLES::

            sage: R.<x> = FunctionField(GF(7)); S.<y> = R[]
            sage: L.<y> = R.extension(y^2 - x^3 - 1); M = L.equation_order()   
            sage: I = M.ideal(y^3); J = M.ideal(y^2)
            sage: Z = I.intersection(J); Z
            Ideal (6*x^6 + 5*x^3 + 6, (6*x^3 + 6)*y) of Order in Function field in y defined by y^2 + 6*x^3 + 6
            sage: y^2 in Z
            False
            sage: y^3 in Z
            True
        """
        if not isinstance(other, FunctionFieldIdeal_module):
            other = self.ring().ideal(other)
        if self.ring() != other.ring():
            raise ValueError, "rings must be the same"
        return FunctionFieldIdeal_module(self.ring(), self.module().intersection(other.module()))

    def __cmp__(self, other):
        """
        Compare self and other.
        
        EXAMPLES::

            sage: R.<x> = FunctionField(GF(7)); S.<y> = R[]
            sage: L.<y> = R.extension(y^2 - x^3 - 1)
            sage: M = L.equation_order()   
            sage: I = M.ideal(y*(y+1)); J = M.ideal((y^2-2)*(y+1))
            sage: I+J == J+I            # indirect test
            True
            sage: I == J
            False
            sage: I < J
            True
            sage: J < I
            False        
        """
        if not isinstance(other, FunctionFieldIdeal_module):
            other = self.ring().ideal(other)
        if self.ring() != other.ring():
            raise ValueError, "rings must be the same"
        return cmp(self.module(), other.module())

    def __invert__(self):
        """
        Return the inverse of this fractional ideal.

        EXAMPLES::

            sage: R.<x> = FunctionField(GF(7)); S.<y> = R[]
            sage: L.<y> = R.extension(y^2 - x^3 - 1); M = L.equation_order()   
            sage: I = M.ideal(y)
            sage: I.__invert__()
            Ideal (1, (6/(x^3 + 1))*y) of Order in Function field in y defined by y^2 + 6*x^3 + 6
            sage: I^(-1)
            Ideal (1, (6/(x^3 + 1))*y) of Order in Function field in y defined by y^2 + 6*x^3 + 6
            sage: I.__invert__() * I
            Ideal (1, 6*y) of Order in Function field in y defined by y^2 + 6*x^3 + 6
        """
        if len(self.gens()) == 0:
            raise ZeroDivisionError

        # NOTE: If  I = (g0, ..., gn), then {x : x*I is in R}
        # is the intersection over i of {x : x*gi is in R}
        # Thus (I + J)^(-1) = I^(-1) intersect J^(-1).
        
        G = self.gens()
        R = self.ring()
        inv = R.ideal(~G[0])
        for g in G[1:]:
            inv = inv.intersection(R.ideal(~g))
        return inv

def ideal_with_gens(R, gens):
    """
    Return fractional ideal in the order of R with given generators
    over R.

    EXAMPLES::

        sage: R.<x> = FunctionField(QQ); S.<y> = R[]
        sage: L.<y> = R.extension(y^2 - x^3 - 1); M = L.equation_order()   
        sage: sage.rings.function_field.function_field_ideal.ideal_with_gens(M, [y])
        Ideal (x^3 + 1, -y) of Order in Function field in y defined by y^2 - x^3 - 1
    """
    K = R.fraction_field()
    return ideal_with_gens_over_base(R, [b*K(g) for b in R.basis() for g in gens])
    

def ideal_with_gens_over_base(R, gens):
    """
    Return fractional ideal in the given order R with given generators
    over the maximal order of the base field.

    EXAMPLES::

        sage: R.<x> = FunctionField(QQ); S.<y> = R[]
        sage: L.<y> = R.extension(y^2 - x^3 - 1); M = L.equation_order()   
        sage: sage.rings.function_field.function_field_ideal.ideal_with_gens_over_base(M, [x^3+1,-y])
        Ideal (x^3 + 1, -y) of Order in Function field in y defined by y^2 - x^3 - 1    
    """
    K = R.fraction_field()
    V, from_V, to_V = K.vector_space()

    # We handle the case of a rational function field separately,
    # since this is the base case and is used, e.g,. internally
    # by the linear algebra Hermite form code. 
    import function_field_order
    if isinstance(R, function_field_order.FunctionFieldOrder_rational):
        try:
            v = R._ring.ideal([x.element() for x in gens]).gens_reduced()
            assert len(v) == 1
            basis = [to_V(v[0])]
            M = V.span_of_basis(basis, check=False, already_echelonized=True, base_ring=R)
        except Exception, msg:
            print msg   # TODO --for debugging
            raise
    else:
        # General case
        S = V.base_field().maximal_order()
        M = V.span([to_V(b) for b in gens], base_ring=S)

    return FunctionFieldIdeal_module(R, M)
