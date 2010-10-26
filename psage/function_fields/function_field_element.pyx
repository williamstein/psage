include "../../ext/stdsage.pxi"

from sage.structure.element cimport FieldElement, RingElement, ModuleElement, Element

def is_FunctionFieldElement(x):
    """
    Return True if x is any type of function field element.
    
    EXAMPLES::

        sage: t = FunctionField(QQ,'t').gen()
        sage: sage.rings.function_field.function_field_element.is_FunctionFieldElement(t)
        True
        sage: sage.rings.function_field.function_field_element.is_FunctionFieldElement(0)
        False
    """
    return isinstance(x, FunctionFieldElement)

cdef class FunctionFieldElement(FieldElement):

    cdef readonly object _x
    cdef readonly object _matrix
    
    """
    The abstract base class for function field elements.
    
    EXAMPLES::

        sage: t = FunctionField(QQ,'t').gen()
        sage: isinstance(t, sage.rings.function_field.function_field_element.FunctionFieldElement)
        True
    """

    cdef FunctionFieldElement _new_c(self):
        cdef FunctionFieldElement x = <FunctionFieldElement>PY_NEW_SAME_TYPE(self)
        x._parent = self._parent
        return x
    
    def _latex_(self):
        """
        EXAMPLES::
        
            sage: K.<t> = FunctionField(QQ)
            sage: latex((t+1)/t)
            \frac{t + 1}{t}
            sage: latex((t+1)/t^67)
            \frac{t + 1}{t^{67}}
            sage: latex((t+1/2)/t^67)
            \frac{t + \frac{1}{2}}{t^{67}}
        """
        return self._x._latex_()
    
    def matrix(self):
        r"""
        Return the matrix of multiplication by self. 

        EXAMPLES::

        A rational function field:
        
            sage: K.<t> = FunctionField(QQ)
            sage: t.matrix()
            [t]
            sage: (1/(t+1)).matrix()
            [1/(t + 1)]

        Now an example in a nontrivial extension of a rational function field::
        
            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<Y> = K.extension(y^2 - x*y + 4*x^3)
            sage: Y.matrix()
            [     0      1]
            [-4*x^3      x]
            sage: Y.matrix().charpoly('Z')
            Z^2 - x*Z + 4*x^3

        An example in a relative extension, where neither function
        field is rational::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<Y> = K.extension(y^2 - x*y + 4*x^3)
            sage: M.<T> = L[]; Z.<alpha> = L.extension(T^3 - Y^2*T + x)
            sage: alpha.matrix()
            [          0           1           0]
            [          0           0           1]
            [         -x x*Y - 4*x^3           0]
        """
        if self._matrix is None:
            # Multiply each power of field generator on the left by this
            # element; make matrix whose rows are the coefficients of the
            # result, and transpose.
            K = self.parent()
            v = []
            x = K.gen()
            a = K(1)
            d = K.degree()
            for n in range(d):
                v += (a*self).list()
                a *= x
            k = K.base_ring()
            import sage.matrix.matrix_space
            M = sage.matrix.matrix_space.MatrixSpace(k, d)
            self._matrix = M(v)
        return self._matrix

    def trace(self):
        """
        Return the trace of this function field element.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<Y> = K.extension(y^2 - x*y + 4*x^3)
            sage: Y.trace()
            x
        """
        return self.matrix().trace()

    def norm(self):
        """
        Return the norm of this function field element.

        EXAMPLES::
        
            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<Y> = K.extension(y^2 - x*y + 4*x^3)
            sage: Y.norm()
            4*x^3


        The norm is relative::
        
            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<Y> = K.extension(y^2 - x*y + 4*x^3)
            sage: M.<T> = L[]; Z.<alpha> = L.extension(T^3 - Y^2*T + x)
            sage: alpha.norm()
            -x
            sage: alpha.norm().parent()
            Function field in Y defined by y^2 - x*y + 4*x^3
        """
        return self.matrix().determinant()

    def characteristic_polynomial(self, *args, **kwds):
        """
        Return the characteristic polynomial of this function field
        element.  Give an optional input string to name the variable
        in the characteristic polynomial.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<Y> = K.extension(y^2 - x*y + 4*x^3)
            sage: M.<T> = L[]; Z.<alpha> = L.extension(T^3 - Y^2*T + x)
            sage: x.characteristic_polynomial('W')
            W - x
            sage: Y.characteristic_polynomial('W')
            W^2 - x*W + 4*x^3
            sage: alpha.characteristic_polynomial('W')
            W^3 + (-x*Y + 4*x^3)*W + x
        """
        return self.matrix().characteristic_polynomial(*args, **kwds)

    charpoly = characteristic_polynomial

    def minimal_polynomial(self, *args, **kwds):
        """
        Return the minimal polynomial of this function field element.
        Give an optional input string to name the variable in the
        characteristic polynomial.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<Y> = K.extension(y^2 - x*y + 4*x^3)
            sage: M.<T> = L[]; Z.<alpha> = L.extension(T^3 - Y^2*T + x)
            sage: x.minimal_polynomial('W')
            W - x
            sage: Y.minimal_polynomial('W')
            W^2 - x*W + 4*x^3
            sage: alpha.minimal_polynomial('W')
            W^3 + (-x*Y + 4*x^3)*W + x
        """
        return self.matrix().minimal_polynomial(*args, **kwds)

    minpoly = minimal_polynomial
    
    def is_integral(self):
        r"""
        Determine if self is integral over the maximal order of the base field. 
        
        EXAMPLES::
        
            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<Y> = K.extension(y^2 - x*y + 4*x^3)
            sage: Y.is_integral()
            True
            sage: (Y/x).is_integral()
            True
            sage: (Y/x)^2 - (Y/x) + 4*x
            0
            sage: (Y/x^2).is_integral()
            False
            sage: f = (Y/x).minimal_polynomial('W'); f
            W^2 - W + 4*x            
        """
        R = self.parent().base_field().maximal_order()
        return all([a in R for a in self.minimal_polynomial()])

cdef class FunctionFieldElement_polymod(FunctionFieldElement):
    """
    EXAMPLES::

        sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<Y> = K.extension(y^2 - x*y + 4*x^3)
        sage: x*Y + 1/x^3
        x*Y + 1/x^3        
    """
    def __init__(self, parent, x):
        """
        EXAMPLES::

             sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<Y> = K.extension(y^2 - x*y + 4*x^3)
             sage: type(Y)
             <type 'sage.rings.function_field.function_field_element.FunctionFieldElement_polymod'>             
        """
        FieldElement.__init__(self, parent)
        self._x = x

    def element(self):
        """
        Return the underlying polynomial that represents this element.
        
        EXAMPLES::
            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<Y> = K.extension(y^2 - x*y + 4*x^3)
            sage: f = Y/x^2 + x/(x^2+1); f
            1/x^2*Y + x/(x^2 + 1)
            sage: f.element()
            1/x^2*y + x/(x^2 + 1)
            sage: type(f.element())
            <class 'sage.rings.polynomial.polynomial_element_generic.Polynomial_generic_dense_field'>        
        """
        return self._x

    def _repr_(self):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<Y> = K.extension(y^2 - x*y + 4*x^3)
            sage: Y._repr_()
            'Y'
        """
        return self._x._repr(name=self.parent().variable_name())
        
    def __nonzero__(self):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<Y> = K.extension(y^2 - x*y + 4*x^3)
            sage: bool(Y)
            True
            sage: bool(L(0))
            False
        """
        return not not self._x

    cdef int _cmp_c_impl(self, Element other) except -2:
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<Y> = K.extension(y^2 - x*y + 4*x^3)
            sage: cmp(L, 0)
            1
            sage: cmp(0, L)
            -1
        """
        cdef int c = cmp(type(self), type(other))
        if c: return c
        cdef FunctionFieldElement left = <FunctionFieldElement>self
        cdef FunctionFieldElement right = <FunctionFieldElement>other
        c = cmp(left._parent, right._parent)
        return c or cmp(left._x, right._x)

    cpdef ModuleElement _add_(self, ModuleElement right):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<Y> = K.extension(y^2 - x*y + 4*x^3)
            sage: (2*Y + x/(1+x^3))  +  (3*Y + 5*x*Y)         # indirect doctest
            (5*x + 5)*Y + x/(x^3 + 1)
        """
        cdef FunctionFieldElement res = self._new_c()
        res._x = self._x + (<FunctionFieldElement>right)._x
        return res

    cpdef ModuleElement _sub_(self, ModuleElement right):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<Y> = K.extension(y^2 - x*y + 4*x^3)
            sage: (2*Y + x/(1+x^3))  -  (3*Y + 5*x*Y)         # indirect doctest
            (-5*x - 1)*Y + x/(x^3 + 1)
        """
        cdef FunctionFieldElement res = self._new_c()
        res._x = self._x - (<FunctionFieldElement>right)._x
        return res

    cpdef RingElement _mul_(self, RingElement right):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<Y> = K.extension(y^2 - x*y + 4*x^3)
            sage: Y  *  (3*Y + 5*x*Y)                          # indirect doctest
            (5*x^2 + 3*x)*Y - 20*x^4 - 12*x^3
        """
        cdef FunctionFieldElement res = self._new_c()
        res._x = (self._x * (<FunctionFieldElement>right)._x) % self._parent.polynomial()
        return res

    cpdef RingElement _div_(self, RingElement right):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<Y> = K.extension(y^2 - x*y + 4*x^3)
            sage: (2*Y + x/(1+x^3))  /  (2*Y + x/(1+x^3))       # indirect doctest
            1
        """
        return self * ~right

    def __invert__(self):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<Y> = K.extension(y^2 - x*y + 4*x^3)
            sage: a = ~(2*Y + 1/x); a                           # indirect doctest
            (-x^2/(8*x^5 + x^2 + 1/2))*Y + (2*x^3 + x)/(16*x^5 + 2*x^2 + 1)
            sage: a*(2*Y + 1/x)
            1
        """
        P = self._parent
        return P(self._x.xgcd(P._polynomial)[1])

    def list(self):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]; L.<Y> = K.extension(y^2 - x*y + 4*x^3)
            sage: a = ~(2*Y + 1/x); a
            (-x^2/(8*x^5 + x^2 + 1/2))*Y + (2*x^3 + x)/(16*x^5 + 2*x^2 + 1)
            sage: a.list()
            [(2*x^3 + x)/(16*x^5 + 2*x^2 + 1), -x^2/(8*x^5 + x^2 + 1/2)]
            sage: (x*Y).list()
            [0, x]
        """
        return self._x.padded_list(self.parent().degree())
        

cdef class FunctionFieldElement_rational(FunctionFieldElement):
    """
    EXAMPLES::
    
        sage: FunctionField(QQ, 't')
        Rational function field in t over Rational Field
    """
    def __init__(self, parent, x):
        """
        EXAMPLES::
        
            sage: FunctionField(QQ,'t').gen()^3
            t^3
        """
        FieldElement.__init__(self, parent)
        self._x = x

    # nonoptimized

    def element(self):
        """
        Return the underlying fraction field element that represents this element.

        EXAMPLES::

            sage: R.<a> = FunctionField(GF(7))
            sage: a.element()
            a
            sage: type(a.element())
            <type 'sage.rings.fraction_field_element.FractionFieldElement'>        
        """
        return self._x

    def list(self):
        """
        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: t.list()
            [t]
        """
        return [self]

    def _repr_(self):
        """
        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: t._repr_()
            't'
        """
        return repr(self._x)
        
    def __nonzero__(self):
        """
        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: bool(t)
            True
            sage: bool(K(0))
            False
            sage: bool(K(1))
            True
        """
        return not not self._x

    cdef int _cmp_c_impl(self, Element other) except -2:
        """
        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: cmp(t, 0)
            1
            sage: cmp(t, t^2)
            -1
        """
        cdef int c = cmp(type(self), type(other))
        if c: return c
        cdef FunctionFieldElement left = <FunctionFieldElement>self
        cdef FunctionFieldElement right = <FunctionFieldElement>other
        c = cmp(left._parent, right._parent)
        return c or cmp(left._x, right._x)

    cpdef ModuleElement _add_(self, ModuleElement right):
        """
        EXAMPLES::
        
            sage: K.<t> = FunctionField(QQ)
            sage: t + (3*t^3)                      # indirect doctest
            3*t^3 + t
        """
        cdef FunctionFieldElement res = self._new_c()
        res._x = self._x + (<FunctionFieldElement>right)._x
        return res

    cpdef ModuleElement _sub_(self, ModuleElement right):
        """
        EXAMPLES::


            sage: K.<t> = FunctionField(QQ)
            sage: t - (3*t^3)                      # indirect doctest
            -3*t^3 + t
        """
        cdef FunctionFieldElement res = self._new_c()
        res._x = self._x - (<FunctionFieldElement>right)._x
        return res

    cpdef RingElement _mul_(self, RingElement right):
        """
        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: (t+1) * (t^2-1)                  # indirect doctest
            t^3 + t^2 - t - 1
        """
        cdef FunctionFieldElement res = self._new_c()
        res._x = self._x * (<FunctionFieldElement>right)._x
        return res

    cpdef RingElement _div_(self, RingElement right):
        """
        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: (t+1) / (t^2 - 1)                # indirect doctest
            1/(t - 1)
        """
        cdef FunctionFieldElement res = self._new_c()
        res._parent = self._parent.fraction_field()
        res._x = self._x / (<FunctionFieldElement>right)._x
        return res

    def numerator(self):
        """
        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: f = (t+1) / (t^2 - 1/3); f
            (t + 1)/(t^2 - 1/3)
            sage: f.numerator()
            t + 1
        """
        return self._x.numerator()

    def denominator(self):
        """
        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: f = (t+1) / (t^2 - 1/3); f
            (t + 1)/(t^2 - 1/3)
            sage: f.denominator()
            t^2 - 1/3
        """
        return self._x.denominator()
    
    def valuation(self, v):
        """
        EXAMPLES::
        
            sage: K.<t> = FunctionField(QQ)
            sage: f = (t-1)^2 * (t+1) / (t^2 - 1/3)^3
            sage: f.valuation(t-1)
            2
            sage: f.valuation(t)
            0
            sage: f.valuation(t^2 - 1/3)
            -3
        """
        R = self._parent._ring
        return self._x.valuation(R(self._parent(v)._x))
    
    def is_square(self):
        """
        Returns whether self is a square.
        
        EXAMPLES::
            
            sage: K.<t> = FunctionField(QQ)
            sage: t.is_square()
            False
            sage: (t^2/4).is_square()
            True
            sage: f = 9 * (t+1)^6 / (t^2 - 2*t + 1); f.is_square()
            True
            
            sage: K.<t> = FunctionField(GF(5))
            sage: (-t^2).is_square()
            True
            sage: (-t^2).sqrt()
            2*t
        """
        return self._x.is_square()
    
    def sqrt(self, all=False):
        """
        Returns the square root of self.
        
        EXMPLES::
        
            sage: K.<t> = FunctionField(QQ)
            sage: f = t^2 - 2 + 1/t^2; f.sqrt()
            (t^2 - 1)/t
            sage: f = t^2; f.sqrt(all=True)
            [t, -t]
            
        TESTS::
        
            sage: K(4/9).sqrt()
            2/3
            sage: K(0).sqrt(all=True)
            [0]
        """
        if all:
            return [self._parent(r) for r in self._x.sqrt(all=True)]
        else:
            return self._parent(self._x.sqrt())

    def factor(self):
        """
        Factor this rational function.
        
        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: f = (t+1) / (t^2 - 1/3)
            sage: f.factor()
            (t + 1) * (t^2 - 1/3)^-1
            sage: (7*f).factor()
            (7) * (t + 1) * (t^2 - 1/3)^-1
            sage: ((7*f).factor()).unit()
            7
            sage: (f^3).factor()
            (t + 1)^3 * (t^2 - 1/3)^-3
        """
        P = self.parent()
        F = self._x.factor()
        from sage.structure.factorization import Factorization
        return Factorization([(P(a),e) for a,e in F], unit=F.unit())

    def inverse_mod(self, I):
        r"""
        Return an inverse of self modulo the integral ideal `I`, if
        defined, i.e., if `I` and self together generate the unit
        ideal.

        EXAMPLES::

            sage: R.<x> = FunctionField(QQ)
            sage: S = R.maximal_order(); I = S.ideal(x^2+1)
            sage: t = S(x+1).inverse_mod(I); t
            -1/2*x + 1/2
            sage: (t*(x+1) - 1) in I
            True            
        """
        f = I.gens()[0]._x
        assert f.denominator() == 1
        assert self._x.denominator() == 1
        return self.parent()(self._x.numerator().inverse_mod(f.numerator()))
