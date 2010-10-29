from siegel_modular_form import SiegelModularForm, SiegelModularForm_class, SMF_DEFAULT_PREC
from sage.algebras.algebra import Algebra
from sage.misc.all import cached_method
from sage.rings.all import ZZ
from sage.structure.factory import UniqueFactory

class SiegelModularFormsAlgebraFactory(UniqueFactory):
    """
    A factory for creating algebras of Siegel modular forms.  It
    handles making sure that they are unique as well as handling
    pickling.  For more details, see
    :class:`~sage.structure.factory.UniqueFactory` and
    :mod:`~sage.modular.siegel.siegel_modular_forms_algebra`.

    EXAMPLES::

        sage: S = SiegelModularFormsAlgebra()
        sage: S.construction()
        (SMFAlg{"Sp(4,Z)", "even", 2, 101}, Integer Ring)
        sage: R.<a, b> = QQ[]
        sage: S = SiegelModularFormsAlgebra(coeff_ring=R)
        sage: S.construction()
        (SMFAlg{"Sp(4,Z)", "even", 2, 101}, Multivariate Polynomial Ring in a, b over Rational Field)
        sage: S is loads(dumps(S))
        True
        sage: TestSuite(S).run()
    """
    def create_key(self, coeff_ring=ZZ, group='Sp(4,Z)', weights='even', degree=2, default_prec=SMF_DEFAULT_PREC):
        """
        Create a key which uniquely defines this algebra of Siegel modular 
        forms.

        TESTS::

            sage: SiegelModularFormsAlgebra.create_key(coeff_ring=QQ, group='Gamma0(5)', weights='all', degree=2, default_prec=41)
            (SMFAlg{"Gamma0(5)", "all", 2, 41}(FractionField(...)), Integer Ring)
        """
        from sage.rings.ring import is_Ring
        if not is_Ring(coeff_ring):
            raise TypeError, 'The coefficient ring must be a ring'
        if not (isinstance(group, str) or group is None):
            raise TypeError, 'The group must be given by a string, or None'
        if not weights in ('even', 'all'):
            raise ValueError, "The argument weights must be 'even' or 'all'" 
        try:
            degree = int(degree)
        except TypeError:
            raise TypeError, 'The degree must be a positive integer'
        if degree < 1:
            raise ValueError, 'The degree must be a positive integer'
        if degree == 1:
            raise ValueError, 'Use ModularForms if you want to work with Siegel modular forms of degree 1'
        if degree > 2:
            raise NotImplementedError, 'Siegel modular forms of degree > 2 are not yet implemented'
        try:
            default_prec = int(default_prec)
        except TypeError:
            raise TypeError, 'The default precision must be a positive integer'
        
        from pushout import SiegelModularFormsAlgebraFunctor
        F = SiegelModularFormsAlgebraFunctor(group=group, weights=weights, degree=degree, default_prec=default_prec) 

        while hasattr(coeff_ring, 'construction'):
            C = coeff_ring.construction()
            if C is None:
                break
            F = F * C[0]
            coeff_ring = C[1]
        return (F, coeff_ring)

    def create_object(self, version, key):
        """
        Return the algebra of Siegel modular forms corresponding to the 
        key ``key``.

        TESTS::

            sage: key = SiegelModularFormsAlgebra.create_key(coeff_ring=QQ, group='Gamma0(5)', weights='all', degree=2, default_prec=41)
            sage: SiegelModularFormsAlgebra.create_object('1.0', key)
            Algebra of Siegel modular forms of degree 2 and all weights on Gamma0(5) over Rational Field
        """
        C, R = key
        from pushout import CompositConstructionFunctor, SiegelModularFormsAlgebraFunctor
        if isinstance(C, CompositConstructionFunctor):
            F = C.all[-1]
            if len(C.all) > 1:
                R = CompositConstructionFunctor(*C.all[:-1])(R)
        else:
            F = C
        if not isinstance(F, SiegelModularFormsAlgebraFunctor):
            raise TypeError, "We expected a SiegelModularFormsAlgebraFunctor, not %s"%type(F)
        return SiegelModularFormsAlgebra_class(coeff_ring=R, group=F._group, weights=F._weights, degree=F._degree, default_prec=F._default_prec)
                       
SiegelModularFormsAlgebra = SiegelModularFormsAlgebraFactory('SiegelModularFormsAlgebra')


class SiegelModularFormsAlgebra_class(Algebra):
    def __init__(self, coeff_ring=ZZ, group='Sp(4,Z)', weights='even', degree=2, default_prec=SMF_DEFAULT_PREC):
        r"""
        Initialize an algebra of Siegel modular forms of degree ``degree`` 
        with coefficients in ``coeff_ring``, on the group ``group``.  
        If ``weights`` is 'even', then only forms of even weights are 
        considered; if ``weights`` is 'all', then all forms are 
        considered.

        EXAMPLES::

            sage: A = SiegelModularFormsAlgebra(QQ)
            sage: B = SiegelModularFormsAlgebra(ZZ)
            sage: A._coerce_map_from_(B)
            True                                                                                                      
            sage: B._coerce_map_from_(A)
            False                                                                                                                                            
            sage: A._coerce_map_from_(ZZ)
            True   
        """
        self.__coeff_ring = coeff_ring
        self.__group = group
        self.__weights = weights
        self.__degree = degree
        self.__default_prec = default_prec
        R = coeff_ring
        from sage.algebras.all import GroupAlgebra
        if isinstance(R, GroupAlgebra):
            R = R.base_ring()
        from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
        if is_PolynomialRing(R):
            self.__base_ring = R.base_ring()
        else:
            self.__base_ring = R
        from sage.categories.all import Algebras
        Algebra.__init__(self, base=self.__base_ring, category=Algebras(self.__base_ring))
    
    @cached_method
    def gens(self, prec=None):
        """
        Return the ring generators of this algebra.

        EXAMPLES::

            sage: S = SiegelModularFormsAlgebra()
            sage: S.gens()
            [Igusa_4, Igusa_6, Igusa_10, Igusa_12]

        TESTS::

            sage: S = SiegelModularFormsAlgebra(group='Gamma0(2)')
            sage: S.gens()
            Traceback (most recent call last):
            ...
            NotImplementedError: Not yet implemented
        """
        if prec is None:
            prec = self.default_prec()
        return _siegel_modular_forms_generators(self, prec=prec)

    def is_commutative(self):
        r"""
        Return True since all our algebras are commutative.

        EXAMPLES:: 

            sage: S = SiegelModularFormsAlgebra(coeff_ring=QQ)
            sage: S.is_commutative()
            True
        """
        return True

    @cached_method
    def ngens(self):
        r"""
        Return the number of generators of this algebra of Siegel modular
        forms.

        EXAMPLES::

            sage: S = SiegelModularFormsAlgebra(coeff_ring=QQ) 
            sage: S.ngens()
            4
        """
        return len(self.gens())

    @cached_method
    def gen(self, i, prec=None):
        r"""
        Return the `i`-th generator of this algebra.

        EXAMPLES::

            sage: S = SiegelModularFormsAlgebra(coeff_ring=QQ)
            sage: S.ngens()
            4
            sage: S.gen(2)
            Igusa_10
        """
        return self.gens(prec=prec)[i]

    def coeff_ring(self):
        r"""
        Return the ring of coefficients of this algebra.

        EXAMPLES::

            sage: S = SiegelModularFormsAlgebra(coeff_ring=QQ)
            sage: S.coeff_ring()
            Rational Field
        """
        return self.__coeff_ring

    def group(self):
        """
        Return the modular group of this algebra.

        EXAMPLES::

            sage: S = SiegelModularFormsAlgebra(group='Gamma0(7)')
            sage: S.group()
            'Gamma0(7)'
        """
        return self.__group

    def weights(self):
        """
        Return 'even' or 'all' depending on whether this algebra
        contains only even weight forms, or forms of all weights.

        EXAMPLES::

            sage: S = SiegelModularFormsAlgebra(weights='even')
            sage: S.weights()
            'even'
            sage: S = SiegelModularFormsAlgebra(weights='all')
            sage: S.weights()
            'all'
        """
        return self.__weights

    def degree(self):
        """
        Return the degree of the modular forms in this algebra.

        EXAMPLES::

            sage: S = SiegelModularFormsAlgebra()
            sage: S.degree()
            2
        """
        return self.__degree

    def default_prec(self):
        """
        Return the default precision for the Fourier expansions of
        modular forms in this algebra.

        EXAMPLES::

            sage: S = SiegelModularFormsAlgebra(default_prec=41)
            sage: S.default_prec()
            41
        """
        return self.__default_prec
    
    def _repr_(self):
        r"""
        Return the plain text representation of this algebra.

        EXAMPLES::

            sage: S = SiegelModularFormsAlgebra(coeff_ring=QQ)
            sage: S._repr_()
            'Algebra of Siegel modular forms of degree 2 and even weights on Sp(4,Z) over Rational Field'
        """
        return 'Algebra of Siegel modular forms of degree %s and %s weights on %s over %s'%(self.__degree, self.__weights, self.__group, self.__coeff_ring) 

    def _latex_(self):
        r"""
        Return the latex representation of this algebra.

        EXAMPLES::
            
            sage: S = SiegelModularFormsAlgebra(coeff_ring=QQ)
            sage: S._latex_()
            '\\texttt{Algebra of Siegel modular forms of degree }2\\texttt{ and even weights on Sp(4,Z) over }\\Bold{Q}'
        """
        from sage.misc.all import latex
        return r'\texttt{Algebra of Siegel modular forms of degree }%s\texttt{ and %s weights on %s over }%s' %(latex(self.__degree), self.__weights, self.__group, latex(self.__coeff_ring))

    def _coerce_map_from_(self, other):
        r"""
        Return True if it is possible to coerce from ``other`` into the
        algebra of Siegel modular forms ``self``.

        EXAMPLES::

            sage: S = SiegelModularFormsAlgebra(coeff_ring=QQ)
            sage: R = SiegelModularFormsAlgebra(coeff_ring=ZZ)
            sage: S._coerce_map_from_(R)
            True
            sage: R._coerce_map_from_(S)
            False
            sage: S._coerce_map_from_(ZZ)
            True
        """
        if self.base_ring().has_coerce_map_from(other):
            return True
        if isinstance(other, SiegelModularFormsAlgebra_class):
            if (self.group() == other.group()) or (other.group() == 'Sp(4,Z)'):
                if self.coeff_ring().has_coerce_map_from(other.coeff_ring()):
                    if self.degree() == other.degree():
                        return True                                               
        return False

    def _element_constructor_(self, x):
        r"""
        Return the element of the algebra ``self`` corresponding to ``x``.

        EXAMPLES::

            sage: S = SiegelModularFormsAlgebra(coeff_ring=QQ)
            sage: B = SiegelModularFormsAlgebra(coeff_ring=ZZ).1
            sage: S(B)
            Igusa_6
            sage: S(1/5)
            1/5
            sage: S(1/5).parent() is S
            True
            sage: S._element_constructor_(2.67)
            Traceback (most recent call last):
            ...
            TypeError: Unable to construct an element of Algebra of Siegel modular forms of degree 2 and even weights on Sp(4,Z) over Rational Field corresponding to 2.67000000000000
            sage: S.base_extend(RR)._element_constructor_(2.67)
            2.67000000000000
        """
        if isinstance(x, (int, long)):
            x = ZZ(x)
        if isinstance(x, float):
            from sage.rings.all import RR
            x = RR(x)
        if isinstance(x, complex):
            from sage.rings.all import CC
            x = CC(x)

        if isinstance(x.parent(), SiegelModularFormsAlgebra_class):
            d = dict((f, self.coeff_ring()(x[f])) for f in x.coeffs())
            return self.element_class(parent=self, weight=x.weight(), coeffs=d, prec=x.prec(), name=x.name())

        R = self.base_ring()
        if R.has_coerce_map_from(x.parent()):
            d = {(0, 0, 0): R(x)}
            from sage.rings.all import infinity
            return self.element_class(parent=self, weight=0, coeffs=d, prec=infinity, name=str(x))
        else:
            raise TypeError, "Unable to construct an element of %s corresponding to %s" %(self, x)

    Element = SiegelModularForm_class

    def _an_element_(self):
        r"""
        Return an element of the algebra ``self``.

        EXAMPLES::

            sage: S = SiegelModularFormsAlgebra(coeff_ring=QQ)
            sage: z = S._an_element_()
            sage: z in S
            True
        """
        return self(0)

    def base_extend(self, R):
        r"""
        Extends the base ring of the algebra ``self`` to ``R``.

        EXAMPLES::

            sage: S = SiegelModularFormsAlgebra(coeff_ring=QQ)
            sage: S.base_extend(RR)
            Algebra of Siegel modular forms of degree 2 and even weights on Sp(4,Z) over Real Field with 53 bits of precision
        """
        #B = self.base_ring()
        S = self.coeff_ring()
        from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
        if is_PolynomialRing(S):
            xS = S.base_extend(R)
        elif R.has_coerce_map_from(S):
            xS = R
        else:
            raise TypeError, "cannot extend to %s" %R
        return SiegelModularFormsAlgebra(coeff_ring=xS, group=self.group(), weights=self.weights(), degree=self.degree(), default_prec=self.default_prec())

    def construction(self):
        r"""
        Return the construction of the algebra ``self``.

        EXAMPLES::

            sage: S = SiegelModularFormsAlgebra(coeff_ring=QQ)
            sage: S.construction()
            (SMFAlg{"Sp(4,Z)", "even", 2, 101}, Rational Field)
        """
        from pushout import SiegelModularFormsAlgebraFunctor
        return SiegelModularFormsAlgebraFunctor(self.group(), self.weights(), self.degree(), self.default_prec()), self.coeff_ring()



def _siegel_modular_forms_generators(parent, prec=None, degree=0):
    r"""
    Compute the four Igusa generators of the ring of Siegel modular forms 
    of degree 2 and level 1 and even weight (this happens if weights='even').  
    If weights = 'all' you get the Siegel modular forms of degree 2, level 1 
    and even and odd weight.

    EXAMPLES::

        sage: A, B, C, D = SiegelModularFormsAlgebra().gens()
        sage: C[(2, 0, 9)]
        -390420
        sage: D[(4, 2, 5)]
        17689760
        sage: A2, B2, C2, D2 = SiegelModularFormsAlgebra().gens(prec=50)
        sage: A[(1, 0, 1)] == A2[(1, 0, 1)]
        True
        sage: A3, B3, C3, D3 = SiegelModularFormsAlgebra().gens(prec=500)
        sage: B2[(2, 1, 3)] == B3[(2, 1, 3)]
        True

    TESTS::

        sage: from sage.modular.siegel.siegel_modular_forms_algebra import _siegel_modular_forms_generators
        sage: S = SiegelModularFormsAlgebra()
        sage: S.gens() == _siegel_modular_forms_generators(S)
        True
    """
    group = parent.group()
    weights = parent.weights()
    if prec is None:
        prec = parent.default_prec()
    if group == 'Sp(4,Z)' and 0 == degree:
        from sage.modular.all import ModularForms
        E4 = ModularForms(1, 4).gen(0)
        M6 = ModularForms(1, 6)
        E6 = ModularForms(1, 6).gen(0)
        M8 = ModularForms(1, 8)
        M10 = ModularForms(1, 10)
        Delta = ModularForms(1, 12).cuspidal_subspace().gen(0)
        M14 = ModularForms(1, 14)
        A = SiegelModularForm(60*E4, M6(0), prec=prec, name='Igusa_4')
        B = SiegelModularForm(-84*E6, M8(0), prec=prec, name='Igusa_6')
        C = SiegelModularForm(M10(0), -Delta, prec=prec, name='Igusa_10')
        D = SiegelModularForm(Delta, M14(0), prec=prec, name='Igusa_12')
        # TODO: the base_ring of A, B, ... should be ZZ
        # Here a hack for now:
        a = [A, B, C, D]
        b = []         
        from sage.rings.all import ZZ
        for F in a:
            c = F.coeffs()
            for f in c.iterkeys():
                c[f] = ZZ(c[f])
            F = parent.element_class(parent=parent, weight=F.weight(), coeffs=c, prec=prec, name=F.name())
            b.append(F)
        if weights == 'even': return b
        if weights == 'all':
            from fastmult import chi35
            coeffs35 = chi35(prec, b[0], b[1], b[2], b[3])
            from sage.groups.all import KleinFourGroup
            G = KleinFourGroup()
            from sage.algebras.all import GroupAlgebra
            R = GroupAlgebra(G)
            det = R(G.gen(0))
            E = parent.element_class(parent=parent, weight=35, coeffs=coeffs35, prec=prec, name='Delta_35')
            E = E*det
            b.append(E)
            return b
        raise ValueError, "weights = '%s': should be 'all' or 'even'" %weights
    if group == 'Sp(4,Z)' and weights == 'even' and 2 == degree:
        b = _siegel_modular_forms_generators(parent=parent)
        c = []
        for F in b:
            i = b.index(F)
            for G in b[(i+1):]:
                c.append(F.satoh_bracket(G))
        return c
    raise NotImplementedError, "Not yet implemented"

