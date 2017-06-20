r"""
Siegel modular forms

Implementation of a class describing (scalar and vector valued) Siegel
modular forms of degree 2 and arbitrary group and weight.


AUTHORS:

- CC: Craig Citro

- MR: Martin Raum

- NR: Nathan Ryan

- NS: Nils-Peter Skoruppa


HISTORY:

- NS: Class SiegelModularForm_class.

- NS: Factory function SiegelModularForm().
      Code parts concerning the Maass lift are partly based on code
      written by David Gruenewald in August 2006 which in turn is
      based on the PARI/GP code by Nils-Peter Skoruppa from 2003.

- CC, NR: _mult_, _reduce_GL using Cython.

- NR: Hecke operators.

- CC: Rewriting of fastmult.spyx for arbitrary base rings.

- NS: Class Morp.

- MR, NS: SiegelModularFormsAlgebra_class, SiegelModularFormsAlgebra  


REFERENCES:

- [Sko] Nils-Peter Skoruppa, ...

- [I-H] Tomoyoshi Ibukiyama and Shuichi Hayashida, ... 
"""

#*****************************************************************************
#  Copyright (C) 2008 Nils-Peter Skoruppa <nils.skoruppa@uni-siegen.de>
#                     
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import cPickle
import urllib
from siegel_modular_form_prec import SiegelModularFormPrecision
from sage.rings.all import ZZ
from sage.structure.element import AlgebraElement

SMF_DEFAULT_PREC = 101

class SiegelModularForm_class(AlgebraElement):
    r"""
    Describes a Siegel modular form of degree 2.
    """
    def __init__(self, parent, weight, coeffs, prec=SMF_DEFAULT_PREC, name=None):
        r"""
        Create a Siegel modular form of degree 2.

        INPUT:

        - ``parent`` -- SiegelModularFormsAlgebra (ambient space).

        - ``weight`` -- the weight of the form.

        - ``coeffs`` -- a dictionary whose keys are triples `(a,b,c)` 
          of integers representing `GL(2,\ZZ)`-reduced quadratic forms 
          (i.e. forms satisfying `0 \le b \le a \le c`).

        - ``name`` -- an optional string giving a name to the form,
          e.g. 'Igusa_4'.  This modifies the text and latex
          representations of the form.

        OUTPUT:

        - a Siegel modular form

        EXAMPLES::

            sage: E10 = ModularForms(1, 10).0
            sage: Delta = CuspForms(1, 12).0
            sage: F = SiegelModularForm(E10, Delta); F
            Siegel modular form on Sp(4,Z) of weight 10

        TESTS::

            sage: TestSuite(F).run() 
        """
        self.__group = parent.group()
        self.__weight = weight
        self.__coeffs = coeffs
        # TODO: check whether we have enough coeffs
        self.__coeff_ring = parent.coeff_ring()
        self.__precision = SiegelModularFormPrecision(prec)
        self.__prec = _normalized_prec(prec)
        self.__name = name
        self.__rep_lists = dict()
        AlgebraElement.__init__(self, parent)

    def base_ring(self):
        """
        Return the ring of coefficients of the Siegel modular form ``self``.

        EXAMPLES::

            sage: S = SiegelModularFormsAlgebra()
            sage: A = S.gen(0)
            sage: A.base_ring()
            Integer Ring
            sage: C = S.gen(2)
            sage: onehalfC = C * (-1/2)
            sage: onehalfC.base_ring()
            Rational Field
            sage: B = S.gen(1)
            sage: AB = A.satoh_bracket(B)
            sage: AB.base_ring()
            Multivariate Polynomial Ring in x, y over Rational Field
        """
        return self.__coeff_ring
    
    def _repr_(self):
        r"""
        Return the plain text representation of ``self``.

        EXAMPLES::

            sage: C = SiegelModularFormsAlgebra().gen(2)
            sage: C._repr_()
            'Igusa_10'
            sage: (C^2)._repr_()
            'Siegel modular form on Sp(4,Z) of weight 20'
        """
        if self.name() is None:
            return 'Siegel modular form on %s of weight %d' % (self.group(), self.weight())
        else:
            return self.name()

    def _latex_(self):
        r"""
        Return the latex representation of ``self``.

        EXAMPLES::

            sage: A = SiegelModularFormsAlgebra().gen(0)
            sage: A._latex_()
            'Igusa_4'
            sage: (A^2)._latex_()
            '\\texttt{Siegel modular form on }Sp(4,Z)\\texttt{ of weight }8'
        """
        if self.name() is None:
            return r'\texttt{Siegel modular form on }%s\texttt{ of weight }%d' % (self.group(), self.weight())
        else:
            return self.name()

    def group(self):
        r"""
        Return the modular group of ``self`` as a string.

        EXAMPLES::

            sage: C = SiegelModularFormsAlgebra().gen(2)
            sage: C.group()
            'Sp(4,Z)'
        """
        return self.__group

    def weight(self):
        r"""
        Return the weight of ``self``.

        EXAMPLES::

            sage: A, B, C, D = SiegelModularFormsAlgebra().gens()
            sage: A.weight()
            4
            sage: B.weight()
            6
            sage: C.weight()
            10
            sage: D.weight()
            12
            sage: (A * B).weight()
            10
            sage: (A*B + C).weight()
            10
       """
        return self.__weight

    def precision(self):
        r"""
        Return the Precision class object that describes the precision of ``self``.

        EXAMPLES::

            sage: A = SiegelModularFormsAlgebra().gen(0)
            sage: A.precision()
            Discriminant precision for Siegel modular form with bound 101
        """
        return self.__precision

    def __getitem__(self, key):
        r"""
        Return the coefficient indexed by ``key`` in the Fourier expansion
        of ``self``.

        INPUT:

        - ``key`` -- a triple of integers `(a, b, c)` defining a quadratic form
        
        OUTPUT:

        - the coefficient of `q^{(a,b,c)}` as an element of the coefficient
          ring of ``self``

        EXAMPLES::

            sage: A, B, C, D = SiegelModularFormsAlgebra().gens()
            sage: A[(0, 0, 10)]
            272160
            sage: A[(0, 0, 25)]
            3780240
            sage: A[(5, 0, 5)]
            97554240
            sage: A[(0, 0, 26)]
            Traceback (most recent call last):
            ...
            ValueError: precision 101 is not high enough to extract coefficient at (0, 0, 26)
            sage: A[(100, 0, 1)]
            Traceback (most recent call last):    
            ...
            ValueError: precision 101 is not high enough to extract coefficient at (100, 0, 1)
            sage: C[(2, 1, 3)]
            2736
            sage: C[(3, 1, 2)]
            2736

        Any coefficient indexed by a quadratic form that is not semi-positive
        definite is by definition zero::

            sage: A[(-1, 1, 1)]
            0
            sage: A[(1, 10, 1)]
            0

        TESTS::

            sage: A[(1/2, 0, 6)]
            Traceback (most recent call last):
            ...
            TypeError: the key (1/2, 0, 6) must be an integral quadratic form
        """
        a, b, c = key  
        try:
            a = ZZ(a)
            b = ZZ(b)
            c = ZZ(c)
        except TypeError:
            raise TypeError, "the key %s must be an integral quadratic form" %str(key)

        if b**2 - 4*a*c > 0  or a < 0 or c < 0:
            return self.base_ring()(0)
        ## otherwise we have a semi-positive definite form
        if self.__coeffs.has_key((a, b, c)):
            return self.__coeffs[(a, b, c)]
        ## otherwise GL2(ZZ)-reduce (a,b,c) and try again
        from fastmult import reduce_GL
        (a0, b0, c0) = reduce_GL(a, b, c)
        if self.precision().is_in_bound((a0, b0, c0)):
##            return get_coeff_with_action( a,b,c,self.coeffs(),self.base_ring())
            return self.__coeffs.get((a0, b0, c0), self.base_ring()(0))
        ## otherwise error - precision is too low 
        raise ValueError, 'precision %s is not high enough to extract coefficient at %s' %(self.prec(), str(key))

    def coeffs(self, disc=None):
        r"""
        Return the dictionary of coefficients of ``self``.  If ``disc`` is 
        specified, return the dictionary of coefficients indexed by 
        semi-positive definite quadratic forms of discriminant ``disc``.

        INPUT:

        - ``disc`` -- optional (default: None); a negative integer giving
          the discriminant of a semi-positive definite quadratic form

        EXAMPLES::

            sage: A, B, C, D = SiegelModularFormsAlgebra().gens(prec=20)
            sage: C.coeffs().has_key((3, 0, 1))
            False
            sage: C.coeffs().has_key((1, 0, 3))
            True
            sage: C.coeffs()[(1, 0, 3)]
            -272
            sage: len(D.coeffs(disc=-16).keys())
            2
            sage: D.coeffs(disc=-16)[(2, 0, 2)]
            17600

        TESTS::

            sage: A.coeffs(disc=-20)
            Traceback (most recent call last):
            ...
            ValueError: precision is not high enough to extract coefficients of discriminant -20
        """
        if disc is None:
            return self.__coeffs
        elif -disc < self.prec():
            if disc < 0 and (disc % 4) in [0, 1]:
                from sage.quadratic_forms.binary_qf import BinaryQF_reduced_representatives
                forms = [f[:] for f in BinaryQF_reduced_representatives(disc)]
                return dict([(Q, self[Q]) for Q in forms])
            else:
                return {}
        else:
            raise ValueError, 'precision is not high enough to extract coefficients of discriminant %s' %disc

    def support(self):
        r"""
        Return the support of ``self``, i.e. the list of reduced quadratic
        forms indexing non-zero Fourier coefficients of ``self``.

        EXAMPLES::

            sage: A, B, C, D = SiegelModularFormsAlgebra().gens(prec=15)
            sage: A.support()
            [(0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 0, 3), (1, 0, 1), (1, 0, 2), (1, 0, 3), (1, 1, 1), (1, 1, 2), (1, 1, 3), (2, 2, 2)]
            sage: (0, 0, 2) in C._keys()
            True
            sage: C[(0, 0, 2)]
            0
            sage: (0, 0, 2) in C.support()
            False
        """
        return [k for k in self._keys() if self[k] != 0]

    def max_disc(self):
        r"""
        Return the largest discriminant corresponding to a non-zero Fourier 
        coefficient of ``self``.

        Note that since these discriminants are all non-positive, this 
        corresponds in some sense to the "smallest" non-zero Fourier 
        coefficient.  It is analogous to the valuation of a power series 
        (=smallest degree of a non-zero term).

        EXAMPLES::

            sage: gens = SiegelModularFormsAlgebra().gens()
            sage: [f.max_disc() for f in gens]
            [0, 0, -3, -3]
        """    
        return max([b**2 - 4*a*c for (a, b, c) in self.support()])

    def _keys(self):
        r"""
        Return the keys of the dictionary of coefficients of ``self``.

        EXAMPLES::

            sage: A, B, C, D = SiegelModularFormsAlgebra().gens(prec=15)
            sage: A._keys()
            [(0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 0, 3), (1, 0, 1), (1, 0, 2), (1, 0, 3), (1, 1, 1), (1, 1, 2), (1, 1, 3), (2, 2, 2)]
        """
        return sorted(self.coeffs().keys())

    def prec(self):
        r"""
        Return the precision of ``self``. 

        EXAMPLES::

            sage: B = SiegelModularFormsAlgebra().gen(1)
            sage: B.prec()
            101
        """
        return self.precision().prec()

    def name(self):
        r"""
        Return the name of the Siegel modular form ``self``, or None if 
        ``self`` does not have a custom name.

        EXAMPLES::

            sage: A, B = SiegelModularFormsAlgebra().gens()[:2]
            sage: A.name()
            'Igusa_4'
            sage: (A*B).name() is None
            True
        """
        return self.__name

    def truncate(self, prec):
        r"""
        Return a Siegel modular form with the same parent as ``self``, but
        with Fourier expansion truncated at the new precision ``prec``.

        EXAMPLES::

            sage: A = SiegelModularFormsAlgebra().gen(0, prec=20); A
            Igusa_4
            sage: A.prec()
            20
            sage: len(A._keys())
            18
            sage: A[(1, 1, 5)]
            1330560
            sage: At = A.truncate(16); At
            Igusa_4
            sage: At.prec()
            16
            sage: len(At._keys())
            14
            sage: At[(1, 1, 5)]
            Traceback (most recent call last):
            ...
            ValueError: precision 16 is not high enough to extract coefficient at (1, 1, 5)
            sage: At[(1, 1, 4)] == A[(1, 1, 4)]
            True

        TESTS::
            sage: A.truncate(30)
            Traceback (most recent call last):
            ...
            ValueError: truncated precision 30 cannot be larger than the current precision 20
        """
        if prec > self.prec():
            raise ValueError, 'truncated precision %s cannot be larger than the current precision %s' %(prec, self.prec())
        else:
            return _SiegelModularForm_from_dict(group=self.group(), weight=self.weight(), coeffs=self.coeffs(), prec=prec, parent=None, name=self.name())

    ################
    ## Arithmetic ##
    ################

    def _add_(left, right):
        r"""
        Return the sum of ``left`` and ``right``. 

        There are some compatibility conditions that the forms ``left``
        and ``right`` must satisfy in order for addition to work:

        - they must have the same weight
        - they must be defined on the same group, or one of them must
          be defined on the whole group 'Sp(4,Z)'

        The precision of the sum of ``left`` and ``right`` is the minimum
        of the precisions of ``left`` and ``right``.

        EXAMPLES::

            sage: A, B, C, D = SiegelModularFormsAlgebra().gens()
            sage: twoA = A + A
            sage: A2 = A * 2
            sage: A2 == twoA
            True
            sage: A*B + C
            Siegel modular form on Sp(4,Z) of weight 10
            sage: E = C._add_(A*B)
            sage: E
            Siegel modular form on Sp(4,Z) of weight 10
            sage: C[(1, 1, 1)]
            1
            sage: (A*B)[(1, 1, 1)]
            57792
            sage: E[(1, 1, 1)]
            57793

            sage: F = SiegelModularForm('Gamma0(2)', 4, A.coeffs())
            sage: F
            Siegel modular form on Gamma0(2) of weight 4
            sage: A + F
            Siegel modular form on Gamma0(2) of weight 4
            sage: F + A == F + A
            True
            sage: (A + F)[(1, 1, 1)] == A[(1, 1, 1)] + F[(1, 1, 1)]
            True
            sage: A1 = SiegelModularFormsAlgebra().gen(0, prec=20)
            sage: (A1 + F).prec()
            20
            sage: (A1 + F)[(1, 1, 1)] == (A + F)[(1, 1, 1)]
            True

        TESTS::

            sage: A + B
            Traceback (most recent call last):
            ...
            ValueError: cannot add Siegel modular forms of different weights

            sage: F = SiegelModularForm('Gamma0(2)', 4, A.coeffs())
            sage: G = SiegelModularForm('Gamma0(5)', 4, A.coeffs())
            sage: F + G
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '+': 'Algebra of Siegel modular forms of degree 2 and even weights on Gamma0(2) over Integer Ring' and 'Algebra of Siegel modular forms of degree 2 and even weights on Gamma0(5) over Integer Ring' 
            sage: Z = A + (-A)
            sage: Z
            Siegel modular form on Sp(4,Z) of weight 4
            sage: Z.prec()
            101
            sage: Z.coeffs()
            {}
        """
        # take the minimum of the precisions and truncate the forms if necessary
        lp = left.prec()
        rp = right.prec()
        if lp < rp:
            prec = lp
            right = right.truncate(prec)
        elif rp < lp:
            prec = rp
            left = left.truncate(prec)
        else:
            prec = lp

        # figure out the group of the sum = intersection of the groups
        # of the two forms
        # right now, we only allow two identical groups, or different
        # groups where one of them is all of 'Sp(4,Z)'
        lg = left.group()
        rg = right.group()
        if lg is None or lg == 'Sp(4,Z)': 
            group = rg
        elif rg is None or rg == 'Sp(4,Z)': 
            group = lg
        elif lg != rg:
            raise NotImplementedError, "addition of forms on different groups not yet implemented"
        else: 
            group = lg
        
        # the sum of two Siegel modular forms is another Siegel modular
        # form only if the weights are equal
        if left.weight() != right.weight():
            raise ValueError, 'cannot add Siegel modular forms of different weights'
        else:
            wt = left.weight()

        # compute the coefficients of the sum
        d = dict()
        for f in set(left._keys() + right._keys()):
            v = left[f] + right[f]
            if not v.is_zero(): d[f] = v
        par = left.parent()
        from sage.modular.siegel.siegel_modular_forms_algebra import SiegelModularFormsAlgebra
        par2 = SiegelModularFormsAlgebra(coeff_ring=par.coeff_ring(), group=group, weights=par.weights(), degree=par.degree())
        return par2.element_class(parent=par, weight=wt, coeffs=d, prec=prec)
        
    def _mul_(left, right):
        r"""
        Return the product of the Siegel modular form ``left`` by the
        element ``right``, which can be a scalar or another Siegel modular
        form.

        The multiplication of two forms on arbitrary groups is not yet
        implemented.  At the moment, the forms must either be on the same
        group, or one of the groups must be the full 'Sp(4,Z)'.

        EXAMPLES::

            sage: A, B, C, D = SiegelModularFormsAlgebra().gens()
            sage: A3 = A * 3
            sage: B2 = 2 * B
            sage: Conehalf = C / 2
            sage: Donehalf = (1/2) * D
            sage: A3[(1, 1, 1)] / A[(1, 1, 1)]
            3
            sage: B2[(2, 1, 3)] / B[(2, 1, 3)]
            2
            sage: Conehalf[(2, 0, 1)] / C[(2, 0, 1)]
            1/2
            sage: Donehalf[(3, 3, 3)] / D[(3, 3, 3)]
            1/2
            sage: E = C._mul_(D)
            sage: E[(0, 0, 5)]
            0
            sage: E[(2, 3, 3)]
            8
        
        TESTS::

            sage: A = SiegelModularFormsAlgebra(default_prec=61).0
            sage: A.prec()
            61
            sage: C = SiegelModularFormsAlgebra(default_prec=81).2
            sage: C.prec()
            81
            sage: A*C == C*A
            True
            sage: (A*C).prec()
            64
        """                                 
        from sage.modular.siegel.siegel_modular_forms_algebra import SiegelModularFormsAlgebra
        par = left.parent()
        from sage.rings.all import infinity
        if left.prec() is infinity:
            tmp = right 
            right = left
            left = tmp

        wt = left.weight() + right.weight()
        if wt % 2:
            weights = 'all'
        else:
            weights = 'even'
        if right.prec() is infinity:
            prec = left.prec()
            c = right[(0, 0, 0)]
            d = dict()
            for f in left.coeffs():
                v = left[f] * c
                if not v.is_zero(): d[f] = v
            par = SiegelModularFormsAlgebra(coeff_ring=par.coeff_ring(), group=par.group(), weights=weights, degree=par.degree())
            return par.element_class(parent=par, weight=left.weight(), coeffs=d, prec=prec)

        lp = left.prec()
        rp = right.prec()
        if lp < rp:
            right = right.truncate(lp)
        elif rp < lp:
            left = left.truncate(rp)
        prec = min(lp - right.max_disc(), rp - left.max_disc())

        if left.group() is None: 
            group = right.group()
        elif right.group() is None: 
            group = left.group()
        elif left.group() != right.group():
            raise NotImplementedError, "multiplication for differing groups not yet implemented"
        else: 
            group = left.group()

        s1, s2, R = left.coeffs(), right.coeffs(), left.base_ring()
        d = dict()
        _prec = SiegelModularFormPrecision(prec)
        if ZZ == R:
            from fastmult import mult_coeff_int
            for x in _prec:
                v = mult_coeff_int(x[0], x[1], x[2], s1, s2)
                if not v.is_zero(): d[x] = v
        elif left.parent().base_ring() == left.parent().coeff_ring():
            from fastmult import mult_coeff_generic
            for x in _prec:
                v = mult_coeff_generic(x[0], x[1], x[2], s1, s2, R)
                if not v.is_zero(): d[x] = v
        else:
            from fastmult import mult_coeff_generic_with_action
            for x in _prec:
                v = mult_coeff_generic_with_action(x[0], x[1], x[2], s1, s2, R)
                if not v.is_zero(): d[x] = v
            
        par = SiegelModularFormsAlgebra(coeff_ring=par.coeff_ring(), group=group, weights=weights, degree=par.degree())
        return par.element_class(parent=par, weight=wt, coeffs=d, prec=prec)

    def _rmul_(self, c):
        r"""
        Right multiplication -- only by scalars.

        EXAMPLES::

            sage: A, B, C, D = SiegelModularFormsAlgebra().gens()
            sage: onehalfC = C._rmul_(-1/2)
            sage: C[(2, 1, 3)]
            2736
            sage: onehalfC[(2, 1, 3)]
            -1368

        TESTS::

            sage: Z = A._rmul_(0)
            sage: Z
            Siegel modular form on Sp(4,Z) of weight 4
            sage: Z.coeffs()
            {}
        """
        d = dict()
        for f in self.coeffs():
            v = self[f] * c
            if not v.is_zero(): d[f] = v
        par = self.parent()
        return par.element_class(parent=par, weight=self.weight(), coeffs=d, prec=self.prec())

    def _lmul_(self, c):
        r"""
        Left multiplication -- only by scalars.
        
        EXAMPLES::

            sage: A = SiegelModularFormsAlgebra(default_prec=51).0
            sage: A3 = A._lmul_(3)
            sage: A[(1, 1, 1)]
            13440
            sage: A3[(1, 1, 1)]
            40320
        """
        return self._rmul_(c)

    def __eq__(left, right):
        r"""
        Return True if the Siegel modular forms ``left`` and ``right`` 
        have the same group and weight, and the same coefficients up to 
        the smaller of the two precisions.

        EXAMPLES::

            sage: A, B, C, D = SiegelModularFormsAlgebra(default_prec=20).gens()
            sage: A2, B2, C2, D2 = SiegelModularFormsAlgebra().gens() 
            sage: A2 == A
            True
            sage: A == A
            True
            sage: A == B
            False

        TESTS::

            sage: E = SiegelModularForm('Gamma0(2)', 6, B.coeffs())
            sage: E == B
            False
            sage: S = SiegelModularFormsAlgebra(default_prec=21)
            sage: S(0) == S(0)
            True
            sage: S(0) == S(1)
            False
        """
        if not isinstance(right, SiegelModularForm_class):
            return False
        if left.group() != right.group():
            return False
        if left.weight() != right.weight():
            return False

        # we compare up to the minimum of the two precisions
        new_prec = min(left.prec(), right.prec())
        # if new_prec is infinity, then left and right are both
        # constant Siegel modular forms
        from sage.rings.all import infinity
        if new_prec is infinity:
            return left[0, 0, 0] == right[0, 0, 0]

        left = left.truncate(new_prec)
        right = right.truncate(new_prec)
        new_smf_prec = SiegelModularFormPrecision(new_prec)
        
        lc = left.coeffs()
        rc = right.coeffs()
        # make coefficients that are implicitly zero be actually zero
        for k in new_smf_prec:
            if not lc.has_key(k):
                lc[k] = 0
            if not rc.has_key(k):
                rc[k] = 0
            
        return lc == rc

    def pickle(self, file, format=1, name=None):
        r"""
        Dump to a file in a portable format.

        EXAMPLES::

            sage: A = SiegelModularFormsAlgebra(default_prec=20).0
            sage: FILE = open('A.sobj', 'w')
            sage: A.pickle(FILE)

        TO DO
            I don't really think I have the right syntax since I couldn't
            load it after doing this
        """
        from sage.rings.all import QQ
        if self.base_ring().fraction_field() == QQ:
            pol = [0, 1]
            coeffs = self.coeffs()
        else:
            pol = self.base_ring().polynomial().list()
            coeffs = dict()
            for k in self.coeffs():
                coeffs[k] = self.coeffs()[k].list()
        if None == name:
            _name = self._repr_()
        f1 =  'format %d: [this string, wt, name, pol., prec., dict. of coefficients]' %1
        data = [f1, self.weight(), _name, pol, self.prec(), coeffs]
        cPickle.dump(data, file)
        return

    def fourier_jacobi_coeff(self, N):
        r"""
        Return the `N`-th Fourier-Jacobi coefficient of the Siegel modular 
        form ``self`` as a dictionary indexed by pairs `(disc,r)` with 
        `disc<0` and `r^2 \equiv disc \bmod 4N`.
        
        EXAMPLES::

            sage: A = SiegelModularFormsAlgebra(default_prec=51).0
            sage: fj = A.fourier_jacobi_coeff(1)
            sage: fj[(-8, 0)]
            181440
        """
        prec = self.prec() 
        NN = 4*N
        keys = [(disc, r) for disc in range(prec) for r in range(NN) if (r**2 - disc)%NN == 0]
        coeff = dict()
        for disc, r in keys:
            (a, b, c) = (-(r**2 - disc)/NN, r, N) 
            coeff[(-disc, r)] = self[(a, b, c)]
        return coeff

    def satoh_bracket(left, right, names = ['x', 'y']):
        r"""
        Return the Satoh bracket [``left``, ``right``], where ``left``
        and ``right`` are scalar-valued Siegel modular forms.

        The result is a vector-valued Siegel modular form whose coefficients
        are polynomials in the two variables specified by ``names``.

        EXAMPLES::

            sage: A, B = SiegelModularFormsAlgebra().gens()[:2]
            sage: AB = A.satoh_bracket(B)
            sage: AB[(1, 0, 1)]
            -20160*x^2 - 40320*y^2
        """
        from sage.rings.all import PolynomialRing
        R = PolynomialRing(ZZ, names)
        x, y = R.gens()
        d_left_coeffs = dict((f, (f[0]*x*x + f[1]*x*y + f[2]*y*y)*left[f])
                            for f in left.coeffs())
        d_left = SiegelModularForm(left.group(), left.weight(),
                                              d_left_coeffs, left.prec())
        d_right_coeffs = dict((f, (f[0]*x*x + f[1]*x*y + f[2]*y*y)*right[f])
                             for f in right.coeffs())
        d_right = SiegelModularForm(right.group(), right.weight(),
                                               d_right_coeffs, right.prec())
        return d_left*right/left.weight() - d_right*left/right.weight()

    #################
    # Hecke operators
    #################

    def hecke_image(self, ell):
        r"""
        Return the Siegel modular form which is the image of ``self`` under
        the Hecke operator T(``ell``).

        EXAMPLES::

            sage: A, B, C, D = SiegelModularFormsAlgebra().gens()
            sage: Ups20 = -1/1785600*A*B*C - 1/1785600*A^2*D + C^2
            sage: TUps20 = Ups20.hecke_image(2)
            sage: TUps20[(1, 1, 1)] / Ups20[(1, 1, 1)]
            -840960
            sage: TUps20
            Siegel modular form on Sp(4,Z) of weight 20
            sage: A.hecke_image(5)
            T(5)(Igusa_4)
        """
        # TODO: compute precision for T(ell)F
        from sage.functions.all import ceil
        try:
            # TODO: I am not sure whether this sets the right prec
            a, b, c = self.prec()
            prec = (ceil(a/ell), ceil(b/ell), ceil(c/ell))
        except TypeError:
            prec = ceil(self.prec()/ell/ell)
            prec = _normalized_prec(prec)
        d = dict((f, self.hecke_coefficient(ell, f)) for f in self._keys() if _is_bounded(f, prec))
        if self.name():
            name = 'T(' + str(ell) + ')(' + self.name() + ')'
        else:
            name = None
        par = self.parent()
        from sage.modular.siegel.siegel_modular_forms_algebra import SiegelModularFormsAlgebra
        par = SiegelModularFormsAlgebra(coeff_ring=par.coeff_ring(), group=par.group(), weights=par.weights(), degree=par.degree())
        return par.element_class(parent=par, weight=self.weight(), coeffs=d, prec=prec, name=name)

    def hecke_coefficient(self, ell, (a, b, c)):
        r"""
        Return the `(a, b, c)`-indexed coefficient of the image of ``self``
        under the Hecke operator T(``ell``).
        
        EXAMPLES::

            sage: A, B, C, D = SiegelModularFormsAlgebra().gens()
            sage: Ups20 = -1/1785600*A*B*C - 1/1785600*A^2*D + C^2
            sage: Ups20.hecke_coefficient(5, (1, 1, 1))
            5813608045/992
            sage: Ups20.hecke_coefficient(5, (1, 0, 1))
            5813608045/248             
            sage: Ups20.hecke_coefficient(5, (1, 0, 1)) / Ups20[(1, 0, 1)]
            -5232247240500
        """
        k = self.weight()
        coeff = 0
        from sage.rings.all import divisors
        from sage.quadratic_forms.binary_qf import BinaryQF
        qf = BinaryQF([a, b, c])
        for t1 in divisors(ell):
            for t2 in divisors(t1):
                cosets = self._get_representatives(ell, t1/t2)
                for V in cosets:                       
                    aprime, bprime, cprime = qf.matrix_action_right(V)[:]
                    if aprime % t1 == 0 and bprime % t2 == 0 and cprime % t2 == 0:
                        try:
                            coeff = coeff + t1**(k-2)*t2**(k-1)*self[(ell*aprime/t1**2, ell*bprime/t1/t2, ell*cprime/t2**2)]
                        except KeyError, msg:
                            raise ValueError, '%s' %(self,msg)
        return coeff

    def _get_representatives(self, ell, t):
        r"""
        A helper function used in hecke_coefficient that computes the right
        coset representatives of `\Gamma^0(t)\SL(2,Z)` where `\Gamma^0(t)` 
        is the subgroup of `SL(2,Z)` where the upper left hand corner is 
        divisible by `t`.
 
        EXAMPLES::

            sage: A = SiegelModularFormsAlgebra().0
            sage: A._get_representatives(5, 3)
            [
            [ 0  1]  [1 0]  [1 1]  [1 2]
            [-1  0], [0 1], [0 1], [0 1]
            ]

        NOTE
            We use the bijection $\Gamma^0(t)\SL(2,Z) \rightarrow P^1(\Z/t\Z)$
            given by $A \mapsto [1:0]A$.
         """
        try:
            return self.__rep_lists[(ell, t)]
        except KeyError:
            from sage.matrix.all import MatrixSpace
            M = MatrixSpace(ZZ, 2, 2)
        if t == 1:
            return [M([1, 0, 0, 1])]
        from sage.modular.all import P1List
        P = list(P1List(t))
        from sage.rings.all import IntegerModRing, xgcd
        ZZt = IntegerModRing(t)
        rep_list = []
        for x, y in P1List(t):
            e, d, c = xgcd(x, y)
            rep = M([x, y, -c, d])
            rep_list.append(rep)
        self.__rep_lists[(ell, t)] = rep_list
        return rep_list



def SiegelModularForm(arg0, arg1=None, arg2=None, prec=None, name=None, hint=None):
    r"""
    Construct a Siegel modular form.

    INPUT:

        Supported formats

        1. SiegelModularForm(f, g):
           creates the Siegel modular form VI(f,g) from the elliptic
           modular forms f, g as in [Sko].
           Shortforms:
               SiegelModularForm(f) for SiegelModularForm(f, 0-form),
               SiegelModularForm(f, 0) for SiegelModularForm(f, 0-form),
               SiegelModularForm(0, f) for SiegelModularForm(0-form, f).

        2. SiegelModularForm(x):
           if x is an element of a commutative ring creates the constant 
           Siegel modular form with value x.

        3. SiegelModularForm(string):
           creates the Siegel modular form pickled at the location
           described by string (url or path_to_file).

        4. SiegelModularForm(qf):
           constructs the degree 2 theta series associated to the given
           quadratic form.

        5. SiegelModularForm(f):
           --- TODO: Borcherds lift

        6. SiegelModularForm(f, g):
           --- TODO: Yoshida lift

        7. SiegelModularForm(repr):
           --- TODO: Lift from Weil representation

        8. SiegelModularForm(char):
           --- TODO: Singular weight

        9. SiegelModularForm(group, weight, dict):
           if group is a string, weight an integer and dict a dictionary,
           creates a Siegel modular form whose coefficients are contained
           in dict.
            
    EXAMPLES::

        sage: M4 = ModularForms(1, 4)                                        
        sage: M6 = ModularForms(1, 6)
        sage: E4 = M4.basis()[0]         
        sage: F = SiegelModularForm(E4, M6(0), 16); F
        Siegel modular form on Sp(4,Z) of weight 4

        sage: url = 'http://sage.math.washington.edu/home/nils/Siegel-Modular-Forms/data/upsilon-forms_20-32_1000/Upsilon_20.p'  # optional -- internet
        sage: H = SiegelModularForm(url); H  # optional - internet
        Upsilon_20
        sage: H[(4, 4, 4)]                      # optional - internet
        248256200704

        sage: url = 'http://sage.math.washington.edu/home/nils/Siegel-Modular-Forms/data/upsilon-forms_20-32_1000/Upsilon_28.p'   # optional -- internet
        sage: L = SiegelModularForm(url); L  # optional - internet
        Upsilon_28
        sage: L.prec()                        # optional - internet
        1001
        sage: L[(3, 3, 3)]                      # optional - internet
        -27352334316369546*a^2 + 3164034718941090318*a + 2949217207771097198880
        sage: L[(3, 3, 3)].parent()             # optional - internet
        Number Field in a with defining polynomial x^3 - x^2 - 294086*x - 59412960

        sage: l = [(0, 0, 0, 0)]*2 + [(1, 0, 0, 0)] + [(0, 1, 0, 0)] + [(0, 0, 0, 1)] + [(0, 0, 1, 1)]
        sage: char_dict = {tuple(l): 1/4}
        sage: F = SiegelModularForm(char_dict, prec=100, name="van Geemen F_7"); F
        van Geemen F_7
        sage: F.coeffs()[(1, 0, 9)]
        -7
 
        sage: Q37 = QuadraticForm(ZZ, 4, [1,0,1,1, 2,2,3, 10,2, 6])
        sage: Q37.level()
        37
        sage: F37 = SiegelModularForm(Q37, prec=100); F37
        Siegel modular form on Gamma0(37) of weight 2
        sage: F37[(2, 1, 3)]
        0
        sage: F37[(2, 1, 7)]
        4
    """
    try:
        from sage.structure.element import py_scalar_to_element
        arg0 = py_scalar_to_element(arg0)
    except TypeError:
        pass

    from sage.modular.modform.element import ModularFormElement
    from sage.rings.all import RingElement
    if isinstance(arg0, ModularFormElement)\
           and isinstance(arg1, ModularFormElement):
        return _SiegelModularForm_as_Maass_spezial_form(arg0, arg1, prec, name)
    if isinstance(arg0, ModularFormElement) \
           and (0 == arg1 or arg1 is None):
        M = ModularForms(1, arg0.weight() + 2)
        return _SiegelModularForm_as_Maass_spezial_form(arg0, M(0), prec, name)
    if 0 == arg0 and isinstance(arg1, ModularFormElement):
        M = ModularForms(1, arg1.weight() - 2)
        return _SiegelModularForm_as_Maass_spezial_form(M(0), arg1, prec, name)
    from sage.quadratic_forms.all import QuadraticForm
    if isinstance(arg0, QuadraticForm):
        return _SiegelModularForm_from_QuadraticForm(arg0, prec, name)
    if isinstance(arg0, RingElement) and arg1 is None:
        from sage.rings.all import infinity
        return _SiegelModularForm_from_dict(group='Sp(4,Z)', weight=0, coeffs={(0, 0, 0): arg0}, prec=infinity)
    if isinstance(arg0, str) and arg1 is None:
        return _SiegelModularForm_from_file(arg0)
    if isinstance(arg0, dict) and arg1 is None:
        return _SiegelModularForm_from_theta_characteristics(arg0, prec=prec, name=name, hint=hint)
    from sage.rings.all import Integer
    if isinstance(arg0, str) and isinstance(arg1, (int, Integer)) and isinstance(arg2, dict):
        return _SiegelModularForm_from_dict(group=arg0, weight=arg1, coeffs=arg2, prec=prec, name=name)
    raise TypeError, "wrong arguments"


def _SiegelModularForm_as_Maass_spezial_form(f, g, prec=SMF_DEFAULT_PREC, name=None):
    """
    Return the Siegel modular form  I(f,g) (Notation as in [Sko]).

    EXAMPLES::
    
        sage: M14 = ModularForms(group=1, weight=14)
        sage: E14 = M14.eisenstein_subspace()
        sage: f = M14.basis()[0]
        sage: S16 = ModularForms(group=1, weight=16).cuspidal_subspace()
        sage: g = S16.basis()[0]
        sage: from sage.modular.siegel.siegel_modular_form import _SiegelModularForm_as_Maass_spezial_form
        sage: IFG = _SiegelModularForm_as_Maass_spezial_form(f, g, prec=100, name=None)
        sage: IFG[(2, 1, 3)]
        -1080946527072
        
    INPUT
        f:   -- modular form of level 1
        g:   -- cusp form of level 1 amd wt = wt of f + 2
        prec -- either a triple (amax,bmac,cmax) or an integer Dmax
    """
    k = f.weight()
    assert(k+2 == g.weight()) | (f==0) | (g==0), "incorrect weights!"
    assert(g.q_expansion(1) == 0), "second argument is not a cusp form" 

    if isinstance(prec, tuple):
        (amax, bmax, cmax) = prec
        amax = min(amax, cmax)
        bmax = min(bmax, amax)
        clean_prec = (amax, bmax, cmax)
        if bmax <= 0:
            # no reduced forms below prec
            return _SiegelModularForm_from_dict(group='Sp(4,Z)', weight=k, coeffs=dict(), prec=0)
        if 1 == amax:
            # here prec = (0,0,>=0)
            Dtop = 0
        else:
            Dtop = 4*(amax-1)*(cmax-1)
    else:
        clean_prec = max(0, prec)
        if 0 == clean_prec:
            # no reduced forms below prec
            return _SiegelModularForm_from_dict(group='Sp(4,Z)', weight=k, coeffs=dict(), prec=0)
        while 0 != clean_prec%4 and 1 != clean_prec%4:
            clean_prec -= 1
        Dtop = clean_prec - 1
    precision = (Dtop+1)//4 + 1
    # TODO: examine error when called with 1 == prec
    if 1 == precision:
        precision = 2
    
    """
    Create the Jacobi form I(f,g) as in [Sko].

    It suffices to construct for all Jacobi forms phi only the part
    sum_{r=0,1;n} c_phi(r^2-4n) q^n zeta^r.
    When, in this code part, we speak of Jacobi form we only mean this part.
    
    We need to compute Ifg = \sum_{r=0,1; n} c(r^2-4n) q^n zeta^r up to
    4n-r^2 <= Dtop, i.e. n < precision
    """
    ## print 'Creating I(f,g)'

    from sage.rings.all import PowerSeriesRing, QQ
    PS = PowerSeriesRing(QQ, name='q')
    q = PS.gens()[0]

    ## Create the quasi Dedekind eta^-6 power series:
    pari_prec = max(1, precision - 1) # next line yields error if 0 == pari_prec
    from sage.libs.pari.gen import pari
    from sage.rings.all import O
    pari.set_series_precision(pari_prec)
    eta_quasi = PS(pari('Vec(eta(q))')) + O(q**precision)
    etapow = eta_quasi**-6

    ## Create the Jacobi forms A=a*etapow and B=b*etapow in stages.
    ## Recall a = sum_{s != r mod 2} s^2*(-1)^r*q^((s^2+r^2-1)/4)*zeta^r
    ##        b = sum_{s != r mod 2}     (-1)^r*q^((s^2+r^2-1)/4)*zeta^r
    ## r, s run over ZZ but with opposite parities.
    ## For r=0, we need s odd, (s^2-1)/4 < precision, with s=2t+1 hence t^2+t < precision.
    ## For r=1, we need s even, s^2/4 < precision, with s=2t hence t^2 < precision.

    from sage.misc.all import xsrange

    a1 = -2*sum((2*t)**2 * q**(t**2)   for t in xsrange(1, precision) if t*t < precision)
    b1 = -2*sum(           q**(t**2)   for t in xsrange(1, precision) if t*t < precision)
    a0 = 2*sum((2*t+1)**2 * q**(t**2+t) for t in xsrange(precision) if t*t +t < precision)
    b0 = 2*sum(             q**(t**2+t) for t in xsrange(precision) if t*t +t < precision)
    b1 = b1 - 1
    
    ## print 'Done'
    ## Form A and B - the Jacobi forms used in [Sko]'s I map.

    (A0, A1, B0, B1) = (a0*etapow, a1*etapow, b0*etapow, b1*etapow)

    ## Calculate the image of the pair of modular forms (f,g) under
    ## [Sko]'s isomorphism I : M_{k} \oplus S_{k+2} -> J_{k,1}.
    
    fderiv = PS(q * f.qexp(precision).derivative())
    (f, g) = (PS(f.qexp(precision)), PS(g.qexp(precision)))

    ## Finally: I(f,g) is given by the formula below:
    Ifg0 = k/2*f*A0 - fderiv*B0 + g*B0 + O(q**precision)
    Ifg1 = k/2*f*A1 - fderiv*B1 + g*B1 + O(q**precision)

    ## For applying the Maass' lifting to genus 2 modular forms.
    ## we put the coefficients og Ifg into a dictionary Chi
    ## so that we can access the coefficient corresponding to 
    ## discriminant D by going Chi[D].

    ## Note: Ifg.list()[i].list()[j] gives the coefficient of q^i*zeta^j
    
    Cphi = {0: 0}  ## initialise dictionary. Value changed in the loop if we have a 'noncusp form'
    qcoeffs0 = Ifg0.list()
    qcoeffs1 = Ifg1.list()
    for i in xsrange(len(qcoeffs0)):
        Cphi[-4*i] = qcoeffs0[i]
        Cphi[1-4*i] = qcoeffs1[i]

    ## the most negative discriminant occurs when i is largest 
    ## and j is zero.  That is, discriminant 0^2-4*i
    ## Note that i < precision.
    maxD = -4*i

    """
    Create the Maass lift F := VI(f,g) as in [Sko].
    """
    
    ## The constant term is given by -Cphi[0]*B_{2k}/(4*k)
    ## (note in [Sko] this coeff has typos).
    ## For nonconstant terms,
    ## The Siegel coefficient of q^n * zeta^r * qdash^m is given 
    ## by the formula  \sum_{ a | gcd(n,r,m) } Cphi[D/a^2] where 
    ## D = r^2-4*n*m is the discriminant.  
    ## Hence in either case the coefficient 
    ## is fully deterimined by the pair (D,gcd(n,r,m)).
    ## Put (D,t) -> \sum_{ a | t } Cphi[D/a^2]
    ## in a dictionary (hash table) maassc.
    maassc = dict();
    ## First calculate maass coefficients corresponding to strictly positive definite matrices:
    from sage.rings.all import is_fundamental_discriminant
    from sage.misc.all import isqrt
    for disc in [d for d in xsrange(maxD, 0) if is_fundamental_discriminant(d)]:
        for s in xsrange(1, isqrt(maxD//disc)+1):
            ## add (disc*s^2,t) as a hash key, for each t that divides s
            for t in s.divisors():
                maassc[(disc*s**2, t)] = sum([a**(k-1)*Cphi[disc*s**2/a**2] for a in t.divisors()])

    ## Compute the coefficients of the Siegel form $F$:
    from sage.rings.all import gcd
    siegelq = dict();                
    if isinstance(prec, tuple):
        ## Note: m>=n>=r, n>=1 implies m>=n>r^2/4n
        for r in xsrange(0, bmax):
            for n in xsrange(max(r, 1), amax):
                for m in xsrange(n, cmax):
                    D=r**2 - 4*m*n
                    g=gcd([n, r, m])
                    siegelq[(n, r, m)] = maassc[(D, g)]
        bound = cmax
    else:
        bound = 0
        for n in xsrange(1, isqrt(Dtop//3)+1):
            for r in xsrange(n + 1):
                bound = max(bound, (Dtop + r*r)//(4*n) + 1)
                for m in xsrange(n, (Dtop + r*r)//(4*n) + 1):
                    D=r**2 - 4*m*n
                    g=gcd([n, r, m])
                    siegelq[(n, r, m)] = maassc[(D, g)]
    
    ## Secondly, deal with the singular part.
    ## Include the coeff corresponding to (0,0,0):
    ## maassc = {(0,0): -bernoulli(k)/(2*k)*Cphi[0]}
    from sage.rings.all import bernoulli
    siegelq[(0, 0, 0)] = -bernoulli(k)/(2*k)*Cphi[0]
    
    ## Calculate the other discriminant-zero maass coefficients:
    from sage.rings.all import sigma
    for i in xsrange(1, bound):
        ## maassc[(0,i)] = sigma(i, k-1) * Cphi[0]
        siegelq[(0, 0, i)] = sigma(i, k-1) * Cphi[0]

    return _SiegelModularForm_from_dict(group='Sp(4,Z)', weight=k, coeffs=siegelq, prec=clean_prec, name=name)


def _SiegelModularForm_from_file(loc):
    r"""
    Initialize an instance of SiegelModularForm_class from a file.

    EXAMPLE::

        sage: from sage.modular.siegel.siegel_modular_form import _SiegelModularForm_from_file
        sage: url = 'http://sage.math.washington.edu/home/nils/Siegel-Modular-Forms/data/upsilon-forms_20-32_1000/Upsilon_20.p'  # optional -- internet
        sage: H = _SiegelModularForm_from_file(url); H  # optional - internet
        Upsilon_20
        sage: H[(4, 4, 4)]                      # optional - internet
        248256200704
    """
    if 'http://' == loc[:7]:
        f = urllib.urlopen(loc)
    else:
        f = open(loc, 'r')
    data = cPickle.load(f)
    f.close()
    fmt, wt, name, pol, prec, coeffs = data
    if len(pol) > 2:
        from sage.rings.all import PolynomialRing, NumberField
        R = PolynomialRing(ZZ, 'x')
        K = NumberField(R(pol), name='a')
        for k in coeffs.iterkeys():
            coeffs[k] = K(coeffs[k])
    return _SiegelModularForm_from_dict(group='Sp(4,Z)', weight=wt, coeffs=coeffs, prec=prec, name=name)


def _SiegelModularForm_from_theta_characteristics(char_dict, prec=SMF_DEFAULT_PREC, name=None, hint=None):
    r"""
    Return the Siegel modular form
    `\sum_{l \in char_dict} \alpha[l] * \prod_i \theta_l[i](8Z)`,
    where `\theta_l[i]` denote the theta constant with characteristic `l`.   

    INPUT
        char_dict - a dictionary whose keys are *tuples* of theta characteristics
                    and whose values are in some ring.
        prec      - a precision for Siegel modular forms
        name      - a string describing this modular form

    EXAMPLES::

        sage: theta_constants = {((1, 1, 0, 0), (0, 0, 1, 1), (1, 1, 0, 0), (0, 0, 1, 1)): 1}
        sage: from sage.modular.siegel.siegel_modular_form import _SiegelModularForm_from_theta_characteristics
        sage: _SiegelModularForm_from_theta_characteristics(theta_constants, prec=100)
        Siegel modular form on None of weight 2
        sage: theta_constants = {((1, 1, 0, 0), (0, 0, 1, 1), (1, 1, 0, 0), (0, 0, 1, 1)): 1}
        sage: S = _SiegelModularForm_from_theta_characteristics(theta_constants, prec=100)
        sage: S.coeffs()[(2, 0, 10)]
        32

    TODO:
        Implement a parameter hint = "cusp_form" to prevent computing singular parts
    """
    coeffs = dict()
    smf_prec = SiegelModularFormPrecision(prec)
    for f in smf_prec:
        if hint is 'cusp_form':
            a, b, c = f
            if 0 == a: continue
        from theta_constant import _compute_theta_char_poly
        val = _compute_theta_char_poly(char_dict, f)
        if val != 0: coeffs[f] = val
    wt = 0
    for l in char_dict: wt = max(wt, len(l)/2) 
    return _SiegelModularForm_from_dict(group=None, weight=wt, coeffs=coeffs, prec=prec, name=name)


def _SiegelModularForm_from_QuadraticForm(Q, prec, name):
    """
    Return the theta series of degree 2 for the quadratic form Q.

    INPUT:

     - ``Q`` - a quadratic form.

     - ``prec`` - an integer.

    EXAMPLE::

        sage: from sage.modular.siegel.siegel_modular_form import _SiegelModularForm_from_QuadraticForm
        sage: Q11 = QuadraticForm(ZZ, 4, [3,0,11,0, 3,0,11, 11,0, 11])
        sage: _SiegelModularForm_from_QuadraticForm(Q11, 100, None)
        Siegel modular form on Gamma0(11) of weight 2
    """
    N = Q.level()
    k = Q.dim()/2
    coeffs = Q.theta_series_degree_2(prec)
    return _SiegelModularForm_from_dict(group='Gamma0(%d)'%N, weight=k, coeffs=coeffs, prec=prec, name=name)


def _SiegelModularForm_borcherds_lift(f, prec=SMF_DEFAULT_PREC, name=None):
    r"""
    Returns Borcherds lift.
 
    EXAMPLES::

        sage: from sage.modular.siegel.siegel_modular_form import _SiegelModularForm_borcherds_lift
        sage: _SiegelModularForm_borcherds_lift(0)
        Traceback (most recent call last):
        ...
        NotImplementedError: Borcherds lift not yet implemented
    """
    raise NotImplementedError, "Borcherds lift not yet implemented"


def _SiegelModularForm_yoshida_lift(f, g, prec=SMF_DEFAULT_PREC, name=None):
    r"""
    Returns the Yoshida lift.

    EXAMPLES::

        sage: from sage.modular.siegel.siegel_modular_form import _SiegelModularForm_yoshida_lift
        sage: _SiegelModularForm_yoshida_lift(0, 0)
        Traceback (most recent call last):
        ...
        NotImplementedError: Yoshida lift not yet implemented
    """
    raise NotImplementedError, 'Yoshida lift not yet implemented'


def _SiegelModularForm_from_weil_representation(gram, prec=SMF_DEFAULT_PREC, name=None):
    r"""
    Returns a form associated to a Weil representation.

    EXAMPLES::

        sage: from sage.modular.siegel.siegel_modular_form import _SiegelModularForm_from_weil_representation
        sage: _SiegelModularForm_from_weil_representation(matrix([[2, 1], [1, 1]]))
        Traceback (most recent call last):
        ...
        NotImplementedError: Weil representation argument not yet implemented    
    """
    raise NotImplementedError, 'Weil representation argument not yet implemented'


def _SiegelModularForm_singular_weight(gram, prec=SMF_DEFAULT_PREC, name=None):
    r"""
    Returns a singular modular form from gram.

    EXAMPLES::

        sage: from sage.modular.siegel.siegel_modular_form import _SiegelModularForm_singular_weight
        sage: _SiegelModularForm_singular_weight(matrix([[2, 1], [1, 1]]))
        Traceback (most recent call last):
        ...
        NotImplementedError: singular form argument not yet implemented
    """
    raise NotImplementedError, 'singular form argument not yet implemented'


def _SiegelModularForm_from_dict(group, weight, coeffs, prec, degree=2, parent=None, name=None):
    r"""
    Create an instance of SiegelModularForm_class(), where the parent
    is computed from the coeffs.

    EXAMPLES::

        sage: d = {(1, 1, 1): 1, (1, 0, 2): 2}
        sage: from sage.modular.siegel.siegel_modular_form import _SiegelModularForm_from_dict
        sage: S = _SiegelModularForm_from_dict('Sp(4,Z)', 2, d, 1); S
        Siegel modular form on Sp(4,Z) of weight 2
        sage: S.coeffs()
        {}
        sage: S = _SiegelModularForm_from_dict('Sp(4,Z)', 2, d, 4); S
        Siegel modular form on Sp(4,Z) of weight 2
        sage: S.coeffs()
        {(1, 1, 1): 1}
        sage: S = _SiegelModularForm_from_dict('Sp(4,Z)', 2, d, 9); S
        Siegel modular form on Sp(4,Z) of weight 2
        sage: S.coeffs()
        {(1, 1, 1): 1, (1, 0, 2): 2}
    """
    if prec is None:
        prec = -min([b**2 - 4*a*c for (a, b, c) in coeffs]) + 1

    smf_prec = SiegelModularFormPrecision(prec)
    coeffs_up_to_prec = {}
    for k in coeffs:
        if smf_prec.is_in_bound(k):
            coeffs_up_to_prec[k] = coeffs[k]                        
    if parent is None:
        from sage.structure.all import Sequence
        from siegel_modular_forms_algebra import SiegelModularFormsAlgebra
        s = Sequence(coeffs_up_to_prec.values())
        if weight % 2:
            weights = 'all'
        else:
            weights = 'even'
        if len(s) == 0:
            coeff_ring = ZZ
        else:
            coeff_ring = s.universe()
        parent = SiegelModularFormsAlgebra(coeff_ring=coeff_ring, group=group, weights=weights, degree=degree)
    
    return parent.element_class(parent=parent, weight=weight, coeffs=coeffs_up_to_prec, prec=prec, name=name)
                                                                    


def _normalized_prec(prec):
    r"""
    Return a normalized prec for instance of class SiegelModularForm_class.

    EXAMPLES::

        sage: from sage.modular.siegel.siegel_modular_form import _normalized_prec
        sage: _normalized_prec((5, 4, 5))
        (5, 4, 5)
        sage: _normalized_prec(101)
        101
        sage: _normalized_prec(infinity)
        +Infinity

    NOTE:
        $(a,b,c)$ is within the precison defined
        by prec, if $a,b,c < floor(a),floor(b),floor(c)%
        respectively $4ac-b^2 < floor(prec)$.
    """
    from sage.rings.all import infinity
    if prec is infinity: return prec
    from sage.functions.all import floor
    if isinstance(prec, tuple):
        a, b, c = prec
        a = min(floor(a), floor(c))
        b = min(floor(b), a)
        if b <= 0:
            return (0, 0, 0)
        return a, b, c
    D = max(0, floor(prec))
    while 0 != D%4 and 1 != D%4:
        D -= 1
    return D


def _is_bounded((a, b, c), prec):
    r"""
    Return True or False accordingly as (a, b, c) is in the
    region defined by prec (see. SiegelModularForm_class for what
    this means).

    EXAMPLES::

        sage: from sage.modular.siegel.siegel_modular_form import _is_bounded
        sage: _is_bounded((2, 0, 2), 15)
        False
        sage: _is_bounded((2, 0, 2), 16)
        False
        sage: _is_bounded((2, 0, 2), 17)
        True

    NOTE
        It is assumed that prec is normalized accordingly to the output
        of _normalized_prec().
    """
    if isinstance(prec, tuple):
        return (a, b, c) < prec        
    else:
        D = 4*a*c - b*b
        return D < prec if D != 0 else c < (prec+1)/4 



def _prec_min(prec1, prec2):
    r"""
    Returns the min of two precs.
    
    EXAMPLES::
  
        sage: from sage.modular.siegel.siegel_modular_form import _prec_min
        sage: _prec_min(100, 200)
        100
        sage: _prec_min((1, 0, 1), (2, 0, 1))
        (1, 0, 1)

    NOTE
       It is assumed that the arguments are normalized according to the output
       of _normalized_prec().
    """
    if type(prec1) != type(prec2):
        raise NotImplementedError, "addition for differing precs not yet implemented"
    return prec1 if prec1 <= prec2 else prec2 


