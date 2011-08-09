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

"""
Weight 2 Hilbert modular forms over F = Q(sqrt(5)).
"""

from sage.misc.cachefunc import cached_method

from sqrt5 import F, O_F

from psage.number_fields.sqrt5 import primes_of_bounded_norm


from sqrt5_fast import IcosiansModP1ModN
from sage.rings.all import is_Ideal, Integer, prime_divisors, QQ, next_prime, ZZ
from tables import ideals_of_norm
from sage.matrix.all import matrix
from sage.structure.all import Sequence

def ideal(X):
    if not is_Ideal(X):
        return O_F.ideal(X)
    return X

def next_prime_of_characteristic_coprime_to(P, I):
    p = next_prime(P.smallest_integer())
    N = ZZ(I.norm())
    while N%p == 0:
        p = next_prime(p)
    return F.primes_above(p)[0]

def next_prime_not_dividing(P, I):
    while True:
        p = P.smallest_integer()
        if p == 1:
            Q = F.ideal(2)
        elif p % 5 in [2,3]: # inert
            Q = F.primes_above(next_prime(p))[0]
        elif p == 5:
            Q = F.ideal(7)
        else: # p split
            A = F.primes_above(p)
            if A[0] == P:
                Q = A[1]
            else:
                Q = F.primes_above(next_prime(p))[0]
        if not Q.divides(I):
            return Q
        else:
            P = Q # try again

class Space(object):
    def __cmp__(self, right):
        if not isinstance(right, Space):
            raise NotImplementedError
        return cmp((self.level(), self.dimension(), self.vector_space()),
                   (right.level(), right.dimension(), right.vector_space()))

    def subspace(self, V):
        raise NotImplementedError
    
    def vector_space(self):
        raise NotImplementedError
    
    def basis(self):
        return self.vector_space().basis()

    def new_subspace(self, p=None):
        """
        Return (p-)new subspace of this space of Hilbert modular forms.

        WARNING: There are known examples where this is still wrong somehow...

        INPUT:
            - p -- None or a prime divisor of the level

        OUTPUT:
            - subspace of this space of Hilbert modular forms

        EXAMPLES::

        We make a space of level a product of 2 split primes and (2)::
        
            sage: from psage.modform.hilbert.sqrt5.hmf import F, HilbertModularForms
            sage: P = F.prime_above(31); Q = F.prime_above(11); R = F.prime_above(2)
            sage: H = HilbertModularForms(P*Q*R); H
            Hilbert modular forms of dimension 32, level 2*a-38 (of norm 1364=2^2*11*31) over QQ(sqrt(5))

        The full new space::
        
            sage: N = H.new_subspace(); N
            Subspace of dimension 22 of Hilbert modular forms of dimension 32, level 2*a-38 (of norm 1364=2^2*11*31) over QQ(sqrt(5))

        The new subspace for each prime divisor of the level::
        
            sage: N_P = H.new_subspace(P); N_P
            Subspace of dimension 31 of Hilbert modular forms of dimension 32, level 2*a-38 (of norm 1364=2^2*11*31) over QQ(sqrt(5))
            sage: N_Q = H.new_subspace(Q); N_Q
            Subspace of dimension 28 of Hilbert modular forms of dimension 32, level 2*a-38 (of norm 1364=2^2*11*31) over QQ(sqrt(5))
            sage: N_R = H.new_subspace(R); N_R
            Subspace of dimension 24 of Hilbert modular forms of dimension 32, level 2*a-38 (of norm 1364=2^2*11*31) over QQ(sqrt(5))
            sage: N_P.intersection(N_Q).intersection(N_R) == N
            True

        """
        V = self.degeneracy_matrix(p).kernel()
        return self.subspace(V)

    def decomposition(self, B, verbose=False):
        """
        Return Hecke decomposition of self using Hecke operators T_p
        coprime to the level with norm(p) <= B.
        """

        # TODO: rewrite to use primes_of_bounded_norm so that we
        # iterate through primes ordered by *norm*, which is
        # potentially vastly faster.  Delete these functions
        # involving characteristic!
        
        p = next_prime_of_characteristic_coprime_to(F.ideal(1), self.level())
        T = self.hecke_matrix(p)
        D = T.decomposition()
        while len([X for X in D if not X[1]]) > 0:
            p = next_prime_of_characteristic_coprime_to(p, self.level())
            if p.norm() > B:
                break
            if verbose: print p.norm()
            T = self.hecke_matrix(p)
            D2 = []
            for X in D:
                if X[1]:
                    D2.append(X)
                else:
                    if verbose: print T.restrict(X[0]).fcp()
                    for Z in T.decomposition_of_subspace(X[0]):
                        D2.append(Z)
            D = D2
        D = [self.subspace(X[0]) for X in D]
        D.sort()
        S = Sequence(D, immutable=True, cr=True, universe=int, check=False)
        return S

    def new_decomposition(self, verbose=False):
        """
        Return complete irreducible Hecke decomposition of new subspace of self.
        """
        V = self.degeneracy_matrix().kernel()
        p = next_prime_of_characteristic_coprime_to(F.ideal(1), self.level())
        T = self.hecke_matrix(p)
        D = T.decomposition_of_subspace(V)
        while len([X for X in D if not X[1]]) > 0:
            p = next_prime_of_characteristic_coprime_to(p, self.level())
            if verbose: print p.norm()
            T = self.hecke_matrix(p)
            D2 = []
            for X in D:
                if X[1]:
                    D2.append(X)
                else:
                    if verbose: print T.restrict(X[0]).fcp()
                    for Z in T.decomposition_of_subspace(X[0]):
                        D2.append(Z)
            D = D2
        D = [self.subspace(X[0]) for X in D]
        D.sort()
        S = Sequence(D, immutable=True, cr=True, universe=int, check=False)
        return S

class HilbertModularForms(Space):
    def __init__(self, level):
        """
        Space of Hilbert modular forms of weight (2,2) over Q(sqrt(5)).
        
        INPUT:
            - level -- an ideal or element of ZZ[(1+sqrt(5))/2].

        TESTS::

            sage: import psage
            sage: H = psage.hilbert.sqrt5.HilbertModularForms(3); H
            Hilbert modular forms of dimension 1, level 3 (of norm 9=3^2) over QQ(sqrt(5))
            sage: loads(dumps(H)) == H
            True
        """
        self._level = ideal(level)
        self._gen = self._level.gens_reduced()[0]
        self._icosians_mod_p1 = IcosiansModP1ModN(self._level)
        self._dimension = self._icosians_mod_p1.cardinality()
        self._vector_space = QQ**self._dimension
        self._hecke_matrices = {}
        self._degeneracy_matrices = {}

    def __repr__(self):
        return "Hilbert modular forms of dimension %s, level %s (of norm %s=%s) over QQ(sqrt(5))"%(
            self._dimension, str(self._gen).replace(' ',''),
            self._level.norm(), str(self._level.norm().factor()).replace(' ',''))

    def intersection(self, M):
        if isinstance(M, HilbertModularForms):
            assert self == M
            return self
        if isinstance(M, HilbertModularFormsSubspace):
            assert self == M.ambient()
            return M
        raise TypeError

    def level(self):
        return self._level
 
    def vector_space(self):
        return self._vector_space

    def weight(self):
        return (Integer(2),Integer(2))

    def dimension(self):
        return self._dimension

    def hecke_matrix(self, n):
        # I'm not using @cached_method, since I want to ensure that
        # the input "n" is properly normalized.  I also want it
        # to be transparent to see which matrices have been computed,
        # to clear the cache, etc.
        n = ideal(n)
        if self._hecke_matrices.has_key(n):
            return self._hecke_matrices[n]
        t = self._icosians_mod_p1.hecke_matrix(n)
        t.set_immutable()
        self._hecke_matrices[n] = t
        return t

    T = hecke_matrix

    def degeneracy_matrix(self, p=None):
        if self.level().is_prime():
            return matrix(QQ, self.dimension(), 0, sparse=True)
        if p is None:
            A = None
            for p in prime_divisors(self._level):
                if A is None:
                    A = self.degeneracy_matrix(p)
                else:
                    A = A.augment(self.degeneracy_matrix(p))
            return A
        p = ideal(p)
        if self._degeneracy_matrices.has_key(p):
            return self._degeneracy_matrices[p]
        d = self._icosians_mod_p1.degeneracy_matrix(p)
        d.set_immutable()
        self._degeneracy_matrices[p] = d
        return d
                
    def __cmp__(self, other):
        if not isinstance(other, HilbertModularForms):
            raise NotImplementedError
        # first sort by norms
        return cmp((self._level.norm(), self._level), (other._level.norm(), other._level))

    def subspace(self, V):
        return HilbertModularFormsSubspace(self, V)
    
    def elliptic_curve_factors(self):
        D = [X for X in self.new_decomposition() if X.dimension() == 1]
        # Have to get rid of the Eisenstein factor
        p = next_prime_of_characteristic_coprime_to(F.ideal(1), self.level())
        while True:
            q = p.residue_field().cardinality() + 1
            E = [A for A in D if A.hecke_matrix(p)[0,0] == q]
            if len(E) == 0:
                break
            elif len(E) == 1:
                D = [A for A in D if A != E[0]]
                break
            else:
                p = next_prime_of_characteristic_coprime_to(p, self.level())
        return Sequence([EllipticCurveFactor(X, number) for number, X in enumerate(D)],
                        immutable=True, cr=True, universe=int, check=False)

class HilbertModularFormsSubspace(Space):
    def __init__(self, H, V):
        assert H.dimension() == V.degree()
        self._H = H
        self._V = V

    def __repr__(self):
        return "Subspace of dimension %s of %s"%(self._V.dimension(), self._H)

    def subspace(self, V):
        raise NotImplementedError
        #return HilbertModularFormsSubspace(self._H, V)

    def intersection(self, M):
        if isinstance(M, HilbertModularForms):
            assert self.ambient() == M
            return self
        if isinstance(M, HilbertModularFormsSubspace):
            assert self.ambient() == M.ambient()
            H = self.ambient()
            V = self.vector_space().intersection(M.vector_space())
            return HilbertModularFormsSubspace(H, V)
        raise TypeError

    def ambient(self):
        return self._H

    def vector_space(self):
        return self._V

    def hecke_matrix(self, n):
        return self._H.hecke_matrix(n).restrict(self._V)
    T = hecke_matrix

    def degeneracy_matrix(self, p):
        return self._H.degeneracy_matrix(p).restrict_domain(self._V)

    def level(self):
        return self._H.level()

    def dimension(self):
        return self._V.dimension()
        

class EllipticCurveFactor(object):
    """
    A subspace of the new subspace of a space of weight 2 Hilbert
    modular forms that (conjecturally) corresponds to an elliptic
    curve.
    """
    def __init__(self, S, number):
        """
        INPUT:
            - S -- subspace of a space of Hilbert modular forms
            - ``number`` -- nonnegative integer indicating some
              ordering among the factors of a given level.
        """
        self._S = S
        self._number = number

    def __repr__(self):
        """
        EXAMPLES::

            sage: from psage.modform.hilbert.sqrt5.hmf import HilbertModularForms, F
            sage: H = HilbertModularForms(F.prime_above(31)).elliptic_curve_factors()[0]
            sage: type(H)
            <class 'psage.modform.hilbert.sqrt5.hmf.EllipticCurveFactor'>
            sage: H.__repr__()
            'Isogeny class of elliptic curves over QQ(sqrt(5)) attached to form number 0 in Hilbert modular forms of dimension 2, level 5*a-2 (of norm 31=31) over QQ(sqrt(5))'
        """
        return "Isogeny class of elliptic curves over QQ(sqrt(5)) attached to form number %s in %s"%(self._number, self._S.ambient())

    def base_field(self):
        """
        Return the base field of this elliptic curve factor.

        OUTPUT:
            - the field Q(sqrt(5))
        
        EXAMPLES::
        
            sage: from psage.modform.hilbert.sqrt5.hmf import HilbertModularForms, F
            sage: H = HilbertModularForms(F.prime_above(31)).elliptic_curve_factors()[0]
            sage: H.base_field()
            Number Field in a with defining polynomial x^2 - x - 1
        """
        return F

    def conductor(self):
        """
        Return the conductor of this elliptic curve factor, which is
        the level of the space of Hilbert modular forms.

        OUTPUT:
            - ideal of the ring of integers of Q(sqrt(5))
        
        EXAMPLES::
        
        """
        return self._S.level()

    def ap(self, P):
        """
        Return the trace of Frobenius at the prime P, for a prime P of
        good reduction.

        INPUT:
            - `P` -- a prime ideal of the ring of integers of Q(sqrt(5)).
            
        OUTPUT:
            - an integer

        EXAMPLES::

            sage: from psage.modform.hilbert.sqrt5.hmf import HilbertModularForms, F
            sage: H = HilbertModularForms(F.primes_above(31)[0]).elliptic_curve_factors()[0]
            sage: H.ap(F.primes_above(11)[0])
            4
            sage: H.ap(F.prime_above(5))
            -2
            sage: H.ap(F.prime_above(7))
            2

        We check that the ap we compute here match with those of a known elliptic curve
        of this conductor::

            sage: a = F.0; E = EllipticCurve(F, [1,a+1,a,a,0])
            sage: E.conductor().norm()
            31
            sage: 11+1 - E.change_ring(F.primes_above(11)[0].residue_field()).cardinality()
            4
            sage: 5+1 - E.change_ring(F.prime_above(5).residue_field()).cardinality()
            -2
            sage: 49+1 - E.change_ring(F.prime_above(7).residue_field()).cardinality()
            2
        """
        if P.divides(self.conductor()):
            if (P*P).divides(self.conductor()):
                # It is 0, because the reduction is additive.
                return ZZ(0)
            else:
                # TODO: It is +1 or -1, but I do not yet know how to
                # compute which without using the L-function.
                return '?'
        else:
            return self._S.hecke_matrix(P)[0,0]

    @cached_method
    def dual_eigenspace(self, B=None):
        """
        Return 1-dimensional subspace of the dual of the ambient space
        with the same system of eigenvalues as self.  This is useful when
        computing a large number of `a_P`.

        If we can't find such a subspace using Hecke operators of norm
        less than B, then we raise a RuntimeError.  This should only happen
        if you set B way too small, or self is actually not new.
        
        INPUT:
            - B -- Integer or None; if None, defaults to a heuristic bound.
        """
        N = self.conductor()
        H = self._S.ambient()
        V = H.vector_space()
        if B is None:
            # TODO: This is a heuristic guess at a "Sturm bound"; it's the same
            # formula as the one over QQ for Gamma_0(N).  I have no idea if this
            # is correct or not yet. It is probably much too large. -- William Stein
            from sage.modular.all import Gamma0
            B = Gamma0(N.norm()).index()//6 + 1
        for P in primes_of_bounded_norm(B+1):
            P = P.sage_ideal()
            if V.dimension() == 1:
                return V
            if not P.divides(N):
                T = H.hecke_matrix(P).transpose()
                V = (T - self.ap(P)).kernel_on(V)
        raise RuntimeError, "unable to isolate 1-dimensional space"

    @cached_method
    def dual_eigenvector(self, B=None):
        # 1. compute dual eigenspace
        E = self.dual_eigenspace(B)
        assert E.dimension() == 1
        # 2. compute normalized integer eigenvector in dual
        return E.basis_matrix()._clear_denom()[0][0]

    def aplist(self, B, dual_bound=None, algorithm='dual'):
        """
        Return list of traces of Frobenius for all primes P of norm
        less than bound.  Use the function
        psage.number_fields.sqrt5.primes_of_bounded_norm(B)
        to get the corresponding primes.

        INPUT:
            - `B` -- a nonnegative integer
            - ``dual_bound`` -- default None; passed to dual_eigenvector function
            - ``algorithm`` -- 'dual' (default) or 'direct'

        OUTPUT:
             - a list of Sage integers

        EXAMPLES::

        We compute the aplists up to B=50::

            sage: from psage.modform.hilbert.sqrt5.hmf import HilbertModularForms, F
            sage: H = HilbertModularForms(F.primes_above(71)[1]).elliptic_curve_factors()[0]
            sage: v = H.aplist(50); v
            [-1, 0, -2, 0, 0, 2, -4, 6, -6, 8, 2, 6, 12, -4]

        This agrees with what we get using an elliptic curve of this
        conductor::
        
            sage: a = F.0; E = EllipticCurve(F, [a,a+1,a,a,0])
            sage: from psage.ellcurve.lseries.aplist_sqrt5 import aplist
            sage: w = aplist(E, 50)
            sage: v == w
            True

        We compare the output from the two algorithms up to norm 75::

            sage: H.aplist(75, algorithm='direct') == H.aplist(75, algorithm='dual')
            True
        """
        primes = [P.sage_ideal() for P in primes_of_bounded_norm(B)]
        if algorithm == 'direct':
            return [self.ap(P) for P in primes]
        elif algorithm == 'dual':
            v = self.dual_eigenvector(dual_bound)
            i = v.nonzero_positions()[0]
            c = v[i]
            I = self._S.ambient()._icosians_mod_p1
            N = self.conductor()
            aplist = []
            for P in primes:
                if P.divides(N):
                    ap = self.ap(P)
                else:
                    ap = (I.hecke_operator_on_basis_element(P,i).dot_product(v))/c
                aplist.append(ap)
            return aplist
        else:
            raise ValueError, "unknown algorithm '%s'"%algorithm

