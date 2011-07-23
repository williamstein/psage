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

from sqrt5 import F, O_F, primes_of_bounded_norm

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
    def _subspace(self, V):
        raise NotImplementedError
    def vector_space(self):
        raise NotImplementedError
    def basis(self):
        return self.vector_space().basis()

    def new_subspace(self, p=None):
        V = self.degeneracy_matrix(p).kernel()
        return self._subspace(V)

    def decomposition(self, B, verbose=False):
        """
        Return Hecke decomposition of self using Hecke operators T_p coprime to the level with norm(p) <= B.
        """
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
        D = [self._subspace(X[0]) for X in D]
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
        D = [self._subspace(X[0]) for X in D]
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
        if ZZ(n.norm()).gcd(ZZ(self.level().norm())) != 1:
            raise NotImplementedError, "norm of n must currently be coprime to the norm of the level"
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

    def _subspace(self, V):
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
        return Sequence([ModularEllipticCurve(X, number) for number, X in enumerate(D)],
                        immutable=True, cr=True, universe=int, check=False)

class HilbertModularFormsSubspace(Space):
    def __init__(self, H, V):
        self._H = H
        self._V = V

    def __repr__(self):
        return "Subspace of dimension %s of %s"%(self._V.dimension(), self._H)

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
        

class ModularEllipticCurve(object):
    def __init__(self, S, number):
        self._S = S
        self._number = number

    def __repr__(self):
        return "Isogeny class of elliptic curves over QQ(sqrt(5)) attached to form number %s in %s"%(self._number, self._S.ambient())

    def conductor(self):
        return self._S.level()

    def ap(self, P):
        if P.smallest_integer().gcd(ZZ(self.conductor().norm())) != 1:
            return '?'
        else:
            return self._S.hecke_matrix(P)[0,0]

    def aplist(self, B):
        return [self.ap(P) for P in primes_of_bounded_norm(F, B)]
