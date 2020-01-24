from __future__ import print_function
from __future__ import absolute_import
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

from builtins import str
from builtins import range
_multiprocess_can_split_ = True

import random
from sage.all import set_random_seed, primes
from .sqrt5 import *
from .sqrt5_fast import *

def prime_powers(B=100, emax=5):
    for p in primes(B+1):
        for P in F.primes_above(p):
            for e in range(1,emax+1):
                yield p, P, e

################################################################
# Tests of Residue Rings
################################################################

def residue_rings(B=100, emax=5):
    for p, P, e in prime_powers(B, emax):
        yield ResidueRing(P, e)

def residue_ring_creation(B=100, emax=5):
    alpha = -F.gen() - 3
    for R in residue_rings(B, emax):
        R.residue_field()
        R._moduli()
        zero = R(0); one = R(1)
        assert R(alpha) - R(alpha) == zero
        assert one - one == zero
        str(R)
        assert R[0] == zero
        assert R[1] == one

def test_residue_ring_1():
    residue_ring_creation(B=100, emax=4)

def test_residue_ring_2():    
    residue_ring_creation(B=35, emax=6)                

def test_residue_ring_3():    
    residue_ring_creation(B=200, emax=2)

def residue_ring_arith(P, e, n=None):
    """
    Test that arithmetic in the residue field is correct.
    """
    set_random_seed(0)
    R = ResidueRing(P, e)
    n0, n1 = R._moduli()
    alpha = F.gen()
    reps = [i+j*alpha for i in range(n0) for j in range(n1)]
    if n is not None:
        reps = [reps[random.randrange(len(reps))] for i in range(n)]
    reps_mod = [R(a) for a in reps]
    for i in range(len(reps)):
        for j in range(len(reps)):
            assert reps_mod[i] * reps_mod[j] == R(reps[i] * reps[j])
            assert reps_mod[i] + reps_mod[j] == R(reps[i] + reps[j])
            assert reps_mod[i] - reps_mod[j] == R(reps[i] - reps[j])

def test_residue_ring_arith2():
    for i in range(1,4):
        residue_ring_arith(F.primes_above(2)[0], i)

def test_residue_ring_arith3():
    for i in range(1,4):
        residue_ring_arith(F.primes_above(2)[0], i, 50)

def test_residue_ring_arith5():
    for i in range(1,4):
        residue_ring_arith(F.primes_above(5)[0], i, 50)

def test_residue_ring_arith11():
    for i in range(1,4):
        residue_ring_arith(F.primes_above(11)[0], i, 50)
        residue_ring_arith(F.primes_above(11)[1], i, 50)

################################################################
# Testing pow function on elements
################################################################
def test_pow():
    def f(B,e):
        for R in residue_rings(B,e):
            for a in R:
                b = a*a
                for i in range(4):
                    assert b == a**(i+2)
                    if a.is_unit():
                        assert ~b == a**(-(i+2))
                    b *= a
    f(100,1)                    
    f(11,3)

################################################################
# Testing is_square in case of degree 1
################################################################
def is_square_0(B=11, emax=2):
    for R in residue_rings(B, emax):
        assert set([x*x for x in R]) == set([x for x in R if x.is_square()])

def test_is_square_0():
    is_square_0(20,2)
    is_square_0(12,3)

def test_is_square_2(emax=6):
    for R in residue_rings(2, emax):
        assert set([x*x for x in R]) == set([x for x in R if x.is_square()])

def test_is_square_5(emax=5):
    P = F.primes_above(5)[0]
    for e in range(1,emax+1):
        R = ResidueRing(P, e)
        assert set([x*x for x in R]) == set([x for x in R if x.is_square()])

################################################################
# Testing sqrt function on elements
################################################################
def test_sqrt1(B=11,emax=2):
    for R in residue_rings(B,emax):
        for a in R:
            b = a*a
            assert b.is_square()
            assert b.sqrt()**2 == b

def test_sqrt2(B=11,emax=2):
    for R in residue_rings(B,emax):
        for a in R:
            if a.is_square():
                assert a.sqrt()**2 == a

def test_sqrt5(emax=4):
    P = F.primes_above(5)[0]
    for e in range(1,emax+1):
        for a in ResidueRing(P, e):
            if a.is_square():
                assert a.sqrt()**2 == a

def test_sqrt5b(emax=4):
    P = F.primes_above(5)[0]
    for e in range(1,emax+1):
        if e % 2: continue
        for a in ResidueRing(P,e):
            b = a*a
            assert b.is_square()
            assert b.sqrt()**2 == b


################################################################
# Test computing local splitting map
################################################################


def test_local_splitting_maps(B=11,e=3, verbose=True):
    ans = []
    for p, P, e in prime_powers(B,e):
        print(p, P, e)
        r = ModN_Reduction(P**e)

def test_local_splitting_maps_1():
    test_local_splitting_maps(100,2)

def test_local_splitting_maps_2():
    test_local_splitting_maps(11,5)

################################################################
# Test computing R^* \ P^1
################################################################


def test_icosians_mod_p1(B=11, e=3):
    for p, P, e in prime_powers(B,e):
        H = IcosiansModP1ModN(P**e)


################################################################
# Testing mod_long functions
################################################################

def test_sqrtmod_long_2(e_max=15):
    for e in range(1, e_max+1):
        R = Integers(2**e)
        X = set([a*a for a in R])
        for b in X:
            c = sqrtmod_long(int(b), 2, e)
            assert c*c == b
    


################################################################
# Test Hecke operators commute
################################################################

def test_hecke_commutes(Bmin=2, Bmax=50):
    from .tables import ideals_of_bounded_norm
    from .hmf import HilbertModularForms, next_prime_of_characteristic_coprime_to
    for N in ideals_of_bounded_norm(Bmax):
        if N.norm() < Bmin: continue
        print(N.norm())
        H = HilbertModularForms(N)
        p = next_prime_of_characteristic_coprime_to(F.ideal(1), N)
        q = next_prime_of_characteristic_coprime_to(p, N)
        print(p, q)
        T_p = H.T(p)
        T_q = H.T(q)
        assert T_p*T_q == T_q*T_p, "Hecke operators T_{%s} and T_{%s} at level %s (of norm %s) don't commute"%(p, q, N, N.norm())
    



################################################################
# New subspaces
################################################################

## def test_compute_new_subspace(Bmin=2, Bmax=200, try_hecke=False):
##     from tables import ideals_of_bounded_norm
##     from hmf import HilbertModularForms, next_prime_of_characteristic_coprime_to
##     for N in ideals_of_bounded_norm(Bmax):
##         if N.norm() < Bmin: continue
##         print N.norm()
##         H = HilbertModularForms(N)
##         NS = H.new_subspace()
##         if try_hecke:
##             p = next_prime_of_characteristic_coprime_to(F.ideal(1), N)
##             NS.hecke_matrix(p)


def test_hecke_invariance_of_new_subspace(Bmin=2, Bmax=300):
    from .tables import ideals_of_bounded_norm
    from .hmf import HilbertModularForms, next_prime_of_characteristic_coprime_to
    for N in ideals_of_bounded_norm(Bmax):
        if N.norm() < Bmin: continue
        print(N.norm())
        H = HilbertModularForms(N)
        NS = H.new_subspace()
        p = next_prime_of_characteristic_coprime_to(F.ideal(1), N)
        q = next_prime_of_characteristic_coprime_to(p, N)
        print(p, q)
        T_p = NS.T(p)
        T_q = NS.T(q)
        assert T_p*T_q == T_q*T_p, "Hecke operators T_{%s} and T_{%s} at level %s (of norm %s) don't commute"%(p, q, N, N.norm())
 
