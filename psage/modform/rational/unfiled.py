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

Miscellaneous code related to rational modular forms, until I find a
better place to put it.

WARNING: Depend on this code at your own risk!  It will get moved to
another module at some point, probably.

"""
from __future__ import print_function
from __future__ import division


def modular_symbols_from_curve(C, N, num_factors=3):
    """
    Find the modular symbols spaces that shoudl correspond to the
    Jacobian of the given hyperelliptic curve, up to the number of
    factors we consider.
    
    INPUT:
        - C -- a hyperelliptic curve over QQ
        - N -- a positive integer
        - num_factors -- number of Euler factors to verify match up; this is
          important, because if, e.g., there is only one factor of degree g(C), we
          don't want to just immediately conclude that Jac(C) = A_f. 
        
    OUTPUT:
        - list of all sign 1 simple modular symbols factor of level N
          that correspond to a simple modular abelian A_f
          that is isogenous to Jac(C).  

    EXAMPLES::

        sage: from psage.modform.rational.unfiled import modular_symbols_from_curve

        sage: R.<x> = ZZ[]
        sage: f = x^7+4*x^6+5*x^5+x^4-3*x^3-2*x^2+1
        sage: C1 = HyperellipticCurve(f)
        sage: modular_symbols_from_curve(C1, 284)
        [Modular Symbols subspace of dimension 3 of Modular Symbols space of dimension 39 for Gamma_0(284) of weight 2 with sign 1 over Rational Field]
        
        sage: f = x^7-7*x^5-11*x^4+5*x^3+18*x^2+4*x-11
        sage: C2 = HyperellipticCurve(f)
        sage: modular_symbols_from_curve(C2, 284)
        [Modular Symbols subspace of dimension 3 of Modular Symbols space of dimension 39 for Gamma_0(284) of weight 2 with sign 1 over Rational Field]
    """
    # We will use the Eichler-Shimura relation and David Harvey's
    # p-adic point counting hypellfrob.  Harvey's code has the
    # constraint:   p > (2*g + 1)*(2*prec - 1).
    # So, find the smallest p not dividing N that satisfies the
    # above constraint, for our given choice of prec.

    f, f2 = C.hyperelliptic_polynomials()
    if f2 != 0:
        raise NotImplementedError("curve must be of the form y^2 = f(x)")
    if f.degree() % 2 == 0:
        raise NotImplementedError("curve must be of the form y^2 = f(x) with f(x) odd")

    prec = 1
    
    g = C.genus()
    B = (2*g + 1)*(2*prec - 1)

    from sage.rings.all import next_prime
    p = B

    # We use that if F(X) is the characteristic polynomial of the
    # Hecke operator T_p, then X^g*F(X+p/X) is the characteristic
    # polynomial of Frob_p, because T_p = Frob_p + p/Frob_p, according
    # to Eichler-Shimura.  Use this to narrow down the factors. 

    from sage.all import ModularSymbols, Integers, get_verbose
    D = ModularSymbols(N,sign=1).cuspidal_subspace().new_subspace().decomposition()
    D = [A for A in D if A.dimension() == g]

    from sage.schemes.hyperelliptic_curves.hypellfrob import hypellfrob

    while num_factors > 0:
        p = next_prime(p)
        while N % p == 0: p = next_prime(p)
        
        R = Integers(p**prec)['X']
        X = R.gen()
        D2 = []
        # Compute the charpoly of Frobenius using hypellfrob
        M = hypellfrob(p, 1, f)
        H = R(M.charpoly())
        
        for A in D:
            F = R(A.hecke_polynomial(p))
            # Compute charpoly of Frobenius from F(X)
            G = R(F.parent()(X**g * F(X + p/X)))
            if get_verbose(): print((p, G, H))
            if G == H:
                D2.append(A)
        D = D2
        num_factors -= 1

    return D
    
