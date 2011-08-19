#################################################################################
#
# (c) Copyright 2011 William Stein
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
Python code for very fast high precision computation of certain
specific modular forms of interest.
"""

import sys
from sage.rings.all import ZZ, QQ, PolynomialRing
from sage.modular.all import eisenstein_series_qexp, ModularForms
from sage.misc.all import cached_function, cputime, prod

from psage.modform.rational.special_fast import (
    _change_ring_ZZ, _evaluate_series_at_power_of_gen,
    Integer_list_to_polynomial)

def degen(h, n):
    """
    Return power series h(q^n) to the same precision as h.

    INPUT:
        - h -- power series over the rational numbers
        - n -- positive integer
    OUTPUT:
        - power series over the rational numbers

    EXAMPLES::

        sage: from psage.modform.rational.special import degen
        sage: R.<q> = QQ[[]]
        sage: f = 2/3 + 3*q + 14*q^2 + O(q^3)
        sage: degen(f,2)
        2/3 + 3*q^2 + O(q^3)
    """
    if n == 1:
        return h
    return _evaluate_series_at_power_of_gen(h,n,True)

#############################################################    


@cached_function
def cached_eisenstein_series_qexp(k, prec, verbose=False):
    """
    Return q-expansion of the weight k level 1 Eisenstein series to
    the requested precision.  The result is cached, so that subsequent
    calls are quick.

    INPUT:
        - k -- even positive integer
        - prec -- positive integer
        - verbose -- bool (default: False); if True, print timing information
    
    OUTPUT:
        - power series over the rational numbers

    EXAMPLES::

        sage: from psage.modform.rational.special import cached_eisenstein_series_qexp
        sage: cached_eisenstein_series_qexp(4, 10)
        1/240 + q + 9*q^2 + 28*q^3 + 73*q^4 + 126*q^5 + 252*q^6 + 344*q^7 + 585*q^8 + 757*q^9 + O(q^10)
        sage: cached_eisenstein_series_qexp(4, 5, verbose=True)
        Computing E_4(q) + O(q^5)... (time = ... seconds)
        1/240 + q + 9*q^2 + 28*q^3 + 73*q^4 + O(q^5)
        sage: cached_eisenstein_series_qexp(4, 5, verbose=True)  # cache used, so no timing printed
        1/240 + q + 9*q^2 + 28*q^3 + 73*q^4 + O(q^5)
    """
    if verbose: print "Computing E_%s(q) + O(q^%s)..."%(k,prec),; sys.stdout.flush(); t = cputime()
    e = eisenstein_series_qexp(k, prec)
    if verbose: print "(time = %.2f seconds)"%cputime(t)
    return e

def eis_qexp(k, t, prec, verbose=False):
    """
    Return the q-expansion of holomorphic level t weight k Eisenstein
    series.  Thus when k=2, this is E2(q) - t*E2(q^t), and is Ek(q^t)
    otherwise.

    INPUT:
        - k -- even positive integer
        - t -- positive integer
        - prec -- positive integer
        - verbose -- bool (default: False); if True, print timing information

    OUTPUT:
        - power series over the rational numbers

    EXAMPLES::

        sage: from psage.modform.rational.special import eis_qexp
        sage: eis_qexp(2, 1, 8)   # output is 0, since holomorphic
        O(q^8)
        sage: eis_qexp(2, 3, 8)
        1/12 + q + 3*q^2 + q^3 + 7*q^4 + 6*q^5 + 3*q^6 + 8*q^7 + O(q^8)
        sage: E2 = eisenstein_series_qexp(2, 8)
        sage: q = E2.parent().0
        sage: E2(q) - 3*E2(q^3)
        1/12 + q + 3*q^2 + q^3 + 7*q^4 + 6*q^5 + 3*q^6 + 8*q^7 + O(q^8)
        sage: eis_qexp(4, 3, 8)
        1/240 + q^3 + 9*q^6 + O(q^8)
        sage: eisenstein_series_qexp(4, 8)(q^3)
        1/240 + q^3 + 9*q^6 + 28*q^9 + 73*q^12 + 126*q^15 + 252*q^18 + 344*q^21 + O(q^24)

    Test verbose::

        sage: eis_qexp(2, 5, 7, verbose=True)
        Computing E_2(q) + O(q^8)... (time = ... seconds)
        1/6 + q + 3*q^2 + 4*q^3 + 7*q^4 + q^5 + 12*q^6 + O(q^7)
    """
    Ek = cached_eisenstein_series_qexp(k, prec+1, verbose=verbose)
    q = Ek.parent().gen()
    if k == 2:
        e = Ek - t*degen(Ek,t)
    else:
        e = degen(Ek,t)
    return e.add_bigoh(prec)

def eisenstein_gens(N, k, prec, verbose=False):
    r"""
    Find spanning list of 'easy' generators for the subspace of
    `M_k(\Gamma_0(N))` spanned by polynomials in Eisenstein series.

    The meaning of 'easy' is that the modular form is a monomial of
    weight `k` in the images of Eisenstein series from level `1` of
    weight dividing `k`.  Note that this list need not span the full
    space `M_k(\Gamma_0(N))`

    INPUT:
        - N -- positive integer
        - k -- even positive integer
        - prec -- positive integer
        - verbose -- bool (default: False); if True, print timing information

    OUTPUT:
        - list of triples (t,w,E), where E is the q-expansion of the
          weight w, holomorphic Eisenstein series of level t.

    EXAMPLES::

        sage: from psage.modform.rational.special import eisenstein_gens
        sage: eisenstein_gens(5,4,6)
        [(1, 2, O(q^6)), (5, 2, 1/6 + q + 3*q^2 + 4*q^3 + 7*q^4 + q^5 + O(q^6)), (1, 4, 1/240 + q + 9*q^2 + 28*q^3 + 73*q^4 + 126*q^5 + O(q^6)), (5, 4, 1/240 + q^5 + O(q^6))]
        sage: eisenstein_gens(11,2,6)
        [(1, 2, O(q^6)), (11, 2, 5/12 + q + 3*q^2 + 4*q^3 + 7*q^4 + 6*q^5 + O(q^6))]

    Test verbose option::

        sage: eisenstein_gens(11,2,9,verbose=True)
        Computing E_2(q) + O(q^10)... (time = 0.00 seconds)
        [(1, 2, O(q^9)), (11, 2, 5/12 + q + 3*q^2 + 4*q^3 + 7*q^4 + 6*q^5 + 12*q^6 + 8*q^7 + 15*q^8 + O(q^9))]
    """
    assert N > 1
    if k % 2 != 0:
        return []
    div = ZZ(N).divisors()
    gens = []
    for w in range(2, k+1, 2):
        d = div if w > 2 else div[1:]
        for t in div:
            gens.append((t, w, eis_qexp(w, t, prec, verbose=verbose)))
    return gens

class EisensteinMonomial(object):
    """
    A monomial in terms of Eisenstein series.
    """
    def __init__(self, v):
        """
        INPUT:
            - v -- a list of triples (k,t,n) that corresponds to E_k(q^t)^n.
        
        EXAMPLES::

            sage: from psage.modform.rational.special import EisensteinMonomial
            sage: e = EisensteinMonomial([(5,4,2), (5,6,3)]); e
            E4(q^5)^2*E6(q^5)^3
            sage: type(e)
            <class 'psage.modform.rational.special.EisensteinMonomial'>
        """
        self._v = v

    def __repr__(self):
        """
        EXAMPLES::
        
            sage: from psage.modform.rational.special import EisensteinMonomial
            sage: e = EisensteinMonomial([(5,4,12), (5,6,3)]); e.__repr__()
            'E4(q^5)^12*E6(q^5)^3'
        """
        return '*'.join(['E%s%s(q^%s)^%s'%(k,'^*' if k==2 else '', t,e) for t,k,e in self._v])
    
    def qexp(self, prec, verbose=False):
        """
        The q-expansion of a monomial in Eisenstein series.

        INPUT:
            - prec -- positive integer
            - verbose -- bool (default: False)

        EXAMPLES::

            sage: from psage.modform.rational.special import EisensteinMonomial
            sage: e = EisensteinMonomial([(5,4,2), (5,6,3)])
            sage: e.qexp(11)
            -1/7374186086400 + 43/307257753600*q^5 - 671/102419251200*q^10 + O(q^11)
            sage: E4 = eisenstein_series_qexp(4,11); q = E4.parent().gen()
            sage: E6 = eisenstein_series_qexp(6,11)
            sage: (E4(q^5)^2 * E6(q^5)^3).add_bigoh(11)
            -1/7374186086400 + 43/307257753600*q^5 - 671/102419251200*q^10 + O(q^11)        
        """
        z = [eis_qexp(k, t, prec, verbose=verbose)**e for t,k,e in self._v]
        if verbose: print "Arithmetic to compute %s +O(q^%s)"%(self, prec); sys.stdout.flush(); t=cputime()
        p = prod(z)
        if verbose: print "(time = %.2f seconds)"%cputime(t)
        return p


def _monomials(v, n, i):
    """
    Used internally for recursively computing monomials function.
    Returns each power of the ith generator times all products not
    involving the ith generator.

    INPUT:
        - `v` -- see docstring for monomials
        - `n` -- see docstring for monomials
        - `i` -- nonnegative integer

    OUTPUT:
        - list

    EXAMPLES::

        sage: from psage.modform.rational.special import _monomials
        sage: R.<x,y> = QQ[]
        sage: _monomials([(x,2),(y,3)], 6, 0)
        [y^2, x^3]
        sage: _monomials([(x,2),(y,3)], 6, 1)
        [x^3, y^2]
    """
    # each power of the ith generator times all products
    # not involving the ith generator.
    if len(v) == 1:
        b, k = v[0]
        if n%k == 0:
            return [b**(n//k)]
        else:
            return []
    else:
        z, k = v[i]
        w = list(v)
        del w[i]
        m = 0
        y = 1
        result = []
        while m <= n:
            for X in monomials(w, n - m):
                result.append(X*y)
            y *= z
            m += k
        return result

def monomials(v, n):
    """
    Return homogeneous monomials of degree exactly n.

    INPUT:
        - v -- list of pairs (x,d), where x is a ring element that
          we view as having degree d.
        - n -- positive integer

    OUTPUT:
        - list of monomials in elements of v

    EXAMPLES::

        sage: from psage.modform.rational.special import monomials
        sage: R.<x,y> = QQ[]
        sage: monomials([(x,2),(y,3)], 6)
        [y^2, x^3]
        sage: [monomials([(x,2),(y,3)], i) for i in [0..6]]
        [[1], [], [x], [y], [x^2], [x*y], [y^2, x^3]]
    """
    if len(v) == 0:
        return []
    return _monomials(v, n, 0)
        

def eisenstein_basis(N, k, verbose=False):
    r"""
    Find spanning list of 'easy' generators for the subspace of
    `M_k(\Gamma_0(N))` generated by level 1 Eisenstein series and
    their images of even integer weights up to `k`.

    INPUT:
        - N -- positive integer
        - k -- positive integer
        - ``verbose`` -- bool (default: False)

    OUTPUT:
        - list of monomials in images of level 1 Eisenstein series
        - prec of q-expansions needed to determine element of
          `M_k(\Gamma_0(N))`.

    EXAMPLES::

        sage: from psage.modform.rational.special import eisenstein_basis
        sage: eisenstein_basis(5,4)
        ([E4(q^5)^1, E4(q^1)^1, E2^*(q^5)^2], 3)
        sage: eisenstein_basis(11,2,verbose=True)  # warning below because of verbose
        Warning -- not enough series.
        ([E2^*(q^11)^1], 2)
        sage: eisenstein_basis(11,2,verbose=False)
        ([E2^*(q^11)^1], 2)
    """
    assert N > 1
    if k % 2 != 0:
        return []
    # Make list E of Eisenstein series, to enough precision to
    # determine them, until we span space.
    M = ModularForms(N, k)
    prec = M.echelon_basis()[-1].valuation() + 1
    
    gens = eisenstein_gens(N, k, prec)
    R = PolynomialRing(ZZ, len(gens), ['E%sq%s'%(g[1],g[0]) for g in gens])
    z = [(R.gen(i), g[1]) for i, g in enumerate(gens)]
    m = monomials(z, k)
    
    A = QQ**prec
    V = A.zero_subspace()
    E = []
    for i, z in enumerate(m):
        d = z.degrees()
        f = prod(g[2]**d[i] for i, g in enumerate(gens) if d[i])
        v = A(f.padded_list(prec))
        if v not in V:
            V = V + A.span([v])
            w = [(gens[i][0],gens[i][1],d[i]) for i in range(len(d)) if d[i]]
            E.append(EisensteinMonomial(w))
            if V.dimension() == M.dimension():
                 return E, prec

    if verbose: print "Warning -- not enough series."
    return E, prec
    

class FastModularForm(object):
    """
    EXAMPLES::

        sage: from psage.modform.rational.special import FastModularForm

    Level 5, weight 4::
    
        sage: f = FastModularForm(CuspForms(5,4).0); f
        (-250/3)*E4(q^5)^1 + (-10/3)*E4(q^1)^1 + (13)*E2^*(q^5)^2
        sage: f.qexp(10)
        q - 4*q^2 + 2*q^3 + 8*q^4 - 5*q^5 - 8*q^6 + 6*q^7 - 23*q^9 + O(q^10)
        sage: parent(f.qexp(10))
        Power Series Ring in q over Integer Ring
        sage: t=cputime(); h=f.qexp(10^5); assert cputime(t)<5   # a little timing test.

    Level 5, weight 6::

        sage: f = CuspForms(5,6).0; g = FastModularForm(f); g
        (521/6)*E6(q^5)^1 + (-1/30)*E6(q^1)^1 + (248)*E2^*(q^5)^1*E4(q^5)^1

    Level 8, weight 4::

        sage: f = CuspForms(8,4).0; g = FastModularForm(f); g
        (-256/3)*E4(q^8)^1 + (4)*E4(q^4)^1 + (1)*E4(q^2)^1 + (-4/3)*E4(q^1)^1 + (4)*E2^*(q^8)^2

    Level 7, weight 4::

        sage: f = CuspForms(7,4).0; g = FastModularForm(f); g
        (-147/2)*E4(q^7)^1 + (-3/2)*E4(q^1)^1 + (5)*E2^*(q^7)^2

    """
    def __init__(self, f, verbose=False):
        """
        EXAMPLES::

        Level 3, weight 6::

            sage: from psage.modform.rational.special import FastModularForm
            sage: f = CuspForms(3,6).0; g = FastModularForm(f); g
            (549/10)*E6(q^3)^1 + (-3/10)*E6(q^1)^1 + (312)*E2^*(q^3)^1*E4(q^3)^1
            sage: type(g)
            <class 'psage.modform.rational.special.FastModularForm'>            
        """
        import sage.modular.modform.element
        if not isinstance(f, sage.modular.modform.element.ModularForm_abstract):
            raise TypeError
        chi = f.character()
        if not chi or not chi.is_trivial():
            raise ValueError, "form must trivial character"
        self._f = f
        
        N = f.level()
        k = f.weight()
        B, prec = eisenstein_basis(N, k, verbose=verbose)

        # Now write f in terms of q-expansions of elements in B.
        V = QQ**prec
        W = V.span_of_basis([V(h.qexp(prec).padded_list()) for h in B])
        self._basis = B
        self._coordinates = W.coordinates(f.qexp(prec).padded_list())
        self._verbose = verbose

        assert self.qexp(prec) == f.qexp(prec), "bug -- q-expansions don't match"
        
    def qexp(self, prec):
        """
        Return the q-expansion of this fast modular form to the given
        precision.

        EXAMPLES::

            sage: from psage.modform.rational.special import FastModularForm
            sage: f = FastModularForm(CuspForms(5,4).0); f
            (-250/3)*E4(q^5)^1 + (-10/3)*E4(q^1)^1 + (13)*E2^*(q^5)^2
            sage: f.qexp(10)
            q - 4*q^2 + 2*q^3 + 8*q^4 - 5*q^5 - 8*q^6 + 6*q^7 - 23*q^9 + O(q^10)
            sage: g = f.modular_form(); g
            q - 4*q^2 + 2*q^3 + 8*q^4 - 5*q^5 + O(q^6)
            sage: g.qexp(10)
            q - 4*q^2 + 2*q^3 + 8*q^4 - 5*q^5 - 8*q^6 + 6*q^7 - 23*q^9 + O(q^10)        
        """
        f = sum(c*self._basis[i].qexp(prec, verbose=self._verbose)
                for i, c in enumerate(self._coordinates) if c)
        return ZZ[['q']](_change_ring_ZZ(f.polynomial()), prec)

    q_expansion = qexp
        
    def modular_form(self):
        """
        Return the underlying modular form object.

        EXAMPLES::

            sage: from psage.modform.rational.special import FastModularForm
            sage: f = FastModularForm(CuspForms(5,4).0); f
            (-250/3)*E4(q^5)^1 + (-10/3)*E4(q^1)^1 + (13)*E2^*(q^5)^2
            sage: f.qexp(10)
            q - 4*q^2 + 2*q^3 + 8*q^4 - 5*q^5 - 8*q^6 + 6*q^7 - 23*q^9 + O(q^10)
            sage: g = f.modular_form(); g.qexp(10)
            q - 4*q^2 + 2*q^3 + 8*q^4 - 5*q^5 - 8*q^6 + 6*q^7 - 23*q^9 + O(q^10)
        """
        return self._f

    def __repr__(self):
        """
        Print representation.

        EXAMPLES::

            sage: from psage.modform.rational.special import FastModularForm
            sage: f = FastModularForm(CuspForms(5,4).0)
            sage: f.__repr__()
            '(-250/3)*E4(q^5)^1 + (-10/3)*E4(q^1)^1 + (13)*E2^*(q^5)^2'
        """
        return ' + '.join('(%s)*%s'%(c, self._basis[i]) for i, c in enumerate(self._coordinates) if c)






#####################################################
# CM forms    
#####################################################

from sage.all import fast_callable, prime_range, SR
from psage.ellcurve.lseries.helper import extend_multiplicatively_generic

def elliptic_cm_form(E, n, prec, aplist_only=False, anlist_only=False):
    """
    Return q-expansion of the CM modular form associated to the n-th
    power of the Grossencharacter associated to the elliptic curve E.

    INPUT:
        - E -- CM elliptic curve
        - n -- positive integer
        - prec -- positive integer
        - aplist_only -- return list only of ap for p prime
        - anlist_only -- return list only of an

    OUTPUT:
        - power series with integer coefficients

    EXAMPLES::

        sage: from psage.modform.rational.special import elliptic_cm_form
        sage: f = CuspForms(121,4).newforms(names='a')[0]; f
        q + 8*q^3 - 8*q^4 + 18*q^5 + O(q^6)
        sage: E = EllipticCurve('121b')
        sage: elliptic_cm_form(E, 3, 7)
        q + 8*q^3 - 8*q^4 + 18*q^5 + O(q^7)
        sage: g = elliptic_cm_form(E, 3, 100)
        sage: g == f.q_expansion(100)
        True
    """
    if not E.has_cm():
        raise ValueError, "E must have CM"
    n = ZZ(n)
    if n <= 0:
        raise ValueError, "n must be positive"

    prec = ZZ(prec)
    if prec <= 0:
        return []
    elif prec <= 1:
        return [ZZ(0)]
    elif prec <= 2:
        return [ZZ(0), ZZ(1)]
    
    # Derive formula for sum of n-th powers of roots
    a,p,T = SR.var('a,p,T')
    roots = (T**2 - a*T + p).roots(multiplicities=False)
    s = sum(alpha**n for alpha in roots).simplify_full()
    
    # Create fast callable expression from formula
    g = fast_callable(s.polynomial(ZZ))

    # Compute aplist for the curve
    v = E.aplist(prec)

    # Use aplist to compute ap values for the CM form attached to n-th
    # power of Grossencharacter.
    P = prime_range(prec)

    if aplist_only:
        # case when we only want the a_p (maybe for computing an
        # L-series via Euler product)
        return [g(ap,p) for ap,p in zip(v,P)]

    # Default cause where we want all a_n
    anlist = [ZZ(0),ZZ(1)] + [None]*(prec-2)
    for ap,p in zip(v,P):
        anlist[p] = g(ap,p)

    # Fill in the prime power a_{p^r} for r >= 2.
    N = E.conductor()
    for p in P:
        prm2 = 1
        prm1 = p
        pr = p*p
        pn = p**n
        e = 1 if N%p else 0
        while pr < prec:
            anlist[pr] = anlist[prm1] * anlist[p]
            if e:
                anlist[pr] -= pn * anlist[prm2]
            prm2 = prm1
            prm1 = pr
            pr *= p

    # fill in a_n with n divisible by at least 2 primes
    extend_multiplicatively_generic(anlist)

    if anlist_only:
        return anlist

    f = Integer_list_to_polynomial(anlist, 'q')
    return ZZ[['q']](f, prec=prec)

    
    
                    
