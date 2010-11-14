# -*- coding: utf-8 -*-
r"""
p-adic L-series of modular Jacobians with ordinary reduction at p


REFERENCES:
        
- [MTT] B. Mazur, J. Tate, and J. Teitelbaum,
  On `p`-adic analogues of the conjectures of Birch and 
  Swinnerton-Dyer, Inventiones mathematicae 84, (1986), 1-48.
   
- [SW] William Stein and Christian Wuthrich, Computations About Tate-Shafarevich Groups
  using Iwasawa theory, preprint 2009.
      
AUTHORS:

- William Stein and Jennifer Balakrishnan (2010-07-01): first version

"""

######################################################################
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
######################################################################

	
from sage.rings.integer_ring import   ZZ
from sage.rings.rational_field import QQ
from sage.rings.padics.factory import Qp, Zp
from sage.rings.infinity import infinity
from sage.rings.all import LaurentSeriesRing, PowerSeriesRing, PolynomialRing, Integers

from sage.rings.integer import Integer
from sage.rings.arith import valuation, binomial, kronecker_symbol, gcd, prime_divisors, valuation

from sage.structure.sage_object import SageObject

from sage.misc.all import verbose, denominator, get_verbose
from sage.databases.cremona import parse_cremona_label
from sage.schemes.elliptic_curves.constructor import EllipticCurve
import sage.rings.arith as arith

from sage.modules.free_module_element import vector
import sage.matrix.all as matrix

# from sage.interfaces.all import gp
from sage.misc.functional import log

from sage.modular.modsym.modsym import ModularSymbols

class pAdicLseries(SageObject):
    r"""
    The `p`-adic L-series of a modular Jacobian

    EXAMPLES:
    An ordinary example::
    

    """
    def __init__(self, J, p, normalize='L_ratio'):
        r"""
        INPUT:
        
        -  ``J`` - modular abelian variety
        -  ``p`` - a prime of good reduction
        -  ``normalize`` - ``'L_ratio'`` (default), ``'period'`` or ``'none'``;
           this is describes the way the modular symbols
           are normalized 
        """
        self._J = J
        self._level = J.level()
        self._p = ZZ(p)
        self._normalize = normalize
        if not self._p.is_prime():
            raise ValueError, "p (=%s) must be a prime"%p
        A = self._J.modular_symbols(sign=1)
        self._modular_symbols_subspace = A   
        v = A.dual_eigenvector()
        self._dual_eigenvector = v
        self._hecke_eigenvalue_field = v[0].parent()
        
    def __cmp__(self,other):
        r"""
        Compare self and other. 
        
        TESTS::
        sage: lp1 = J0(23)[0].padic_lseries(5)
	    sage: lp2 = J0(23)[0].padic_lseries(7)
	    sage: lp3 = J0(29)[0].padic_lseries(5)
	    sage: lp1 == lp1
	    True
	    sage: lp1 == lp2
	    False
	    sage: lp1 == lp3
	    False

        """
        c = cmp(type(self), type(other))
        if c: 
            return c
        return cmp((self._J, self._p), (other._J, other._p))               

    def modular_symbols_subspace(self):
        """
	"""
        return self._modular_symbols_subspace
    
    def hecke_eigenvalue_field(self):
        """
	The field of definition of the dual eigenform.
	"""
        return self._hecke_eigenvalue_field

    def psi(self):
    	"""
        The embedding $\Q(\alpha) \into \Q_p(a)$ sending $\alpha \mapsto a$.
	"""
        K_f = self._hecke_eigenvalue_field
        p = self._p
#        kbar = K_f.residue_field(p)
        Q = Qp(p)

        ###split case
        pOK = K_f.factor(p)
        if (len(pOK) == 2 and pOK[0][1] == 1):
            R = Q['x']
            r1, r2 = R(K_f.defining_polynomial()).roots()
            psi1 = K_f.hom([r1[0]])
            psi2 = K_f.hom([r2[0]])
            return [psi1, psi2]
        else:
            F = Q.extension(K_f.defining_polynomial(),names='a')
            a = F.gen()
            psi = self._psis = [K_f.hom([a])]
            return psi
	
    def modular_symbol(self,r):
        """
        Compute the modular symbol at the cusp $r$.
	"""
        v = self._dual_eigenvector

        try:
            psis = self._psis
        except AttributeError:
            psis = self._psis = self.psi()
                                        
        # TODO: rewrite this function to be a separate Cython class
        # that just does reducing [r,infinity] to a rational using
        # exactly the data that ones needs to do this (not going
        # through modular symbols), and it'll make this massively
        # faster.
        if len(psis) == 1:
            psi = psis[0]
            M = self.modular_symbols_subspace().ambient()
            s = M([r,infinity])
            return psi(v.dot_product(s.element()))
        else:
	    raise NotImplementedError

#############stuff for p-adic BSD, should be moved######

    def regulator(self):
        """
        """
        list = [23, 29, 31, 35, 39, 63,65, 87, 117, 135, 175, 189]
        if (self.abelian_variety().level() in list):
            return 1
        else:
            raise NotImplementedError

    def rank(self):
        """
        """
        list = [23, 29, 31, 35,39, 63,65, 87, 117, 135, 175, 189]
        if (self.abelian_variety().level() in list):
            return 1
        else:
            raise NotImplementedError

    def tamagawa_prod(self):
        """
        """
        list = [23, 29, 31, 35,39, 63,65, 87, 117, 125, 133, 135, 175, 189]
        A = self.abelian_variety()
        if A.dimension() != 2:
            raise NotImplementedError
        lev = self._level
        if lev not in list:
            raise NotImplementedError
        elif lev == 23:
            return 11
        elif lev == 29:
            return 7
        elif lev == 31:
            return 5
        elif lev == 35:
            return 32
        elif lev == 39:
            return 28
        elif lev == 63:
            return 6
        elif self.label() == '65a(1,65)':
            return 7
        elif self.label() == '65b(1,65)':
            return 3
        elif lev == 87:
            return 5
        elif self.label() == '117a(1,117)':
            return 12
        elif self.label() == '117b(1,117)':
            return 4
        elif self.label() == '125b(1,125)':
            return 5
        elif self.label() == '133a(1,133)':
            return 5
        elif lev == 135:
            return 3
        elif lev == 175:
            return 5
        elif lev == 189:
            return 3

    def torsion_order(self):
        """
        """
        list = [23, 29, 31, 35,39, 63,65, 87, 117, 125, 133, 135, 175, 189]
        A = self.abelian_variety()
        if A.dimension() != 2:
            raise NotImplementedError
        lev = self._level
        if lev not in list:
            raise NotImplementedError
        elif lev == 23:
            return 11
        elif lev == 29:
            return 7
        elif lev == 31:
            return 5
        elif lev == 35:
            return 16
        elif lev == 39:
            return 28
        elif lev == 63:
            return 6
        elif self.label() == '65a(1,65)':
            return 14
        elif self.label() == '65b(1,65)':
            return 6
        elif lev == 87:
            return 5
        elif self.label() == '117a(1,117)':
            return 6
        elif self.label() == '117b(1,117)':
            return 2
        elif self.label() == '125b(1,125)':
            return 5
        elif self.label() == '133a(1,133)':
            return 5
        elif lev == 135:
            return 3
        elif lev == 175:
            return 5
        elif lev == 189:
            return 3

    def sha(self):
    	"""
	"""
        list = [23, 29, 31, 35,39, 63,65, 87, 117, 125, 133, 135, 175, 189]
        A = self.abelian_variety()
        if A.dimension() != 2:
            raise NotImplementedError
        lev = self._level
        if lev not in list:
            raise NotImplementedError
        elif lev == 23:
            return 1
        elif lev == 29:
            return 1
        elif lev == 31:
            return 1
        elif lev == 35:
            return 1
        elif lev == 39:
            return 1
        elif lev == 63:
            return 1
	elif self.label() == '65a(1,65)':
            return 2
	elif self.label() == '65b(1,65)':
            return 2
	elif lev == 87:
            return 1
	elif self.label() == '117a(1,117)':
            return 1
	elif self.label() == '117b(1,117)':
	    return 1
	elif self.label() == '125b(1,125)':
            return 4
	elif self.label() == '133a(1,133)':
            return 4
	elif lev == 135:
            return 1
	elif lev == 175:
            return 1
	elif lev == 189:
            return 1

    def rhs(self):
        """
	"""
        list = [23, 29, 31, 35,39, 63,65, 87, 117, 125, 133, 135, 175, 189]
        lev = self.abelian_variety().level()
        if lev not in list:
            raise NotImplementedError
        else:
            try:
   	        eps = (1-1/self.alpha()).norm()**2
            except AttributeError:
                eps = (1-1/self.alpha())**4
	    return eps*(self.tamagawa_prod()*self.sha())/(self.torsion_order()**2)

        




                

#######################################################

    def abelian_variety(self):
        r"""
        Return the abelian variety to which this `p`-adic L-series is associated.
        
        EXAMPLES::
        
	    sage: L = J0(23)[0].padic_lseries(5)
            sage: L.abelian_variety()
	    Simple abelian variety J0(23) of dimension 2

            """
        return self._J

    def prime(self):
        r"""
        Returns the prime `p` as in 'p-adic L-function'.        
        
        EXAMPLES::
        
	    sage: L = J0(23)[0].padic_lseries(5)
	    sage: L.prime()
	    5

        """
        return self._p
		
    def _repr_(self):
        r"""
        Return print representation.

        EXAMPLES::
        
        """
        s = "%s-adic L-series of %s"%(self._p, self._J)
        if not self._normalize == 'L_ratio':
            s += ' (not normalized)'
        return s
   	
    def ap(self):
        """
	Return the Hecke eigenvalue $a_p$.

        EXAMPLES::

            sage: J = J0(23)[0]
	    sage: for p in prime_range(5,30):
	    ....:     L = J.padic_lseries(p)
	    ....:     p, L.ap()
	    ....:     
	    (5, 2*alpha)
	    (7, 2*alpha + 2)
	    (11, -2*alpha - 4)
	    (13, 3)
	    (17, -2*alpha + 2)
	    (19, -2)
	    (23, 1)
	    (29, -3)

	"""
        try:
            A = self._modular_symbols_subspace
        except AttributeError:
            A = self._modular_symbols_subspace = self.modular_symbols_subspace()
        a_p = self._ap = A.eigenvalue(self._p)
	return a_p 
 
    def is_ordinary(self):
    	"""
	Check if $p$ is an ordinary prime.
	"""
        try:
    	    K_f = self._hecke_eigenvalue_field
        except AttributeError:
            K_f = self._hecke_eigenvalue_field = self.hecke_eigenvalue_field()
        try:
            a_p = self._ap
        except AttributeError:
            a_p = self._ap = self.ap()
        frak_p = [x[0] for x in K_f.factor(self._p)]
        not_in_p = [x for x in frak_p if a_p not in frak_p]
        if len(not_in_p) == 0:
            return False
        else: 
            return True
        
    def measure(self,a,n,prec,quadratic_twist=+1,alpha=[]):
        """
	"""
        if quadratic_twist > 0:
            s = +1
        else:
            s = -1
        try:
            p, alpha, z, w, f = self.__measure_data[(n,prec,s)]
        except (KeyError, AttributeError):
            if not hasattr(self, '__measure_data'):
                self.__measure_data = {}
            p = self._p
            z = 1/(alpha**n)
            w = p**(n-1)
            f = self.modular_symbol

            self.__measure_data[(n,prec,s)] = (p,alpha,z,w,f)

        if quadratic_twist == 1:
            return z * f(a/(p*w)) - (z/alpha) * f(a/w)
        else:
            D = quadratic_twist
            chip = kronecker_symbol(D,p)
            if self._E.conductor() % p == 0:
                mu = chip**n * z * sum([kronecker_symbol(D,u) * f(a/(p*w)+ZZ(u)/D) for u in range(1,abs(D))])
            else:
                mu = chip**n * sum([kronecker_symbol(D,u) *(z * f(a/(p*w)+ZZ(u)/D) - chip *(z/alpha)* f(a/w+ZZ(u)/D)) for u in range(1,abs(D))])
            return s*mu
            

#    def measure(self, a, n, prec, quadratic_twist=+1 ):
#        r"""
#        Return the measure on `\ZZ_p^{\times}` defined by
#        
#           `\mu_{J,\alpha}^+ ( a + p^n \ZZ_p  ) =
#           \frac{1}{\alpha^n} \left [\frac{a}{p^n}\right]^{+} -
#           \frac{1}{\alpha^{n+1}} \left[\frac{a}{p^{n-1}}\right]^{+}`
#           
#        where `[\cdot]^{+}` is the modular symbol. This is used to define
#        this `p`-adic L-function (at least when the reduction is good).
#
#        The optional argument ``quadratic_twist`` replaces `J` by the twist in the above formula,
#        but the twisted modular symbol is computed using a sum over modular symbols of `J`
#        rather then finding the modular symbols for the twist.
#        
#        Note that the normalisation is not correct at this 
#        stage: use  ``_quotient_of periods`` and ``_quotient_of periods_to_twist`` 
#        to correct.
#
#        Note also that this function does not check if the condition
#        on the ``quadratic_twist=D`` is satisfied. So the result will only
#        be correct if for each prime `\ell` dividing `D`, we have
#        `ord_{\ell}(N)<= ord_{\ell}(D)`, where `N` is the level
#
#        INPUT:
#        
#        -  ``a`` - an integer
#        
#        -  ``n`` - a non-negative integer
#        
#        -  ``prec`` - an integer
#        
#        -  ``quadratic_twist`` (default = 1) - a fundamental discriminant of a quadratic field,
#           should be coprime to the level of `J`
#
#        EXAMPLES::
#        
#
#        """
#        
#        if quadratic_twist > 0:
#            s = +1
#        else:
#            s = -1
#        try:
#            p, alpha, z, w, f = self.__measure_data[(n,prec,s)]
#        except (KeyError, AttributeError):
#            if not hasattr(self, '__measure_data'):
#                self.__measure_data = {}
#            p = self._p
#
#            alpha = self.alpha(prec=prec)
#            print alpha
#        try:
#            if len(alpha)==1:
#                alpha = alpha[0]
#                z = 1/(alpha**n)
#                w = p**(n-1)
#                f = self.modular_symbol
  
#                self.__measure_data[(n,prec,s)] = (p,alpha,z,w,f)
         
#                if quadratic_twist == 1:
#                    return z * f(a/(p*w)) - (z/alpha) * f(a/w)
#                else:
#                    D = quadratic_twist
#                    chip = kronecker_symbol(D,p)
#                    if self._E.conductor() % p == 0:
#                        mu = chip**n * z * sum([kronecker_symbol(D,u) * f(a/(p*w)+ZZ(u)/D) for u in range(1,abs(D))])
#                    else:
#                        mu = chip**n * sum([kronecker_symbol(D,u) *(z * f(a/(p*w)+ZZ(u)/D) - chip *(z/alpha)* f(a/w+ZZ(u)/D)) for u in range(1,abs(D))])
#                    return s*mu
#            else:
#	        meas = []
#                for a in alpha:
#                    z = 1/(alpha**n)
#                    w = p**(n-1)
#                    f = self.modular_symbol
# 
#                    self.__measure_data[(n,prec,s)] = (p,alpha,z,w,f)
# 
#                    if quadratic_twist == 1:
#                        return z * f(a/(p*w)) - (z/alpha) * f(a/w)
#                    else:
#                        D = quadratic_twist
#                        chip = kronecker_symbol(D,p)
#                        if self._E.conductor() % p == 0:
#                            mu = chip**n * z * sum([kronecker_symbol(D,u) * f(a/(p*w)+ZZ(u)/D) for u in range(1,abs(D))])
#                        else:
#                            mu = chip**n * sum([kronecker_symbol(D,u) *(z * f(a/(p*w)+ZZ(u)/D) - chip *(z/alpha)* f(a/w+ZZ(u)/D)) for u in range(1,abs(D))])
#                        meas = meas + [s*mu]
#		    return meas

#        except TypeError:
#	    print "alpha = %s"%alpha
#            z = 1/(alpha**n)
#            w = p**(n-1)
#            f = self.modular_symbol

#            self.__measure_data[(n,prec,s)] = (p,alpha,z,w,f)

#            if quadratic_twist == 1:
#                return z * f(a/(p*w)) - (z/alpha) * f(a/w)
#            else:
#                D = quadratic_twist
#                chip = kronecker_symbol(D,p)
#                if self._E.conductor() % p == 0:
#                    mu = chip**n * z * sum([kronecker_symbol(D,u) * f(a/(p*w)+ZZ(u)/D) for u in range(1,abs(D))])
#                else:
#                    mu = chip**n * sum([kronecker_symbol(D,u) *(z * f(a/(p*w)+ZZ(u)/D) - chip *(z/alpha)* f(a/w+ZZ(u)/D)) for u in range(1,abs(D))])
#		return s*mu

        
    def alpha(self, prec=20):
        r"""
        Return a `p`-adic root `\alpha` of the polynomial `x^2 - a_p x
        + p` with `ord_p(\alpha) < 1`.  In the ordinary case this is
        just the unit root.

        INPUT:
        -  ``prec`` - positive integer, the `p`-adic precision of the root.

        EXAMPLES:

        """
        try:
            return self._alpha[prec]
        except AttributeError:
            self._alpha = {}
        except KeyError:
            pass
       
        J = self._J
        p = self._p
        Q = Qp(p)
        try:
            a_p = self._ap
        except AttributeError:
            a_p = self._ap = self.ap() 

        try:
            psis = self._psis
        except AttributeError:
            psis = self._psis = self.psi()

        K_f = self.hecke_eigenvalue_field()
        if len(psis) == 1:
   	    F = Q.extension(K_f.defining_polynomial(),names='a')
            a = F.gen()
            G = K_f.embeddings(K_f)
            if G[0](K_f.gen()) == K_f.gen():
                conj_map = G[1]
            else: 
                conj_map = G[0]
            v = self._dual_eigenvector
            v_conj = vector(conj_map(a) for a in v)
            a_p_conj = conj_map(a_p)
            R = F['x']
            x = R.gen()
	    psi = psis[0]
            a_p_padic = psi(a_p)
	    a_p_conj_padic = psi(a_p_conj)
            f = x**2 - (a_p_padic)*x + p
            fconj = x**2 - (a_p_conj_padic)*x + p
            norm_f = f*fconj
            norm_f_basefield = norm_f.change_ring(Q)
            FF = norm_f_basefield().factor()
            root0 = -f.gcd(FF[0][0])[0]
            root1 = -f.gcd(FF[1][0])[0]
            if root0.valuation() < 1:
                padic_lseries_alpha = [root0]
            else:
                padic_lseries_alpha = [root1]

        else:
            a_p_conj_padic = []
            a_p_padic = []
            for psi in psis:
                a_p_padic = a_p_padic + [psi(a_p)]

            R = Q['x']
            x = R.gen()
            padic_lseries_alpha = []
            for aps in a_p_padic:
                f = R(x**2 - aps*x + p)
                roots = f.roots()
                root0 = roots[0][0]
                root1 = roots[1][0]
                if root0.valuation() < 1:
                    padic_lseries_alpha = padic_lseries_alpha + [root0]
                else:
                    padic_lseries_alpha = padic_lseries_alpha + [root1]
                
        return padic_lseries_alpha  

    def order_of_vanishing(self):
        r"""
        Return the order of vanishing of this `p`-adic L-series.

        The output of this function is provably correct, due to a
        theorem of Kato [Ka].  This function will terminate if and only if
        the Mazur-Tate-Teitelbaum analogue [MTT] of the BSD conjecture about
        the rank of the curve is true and the subgroup of elements of
        `p`-power order in the Shafarevich-Tate group of this curve is
        finite.  I.e. if this function terminates (with no errors!),
        then you may conclude that the `p`-adic BSD rank conjecture is
        true and that the `p`-part of Sha is finite.

        NOTE: currently `p` must be a prime of good ordinary reduction.
        
        REFERENCES:
        
        - [MTT] B. Mazur, J. Tate, and J. Teitelbaum,
          On `p`-adic analogues of the conjectures of Birch and 
          Swinnerton-Dyer, Inventiones mathematicae 84, (1986), 1-48.

        - [Ka] Kayuza Kato, `p`-adic Hodge theory and values of zeta functions of modular
          forms, Cohomologies `p`-adiques et applications arithmetiques III,
          Asterisque vol 295, SMF, Paris, 2004.
        
        EXAMPLES::
        
        """
        try:
            return self.__ord
        except AttributeError:
            pass
        
        if not self.is_ordinary():
            raise NotImplementedError
        E = self.elliptic_curve()
        if not E.is_good(self.prime()):
            raise ValueError, "prime must be of good reduction"
        r = E.rank()
        n = 1
        while True:
            f = self.series(n)
            v = f.valuation()
            if v < r:
                raise RuntimeError, "while computing p-adic order of vanishing, got a contradiction: the curve is %s, the curve has rank %s, but the p-adic L-series vanishes to order <= %s"%(E, r, v)
            if v == r:
                self.__ord = v
                return v
            n += 1

    def _c_bounds(self, n):
        raise NotImplementedError

    def _prec_bounds(self, n,prec):
        raise NotImplementedError

    def teichmuller(self, prec):
        r"""
        Return Teichmuller lifts to the given precision.
        
        INPUT:
        
        - ``prec`` - a positive integer.
            
        OUTPUT:
        
        - a list of `p`-adic numbers, the cached Teichmuller lifts 

        EXAMPLES::
        
            sage: L = EllipticCurve('11a').padic_lseries(7)
            sage: L.teichmuller(1)
            [0, 1, 2, 3, 4, 5, 6]
            sage: L.teichmuller(2)
            [0, 1, 30, 31, 18, 19, 48]
        """
        p = self._p
        K = Qp(p, prec, print_mode='series')        
        return [Integer(0)] + \
               [a.residue(prec).lift() for a in K.teichmuller_system()]

    def _e_bounds(self, n, prec):
        p = self._p
        prec = max(2,prec)
        R = PowerSeriesRing(ZZ,'T',prec+1)
        T = R(R.gen(),prec +1)
        w = (1+T)**(p**n) - 1
        return [infinity] + [valuation(w[j],p) for j in range(1,min(w.degree()+1,prec))]
        
    def _get_series_from_cache(self, n, prec, D):
        try:
            return self.__series[(n,prec,D)]
        except AttributeError:
            self.__series = {}
        except KeyError:
            for _n, _prec, _D in self.__series.keys():
                if _n == n and _D == D and _prec >= prec:  
                    return self.__series[(_n,_prec,_D)].add_bigoh(prec)
        return None

    def _set_series_in_cache(self, n, prec, D, f):
        self.__series[(n,prec,D)] = f
  

    def _quotient_of_periods_to_twist(self,D):
        r"""
        For a fundamental discriminant `D` of a quadratic number field this computes the constant `\eta` such that
        `\sqrt{D}\cdot\Omega_{E_D}^{+} =\eta\cdot \Omega_E^{sign(D)}`. As in [MTT]_ page 40.
        This is either 1 or 2 unless the condition on the twist is not satisfied, e.g. if we are 'twisting back'
        to a semi-stable curve.
        
        REFERENCES:
        
        - [MTT] B. Mazur, J. Tate, and J. Teitelbaum,
          On `p`-adic analogues of the conjectures of Birch and 
          Swinnerton-Dyer, Invertiones mathematicae 84, (1986), 1-48.        
        
        .. note: No check on precision is made, so this may fail for huge `D`.
        
        EXAMPLES::
        
        """
        from sage.functions.all import sqrt
        # This funciton does not depend on p and could be moved out of this file but it is needed only here
        
        # Note that the number of real components does not change by twisting.
        if D == 1:
            return 1
        if D > 1:
            Et = self._E.quadratic_twist(D)
            qt = Et.period_lattice().basis()[0]/self._E.period_lattice().basis()[0]
            qt *= sqrt(qt.parent()(D))
        else:
            Et = self._E.quadratic_twist(D)
            qt = Et.period_lattice().basis()[0]/self._E.period_lattice().basis()[1].imag()
            qt *= sqrt(qt.parent()(-D))
        verbose('the real approximation is %s'%qt)
        # we know from MTT that the result has a denominator 1
        return QQ(int(round(8*qt)))/8
    

class pAdicLseriesOrdinary(pAdicLseries):
    """
    EXAMPLE:

        sage: from psage.modform.rational.padic_lseries import *
        sage: D = J0(188).decomposition()
        sage: A = D[-2]
        sage: L = pAdicLseriesOrdinary(A, 7)
        sage: L.series()
        O(7^20) + ((2*a + 4) + (6*a + 2)*7 + 6*7^2 + 4*a*7^3 + (4*a + 5)*7^4 + (4*a + 5)*7^5 + (6*a + 3)*7^6 + (2*a + 1)*7^7 + (5*a + 2)*7^8 + 5*7^9 + (2*a + 2)*7^10 + (5*a + 4)*7^11 + 3*a*7^12 + (5*a + 4)*7^13 + 5*a*7^14 + (3*a + 6)*7^15 + (5*a + 6)*7^16 + (6*a + 4)*7^17 + 5*a*7^18 + (3*a + 5)*7^19 + O(7^20))*T + ((5*a + 3) + (a + 6)*7 + (2*a + 1)*7^2 + (3*a + 2)*7^3 + (4*a + 2)*7^4 + 4*a*7^5 + (2*a + 6)*7^6 + 3*7^7 + (3*a + 5)*7^8 + (5*a + 2)*7^9 + (a + 3)*7^10 + 6*a*7^11 + (5*a + 5)*7^12 + (6*a + 6)*7^13 + (3*a + 4)*7^14 + (2*a + 4)*7^15 + (3*a + 6)*7^16 + (6*a + 1)*7^17 + 5*7^18 + (6*a + 5)*7^19 + O(7^20))*T^2 + ((3*a + 6) + (6*a + 6)*7 + 5*7^2 + 6*a*7^3 + (a + 4)*7^4 + (3*a + 2)*7^5 + 3*a*7^6 + (4*a + 6)*7^7 + (5*a + 2)*7^8 + 6*a*7^9 + (3*a + 1)*7^10 + (5*a + 3)*7^11 + (a + 1)*7^12 + (a + 4)*7^13 + (6*a + 2)*7^14 + (a + 2)*7^15 + (3*a + 1)*7^17 + (2*a + 5)*7^18 + (a + 3)*7^19 + O(7^20))*T^3 + ((3*a + 6) + (4*a + 3)*7 + (2*a + 3)*7^2 + (2*a + 2)*7^3 + 5*a*7^4 + 3*a*7^5 + (2*a + 5)*7^6 + (6*a + 1)*7^7 + (a + 2)*7^8 + (a + 1)*7^9 + (4*a + 5)*7^10 + 5*a*7^11 + (3*a + 4)*7^12 + (5*a + 6)*7^13 + (a + 5)*7^14 + 5*7^15 + 4*a*7^16 + (a + 3)*7^17 + (6*a + 2)*7^18 + (2*a + 3)*7^19 + O(7^20))*T^4 + O(T^5)
    """
    def series(self, n=2, quadratic_twist=+1, prec=5):
        r"""
        Returns the `n`-th approximation to the `p`-adic L-series as
        a power series in `T` (corresponding to `\gamma-1` with
        `\gamma=1+p` as a generator of `1+p\ZZ_p`).  Each
        coefficient is a `p`-adic number whose precision is provably
        correct.
        
        Here the normalization of the `p`-adic L-series is chosen
        such that `L_p(J,1) = (1-1/\alpha)^2 L(J,1)/\Omega_J`
        where `\alpha` is the unit root 

        INPUT:
        
        -  ``n`` - (default: 2) a positive integer
        -  ``quadratic_twist`` - (default: +1) a fundamental discriminant
           of a quadratic field, coprime to the 
           conductor of the curve
        -  ``prec`` - (default: 5) maximal number of terms of the series
           to compute; to compute as many as possible just
           give a very large number for ``prec``; the result will
           still be correct.

        ALIAS: power_series is identical to series.

        EXAMPLES:

	    sage: J = J0(188)[0]
	    sage: p = 7
	    sage: L = J.padic_lseries(p)
	    sage: L.is_ordinary()
	    True
	    sage: f = L.series(2)
	    sage: f[0]
	    O(7^20)
	    sage: f[1].norm()
	    3 + 4*7 + 3*7^2 + 6*7^3 + 5*7^4 + 5*7^5 + 6*7^6 + 4*7^7 + 5*7^8 + 7^10 + 5*7^11 + 4*7^13 + 4*7^14 + 5*7^15 + 2*7^16 + 5*7^17 + 7^18 + 7^19 + O(7^20)

        """
        n = ZZ(n)
        if n < 1:
            raise ValueError, "n (=%s) must be a positive integer"%n
        if not self.is_ordinary():
            raise ValueError, "p (=%s) must be an ordinary prime"%p
        # check if the conditions on quadratic_twist are satisfied
        D = ZZ(quadratic_twist)
        if D != 1:
            if D % 4 == 0:
                d = D//4
                if not d.is_squarefree() or d % 4 == 1:
                    raise ValueError, "quadratic_twist (=%s) must be a fundamental discriminant of a quadratic field"%D
            else:
                if not D.is_squarefree() or D % 4 != 1:
                    raise ValueError, "quadratic_twist (=%s) must be a fundamental discriminant of a quadratic field"%D
            if gcd(D,self._p) != 1:
                raise ValueError, "quadratic twist (=%s) must be coprime to p (=%s) "%(D,self._p)
            if gcd(D,self._E.conductor())!= 1:
                for ell in prime_divisors(D):
                    if valuation(self._E.conductor(),ell) > valuation(D,ell) :
                        raise ValueError, "can not twist a curve of conductor (=%s) by the quadratic twist (=%s)."%(self._E.conductor(),D)
                    
            
        p = self._p
        if p == 2 and self._normalize :
            print 'Warning : For p=2 the normalization might not be correct !'
        #verbose("computing L-series for p=%s, n=%s, and prec=%s"%(p,n,prec))
        
#        bounds = self._prec_bounds(n,prec)
#        padic_prec = max(bounds[1:]) + 5
        padic_prec = 10
#        verbose("using p-adic precision of %s"%padic_prec)
        
        res_series_prec = min(p**(n-1), prec)
        verbose("using series precision of %s"%res_series_prec)
        
        ans = self._get_series_from_cache(n, res_series_prec,D)
        if not ans is None:
            verbose("found series in cache")
            return ans
 
        K = QQ
        gamma = K(1 + p)
        R = PowerSeriesRing(K,'T',res_series_prec)
        T = R(R.gen(),res_series_prec )
        L = R(0) 
        one_plus_T_factor = R(1) 
        gamma_power = K(1)
        teich = self.teichmuller(padic_prec)
        p_power = p**(n-1)

        verbose("Now iterating over %s summands"%((p-1)*p_power))
        verbose_level = get_verbose()
        count_verb = 0
        alphas = self.alpha()
        Lprod = []
        for alpha in alphas:
            for j in range(p_power):
                s = K(0)
                if verbose_level >= 2 and j/p_power*100 > count_verb + 3:
                    verbose("%.2f percent done"%(float(j)/p_power*100))
                    count_verb += 3
                for a in range(1,p):
                    b = teich[a] * gamma_power
                    s = s + self.measure(b, n, padic_prec,D,alpha)
                L += s * one_plus_T_factor
                one_plus_T_factor *= 1+T
                gamma_power *= gamma
            Lprod = Lprod + [L]
            if len(Lprod)==1:
                return Lprod[0]
            else:
                return Lprod[0]*Lprod[1]

    power_series = series

