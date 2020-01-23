"""
Triple product L-series.

WARNING: The code in here gets the algorithm down, but is (1) pure
Python, and (2) runs into a limit since in order to actually use it in
even the first example, one needs to send a huge number of
coefficients to GP/PARI over a ptty, which results in a crash.
So currently this code does not actually work in a single case!
To fix it:

   * find a better way to send a large number of coefficients to pari
   
   * make a version of the code that instead uses lcalc

Note that this code still has some value, since it can be used to
compute the Dirichlet series coefficients of the triple product.
"""

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


from builtins import range
from builtins import object
from sage.all import ZZ, RDF, CDF
R_cdf = CDF['x']

pi = RDF.pi()

def quad_roots(a, p):
    """
    Return the two complex roots of X^2 - a*X + p.
    
    INPUT:
        - `a` -- a real number
        - `p` -- a prime number
    OUTPUT:
        - 2-tuple of complex numbers

    EXAMPLES::
    """
    t = R_cdf([p, -a, 1]).roots()
    return (t[0][0], t[1][0])

class TripleProductLseries(object):
    """
    A triple product `L`-series attached to three newforms of level `N`.
    """
    def __init__(self, N, f, g, h):
        """
        INPUT:
            - N -- squarefree integer
            - f -- an object such that for n>=0, we have
              f[n] = a_n(f) = n-th Fourier coefficient of
              a newform on Gamma_0(N).
            - g -- like f
            - h -- like f
        """
        self._N = ZZ(N)
        if not (self._N.is_squarefree() and self._N > 0):
            raise ValueError("N (={0}) must be a squarefree positive integer".format(self._N))
        self._newforms = (f,g,h)
        self._gen = RDF['X'].gen()
        self._genC = CDF['X'].gen()
        self._series = RDF[['X']]

    def __repr__(self):
        """
        Return text representation of this triple product `L`-function.
        """
        return "Triple product L-function L(s,f,g,h) of three newforms on Gamma_0(%s)"%self._N

    def level(self):
        """
        Return the common level `N` of the newforms in the triple product.

        OUTPUT:
            - Integer
        
        EXAMPLES::
        """
        return self._N

    def newforms(self):
        """
        Return 3-tuple (f,g,h) of the data that defines the newforms
        in the triple product.

        OUTPUT:
            - 3-tuple
        
        EXAMPLES::
        """
        return self._newforms

    def _local_series(self, p, prec):
        """
        Return power series in `X` (which you should think of as `p^{-s}`)
        that is the expansion to precision prec of the local factor at `p`
        of this `L`-series.

        INPUT:
            - p -- prime
            - prec -- positive integer

        OUTPUT:
            - power series that ends in a term ``O(X^prec)``.
        """
        f = self._series(self.charpoly(p), prec)
        return f**(-1)

    def dirichlet_series(self, prec, eps=1e-10):
        """
        Return the Dirichlet series representation of self, up to the given
        precision.

        INPUT:
           - prec -- positive integer
           - eps -- None or a positive real; any coefficient with absolute
             value less than eps is set to 0.
        """
        coeffs = self.dirichlet_series_coeffs(prec, eps)
        return DirichletSeries(coeffs, 's')
            
    def dirichlet_series_coeffs(self, prec, eps=1e-10):
        """
        Return the coefficients of the Dirichlet series representation
        of self, up to the given precision.

        INPUT:
           - prec -- positive integer
           - eps -- None or a positive real; any coefficient with absolute
             value less than eps is set to 0.
        """
        # Use multiplicativity to compute the Dirichlet series
        # coefficients, then make a DirichletSeries object.
        zero = RDF(0)
        coeffs = [RDF(0),RDF(1)] + [None]*(prec-2)

        from sage.all import log, floor   # TODO: slow
        
        # prime-power indexed coefficients
        for p in prime_range(2, prec):
            B = floor(log(prec, p)) + 1
            series = self._local_series(p, B)
            p_pow = p
            for i in range(1, B):
                coeffs[p_pow] = series[i] if (eps is None or abs(series[i])>eps) else zero
                p_pow *= p

        # non-prime-powers
        from sage.all import factor
        for n in range(2, prec):
            if coeffs[n] is None:
                a = prod(coeffs[p**e] for p, e in factor(n))
                coeffs[n] = a if (eps is None or abs(a) > eps) else zero

        return coeffs

    def charpoly(self, p):
        """
        Return the denominator as a polynomial in `X` (=`p^{-s}`) of the local
        factor at `p` of this `L`-series.

        The degree of the polynomial is ???? [[todo]].

        INPUT:
            - `p` -- prime

        OUTPUT:
            - polynomial in `X`

        EXAMPLES::

        """
        if self._N % p == 0:
            return self._charpoly_bad(p)
        else:
            return self._charpoly_good(p)

    def _charpoly_good(self, p):
        """
        Internal function that returns the local charpoly at a good prime.
        
        INPUT:
            - `p` -- prime

        OUTPUT:
            - polynomial in `X`

        EXAMPLES::

        """
        Y = self._genC
        a = [quad_roots(f[p], p) for f in self._newforms]
        L = 1
        for n in range(8):
            d = ZZ(n).digits(2)
            d = [0]*(3-len(d)) + d
            L *= 1 - prod(a[i][d[i]] for i in range(3))*Y
        return self._gen.parent()([x.real_part() for x in L])

    def _charpoly_bad(self, p):
        """
        Internal function that returns the local charpoly at a bad prime.
        
        INPUT:
            - `p` -- prime

        OUTPUT:
            - polynomial in `X`

        EXAMPLES::

        """
        X = self._gen
        a_p, b_p, c_p = [f[p] for f in self._newforms]
        return (1 - a_p*b_p*c_p * X) * (1 - a_p*b_p*c_p*p*X)**2

    def epsilon(self, p=None):
        """
        Return the local or global root number of this triple product
        L-function.

        INPUT:
            - p -- None (default) or a prime divisor of the level
        """
        if p is None:
            # Right below equation (1.11) in [Gross-Kudla]
            return -prod(self.epsilon(p) for p in self.level().prime_divisors())
        else:
            if not ZZ(p).is_prime():
                raise ValueError("p must be prime")
            if self.level() % p != 0:
                raise ValueError("p must divide the level")
            # Equation (1.3) in [Gross-Kudla]
            a_p, b_p, c_p = [f[p] for f in self._newforms]
            return -a_p*b_p*c_p

    def dokchitser(self, prec):
        # NOTE: In order to get the Dokchitser parameters of an L-function,
        # it is useful to know that
        #
        #         gamma(s) = 2^s*gamma(s/2)*gamma((s+1)/2) / (2*sqrt(pi))
        #
        conductor = self.level()**10
        gammaV = [-1,-1,-1,0,0,0,0,1]
        weight = 4
        eps = self.epsilon()
        poles = []
        residues = []
        
        from sage.lfunctions.dokchitser import Dokchitser
        L = Dokchitser(conductor = conductor,
                       gammaV = gammaV,
                       weight = weight,
                       eps = eps,
                       poles = poles, residues=[])
        #s = 'v=%s; a(k)=if(k>%s,0,v[k])'%(self.dirichlet_series_coeffs(prec), prec)
        s = 'v=%s; a(k)=v[k]'%(self.dirichlet_series_coeffs(prec))
        L.init_coeffs('a(k)', pari_precode=s)
        return L
        
        

class DirichletSeries(object):
    """
    A Dirichlet series.
    """
    def __init__(self, coeffs, variable='s'):
        """
        INPUT:
            - ``coeffs`` -- a list of the coefficients of the Dirichlet series
              such that the coefficient `a_n` in the sum `a_n/n^s` is ``coeffs[n]``.
            - ``variable`` - a string
        """
        self._coeffs = coeffs
        self._variable = variable

    def __repr__(self):
        if self._coeffs[0]._is_atomic():
            v = ['%s/%s^%s'%(self._coeffs[n], n, self._variable) for
                 n in range(1,len(self._coeffs)) if self._coeffs[n]]
            s = ' + '.join(v)
            s = s.replace(' + -',' - ')
        else:
            v = ['(%s)/%s^%s'%(self._coeffs[n], n, self._variable) for
                 n in range(1,len(self._coeffs)) if self._coeffs[n]]
            s = ' + '.join(v)
        return s

    def __getitem__(self, *args):
        return self._coeffs.__getitem__(*args)

                          
