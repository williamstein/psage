r"""
Generator functions for the fourier expansion of Siegel modular forms. 

AUTHORS:

- Nils-Peter Skorupa 
  Factory function SiegelModularForm().
  Code parts concerning the Maass lift are partly based on code
  written by David Gruenewald in August 2006 which in turn is
  based on the PARI/GP code by Nils-Peter Skoruppa from 2003.
- Martin Raum (2009 - 07 - 28) Port to new framework.
- Martin Raum (2010 - 03 - 16) Implement the vector valued case.

REFERENCES:

- [Sko] Nils-Peter Skoruppa, ...
- [I-H] Tomoyoshi Ibukiyama and Shuichi Hayashida, ... 

"""

#===============================================================================
# 
# Copyright (C) 2009 Martin Raum
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================

from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_lazyelement import EquivariantMonoidPowerSeries_lazy
from psage.modform.paramodularforms.siegelmodularformg2_fourierexpansion import SiegelModularFormG2FourierExpansionRing
from sage.combinat.partition import number_of_partitions
from sage.misc.functional import isqrt
from sage.misc.latex import latex 
from sage.modular.modform.constructor import ModularForms
from sage.modular.modform.element import ModularFormElement
from sage.rings.arith import bernoulli, sigma
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.rational_field import QQ
from sage.structure.sage_object import SageObject
from psage.modform.paramodularforms.siegelmodularformg2_fourierexpansion import SiegelModularFormG2Filter_discriminant
from psage.modform.paramodularforms import siegelmodularformg2_misc_cython

#===============================================================================
# SiegelModularFormG2MaassLift
#===============================================================================

def SiegelModularFormG2MaassLift(f, g, precision, is_integral = True, weight = None) :    
    factory = SiegelModularFormG2Factory(precision)
    if is_integral :
        expansion_ring = SiegelModularFormG2FourierExpansionRing(ZZ)
    else :
        expansion_ring = SiegelModularFormG2FourierExpansionRing(QQ)

    if f == 0 :
        fexp = lambda p : 0
    if hasattr(f, "qexp") :
        fexp = lambda p : f.qexp(p)
    else :
        fexp = f
        
    if g == 0 :
        gexp = lambda p : 0
    if hasattr(g, "qexp") :
        gexp = lambda p : g.qexp(p)
    else :
        gexp = g
    
    if weight is None :
        try :
            weight = f.weight()
        except AttributeError :
            weight = g.weight() -2
    
    coefficients_factory = DelayedFactory_SMFG2_masslift( factory, fexp, gexp, weight )

    return EquivariantMonoidPowerSeries_lazy(expansion_ring, precision, coefficients_factory.getcoeff)

#===============================================================================
# DelayedFactory_SMFG2_masslift
#===============================================================================

class DelayedFactory_SMFG2_masslift :
    def __init__(self, factory, f, g, weight) :
        self.__factory = factory
        self.__f = f
        self.__g = g
        self.__weight = weight
        
        self.__series = None
    
    def getcoeff(self, key, **kwds) :
        (_, k) = key
        # for speed we ignore the character 
        if self.__series is None :
            self.__series = \
             self.__factory.maass_form( self.__f, self.__g, self.__weight,
                                        is_integral = True)
             
        return self.__series[k]

#===============================================================================
# SiegelModularFormG2Factory
#===============================================================================

_siegel_modular_form_g2_factory_cache = dict()

def SiegelModularFormG2Factory(precision) :
    if not isinstance(precision, SiegelModularFormG2Filter_discriminant) :
        precision = SiegelModularFormG2Filter_discriminant(precision)
    
    global _siegel_modular_form_g2_factory_cache
    try :
        return _siegel_modular_form_g2_factory_cache[precision]
    except KeyError :
        factory = SiegelModularFormG2Factory_class(precision)
        _siegel_modular_form_g2_factory_cache[precision] = factory
        
        return factory
    
#===============================================================================
# SiegelModularFormG2Factory_class
#===============================================================================

class SiegelModularFormG2Factory_class( SageObject ) :
    r"""
    A class producing dictionarys saving fourier expansions of Siegel
    modular forms of genus `2`.
    """
    
    def __init__(self, precision) :
        r"""
        INPUT:
        
        - ``prec``  -- a SiegelModularFormPrec, which will be the precision of
                       every object returned by an instance of this class.           
        """
        self.__precision = precision
        
        ## Conversion of power series is not expensive but powers of interger
        ## series are much cheaper then powers of rational series
        self._power_series_ring_ZZ = PowerSeriesRing(ZZ, 'q')
        self._power_series_ring = PowerSeriesRing(QQ, 'q')
    
    def power_series_ring(self) :
        return self._power_series_ring
    
    def integral_power_series_ring(self) :
        return self._power_series_ring_ZZ
    
    def _get_maass_form_qexp_prec(self) :
        r"""
        Return the precision of eta, needed to calculate the Maass lift in
        self.as_maass_spezial_form 
        """
        return ( self.__precision.discriminant() + 1)//4 + 1
    
    def _divisor_dict(self) :
        r"""
        Return a dictionary of assigning to each `k < n` a list of its divisors.
        
        INPUT:
        
        - `n`    -- a positive integer
        """
        try :
            return self.__divisor_dict
        except AttributeError :
            self.__divisor_dict = siegelmodularformg2_misc_cython.divisor_dict(self.__precision.discriminant())

            return self.__divisor_dict
        
    def _eta_power(self) :
        try :
            return self.__eta_power
        except AttributeError :
            qexp_prec = self._get_maass_form_qexp_prec()
        
            self.__eta_power = self.integral_power_series_ring() \
                 ([number_of_partitions(n) for n in xrange(qexp_prec)], qexp_prec)**6
                 
            return self.__eta_power

    def precision(self) :
        return self.__precision
    
    def _negative_fundamental_discriminants(self) :
        r"""
        Return a list of all negative fundamental discriminants.
        """
        try :
            return self.__negative_fundamental_discriminants
        except AttributeError :
            self.__negative_fundamental_discriminants = \
             siegelmodularformg2_misc_cython.negative_fundamental_discriminants(
                                              self.__precision.discriminant() )

            return self.__negative_fundamental_discriminants
    
    def maass_form( self, f, g, k = None, is_integral = False) :
        r"""
        Return the Siegel modular form `I(f,g)` (Notation as in [Sko]).
    
        INPUT:
        - `f`              -- modular form of level `1`
        - `g`              -- cusp form of level `1` and weight = ``weight of f + 2``
        - ``is_integral``  -- ``True`` if the result is garanteed to have integer
                              coefficients
        """
        
        ## we introduce an abbreviations
        if is_integral :
            PS = self.integral_power_series_ring()
        else :
            PS = self.power_series_ring()
        
        fismodular = isinstance(f, ModularFormElement)
        gismodular = isinstance(g, ModularFormElement)
    
        ## We only check the arguments if f and g are ModularFormElements.
        ## Otherwise we trust in the user 
        if fismodular and gismodular :
            assert( f.weight() + 2 == g.weight() | (f==0) | (g==0)), \
                    "incorrect weights!"
            assert( g.q_expansion(1) == 0), "second argument is not a cusp form"

        qexp_prec = self._get_maass_form_qexp_prec()
        if qexp_prec is None : # there are no forms below prec
            return dict()

        if fismodular :
            k = f.weight()
            if f == f.parent()(0) :
                f = PS(0, qexp_prec)
            else :
                f = PS(f.qexp(qexp_prec), qexp_prec)
        elif f == 0 :
            f = PS(0, qexp_prec)
        else :
            f = PS(f(qexp_prec), qexp_prec)
        
        if gismodular :
            k = g.weight() - 2
            if g == g.parent()(0) :
                g = PS(0, qexp_prec)
            else :
                g = PS(g.qexp(qexp_prec), qexp_prec)
        elif g == 0 :
            g = PS(0, qexp_prec)
        else :
            g = PS(g(qexp_prec), qexp_prec)
                
        if k is None :
            raise ValueError, "if neither f nor g are not ModularFormElements " + \
                              "you must pass k"
            
        fderiv = f.derivative().shift(1)
        f *= Integer(k/2)
        gfderiv = g - fderiv

        ## Form A and B - the Jacobi forms used in [Sko]'s I map.
        ## This is not necessary if we multiply Ifg0 and Ifg1 by etapow
        # (A0,A1,B0,B1) = (a0*etapow, a1*etapow, b0*etapow, b1*etapow)
    
        ## Calculate the image of the pair of modular forms (f,g) under
        ## [Sko]'s isomorphism I : M_{k} \oplus S_{k+2} -> J_{k,1}.
        
        # Multiplication of big polynomials may take > 60 GB, so wie have
        # to do it in small parts; This is only implemented for integral
        # coefficients.

        """
        Create the Jacobi form I(f,g) as in [Sko].
    
        It suffices to construct for all Jacobi forms phi only the part
        sum_{r=0,1;n} c_phi(r^2-4n) q^n zeta^r.
        When, in this code part, we speak of Jacobi form we only mean this part.
        
        We need to compute Ifg = \sum_{r=0,1; n} c(r^2-4n) q^n zeta^r up to
        4n-r^2 <= Dtop, i.e. n < precision
        """

        ## Create the Jacobi forms A=a*etapow and B=b*etapow in stages.
        ## Recall a = sum_{s != r mod 2} s^2*(-1)^r*q^((s^2+r^2-1)/4)*zeta^r
        ##        b = sum_{s != r mod 2}     (-1)^r*q^((s^2+r^2-1)/4)*zeta^r
        ## r, s run over ZZ but with opposite parities.
        ## For r=0, we need s odd, (s^2-1)/4 < precision, with s=2t+1 hence t^2+t < precision.
        ## For r=1, we need s even, s^2/4 < precision, with s=2t hence t^2 < precision.
    
        ## we use a slightly overestimated ab_prec 
        
        ab_prec = isqrt(qexp_prec + 1)
        a1dict = dict(); a0dict = dict()
        b1dict = dict(); b0dict = dict()
    
        for t in xrange(1, ab_prec + 1) :
            tmp = t**2
            a1dict[tmp] = -8*tmp
            b1dict[tmp] = -2
        
            tmp += t
            a0dict[tmp] = 8*tmp + 2
            b0dict[tmp] = 2
        b1dict[0] = -1
        a0dict[0] = 2; b0dict[0] = 2 
        
        a1 = PS(a1dict); b1 = PS(b1dict)
        a0 = PS(a0dict); b0 = PS(b0dict)

        ## Finally: I(f,g) is given by the formula below:
        ## We multiply by etapow explecitely and save two multiplications
        # Ifg0 = k/2*f*A0 - fderiv*B0 + g*B0 + O(q^precision)
        # Ifg1 = k/2*f*A1 - fderiv*B1 + g*B1 + O(q^precision)
        Ifg0 = (self._eta_power() * (f*a0 + gfderiv*b0)).list()
        Ifg1 = (self._eta_power() * (f*a1 + gfderiv*b1)).list()

        if len(Ifg0) < qexp_prec :
            Ifg0 += [0]*(qexp_prec - len(Ifg0))
        if len(Ifg1) < qexp_prec :
            Ifg1 += [0]*(qexp_prec - len(Ifg1))
        
        ## For applying the Maass' lifting to genus 2 modular forms.
        ## we put the coefficients of Ifg into a dictionary Chi
        ## so that we can access the coefficient corresponding to 
        ## discriminant D by going Chi[D].
        
        Cphi = dict([(0,0)])
        for i in xrange(qexp_prec) :
            Cphi[-4*i] = Ifg0[i]
            Cphi[1-4*i] = Ifg1[i]

        del Ifg0[:], Ifg1[:]

        """
        Create the Maas lift F := VI(f,g) as in [Sko].
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

        maass_coeffs = dict()
        divisor_dict = self._divisor_dict()

        ## First calculate maass coefficients corresponding to strictly positive definite matrices:        
        for disc in self._negative_fundamental_discriminants() :
            for s in xrange(1, isqrt((-self.__precision.discriminant()) // disc) + 1) :
                ## add (disc*s^2,t) as a hash key, for each t that divides s
                for t in divisor_dict[s] :
                    maass_coeffs[(disc * s**2,t)] = \
                       sum( a**(k-1) * Cphi[disc * s**2 / a**2] 
                            for a in divisor_dict[t] )

        ## Compute the coefficients of the Siegel form $F$:
        siegel_coeffs = dict()
        for (n,r,m), g in self.__precision.iter_positive_forms_with_content() :
            siegel_coeffs[(n,r,m)] = maass_coeffs[(r**2 - 4*m*n, g)]

        ## Secondly, deal with the singular part.
        ## Include the coeff corresponding to (0,0,0):
        ## maass_coeffs = {(0,0): -bernoulli(k)/(2*k)*Cphi[0]}
        siegel_coeffs[(0,0,0)] = -bernoulli(k)/(2*k)*Cphi[0]
        if is_integral :
            siegel_coeffs[(0,0,0)] = Integer(siegel_coeffs[(0,0,0)])
        
        ## Calculate the other discriminant-zero maass coefficients.
        ## Since sigma is quite cheap it is faster to estimate the bound and
        ## save the time for repeated calculation
        for i in xrange(1, self.__precision._indefinite_content_bound()) :
            ## maass_coeffs[(0,i)] = sigma(i, k-1) * Cphi[0]
            siegel_coeffs[(0,0,i)] = sigma(i, k-1) * Cphi[0]

        return siegel_coeffs

    def maass_eisensteinseries(self, k) :
        if not isinstance(k, (int, Integer)) :
            raise TypeError, "k must be an integer"
        if k % 2 != 0 or k < 4 :
            raise ValueError, "k must be even and greater than 4"
        
        return self.maass_form(ModularForms(1,k).eisenstein_subspace().gen(0), 0)
  
    ## TODO: Rework these functions
    def from_gram_matrix( self, gram, name = None) :
        raise NotImplementedError, "matrix argument not yet implemented"

    def _SiegelModularForm_borcherds_lift( self, f, prec = 100, name = None):
    
        raise NotImplementedError, "matrix argument not yet implemented"

    def _SiegelModularForm_yoshida_lift( self, f, g, prec = 100, name = None):
    
        raise NotImplementedError, "matrix argument not yet implemented"

    def _SiegelModularForm_from_weil_representation( self, gram, prec = 100, name = None):
    
        raise NotImplementedError, "matrix argument not yet implemented"

    def _SiegelModularForm_singular_weight( self, gram, prec = 100, name = None):
    
        raise NotImplementedError, "matrix argument not yet implemented"

    def _repr_(self) :
        return "Factory class for Siegel modular forms of genus 2 with precision %s" % \
               repr(self.__precision.discriminant())
    
    def _latex_(self) :
        return "Factory class for Siegel modular forms of genus $2$ with precision %s" % \
               latex(self.__precision.discriminant())
