r"""
The Hecke action on Siegel modular forms of genus 2.

AUTHORS:

- Martin Raum (2009 - 08 - 28) Initial version based on code by Nathan Ryan.
"""
from __future__ import division

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

from sage.misc.cachefunc import cached_method
from sage.misc.latex import latex
from sage.modular.modsym.p1list import P1List
from sage.rings.integer import Integer
from sage.structure.sage_object import SageObject
from psage.modform.paramodularforms import siegelmodularformg2_misc_cython

_siegelmodularformg2_heckeoperator_cache = dict()

#===============================================================================
# SiegelModularFormG2FourierExpansionHeckeAction
#===============================================================================

def SiegelModularFormG2FourierExpansionHeckeAction( n ) :
    global _siegelmodularformg2_heckeoperator_cache
    
    try :
        return _siegelmodularformg2_heckeoperator_cache[n]
    except KeyError :
        A = SiegelModularFormG2FourierExpansionHeckeAction_class(n)
        _siegelmodularformg2_heckeoperator_cache[n] = A
        
        return A

#===============================================================================
# SiegelModularFormG2FourierExpansionHeckeAction_class
#===============================================================================

class SiegelModularFormG2FourierExpansionHeckeAction_class (SageObject):
    r"""
    Implements a Hecke operator acting on Siegel modular.
    """

    def __init__(self, n):
        self.__l = n

        self.__l_divisors = siegelmodularformg2_misc_cython.divisor_dict(n + 1)

    def eval(self, expansion, weight = None) :
        if weight is None :
            try :
                weight = expansion.weight()
            except AttributeError :
                raise ValueError("weight must be defined for the Hecke action")
            
        precision = expansion.precision()
        if precision.is_infinite() :
            precision = expansion._bounding_precision()
        else :
            precision = precision._hecke_operator(self.__l)
        characters = expansion.non_zero_components()
        
        hecke_expansion = dict()
        for ch in characters :
            hecke_expansion[ch] = dict( (k, self.hecke_coeff(expansion, ch, k, weight)) for k in precision )
        
        result = expansion.parent()._element_constructor_(hecke_expansion)
        result._set_precision(expansion.precision()._hecke_operator(self.__l))
        
        return result
        
    # TODO : reimplement in cython
    def hecke_coeff(self, expansion, ch, a_b_c, k) :
        r"""
        Computes the coefficient indexed by `(a, b, c)` of `T(\ell) (F)`.
        """
        (a,b,c) = a_b_c
        character_eval = expansion.parent()._character_eval_function()
        
        ell = self.__l
        coeff = 0
        for t1 in self.__l_divisors[ell]:
            for t2 in self.__l_divisors[t1]:
                for V in self.get_representatives(t1/t2):
                    aprime, bprime, cprime = self.change_variables(V,(a,b,c))
                    if aprime % t1 == 0 and bprime % t2 == 0 and cprime % t2 == 0:
                        try:
                            coeff = coeff + character_eval(V, ch) * t1**(k-2)*t2**(k-1)*expansion[( ch, ((ell*aprime) //t1**2,
                                                                                                         (ell*bprime) //t1//t2,
                                                                                                         (ell*cprime) //t2**2) )]
                        except KeyError as msg:
                            raise ValueError('{0},{1}'.format(expansion,msg))
        return coeff

    @cached_method
    def get_representatives( self, t) :
        r"""
        A helper function used in hecke_coeff that computes the right
        coset representatives of $\Gamma^0(t)\SL(2,Z)$ where
        $\Gamma^0(t)$ is the subgroup of $SL(2,Z)$ where the upper right hand
        corner is divisible by $t$.

        NOTE:
        
            We use the bijection $\Gamma^0(t)\SL(2,Z) \rightarrow P^1(\Z/t\Z)$
            given by $A \mapsto [1:0]A$.
        """
        if t == 1 : return [(1,0,0,1)]
        
        rep_list = []
        
        for x,y in P1List(t):
            ## we calculate a pair c,d satisfying a minimality condition
            ## to make later multiplications cheaper
            (_, d, c) = Integer(x)._xgcd(Integer(y), minimal=True)
            rep_list.append((x,y,-c,d))
                
        return rep_list

    @staticmethod
    def change_variables(a_b_c_d, n_r_m):
        r"""
        A helper function used in hecke_coeff that computes the
        quadratic form [nprime,rprime,mprime] given by $Q((x,y) V)$
        where $Q=[n,r,m]$ and $V$ is a 2 by 2 matrix given by `(a, b, c, d)`
        """        
        (a,b,c,d) = a_b_c_d
        (n,r,m) = n_r_m
        return ( n*a**2 + r*a*b + m*b**2, 2*(n*a*c + m*b*d) + r*(a*d + c*b), \
                 n*c**2 + r*c*d + m*d**2 )

    def _repr_(self) :
        return '%s-th Hecke operator for Siegel modular forms' % self.__l

    def _latex_(self) :
        return latex(self.__n) + '-th Hecke operator for Siegel modular forms'
