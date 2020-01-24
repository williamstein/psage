r"""
Classes describing the Fourier expansion of Siegel modular forms of genus `4`.

AUTHOR :
    -- Martin Raum (2009 - 05 - 19) Initial version
"""

#===============================================================================
# 
# Copyright (C) 2010 Martin Raum
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

from builtins import map
from builtins import range
from copy import copy
from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_basicmonoids \
              import TrivialCharacterMonoid, TrivialRepresentation
from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ring \
              import EquivariantMonoidPowerSeriesRing
from operator import xor
from psage.modform.paramodularforms.siegelmodularformgn_fourierexpansion import SiegelModularFormGnFilter_diagonal_lll,\
    SiegelModularFormGnIndices_diagonal_lll
from sage.matrix.constructor import diagonal_matrix, matrix, zero_matrix, identity_matrix
from sage.misc.functional import isqrt
from sage.misc.latex import latex
from sage.rings.infinity import infinity
from sage.rings.integer_ring import ZZ
from functools import reduce

#===============================================================================
# SiegelModularFormG4Indices_diagonal_lll
#===============================================================================

class SiegelModularFormG4Indices_diagonal_lll ( SiegelModularFormGnIndices_diagonal_lll ) :
    r"""
    All positive definite quadratic forms filtered by their discriminant.
    """
    def __init__(self, reduced = True) :
        SiegelModularFormGnIndices_diagonal_lll.__init__(self, 4, reduced)

    def monoid(self) :
        return SiegelModularFormG4Indices_diagonal_lll(False) 
    
    def filter(self, bound) :
        return SiegelModularFormG4Filter_diagonal_lll(bound, self.is_reduced())
        
    def filter_all(self) :
        return SiegelModularFormG4Filter_diagonal_lll(infinity, self.is_reduced())
        
    def _reduction_function(self) :
        return self.reduce

    def decompositions(self, t) :
        ## We find all decompositions t1 + t2 = t
        ## We first find possible upper left matrices
        
        sub2 = list()
        for a0 in range(0, t[0,0] + 1, 2) :
            for a1 in range(0, t[1,1] + 1, 2) :
                # obstruction for t1[0,1]
                B1 = isqrt(a0 * a1)
                # obstruction for t2[0,1]
                B2 = isqrt((t[0,0] - a0) * (t[1,1] - a1))
                
                for b01 in range(max(-B1, t[0,1] - B2), min(B1, t[0,1] + B2) + 1) :
                    sub2.append((a0,a1,b01))
        
        sub3 = list()
        for (a0, a1, b01) in sub2 :
            sub3s = list()
            for a2 in range(0, t[2,2] + 1, 2) :
                # obstruction for t1[0,2]
                B1 = isqrt(a0 * a2)
                # obstruction for t2[0,2]
                B2 = isqrt((t[0,0] - a0) * (t[2,2] - a2))
                
                for b02 in range(max(-B1, t[0,2] - B2), min(B1, t[0,2] + B2) + 1) :
                    # obstruction for t1[1,2]
                    B3 = isqrt(a1 * a2)
                    # obstruction for t2[1,2]
                    B4 = isqrt((t[1,1] - a1) * (t[2,2] - a2))
                    
                    for b12 in range(max(-B3, t[1,2] - B4), min(B3, t[1,2] + B4) + 1) :
                        # obstruction for the minor [0,1,2] of t1
                        if a0*a1*a2 - a0*b12**2 + 2*b01*b12*b02 - b01**2*a2 - a1*b02**2 < 0 :
                            continue
                        # obstruction for the minor [0,1,2] of t2
                        if  (t[0,0] - a0)*(t[1,1] - a1)*(t[2,2] - a2) - (t[0,0] - a0)*(t[1,2] - b12)**2 \
                           + 2*(t[0,1] - b01)*(t[1,2] - b12)*(t[0,2] - b02) - (t[0,1] - b01)**2*(t[2,2] - a2) \
                           - (t[1,1] - a1)*(t[0,2] - b02)**2 < 0 :
                            continue
                        sub3s.append((a2, b02, b12))
            sub3.append((a0, a1, b01, sub3s))

        for (a0,a1,b01, sub3s) in sub3 :
            for (a2, b02, b12) in sub3s :
                for a3 in range(0, t[3,3] + 1, 2) :
                    # obstruction for t1[0,3]
                    B1 = isqrt(a0 * a3)
                    # obstruction for t2[0,3]
                    B2 = isqrt((t[0,0] - a0) * (t[3,3] - a3))
                
                    for b03 in range(max(-B1, t[0,3] - B2), min(B1, t[0,3] + B2) + 1) :
                        # obstruction for t1[1,3]
                        B3 = isqrt(a1 * a3)
                        # obstruction for t2[1,3]
                        B4 = isqrt((t[1,1] - a1) * (t[3,3] - a3))

                        for b13 in range(max(-B3, t[1,3] - B4), min(B3, t[1,3] + B4) + 1) :
                            # obstruction for the minor [0,1,3] of t1
                            if a0*a1*a3 - a0*b13**2 + 2*b01*b13*b03 - b01**2*a3 - a1*b03**2 < 0 :
                                continue
                            # obstruction for the minor [0,1,3] of t2
                            if  (t[0,0] - a0)*(t[1,1] - a1)*(t[3,3] - a3) - (t[0,0] - a0)*(t[1,3] - b13)**2 \
                               + 2*(t[0,1] - b01)*(t[1,3] - b13)*(t[0,3] - b03) - (t[0,1] - b01)**2*(t[3,3] - a3) \
                               - (t[1,1] - a1)*(t[0,3] - b03)**2 < 0 :
                                continue
                            
                            # obstruction for t1[2,3]
                            B3 = isqrt(a2 * a3)
                            # obstruction for t2[2,3]
                            B4 = isqrt((t[2,2] - a2) * (t[3,3] - a3))
                            
                            for b23 in range(max(-B3, t[2,3] - B4), min(B3, t[2,3] + B4) + 1) :
                                # obstruction for the minor [0,2,3] of t1
                                if a0*a2*a3 - a0*b23**2 + 2*b02*b23*b03 - b02**2*a3 - a2*b03**2 < 0 :
                                    continue
                                # obstruction for the minor [0,2,3] of t2
                                if  (t[0,0] - a0)*(t[2,2] - a2)*(t[3,3] - a3) - (t[0,0] - a0)*(t[2,3] - b23)**2 \
                                   + 2*(t[0,2] - b02)*(t[2,3] - b23)*(t[0,3] - b03) - (t[0,2] - b02)**2*(t[3,3] - a3) \
                                   - (t[2,2] - a2)*(t[0,3] - b03)**2 < 0 :
                                    continue
                                
                                # obstruction for the minor [1,2,3] of t1
                                if a1*a2*a3 - a1*b23**2 + 2*b12*b23*b13 - b12**2*a3 - a2*b13**2 < 0 :
                                    continue
                                # obstruction for the minor [1,2,3] of t2
                                if  (t[1,1] - a1)*(t[2,2] - a2)*(t[3,3] - a3) - (t[1,1] - a1)*(t[2,3] - b23)**2 \
                                   + 2*(t[1,2] - b12)*(t[2,3] - b23)*(t[1,3] - b13) - (t[1,2] - b12)**2*(t[3,3] - a3) \
                                   - (t[2,2] - a2)*(t[1,3] - b13)**2 < 0 :
                                    continue

                                t1 = matrix(ZZ, 4, [a0, b01, b02, b03, b01, a1, b12, b13, b02, b12, a2, b23, b03, b13, b23, a3], check = False)
                                if t1.det() < 0 :
                                    continue
                                t2 = t - t1
                                if t2.det() < 0 :
                                    continue
                                
                                t1.set_immutable()
                                t2.set_immutable()
                                
                                yield (t1, t2)

        raise StopIteration

    def __hash__(self) :
        return hash(self.is_reduced())
    
    def _repr_(self) :
        if self.is_reduced() :
            return "Reduced quadratic forms of rank 4 over ZZ"
        else :
            return "Quadratic forms over of rank 4 ZZ"
            
    def _latex_(self) :
        if self.is_reduced() :
            return "Reduced quadratic forms of rank $4$ over $\mathbb{Z}$"
        else :
            return "Quadratic forms over of rank $4$ over $\mathbb{Z}$"

#===============================================================================
# SiegelModularFormG4Filter_diagonal_lll
#===============================================================================

class SiegelModularFormG4Filter_diagonal_lll ( SiegelModularFormGnFilter_diagonal_lll ) :
    def __init__(self, bound, reduced = True) :
        SiegelModularFormGnFilter_diagonal_lll.__init__(self, 4, bound, reduced)

    def filter_all(self) :
        return SiegelModularFormG4Filter_diagonal_lll(infinity, self.is_reduced())
    
    def zero_filter(self) :
        return SiegelModularFormG4Filter_diagonal_lll(0, self.is_reduced()) 

    def _calc_iter_reduced_sub2(self) :
        try :
            return self.__iter_reduced_sub2
        except AttributeError :
            pass
        
        sub2 = list()
        for a0 in range(2, 2 * self.index(), 2) :
            for a1 in range(a0, 2 * self.index(), 2) :
                # obstruction for t[0,1]
                B1 = isqrt(a0 * a1 - 1)
                
                for b01 in range(0, min(B1, a0 // 2) + 1) :
                    sub2.append((a0,a1,-b01))
        
        self.__iter_reduced_sub2 = sub2

        return sub2

    def _calc_iter_reduced_sub3(self) :
        try :
            return self.__iter_reduced_sub3
        except AttributeError :
            pass
        
        sub2 = self._calc_iter_reduced_sub2()
        
        sub3 = list()
        for (a0, a1, b01) in sub2 :
            sub3s = list()
            for a2 in range(2, 2 * self.index(), 2) :
                # obstruction for t[0,2]
                B1 = isqrt(a0 * a2 - 1)
                
                for b02 in range(-B1, B1 + 1) :
                    # obstruction for t[1,2]
                    B3 = isqrt(a1 * a2 - 1)
                    
                    for b12 in range(-B3, B3 + 1) :
                        # obstruction for the minor [0,1,2] of t
                        if a0*a1*a2 - a0*b12**2 + 2*b01*b12*b02 - a2*b01**2 - a1*b02**2 <= 0 :
                            continue
                        
                        t = matrix(ZZ, 3, [a0, b01, b02, b01, a1, b12, b02, b12, a2])
                        u = t.LLL_gram()
                        if u.transpose() * t * u != t :
                            continue
                        
                        sub3s.append((a2, b02, b12))
                        
            sub3.append((a0, a1, b01, sub3s))
        
        self.__iter_reduced_sub3 = sub3
        
        return sub3

    def _calc_iter_reduced_sub4(self) :
        try :
            return self.__iter_reduced_sub4
        except AttributeError :
            pass
        
        sub3 = self._calc_iter_reduced_sub3()
        
        sub4 = list()
        for (a0,a1,b01,sub3s) in sub3 :
            sub4s = list()
            for (a2, b02, b12) in sub3s :
                sub4ss = list()
                for a3 in range(2, 2 * self.index() + 1, 2) :
                    # obstruction for t[0,3]
                    B1 = isqrt(a0 * a3 - 1)
                
                    for b03 in range(-B1, B1 + 1) :
                        # obstruction for t[1,3]
                        B3 = isqrt(a1 * a3 - 1)

                        for b13 in range(-B3, B3 + 1) :
                            # obstruction for the minor [0,1,3] of t
                            if a0*a1*a3 - a0*b12**2 + 2*b01*b13*b03 - b01**2*a3 - a1*b03**2 <= 0 :
                                continue
                            
                            # obstruction for t[2,3]
                            B3 = isqrt(a2 * a3 - 1)
                            
                            for b23 in range(-B3, B3 + 1) :
                                # obstruction for the minor [0,2,3] of t
                                if a0*a2*a3 - a0*b23**2 + 2*b02*b23*b03 - b02**2*a3 - a2*b03**2 <= 0 :
                                    continue
                                
                                # obstruction for the minor [1,2,3] of t
                                if a1*a2*a3 - a1*b23**2 + 2*b12*b23*b13 - b12**2*a3 - a2*b13**2 <= 0 :
                                    continue

                                t = matrix(ZZ, 4, [a0, b01, b02, b03, b01, a1, b12, b13, b02, b12, a2, b23, b03, b13, b23, a3], check = False)
                                if t.det() <= 0 :
                                    continue
                                
                                u = t.LLL_gram()
                                if u.transpose() * t * u != t :
                                    continue
                                                                   
                                sub4ss.append((a3, b03, b13, b23))
                sub4s.append((a2, b02, b12, sub4ss))
            sub4.append((a0, a1, b01, sub4s))
        
        self.__iter_reduced_sub4 = sub4

        return sub4

    def __iter__(self) :
        if self.index() is infinity :
            raise ValueError("infinity is not a true filter index")

        if self.is_reduced() :
            ## We only iterate positive definite matrices
            ## and later build the semidefinite ones
            ## We first find possible upper left matrices
            
            sub2 = self._calc_iter_reduced_sub2()
            sub3 = self._calc_iter_reduced_sub3()
            sub4 = self._calc_iter_reduced_sub4()

                
            t = zero_matrix(ZZ, 4)
            t.set_immutable()
            yield t
            
            for a0 in range(2, 2 * self.index(), 2) :
                t = zero_matrix(ZZ, 4)
                t[3,3] = a0
                t.set_immutable()
                
                yield t
                
            for (a0, a1, b01) in sub2 :
                t = zero_matrix(ZZ, 4)
                t[2,2] = a0
                t[3,3] = a1
                t[2,3] = b01
                t[3,2] = b01
                t.set_immutable()
                
                yield t
                
            for (a0, a1, b01, sub3s) in sub3 :
                t = zero_matrix(ZZ, 4)
                t[1,1] = a0
                t[2,2] = a1
                t[1,2] = b01
                t[2,1] = b01
                
                for (a2, b02, b12) in sub3s :
                    ts = copy(t)
                    ts[3,3] = a2
                    ts[1,3] = b02
                    ts[3,1] = b02
                    ts[2,3] = b12
                    ts[3,2] = b12
                    ts.set_immutable()
                    
                    yield ts
                    
            for (a0, a1, b01, sub4s) in sub4 :
                t = zero_matrix(ZZ, 4)
                t[0,0] = a0
                t[1,1] = a1
                t[0,1] = b01
                t[1,0] = b01
                
                for (a2, b02, b12, sub4ss) in sub4s :
                    ts = copy(t)
                    ts[2,2] = a2
                    ts[0,2] = b02
                    ts[2,0] = b02
                    ts[1,2] = b12
                    ts[2,1] = b12
                    
                    for (a3, b03, b13, b23) in sub4ss :
                        tss = copy(ts)
                        tss[3,3] = a3
                        tss[0,1] = b01
                        tss[1,0] = b01
                        tss[0,2] = b02
                        tss[2,0] = b02
                        tss[0,3] = b03
                        tss[3,0] = b03
                        tss.set_immutable()
                        
                        yield tss

        #! if self.is_reduced()
        else :
            ## We first find possible upper left matrices
            
            sub2 = list()
            for a0 in range(2 * self.index(), 2) :
                for a1 in range(2 * self.index(), 2) :
                    # obstruction for t[0,1]
                    B1 = isqrt(a0 * a1)
                    
                    for b01 in range(-B1, B1 + 1) :
                        sub2.append((a0,a1,b01))
                        
            sub3 = list()
            for (a0, a1, b01) in sub2 :
                sub3s = list()
                for a2 in range(2 * self.index(), 2) :
                    # obstruction for t[0,2]
                    B1 = isqrt(a0 * a2)
                    
                    for b02 in range(-B1, B1 + 1) :
                        # obstruction for t[1,2]
                        B3 = isqrt(a1 * a2)
                        
                        for b12 in range(-B3, B3 + 1) :
                            # obstruction for the minor [0,1,2] of t
                            if a0*a1*a2 - a0*b12**2 + 2*b01*b12*b02 - b01**2*a2 - a1*b02**2 < 0 :
                                continue
                            sub3s.append((a2, b02, b12))
                sub3.append((a0, a1, b01, sub3s))
    
            for (a0,a1,b01, sub3s) in sub3 :
                for (a2, b02, b12) in sub3s :
                    for a3 in range(2 * self.index(), 2) :
                        # obstruction for t[0,3]
                        B1 = isqrt(a0 * a3)
                    
                        for b03 in range(-B1, B1 + 1) :
                            # obstruction for t[1,3]
                            B3 = isqrt(a1 * a3)
    
                            for b13 in range(-B3, B3 + 1) :
                                # obstruction for the minor [0,1,3] of t
                                if a0*a1*a3 - a0*b13**2 + 2*b01*b13*b03 - b01**2*a3 - a1*b03**2 < 0 :
                                    continue
                                
                                # obstruction for t[2,3]
                                B3 = isqrt(a2 * a3)
                                
                                for b23 in range(-B3, B3 + 1) :
                                    # obstruction for the minor [0,2,3] of t
                                    if a0*a2*a3 - a0*b23**2 + 2*b02*b23*b03 - b02**2*a3 - a2*b03**2 < 0 :
                                        continue
                                    
                                    # obstruction for the minor [1,2,3] of t
                                    if a1*a2*a3 - a1*b23**2 + 2*b12*b23*b13 - b12**2*a3 - a2*b13**2 < 0 :
                                        continue
    
                                    t = matrix(ZZ, 4, [a0, b01, b02, b03, b01, a1, b12, b13, b02, b12, a2, b23, b03, b13, b23, a3], check = False)
                                    if t.det() < 0 :
                                        continue
                                    
                                    t.set_immutable()
                                    
                                    yield t

        raise StopIteration
    
    def iter_positive_forms(self) :
        if self.is_reduced() :
            sub4 = self._calc_iter_reduced_sub4()
            
            for (a0, a1, b01, sub4s) in sub4 :
                t = zero_matrix(ZZ, 4)
                t[0,0] = a0
                t[1,1] = a1
                t[0,1] = b01
                t[1,0] = b01
                
                for (a2, b02, b12, sub4ss) in sub4s :
                    ts = copy(t)
                    ts[2,2] = a2
                    ts[0,2] = b02
                    ts[2,0] = b02
                    ts[1,2] = b12
                    ts[2,1] = b12
                    
                    for (a3, b03, b13, b23) in sub4ss :
                        tss = copy(ts)
                        tss[3,3] = a3
                        tss[0,1] = b01
                        tss[1,0] = b01
                        tss[0,2] = b02
                        tss[2,0] = b02
                        tss[0,3] = b03
                        tss[3,0] = b03
                        tss.set_immutable()
                        
                        yield tss
        else :
            raise NotImplementedError

    def __hash__(self) :
        return reduce(xor, list(map(hash, [self.is_reduced(), self.index()])))
                   
    def _repr_(self) :
        return "Diagonal filter (%s)" % self.index()
    
    def _latex_(self) :
        return "Diagonal filter (%s)" % latex(self.index())

def SiegelModularFormG4FourierExpansionRing(K, genus) :

    R = EquivariantMonoidPowerSeriesRing(
         SiegelModularFormG4Indices_diagonal_lll(genus),
         TrivialCharacterMonoid("GL(%s,ZZ)" % (genus,), ZZ),
         TrivialRepresentation("GL(%s,ZZ)" % (genus,), K) )

    return R
