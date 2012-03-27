r"""
Elements of rings of paramodular forms of degree 2.

AUTHORS:

- Martin Raum (2010 - 05 - 06) Initial version.
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

from psage.modform.fourier_expansion_framework.modularforms.modularform_element import ModularForm_generic
from psage.modform.paramodularforms.paramodularformd2_fourierexpansion import ParamodularFormD2Filter_trace
from psage.modform.paramodularforms.paramodularformd2_fourierexpansion_cython import apply_GL_to_form
from sage.functions.other import ceil
from sage.modular.modsym.p1list import P1List
from sage.rings.all import Integer
from sage.rings.complex_field import is_ComplexField
from sage.structure.sequence import Sequence

class ParamodularFormD2_generic :

    def is_cusp_form(self) :
        raise NotImplementedError
        
class ParamodularFormD2_classical ( ParamodularFormD2_generic, ModularForm_generic ) :

    def is_gritsenko_form(self) :
        raise NotImplementedError    

    def atkin_lehner_eigenvalue_numerical(self, (tau1, z, tau2)) :
        from sage.libs.mpmath import mp
        from sage.libs.mpmath.mp import exp, pi
        from sage.libs.mpmath.mp import j as i
        
        if not Integer(self.__level()).is_prime() :
            raise ValueError, "The Atkin Lehner involution is only unique if the level is a prime"
        
        precision = ParamodularFormD2Filter_trace(self.precision())
        
        s = Sequence([tau1, z, tau2])
        if not is_ComplexField(s) :
            mp_precision = 30
        else :
            mp_precision = ceil(3.33 * s.universe().precision() ) 
        mp.dps = mp_precision
                
        (tau1, z, tau2) = tuple(s)
        (tau1p, zp, tau2p) = (self.level()*tau1, self.level()*z, self.level()*tau2)
        
        (e_tau1, e_z, e_tau2) = (exp(2 * pi * i * tau1), exp(2 * pi * i * z), exp(2 * pi * i * tau2))
        (e_tau1p, e_zp, e_tau2p) = (exp(2 * pi * i * tau1p), exp(2 * pi * i * zp), exp(2 * pi * i * tau2p))
        
        self_value = s.universe().zero()
        trans_value = s.universe().zero()
        
        for k in precision :
            (a,b,c) = apply_GL_to_form(self._P1List()(k[1]), k[0])
            
            self_value = self_value + self[k] * (e_tau1**a * e_z**b * e_tau2**c)
            trans_value = trans_value + self[k] * (e_tau1p**a * e_zp**b * e_tau2p**c)
        
        return trans_value / self_value
    
    def schmidt_t5_eigenvalue_numerical(self, (tau1, z, tau2)) :
        from sage.libs.mpmath import mp
        from sage.libs.mpmath.mp import exp, pi
        from sage.libs.mpmath.mp import j as i
        
        if not Integer(self.__level()).is_prime() :
            raise ValueError, "T_5 is only unique if the level is a prime"
        
        precision = ParamodularFormD2Filter_trace(self.precision())
        
        s = Sequence([tau1, z, tau2])
        if not is_ComplexField(s) :
            mp_precision = 30
        else :
            mp_precision = ceil(3.33 * s.universe().precision())
        mp.dps = mp_precision
        
        p1list = P1List(self.level())
        
        ## Prepare the operation for d_1(N)
        ## We have to invert the lifts since we will later use apply_GL_to_form
        d1_matrices = [p1list.lift_to_sl2z(i) for i in xrange(len(p1list))]
        d1_matrices = map(lambda (a,b,c,d): (d,-b,-c,a), d1_matrices)
        
        ## Prepare the evaluation points corresponding to d_02(N)
        d2_points = list()
        for i in xrange(len(p1list())) :
            (a, b, c, d) = p1list.lift_to_sl2z(i)
            tau1p = (a * tau1 + b) / (c * tau1 + d)
            zp = z / (c * tau1 + d)
            tau2p = tau2 - c * z**2 / (c * tau1 + d)
            
            (e_tau1p, e_zp, e_tau2p) = (exp(2 * pi * i * tau1p), exp(2 * pi * i * zp), exp(2 * pi * i * tau2p))
            d2_points.append((e_tau1p, e_zp, e_tau2p))
            
        (e_tau1, e_z, e_tau2) = (exp(2 * pi * i * tau1), exp(2 * pi * i * z), exp(2 * pi * i * tau2))
            
        self_value = s.universe().zero()
        trans_value = s.universe().zero()
        
        for k in precision :
            (a,b,c) = apply_GL_to_form(self._P1List()(k[1]), k[0])
            
            self_value = self_value + self[k] * e_tau1**a * e_z**b * e_tau2**c
            for m in d1_matrices :
                (ap, bp, cp) = apply_GL_to_form(m, (a,b,c))
                
                for (e_tau1p, e_zp, e_tau2p) in d2_points :
                    trans_value = trans_value + self[((ap,bp,cp),0)] * e_tau1p**ap * e_zp**bp * e_tau2p**cp
        
        return trans_value / self_value
        