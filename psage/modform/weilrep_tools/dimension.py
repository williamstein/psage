r"""
Class representing a space of vector valued modular forms
transforming with the Weil representation of a discriminant module.
...

AUTHORS:

- Stephan Ehlen (2012-11-12): initial version
...

"""

#*****************************************************************************
#       Copyright (C) 2012 Stephan Ehlen <ehlen@mathematik.tu-darmstadt.de>
#
#  Distributed under the terms of the GNU General Public License (GPLv2)
#  as published by the Free Software Foundation
#                  http://www.gnu.org/licenses/
#*****************************************************************************

#from psage.modules import *
from sage.all import SageObject, Integer, RR, is_odd, next_prime, floor, RealField, ZZ, ceil, log, ComplexField, real, sqrt, exp
import sys

try:
    from psage.modules.finite_quadratic_module import FiniteQuadraticModule
    from psage.modules.invariants import cython_invariants_dim
except ImportError:
    raise

try:
    from sage_code.fqm.genus_symbol import GenusSymbol
except ImportError:
    GENUS_SYMBOLS = False

def BB(x):
    RF=RealField(100)
    x = RF(x)
    onehalf = RF(1)/2
    return x - onehalf*(floor(x)-floor(-x))

class VectorValuedModularForms(SageObject):
    r"""
    Class representing a space of vector valued modular forms
    transforming with the Weil representation of a discriminant module.


    EXAMPLES:

        As input we use any input that one can use for a finite quadratic module.

        (For instance a genus symbol.)

        ::

            sage: V=VectorValuedModularForms('2_7^+1.3^-2')
            sage: V
            Vector valued modular forms for the Weil representation corresponding to: 
            Finite quadratic module in 3 generators:
            gens: e0, e1, e2
            form: 3/4*x0^2 + 2/3*x1^2 + 1/3*x2^2

    SEE ALSO:

        :class:`psage.modules.finite_quadratic_modules.FiniteQuadraticModule`
    """

    def __init__(self, A, use_genus_symbols = False, aniso_formula = False, use_reduction = False):
        self._use_reduction = use_reduction
        if not GENUS_SYMBOLS:
            use_genus_symbols = False
        if use_genus_symbols:
            if isinstance(A, str):
                g = GenusSymbol(A)
            else:
                try:
                    g = GenusSymbol(A.jordan_decomposition().genus_symbol())
                except:
                    raise ValueError
            self._g = g
            n2 = self._n2 = g.torsion(2)
            self._v2 = g.two_torsion_values()
            self._M = None
            self._aniso_formula = aniso_formula
        else:
            self._M = FiniteQuadraticModule(A)
            self._g = None
            self._level = self._M.level()
            self._aniso_formula = False
            if is_odd(self._M.order()):
                self._n2 = n2 = 1
                self._v2 = {0: 1}
            else:
                self._M2 = M2 = self._M.kernel_subgroup(2).as_ambient()[0]
                self._n2 = n2 = self._M2.order()
                self._v2 = self._M2.values()

        if use_genus_symbols:
            self._signature = g.signature()
            m = g.order()
        else:
            self._signature = self._M.signature()
            m = self._M.order()
            
        self._m = m
        d = Integer(1)/Integer(2)*(m+n2) # |discriminant group/{+/-1}|
        self._d = d
        self._alpha3 = None
        self._alpha4 = None

    def __repr__(self):
        return "Vector valued modular forms for the Weil representation corresponding to: \n" + self._M.__repr__()

    def finite_quadratic_module(self):
        return self._M

    def dimension(self,k,ignore=False, debug = 0):
        if k < 2 and not ignore:
            raise NotImplementedError("k has to >= 2")
        s = self._signature
        if not (2*k in ZZ):
            raise ValueError("k has to be integral or half-integral")
        if (2*k+s)%4 != 0 and not ignore:
            raise NotImplementedError("2k has to be congruent to -signature mod 4")
        if self._v2.has_key(0):
            v2 = self._v2[0]
        else:
            v2 = 1

        if self._g != None:
            if not self._aniso_formula:
                vals = self._g.values()
            #else:
                #print "using aniso_formula"
            M = self._g
        else:
            vals = self._M.values()
            M = self._M
            
        prec = ceil(max(log(M.order(),2),52)+1)
        #print prec
        RR = RealField(prec)
        CC = ComplexField(prec)
        d = self._d
        m = self._m
        if debug > 0: print d,m
            
        if self._alpha3 == None:
            if self._aniso_formula:
                self._alpha4 = 1
                self._alpha3 = -sum([BB(a)*mm for a,mm in self._v2.iteritems() if a != 0])
                #print self._alpha3
                self._alpha3 += Integer(d) - Integer(1) - self._g.a5prime_formula()
                #print self._alpha3, self._g.a5prime_formula()
                self._alpha3 = self._alpha3/RR(2)
            else:
                self._alpha3 = sum([(1-a)*mm for a,mm in self._v2.iteritems() if a != 0])
                #print self._alpha3
                self._alpha3 += sum([(1-a)*mm for a,mm in vals.iteritems() if a != 0])
                #print self._alpha3
                self._alpha3 = self._alpha3 / Integer(2)
                self._alpha4 = 1/Integer(2)*(vals[0]+v2) # the codimension of SkL in MkL
        alpha3 = self._alpha3
        alpha4 = self._alpha4
        if debug > 0: print alpha3, alpha4
        g1=M.char_invariant(1)
        g1=CC(g1[0]*g1[1])
        #print g1
        g2=M.char_invariant(2)
        g2=RR(real(g2[0]*g2[1]))
        if debug > 0: print g2, g2.parent()
        g3=M.char_invariant(-3)
        g3=CC(g3[0]*g3[1])
        if debug > 0: print RR(d) / RR(4), sqrt(RR(m)) / RR(4), CC(exp(2 * CC.pi() * CC(0,1) * (2 * k + s) / Integer(8)))
        alpha1 = RR(d) / RR(4) - (sqrt(RR(m)) / RR(4)  * CC(exp(2 * CC.pi() * CC(0,1) * (2 * k + s) / Integer(8))) * g2)
        if debug > 0: print alpha1
        alpha2 = RR(d) / RR(3) + sqrt(RR(m)) / (3 * sqrt(RR(3))) * real(exp(CC(2 * CC.pi() * CC(0,1) * (4 * k + 3 * s - 10) / 24)) * (g1+g3))
        if debug > 0: print alpha1, alpha2, g1, g2, g3, d, k, s
        dim = real(d + (d * k / Integer(12)) - alpha1 - alpha2 - alpha3)
        if debug > 0:
            print dim
        if abs(dim-round(dim)) > 1e-6:
            raise RuntimeError("Error ({0}) too large in dimension formula".format(abs(dim-round(dim))))
        dim = round(dim)
        if k >=2 and dim < 0:
            raise RuntimeError("Negative dimension!")
        return dim

    def dimension_cusp_forms(self, k, ignore=False, no_inv = False, test_positive = False, proof = False):
        if k == Integer(3)/2:
            dim = self.dimension(k, True) - self._alpha4
            if not test_positive or dim <= 0:
                if self._M == None:
                    self._M = self._g.finite_quadratic_module()
                corr = weight_one_half_dim(self._M, self._use_reduction, proof = proof)
                #print "weight one half: {0}".format(corr)
                dim += corr
        else:
            dim = self.dimension(k,ignore)-self._alpha4
        if k == Integer(2) and not no_inv:
            if self._M == None:
                self._M = self._g.finite_quadratic_module()
            if self._M.level() == 1:
                return dim + 1
            dinv = cython_invariants_dim(self._M,self._use_reduction)
            dim = dim + dinv
        if dim < 0:
            raise RuntimeError("Negative dimension!")
        return dim
        
def test_real_quadratic(minp=1,maxp=100,minwt=2,maxwt=1000):
    for p in prime_range(minp,maxp):
        if p%4==1:
            print "p = ", p
            gram=Matrix(ZZ,2,2,[2,1,1,(1-p)/2])
            M=VectorValuedModularForms(gram)
            if is_odd(minwt):
                minwt=minwt+1
            for kk in range(minwt,round(maxwt/2-minwt)):
                k = minwt+2*kk
                if M.dimension_cusp_forms(k)-dimension_cusp_forms(kronecker_character(p),k)/2 != 0:
                    print "ERROR: ", k, M.dimension_cusp_forms(k), dimension_cusp_forms(kronecker_character(p),k)/2
                    return False
    return true

#sys.path.append('/home/stroemberg/Programming/Sage/sage-add-ons3/nils')
#from jacobiforms.all import *

def test_jacobi(index=1,minwt=4,maxwt=100,eps=-1):
    m=index
    gram=Matrix(ZZ,1,1,[-eps*2*m])
    M=VectorValuedModularForms(gram)
    if is_odd(minwt):
        minwt=minwt+1
    for kk in range(0,round(maxwt/2-minwt)+2):
        k = minwt+2*kk+(1+eps)/2
        print "Testing k = ", k
        if eps==-1:
            dimJ=dimension_jac_forms(k,m,-1)
            dimV=M.dimension(k-1/2)
        else:
            dimJ=dimension_jac_cusp_forms(k,m,1)
            dimV=M.dimension_cusp_forms(k-1/2)
        if dimJ-dimV != 0:
            print "ERROR: ", k, dimJ, dimV
            return False
    return true
