# -*- coding: utf-8 -*-
#*****************************************************************************
#  Copyright (C) 2010 Fredrik Strömberg <stroemberg@mathematik.tu-darmstadt.de>,
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
#*****************************************************************************

r"""

This file implements:
- Vector-valued harmonic weak Maass forms for the above Weil representation.

AUTHORS:

- Fredrik Strömberg

EXAMPLES::

# Construct the vector-valued Harmonic weak maass form corresponding to the holomorphic modular form of weight 2 on Gamma0(11)

    sage: WR = WeilRepDiscriminantForm(11,dual=True);WR
    Dual of Weil representation of the discriminant form given by ZZ/22ZZ with quadratic form Q(x)=11*x**2 mod 1.
    sage: M=VVHarmonicWeakMaassForms(WR,0.5,20);M
    Space of Vector-Valued harmonic weak Maass forms on Modular Group SL(2,Z) of weight 1/2  and values in CC[ZZ/22ZZ].
    Representation is Dual of Weil representation of the discriminant form given by ZZ/22ZZ with quadratic form Q(x)=11*x**2 mod 1.
    
## Vector-valued weakly holomorphic modular form who is a generating function of the partition function.
    sage: H=VVHarmonicWeakMaassForms(-6,-1/2,holomorphic=True)
    sage: PP={(1/12,0):1,(5/12,0):-1}
    sage: F=H.get_element(PP,maxC=10)




"""
from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

from builtins import map
from builtins import str
from builtins import range
import mpmath as mpmath
# ugly fix to support pickling of mpmath matrices
import mpmath.matrices.matrices
mpmath.matrices.matrices.matrix = mpmath.matrix
import random
import tempfile,os
from sage.all import Parent,SageObject,Integer,Rational,SL2Z,QQ,ZZ,CC,RR,Newform,sign,Newforms,RealField,save,load,\
    latex,log_b,vector,dimension_new_cusp_forms,dimension_cusp_forms,is_square,kronecker,numerator,denominator,\
    cached_function,is_odd,MatrixSpace,is_fundamental_discriminant
from sage.rings.complex_mpc import MPComplexField,MPComplexNumber

from psage.modform.arithgroup.mysubgroup import *
from .automorphic_forms import *
from .weil_rep_simple import *
from .vv_harmonic_weak_maass_forms_alg import *
from psage.modules.vector_complex_dense import *
from psage.matrix.matrix_complex_dense import *




mp1=mpmath.mpf(1)
mp2=mpmath.mpf(2)
mp4=mpmath.mpf(4)
mppi=mpmath.mp.pi()
mpsqrtpi=mpmath.sqrt(mpmath.mp.pi())


class VVHarmonicWeakMaassForms(AutomorphicFormSpace):
    r"""
    Space of vector-valued harmonic weak Maass forms for the Weil representation of a finite quadratic module.

    """

    def __init__(self,WM,k=QQ(1)/QQ(2),holomorphic=False,dprec=15,sym=True,dual=False,verbose=0,**kwds):
        r"""
        Create a space of vector-valued harmonic weak Maass forms for the Weil representation of a finite quadratic module.

        INPUT:

            -''WM'' -- Weil representation (or other finite dimensional representation)

            -''k'' --  real

            -''dprec''-- integer (default 15)

            -''sym'' -- logical (default True)

            -''dr'' -- logical (default None) = True if we force usage of dual representation and False if we force standard rep.

        EXAMPLES::

            sage: WM=WeilRepDiscriminantForm(1)
            sage: M=VVHarmonicWeakMaassForms(WM,0.5,20);M
            Space of Vector-Valued harmonic weak Maass forms on Modular Group SL(2,Z) of weight 1/2  and dimension 0.
            Representation is Weil representation of the discriminant form given by ZZ/2ZZ with quadratic form Q(x)=1*x**2 mod 1.

        A ValueError is raised if we try to construct a space consisting of the zero function alone

            sage: WM=WeilRepDiscriminantForm(1,dual=True)
            sage: M=VVHarmonicWeakMaassForms(WM,0.5,20)
            ....
            ValueError: Space only contains the zero function! Change weight (1/2) or representation (dual rep.)

        """
        # If WR is an integer we construct the corr. Weil rep.
        if isinstance(WM,WeilRepMultiplier):
            self.WM=WM
        elif is_int(WM) or isinstance(WM,WeilRepDiscriminantForm):
            if is_int(WM):
                N=WM
            else:
                N = WM.N
            WR = WeilModule(N)
            self.WM=WeilRepMultiplier(WR,weight=QQ(k),dual=dual)
        else:
            raise ValueError("Got invalid multiplier: {0}".format(WM))
        self.group=MySubgroup(self.WM.group())
        self.weight=k
        self._weight_rat=QQ(RR(k))
        self.dprec=dprec
        self._pullback_prec = ceil(dprec*3.32)
        self._setupV_prec = ceil(dprec*3.32)
        self._verbose=verbose
        if not holomorphic:
            self._harmonic = True
            self._weak_holomorphic=False
        else:
            self._weak_holomorphic=True
        #self.prec=prec
        # if we want to force dual representation and try to use
        # a Weil representation WM which is not dual we change WM
        if isinstance(WM,WeilRepMultiplier):
            multiplier = WM
        else:
            multiplier = WeilRepMultiplier(WM,dual=dual,weight=QQ(k))
        self.group=multiplier.group()       
        self._is_dual_rep=multiplier.is_dual()
        
        AutomorphicFormSpace.__init__(self,self.group,weight=self.weight,multiplier=multiplier,holomorphic=holomorphic,weak=True,cuspidal=False,dprec=self.dprec,verbose=verbose)        

        self._index_set = multiplier.D()
        self._rank=multiplier.rank()
        if len(self._index_set)==0:
            if self._is_dual_rep:
                rep='dual rep.'
            else:
                rep='reg. rep.'
            raise ValueError("Space only contains the zero function! Change weight ({0}) or representation ({1})".format(self._weight_rat,rep))
        self.members=list() 
        #if hasattr(self.WM.WR
        #self._dimension_cusp_forms = 

    def index_set(self):
        return self._index_set

    def rank(self):
        return self._rank

    def __reduce__(self):
        r""" Used for pickling.


        EXAMPLES::

        sage: M=VVHarmonicWeakMaassForms(WR,0.5,100)
        sage: save(M,"M.sobj")

        """
        return(VVHarmonicWeakMaassForms,(self.multiplier(),self.weight,self.dprec,self._sym_type,self._is_dual_rep))

    def _cmp_(self,other):
        r""" Compare self to other"""
        if(not isinstance(other,VVHarmonicWeakMaassForms)):
            return False
        eq=(self.multiplier() == other.WR) and (self.weight_rat==other.weight_rat)
        eq = eq and (self.prec==other.prec) and (self._sym_type==other._sym_type)
        eq = eq and (self._is_dual_rep==other._is_dual_rep)
        return eq
    
    def _repr_(self):
        r""" Return string representation of self.

        EXAMPLES::
        
            sage: WR=WeilRepDiscriminantForm(1,dual=False)
            sage: M=VVHarmonicWeakMaassForms(WR,0.5,20);M
            Space of Vector-Valued harmonic weak Maass forms on Modular Group SL(2,Z) of weight 1/2  and dimension 1.
            Representation is Weil representation of the discriminant form given by ZZ/2ZZ with quadratic form Q(x)=1*x**2 mod 1.
            

        """
        s="Space of Vector-Valued harmonic weak Maass forms"
        s+=" on "+str(self.multiplier().group)+" of weight "+str(self._weight_rat)+" "
        s+=" and values in CC[ZZ/"+str(2*self.multiplier().N)+"ZZ]."
        s+="\nRepresentation is "+str(self.multiplier())
        return s

    def _latex_(self):
        r""" Return LaTeX string representation of self.

        EXAMPLES::
        
            sage: WR=WeilRepDiscriminantForm(1,dual=False)
            sage: M=VVHarmonicWeakMaassForms(WR,0.5,20)

            

        """
        p=self._weight_rat.numer()
        q=self._weight_rat.denom()
        old=s="\\begin{verbatim}\\end{verbatim}"
        new=""
  ##       s="\\text{Space of Vector-Valued harmonic weak Maass forms on }"
##         s+=latex(self.multiplier().group)+" \\text{ of weight } \\frac{"+str(p)+"}{"+str(q)+"}"
##         s+="\\text{and values in } \\mathbb{C}\\left[\\mathbb{Z}/"+latex(2*self.multiplier().N)+"\\mathbb{Z}\\right]\\text{.}"
##         s+="$ \\text{ The representation is }"+latex(self.multiplier())+"\\text{.}"
        s="\\begin{verbatim}\\end{verbatim}"
        s+=" Space of Vector-Valued harmonic weak Maass forms on $"
        s+=latex(self.multiplier().group)+"$  of weight  $\\frac{"+str(p)+"}{"+str(q)+"}$"
        s+="and values in $\\mathbb{C}\\left[\\mathbb{Z}/"+latex(2*self.multiplier().N)+"\\mathbb{Z}\\right]$. "
        s+="The representation is "+self.multiplier()._latex_().replace(old,new)+"."
        
        return s

    def sym_type(self):
        return self._sym_type
    
    def an_element(self):
        pp = self.smallest_pp()
        F = self.get_element(pp,maxC=10)
        return F


    def get_element(self,PP=None,prec=None,maxD=None,maxC=None,cusp_form=False,ef=True,mp=1E-8,M0_set=None,Y_set=None,SetCs={},use_mpmath=False,constant_term="pp"):
        r"""
        Get an element of the space of harmonic weak Maass form
        with specified principal part.

        INPUT:

        -''PP'' -- principal part
        -''prec'' -- number of digits desired
        -''maxD'' -- number of discriminants desired (use without prec, will override any value of prec and set it automatically)
        -''maxC'' -- number of coefficients desired for each component (use without prec, will set prec automatically)
        -''cusp_form'' -- True if we force c(h,n)=0 with n+q(h)=0. 

        -``use_mpmath`` -- use mpmath package for multiprecision instead of mpc/mpfr

        -``constant_term`` -- set to "pp" if the constant term is included in the principal part (default) and otherwise to "var" if constant terms are treated as variables.

        EXAMPLES::

            sage: WR=WeilRepDiscriminantForm(1,dual=False)
            sage: M=VVHarmonicWeakMaassForms(WR,0.5,20)
            sage: PP={(1/2,-1):1}
            sage: C=M.get_element(PP,12)
            sage: C[0][0][1]
            mpc(real='53503.999999999985', imag='-5.0861876730364994e-12')
            sage: abs(C[0][0][1]-mpmath.mpf(53504))
            mpf('1.5415172456347141e-11')

        """
        # There are three permitted formats for the principal part:  
        # P = dict or list of dicts of the form
        # i)  {(r/2N,m):c}
        # ii) {(r,m):c}
        # iii){D:c}
        # if the principal part contains
        # i)   c*e(m+Q(r/2N))
        # ii)  c*e(m+Q(r))
        #iii)  c*e(D/4N)
        minprec=mp
        #M0_set=None
        if PP==None:
            PP=self.smallest_pp() ## Use the smallest discriminant possible
        if prec == None and maxD == None and maxC == None:
            raise ValueError("Need either desired precision, number of discriminants or number of coefficients!")
        elif maxD != None:
            # Need the correct M0 to use for getting this many coefficients. 
            prec = None
            M0_set=ceil(maxD/self.rank()) # should be approximately good
            for m in range(1,maxD*self.rank()+1):
                min_d=maxD
                ## We need to make sure that all discriminants below maxD are accounted for
                for r in self.index_set():
                    D=D_from_rn(self.multiplier(),(r,m))
                    if D < min_d:
                        min_d=D
                if min_d >= maxD:
                    M0_set=m
                    break
        elif maxC != None:
            prec=None
            M0_set=maxC
            
        if not isinstance(PP,(dict,list)):
            raise TypeError("Need principal part in form of dictionary! Got:{0}".format(PP))
        if isinstance(PP,dict):
            Ptmp = [PP]
        else:
            Ptmp = PP
        type='vector'
        P0 = list()
        for p in Ptmp:
            d={}
            for t in p:
                if isinstance(t,tuple):
                    (r,m) = t
                if isinstance(t,(int,Integer)):
                    if not self.is_Heegner_disc(t):
                        raise ValueError("Need discriminant satisfying Heegner condition, i.e. square mod {0}. Got: {1}".format(self.multiplier().level(),t))
                    (r,m) = rn_from_D(self.multiplier(),t)
                d[(r,m)]=p[t]
                # Also check that principal part adheres to the symmetry if present
                if self._sym_type!=0:
                    minus_r = self.multiplier().weil_module().neg_index(r)
                    if (minus_r,m) in p:
                        if p[(minus_r,m)]!=self.sym_type*p[(r,m)]:
                            raise ValueError("Need symmetric principal part! Got:{0}".format(PP))
            P0.append(d)
        P=P0[0]  ### More than one principal part is not implemented yet...
        if self._verbose > 0:
            print("P={0}".format(P))
            print("M0_set={0}".format(M0_set))
        F=VVHarmonicWeakMaassFormElement(self,P)
        if self._verbose>0:
            sys.stdout.flush()
        if M0_set != None and M0_set != 0:
            if self._verbose > 0:
                print("M0_set={0}".format(M0_set))
            mpmath.mp.dps = 53 ## Do the estimate in low precision
            Y=mpmath.mpf(0.75); M=M0_set
            [er1,er2]=F.get_error_estimates(Y,M0_set)
            prec=max(er1,er2)
        elif prec != None:
            [Y,M]=self.get_Y_and_M(P,self.weight,prec)
        else:
            raise ValueError("Could not deicide number of coefficients to compute from input!")
        if minprec != None and prec > minprec:
            prec=minprec
            [Y,M]=self.get_Y_and_M(P,self.weight,minprec)
        Q=M+50
        if(Y_set!=None):  
            if Y_set>0 and Y_set < 0.866:
                Y=Y_set
        if self._verbose > 1:
            print("prec=",prec)
            print("Y=",Y)
            print("M=",M)
            print("Q=",Q)
        dold=mpmath.mp.dps
        mpmath.mp.dps=max(self._dprec,prec+10)
        if self._verbose > 0:
            print("using {0} digits".format(mpmath.mp.dps))
            print("P={0}".format(P))
            print("ef={0}".format(ef))
	RF = RealField(self._prec)
        if(ef):
            if use_mpmath:
                Y = mpmath.mp.mpf(Y)
                W=vv_harmonic_wmwf_setupV_ef(self,P,Y,M,Q,mpmath.mp.mpf(self.weight))
            else:
                Y = RF(Y)
                mpmath.mp.dps = self._setupV_prec
                W=vv_harmonic_wmwf_setupV_mpc2(self,P,Y,M,Q)
                mpmath.mp.dps = self._prec
        else:
            W=vv_harmonic_wmwf_setupV(self,P,Y,M,Q,self.weight,self._sym_type,verbose=self._verbose)
        W['space']=self
        W['PP']=P
        s="tmpWN"+str(self.multiplier().N)+"-"
        dirv=["/local/stroemberg/tmp/","/tmp","."]
        try:
            for dirs in dirv:
                try:
                    st=os.stat(dirs)
                    [f,tmpfilename] = tempfile.mkstemp(prefix=s,dir=dirs,suffix='.sobj')
                    os.close(f)
                    raise StopIteration()
                except OSError:
                    continue
        except StopIteration:
            pass
        if self._verbose > 0:
            print("tmpfilename={0}".format(tmpfilename))
        try:
            save(W,tmpfilename)
        except MemoryError:
            print("Could not save to file!")
            pass
        #if(PP.has_key((0,0))):
        #    N=self.set_norm(P=PP,C=SetCs)
        #else:
        N=self.set_norm(P=PP,C=SetCs)
        if self._verbose>0:
            print("N = {0}".format(N))
            print("V.parent={0}".format(W['V'].parent()))
        #return (W,N)
        if use_mpmath:
            C=solve_system_for_vv_harmonic_weak_Maass_waveforms(W,N,deb=False)
        else:
            C=solve_system_for_vv_harmonic_weak_Maass_waveforms_new(self,W,N)
        if self._verbose > 0:
            print("C001={0}".format(C[0][0][1]))
        D=self.D
        CC=dict()
        # we need to shift the indices of C to the (semi-)correct indices
        # i.e. for effficiency we use integers instead of rationals
	#return C
        F=list()
        for j in C.keys():
            CC[j]=dict()
            for k in C[j].keys():
                #print "k=",k
                CC[j][self.D_as_int[k]]=C[j][k]
            F.append(VVHarmonicWeakMaassFormElement(self,P,CC[j],prec))
        mpmath.mp.dps=dold
        if(len(F)==1):
            return F[0]
        else:
            return F
        #return P
        # nw we have a 

    def smallest_pp(self,sgn=-1,n=1):
        r""" Returns the smallest valid principal part.
        INPUT:
        -''sgn'' -- integer, if we want a negative or a positive power of Q in the principal part
                      (technically speaking it is not really the principal part ifthe power is > 0)
        -''n''  -- integer. seek the n-th smallest discriminant 
        """
        res=list()
        sz=ceil(QQ(n)/QQ(len(self.D)))+1
        if self._is_dual_rep:
            m_start=0
        else:
            m_start=min(sgn,0)
        if(n>1):
            m_stop=sz
        else:
            m_stop=m_start+1
        rmin=1; xmin=self.D_as_int[0]
        if self._verbose > 0:
            print("m_start={0}".format(m_start))
        for m in range(m_start,m_stop):
            #m=m*sgn
            for x in self.multiplier().Qv:
                y=QQ(x+m)
                D=ZZ(y*self.multiplier().level())
                if self._verbose > 0:
                    print(m,self.multiplier().Qv.index(x),D,y)
                if D==0 and sgn!=0: # don't count the zero unless we want zero
                    continue
                if res.count(D) == 0:
                    res.append(D)
                if abs(y) < xmin and n == 1:
                    j = self.multiplier().Qv.index(x)
                    if j in self.D(): #_as_int):
                        rmin=j
                        xmin=x                
            if(len(res)>=n):
                break
        if(sgn==-1):
            res.sort(reverse=True)
        else:
            res.sort(reverse=True)
        if self._verbose > 0:
            print(res)
        # make unique
        return {res[n-1]:1}
                
    def next_heegner_disc(self,n=0,sgn=-1,fd=False):
        r""" Returns the smallest (in absolute value) discriminant greater than n which satisfies the Heegner condition,
        i.e. which appears as an index of Fourier coefficients of forms in M.
        INPUT:
        -'n' -- integer (default 0), look at discriminants |D| > n
        -'sgn' -- integer (default -1), look at negative discriminants with sign = sgn
        -'fd'' -- logical(default False), if True, look at fundamental discriminats        
        """
        for D in range(n+1,2*n+100):
            DD=D*sgn
            if(self.is_Heegner_disc(DD)):
                return DD
        raise ArithmeticError(" COuld not find any Heegner discriminat > {0} !".format(n))
            

    def is_Heegner_disc(self,D):
        r""" Returns true is \pm D is appears as an index of a Fourier coefficients of M,
        i.e. in the simplest case if it is a square mod 4N
        INPUT:
        -''D'' -- integer
        """ 
        Dr = D % self.multiplier().level()
        if self.multiplier().is_dual():
            Dr = -Dr
        if Dr in self.multiplier().Qv_times_level:
            return True
        else:
            return False

    def get_Y_and_M(self,PP,weight,prec,Yin=None):
        r"""
        Find a good Y and M for computing coefficents with precison 10^-prec
        
        """
        # generalized_level
        if(Yin!=None):
            Y0=Yin
        else:
            Y0=min(self.group.minimal_height(),0.5)
        Cmax=1
        Kmax=0
        Cmax=max(PP.values())
        for t in PP.keys():
            if isinstance(t,tuple):
                (c,l) = t
            elif isinstance(t,(int,Integer)):
                (c,l)=rn_from_D(self.multiplier(),t)
            else:
                raise ValueError("Incorrect principal part: t={0}".format(t))
            if c in self.multiplier().D:
                tmp=l+self.multiplier().Qv[self.multiplier().D.index(c)]
            elif c in range(len(self.multiplier().Qv)):
                tmp=l+self.multiplier().Qv[c]
            else:
                    raise ValueError("Incorrect principal part: c,l={0},{1}".format(c,l))
            if self._verbose>0:
                print("tmp={0}".format(tmp))
            if abs(tmp)>Kmax:
                Kmax=abs(tmp)
            #x
        # then get corresponding M
        #print "Kmax=",Kmax
        #print "Cmax=",Cmax
        M0=self.get_M(Y0,Kmax,Cmax,prec)
        return [Y0,M0]
    
    def get_M(self,Y,K0,K1,prec):
        r""" Computes truncation point for Harmonic Maass waveform.
        """
        
        # # Use low precision
        dold=mpmath.mp.dps
        # mpmath.mp.dps=int(mpmath.ceil(abs(mpmath.log10(eps))))+5
        mpmath.mp.dps=max(dold,prec)+5
        twopi=2*mpmath.pi()
        twopiy=twopi*mpmath.mpf(Y)
        # an extra two for the accumulation of errors
        eps=mpmath.mpf(10)**mpmath.mpf(-prec)
        minm=max(10,abs(int(1-self.weight)+1)/(2*mpmath.pi()*Y))
        #print "K0=",K0
        #print "K1=",K1
        [Cp0,Cp1]=self.get_Cp(K0)
        Cm=self.get_Cm(K0,K1)
        #print "Cp0,Cp1,Cm=",mppr(Cp0),mppr(Cp1),mppr(Cm)
        fak=len(self.multiplier().D)
        try:
            for m in range(minm,minm+10000):
                errest1=fak*self.err_est_vv_hwmf_pos(Y,m,Cp0,Cp1)
                errest2=fak*self.err_est_vv_hwmf_neg(Y,m,Cm)
                #print "er1+(",m,")=",mppr(errest1)
                #print "er2-(",m,")=",mppr(errest2)
                if(max(abs(errest1),abs(errest2))<eps):
                    raise StopIteration()
            raise ArithmeticError("Could not find M<%s such that error bound in truncation is <{0}! and Y,K0,K1={1},{2},{3} \n err+={4} \n err-={5}".format(m,eps,mppr(Y),K0,K1,mppr(errest1),mppr(errest2)))
        except StopIteration:
            if(self._verbose > 2):
                print("er +={0}".format(errest1))
                print("er -={0}".format(errest2))
                print("m={0}".format(m))
                print("Y={0}".format(Y))
        mpmath.mp.dps=dold
        return m


    def get_sym_type(self):
        r"""
        Calculate the symmetry type (even/odd) for the combination
        of representation and weight.
        """
        
        t=self.weight-0.5
        if(not is_int(t)):
            raise ValueError("Need half-integral value of weight! Got k={0}, k-1/2={1}".format(self.weight,t))
        ti=Integer(float(t))
        if is_odd(ti):
            sym_type=-1
        else:
            sym_type=1
        if self._is_dual_rep:
            sym_type=-sym_type
        return sym_type

    

    def get_Cp(self,K0):
        r"""
        Set constants for the error estimates.
        """
        #if(self.weight>=1.5):
        #    raise ValueError," Error bounds only accurate for k<1.5! got k=%s" % self.weight
        mp2=mpmath.mpf(2)
        twominusk=mp2-self.weight        
        tmp=mpmath.mpf(len(self.multiplier().D))
        tmp0=mpmath.sqrt(tmp)+mpmath.mpf(1)
        tmp1=mpmath.pi()*mpmath.mpf(4)
        Cp1=tmp1*mpmath.sqrt(abs(K0))
        tmp1=mpmath.power(tmp1,twominusk)
        tmp2=mpmath.besseli(1-self.weight,1.0)
        tmp3=mpmath.zeta(twominusk)
        if(K0==0):
            tmp4=1
        else:
            tmp4=mpmath.power(K0,1-self.weight)
        Cp0=tmp0*tmp1*tmp2*tmp3*tmp4
        return [Cp0,Cp1]

    def get_Cm(self,K0,K1):
        r""" Constant in error bound for negative part.
        """
        #if(self.weight>=1.5):
        #    raise ValueError," Error bounds only accurate for k<1.5! got k=%s" % self.weight
        twominusk=mp2-self.weight
        tmp=mpmath.mpf(len(self.multiplier().D))
        tmp1=mppi*mp2
        tmp1=mpmath.power(tmp1,twominusk)
        tmp3=mpmath.zeta(twominusk)
        if(K0==0):
            tmp4=1
        else:
            tmp4=mpmath.power(K0,1-self.weight)
        g1=mpmath.gamma(1-self.weight)
        g2=mpmath.gamma(2-self.weight)

        Cm=mp2/g1+mp4*tmp1/g1/g2*tmp*tmp3*tmp4
        return Cm


    def err_est_vv_hwmf_pos(self,Y,m,Cp0,Cp1):
        r""" Error estimate. See paper...
        """
        #if(self.weight>=1.5):
        #    raise ValueError," Error bounds only accurate for k<1.5! got k=%s" % weight


        twopiY=mpmath.pi()*mp2*Y
        fourpiY=mp2*twopiY
        tmp=mpmath.mpf(len(self.multiplier().D))        
        etmp1=mpmath.sqrt(mpmath.mpf(m))-Cp1/fourpiY
        etmp2=mpmath.exp(-twopiY*etmp1**2)
        etmp3=mp2+mpsqrtpi*Cp1/mp2/mpmath.sqrt(twopiY)
        etmp4=Cp0/twopiY
        err_pos=tmp*etmp4*etmp3*etmp2
        return err_pos
        

    def err_est_vv_hwmf_neg(self,Y,m,Cm):
        r""" Errorbound for negative part.
        """
        #print "Cm=",Cm
        #if(self.weight>=1.5):
        #    raise ValueError," Error bounds only accurate for k<1.5! got k=%s" % weight
        etmp1=abs(mpmath.mpf(1-self.weight))
        twopiY=mp2*Y*mpmath.pi()
        tmp=mpmath.mpf(len(self.multiplier().D))
        if(self.weight>0):
            etmp2=mp2*mpmath.mpf(m-1)*twopiY
            etmp2=mpmath.power(etmp2,-self.weight)
            etmp3=mpmath.exp(-twopiY*mpmath.mpf(m))
            etmp4=mp1/(mp1-mpmath.exp(-twopiY))
            err_neg=Cm*tmp*etmp1*etmp2*etmp3*etmp4
        else:
            etmp2=mpmath.power(mp2,-self.weight)
            etmp3=mpmath.power(twopiY,-self.weight-mp1)
            etmp4=mpmath.power(mpmath.mpf(m),-self.weight)
            etmp5=mpmath.exp(-twopiY*mpmath.mpf(m))
            err_neg=Cm*tmp*etmp1**2*etmp2*etmp3*etmp4*etmp5
        return err_neg


    #def set_norm_vv_harmonic_weak_maass_forms(WR,cusp_form=True,holomorphic=True,SetCs=None):
    def set_norm(self,P={},C={},c_t="pp"):
        r"""
        Set the normalization dictionary corresponding to self and computation of a form
        with principal part P and set fourier coefficients in C
        """
        N=dict()
        if isinstance(P,list):
            Pl=P
        else:
            Pl=[P]
        if isinstance(C,list):
            Cl=C
        else:
            Cl=[C]        
        if len(Pl)>0:
            N['comp_dim']=len(Pl)
        if len(Cl) > 0:
            if len(Cl) != len(Pl):
                raise ValueError("Need same number of principal parts and coefficients to set!")
            keys = list(Cl[0].keys())
            for j in range(1,N['comp_dim']):
                if list(Cl[j].keys()) != keys:
                    raise ValueError("Need to set the same coefficients! (or call the method more than once)")
        else:
            Cl=[]
            for j in range(N['comp_dim']):
                Cl.append(C)
        N['Vals']=list()
        N['Vals']=list()
        N['SetCs']=list()
        for i in range(N['comp_dim']):
            N['Vals'].append({})
            N['Vals'].append({})
            N['SetCs'].append([])
            for j in range(len(self.multiplier().D)):
                a=self.multiplier().D[j]
                x=self.multiplier().Qv[j]
                #N['Vals'][i][(0,j)]=dict()
                if x==0:
                    if c_t=="pp":
                        #N['comp_dim']=N['comp_dim']+1
                        N['SetCs'][i].append((j,0))
                        if (0,j) in Pl[i]:
                            N['Vals'][i][(j,0)]=Pl[i][(j,0)]
                        else:
                            N['Vals'][i][(j,0)]=0 #P[(0,0)]
                elif x<0 and self._holomorphic:
                    N['SetCs'][i].append((j,0))
                    N['Vals'][i][(j,0)]=0 #P[(0,0)]            
            
            for (r,n) in Cl[i].keys():
                if(N['SetCs'][i].count((r,n))==0):
                    N['SetCs'][i].append((r,n))
                    N['Vals'][i][(r,n)]=Cl[i][(r,n)] 
        return N







    # def rn_from_D(self,D):
    #     r""" Find the pair(s) (r,n) s.t. +-D/4N= n +- q(r) for D in D 

    #     INPUT:
    #     -''D'' -- integer or list of integers

    #     OUTPUT:
    #     -''t'' -- tuple (r,n) or list of tuples



    #     """
    #     if(isinstance(D,list)):
    #         lout=list()
    #         for DD in D: 
    #             t=self._one_rn_from_D(DD)
    #             if(t<>None):
    #                 lout.append(t)
    #         return lout
    #     else:
    #         return self._one_rn_from_D(D)


    # def _one_rn_from_D(self,D):
    #     r""" Find the (r,n) s.t. +-D/4N= n +- q(r)
    #     """            
    #     Dv=QQ(D)/QQ(self.multiplier().level())
    #     sig=1
    #     if self.multiplier()._is_dual_rep:
    #         sig=-1
    #     for r in self.multiplier().D:
    #         x=self.multiplier().Q(r)
    #         if(is_int(Dv-x)):
    #             rr=self.multiplier().D.index(r)
    #             n=sig*int(Dv-x)
    #             return (rr,n)
    #     return None


    # def D_from_rn(self,t):
    #     r""" Find the D s.t. +-D/4N= n +- q(r)
    #     """
    #     if(isinstance(t,list)):
    #         lout=list()
    #         for (r,n) in t:
    #             D=self._one_D_from_rn((r,n))
    #             if(D<>None):
    #                 lout.append(D)
    #         return lout
    #     else:
    #         return self._one_D_from_rn(t)



    # def _one_D_from_rn(self,t):
    #     r""" Find the D s.t. +-D/4N= n +- q(r)
    #     """
    #     #print "t=",t,type(t)
    #     if(not isinstance(t,tuple)):
    #         raise TypeError,"Need a tuple of integers! Got:%s" % t
    #     (r,n)=t
    #     sig=1
    #     #print "r=",r
    #     if(r in self.multiplier().D):
    #         x=self.multiplier().Q(r)
    #     elif(r in self.multiplier().D_as_integers):
    #         x=self.multiplier().Q(self.multiplier().D[r])
    #     else:
    #         raise TypeError,"Need (r,n) in proper format forcoefficients! I.e. n integer and r in D or integer!"
    #     #print "x=",x
    #     if self.multiplier()._is_dual_rep:
    #         sig=-1
    #     D=sig*self.multiplier().level()*(n+sig*x)
    #     return D


def dist_from_int(x):
    r"""
    Return the distance to the closest integer to x, as well as the closest integer
    """
    m1=floor(x); m2=ceil(x)
    er1=abs(m1-x); er2=abs(m2-x)
    #print "er1=",er1
    #print "er2=",er2
    if(er1<er2):
        return [er1,m1]
    elif(er2<er1):
        return [er2,m2]
    else:
        random.seed()
        d=random.getrandbits(1)
        #print "d=",d
        if(d==0):
            return [er1,m1]
        else:
            return [er2,m2]



def norm_sci_pretty_print(c,nd=0,es='e',latex_pow=False):
    x=c.real; y=c.imag
    if(abs(x)<1E-5):
        s = sci_pretty_print(x,2,es,latex_pow)
    else:
        s = sci_pretty_print(x,nd,es,latex_pow)
    if(y>0):
        p="+"
    else:
        p=""
    if(abs(y)<1E-5):
        s=s+p+sci_pretty_print(y,2,es,latex_pow)+"i"
    else:
        s=s+p+sci_pretty_print(y,nd,es,latex_pow)+"i"
    return s

def sci_pretty_print(s,nd=0,es='e',latex_pow=False):
    r""" Take a string representation of a number and returns it in scientific notation to desired number of digits.
    """
    # remove leading sign #
    #x=mpmath.mp.mpf(s)
    #if(abs(x)<1):
    #    raise NotImplementedError," Only implemented for |x|>1!!  got x=%s" %x
    if(not isinstance(s,str)):
        s=str(s)
        
    s=s.replace("(","")
    s=s.replace(")","")
    s=s.strip()
    if s.count("I")+s.count("i")+s.count("j") > 0:
        # Get a default complex notation
        s=s.replace("*","")
        s=s.replace("I","i")
        s=s.replace("j","i")
        ## We have to find imaginary and real parts
        l=s.split("+")                
        if(len(l)>1):
            (s1,s2)=l
        else:
            (s1,s2)=s.split("-")                
            s2="-"+s2
        if(s1.count("i")>0): # put imaginary part in standard form
            ims=s1.strip("i").lstrip(); res=s2             
        else:
            ims=s2.strip("i").lstrip(); res=s1 
        if(ims==""): ims="1"
        sres=sci_pretty_print(res,nd,es,latex_pow)    
        sims=sci_pretty_print(ims,nd,es,latex_pow)    
        if(sres=="0"): sres=""
        if(sims=="0"):
            sims=""
        else:
            sims=sims+"i"

        if sims.count("-") > 0:
            return sres+" "+sims.replace(" -"," - ")
        elif sims != "" and sres!="":
            return sres+" + "+sims
        elif sres != "":
            return sres
        elif sims != "":
            return sims
        else:
            raise ValueError("Could not find pretty print for s={0} ".format(s))
    s=s.strip()    
    if len(s.replace(".","").strip("0")) == 0:
        return "0"
    if s.count(".") == 0:
        s=s+".0"
    if s[0] == '-':
        ss=s.strip("-")
        ss=sci_pretty_print(ss,nd,es,latex_pow)    
        return "-"+ss

    l=s.split(".")
    if(len(l)>1):
        (sint,sdigs)=l
    elif(len(s)<nd):
        return s
    elif(len(l)>0):
        sint=l[0]
        sdigs=""
    else:
        raise ValueError(" Can not do pretty print for s={0}".format(s))

    if sdigs.count("e") > 0:
        l=sdigs.split("e")
        sdigs=l[0]
        ex=int(l[1])
    else:
        ex=0
    if len(sint) == 1 and sint == "0":
        # find the number of leading zeros
        sss=sdigs.lstrip("0")        
        nz=len(sdigs)-len(sss)+1
        if(nz<10):
            ex="-0"+str(nz)
        else:
            ex=str(-nz)
        # Fix correct rounding
        rest=sss[nd:len(sss)]
        ix=nd-1
        if len(rest)>0:
            if int(rest) < 5*10**(len(rest)-1):
                # print " round < : since "+rest+"<"+str(5*10**(len(rest)-1))
                d=int(sss[ix])
            elif int(rest) > 5*10**(len(rest)-1):
                # print " round > : since "+rest+">"+str(5*10**(len(rest)-1))
                d=int(sss[ix])+1
            else:
                # if we have an exact half we round randomly 
                random.seed()
                d=int(sss[ix])+int(random.getrandbits(1))            
            if(d<10):
                ssdigs=sss[0:ix]+str(d) # We account for the leading digit too
            else:
                ssdigs=sss[0:ix-1]+str(int(sss[ix-1])+1)+str(d-10) # We account for the leading digit too
            if(latex_pow):
                return ssdigs[0]+"."+ssdigs[1:nd]+es+"10^{"+ex+"}"
            else:
                return ssdigs[0]+"."+ssdigs[1:nd]+es+ex            
    ex=int(ex)+len(sint)-1
    if(abs(ex)<10):
        ex="0"+str(ex)
    else:
        ex=str(ex)
    #ssdigs=sint[1:len(sint)]+sdigs
    # cut away to nd digits
    if(nd>0):
        #ssdigs=sdigs[0:nd-1] # We acount the leading digit too
        # Try to do correct rounding        
        rest=sdigs[nd-len(sint):len(sdigs)]
        #print "sdigs=",sdigs," nd=",nd
        #print "rest=",rest
        ix=nd-len(sint)-1
        if(len(rest)>0):
            if(int(rest) < 5*10**(len(rest)-1)):
                # print " round < : since "+rest+"<"+str(5*10**(len(rest)-1))
                d=int(sdigs[ix])
            elif(int(rest) > 5*10**(len(rest)-1)):
                # print " round > : since "+rest+">"+str(5*10**(len(rest)-1))
                d=int(sdigs[ix])+1
            else:
                # if we have an exact half we round randomly 
                random.seed()
                d=int(sdigs[ix])+int(random.getrandbits(1))            
            if(d<10):
                ssdigs=sdigs[0:ix]+str(d) # We account for the leading digit too
            else:
                if(ix>0):
                    ssdigs=sdigs[0:ix-1]+str(int(sdigs[ix-1])+1)+str(d-10) # We account for the leading digit too
                else:
                    ssdigs=str(d-10) # We account for the leading digit too
                    if(len(sint)==1):
                        sint=str(int(sint)+1)
                    else:
                        ll=len(sint); stmp=sint;
                        sint=stmp[1:ll-1]+str(int(stmp[ll-1])+1)
        else:
            ssdigs=sdigs[0:ix] 
        #print "rest=",rest,len(rest)
        ssdigs=sint[1:len(sint)]+ssdigs    
    else:
        ssdigs=sint[1:len(sint)]+sdigs    

    if latex_pow:
        res=sint[0]+"."+ssdigs+" \cdot 10^{"+ex+"}"
        return res 
    else:
        return sint[0]+"."+ssdigs+es+ex

## def dist_from_int(x):
##     r""" Compute distance from x to the nearest integer.
##     """
##     d1=abs(x-mpmath.mp.floor(x))
##     d2=abs(x-mpmath.mp.ceil(x))
##     if(d1<d2):
##         y=mpmath.mp.floor(x)
##         d=d1
##     else:
##         y=mpmath.mp.ceil(x)
##         d=d2
##     return RR(d),y

## def is_int(q):
##     r"""
##     Find out if the rational number q is an integer.
##     """
##     try:
##         if(isinstance(q,sage.rings.integer.Integer) or isinstance(q,int)):
##             return True
##         if(isinstance(q,sage.rings.rational.Rational)):
##             n=q.denominator()
##             if(n==1):
##                 return True
##         if(floor(q)==ceil(q)):
##             return True
##     except TypeError:
##         pass
##     return False
        

def solve_system_for_vv_harmonic_weak_Maass_waveforms(W,N=None,deb=False,gr=False,cn=False):
    r"""
    Solve the linear system to obtain the Fourier coefficients of Maass forms

    INPUT:

    - ``W`` --   (system) dictionary
        - ``W['Ms']``  -- M start
        - ``W['Mf']``  -- M stop
        - ``W['nc']``  -- number of cusps
        - ``W['V']``   -- matrix of size ((Ms-Mf+1)*nc)**2
        - ``W['RHS']`` -- right hand side (for inhomogeneous system) matrix of size ((Ms-Mf+1)*nc)*(dim)
    - ``N`` -- normalisation (dictionary, output from the set_norm_for_maass function, default None)
        - if N=None we assume that the solution is uniquely determined by the prinicpal part (in the right hand side)
        - ``N['SetCs']``   -- Which coefficients are set
        - ``N['Vals'] ``   -- To which values are these coefficients set
        - ``N['comp_dim']``-- How large is the assumed dimension of the solution space
        - ``N['num_set']`` -- Number of coefficients which are set
        

        
    - ``deb`` -- print debugging information (default False)
    - ''gr''  -- only return the reduced matrix and right hand side. do not perform the solving .
    - ''cn''  -- logical (default False) set to True to compute the max norm of V^-1
    OUTPUT:
    
    - ``C`` -- Fourier coefficients

    EXAMPLES::

        sage: G=MySubgroup(Gamma0(1))
        sage: mpmath.mp.dps=20
        sage: R=mpmath.mpf(9.533695261353557554344235235928770323821256395107251982375790464135348991298347781769255509975435366)
        sage: Y=mpmath.mpf(0.5)
        sage: W=setup_matrix_for_Maass_waveforms(G,R,Y,12,22)
        sage: N=set_norm_maass(1)
        sage: C=solve_system_for_Maass_waveforms(W,N)
        sage: C[0][2]*C[0][3]-C[0][6]
        mpc(real='-1.8055426724989656270259e-14', imag='1.6658248366482944572967e-19')

    If M is too large and the precision is not high enough the matrix might be numerically singular

        W=setup_matrix_for_Maass_waveforms(G,R,Y,20,40)  
        sage: C=solve_system_for_Maass_waveforms(W,N)
        Traceback (most recent call last)
        ...
        ZeroDivisionError: Need higher precision! Use > 23 digits!

    Increasing the precision helps
    
        sage: mpmath.mp.dps=25
        sage: R=mpmath.mpf(9.533695261353557554344235235928770323821256395107251982375790464135348991298347781769255509975435366)
        sage: C=solve_system_for_Maass_waveforms(W,N)
        sage: C[0][2]*C[0][3]-C[0][6]
        mpc(real='3.780824715556976438911480324e-25', imag='2.114746048869188750991752872e-99')

        
        

    """
    V=W['V']
    Ms=W['Ms']
    Mf=W['Mf']
    H = W['space']
    WR= H.multiplier()
    #nc=W['nc']
    Ml=Mf-Ms+1
    get_reduced_matrix=gr
    comp_norm=cn
    if(W['sym_type']==1):
        setD=list(range(0,WR.N+1))  # 0,1,...,N
    elif(W['sym_type']==-1):
        setD=list(range(1,WR.N))    # 1,2,...,N-1  (since -0=0 and -N=N)
    nc=len(setD)
    if V.cols != Ml*nc or V.rows != Ml*nc:
        raise Exception(" Wrong dimension of input matrix!")
    if N == None:
        SetCs=[]
        Vals=dict();Vals[0]=dict()
        for b in WR.D:
            if(WR.Q(b)==0):
                Vals[0][b]=dict()
                SetCs.append((b,0))
                Vals[0][b][0]=0
        comp_dim=1
    else:
        SetCs=N['SetCs']
        Vals=N['Vals']
        comp_dim=N['comp_dim']
    #for a in Vals.keys():
    #    Vals[a]=mpmath.mp.mpc(Vals[a])
    if(Ms<0):
        use_sym=0
    else:
        use_sym=1
    if(use_sym==1 and len(SetCs)>0):
        num_set=len(SetCs)-1
    else:
        num_set=len(SetCs)
    t=V[0,0]
    if(isinstance(t,float)):
        mpmath_ctx=mpmath.fp
    else:  
        mpmath_ctx=mpmath.mp
    #print "mpmath_ctx=",mpmath_ctx
    #use_symmetry=False
    RHS=mpmath_ctx.matrix(int(Ml*nc-num_set),int(comp_dim))
    if 'RHS' in W:
        if W['RHS'].cols!=comp_dim:
            raise ValueError("Incorrect number of right hand sides!")

    LHS=mpmath_ctx.matrix(int(Ml*nc-num_set),int(Ml*nc-num_set))
    roffs=0
    if(deb):
        print("Ml={0}".format(Ml))
        print("num_set={0}".format(num_set))
        print("SetCs={0}".format(SetCs))
        print("Vals={0}".format(Vals))
        print("V.rows={0}".format(V.rows))
        print("V.cols={0}".format(V.cols))
        print("LHS.rows={0}".format(LHS.rows))
        print("LHS.cols={0}".format(LHS.cols))
        print("RHS.rows={0}".format(RHS.rows))
        print("RHS.cols={0}".format(RHS.cols))
        print("use_sym={0}".format(use_sym))

    if V.rows != nc*Ml:
        raise ArithmeticError(" Matrix does not have correct size!")
    if len(SetCs)>0:
        for a in range(nc):
            for n in range(Ms,Mf+1):
                r=a*Ml+n-Ms
                if SetCs.count((a,n))>0:
                    print(" set row a,n={0}, {1}".format(a,n))
                    roffs=roffs+1
                    continue
                for fn_j in range(comp_dim):
                    if 'RHS' in W:
                        RHS[r-roffs,fn_j]=-W['RHS'][r,fn_j]
                    for (i,cset) in SetCs:
                        v=Vals[fn_j][i][cset]
                        if mpmath_ctx == mpmath.mp:
                            tmp=mpmath_ctx.mpmathify(v)
                        elif isinstance(v,float):
                            tmp=mpmath_ctx.mpf(v)
                        else:
                            tmp=mpmath_ctx.mpc(v)
                        tmp=tmp*V[r,i*Ml+cset-Ms]
                    RHS[r-roffs,fn_j]=RHS[r-roffs,fn_j]-tmp
                coffs=0
                for b in range(nc):
                    for l in range(Ms,Mf+1):
                        k=b*Ml+l-Ms
                        if(SetCs.count((b,l))>0):
                            #print " set col b,l=",b,l
                            coffs=coffs+1
                            continue
                        LHS[r-roffs,k-coffs]=V[r,k]
            #print "LHS[",r,k,"]=",LHS[r-roffs,k-coffs]
    else:
        LHS=V
        RHS=-W['RHS']
    if(get_reduced_matrix):
        return [LHS,RHS]
    dpold=mpmath.mp.dps    
    maxit=100;i=0
    done=False
    while (not done and i<=maxit):
        try:
            A, p = mpmath_ctx.LU_decomp(LHS)
            done=True
        except ZeroDivisionError:
            try:
                sinf=smallest_inf_norm_mpmath(LHS)
                t=int(mpmath_ctx.ceil(-mpmath_ctx.log10(sinf)))
            except ValueError:
                print("Warning: Got smallest inf. norm={0}".format(sinf))
                t = mpmath.mp.dps+10
            mpmath.mp.dps=t+5*i; i=i+1
            print("raising number of digits to:{0}".format(mpmath.mp.dps))
            # raise ZeroDivisionError,"Need higher precision! Use > %s digits!" % t
    if(i>=maxit):
        raise ZeroDivisionError("Can not raise precision enough to solve system! Should need > {0} digits! and {1} digits was not enough!".format(t,mpmath.mp.dps))
    # Use the LU-decomposition to compute the inf- norm of A^-1
    # Note that Ax=LUx=y and we get columns of A^-1 by solving for y1=(1,0,...), y2=(0,1,0,...) etc.
    if(comp_norm):
        max_norm=0
        for j in range(LHS.rows):
            y=mpmath_ctx.matrix(LHS.rows,int(1)); y[j,0]=1            
            b = mpmath_ctx.L_solve(A,y, p)
            TMP = mpmath_ctx.U_solve(A, b)
            tmpnorm=max(list(map(abs,TMP)))
            if(tmpnorm>max_norm):
                max_norm=tmpnorm
        print("max norm of V^-1={0}".format(max_norm))
    mpmath.mp.dps=dpold
    X=dict()
    for fn_j in range(comp_dim):
        X[fn_j] = dict() #mpmath.matrix(int(Ml),int(1))
        print("len(B)={0}".format(len(RHS.column(fn_j))))
        b = mpmath_ctx.L_solve(A, RHS.column(fn_j), p)
        #return b
        TMP = mpmath_ctx.U_solve(A, b)
        roffs=0
        res = mpmath_ctx.norm(mpmath_ctx.residual(LHS, TMP, RHS.column(fn_j)))
        print("res({0}={1}".format(fn_j,res))
        for ai in range(nc):
            X[fn_j][ai]=dict()
            for n in range(Ms,Mf+1):
                ni=Ml*ai+n-Ms
                if(SetCs.count((ai,n))>0):
                    roffs=roffs+1
                    # print "X[",fn_j,",",n,",Vals[fn_j][n]
                    # print "X[",fn_j,",",n,",Vals[fn_j][n]
                    X[fn_j][ai][n]=mpmath_ctx.mpc(Vals[fn_j][ai][n])
                    continue
                #print "roffs,n,ni=",roffs,n,ni-roffs
                #print "TMP=",TMP[ni-roffs,0]
                X[fn_j][ai][n]=TMP[ni-roffs,0]
    # return x
    return X

def solve_system_for_vv_harmonic_weak_Maass_waveforms_new(H,W,N=None,gr=False,cn=False,verbose=0):
    r"""
    Solve the linear system to obtain the Fourier coefficients of Maass forms

    INPUT:

    - ``H`` -- Space of vector-valued modular forms
    - ``W`` --   (system) dictionary
        - ``W['Ms']``  -- M start
        - ``W['Mf']``  -- M stop
        - ``W['nc']``  -- number of cusps
        - ``W['V']``   -- matrix of size ((Ms-Mf+1)*nc)**2
        - ``W['RHS']`` -- right hand side (for inhomogeneous system) matrix of size ((Ms-Mf+1)*nc)*(dim)
    - ``N`` -- normalisation (dictionary, output from the set_norm_for_maass function, default None)
        - if N=None we assume that the solution is uniquely determined by the prinicpal part (in the right hand side)
        - ``N['SetCs']``   -- Which coefficients are set
        - ``N['Vals'] ``   -- To which values are these coefficients set
        - ``N['comp_dim']``-- How large is the assumed dimension of the solution space
        - ``N['num_set']`` -- Number of coefficients which are set
        

        
    - ''gr''  -- only return the reduced matrix and right hand side. do not perform the solving .
    - ''cn''  -- logical (default False) set to True to compute the max norm of V^-1
    OUTPUT:
    
    - ``C`` -- Fourier coefficients

    EXAMPLES::

        sage: G=MySubgroup(Gamma0(1))
        sage: mpmath.mp.dps=20
        sage: R=mpmath.mpf(9.533695261353557554344235235928770323821256395107251982375790464135348991298347781769255509975435366)
        sage: Y=mpmath.mpf(0.5)
        sage: W=setup_matrix_for_Maass_waveforms(G,R,Y,12,22)
        sage: N=set_norm_maass(1)
        sage: C=solve_system_for_Maass_waveforms(W,N)
        sage: C[0][2]*C[0][3]-C[0][6]
        mpc(real='-1.8055426724989656270259e-14', imag='1.6658248366482944572967e-19')

    If M is too large and the precision is not high enough the matrix might be numerically singular

        W=setup_matrix_for_Maass_waveforms(G,R,Y,20,40)  
        sage: C=solve_system_for_Maass_waveforms(W,N)
        Traceback (most recent call last)
        ...
        ZeroDivisionError: Need higher precision! Use > 23 digits!

    Increasing the precision helps
    
        sage: mpmath.mp.dps=25
        sage: R=mpmath.mpf(9.533695261353557554344235235928770323821256395107251982375790464135348991298347781769255509975435366)
        sage: C=solve_system_for_Maass_waveforms(W,N)
        sage: C[0][2]*C[0][3]-C[0][6]
        mpc(real='3.780824715556976438911480324e-25', imag='2.114746048869188750991752872e-99')

    """
    V=W['V']
    Ms=W['Ms']
    Mf=W['Mf']
    #WR=W['WR']
    X = H.multiplier()
    WM = H.multiplier().weil_module()
    nc=W['nc']
    Ml=Mf-Ms+1
    get_reduced_matrix=gr
    verbose = H._verbose
    comp_norm=cn
    r = W['r'] #WM.rank()
    if verbose>0:
        print("Ml={0}".format(Ml))
        print("nv={0}".format(nc))
    #if W['sym_type']==1:
    #    setD=range(0,WM.rank()/2+1)  # 0,1,...,N
    #    #setD=range(0,WR.N+1)  # 0,1,...,N
    #elif W['sym_type']==-1:
    #    setD=range(1,WM.rank()/2)    # 1,2,...,N-1  (since -0=0 and -N=N)
    setD = X.D()
    Ds =setD[0]
    num_comp=nc*r
    if num_comp != len(setD):
        raise Exception(" Wrong dimension of input matrix! r={0},nc={1},setD={2}".format(r,nc,setD))
    if V.ncols()!=Ml*num_comp or V.nrows()!=Ml*num_comp:
        raise Exception(" Wrong dimension of input matrix!")
    if N == None:
        SetCs=[[]]
        Vals=[{}] 
        for b in WM.basis():
            if WM.Q(b) == 0:
                SetCs[0].append((b,0))
                Vals[0][(b,0)]=0
        comp_dim=1
    else:
        SetCs=N['SetCs']
        Vals=N['Vals']
        comp_dim=N['comp_dim']
   

    #for a in Vals.keys():
    #    Vals[a]=mpmath.mp.mpc(Vals[a])
    if Ms<0:
        use_sym=0
    else:
        use_sym=1
    if verbose>0:
        print("Before SetCs={0}".format(SetCs))
    if use_sym == 1:
        if len(SetCs[0]) > 0:
            num_set=0 #
            for (r,m) in SetCs[0]:
                if r in X.D(): #)_as_int:
                    num_set = num_set + 1
    else:
        num_set=len(SetCs[0])
    t=V[0,0]
    CF=MPComplexField(H._prec)
    MS = MatrixSpace(CF,int(Ml*num_comp-num_set),int(comp_dim))
    RHS=Matrix_complex_dense(MS,0,True,True)
    #if W.has_key('RHS'):
    #    if(W['RHS'].ncols()<>comp_dim):
    #        raise ValueError,"Inum_comporrect number of right hand sides!"
    MS = MatrixSpace(CF,int(Ml*num_comp-num_set),int(Ml*num_comp-num_set))
    LHS=Matrix_complex_dense(MS,0,True,True)
    if verbose>0:
        print("Ml={0}".format(Ml))
        print("num_set={0}".format(num_set))
        print("num_comp={0}".format(num_comp))
        print("SetCs={0}".format(SetCs))
        print("setD={0}".format(setD))
        print("Vals={0}".format(Vals))
        print("V.rows={0}".format(V.nrows()))
        print("V.cols={0}".format(V.ncols()))
        print("LHS.rows={0}".format(LHS.nrows()))
        print("LHS.cols={0}".format(LHS.ncols()))
        print("RHS.rows={0}".format(RHS.nrows()))
        print("RHS.cols={0}".format(RHS.ncols()))
        print("use_sym={0}".format(use_sym))
        print("N={0}".format(N))

    num_rhs=0
    if 'RHS' in W:
        num_rhs=W['RHS'].ncols()
    else:
        num_rhs = 1
    if num_rhs!=1 and num_rhs!=comp_dim:
        raise ValueError("Need same number of right hand sides (or just one) as the number of set coefficients!")

    if V.nrows() != num_comp*Ml:
        raise ArithmeticError(" Matrix does not have correct size!")
    roffs=0
    for a in range(num_comp):
        for n in range(Ms,Mf+1):
            r=a*Ml+n-Ms
            if SetCs[0].count((a+Ds,n))>0  and a+Ds in setD:
                #print " set row a,n=",a,n
                roffs=roffs+1
                continue
            for fn_j in range(comp_dim):
                #ztmp=CF(ztmp.real(),ztmp.imag())
                if verbose>3:
                    print("fn_j={0}".format(fn_j))
                    print("r,roffs={0},{1},{2}".format(r,roffs,r-roffs))
                if num_rhs==comp_dim:
                    rhs_j =fn_j
                else:
                    rhs_j = 0
                if 'RHS' in W:
                    RHS[r-roffs,fn_j]=-W['RHS'][r,rhs_j]
                else:
                    RHS[r-roffs,fn_j]=0
                for (i,cset) in SetCs[fn_j]:
                    if i in setD:
                        v=Vals[fn_j][(i,cset)]
                        #k = i*Ml+cset-Ms
                        k = (i-Ds)*Ml+cset-Ms
                        if verbose>3:
                            print("fi,cset={0},{1}".format(i,cset))
                            print("r,k={0},{1}".format(r,k))
                        tmp = CF(v)*V[r,k]
                        RHS[r-roffs,fn_j]=RHS[r-roffs,fn_j]- tmp
            coffs=0
            for b in range(num_comp):
                for l in range(Ms,Mf+1):
                    k=b*Ml+l-Ms
                    if SetCs[0].count((b+Ds,l))>0 and b+Ds in setD:
                        coffs=coffs+1
                        continue
                    if verbose>3:
                        print("r,k={0},{1}".format(r,k))
                        print("r-roffs,k-coffs={0},{1}".format(r-roffs,k-coffs))
                    LHS[r-roffs,k-coffs]=V[r,k]
        #print "LHS[",r,k,"]=",LHS[r-roffs,k-coffs]
    #else:
    #	LHS = V
    #	RHS = -W['RHS']
    if get_reduced_matrix:
        return [LHS,RHS]
    maxit=100;i=0
    done=False
    dps0=CF.prec()
    while (not done and i<=maxit):
        try:
            Q,R = LHS.qr_decomposition()
            done=True
        except ZeroDivisionError:
            t=int(ceil(-log_b(smallest_inf_norm(LHS),10)))
            dps=t+5*i; i=i+1
            print("raising number of digits to:{0}".format(dps))
            LHS.set_prec(dps)
            # raise ZeroDivisionError,"Need higher precision! Use > %s digits!" % t
    if i>=maxit:
        raise ZeroDivisionError("Can not raise precision enough to solve system! Should need > {0} digits! and {1} digits was not enough!".format(t,dps))
    if comp_norm:
        max_norm=LHS.norm()
        for j in range(LHS.rows):
            #y=mpmath_ctx.matrix(LHS.rows,int(1)); y[j,0]=1
            y = Vector_complex_dense(vector(F,LHS.rows).parent(),0)
            y[j]=1
            TMP = RHS.solve(b) #pmath_ctx.U_solve(A, b)
            tmpnorm=max(list(map(abs,TMP)))
            if(tmpnorm>max_norm):
                max_norm=tmpnorm
        print("max norm of V^-1={0}".format(max_norm))
    X=dict()
    for fn_j in range(comp_dim):
        X[fn_j] = dict() #mpmath.matrix(int(Ml),int(1))
        #b = mpmath_ctx.L_solve(A, RHS.column(fn_j), p)
        v = RHS.column(fn_j)
        print("len(B)={0}".format(len(v)))
        TMP = LHS.solve(v)
        #TMP = mpmath_ctx.U_solve(A, b)
        roffs=0
        #res = mpmath_ctx.norm(mpmath_ctx.residual(LHS, TMP, RHS.column(fn_j)))
        res = (LHS*TMP-v).norm()
        print("res(",fn_j,")={0}".format(res))
        for ai in range(num_comp):
            X[fn_j][ai]=dict()
            for n in range(Ms,Mf+1):
                ni=Ml*ai+n-Ms
                if SetCs[fn_j].count((ai+Ds,n))>0 and ai+Ds in setD:
                    roffs=roffs+1
                    X[fn_j][ai][n]=CF(Vals[fn_j][(ai+Ds,n)])
                    continue
                #if verbose>0:
                #    print "ni-roffs=",ni-roffs
                X[fn_j][ai][n]=TMP[ni-roffs]
    return X




def vv_harmonic_wmwf_phase2_tst1(M,PP,C,Ns,Is=None,prec=20,Yin=None):
    try:
        CC=vv_harmonic_wmwf_phase2_1(M,PP,C,Ns,Is,prec,Yin)
        return CC
    except KeyboardInterrupt:
        print("Stopping!")

def vv_harmonic_wmwf_phase2_tst2(M,PP,C,Ns,Is=None,prec=20,Yin=None,do_save=False):
    try:
        CC=vv_harmonic_wmwf_phase2_2(M,PP,C,Ns,Is,prec,Yin,do_save)
        return CC
    except KeyboardInterrupt:
        print("Stopping!")



def vv_harmonic_wmwf_phase2_1(M,PP,C,Ns,Is=None,prec=20,Yin=None):
    r"""
    Phase 2 for vector-valued harmonic weak Maass forms.
    """
    WR=M.WR;
    kappa=M.weight
    D=WR.D  
    Dsym=M.D # the symmetrized index set
    if len(Dsym) != len(C.keys()):
        raise ValueError("Got incompatible coefficient vector! indices={0}".format(list(C.keys())))
    
    #we only use symmetrized values
    if(Is==None):
        Is=[0,len(D)]
        
    N=WR.N
    t=-1
    sym_type=mpmath.mpf(M._sym_type)
    verbose=M._verbose
    #ndig=12
    eps=mpmath.power(mpmath.mpf(10),mpmath.mpf(-prec))
    if verbose > 0:
        print("eps={0}".format(eps))
        print("Yin={0}".format(Yin))
    betai=dict();mbeta=dict(); mbetai=dict(); mm=dict(); mmm=dict()
    mptwo=mpmath.mp.mpf(2); mpfour=mpmath.mp.mpf(4)
    twopi=mptwo*mpmath.pi(); twopii=mpmath.mp.mpc(0,twopi)
    fourpi=mptwo*twopi; mp0=mpmath.mpf(0)
    weight=mpmath.mp.mpf(kappa); weight_over_two=weight/mpmath.mp.mpf(2)
    K0=0; K1=0
    for (beta,m) in PP:
        #if( (not beta in D) or (not 1-beta in D)):
        if( (not beta in D) and (not 1-beta in D)):
            raise Exception("Need beta={0} in D={1}".format(beta,D))
        betai[beta]=D.index(beta)
        mbeta[beta]=1-beta
        mbetai[beta]=D.index(1-beta)
        mm[(beta,m)]=(m+WR.Qv[betai[beta]])
        mmm[(beta,m)]=mpmath.mp.mpf(mm[(beta,m)])
        if verbose > 0:
            print("beta,m={0}, {1}".format(beta,m))
            print("mm={0}".format(mm[(beta,m)]))
            #print "-beta=",minus_beta
            #print "mm=",mm[(beta,m)]
        if mm[(beta,m)]>t:
            t=mm
        if abs(mm[(beta,m)])<K0:
            K0=abs(mm[(beta,m)])
        if abs(PP[(beta,m)])>K1:
            K1=abs(PP[(beta,m)])
        # One way to specify the principal part
        # is to only specify half and let the rest be decided
        # by the symmetry. If we specify the rest it must match
        if (mbeta[beta],m) in PP and (beta,m) in PP:
            test=abs(PP[(beta,m)]-sym_type*PP[(mbeta[beta],m)])
            if test>0: # and not test.ae(mp0):
                raise ValueError("The principal part has not correct symmetry: type={0}, PP={1}".format(sym_type,PP))
        else:
            pass
    abD=len(WR.D)
   
    if(Yin==None):
        Y0=mpmath.mp.mpf(0.5)
    else:
        Y0=mpmath.mp.mpf(Yin)
    kint=mpmath.mp.mpf(1-weight)
    sym_type=M.get_sym_type()
    NA=Ns[0]; NB=Ns[1]
    if(sym_type==1):
        Dstart=int(0); Dfinish=int(WR.N) # 0,1,...,N
    elif(sym_type==-1):
        Dstart=int(1); Dfinish=int(WR.N-1) # 1,2,...,N-1  (since -0=0 and -N=N)
    IA=int(max(Is[0],Dstart)); IB=int(min(Is[1],Dfinish))
    Ms=int(min(C[Dstart].keys())); Mf=int(max(C[Dstart].keys())); Ml=int(Mf-Ms+1)
    #Ms=int(min(C[D[Dstart]].keys())); Mf=int(max(C[D[Dstart]].keys())); Ml=int(Mf-Ms+1)
    if verbose > 0:
        print("Ms,Mf,Ml={0}, {1}, {2}".format(Ms,Mf,Ml))
    K1=K1*2*N
    NAA=NA; IAA=IA
    numys=2
    # have    Y=mpmath.mp.mpf(Y_in)
    # have to find suitable Q for the given Y
    if verbose>0:
        print("dps={0}".format(mpmath.mp.dps))
    ctmp=dict(); ctmp_neg=dict()
    Cout=dict()
    for bi in range(IA,IB+1):
        Cout[bi]=dict()
    stw=str(weight)[0:5]
    Qadd=0; Yfak=mpmath.mpf(0.95); Yvold=list(range(2))
    Xm=dict();Xpb=dict();Ypb=dict(); Cv=dict()
    Q=dict(); Qs=dict();Qf=dict(); QQ=dict()
    for yloop in range(1000):
        Yv=[Y0,Yfak*Y0]
        for i in range(2):
            Q[i]=M.get_M(Yv[i],K0,K1,prec)+Qadd
            Qs[i]=1-Q[i]; Qf[i]=Q[i]; QQ[i]=mpmath.mp.mpf(1)/mpmath.mp.mpf(2*Q[i])
            if verbose>0:
                print("Yv[{0}]={1},{2}".format(i,mppr(Yv[0]),mppr(Yv[1])))
                print("Q(Y)[{0}]={1}".format(i,Q[i]))
                #print "1/2Q[",i,"]=",mppr(QQ[i])
        # Recall that the first Y-value is always the larger
        if verbose > 1:
            print("Yvold={0},{1}".format(mppr(Yvold[0]),mppr(Yvold[1])))
        if Yv[0].ae(Yvold[1]):
            if verbose > 1:
                print("do not evaluate for Yv={0}".format(Yv[0]))

            [Xm[0],Xpb[0],Ypb[0],Cv[0]]=[Xm[1],Xpb[1],Ypb[1],Cv[1]]
            [Xm[1],Xpb[1],Ypb[1],Cv[1]]=pullback_pts_weil_rep(WR,Q[1],Yv[1],weight,Dstart,Dfinish)
        else:
            for i in range(2):
                [Xm[i],Xpb[i],Ypb[i],Cv[i]]=pullback_pts_weil_rep(WR,Q[i],Yv[i],weight,Dstart,Dfinish)
        Yvold=Yv; Zipb=dict()
        for yj in range(0,numys):
            Zipb[yj]=dict()
            for j in range(Qs[yj],Qf[yj]+1):
                Zipb[yj][j]=mpmath.mp.mpc(-Ypb[yj][j],Xpb[yj][j])

        gamma_fak=dict()
        for yj in range(0,numys):
            for bi in range(IA,IB+1):
                for l in range(Ms,Mf+1):
                    lr=mpmath.mp.mpf(l+WR.Qv[bi])
                    if(lr<0):
                        lrtwo=lr*mptwo
                        for j in range(Qs[yj],Qf[yj]+1):
                            gamma_fak[yj,bi,l,j]=mpmath.gammainc(kint,abs(lrtwo)*Ypb[yj][j])*mpmath.mp.exp(-lr*Ypb[yj][j])
                    else:
                        for j in range(Qs[yj],Qf[yj]+1):
                            gamma_fak[yj,bi,l,j]=mpmath.mp.exp(-lr*Ypb[yj][j])
        if verbose > 0:
            print("Got pullback points!")
            print("dps={0}".format(mpmath.mp.dps))
            print("NAA={0}".format(NAA))
        # If we want to do negative coefficients too we save time by computing simultaneously
        do_neg=True
        try:
            for n in range(NAA,NB+1):
                if verbose > 0:
                    print("n=",n)
                for ai in range(IAA,IB+1):
                    if verbose > 0:
                        print("ai={0}".format(ai))
                    mai=-ai % abD
                    nr=mpmath.mp.mpf(n+WR.Qv[ai])
                    nrtwo=mp2*nr
                    nri=mpmath.mp.mpc(0,nr)
                    if(do_neg):
                        nrm=mpmath.mp.mpf(-n+WR.Qv[ai])
                        nrmi=mpmath.mp.mpc(0,nrm)
                        nrmtwo=mp2*nrm
                    for yj in range(0,numys):
                        Y=Yv[yj]*twopi
                        summa=mp0;
                        summa_neg=mp0
                        #print "IA,IB=",IA,IB
                        fak=dict()
                        for j in range(Qs[yj],Qf[yj]+1):
                            fak[j]=mpmath.mp.exp(-nri*Xm[yj][j])
                        if(do_neg):
                            fak_neg=dict()
                            for j in range(Qs[yj],Qf[yj]+1):
                                fak_neg[j]=mpmath.mp.exp(-nrmi*Xm[yj][j])
                        for bi in range(IA,IB+1):
                            mbi=-bi % abD
                            #print "mbi=",mbi
                            for l in range(Ms,Mf+1):
                                if C[bi][l]==0 or abs(C[bi][l]).ae(mp0):
                                    if verbose > 0:
                                        print("Skip coeff {0}, {1}".format(bi,l))
                                    #continue
                                lr=mpmath.mp.mpf(l+WR.Qv[bi])
                                ilr=mpmath.mp.mpc(0,lr)
                                Vtmp=mp0;
                                Vtmp_neg=mp0
                                for j in range(Qs[yj],Qf[yj]+1):
                                    if mbi != bi:
                                        ch=Cv[yj][j][ai,bi]+sym_type*Cv[yj][j][ai,mbi]
                                    else:
                                        ch=Cv[yj][j][ai,bi]
                                    if ch==0 or ch.ae(mp0):
                                        continue
                                    tmp=(ch*mpmath.exp(ilr*Xpb[yj][j]))*gamma_fak[yj,bi,l,j]
                                    Vtmp=Vtmp+tmp*fak[j]
                                    if(do_neg):
                                        Vtmp_neg=Vtmp_neg+tmp*fak_neg[j]
                                summa=summa+Vtmp*C[bi][l]
                                if(do_neg):
                                    summa_neg=summa_neg+Vtmp_neg*C[bi][l]
                        if verbose > 1:
                            print("summa({0})={1}".format(yj,summa))
                            if(do_neg):
                                print("summa_neg({0})={1}".format(yj,summa_neg))
                        wsumma=mp0;
                        wsumma_neg=mp0
                        for (beta,m) in PP:
                            app=mpmath.mp.mpf(PP[(beta,m)])
                            lr=mpmath.mp.mpf(m+WR.Qv[betai[beta]])
                            tmpsumma=mp0;
                            tmpsumma_neg=mp0
                            for j in range(Qs[yj],Qf[yj]+1):
                                if betai[beta] != mbetai[beta]:
                                    ch=Cv[yj][j][ai,betai[beta]]+sym_type*Cv[yj][j][ai,mbetai[beta]]
                                else:
                                    ch=Cv[yj][j][ai,betai[beta]]
                                if ch==0 or ch.ae(mp0):
                                    continue
                                tmp=ch*mpmath.exp(lr*Zipb[yj][j])
                                tmpsumma=tmpsumma+tmp*fak[j]
                                if do_neg:
                                    tmpsumma_neg=tmpsumma_neg+tmp*fak_neg[j]
                            wsumma=wsumma+app*tmpsumma
                            if do_neg:
                                wsumma_neg=wsumma_neg+app*tmpsumma_neg
                        if verbose > 1:
                            print("wsumma({0})={1}".format(yj,wsumma))
                            if(do_neg):
                                print("wsumma_neg({0})={1}".format(yj,wsumma_neg))
                        sumtmp=(summa+wsumma)*QQ[yj]
                        if ((D[ai],n) in PP)>0:
                            sum_tmp=sum_tmp-PP[(D[ai],n)]*mpmath.mp.exp(-nr*Y)
                        lhs=mpmath.mp.exp(nr*Y)
                        if verbose > 0:
                            print("exp(2pinY)={0}".format(mppr(lhs)))
                        ctmp[yj]=sumtmp*lhs
                        if(do_neg):
                            sumtmp_neg=(summa_neg+wsumma_neg)*QQ[yj]
                            if (D[ai],-n) in PP:
                                sumtmp_neg=sumtmp_neg-PP[(D[ai],-n)]*mpmath.mp.exp(-nrm*Y)
                            lhs=mpmath.gammainc(kint,abs(nrmtwo)*Y)*mpmath.mp.exp(-nrm*Y)
                            if verbose > 0:
                                print("Gamma(1-k,4pinY)={0}".format(mppr(lhs)))
                            ctmp_neg[yj]=sumtmp_neg/lhs
                    # end for yj
                    if(verbose>-1):
                        print("C1[{0},{1}]={2}".format(n,ai,ctmp[0].real))
                        print("C2[{0},{1}]={2}".format(n,ai,ctmp[1].real))
                        if(do_neg):
                            print("C1[{0},{1}]={2}".format(-n,ai,ctmp_neg[0].real))
                            print("C2[{0},{1}]={2}".format(-n,ai,ctmp_neg[1].real))
                    if(do_neg):
                        err_pos=abs(ctmp[1]-ctmp[0])
                        err_neg=abs(ctmp_neg[1]-ctmp_neg[0])
                        err=err_pos # max(err_pos,err_neg)
                        if(verbose>-1):
                            print("err_pos={0}".format(mppr(err_pos)))
                            print("err_neg={0}".format(mppr(err_neg)))
                    else:
                        err=abs(ctmp[1]-ctmp[0])
                        if(verbose>-1):
                            print("err={0}".format(mppr(err)))
                    if verbose > 0:
                        if list(C.keys()).count(ai)>0:
                            if n in C[ai]:
                                print("Cin({0},{1})={2}".format(ai,n,C[ai][n].real))
                                if -n in C[ai]:
                                    print("Cin({0},{1})={2}".format(ai,-n,C[ai][-n].real))
                    sys.stdout.flush()
                    if(err>eps):

                        # Have to modify
                        if verbose>0:
                            print(" Need to decrease Y!")
                        Y0=Yv[0]*Yfak
                        #Qadd=Qadd+10
                        #Yv[0]=Y0; Yv[1]=Y0*mpmath.mp.mpf(0.95)
                        NAA=n; IAA=ai
                        raise StopIteration()
                    else:
                        Cout[ai][n]=ctmp[1]
                        if(do_neg):
                            Cout[ai][-n]=ctmp_neg[1]
                        if verbose > 0:
                            print("OK! av={0}".format((ctmp[1]+ctmp[0])/mpmath.mpf(2)))
                        # If we are in the range of the used C's we update
                        #if(C.keys().count(ai)>0):
                        #    if(C[ai].keys().count(n)>0):
                        #        C[ai][n]=ctmp[1]
                        #    if(do_neg and C[ai].keys().count(-n)>0):
                        #        C[ai][-n]=ctmp_neg[1]
                        #continue
                # end for ai
            #print "n at end=",n
            # end for n
            return Cout
        except StopIteration:
            if verbose > 0:
                print("Iteration stopped!")
            continue
    raise StopIteration()





class VVHarmonicWeakMaassFormElement(AutomorphicFormElement):
    r"""
    A harmonic weak Maass form.
    """
    def __init__(self,M,principal_part=None,C=None,prec=53,Lv=None):
        r"""
        Initialize a harmonic weak Maass form element.
        INPUT:

        -''M'' -- space of harmonic weak Maass forms

        -''PP''-- Principal part

        -''C''-- Fourier coefficients

        -''prec'' -- integer, precision (if given by construction, default None)

        EXAMPLES::

        
            sage: WR=WeilRepDiscriminantForm(11,dual=True)
            sage: M=VVHarmonicWeakMaassForms(WR,0.5,100)
            sage: PP={(7/22,0):1}
            sage: F=M.get_element(PP,12);F
            Element of Space of Vector-Valued harmonic weak Maass forms on Modular Group SL(2,Z) of weight 1/2  and dimension 10.
            Representation is Dual of Weil representation of the discriminant form given by ZZ/22ZZ with quadratic form Q(x)=11*x**2 mod 1. with principal part: q^-5/44
        


        """

        if(isinstance(M,type(Newforms(1,12)[0]))):
            self._init_from_newform_(M,principal_part,C,prec)
            return
        if C != None:
            if M.dim != len(C.keys()):
                raise ValueError("Coefficient vector of wrong format! Got length={0}".format(len(C)))
    	## We inherit symmetry from the space
	    self._sym_type = M._sym_type
        self._class_name = "VVHarmonicWeakMaassFormElement"
        AutomorphicFormElement.__init__(self,M,C=C,prec=prec,principal_part=principal_part)

        self._verbose = self._space._verbose
        self.maxdigs=prec # the number of digits needed to be displayed to print all digits of the coefficients
        if Lv != None:
            self._Lv=Lv
        else:
            self._Lv=dict()
        self.Cp0=0; self.Cp1=0; self.Cm=0        
        ## We also find the space corresponding to this space via inverse of xi_k and Shimura corr.
        t=QQ(RR(M.weight))-QQ(1)/QQ(2)
        if(is_int(t)):
            k=QQ(3)-QQ(2)*QQ(RR(M.weight))
            if self._verbose > 0:
                print("k={0}".format(k))
            self.shim_corr=Newforms(self._space.WR.N,k,names='a')  #.new_subspace()
        else:
            self.shim_corr=None


   

    def __reduce__(self):
        r"""
        """
        #return(HarmonicWeakMaassFormElement,(self.space,self.principal_part.items(),self._coeffs,self.prec))
        return(VVHarmonicWeakMaassFormElement,(self._space,self._principal_part,self._coeffs,self.prec,self._Lv))

    def _repr_(self):
        r""" Return string representation of self.

        EXAMPLES:

            sage: WR=WeilRepDiscriminantForm(11,dual=True)
            sage: M=VVHarmonicWeakMaassForms(WR,0.5,100)
            sage: PP={(7/22,0):1}
            sage: F=M.get_element(PP,12);F
            Element of Space of Vector-Valued harmonic weak Maass forms on Modular Group SL(2,Z) of weight 1/2  and dimension 10.
            Representation is Dual of Weil representation of the discriminant form given by ZZ/22ZZ with quadratic form Q(x)=11*x**2 mod 1. with principal part: q^-5/44
            
        
        """
        
        s="Element of "+str(self._space)+" with principal part: "
        WR=self._space.multiplier()
        sp=""
        for (b,m) in self._principal_part:
            a=self._principal_part[(b,m)]
            if(a!=0):
                x=QQ(m+WR.Qv[WR.D.index(b)])
                if(a!=1):
                    if(a>0 and len(sp)>0):
                        ast="+"+str(a)
                    else:
                        ast=str(a)
                    sp=sp+ast+"q^{"+str(x)+"}"
                else:
                    sp=sp+"q^{"+str(x)+"}"
        s=s+sp
        return s

    def _latex_(self):
        r""" Return LaTeX string representation of self.

        EXAMPLES:


            sage: WR=WeilRepDiscriminantForm(11,dual=True)
            sage: M=VVHarmonicWeakMaassForms(WR,0.5,100)
            sage: PP={(7/22,0):1}
            sage: F=M.get_element(PP,12);F
            Element of Space of Vector-Valued harmonic weak Maass forms on Modular Group SL(2,Z) of weight 1/2  and dimension 10.
            Representation is Dual of Weil representation of the discriminant form given by ZZ/22ZZ with quadratic form Q(x)=11*x**2 mod 1. with principal part: q^-5/44
            
        
        """
        old=s="\\begin{verbatim}\\end{verbatim}"
        new=""
        s="\\begin{verbatim}\\end{verbatim}"
        s+="Element of "+self._space._latex_().replace(old,new)+" With principal part "
        WR=self._space.multiplier()
        # If we have more than one non-zero element in the principal part we have to
        # addd a + between terms
        sp=""
        for (b,m) in self._principal_part:
            a=self._principal_part[(b,m)]
            if a!=0:
                x=QQ(m+WR.Qv[WR.D.index(b)])            
                if a!=1:
                    if a>0 and len(sp)>0:
                        ast="+"+str(a)
                    else:
                        ast=str(a)
                    sp=sp+ast+"q^{"+str(x)+"}"
                else:
                    sp=sp+"q^{"+str(x)+"}"
        s=s+sp+"$."
        return s


    def _init_from_newform_(self,g,PP,C,prec):
        r""" Init a function from a newform (without computing Fourier coefficients). 

        """
        k=g.weight()
        N=g.level()
        if(k % 4==0):
            ep=g.atkin_lehner_eigenvalue()
        else:
            ep=-g.atkin_lehner_eigenvalue()
        kappa=mpmath.mpf(3-k)/mp2
        if(ep==-1):
            WR=WeilRepDiscriminantForm(N,False)
        elif ep==1:
            WR=WeilRepDiscriminantForm(N,True)
        else:
            raise ValueError(" Sign of functional equation must be 1 or -1! Got:{0}".format(ep))
        #print "kappa=",kappa
        #print "WR=",WR
        M=VVHarmonicWeakMaassForms(WR,kappa,prec)
        #print "M=",M
        self._space=M
        self.prec=prec
        self.coeff=dict()
        # We want a Harmonic weak maass form corresponding to the form g
        # that means we need to avoid any other cuspforms as well as
        # theta series...
        # If there are no oldforms we are happy
        if dimension_new_cusp_forms(N,k)==dimension_cusp_forms(N,k):
            # and to avoid theta series we need to avoid square discriminants
            # in the principal part
            if M._is_dual_rep:
                nset=[0,-1]
            else:
                nset=[-1]
            try:
                for n in nset:
                    for r in WR.D(): #)_as_integers:
                        D=M.D_from_rn((r,n))
                        if not is_square(D):
                            PP={(WR.D[r],n):1}
                            self._principal_part=PP
                            #print "PP=",PP,"is ok!"
                            raise StopIteration()
            except StopIteration:
                pass
        #if(C<>None and C>0):
            
    def get_principal_part(self,str=False,disc=False):
        r""" Return principal part of self.

        INPUT:

        -disc -- logical (default False) if True, return the principal part as discriminants

        EXAMPLES:

            sage: WR=WeilRepDiscriminantForm(11,dual=True)
            sage: M=VVHarmonicWeakMaassForms(WR,0.5,100)
            sage: PP={(7/22,0):1}
            sage: F=M.get_element(PP,12);F
            Element of Space of Vector-Valued harmonic weak Maass forms on Modular Group SL(2,Z) of weight 1/2  and dimension 10.
            Representation is Dual of Weil representation of the discriminant form given by ZZ/22ZZ with quadratic form Q(x)=11*x**2 mod 1. with principal part: q^-5/44       """
        if(not disc and not str):
            return self._principal_part
        
        L=list()
        for (r,n) in self._principal_part:
            D=self._space.D_from_rn((r,n))
            L.append((D,self._principal_part[(r,n)]))
        if(disc):
            return L
        if(str):
            if(len(L)==1 and L[0][1]==1):
                return L[0][0]
            else:
                return L


    def pairing(self,g,t=1):
        r""" Compute the bilinear pairing {g,self} for a holomorphic form g 
        """
        # Have to make test of g!!!
        #raise NotImplementedError," Need a proper class of holomorphic vector-valued forms!"
        # # We need to obtain the coefficients of the inverse Shimura rep. first
        # # if we want to apply this to a scalar holomorphic form g
        # First locate the maximum discriminant we need
        Dmax=0
        PP=self._principal_part
        for (r,n) in PP:
            D=self._space.D_from_rn((r,n))
            if(abs(D)>abs(Dmax)):
                Dmax=D
        sig=sign(Dmax)
        Dmax=10 # abs(Dmax)
        print("Dmax={0}".format(sig*Dmax))
        #t=1 # if this doesn't work we have to choose another t
        # I also assume that g and G have trivial characters
        syst=matrix(ZZ,Dmax,Dmax)
        rhs=matrix(ZZ,Dmax,1)
        k=Integer(g.weight()/Integer(2))
        for n in range(1,Dmax+1):            
            rhs[n-1,0]=g[n]
            for d in range(1,Dmax+1):
                if(n % d !=0):
                    continue
                ## I use the character d -> (4N / d) 
                #chi=(kronecker(-1,d)**k)*kronecker(t,d)*kronecker(d,F.space.WR.level())
                chi=kronecker(t,d) #*kronecker(d,F.space.WR.level)
                am=chi*d**(k-1)
                #print "am[",n,d,"]=",am
                syst[n-1,n/d-1]=am
        X=syst.solve_right(rhs)
        C=dict()
        for j in range(1,Dmax+1):
            C[t*j**2]=X[j-1,0]
            print("C[{0}={1}".format(t*j**2,X[j-1,0]))
        return C


        summa=0

        PP=self._principal_part
        for (r,n) in PP:
            summa=summa+PP[(r,n)]*g.coeff((r,n))
        return summa
    
    def compute_coefficients(self,nrange,irange=None,prec=10,ef=True,Qadd=0):
        r""" Compute a list of coeficients.
        INPUT:

        - nrange -- range of integers
        - irange --  range of integers
        - prec   -- integer
        - ef     -- logical
        - Qadd   -- integer   (take more sample points along the horocycle)

        OUTPUT:
        Coefficients C(i,n) for i  in irange and n in nrange.

        EXAMPLES:

            sage: WR=WeilRepDiscriminantForm(11,dual=True)
            sage: M=VVHarmonicWeakMaassForms(WR,0.5,100)
        """
        # we first need an initial set of coefficients
        C=self._coeffs; P=self._principal_part; M=self._space
        WR=M.multiplier(); weight=M.weight
        if(self.prec>= prec or len(C)>0):
            # presumable we already have good coefficients
            pass
        else:
            # Need initial set first
            print("Computing initial set of coefficients!")
            self.prec=prec
            [Y,M0]=self._space.get_Y_and_M(P,weight,prec)
            Q=M0+10
            W=vv_harmonic_wmwf_setupV(WR,P,Y,M0,Q,weight,self._space._sym_type,verbose=self._space._verbose)
            if (0,0) in P:
                N = self._space.set_norm()
                # N=set_norm_vv_harmonic_weak_maass_forms(WR,cusp_form=True,holomorphic=self._holomorphic)
            else:
                N = self._space.set_norm()
                # N=set_norm_vv_harmonic_weak_maass_forms(WR,cusp_form=False,holomorphic=self._holomorphic)
            C=solve_system_for_vv_harmonic_weak_Maass_waveforms(W,N,verbose=self._verbose)
            
        # endif
        # check if we  have all coefficients we wanted
        maxc=max(C[list(C.keys())[0]].keys())
        if maxc >= max(nrange):
            print("Have all we need!")
            pass
        else:
            # we do not have all coefficients we need
            print("Need to compute more!!")
            Ns=nrange # [maxc,max(nrange)]
            if irange!=None:
                Is=irange
            else:
                Is=[min(M.D()),max(M.D())]

            # Try to find good Y
            # Recall that the error in the negative part is usually smaller than in the positive part
            M_minus=abs(min(self._coeffs[list(self._coeffs.keys())[0]]))
            M_plus=abs(max(self._coeffs[list(self._coeffs.keys())[0]]))
            # Assume we computed these coefficients at (almost) the highest horocycle
            Y0=mpmath.sqrt(3)/mpmath.mpf(2)*mpmath.mpf(0.995)
            [err_minus,err_plus]=self.get_error_estimates(Y0,M_minus,M_plus)
            kint=mpmath.mp.mpf(1-self._space.weight)
            print("original:")
            print("err-={0}".format(err_minus))
            print("err+={0}".format(err_plus))
            Y0=mpmath.mpf(0.5)
            Yin=Y0
            for j in range(5000):
                Y=Y0*mpmath.power(mpmath.mpf(0.99),j)
                t=mpmath.pi()*2*Y*abs(Ns[0])
                tmp1=mpmath.exp(t)
                err1=err_plus*tmp1
                #print "err+=",err1
                tmp2=mpmath.gammainc(kint,2*t)
                err2=err_plus*mpmath.exp(-t)/tmp2
                #print "err-=",err2
                if(max(err1,err2)<mpmath.power(10,-prec)):
                    Yin=Y
                    break
                #t=max(1.0,abs(mpmath.log10(prec)-mpmath.log10(self.prec)))
                #Yin=t/mpmath.mpf(Ns[0]+Ns[1])*mpmath.mpf(2.0) ## This should be good on average
            #Yin=Yin*mpmath.mpf(0.2)
            print("err={0}".format(max(err1,err2)))
            print("Yin={0}".format(Yin))
            sys.stdout.flush()
            #Qadd=40
            try:
                if(ef):
                    CC=vv_harmonic_wmwf_phase2_2_ef(self,Ns,Is,prec,Yin,Qadd_in=Qadd)
                else:
                    CC=vv_harmonic_wmwf_phase2_2(M,P,C,Ns,Is,prec,Yin)
                for x in CC.keys():
                    C[x]=CC[x]
            except KeyboardInterrupt:
                print("Manually stopping...")


    def get_error_estimates(self,Y,M1,M2=None):
        r""" Compute the constants needed to make error estimates.
        """
        # First K0 and K1
        Mminus=M1
        if M2==None:
            Mplus=M1
        else:
            Mplus=M2
        if self.Cp0 != 0 and self.Cp1 != 0 and self.Cm != 0:
            Cp0=self.Cp0; Cp1=self.Cp1; Cm=self.Cm
        else:
            PP=self.principal_part()
            Cmax=max(PP.values());Kmax=0
            for t in PP.keys():
                if isinstance(t,tuple):
                    (c,l) = t
                elif isinstance(t,(int,Integer)):
                    (c,l)=rn_from_D(self._space.multiplier(),t)
                else:
                    raise ValueError("Incorrect principal part: t={0}".format(t))
                if c in self._space.multiplier().D():
                    tmp=l+self._space.multiplier().Qv[self._space.index_set().index(c)]
                elif c in range(len(self._space.multiplier().Qv)):
                    tmp=l+self._space.multiplier().Qv[c]
                else:
                    raise ValueError("Incorrect principal part: c,l={0},{1}".format(c,l))
                if(abs(tmp)>Kmax):
                    Kmax=abs(tmp)
            [Cp0,Cp1]=self._space.get_Cp(Cmax)
            Cm=self._space.get_Cm(Kmax,Cmax)
            self.Cp0=Cp0; self.Cp1=Cp1; self.Cm=Cm

        fak=len(self._space.index_set())
        #print "Cp0,Cp1,Cm=",Cp0,Cp1,Cm
        #print "fak=",fak

        er1=fak*self._space.err_est_vv_hwmf_neg(Y,Mminus,Cm)
        er2=fak*self._space.err_est_vv_hwmf_pos(Y,Mplus,Cp0,Cp1)
        return [er1,er2]
        
    def get_coefficient(self,L,n=None):
        r""" Return a coefficient or a list of coeficients.


        EXAMPLES:

            sage: WR=WeilRepDiscriminantForm(11,dual=True)
            sage: M=VVHarmonicWeakMaassForms(WR,0.5,100)
            sage: F=M.get_element({(7/22,0):1},12)
        """
        if(isinstance(L,list)):
            l=list()
            for t in L:
                if isinstance(t,tuple):
                    tt=t
                elif is_int(t):
                    tt=rn_from_D(self._space.multiplier(),t)
                if tt != None:
                    c=self.get_one_coefficient(tt[0],tt[1])
                    l.append(c)
            return l
        elif is_int(n):
            return self.get_one_coefficient(L,n)
        elif is_int(L):
            tt=rn_from_D(self._space.multiplier(),L)
            if tt != None:
                return self.get_one_coefficient(tt[0],tt[1])
        else:
            raise ValueError("Incorrect keys for coefficents: L,n={0}, {1}".format(L,n))

    def C(self,r,n=None):
        r"""
        Alias to get_coefficient.
        """
        return self.get_coefficient(r,n)
            
    def get_one_coefficient(self,r,n):
        r""" Return coefficient c(r,n) if it exists
        """
        c=None
        if not r in self._space.multiplier().D(): # and not r in self._space.WR.D_as_integers:
            raise ValueError("Component r={0} is not valid!".format(r))
        # We also see if the coefficient can be provided via symmetry
        # If r is in D we swap it to its index
        if r in self._space.multiplier().D():
            rr=self._space.multiplier().D().index(r)
        else:
            rr=r
        minus_rr=(len(self._space.multiplier().D())-rr) % len(self._space.multiplier().D())
        #print "rr=",rr
        #print "-rr=",minus_rr
        if self._space._is_dual_rep:
            if rr == minus_rr:
                return 0
        if rr in self._coeffs:
            if n in self._coeffs[rr]:
                c=self._coeffs[rr][n]        
        elif minus_rr in self._coeffs:
            if n in self._coeffs[minus_rr]:
                c=self._coeffs[minus_rr][n]
        return c

    def add_coefficients_from_file(self,file=None,overwrite=False,nmax=None):
        r"""
        Add a set of coefficients from a file of the format:
        r n C(r,n)
        """
        C=dict()
        f=open(file,'r')
        i=0
        md=mpmath.mp.dps
        mpold=mpmath.mp.dps
        for line in f:
            i=i+1
            if nmax != None and i > nmax:
                break
            #return line
            l=line.strip().split()
            if len(l) < 2:
                continue
            #print "l=",l
            r=int(l[0])
            n=int(l[1])
            #print "r=",r,"n=",n,self._coeffs.keys(),self._coeffs.keys().count(r)
            if(list(self._coeffs.keys()).count(r)==0):
                continue
            cs="".join(l[2:len(l)])
            #print cs
            ## see if the string is given in arprec format: 10^a x b
            if(cs.count("^")==1 and cs.count("x")==1):
                s=cs.split("^")[1]
                a=s.split("x")
                cs=a[1]+"E"+a[0]                
            cs=cs.replace(" ","")
            if(len(cs)>md):
                md=len(cs)
            ##mpmath.mp.dps=md
            #if(r==1 and n==79):
            #    print "l=",l
            #    print "c(",r,n,")=",cs
            c=mpmath.mpf(cs)
            #mpmath.mp.dps=mpold
            #print "c(",r,n,")=",c
            C[(r,n)]=c
        #return C
        self.add_coefficients(C,overwrite)
            
    def add_coefficients(self,L,overwrite=False):
        r""" Add one or more coefficients to self.

        INPUT:

        -''L'' -- dictionary of pairs of indices and coefficients
        -''overwrite'' -- logical, set to True if we want to overwrite present coefficients
        
        """
        if not isinstance(L,dict):
            raise ValueError("Call with dictionary as argument!")

        for p in L.keys():
            c=mpmath.mpmathify(L[p])
            #print "c=",c
            cd=ceil(mpmath.log10(abs(c)))
            if(cd>self.maxdigs):
                self.maxdigs=cd
            #print "p=",p
            if(is_int(p)):
                (r,n)=rn_from_D(self._space.WR,p)
            elif(isinstance(p,tuple)):
                (r,n)=p
            if r in self._coeffs:
                if n in self._coeffs[r]:
                    c_old=self._coeffs[r][n]
                    ## Try to determine (heuristically) if the new coefficient is really better
                    d1=dist_from_int(c)[0]
                    d2=dist_from_int(c_old)[0]
                    if(overwrite):
                        self._coeffs[r][n]=c
                else:
                    self._coeffs[r][n]=c
            else:
                # see if it is a possible index at all
                if not r < 0 or r > self._space.multiplier().ambient_rank():
                    raise ValueError("Key {0} corr to (r,n)=({1},{2}) is invalid for the current space!".format(p,r,n))
                elif r not in self._space.multiplier().D():
                    if self._space._sym_type==-1 and (r==0 or r==self._space.multiplier().N):
                        # Should be 0 by symmetry
                        if abs(c) > 10**(1-self.prec):
                            raise ValueError("Coefficient should be zero by symmetry. Got c({0},{1})={2}!".format(r,n,c))
                        else:
                            self._coeffs[r][n]=0
                    else:
                        # is equal to +- c(-r,n)
                        mr=2*self.multiplier().N-r
                        if mr in self._coeffs:
                            if n in self._coeffs[mr]:
                                c_old=self._coeffs[mr][n]
                                if abs(c-self._space.multiplier()._sym_type*c_old) > 10**(1-self.prec):
                                    st="Might add an erronous coefficient! Got c({0},{1})={2}. ".format(r,n,c)
                                    st+="From previous coefficients should have {0}".format(self._space._sym_type*c_old)
                                    raise ValueError(st)
                                if overwrite:
                                    self._coeffs[mr][n]=c
                else:
                    raise ValueError("Coefficient should be zero by symmetry. Got c({0},{1})={2}!" .format(r,n,c))


    def list_coefficients(self,format='components',fd=True,pos=True,neg=True,printimag=False,norm_neg=True,nmin=0,nmax=0,latex=False,nd=0,Lv=False,prime=False):
        r""" List all coefficients C^{+}(Delta} corresponding to fundamental discriminants up to Dmax.

        INPUT:

        
            -``format`` -- string.
                        == components (default) means that we list coefficients component-wise
                        == disc menas that we list according to discriminant
            -``fd`` -- logical (default True) if set to True only prints coefficients given by fundamental discriminants                    
            -``max`` -- integer (default 0) the largest coefficient (either max n or max D). If 0 we list all coefficients we got

            -``pos``  -- logical (default True) if set to true prints positive coefficients

            -``neg``  -- logical (default True) if set to true prints negative coefficients

            -``norm_neg`` -- logical( default True) if we want to normalize the negative coefficients by dividing c(D) with some non-zero c(j) and sqrt(D)

            -``printimag``  -- logical (default False) if set to True prints imaginary parts (which otherwise are assumed zero)

            -``max``  -- integer (default 0) if >0 denotes the maximum index of coefficients to print

            -''Lvals'' -- logical (default False) if True add values from self._Lv in a third column

            -''prime'' -- only list coefficients of discriminant relatively prime to the level

        """        
        M=self._space; WR=M.WR
        C=self._coeffs
        if format[0]=="C" or format[0]=="c":
            self._list_coefficients_by_components(fd,pos,neg,printimag,norm_neg,nmin,nmax,latex,nd,Lv,prime)
        else:
            self._list_coefficients_by_discriminant(fd,pos,neg,printimag,norm_neg,nmin,nmax,latex,nd,Lv,prime)


    def _list_coefficients_by_components(self,fd=True,pos=True,neg=True,printimag=False,norm_neg=True,nmin=0,nmax=0,latex=False,nd=0,Lvals=False,prime=False):
        r""" List all coefficients C^{+}(Delta} corresponding to fundamental discriminants up to Dmax.

        INPUT:

        -``fd`` -- logical (default True) if set to True only prints coefficients given by fundamental discriminants
        -``pos``  -- logical (default True) if set to True prints positive coefficients
        -``neg``  -- logical (default True) if set to True prints negative coefficients
        -``printimag``  -- logical (default False) if set to True prints imaginary parts (which otherwise are assumed zero)
        -``norm_neg`` -- logical( default True) if we want to normalize the negative coefficients by dividing c(D) with some non-zero c(j) and sqrt(D)
        -``nmax``  -- integer (default 0) if >0 denotes the maximum index of coefficients to print

        -''Lvals'' logical (default False) if True add values from self._Lv in a third column
        -''prime'' -- only list coefficients of discriminant relatively prime to the level
        
        """
        sig=1
        if(self._space.WR.is_dual()):
            sig=-1
        maxi=max(self._coeffs.keys())
        w1=len(str(maxi))
        w2=max(list(map(len,str(self._space.WR.D()).split())))
        maxn=max(self._coeffs[list(self._coeffs.keys())[0]].keys())
        w3=len(str(maxn))+1
        C=self._coeffs
        mp0=mpmath.mpf(0)
        mpold=mpmath.mp.dps
        N=self._space.WR.N
        if(mpmath.mp.dps < self.maxdigs):
            mpmath.mp.dps=self.maxdigs
        if norm_neg:
            cnorm=0
            tnorm=(0,0)
            for j in range(1,100):
                t=rn_from_D(self.space.WR,-j*sig)
                if(t==None):
                    continue
                if(t[1]+self._space.WR.Qv[t[0]]>=0):
                    continue
                c1=self.get_coefficient(t[0],t[1])
                if c1==None:
                    continue
                if abs(c1)>self._prec:
                    cnorm=c1
                    tnorm=t
                    print("c1=c({0})={1}".format(tnorm,cnorm))
                    break
        for r in C.keys():
            for n in range(min(C[r].keys()),max(C[r].keys())+1):
                if nmin > 0 and abs(n) < nmin:
                    continue
                if nmax > 0 and abs(n) > nmax:
                    continue
                nn=n+self._space.WR.Qv[r]
                if not neg and nn < 0:
                    continue
                if not pos and nn >= 0:
                    continue
                D=self._space.D_from_rn((r,n))
                if(fd):
                    if fd and not is_fundamental_discriminant(D) and D != 1:
                        continue
                if prime and gcd(D,N)>1:
                    continue
                c=self.get_coefficient(r,n)
                cs=""
                if c != 0 and c != None:
                    if(nn>=0): ss="+"
                    if(nn<0):
                        ss="-"
                        if(norm_neg):
                            #print "r,n=",r,n
                            #print "cnorm=",cnorm
                            #print "tnorm=",tnorm
                            D=self._space.D_from_rn((r,n))
                            if ((r,n)!= tnorm) and cnorm != 0:
                                c=c/cnorm*mpmath.sqrt(mpmath.mpf(abs(D)))
                    if c.real() >= 0: cs=" "
                    if printimag == False:
                        if nd > 0:
                            cs=str(c.real()).strip()
                            cs=sci_pretty_print(cs,nd,latex_pow=latex)
                        else:
                            cs=str(c.real())
                    else:
                        cs=cs+str(c)
                    if Lvals and list(self._Lv.keys()).count(D) == 1:
                        ls="\t"+str(self._Lv[D])
                    else:
                        if latex:
                            ls="\\\\ \n"
                        else:
                            ls=""
                    if latex:
                        D=self._space.WR.D()[r]
                        if(is_int(D)):
                            p=numerator(D); q=denominator(D)                        
                            sr="\\frac{"+str(p)+"}{"+str(q)+"}"
                        else:
                            sr=str(D)
                        ss=""
                        print("$C{0}({1},{2}) $ & $ {3} $ {4}".format(ss,sr.ljust(w1),str(n).ljust(w3),cs,ls))
                    else:
                        print("C^{0}[{1}][{2}] = {3}".format(ss,str(r).ljust(w1),str(n).ljust(w3),cs+ls))
        mpmath.mp.dps=mpold


                
    def _list_coefficients_by_discriminant(self,fd=True,pos=True,neg=True,printimag=False,norm_neg=True,dmin=0,dmax=0,latex=False,nd=0,Lvals=False,prime=False):
        r""" List all coefficients C^{+}(Delta} corresponding to fundamental discriminants up to Dmax.

        INPUT:

        -``fd`` -- logical (default True) if set to True only prints coefficients given by fundamental discriminants

        -``pos``  -- logical (default True) if set to True prints positive coefficients

        -``neg``  -- logical (default True) if set to True prints negative coefficients

        -``printimag``  -- logical (default False) if set to True prints imaginary parts (which otherwise are assumed zero)

        -``norm_neg`` -- logical( default True) if we want to normalize the negative coefficients by dividing c(D) with some non-zero c(j) and sqrt(D)

        -``dmax``  -- integer (default 0) if >0 denotes the maximum index of coefficients to print

        -''latex'' logical (default False) set to true if you want a latex table

        -''Lvals'' logical (default False) if True add values from self._Lv in a third column

        -''prime'' -- only list coefficients of discriminant relatively prime to the level
         
        """
        sig=1
        S="$"
        if(self._space.WR.is_dual()):
            sig=-1
        maxn=max(self._coeffs[list(self._coeffs.keys())[0]].keys())
        maxD=self._space.WR.level()*(maxn+1)
        N=self._space.WR.N
        if(dmax>0):
            w1=len(str(dmax))+1
        else:
            w1=len(str(maxD))+1
        w2=max(list(map(len,str(self._space.WR.D()).split())))
        w3=len(str(maxn))+1
        mp0=mpmath.mpf(0)
        mpold=mpmath.mp.dps
        if(mpmath.mp.dps < self.maxdigs):
            mpmath.mp.dps=self.maxdigs
        if(norm_neg and neg):
            cnorm=0
            tnorm=(0,0)
            for j in range(1,100):
                t=rn_from_D(self.space.WR,-j*sig)
                if(t==None):
                    continue
                if(t[1]+self._space.WR.Qv[t[0]]>=0):
                    continue
                c1=self.get_coefficient(t[0],t[1])
                if(c1 == None):
                    continue
                #print "c1 =",c1
                # If the first coefficient is zero to the precision we assume we shouldn't divide by it
                if(abs(c1)>mpmath.power(10,-self.prec)):
                    cnorm=c1*mpmath.sqrt(j)
                    tnorm=t
                    print("c1=c({0})=c({1})={2}".format(tnorm,-j*sig,cnorm))
                    break
                
        for sn in [1,-1]:
            for D in range(1,maxD):
                #print "D=",D
                if(dmin>0 and abs(D)<dmin):
                    continue
                if dmax > 0 and abs(D) > dmax:
                    continue
                DD=sig*D*sn
                # print "D=",D,is_fundamental_discriminant(D)
                if fd and not is_fundamental_discriminant(DD) and DD != 1:
                    # print "D=",D,is_fundamental_discriminant(D)
                    continue
                if prime and gcd(D,N) > 1:
                    continue
                t=rn_from_D(self._space.WR,DD)
                if t == None:
                    continue
                else:
                    (r,n)=t
                #print " DD=",DD,t
                nn=n+self._space.WR.Qv[r]
                if(not pos and nn>=0):
                    continue
                if(not neg and nn<0):
                    continue
            
                c=self.get_coefficient(r,n)
                cs=""
                erms="";erm=10
                if c != 0 and c != None:
                    if nn >= 0: ss="+"
                    if nn < 0:
                        ss="-"
                        if(norm_neg):
                            if ((r,n) != tnorm) and cnorm != 0:
                                c=c/cnorm*mpmath.sqrt(mpmath.mpf(abs(D)))
                            x=c.real(); x1=floor(x); x2=ceil(x); er1=abs(x1-x); er2=abs(x2-x)                        
                            erm=min(er1,er2);erms=sci_pretty_print(erm,2,latex_pow=latex)
                    if(erm<0.001):
                        if(er1<er2):
                            cs=str(x1)
                        else:
                            cs=str(x2)
                    elif(printimag==False):
                        if(nd>0):
                            cs=str(c.real()).strip()
                            cs=sci_pretty_print(cs,nd,latex_pow=latex)
                        else:
                            cs=str(c.real())
                    else:
                        if(nd>0):
                            cs=str(c).strip()
                            cs=sci_pretty_print(cs,nd,latex_pow=latex)
                        else:
                            cs=str(c)                        
                    if(c.real()>=0 and latex):
                        cs="\hphantom{-}"+cs
                    elif(c.real()>=0):
                        cs=" "+cs 
                    if(latex):
                        O=" & "
                        if(Lvals and list(self._Lv.keys()).count(DD)==1):
                            ls="&"+S+sci_pretty_print(self._Lv[DD],nd,latex_pow=latex)+S
                        else:
                            ls=""
                        if(len(erms)==0):
                            s= S+str(DD).center(w1)+S+"&"+S+cs+S+ls+"\\\\"
                        else:
                            s= S+str(DD).center(w1)+S+"&"+S+cs+S+ls+O+S+erms+S+"\\\\"
                    else:
                        if(Lvals and list(self._Lv.keys()).count(DD)==1):
                            ls="\t"+sci_pretty_print(self._Lv[DD],nd)
                        else:
                            ls=""
                        if(len(erms)==0):
                            s= "C^"+ss+"["+str(DD).center(w1)+"] = "+cs+ls
                        else:
                            s= "C^"+ss+"["+str(DD).center(w1)+"] = "+cs+ls+" "+erms+"\n"
                        #                s=s+str(self._space.WR.D[r]).ljust(w2)+","+str(n).ljust(w3)+"] = "+cs
                    print(s)
        mpmath.mp.dps=mpold

        
    def add_lderiv_vals(self,file=None):
        r""" Reads values of the corresponding central derivatives: L'(G,chi_D,3/2-k) ""
        """
        Lv=dict()
        f=open(file,'r')
        for line in f:
            #print "line=",line
            #print "split=",line.split(" ")
            l=line.rsplit(" ",1)
            if(len(l)==1):
                l=line.rsplit("\t",1)
            D=l[0].strip(); x=l[1].strip()
            Lv[int(D)]=float(x)
        self._Lv=Lv

    def find_vanishing_lderivs(self,do_print=True,latex=True,nd=50):
        r""" Returns a list with all discriminats for which it is likely that L'(G,chi_D,3/2-k)=0
        """
        res=list()
        if(latex):
            S=" $ "; O=" & "
        else:
            S=" "; O=" "
        if(len(list(self._Lv.keys()))==0):
            return res
        L=list(self._Lv.keys()); L.sort(); L.reverse()
        s=""; sc=""
        ## increase mpmath.mp.dps to print all relevant digits
        mpold=mpmath.mp.dps
        mpmath.mp.dps=self.maxdigs
        for DD in L:
            x=self._Lv[DD]
            if(abs(x)<1E-10):
                #res.append((DD,x))
                res.append(DD)
                s=s+S+str(DD)+S+O+S+sci_pretty_print(self._Lv[DD],nd,latex_pow=latex)+S+"\\\\ \n"
                c=self.get_coefficient(DD)
                if c != None:
                    x=c.real(); x1=floor(x); x2=ceil(x); er1=abs(x1-x); er2=abs(x2-x)
                    erm=min(er1,er2)
                    print("erm({0})={1}".format(DD,erm))
                    erms=sci_pretty_print(erm,2,latex_pow=latex)
                    if(er1<er2):
                        xi=x1;
                    else:
                        xi=x2
                    #sc=sc+S+str(DD)+S+"\t"+O+S+sci_pretty_print(c.real,nd,latex_pow=latex)+"\\ \n"
                    sc=sc+S+str(DD)+S+O+S+str(xi)+S+O+S+erms+S+"\\\\ \n"
                else:
                    sc=sc+S+str(DD)+S+O+S+" "+S+O+S+" "+S+"\\\\ \n"
        print(s)
        print(sc)
        mpmath.mp.dps=mpold
        return res





### Routines to switch between D's and (r,n)'s
@cached_function
def rn_from_D(WR,D):
    r""" Find the pair(s) (r,n) s.t. +-D/4N= n +- q(r) for D in D 

    INPUT:
    -''D'' -- integer or list of integers

    OUTPUT:
    -''t'' -- tuple (r,n) or list of tuples



    """
    if(isinstance(D,list)):
        lout=list()
        for DD in D: 
            t=one_rn_from_D(WR,DD)
            if t != None:
                lout.append(t)
        return lout
    else:
        return one_rn_from_D(WR,D)

@cached_function
def one_rn_from_D(WR,D):
    r""" Find the (r,n) s.t. +-D/4N= n +- q(r)
    """            
    Dv=QQ(D)/QQ(WR.level())
    sig=1
    if WR.is_dual():
        sig=-1
    for r in WR.D():
        x=WR.Qv[r]
        nn = Dv-sig*x
        print("D/N -{0} Q({1})={2}".format(sig,r,x))
        if is_int(nn):
            rr=WR.D().index(r)
            n=int(Dv-sig*x)
            return (rr,n)
    return None

@cached_function
def D_from_rn(WR,t):
    r""" Find the D s.t. +-D/Level = n +- q(r)
    """
    if isinstance(t,list):
        lout=list()
        for (r,n) in t:
            D=_one_D_from_rn(WR,(r,n))
            if D!=None:
                lout.append(D)
        return lout
    else:
        return _one_D_from_rn(WR,t)


@cached_function
def _one_D_from_rn(WR,t):
    r""" Find the D s.t. +-D/Level = n +- q(r)
    """
    #print "t=",t,type(t)
    if not isinstance(t,tuple):
        raise TypeError("Need a tuple of integers! Got:{0}".format(t))
    (r,n)=t
    sig=1
    #print "r=",r
    if r in WR.D():
        x=WR.Qv[r]
    #elif r in WR.D_as_integers:
    #    x=WR.Q(WR.D[r])
    else:
        raise TypeError("Need (r,n) in proper format forcoefficients! I.e. n integer and r in D or integer!")
    #print "x=",x
    if WR.is_dual():
        sig=-1
    D=sig*WR.level()*(n+sig*x)
    return D





def smallest_inf_norm_mpmath(V):
    r"""
    Computes the smallest of the supremum norms of the columns of a matrix.

    INPUT:

        - ``V`` -- matrix (real/complex)

    OUTPUT:

        - ``t`` -- minimum of supremum norms of the columns of V

    EXAMPLE::


        sage: A=mpmath.matrix([['0.1','0.3','1.0'],['7.1','5.5','4.8'],['3.2','4.4','5.6']])
        sage: smallest_inf_norm(A)
        mpf('5.5')

    
    """
    minc=mpmath.mpf(100)
    mi=0
    for j in range(V.cols):
        maxr=mpmath.mpf(0)
        for k in range(V.rows):
            t=abs(V[k,j])
            if(t>maxr):
                maxr=t
        if(maxr<minc):
            minc=maxr
            mi=j
    return minc

def matrix_norm(V,nt='max'):
    r"""
    Computes the smallest of the supremum norms of the columns of a matrix.

    INPUT:

        - ``V`` -- matrix (real/complex)
        -''nt'' -- string (default 'max') type of matrix norm

    OUTPUT:

        - ``t`` -- norm of V

    EXAMPLE::


        sage: A=mpmath.matrix([['0.1','0.3','1.0'],['7.1','5.5','4.8'],['3.2','4.4','5.6']])
        sage: smallest_inf_norm(A)
        mpf('5.5')

    
    """
    t=0
    ix=(0,0)
    if(nt=='max'):
        for j in range(V.cols):
            for k in range(V.rows):
                if(abs(V[k,j])>t):
                    t=abs(V[k,j])
                    ix=(k,j)
    elif(nt=='inf' or nt=='Inf'):
        # Max row sum 
        for n in range(V.rows):
            summa=0
            for l in range(V.cols):
                summa=summa+abs(V[n,l])
            if(summa>t):
                t=summa
    elif(nt=='L1' or nt=='l1'):
        # max column sum
        for l in range(V.cols):
            summa=0
            for n in range(V.rows):
                summa=summa+abs(V[n,l])
            if(summa>t):
                t=summa
    else:
        raise NotImplementedError
    return t,ix

def is_int(q):
    r"""
    Find out if the rational number q is an integer.

    INPUT:
    -''q'' -- integer/rational/real

    OUTPUT:
    - logical -- True if q is an integer otherwise False

    EXAMPLES::

        sage: is_int(1)
        True
        sage: is_int(float(1.0))
        True
        sage: is_int(RR(1.0))   
        True
        sage: is_int(6/3) 
        True
        sage: is_int(6/4)
        False
        sage: is_int(1.5)    
        False
        sage: is_int(Gamma0(1))    
        False
        
    """
    if(isinstance(q,Integer) or isinstance(q,int)):
        return True
    if(isinstance(q,Rational)):
        n=q.denominator()
        if(n==1):
            return True
    if(isinstance(q,tuple)):
        return False
    try:
        if(floor(q)==ceil(q)):
            return True
    except:
        pass
    return False




def get_list_of_forms(Nmin,Nmax,compute=False):
    r""" compute a whole list of various forms."""
    FL=dict()
    try:
        for N in range(Nmin,Nmax):
            print("N={0}".format(N))
            M=VVHarmonicWeakMaassForms(int(N),-0.5,75,dr=False,verbose=1)
            
            print("minPP={0}".format(M.smallest_pp()))
            print("Compute form on {0}".format(M))
            #s="FN"+str(N)+"-DR-D"+str(F.get_principal_part(str=True))+".sobj"
            s="FN"+str(N)+"-DR-D"+str(list(M.smallest_pp().keys())[0])+".sobj"
            print("trying :{0}".format(s))
            try:
                F=load(s)
            except IOError:
                # the file did not exist
                if(compute):
                    F=M.get_element(maxD=500)
                    save(F,s)
                else:
                    continue
            FL[N]=F
    except KeyboardInterrupt:
        pass
    return FL

def check_relevant_forms(L):
    r"""
    Filter out those forms which might correspond to newforms on
    Gamma0(N) 
    """
    L2=list()
    for F in L.values():
        #needed_ev=(
        S=F.shim_corr
        #print "S=",S
        ok_ev=0
        for g in S:
            if g.atkin_lehner_eigenvalue() == -1:
                ok_ev=ok_ev+1
        if ok_ev > 0:
            print("Number of ok forms on ",F.space.WR.N," :",ok_ev)
            F.list_coefficents('D',fd=True,neg=False,nmin=0,nmax=1000,latex=False,nd=50,prime=True)
            L2.append(F)
    return L2
#__main__.VVHarmonicWeakMaassForms=VVHarmonicWeakMaassForms
#__main__.WeilRepDiscriminantForm=WeilRepDiscriminantForm
#__main__.HarmonicWeakMaassFormElement=HarmonicWeakMaassFormElement

