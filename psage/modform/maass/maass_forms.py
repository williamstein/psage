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


#from sage.all_cmdline import *   # import sage library
#import mpmath as mpmath
#import mpmath as mpmath
from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.structure.sage_object import SageObject,cPickle
#from mpmath import mpf
from sage.functions.all import ln,sqrt,floor
from sage.rings.arith import divisors,gcd
from sage.modular.dirichlet import DirichletGroup
from sage.rings.all import RR
# needed for attaching the file locally
##import __main__
from sage.modular.arithgroup.all import Gamma0

from maass_forms_alg import *
r"""
Maass waveforms for subgroups of the modular group

AUTHORS:

 - Fredrik Strömberg (March 2010)

EXAMPLES::


 ?


 TODO:
   - Nontrivial multiplier systems and weights
   - improve eigenvalue finding alorithms


"""


#load "maass_forms_alg.pyx"


class MaassWaveForms (Parent):
    r"""
    Describes a space of Maass waveforms
    """
    def __init__(self,G,prec=500,verbose=None):
        r"""
        Creates an ambient space of Maass waveforms
        (i.e. there are apriori no members).

        INPUT:

            - `'G``    -- subgroup of the modular group
            - ``prec`` -- default working precision in bits 

        EXAMPLES::


            sage: S=MaassWaveForms(Gamma0(1)); S
            Space of Maass waveforms on the group G:
            Arithmetic Subgroup of PSL2(Z) with index 1. Given by: 
            perm(S)=()
            perm(ST)=()
            Constructed from G=Modular Group SL(2,Z)

        """
        if(not isinstance(G,sage.modular.maass.mysubgroup.MySubgroup)):
            if(isinstance(G,int) or isinstance(G,sage.rings.integer.Integer)):
                self._G=MySubgroup(Gamma0(G))
            else:
                try:
                    self._G=MySubgroup(G)
                except TypeError:
                    raise TypeError,"Incorrect input!! Need subgroup of PSL2Z! Got :%s" %(G)
        else:
            self._G=G

        self.verbose=verbose
        self.prec=prec
        self._Weyl_law_const=self._Weyl_law_consts()
        maass_forms=dict() # list of members 
        if(verbose==None):
            self._verbose=0 # output debugging stuff or not
        else:
            self._verbose=verbose
        self._symmetry=None

    def _repr_(self):
        r"""
        Return the string representation of self.

        EXAMPLES::


            sage: M=MaassWaveForms(MySubgroup(Gamma0(1)));M
            Space of Maass waveforms on the group G:
            Arithmetic Subgroup of PSL2(Z) with index 1. Given by:
                perm(S)=()
                perm(ST)=()
            Constructed from G=Modular Group SL(2,Z)


        """
        s="Space of Maass waveforms on the group G:\n"+str(self._G)
        return s

    def __reduce__(self):
        r""" Used for pickling.
        """
        return(MaassWaveForms,(self._G,self.prec,self.verbose))

    def __cmp__(self,other):
        r""" Compare self to other
        """
        if(self._G <> other._G or self.prec<>other.prec):
            return False
        else:
            return True
        
    def get_element_in_range(self,R1,R2,ST=None,neps=10):
        r""" Finds element of the space self with R in the interval R1 and R2

        INPUT:

        - ``R1`` -- lower bound (real)
        - ``R1`` -- upper bound (real)
        
        """
        if(ST <>None):
            ST0=ST
        else:
            ST0=self._ST
        l=self.split_interval(R1,R2)
        if(self._verbose>1):
            print "Split into intervals:"
            for [r1,r2,y] in l:
                print "[",r1,",",r2,"]:",y
        Rl=list()
        for [r1,r2,y] in l:
            [R,er]=find_single_ev(self,r1,r2,Yset=y,ST=ST0,neps=neps)
            Rl.append([R,er])
        print "R=",R
        print "er=",er


    def _Weyl_law_consts(self):
        r"""
        Compute constants for the Weyl law on self._G

        OUTPUT:

        - tuple of real numbers

        EXAMPLES::


            sage: M=MaassWaweForms(MySubgroup(Gamma0(1))
            sage: M._Weyl_law_consts  
            (0, 2/pi, (log(pi) - log(2) + 2)/pi, 0, -2)
        """
        import mpmath.fp
        pi=mpmath.fp.pi
        ix=Integer(self._G.index())
        nc=Integer(len(self._G.cusps()))
        if(self._G.is_congruence()):
            lvl=Integer(self._G.level())
        else:
            lvl=0
        n2=Integer(self._G.nu2())
        n3=Integer(self._G.nu3())
        c1=ix/Integer(12)
        c2=Integer(2)*nc/pi
        c3=nc*(Integer(2)-ln(Integer(2))+ln(pi))/pi
        if(lvl<>0):
            A=1
            for q in divisors(lvl):
                num_prim_dc=0
                DG=DirichletGroup(q)
                for chi in DG.list():
                    if(chi.is_primitive()):
                        num_prim_dc=num_prim_dc+1
                for m in divisors(lvl):
                    if(lvl % (m*q) == 0   and m % q ==0 ): 
                        fak=(q*lvl)/gcd(m,lvl/m)
                        A=A*Integer(fak)**num_prim_dc        
            c4=-ln(A)/pi
        else:
            c4=Integer(0)
        # constant term
        c5=-ix/144+n2/8+n3*2/9-nc/4-1
        return (c1,c2,c3,c4,c5)

    def Weyl_law_N(self,T,T1=None):
        r"""
        The counting function for this space. N(T)=#{disc. ev.<=T}
        
        INPUT:
        
        -  ``T`` -- double


        EXAMPLES::

            sage: M=MaassWaveForms(MySubgroup(Gamma0(1))
            sage: M.Weyl_law_N(10)
            0.572841337202191
            
        """
        (c1,c2,c3,c4,c5)=self._Weyl_law_const
        cc1=RR(c1); cc2=RR(c2); cc3=RR(c3); cc4=RR(c4); cc5=RR(c5)
        #print "c1,c2,c3,c4,c5=",cc1,cc2,cc3,cc4,cc5
        t=sqrt(T*T+0.25)
        try: 
            lnt=ln(t)
        except TypeError:
            lnt=mpmath.ln(t)
        #print "t,ln(t)=",t,lnt
        NT=cc1*t*t-cc2*t*lnt+cc3*t+cc4*t+cc5
        if(T1<>None):
            t=sqrt(T1*T1+0.25)
            NT1=cc1*(T1*T1+0.25)-cc2*t*ln(t)+cc3*t+cc4*t+cc5
            return RR(abs(NT1-NT))
        else:
            return RR(NT)

    def next_eigenvalue(self,R):
        r"""
        An estimate of where the next eigenvlue will be, i.e. the smallest R1>R so that N(R1)-N(R)>=1

        INPUT:
        - ``R`` -- real > 0

        OUTPUT:
        - real > R

        EXAMPLES::

            sage: M.next_eigenvalue(10.0)
            12.2500000000000

        
        """
        #cdef nmax
        N=self.Weyl_law_N(R)
        try:
            for j in range(1,10000):
                R1=R+j*RR(j)/100.0
                N1=self.Weyl_law_N(R1)
                if(N1-N >= 1.0):
                    raise StopIteration()
        except StopIteration:
            return R1
        else:
            raise ArithmeticError,"Could not find next eigenvalue! in interval: [%s,%s]" %(R,R1)
        
    def Weyl_law_Np(self,T,T1=None):
        r"""
        The derviative of the counting function for this space. N(T)=#{disc. ev.<=T}
        
        INPUT:
        
        -  ``T`` -- double


        EXAMPLES::

            sage: M=MaassWaweForms(MySubgroup(Gamma0(1))
            sage: M.Weyl_law_Np(10)
        
        """
        (c1,c2,c3,c4,c5)=self._Weyl_law_const
        cc1=RR(c1); cc2=RR(c2); cc3=RR(c3); cc4=RR(c4); cc5=RR(c5)
        #print "c1,c2,c3,c4,c5=",c1,c2,c3,c4,c5
        NpT=2.0*cc1*T-cc2*(ln(T)+1.0)+cc3+cc4
        return RR(NpT)


    #### Split an interv
    def split_interval(self,R1,R2):
        r"""
        Split an interval into pieces, each containing (on average) at most one
        eigenvalue as well as a 0<Y<Y0 s.t. K_IR(Y) has no zero here

        INPUT:
        
            - ''R1'' -- real
            - ''R2'' -- real

        OUPUT:

            - list of triplets (r1,r2,y) where [r1,r2] does not contain a zero of K_ir(y)

        EXAMPLES::


            sage: M._next_kbes
            sage: l=M.split_interval(9.0,11.0)
            sage: print l[0],'\n',l[1],'\n',l[2],'\n',l[3]
            (9.00000000000000, 9.9203192604549457, 0.86169527676551638)
            (9.9203192704549465, 10.135716354265259, 0.86083358148875089)
            (10.13571636426526, 10.903681677771321, 0.86083358148875089)
            (10.903681687771321, 11.0000000000000, 0.85997274790726208)
            
        """
        # It is enough to work with double precision
        base=mpmath.fp
        pi=base.pi
        # Locate all zeros of K_IR(Y0) first
        #def f(r):
        #    ir=base.mpc(0 ,r)
        #    return base.besselk(ir,Y0)
        # First we find the next zero 
        # First split into intervals having at most one zero
        ivs=list()
        rnew=R1; rold=R1
        while (rnew < R2):
            rnew=min(R2,self.next_eigenvalue(rold))
            if( abs(rold-rnew)==0.0):
                print "ivs=",ivs
                exit
            iv=(rold,rnew)
            ivs.append(iv)
            rold=rnew

        # We now need to split these intervals into pieces with at most one zero of the K-Bessel function
        Y00=base.mpf(0.995)*base.sqrt(base.mpf(3))/base.mpf(2 *self._G._level)
        new_ivs=list()
        for (r1,r2) in ivs:
            print "r1,r2=",r1,r2
            Y0=Y00; r11=r1
            i=0
            while(r11 < r2 and i<1000):
                t=self._next_kbessel_zero(r11,r2,Y0*pi);i=i+1
                print "t=",t
                iv=(r11,t,Y0); new_ivs.append(iv)
                # must find Y0 s.t. |besselk(it,Y0)| is large enough
                Y1=Y0
                k=base.besselk(base.mpc(0,t),Y1).real*mpmath.exp(t*0.5*base.pi)
                j=0
                while(j<1000 and abs(k)<1e-3):
                    Y1=Y1*0.999;j=j+1
                    k=base.besselk(base.mpc(0,t),Y1).real*mpmath.exp(t*0.5*base.pi)
                Y0=Y1
                r11=t+1E-08
        return new_ivs


    def _next_kbessel_zero(self,r1,r2,y):
        r"""
        The first zero after r1 i the interval [r1,r2] of K_ir(y),K_ir(2y)
        
        INPUT:

            - ´´r1´´ -- double
            - ´´r2´´ -- double
            - ´´y´´  -- double
    
        OUTPUT:

            - ''double''

        EXAMPLES::

        
            sage: M=MaassWaveForms(MySubgroup(Gamma0(1)))
            sage: Y0=0.995*sqrt(3.0)/2.0
            sage: M._next_kbessel_zero(9.0,15.0,Y0)
            9.9203192604549439
            sage: M._next_kbessel_zero(9.921,15.0,Y0)
            10.139781183668587


        CAVEAT:
            The rootfinding algorithm is not very sophisticated and might miss roots
            
        """
        base=mpmath.fp
        h=(r2-r1)/500.0
        t1=-1.0; t2=-1.0
        r0=base.mpf(r1)
        kd0=my_mpmath_kbes_diff_r(r0,y,base)
        #print "r0,y=",r0,y,kd0
        while(t1<r1 and r0<r2):

            # Let us first find a R-value for which the derivative changed sign
            kd=my_mpmath_kbes_diff_r(r0,y,base)
            i=0
            while(kd*kd0>0 and i<500 and r0<r2):
                i=i+1
                r0=r0+h
                kd=my_mpmath_kbes_diff_r(r0,y,base)
                #print "r0,kd=",r0,kd
                #print "kd*kd0=",kd*kd0
            #print "-r0,y,kd=",r0,y,kd
            #t1=base.findroot(lambda x :  base.besselk(base.mpc(0,x),base.mpf(y),verbose=True).real,r0)
            try:
                t1=base.findroot(lambda x :  my_mpmath_kbes(x,y,base),r0)
            except ValueError:
                t1=base.findroot(lambda x :  my_mpmath_kbes(x,y,mpmath.mp),r0)
            r0=r0+h
        if(r0>=r2 or t1>=r2):
            t1=r2
        r0=r1
        kd0=my_mpmath_kbes_diff_r(r0,y,base)
        while(t2<r1 and r0<r2):
            kd=my_mpmath_kbes_diff_r(r0,y,base)
            i=0
            while(kd*kd0>0 and i<500 and r0<r2):
                i=i+1
                r0=r0+h
                kd=my_mpmath_kbes_diff_r(r0,y,base)
            try:
                t2=base.findroot(lambda x :  my_mpmath_kbes(x,2*y,base),r0)
            except ValueError:
                t2=base.findroot(lambda x :  my_mpmath_kbes(x,2*y,mpmath.mp),r0)
            #t2=base.findroot(lambda x :  base.besselk(base.mpc(0,x),base.mpf(2*y),verbose=True).real,r0)
            r0=r0+h
        if(r0>=r2 or t2>=r2):
            t2=r2
            #print "zero(besselk,y1,y2)(",r1,r2,")=",t1,t2
        t=min(min(max(r1,t1),max(r1,t2)),r2)
        return t



def my_mpmath_kbes(r,x,mp_ctx=None):
    r"""Scaled K-Bessel function with

    INPUT:

        - ''r'' -- real
        - ''x'' -- real
        - ''mp_ctx'' -- mpmath context (default None)

    OUTPUT:

        - real -- K_ir(x)*exp(pi*r/2)

    EXAMPLES::


        sage: my_mpmath_kbes(9.0,1.0)
        mpf('-0.71962866121965863')
        sage: my_mpmath_kbes(9.0,1.0,mpmath.fp)
        -0.71962866121967572

    

    """
    import mpmath.mp,mpmath.fp
    if(mp_ctx==None):
        mp_ctx=mpmath.mp
    if(mp_ctx==mpmath.mp):
        pi=mpmath.mp.pi()
    else:
        pi=mpmath.fp.pi
    try:
        k=mp_ctx.besselk(ctx.mpc(0,r),ctx.mpf(x))
        f=k*mp_ctx.exp(r*ctx.mpf(0.5)*pi)
    except OverflowError:
        k=mp_cyx.besselk(mp_ctx.mpc(0,r),mp_ctx.mpf(x))
        f=k*mp_ctx.exp(r*mp_ctx.mpf(0.5)*pi)
    return f.real

def my_mpmath_kbes_diff_r(r,x,mp_ctx=None):
    r"""
    Derivative with respect to R of the scaled K-Bessel function.

    INPUT:

        - ''r'' -- real
        - ''x'' -- real
        - ''ctx'' -- mpmath context (default mpmath.mp)

    OUTPUT:

        - real -- K_ir(x)*exp(pi*r/2)

    EXAMPLES::


        sage: my_mpmath_kbes_diff_r(9.45,0.861695276766 ,mpmath.fp)
        -0.31374673969963851
        sage: my_mpmath_kbes_diff_r(9.4,0.861695276766 ,mpmath.fp)
        0.074219541623676832


    """
    import mpmath.mp,mpmath.fp
    if(mp_ctx==None):
        mp_ctx=mpmath.mp
    if(mp_ctx==mpmath.mp):
        pi=mpmath.mp.pi()
    else:
        pi=mpmath.fp.pi
    try:
        k=mp_ctx.besselk(ctx.mpc(0,r),ctx.mpf(x))
        f=k*mp_ctx.exp(r*ctx.mpf(0.5)*pi)
    except OverflowError:
        k=mp_ctx.besselk(mp_ctx.mpc(0,r),mp_ctx.mpf(x))
        f=k*mp_ctx.exp(r*mp_ctx.mpf(0.5)*pi)
    f1=f.real
    try:
        h=mp_ctx.mpf(1e-8)
        k=mp_ctx.besselk(mp_ctx.mpc(0,r+h),mp_ctx.mpf(x))
        f=k*mp_ctx.exp((r+h)*mp_ctx.mpf(0.5)*pi)
    except OverflowError:
        h=mp_ctx.mpf(1e-8)
        k=mp_ctx.besselk(mp_ctx.mpc(0,r+h),mp_ctx.mpf(x))
        f=k*mp_ctx.exp((r+h)*mp_ctx.mpf(0.5)*pi)
    f2=f.real
    diff=(f2-f1)/h
    return diff

#class MaassWaveformElement (SageObject):
class MaassWaveformElement (Parent):
    r"""
    An element of a space of Maass waveforms

    EXAMPLES::

        sage: G=MySubgroup(Gamma0(1))
        sage: R=mpmath.mpf(9.53369526135355755434423523592877032382125639510725198237579046413534)
        sage: F=MaassWaveformElement(G,R)    
        Maass waveform with parameter R=9.5336952613536
        in Space of Maass waveforms on the group G:
        Arithmetic Subgroup of PSL2(Z) with index 1. Given by:
            perm(S)=()
            perm(ST)=()
        Constructed from G=Modular Group SL(2,Z)
        sage: G=MySubgroup(Gamma0(4))
        sage: R=mpmath.mpf(3.70330780121891)
        sage: F=MaassWaveformElement(G,R);F
        Maass waveform with parameter R=3.70330780121891
        Member of the Space of Maass waveforms on the group G:
        Arithmetic Subgroup of PSL2(Z) with index 6. Given by: 
        perm(S)=(1,2)(3,4)(5,6)
        perm(ST)=(1,3,2)(4,5,6)
        Constructed from G=Congruence Subgroup Gamma0(4)
        sage: F.C(0,-1)
        mpc(real='-1.0000000000014575', imag='5.4476887980094281e-13')
        sage: F.C(0,15)-F.C(0,5)*F.C(0,3)
        mpc(real='-5.938532886679327e-8', imag='-1.0564743382278074e-8')
        sage: F.C(0,3)
        mpc(real='0.53844676975670527', imag='-2.5525466782958545e-13')
        sage: F.C(1,3)
        mpc(real='-0.53844676975666916', imag='2.4484251009604091e-13')
        sage: F.C(2,3)
        mpc(real='-0.53844676975695485', imag='3.3624257152434837e-13')
        
        

        
    
    """
    def __init__(self,G,R,C=None,nd=12,sym_type=None,verbose=0):
        r"""
        Construct a Maass waveform on thegroup G with spectral parameter R and coefficients C
        
        INPUT:
        - ``G`` -- Group
        - ``R`` -- Spectral parameter
        - ``C`` -- Fourier coefficients (default None)
        - ``nd``-- Number of desired digits (default 15)
        """
        import mpmath.mp,mpmath.fp
        self._space= MaassWaveForms (G)
        self._group=self._space._G
        self._R=R
        self._nd=nd
        self._sym_type=sym_type
        self._space._verbose=verbose
        if(nd>15):
            mpmath.mp.dps=nd
            self.mp_ctx=mpmath.mp
        else:
            self.mp_ctx=mpmath.fp
        ## We use the Fourier coefficients to verify whether we really have an eigenvalue
        if(C<>None):
            self.coeffs=C
            er=self.test(C)
        else:
            #print "Need to compute coefficients!"
            # We compute a set of Fourier coefficients
            #(Y,M)=find_Y_and_M(G,R)
            #Q=M+10
            self.coeffs=Maassform_coeffs(self._space,R,ndigs=self._nd)[0]


    def _repr_(self):
        r""" Returns string representation of self.

        EXAMPLES::


            sage: R=mpmath.mpf(9.53369526135355755434423523592877032382125639510725198237579046413534)
            sage: F=MaassWaveformElement(Gamma0(1),R,nd=50);F
            Maass waveform with parameter R=9.5336952613535575543442352359287703238212563951073
            Member of the Space of Maass waveforms on the group G:
            Arithmetic Subgroup of PSL2(Z) with index 1.Given by
                perm(S)=()
                perm(ST)=()
            Constructed from G=Modular Group SL(2,Z)

        
        """
        s="Maass waveform with parameter R="+str(self._R)+"\nMember of the "+str(self._space)
        return s


    def __reduce__(self):
        r""" Used for pickling.
        """
        return(MaassWaveformElement,(self._group,self._R,self.coeffs,self._nd))

    def group():
        r"""
        Return self._group 
        """
        return self._group

    def eigenvalue():
        return self._R  #eigenvalue

    def C(self,i,j=None):
        r"""
        Return the coefficient C(i,j) i.e. coefficient nr. j at cusp i (or if j=none C(1,i))


        EXAMPLES::

        
            sage: G=MySubgroup(Gamma0(1))
            sage: R=mpmath.mpf(9.53369526135355755434423523592877032382125639510725198237579046413534)
            sage: F=MaassWaveformElement(G,R)    
            sage: F.C(2)
            mpc(real='-1.068333551223568', imag='2.5371356217909904e-17')
            sage: F.C(3)
            mpc(real='-0.45619735450601293', imag='-7.4209294760716175e-16')
            sage: F.C(2)*F.C(3)-F.C(6)
            mpc(real='8.2470016583210667e-8', imag='1.6951583479643061e-9')


        """
        if(j==None):
            cusp=0
            j=i
        else:
            cusp=i
        if(not self.coeffs.has_key(cusp)):
                raise ValueError," Need a valid index of a cusp as first argument! I.e in %s" %self.coeffs.keys()
        if(not self.coeffs[cusp].has_key(j)):
            return None
        return self.coeffs[cusp][j]
        
    def test(self,do_twoy=False,up_to_M=0):
        r""" Return the number of digits we believe are correct (at least) 

        EXAMPLES::

        
            sage: G=MySubgroup(Gamma0(1))
            sage: R=mpmath.mpf(9.53369526135355755434423523592877032382125639510725198237579046413534)
            sage: F=MaassWaveformElement(G,R)    
            sage: F.test()
            7



        """
        # If we have a Gamma_0(N) we can use Hecke operators
        if(not do_twoy and (isinstance(self._space,sage.modular.arithgroup.congroup_gamma0.Gamma0_class) or \
            self._space._G.is_congruence)):
            # isinstance(self._space,sage.modular.arithgroup.congroup_sl2z.SL2Z_class_with_category))):
            if(self._space._verbose>1):
                print "Check Hecke relations!"
            er=test_Hecke_relations(2,3,self.coeffs)
            d=floor(-mpmath.log10(er))
            if(self._space._verbose>1):
                print "Hecke is ok up to ",d,"digits!"
            return d
        else:
            # Test two Y's
            nd=self._nd+5
            [M0,Y0]=find_Y_and_M(G,R,nd)
            Y1=Y0*0.95
            C1=Maassform_coeffs(self,R,Mset=M0,Yset=Y0 ,ndigs=nd )[0]
            C2=Maassform_coeffs(self,R,Mset=M0,Yset=Y1 ,ndigs=nd )[0]
            er=mpmath.mpf(0)
            for j in range(2,max(M0/2,up_to_M)):
                t=abs(C1[j]-C2[j])
                print "|C1-C2[",j,"]|=|",C1[j],'-',C2[j],'|=',t
                if(t>er):
                    er=t
            d=floor(-mpmath.log10(er))
            print "Hecke is ok up to ",d,"digits!"
            return d

def Maassform_coeffs(S,R,Mset=0 ,ST=None ,ndigs=12,twoy=None):
    r"""
    Compute M Fourier coefficients (at each cusp) of a Maass (cusp)
    waveform with eigenvalue R for the group G.


    INPUT:
    
    - ''S''    -- space of Maass waveforms
    - ''R''    -- Real number with precision prec
    - ''Mset'' -- integer : number of desired coefficients
    - ''ST''   -- set symmetry
    - ''Yset'' -- real
    - ''ndigs''-- integer


    OUTPUT:

    -''D'' -- dictionary of Fourier coefficients

    EXAMPLES:


        sage: R=mpmath.mpf(9.53369526135355755434423523592877032382125639510725198237579046413534)
        sage:         sage: M=MaassWaveForms(Gamma0(1))
        sage:         sage: C=Maassform_coeffs(M,R)


    """
    G=S._G
    if(S._verbose>1):
        print "S=",S
    [YY,M0]=find_Y_and_M(G,R,ndigs)
    if(ndigs>=15):
        Y=mpmath.mp.mpf(YY)
    else:
        Y=YY
    #print "M,Y=",M0,Y0
    #if(Yset<>None  and Yset<=Y0):
    #    Y=Yset
    #else:
    #    Y=Y0
    if(Mset>=M0):
        M=Mset
    else:
        M=M0
    if(ST<>None):
        if(ST.has_key('sym_type')):
            sym_type=ST['sym_type']
        else:
            sym_type=False
    else:
        sym_type=None
    Q=M+10 
    dold=mpmath.mp.dps
    ###mpmath.mp.dps=ndigs+2  # We work with low precision initially
    if(S._verbose>1):
        print "R,Y,M,Q=",R,Y,M,Q
    #print "sym_type=",sym_type
    #print "nd=",mpmath.mp.dps
    X = coefficients_for_Maass_waveforms(S,R,Y,M,Q,ndigs,cuspidal=True,sym_type=sym_type)
    return X



def coefficients_for_Maass_waveforms(S,R,Y,M,Q,ndigs,cuspidal=True,sym_type=None):
    r"""
    Compute coefficients of a Maass waveform given a specific M and Y.
     INPUT:
    
    - ''S''    -- Space of Maass waveforms
    - ''R''    -- real : 1/4+R*R is the eigenvalue
    - ''Y''    -- real number > 9
    - ''M''    -- integer
    - ''Q''    -- integer        
    - ''cuspidal''-- logical (default True)
    - ''sym_type'' -- integer (default None)
    - ''ndigs''-- integer : desired precision

    OUTPUT:

    -''D'' -- dictionary of Fourier coefficients


    EXAMPLES::

        sage: S=MaassWaveForms(Gamma0(1))
        sage: R=mpmath.mpf(9.53369526135355755434423523592877032382125639510725198237579046413534)
        sage: Y=mpmath.mpf(0.85)
        sage: C=coefficients_for_Maass_waveforms(S,R,Y,10,20,12)

    

    """
    G=S._G
    if(ndigs<=12):
        W=setup_matrix_for_Maass_waveforms(G,R,Y,M,Q,cuspidal=True,sym_type=sym_type,low_prec=True)
    else:
        W=setup_matrix_for_Maass_waveforms(G,R,Y,M,Q,cuspidal=True,sym_type=sym_type)
    dim=1
    N=set_norm_maass(dim,cuspidal=True)
    #return [W, N]
    dold=mpmath.mp.dps
    mpmath.mp.dps=max(dold,50)  # We work with low precision initially
    if(S._verbose>1):
        deb=True
    else:
        deb=False
    done=False; j=0
    while(done==False and j<=10):
        if(S._verbose>1):
            print "Trying to solve with prec=",mpmath.mp.dps
        try:
            X=solve_system_for_Maass_waveforms(W,N,deb=deb)
        except ZeroDivisionError:
            pass
        if(is_Integer(X) or isinstance(X,int)):
            mpmath.mp.dps=X+5
        elif(isinstance(X,dict)):
            done=True
        else:
            raise ArithmeticError," Could not solve system!"
        j=j+1
    if(S._verbose>1):
        for m in range(dim):
            print "Function nr. ",m+1 
            if(sym_type==None):
                for n in range(M,1 ,-1 ): 
                    print "C[",n,"]=",X[m][n]
                for n in range(M): 
                    print "C[",n,"]=",X[m][n]
            else:
                for n in range(1,M+1): 
                    print "C[",n,"]=",X[m][n]
    #print "C2=",X[0][2]
    #print "c2c3-c6=",X[0][2]*X[0][3]-X[0][6]
    mpmath.mp.dps=dold
    return X

def verify_eigenvalue(S,R,nd=10,ST=None,method='TwoY'):
    r""" Verify an eigenvalue and give an estimate of the error.

    INPUT:
    -''S'' -- Space of Maass waveforms
    -''R'' -- real: (tentative) eigenvalue = 1/4+R**2 
    -''nd''-- integer : number of digits we try to get (at a minimum)
    """
    C=Maassform_coeffs(S,R,ST=ST ,ndigs=nd)





def find_single_ev(S,R1in,R2in,Yset=None,ST=None,neps=10,method='TwoY',verbose=0):
    r""" Locate a single eigenvalue on G between R1 and R2

    INPUT:(tentative)

    - ''S''    -- space of Maass waveforms
    - ''R1in'' -- real
    - ''R1in'' -- real
    - ''Yset'' -- real (use this value of Y to compute coefficients)
    - ''ST'''  -- dictionary giving symmetry
    - ''neps'' -- number of desired digits

    OUPUT:
    
    - ''R'' --

    
    """
    G=S._G
    jmax=1000  # maximal number of interation
    if(neps>=15):
        R1=mpmath.mp.mpf(R1in);R3=mpmath.mp.mpf(R2in)
        print "mpmath.mp.dps=",mpmath.mp.dps
        print "R1=",R1,type(R1)
        print "R3=",R3,type(R3)
    else:
        R1=mpmath.fp.mpf(R1in);R3=mpmath.fp.mpf(R2in)
    if(Yset==None):
        [Y,M]=find_Y_and_M(G,R1,neps)
    else:
        [Y,M]=find_Y_and_M(G,R1,neps,Yset=Yset)
    Y1=Y; Y2=mpmath.mpf(0.995)*Y1
    tol=mpmath.mpf(10)**mpmath.mpf(-neps)
    dold=mpmath.mp.dps
    mpmath.mp.dps=neps+3  # We work with low precision initially
    h=dict()
    signs=dict();diffs=dict()
    c=dict(); h=dict()
    c[1]=2 ; c[2 ]=3 ; c[3 ]=4 
    #met='Hecke'

    met=method
    [diffs[1 ],h[1 ]]=functional(S,R1,M,Y1,Y2,signs,c,ST,first_time=True,method=met,ndigs=neps)
    [diffs[3 ],h[3 ]]=functional(S,R3,M,Y1,Y2,signs,c,ST,first_time=True,method=met,ndigs=neps)
    if(S._verbose>1):
        print "diffs: met=",met
        print "R1=",R1
        print "R3=",R3
        for n in list(c.keys()): #.sort():
            for j in list(diffs.keys()): #.sort():
                print "diff[",j,c[n],"]=",diffs[j][n]
    # Sset signs and check zeros
    if(met=='Hecke'):
        if(h[1 ]*h[3]>mpmath.eps()):
                # We do not have a signchange
                return [0 ,0 ]
    else:
        var=0.0
        for j in range(1 ,3 +1 ):
            var+=abs(diffs[1 ][j])+abs(diffs[3 ][j])
        print "var=",var
        for j in range(1 ,3 +1 ):
            signs[j]=1 
            if(diffs[1 ][j]*diffs[3 ][j]>mpmath.eps()):
                # If we do not have a signchange
                # and the absolute values are relatively large
                # there is probably no zero here
                if(abs(diffs[1][j])+abs(diffs[3][j]) > 0.01*var):
                    return [0 ,0 ]
            elif(diffs[1 ][j]>0 ):
                signs[j]=-1         
            # Recompute functionals using the signs
        if(S._verbose>1):
            print "h1=",h
            print "diffs1=",diffs
            print "signs=",signs
        for k in [1,3]:
            h[k]=0 
            for j in range(1,3+1):
                h[k]=h[k]+signs[j]*diffs[k][j]
    Rnew=prediction(h[1 ],h[3 ],R1,R3)
    if(S._verbose>1):
        print "h=",h
        print "Rnew=",Rnew
    [diffs[2],h[2]]=functional(S,Rnew,M,Y1,Y2,signs,c,ST,first_time=False,method=met,ndigs=neps)
    zero_in=is_zero_in(h)
    if(zero_in == -1 ):
        R3=Rnew; h[3]=h[2 ]; diffs[3 ]=diffs[2 ]; errest=abs(Rnew-R1)
    else:
        R1=Rnew; h[1 ]=h[2 ]; diffs[1 ]=diffs[2 ]; errest=abs(Rnew-R3)
    step=0
    for j in range(100):
        Rnew=prediction(h[1 ],h[3 ],R1,R3)
        errest=max(abs(Rnew-R1),abs(Rnew-R3))
        if(S._verbose>1):
            print "R1,R3,Rnew,errest=",R1,R3,Rnew,errest
        if(errest<tol):
            return [Rnew,errest]
        [diffs[2 ],h[2 ]]=functional(S,Rnew,M,Y1,Y2,signs,c,ST,first_time=False,method=met,ndigs=neps)
        zero_in=is_zero_in(h)
        if(zero_in==0):
            return [Rnew,errest]
        elif(zero_in not in [1,-1]):
            raise StopIteration()
        if(zero_in==-1):
            stepz=abs(Rnew-R3)
            R3=Rnew; h[3 ]=h[2 ]; diffs[3 ]=diffs[2 ]; errest=abs(Rnew-R1)
        elif(zero_in==1):
            stepz=abs(Rnew-R1)
            R1=Rnew; h[1 ]=h[2 ]; diffs[1 ]=diffs[2 ]; errest=abs(Rnew-R3)
        # If we have gone in the same direction too many times we need to modify our approach
        step=step+zero_in
        if(S._verbose>1):
            print "step=",step
        if(step>2):    # Have gone too many times to the left
            Rtest=Rnew + mpmath.mpf(0.5)*stepz  # Need to test if this modified R3 work:
            if(S._verbose>1):
                print "Rtest(R)=",Rtest
            [diffs[2 ],h[2 ]]=functional(S,Rtest,M,Y1,Y2,signs,c,ST,False,met,neps)
            if(is_zero_in(h) ==-1): # all is ok
                R3=Rtest; h[3]=h[2]; diffs[3]=diffs[2]; step=step-1
            else: # Test another one
                Rtest=Rnew + mpmath.mpf(0.5)*abs(R1-R3)  # Need to test if this modified R3 work:
                if(S._verbose>1):
                    print "Rtest(R)=",Rtest
                [diffs[2 ],h[2 ]]=functional(G,Rtest,M,Y1,Y2,signs,c,ST,False,met,neps)
                if(is_zero_in(h) ==-1): # all is ok
                    R3=Rtest; h[3]=h[2]; diffs[3]=diffs[2]; step=step-1
        elif(step<-2): # Have gone too many times to the right
            Rtest=Rnew - mpmath.mpf(0.5)*stepz
            if(S._verbose>1):
                print "Rtest(L)=",Rtest
            [diffs[2 ],h[2 ]]=functional(S,Rtest,M,Y1,Y2,signs,c,ST,False,met,neps)
            if(is_zero_in(h) == 1): # all is ok
                R1=Rtest; h[1]=h[2]; diffs[1]=diffs[2]; step=step+1
            else:
                Rtest=Rnew - mpmath.mpf(0.5)*abs(R3-R1)
                if(S._verbose>1):
                    print "Rtest(L)=",Rtest
                [diffs[2 ],h[2 ]]=functional(S,Rtest,M,Y1,Y2,signs,c,ST,False,met,neps)
                if(is_zero_in(h) == 1): # all is ok
                    R1=Rtest; h[1]=h[2]; diffs[1]=diffs[2]; step=step+1
                
####
            
def is_zero_in(h):
    r"""
    Tells which interval contains changes of sign.
    

    INPUT:
    
    - ''h'' -- dictionary h[j]=[h1,h2,h3]

    OUTPUT:
    
    - integer  (-1 0 1)

    EXAMPLES::
    
    
    """
    zi=dict(); i=0
    if(h[1]*h[2] < 0):
        zi[-1]=1; i=-1
    if(h[3]*h[2] < 0):
        zi[1]=1; i=1
    if(zi.values().count(1) >1 ): # need to split
        return -2 
    #s="Neeed to split! Not implemented!"
    #raise ValueError,s
    return i

def prediction(f0,f1,x0,x1):
    r"""
    Predict zero using the secant method.

    INPUT:

        - ''f0'' -- real
        - ''f1'' -- real
        - ''x0'' -- real
        - ''x1'' -- real

    OUTPUT:

        - real


    EXAMPLES::


        sage: prediction(-1,1,9,10)
        19/2
        sage: prediction(-1.0,1.0,9.0,10.0)
        9.50000000000000

    
    """
    xnew=x0-f0*(x1-x0)/(f1-f0)
    if(xnew<x0 or xnew>x1):
        st= "Secant method ended up outside interval! \n"
        st+="input: f0,f1,x0,x1=%s,%s,%s,%s \n xnew=%s"
        raise ValueError,st%(f0,f1,x0,x1,xnew)
    return xnew

def prediction_newton(x,f,df):
    r"""
    Predict zero using the secant method.

    INPUT:

        - ''x''  -- real
        - ''f''  -- real, f(x)
        - ''df'' -- real, f'(x)

    OUTPUT:

        - real


    EXAMPLES::


        sage: prediction_newton(1,9,10)
        19/2
        sage: prediction_newton(1.0,9.0,10.0)
        9.50000000000000

    
    """
    if(df==0.0):
        st= "Newtons method failed! \n"
        st+="f'(x)=0. input: f,df,x=%s,%s,%s"
        raise ValueError,st%(f,df,x)
    xnew=x-f/df
    return xnew


def functional(S,r,M,Y1,Y2,signs,c,ST,first_time=False,method='Hecke',ndigs=12):
    r"""
    Computes the functional we use as an indicator of an eigenvalue.

    INPUT:

        -''S'' -- space of Maass waveforms
        -''r'' -- real
        -''M'' -- integer
        -''Y1''-- real
        -''Y2''-- real
        -''signs'' -- dict
        -''c'' -- set which coefficients to use
        -''ST''-- normalization
        -''first_time'' --
        -''method'' -- Hecke/two y
        -''ndigs'' -- integer (number of digits wanted)

    OUTPUT:

    -  list of real values
    
    
    
    """
    diffsx=dict()
    h=0
    if(S._verbose>1):
        print "r,Y1,Y2=",r,Y1,Y2
    if(ndigs>=15 and not isinstance(r,sage.libs.mpmath.ext_main.mpf)):
        raise TypeError,"Need mpmath input! got r=",r
    C1=Maassform_coeffs(S,r,M,ST,Yset=Y1,ndigs=ndigs)
    if(method=='TwoY'):
        C2=Maassform_coeffs(S,r,M,ST,Yset=Y2,ndigs=ndigs)
        for j in range(1 ,4 ):
            diffsx[j]=(C1[0 ][c[j]]-C2[0 ][c[j]]).real
            if(not first_time and signs.keys().count(j)>0 ):
                h=h+signs[j]*diffsx[j]
            else:
                h=h+abs(diffsx[j])
        return [diffsx,h]
    elif(method=='Hecke'):
        h=(C1[0 ][2 ]*C1[0 ][3 ]-C1[0 ][6 ]).real
        diffsx[1 ]=0 ;diffsx[1 ]=0 ;diffsx[3 ]=0 ;
        #print "c2c3-c6=",h
        return [diffsx,h]

def find_Y_and_M(G,R,ndigs=12,Yset=None,Mset=None):
    r"""
    Compute a good value of M and Y for Maass forms on G

    INPUT:

    - ''G'' -- group
    - ''R'' -- real
    - ''ndigs'' -- integer (number of desired digits of precision)
    - ''Yset'' -- real (default None) if set we return M corr. to this Y
    - ''Mset'' -- integer (default None) if set we return Y corr. to this M

    OUTPUT:

    - [Y,M] -- good values of Y (real) and M (integer)

    EXAMPLES::

    

    TODO:
    Better and more effective bound
    """
    l=G._level
    if(Mset <> None):
        # then we get Y corr. to this M
        Y0=mpmath.sqrt(3.0)/mpmath.mpf(2*l)
        
    if(Yset==None):
        Y0=mpmath.sqrt(3.0)/mpmath.mpf(2*l)
        Y=mpmath.mpf(0.95*Y0)
    else:
        Y=Yset
    #print "Y=",Y,"Yset=",Yset
    IR=mpmath.mpc(0,R)
    eps=mpmath.mpf(10 **-ndigs)
    twopiY=mpmath.pi()*Y*mpmath.mpf(2)

    M0=get_M_for_maass(R,Y,eps) 
    if(M0<10):
        M0=10
    ## Do this in low precision
    dold=mpmath.mp.dps
    #print "Start M=",M0
    #print "dold=",dold
    #mpmath.mp.dps=100
    try:
        for n in range(M0,10000 ):
            X=mpmath.pi()*Y*mpmath.mpf(2 *n)
            #print "X,IR=",X,IR
            test=mpmath.fp.besselk(IR,X)
            if(abs(test)<eps):
                raise StopIteration()
    except StopIteration:
        M=n
    else:
        M=n
        raise Exception,"Error: Did not get small enough error:=M=%s gave err=%s" % (M,test)
    mpmath.mp.dps=dold
    return [Y,M]







def testing_kbes(Rt,Xt):
    [R0,R1,NR]=Rt
    [X0,X1,NX]=Xt
    NRr=mpmath.mpf(NR)
    NXr=mpmath.mpf(NX)
    for j in range(1,NR):
        rj=mpmath.mpf(j)
        R=R0+R1*rj/NRr
        iR=mpmath.mpc(0,R)
        for k in range(1,NX):
            rk=mpmath.mpf(k)
            x=X0+X1*rk/NXr
            print "r,x=",R,x
            if(x>R):
                print "kbes_asymp="
                timeit( "kbes_asymp(R,x)",repeat=1)
            else:
                print "kbes_rec="
                timeit( "kbes_rec(R,x)",repeat=1)
            print "mpmath.besselk="
            timeit("mpmath.besselk(iR,x)",repeat=1)
            

            #print "t1(",R,x,")=",t1
            #print "t2(",R,x,")=",t2
            if(R<15.0):
                if(x<0.3 *R):
                    print "Case 1"
                elif(x<=max(10.0 +1.2*R,2 *R)):
                    print "Case 2"
            elif(R>20  and x>4 *R):
                print "Case 3"
            else:
                print "Case 4"


def test_Hecke_relations(a,b,C):
    r"""Testing Hecke relations for the Fourier coefficients in C



    INPUT:
    -''C'' -- dictionary of complex (Fourier coefficients)
    -''a'' -- integer
    -''b'' -- integer

    OUTPUT:
    -''diff'' -- real : |C(a)C(b)-C(ab)| if (a,b)=1

    EXAMPLE::

    
    sage: S=MaassWaveForms(Gamma0(1))
    sage: R=mpmath.mpf(9.53369526135355755434423523592877032382125639510725198237579046413534)
    sage: Y=mpmath.mpf(0.85)
    sage: C=coefficients_for_Maass_waveforms(S,R,Y,10,20,12)
    sage: d=test_Hecke_relations(C,2,3); mppr(d)
    '9.29e-8'
    sage: C=coefficients_for_Maass_waveforms(S,R,Y,30,50,20)
    sage: d=test_Hecke_relations(C,2,3); mppr(d)
    '3.83e-43'
    
    
    """
    c=gcd(Integer(a),Integer(b))
    lhs=C[0][a]*C[0][b]
    rhs=mpmath.mpf(0)
    for d in divisors(c):
        rhs=rhs+C[0][Integer(a*b/d/d)]
    return abs(rhs-lhs)


def test_Hecke_relations_all(C):
    r"""
    Test all possible Hecke relations.


    EXAMPLE::

    
        sage: S=MaassWaveForms(Gamma0(1))
        sage: mpmath.mp.dps=100
        sage: R=mpmath.mpf(9.53369526135355755434423523592877032382125639510725198237579046413534899129834778176925550997543536649304476785828585450706066844381418681978063450078510030977880577576)
        sage: Y=mpmath.mpf(0.85)
        sage: C=coefficients_for_Maass_waveforms(S,R,Y,30,50,20)
        sage: test=test_Hecke_relations_all(C); test  
        {4: '9.79e-68', 6: '4.11e-63', 9: '4.210e-56', 10: '9.47e-54', 14: '2.110e-44', 15: '4.79e-42', 21: '4.78e-28', 22: '1.02e-25', 25: '9.72e-19', 26: '2.06e-16'}
        

    We can see how the base precision used affects the coefficients

        sage: mpmath.mp.dps=50
        sage: C=coefficients_for_Maass_waveforms(S,R,Y,30,50,20)
        sage: test=test_Hecke_relations_all(C); test               
        sage: test=test_Hecke_relations_all(C); test  
{4: '1.83e-48', 6: '4.75e-43', 9: '3.21e-36', 10: '1.24e-33', 14: '4.41e-25', 15: '1.53e-23', 21: '6.41e-8', 22: '4.14e-6', 25: '91.455', 26: '12591.0'}

    

    """
    N=max(C[0].keys())
    test=dict()
    for a in prime_range(N):
        for b in prime_range(N):
            if(a*b <= N):
                test[a*b]=mppr(test_Hecke_relations(C,a,b))
    return test

def set_norm_maass(k,cuspidal=True):
    r""" Set normalization for computing maass forms.

    INPUT:

     - ``k`` -- dimenson
     - ``cuspidal`` -- cuspidal maass waveforms (default=True)


     OUTPUT:

     - ``N`` -- normalization (dictionary)
     -- N['comp_dim'] -- dimension of space we are computing in 
     -- N['SetCs']``     -- which coefficients are set
     -- N['Vals']     -- values of set coefficients     

     EXAMPLES::

         sage: set_norm_maass(1)
         {'Vals': {0: {0: 0, 1: 1}}, 'comp_dim': 1, 'num_set': 2, 'SetCs': [0, 1]}
         sage: set_norm_maass(1,cuspidal=False)
         {'Vals': {0: {0: 1}}, 'comp_dim': 1, 'num_set': 1, 'SetCs': [0]}
         sage: set_norm_maass(2)
         {'Vals': {0: {0: 0, 1: 1, 2: 0}, 1: {0: 0, 1: 0, 2: 1}}, 'comp_dim': 2, 'num_set': 3, 'SetCs': [0, 1, 2]}
         
    """
    C=dict()
    Vals=dict()
    #  set coeffs c(0),c(1),...,c(k-1) if not cuspidal
    #  set coeffs c(0)=0,c(1),...,c(k) if cuspidal 
    if(cuspidal and k>0):
        SetCs=range(0,k+1)
    else:
        SetCs=range(0,k)
    #if(cuspidal):  # have to set other cusps too
    #    for i  in range(1,len(G.cusps())+1):
    #       SetCs.append(0+Ml*i)
    if(cuspidal):
        C['cuspidal']=True
    else:
        C['cuspidal']=False
    for j in range(k):
        Vals[j]=dict()
        for n in SetCs:
            Vals[j][n]=0
        if(cuspidal):
            Vals[j][j+1]=1
        else:
            Vals[j][j]=1
    C['comp_dim']=k
    C['SetCs']=SetCs
    C['Vals']=Vals
    return C

def solve_system_for_Maass_waveforms(W,N,deb=False):
    r"""
    Solve the linear system to obtain the Fourier coefficients of Maass forms

    INPUT:

    - ``W`` --   (system) dictionary

        - ``W['Ms']``  -- M start
        - ``W['Mf']``  -- M stop
        - ``W['nc']``  -- number of cusps
        - ``W['V']``   -- matrix of size ((Ms-Mf+1)*nc)**2
        - ``W['RHS']`` -- right hand side (for inhomogeneous system) matrix of size ((Ms-Mf+1)*nc)*(dim)

    - ``N`` -- normalisation (dictionary, output from the set_norm_for_maass function)

        - ``N['SetCs']``   -- Which coefficients are set
        - ``N['Vals'] ``   -- To which values are these coefficients set
        - ``N['comp_dim']``-- How large is the assumed dimension of the solution space
        
    - ``deb`` -- print debugging information (default False)

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
    import mpmath.mp,mpmath.fp
    V=W['V']
    Ms=W['Ms']
    Mf=W['Mf']
    nc=W['nc']
    Ml=Mf-Ms+1
    if(V.cols<>Ml*nc or V.rows<>Ml*nc):
        raise Exception," Wrong dimension of input matrix!"
    SetCs=N['SetCs']
    Vals=N['Vals']
    comp_dim=N['comp_dim']
    if(N['cuspidal']):
        for i in range(1,nc):
            if(SetCs.count(i*Ml)==0):
                SetCs.append(i*Ml)
            for fn_j in range(comp_dim):
                Vals[fn_j][i*Ml]=0

    if(Ms<0):
        use_sym=0
    else:
        use_sym=1
    if(use_sym==1 and SetCs.count(0)>0):
        num_set=len(N['SetCs'])-1
    else:
        num_set=len(N['SetCs'])
    t=V[0,0]
    if(isinstance(t,float)):
        mpmath_ctx=mpmath.fp
    else:  
        mpmath_ctx=mpmath.mp
    if(W.has_key('RHS')):
        RHS=W['RHS']
    else:
        RHS=mpmath_ctx.matrix(int(Ml*nc-num_set),int(comp_dim))
    LHS=mpmath_ctx.matrix(int(Ml*nc-num_set),int(Ml*nc-num_set))
    roffs=0

    if(deb):
        print "num_set,use_sym=",num_set,use_sym
        print "SetCs,Vals=",SetCs,Vals
        print "V.rows,cols=",V.rows,V.cols
        print "LHS.rows,cols=",LHS.rows,LHS.cols
        print "RHS.rows,cols=",RHS.rows,RHS.cols
    for r in range(V.rows):
        cr=r+Ms
        if(SetCs.count(r+Ms)>0):
            roffs=roffs+1
            continue
        for fn_j in range(comp_dim):
            RHS[r-roffs,fn_j]=mpmath_ctx.mpf(0)
            for cset in SetCs:
                v=Vals[fn_j][cset]
                if(mpmath_ctx==mpmath.mp):
                    tmp=mpmath_ctx.mpmathify(v)
                elif(isinstance(v,float)):
                    tmp=mpmath_ctx.mpf(v)
                else:
                    tmp=mpmath_ctx.mpc(v)
                #print "tmp=",tmp
                #print "V[",r,cset-Ms,"]=",V[r,cset-Ms]
                tmp=tmp*V[r,cset-Ms]
                RHS[r-roffs,fn_j]=RHS[r-roffs,fn_j]-tmp
        coffs=0
        for k in range(V.cols):            
            if(SetCs.count(k+Ms)>0):
                coffs=coffs+1
                continue
            # print "r,k=",r,k
            #print "roffs,coffs=",roffs,coffs
            #print "r-roffs,k-coffs=",r-roffs,k-coffs
            LHS[r-roffs,k-coffs]=V[r,k]
            #print "LHS[",r,k,"]=",LHS[r-roffs,k-coffs]
    #print "RHS="
    #for j in range(RHS.rows):
    #    print j,RHS[j,0]#
    try:
        A, p = mpmath_ctx.LU_decomp(LHS)
    except ZeroDivisionError:
        t1=smallest_inf_norm(LHS)
        print "n=",smallest_inf_norm(LHS)
        t2=mpmath_ctx.log10(smallest_inf_norm(LHS))
        t3=mpmath_ctx.ceil(-t2)
        isinf=False
        if(isinstance(t3,float)):
            isinf = (t3 == float(infinity))
        if(isinstance(t3,sage.libs.mpmath.ext_main.mpf)):
            isinf = ((t3.ae(mpmath.inf)) or t3==mpmath.inf)
        if(isinstance(t3,sage.rings.real_mpfr.RealLiteral)):
            isinf = t3.is_infinity()
        if(isinf):
            raise ValueError, " element in LHS is infinity! t3=%s" %t3
        #print "LHS="
        #for j in range(LHS.rows):
        #    print j,LHS.column(j)
        #print "t3=",t3
        t=int(t3)
        #t=int(mpmath_ctx.ceil(-mpmath_ctx.log10(smallest_inf_norm(LHS))))
        #raise ZeroDivisionError,"Need higher precision! Use > %s digits!" % t
        return t
    X=dict()
    for fn_j in range(comp_dim):
        X[fn_j] = dict() #mpmath.matrix(int(Ml),int(1))
        b = mpmath_ctx.L_solve(A, RHS.column(fn_j), p)
        TMP = mpmath_ctx.U_solve(A, b)
        roffs=0
        res = mpmath_ctx.norm(mpmath_ctx.residual(LHS, TMP, RHS.column(fn_j)))
        #print "res(",fn_j,")=",res
        for i in range(nc):
            X[fn_j][i]=dict()
        for n in range(Ml):
            if(SetCs.count(n+Ms)>0):
                roffs=roffs+1
                #print "X[",fn_j,",",n,",Vals[fn_j][n]
                X[fn_j][0][n+Ms]=Vals[fn_j][n+Ms]
                continue
            X[fn_j][0][n+Ms]=TMP[n-roffs,0]
            #print "X[",fn_j,",",n+Ms,"=",TMP[n-roffs,0]
        for i in range(1,nc):
            for n in range(Ml):
                if(SetCs.count(n+Ms+i*Ml)>0):
                    roffs=roffs+1
                    # print "X[",fn_j,",",n,",Vals[fn_j][n]
                    X[fn_j][i][n+Ms]=Vals[fn_j][n+Ms+i*Ml]
                    continue
                X[fn_j][i][n+Ms]=TMP[n+i*Ml-roffs,0]
    # return x
    return X




# Stuff needed for attaching the file locally
#__main__.MaassWaveForms=MaassWaveForms
#__main__.MaassWaveformElement=MaassWaveformElement
