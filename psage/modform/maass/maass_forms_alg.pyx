# -*- coding=utf-8 -*-
#*****************************************************************************
#  Copyright (C) 2010 Fredrik Strömberg <stroemberg@mathematik.tu-darmstadt.de>
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

### Cython includes 
include 'stdsage.pxi'
include 'cdefs.pxi'
include 'interrupt.pxi'
import sage.structure.element
cimport sage.structure.element
from sage.structure.element cimport Element, ModuleElement, RingElement
from sage.modular.cusps import Cusp
from sage.rings.infinity import infinity
from sage.rings.integer import Integer,is_Integer
from sage.rings.complex_double import CDF
from numpy import array
import numpy as np
cimport numpy as cnp
#cimport cython
import cython
cdef extern from "math.h":
    double fabs(double)
    double fmax(double,double)
    int ceil(double) 

## Try to use mpmath for high precision computations
# MPMATH includes
import sage.libs.mpmath.all as mpmath
#from sage.libs.mpmath.ext_main import *

from sage.libs.mpfr cimport *


from sage.functions.all import ceil as pceil
from sage.modular.arithgroup.congroup_sl2z import SL2Z
from mysubgroup import MySubgroup
from lpkbessel import besselk_dp

mpftype=sage.libs.mpmath.ext_main.mpf
mp0=mpmath.mpf(0)
mp5=mpmath.mpf(5)
mp8=mpmath.mpf(8)
mp4=mpmath.mpf(4)
mp1=mpmath.mpf(1)
mp2=mpmath.mpf(2)
mp24=mpmath.mpf(24)
mp36=mpmath.mpf(36)



def pullback_to_G(G,x,y,mp_ctx=mpmath.mp):
    r""" Mpmath version of general pullback alg

    INPUT:
    
     - ``x`` -- real
     - ``y`` -- real > 0
     - ``G`` -- MySubgroup 

    EXAMPLES::

    
        sage: G=MySubgroup(SL2Z)
        sage: [x,y,A]=pullback_to_G_mp(G,mpmath.mpf(0.1),mpmath.mpf(0.1))
        sage: x,y
        (mpf('6.9388939039072274e-16'), mpf('4.9999999999999991'))
        sage: A
        [ 5 -1]
        [ 1  0]
        sage: G=MySubgroup(Gamma0(3))
        sage: [x,y,A]=pullback_to_G_mp(G,mpmath.mpf(0.1),mpmath.mpf(0.1))
        sage: x,y
        (mpf('-0.038461538461538491'), mpf('0.19230769230769232'))
        sage: A
        [-1  0]
        [ 6 -1]

    

    """
    try:
        reps=G._coset_reps
    except AttributeError:
        reps=list(G.coset_reps())

    if(mp_ctx==mpmath.mp):
        try:
            x.ae; y.ae
        except AttributeError:
            raise Exception,"Need x,y in mpmath format!"
        [x1,y1,[a,b,c,d]]=pullback_to_psl2z_mp(x,y)
    elif(mp_ctx==mpmath.fp):
        [x1,y1,[a,b,c,d]]=pullback_to_psl2z_dble(x,y)
    else:
        raise NotImplementedError," Only implemented for mpmath.mp and mpmath.fp!"
    
    A=SL2Z([a,b,c,d])
    try:
        for V in reps:
            B=V*A
            if(B in G):
                raise StopIteration
    except StopIteration:            
        pass
    else:
        raise Exception,"Did not find coset rep. for A=%s" % A
    [xpb,ypb]=apply_sl2_map(x,y,B,mp_ctx=mp_ctx)
    return [xpb,ypb,B]
    



        
def apply_gl2_map(x,y,A):
    r"""
    Apply the matrix A in SL(2,R) to the point z=xin+I*yin 

    INPUT:

     - ``x`` --  real number
     - ``y`` --  real number >0
     - ``A`` --  2x2 Real matrix with determinant one

     OUTPUT:

     - [u,v] -- pair of real points u,v>0

     EXAMPLES::


         sage: A=matrix(RR,[[1,2],[3,4]])
         sage: apply_gl2_map(mpmath.mpf(0.1),mpmath.mpf(0.1),A)
         [mpf('0.48762109795479019'), mpf('-0.010764262648008612')]

    
    """
    try:
        x.ae; y.ae
    except AttributeError:
        raise Exception,"Need x,y in mpmath format!"
    a=A[0,0]; b=A[0,1]; c=A[1,0]; d=A[1,1]
    aa=mpmath.mpf(a);bb=mpmath.mpf(b);cc=mpmath.mpf(c);dd=mpmath.mpf(d)
    tmp1=cc*x+dd
    tmp2=cc*y
    den=tmp1*tmp1+tmp2*tmp2
    u=(aa*cc*(x*x+y*y)+bb*dd+x*(aa*dd+bb*cc))/den
    det=aa*dd-bb*cc
    v=y*det/den
    return [u,v]

def apply_sl2_map(x,y,A,mp_ctx=mpmath.mp,inv=False):
    r"""
    Apply the matrix A in SL(2,R) to the point z=xin+I*yin 

    INPUT:
    
    - ``x`` --  real number
    - ``y`` --  real number >0
    - ``A`` -- 2x2 Real matrix with determinant one
    - ``inv`` -- (default False) logical
           = True  => apply the inverse of A
           = False => apply A

   OUTPUT:

     - [u,v] -- pair of real points u,v>0

    EXAMPLES::


        sage: apply_sl2_map(mpmath.mpf(0.1),mpmath.mpf(0.1),A)
        [mpf('0.69893607327370766'), mpf('1.3806471921778053e-5')]


    

    """
   
    if(not (test_ctx(x,mp_ctx) or test_ctx(y,mp_ctx))):
        raise TypeError,"Need x,y of type : %s" %(type(mp_ctx.mpf(1)))
    if(inv):
        d=A[0,0]; b=-A[0,1]; c=-A[1,0]; a=A[1,1]
    else:
        a=A[0,0]; b=A[0,1]; c=A[1,0]; d=A[1,1]
        
    aa=mp_ctx.mpf(a);bb=mp_ctx.mpf(b);cc=mp_ctx.mpf(c);dd=mp_ctx.mpf(d)
    tmp1=cc*x+dd
    tmp2=cc*y
    den=tmp1*tmp1+tmp2*tmp2
    u=(aa*cc*(x*x+y*y)+bb*dd+x*(aa*dd+bb*cc))/den
    v=y/den
    return [u,v]

def test_ctx(x,mp_ctx):
    r""" Check that x is of type mp_ctx

    INPUT:

        - ``x`` -- real
        - ``mp_ctx`` -- mpmath context

    OUTPUT:
    
        - ``test`` -- True if x has type mp_ctx otherwise False 

    EXAMPLES::


        sage: test_ctx(mpmath.mp.mpf(0.5),mpmath.fp)
        False
        sage: test_ctx(mpmath.mp.mpf(0.5),mpmath.mp)
        True
    

    """
    if(not ( isinstance(mp_ctx,mpmath.ctx_mp.MPContext) or isinstance(mp_ctx,mpmath.ctx_fp.FPContext))):
        raise TypeError," Need mp_ctx to be a mpmath mp or fp context. Got: mp_ctx of type %s" % type(mp_ctx)
    t=type(mp_ctx.mpf(0.1))
    if(isinstance(x,t)):
        return True
    else:
        return False

    
def normalize_point_to_cusp(G,x,y,cu,mp_ctx=mpmath.mp,inv=False):
    r"""
    Compute the normalized point with respect to the cusp cu

    INPUT:

    - ``G``  -- MySubgroup
    - ``x``  -- real
    - ``y``  -- real > 0
    - ``cu`` -- Cusp of G
    - ``inv``-- (default False) logical
              - True  => apply the inverse of the normalizer
              - False => apply the normalizer

    OUTPUT:

    - [u,v] -- pair of reals, v > 0


    EXAMPLES::


        sage: G=MySubgroup(Gamma0(3))
        sage: x=mpmath.mpf(0.1); y=mpmath.mpf(0.1)
        sage: normalize_point_to_cusp(G,x,y,Cusp(infinity))
        [mpf('0.10000000000000001'), mpf('0.10000000000000001')]
        sage: normalize_point_to_cusp(G,x,y,Cusp(0))
        [mpf('-1.6666666666666665'), mpf('1.6666666666666665')]
    
    
    """
    if(mp_ctx==mpmath.mp):
        try:
            x.ae; y.ae
        except AttributeError:
            raise Exception,"Need x,y in mpmath format!"
    if(cu==Cusp(infinity) and G.cusp_width(cu)==1): # Inf is the first cusp
        return [x,y]
    wi=mp_ctx.mpf(G.cusp_width(cu))
    if(inv):
        [d,b,c,a]=G.cusp_normalizer(cu)
        b=-b; c=-c
    else:
        [a,b,c,d]=G.cusp_normalizer(cu)
        x=x*wi; y=y*wi
    aa=mp_ctx.mpf(a);bb=mp_ctx.mpf(b);cc=mp_ctx.mpf(c);dd=mp_ctx.mpf(d)
    tmp1=cc*x+dd
    tmp2=cc*y
    den=tmp1*tmp1+tmp2*tmp2
    u=(aa*cc*(x*x+y*y)+bb*dd+x*(aa*dd+bb*cc))/den
    v=y/den
    if(inv):
        u=u/wi
        v=v/wi
    return [u,v]

def j_fak(c,d,x,y,k,unitary=True):
    r"""
    Computes the argument (with principal branch)  of the
    automorphy factor j_A(z;k) defined by:
    For A=(a b ; c d ) we have j_A(z)=(cz+d)^k=|cz+d|^k*exp(i*k*Arg(cz+d))

    INPUT:

    - ``c`` -- integer or real
    - ``d`` -- integer or real
    - ``x`` -- real (mpmath.mpf)
    - ``y`` -- real (mpmath.mpf) > 0
    - ``k`` -- real (mpmath.mpf)
    - ``unitary`` -- (default True) logical
               = True  => use the unitary version (i.e. for Maass forms)
               = False => use the non-unitary version (i.e. for holomorphic forms)

    OUTPUT:

    - ``t``     -- real: t = k*Arg(cz+f) (if unitary==True)
    - ``[r,t]`` -- [real,real]: r=|cz+d|^k, t=k*Arg(cz+f) (if unitary==False)
    EXAMPLES::

        sage: j_fak(3,2,mpmath.mpf(0.1),mpmath.mpf(0.1),mpmath.mpf(4),True)
        mpf('0.51881014862364842780456450654158898199306583839137025')
        sage: [u,v]=j_fak(3,2,mpmath.mpf(0.1),mpmath.mpf(0.1),mpmath.mpf(4),False)
        sage: u
        mpf('28.944399999999999999999999999999999999999999999999763')
        sage: v
        mpf('0.51881014862364842780456450654158898199306583839137025')


    

    """
    if(not isinstance(c,sage.libs.mpmath.ext_main.mpf)):
        cr=mpmath.mpf(c); dr=mpmath.mpf(d)        
    else:
        cr=c; dr=d
    try:
        x.ae; y.ae; k.ae
    except AttributeError:
        print "x,y,k=",x,y,k
        raise TypeError,"Need x,y,c,d,k in mpmath format!"
    if( (cr.ae(0) and dr.ae(1)) or (k==0) or (k.ae(mpmath.mpf(0)))):
        if(not unitary):
            return (mpmath.mp.mpf(1),mpmath.mp.mpf(0))
        else:
            return mpmath.mp.mpf(0)
    tmp=mpmath.mpc(cr*x+dr,cr*y)
    tmparg=mpmath.arg(tmp)
    # Consider principal branch of the argument
    tmparg=tmparg*mpmath.mpf(k)
    #while tmparg <= -mpmath.pi():
    #    tmparg=tmparg+mptwopi
    #while tmparg > mpmath.pi():
    #    tmparg=tmparg-mptwopi
    if(not unitary):
        tmpabs=mpmath.power(mpmath.absmax(tmp),k)
        return (tmpabs,tmparg)
    else:
        return tmparg

    
def pullback_pts(G,Qs,Qf,Y,weight=None,holo=False,deb=False,mp_ctx=mpmath.mp):
    r""" Computes a whole array of pullbacked points-

    INPUT:

    - ``G``  -- MySubgroup
    - ``Qs`` -- integer
    - ``Qf`` -- integer
    - ``Y``  -- real (mpmath) 
    - ``weight`` -- (optional) real
    - ``holo``   -- (False) logical
    - ``deb``    -- (False) logical
    -``mp_x``  -- Context for mpmath (default mp)
    
    OUTPUT:

    - ``pb`` -- dictonary with entries:
       - 'xm'   -- real[Qf-Qs+1]  : x_m=2*pi*(1-2*m)/2(Qf-Qs+1)
       - 'xpb'  -- real[0:nc,0:nc,Qf-Qs+1]  : real part of pullback of x_m+iY
       - 'ypb'  -- real[0:nc,0:nc,Qf-Qs+1]  : imag. part of pullback of x_m+iY
       - 'cvec' -- complex[0:nc,0:nc,Qf-Qs+1] : mult. factor

    EXAMPLES::


        sage: G=MySubgroup(Gamma0(2))
        sage: [Xm,Xpb,Ypb,Cv]=pullback_pts(G,1,3,mpmath.mpf(0.1))   

    Note that we need to have $Y<Y_{0}$ so that all pullbacked points
    are higher up than Y.

        sage: pullback_pts(G,1,10,mpmath.mpf(0.44))
        Traceback (most recent call last):
        ...
        ArithmeticError: Need smaller value of Y

    
    """
    if(not isinstance(G,sage.modular.maass.mysubgroup.MySubgroup)):
        GG=MySubgroup(G)
    else:
        GG=G
    mpmath_ctx=mp_ctx
    #mpmath.mp.dps=500
    if(mpmath_ctx == mpmath.fp):
        pi=mpmath.pi
    elif(mpmath_ctx == mpmath.mp):
        pi=mpmath.pi()
        try:
            Y.ae
        except AttributeError:
            raise TypeError,"Need Y in mpmath format for pullback!"
    else:
        raise NotImplementedError," Only implemented for mpmath.mp and mpmath.fp!"
    
    twopi=mpmath_ctx.mpf(2)*pi
    twopii=mpmath_ctx.mpc(0,twopi)
    mp_i=mpmath_ctx.mpc(0,1)
    Xm=dict()  #vector(RF,2*Q)
    Xpb=dict()  #vector(RF,2*Q)
    Ypb=dict()  #vector(RF,2*Q)
    Cvec=dict()
    if(holo and weight==None):
        raise Exception,"Must support a weight for holomorphic forms!"
    elif(holo):
        try:
            weight.ae
        except AttributeError:
            raise TypeError,"Need weight in mpmath format for pullback!"
    #Qs=1-Q; Qf=Q
    if(deb):
        fp0=open("xm.txt","w")
        fp1=open("xpb.txt","w")
        fp2=open("ypb.txt","w")
        fp3=open("cv.txt","w")
    twoQl=mpmath_ctx.mpf(2*(Qf+1-Qs))
    if(Qs<0):
        for j in range(Qs,Qf+1): #1-Q,Q):
            Xm[j]=mpmath_ctx.mpf(2*j-1)/twoQl        
            if(deb):
                s=str(Xm[j])+"\n"
                fp0.write(s)
    else:
        for j in range(Qs,Qf+1): #1-Q,Q):
            Xm[j]=mpmath_ctx.mpf(2*j-1)/twoQl        
            #print "Xm[",j,"]=",Xm[j]
            if(deb):
                s=str(Xm[j])+"\n"
                fp0.write(s)
    for ci in GG.cusps():
        cii=GG.cusps().index(ci)
        if(deb):
            print "cusp =",ci
        swi=mpmath_ctx.sqrt(mpmath_ctx.mpf(GG.cusp_width(ci)))
        for j in range(Qs,Qf+1): #1-Q,Q):
            [x,y]   = normalize_point_to_cusp(GG,Xm[j],Y,ci,mp_ctx=mpmath_ctx)
            [x1,y1,Tj] =  pullback_to_G(GG,x,y,mp_ctx=mpmath_ctx)
            v = GG.closest_vertex(x1,y1)
            cj= GG._cusp_representative[v]
            cjj=GG._cusps.index(cj)
            swj=mpmath_ctx.sqrt(mpmath_ctx.mpf(GG._cusp_width[cj][0]))
            U = GG._vertex_map[v]
            if(U<>SL2Z.gens()[0]**4):
                [x2,y2] = apply_sl2_map(x1,y1,U)
            else:
                x2=x1; y2=y1;
            [x3,y3] = normalize_point_to_cusp(GG,x2,y2,cj,inv=True,mp_ctx=mpmath_ctx)
            #Xpb[cii,cjj,j]=mpmath_ctx.mpc(0,x3*twopi)
            Xpb[cii,cjj,j]=x3*twopi
            # Recall that Ypb must be greater than Y otherwise we need a better Y
            if(y3>Y):
                Ypb[cii,cjj,j]=y3
            else:
                raise ArithmeticError,"Need smaller value of Y" 
            if(deb):
                print "ci,cj=",ci,cj
                print "Xm=",Xm[j]
                print "x,y=",x,"\n",y
                print "x1,y1=",x1,"\n",y1
                print "x2,y2=",x2,"\n",y2
                print "x3,y3=",x3,"\n",y3
                print "Xpb=",Xpb[cii,cjj,j]
                print "Ypb=",Ypb[cii,cjj,j]
            # We also get the multiplier if we need it
            if(holo):
                c = mpmath_ctx.mpf(GG._cusp_normalizer[ci][1,0])*swi
                d = mpmath_ctx.mpf(GG._cusp_normalizer[ci][1,1])/swi
                m1=j_fak(c,d,Xm[j],Y,-weight)
                c = mpmath_ctx.mpf(GG._cusp_normalizer[cj][1,0])*swj
                d = mpmath_ctx.mpf(GG._cusp_normalizer[cj][1,1])/swj
                m2=j_fak(c,d,x3,y3,weight)
                A=(U*Tj)
                #A=A**-1
                c=mpmath_ctx.mpf(-A[1,0]); d=mpmath_ctx.mpf(A[0,0])
                m3=j_fak(c,d,x2,y2,weight)
                tmp=m1+m2+m3
                Cvec[cii,cjj,j]=tmp # mpmath_ctx.mpc(0,tmp)
            if(deb):
                fp1.write(str(x3)+"\n")
                fp2.write(str(y3)+"\n")
                fp3.write(str(tmp)+"\n")

    for j in range(Qs,Qf+1): #1-Q,Q):
        #        Xmi[j]=mpmath_ctx.mpc(0,Xm[j]*twopi)
        Xm[j]=Xm[j]*twopi
    if(deb):
        fp1.close()
        fp2.close()
        fp3.close()
    pb=dict()
    pb['xm']=Xm; pb['xpb']=Xpb; pb['ypb']=Ypb; pb['cvec']=Cvec
    return pb
    #return [Xm,Xpb,Ypb,Cvec]

def setup_matrix_for_Maass_waveforms(G,R,Y,int M,int Q,cuspidal=True,sym_type=None,low_prec=False):
    r"""
    Set up the matrix for the system of equations giving the Fourier coefficients of the Maass waveforms.

    INPUT:
    
    - ``G`` -- MySubgroup
    - ``R`` -- real (mpmath)
    - ``Y`` -- real (mpmath)
    - ``M`` -- integer
    - ``Q`` -- integer > M
    - ``cuspidal`` -- (False) logical : Assume cusp forms
    - ``sym_type``  -- (None) integer

          -  0 -- even / sine
          -  1 -- odd  / cosine

    - ``low_prec`` -- (False) : Use low (double) precision K-Bessel

    OUTPUT:
    
    - ``W`` -- dictionary
    - ``W['Mf']`` -- M start
    - ``W['nc']`` -- number of cusps
    - ``W['Ms']`` -- M stop
    - ``W['V']``  -- ((Ms-Mf+1)*num_cusps)**2 real number
    
    
    EXAMPLES::

        sage: R=mpmath.mpf(9.533695261353557554344235235928770323821256395107251982375790464135348991298347781769255509975435)
        sage: Y=mpmath.mpf(0.1)
        sage: G=MySubgroup(Gamma0(1))
        sage: W=setup_matrix_for_Maass_waveforms(G,R,Y,10,20)
        sage: W.keys()
        ['Mf', 'nc', 'Ms', 'V']
        sage: print W['Ms'],W['Mf'],W['nc'];
        -10 10 1
        sage: V.rows,V.cols
        (21, 21)
        sage: N=set_norm_maass(1)
        sage: C=solve_system_for_Maass_waveforms(W['V'],N)
        sage: c2.real
        mpf('-1.0683335512235690597444640854010889121099661684073309')
        sage: C[0][2].real
        mpf('-1.0683335512235690597444640854010889121099661684073309')
        sage: C[0][3].real
        mpf('-0.4561973545061180270983019925431624077983243678097305')
        sage: C[0][6].real
        mpf('0.48737093979829448482994037031228389306060715149700209')
        sage: abs(C[0][6]-C[0][2]*C[0][3])
        mpf('0.00000000000002405188994758679782568437503269848800837')
        
    

    """
    #print "sym_type,R=",sym_type,R
    cdef int l,j,icusp,jcusp,n,ni,li,iMl,iMs,iMf,iQs,iQf,iQl,s,nc
    cdef double RR,YY
    if(low_prec==True):
        RR=<double>R
        YY=<double>Y
        if(sym_type<>None):
            W=setup_matrix_for_Maass_waveforms_np_dble(G,RR,YY,int(M),int(Q),int(cuspidal),int(sym_type))
        elif(sym_type==None):
            W=setup_matrix_for_Maass_waveforms_np_cplx(G,RR,YY,int(M),int(Q),int(cuspidal))
        return W
    if(Q < M):
        raise ValueError," Need integers Q>M. Got M=%s, Q=%s" %(M,Q)
    mpmath_ctx=mpmath.mp
    if( not test_ctx(R,mpmath_ctx) or not test_ctx(Y,mpmath_ctx)):
        raise TypeError," Need Mpmath reals as input!"
    if(not isinstance(G,sage.modular.maass.mysubgroup.MySubgroup)):
        GG=MySubgroup(G)
    else:
        GG=G
    pi=mpmath.mp.pi()
    mp2=mpmath_ctx.mpf(2)
    twopi=mp2*pi
    IR=mpmath_ctx.mpc(0,R)
    nc=int(GG.ncusps())
    if(Q<M):
        Q=M+20
    if(sym_type<>None): 
        Ms=1; Mf=M; Ml=Mf-Ms+1
        Qs=1; Qf=Q;
        Qfak=mpmath_ctx.mpf(Q)/mp2
    else:
        Ms=-M; Mf=M; Ml=Mf-Ms+1
        Qs=1-Q; Qf=Q
        Qfak=mpmath_ctx.mpf(2*Q)
    iMl=int(Ml); iMs=int(Ms); iMf=int(Mf)
    iQf=int(Qf); iQs=int(Qs); iQl=int(Qf-Qs+1)
    # Pullpack points
    pb=pullback_pts(GG,Qs,Qf,Y,mp_ctx=mpmath_ctx)
    #[Xm,Xpb,Ypb,Cv]=
    s=int(nc*iMl)
    if(sym_type<>None):
        V=mpmath_ctx.matrix(int(s),int(s),force_type=mpmath_ctx.mpf)
    else:
        V=mpmath_ctx.matrix(int(s),int(s),force_type=mpmath_ctx.mpc)
    besarg_old=-1
    nvec=list()
    #print "keys=",Xpb.keys()
    for l in range(iMs,iMf+1):
        nvec.append(mpmath_ctx.mpf(l))
    if(sym_type==1):
        fun=mpmath_ctx.sin
    elif(sym_type==0):
        fun=mpmath_ctx.cos
    ef1=dict(); ef2=dict()
    for n from iMs <=n <= iMf:        
        nr=nvec[n-iMs]
        inr=mpmath_ctx.mpc(0,nr)
        for j from iQs <= j <= iQf:
            if(sym_type<>None):
                argm=nr*pb['xm'][j]
                ef2[n,j]=fun(argm)
            else:
                iargm=inr*pb['xm'][j]
                ef2[n,j]=mpmath_ctx.exp(-iargm)
            for jcusp in range(nc):
                for icusp in range(nc):
                    if(not pb['xpb'].has_key((icusp,jcusp,j))):
                        continue
                    if(sym_type<>None):
                        argpb=nr*pb['xpb'][icusp,jcusp,j]
                        ef1[j,icusp,jcusp,n]=fun(argpb)
                    else:
                        iargpb=inr*pb['xpb'][icusp,jcusp,j]
                        #if(n==1):
                        #    print "Xpb[",icusp,jcusp,j,"]=",Xpb[icusp,jcusp,j]
                        #    print "iargpb=",iargpb
                        ef1[j,icusp,jcusp,n]=mpmath_ctx.exp(iargpb)
                        
    fak=mpmath_ctx.exp(pi*mpmath_ctx.mpf(R)/mp2)
    for l from  iMs <= l <= iMf:
        if(cuspidal and l==0):
            continue
        lr=nvec[l-iMs]*twopi
        for j from iQs <= j <= iQf:
            for jcusp from int(0) <= jcusp < int(nc):
                lj=iMl*jcusp+l-iMs
                for icusp from int(0) <= icusp < int(nc):
                    if(not pb['ypb'].has_key((icusp,jcusp,j))):
                        continue
                    ypb=pb['ypb'][icusp,jcusp,j]
                    if(ypb==0):
                        continue
                    besarg=abs(lr)*ypb
                    #kbes=mpmath_ctx.sqrt(ypb)*(mpmath_ctx.besselk(IR,besarg).real)*fak
                    kbes=mpmath_ctx.sqrt(ypb)*(mpmath_ctx.besselk(IR,besarg).real)*fak
                    ckbes=kbes*ef1[j,icusp,jcusp,l]
                    for n from iMs <= n <= iMf:
                        if(cuspidal and n==0):
                            continue
                        ni=iMl*icusp+n-iMs
                        tmpV=ckbes*ef2[n,j]
                        V[ni,lj]=V[ni,lj]+tmpV
##                         if(ni==16 and lj==16):
##                             print "besarg(",j,")=",besarg
##                             print "ef2(",j,")=",ef2[n,j]
##                             print "ef1(",j,")=",ef1[j,icusp,jcusp,l]
##                             print "kbes(",j,")=",kbes
##                             print "ckbes(",j,")=",ckbes
##                             print "tmpV(",j,")=",tmpV
##                             print "V[16,16](",j,")=",V[ni,lj]


#    print "V[16,16]=",V[16,16]
    for n from 0<=n< V.rows:
        for l from 0 <= l < V.cols:        
            V[n,l]=V[n,l]/Qfak

#    print "V[16,16]=",V[16,16]
    #fak=fak*mpmath_ctx.sqrt(Y)
    twopiY=Y*twopi
#    print "sqrtY=",mpmath_ctx.sqrt(Y)
    for n from iMs<=n <=iMf:
        if(n==0 and cuspidal):
            continue
        nr=mpmath_ctx.mpf(abs(n))
        nrY2pi=nr*twopiY
        for icusp from int(0)<=icusp<nc:
            ni=iMl*icusp+n-iMs
            kbes=(mpmath_ctx.besselk(IR,nrY2pi).real)*fak
            #if(ni==16):
            #    print "arg=",nrY2pi
            #    print "kbes(",ni,")=",kbes
            V[ni,ni]=V[ni,ni]-kbes*mpmath_ctx.sqrt(Y)
    #print "V[16,16]=",V[16,16]
    W=dict()
    W['V']=V
    W['Ms']=Ms
    W['Mf']=Mf
    W['nc']=nc
    return W




DTYPE = np.float
ctypedef cnp.float_t DTYPE_t
CTYPE = np.complex128
ctypedef cnp.complex128_t CTYPE_t

@cython.boundscheck(False)
def setup_matrix_for_Maass_waveforms_np_dble(G,double R,double Y,int M,int Q,int cuspidal=-1, int sym_type=0):
    r"""
    Set up the matrix for the system of equations giving the Fourier coefficients of the Maass waveforms.

    INPUT:
    
    - ``G`` -- Subgroup of SL2Z (MySubgroup)
    - ``R`` -- real (mpmath)
    - ``Y`` -- real (mpmath)
    - ``M`` -- integer
    - ``Q`` -- integer > M
    - ``sym_type``  -- (None) integer

          -  0 -- even / sine
          -  1 -- odd  / cosine

    - ``cuspidal`` -- (False) logical : Assume cusp forms

    OUTPUT:
    
    - ``W`` -- dictionary
    - ``W['Mf']`` -- M start
    - ``W['nc']`` -- number of cusps
    - ``W['Ms']`` -- M stop
    - ``W['V']``  -- ((Ms-Mf+1)*num_cusps)**2 real number
    
    
    EXAMPLES::

        sage: R=9.533695261353557554344235235928770323821256395107
        sage: Y=0.5
        sage: G=MySubgroup(Gamma0(1))
        sage: W=setup_matrix_for_Maass_waveforms_np_dble(G,R,Y,10,20,1)
        sage: W.keys()
        ['Mf', 'nc', 'Ms', 'V']
        sage: print W['Ms'],W['Mf'],W['nc'];
        -10 10 1
        sage: V.rows,V.cols
        (21, 21)
        sage: N=set_norm_maass(1)
        sage: C=solve_system_for_Maass_waveforms(W['V'],N)
        sage: c2.real
        mpf('-1.0683335512235690597444640854010889121099661684073309')
        sage: C[0][2].real
        mpf('-1.0683335512235690597444640854010889121099661684073309')
        sage: C[0][3].real
        mpf('-0.4561973545061180270983019925431624077983243678097305')
        sage: C[0][6].real
        mpf('0.48737093979829448482994037031228389306060715149700209')
        sage: abs(C[0][6]-C[0][2]*C[0][3])
        mpf('0.00000000000002405188994758679782568437503269848800837')
        
    
    """
    cdef int l,j,icusp,jcusp,n,ni,li,iMl,iMs,iMf,iQs,iQf,iQl,s,nc
    cdef DTYPE_t kbes,tmpV,rtmpV,sqrtY,Y2pi,nrY2pi,argm,argpb,twopi,Qfak
    mpmath_ctx=mpmath.fp
    pi=mpmath_ctx.pi
    sqrtY=mpmath_ctx.sqrt(Y)    
    two=mpmath_ctx.mpf(2)
    Y2pi=Y*pi*two
    twopi=two*pi
    IR=mpmath_ctx.mpc(0,R)
    if(not isinstance(G,sage.modular.maass.mysubgroup.MySubgroup)):
        GG=MySubgroup(G)
    else:
        GG=G
    nc=int(GG.ncusps())
    if(Q<M):
        Q=M+20
    # Recall that we use symmetrization here
    Ms=1; Mf=M; Ml=Mf-Ms+1
    Qs=1; Qf=Q
    Qfak=mpmath_ctx.mpf(Q)/two
    iMl=int(Ml); iMs=int(Ms); iMf=int(Mf)
    iQf=int(Qf); iQs=int(Qs); iQl=int(Qf-Qs+1)
    # Pullpack points
    pb=pullback_pts(GG,Qs,Qf,Y,mp_ctx=mpmath_ctx)
    s=int(nc*iMl)
    cdef nvec=np.arange(iMs,iMf+1,dtype=DTYPE)
    cdef cnp.ndarray[DTYPE_t,ndim=2] V=np.zeros([int(s), int(s)], dtype=DTYPE)
    cdef cnp.ndarray[DTYPE_t,ndim=4] ef1=np.zeros([int(iQl), int(nc),int(nc),int(iMl)], dtype=DTYPE)         
    cdef cnp.ndarray[DTYPE_t,ndim=2] ef2=np.zeros([int(iMl), int(iQl)], dtype=DTYPE)
    if(sym_type==1):
        fun=mpmath_ctx.sin
    elif(sym_type==0):
        fun=mpmath_ctx.cos
    else:
        raise ValueError, "Only use together with symmetry= 0 or 1 ! Got:%s"%(sym_type)
    for n in range(iMs,iMf+1):        
        nr=nvec[n-iMs]
        for j in range(Qs,Qf+1):
            argm=nr*pb['xm'][j]; ef2[n-iMs,j-iQs]=fun(argm)
            for jcusp from 0<=jcusp <nc:
                for icusp from 0<=icusp <nc:
                    if(not pb['xpb'].has_key((icusp,jcusp,j))):
                        continue
                    argpb=nr*pb['xpb'][icusp,jcusp,j]; ef1[j-iQs,icusp,jcusp,n-iMs]=fun(argpb)

    for l from  iMs <= l <= iMf:
        #if(l==0 and cuspidal==1):
        #    continue
        lr=nvec[l-iMs]*twopi
        for j from iQs <= j <= iQf:
            for jcusp from int(0) <= jcusp < int(nc):
                lj=iMl*jcusp+l-iMs
                for icusp from int(0) <= icusp < int(nc):
                    if(not pb['ypb'].has_key((icusp,jcusp,j))):
                        continue
                    ypb=pb['ypb'][icusp,jcusp,j]
                    if(ypb==0):
                        continue
                    besarg=lr*ypb
                    kbes=mpmath_ctx.sqrt(ypb)*besselk_dp(R,besarg,pref=1)
                    kbes=kbes*ef1[int(j-iQs),int(icusp),int(jcusp),int(l-iMs)]
                    for n from iMs <= n <= iMf:
                        #if(cuspidal==1 and n==0):
                        #    continue
                        ni=iMl*icusp+n-iMs
                        rtmpV=kbes*ef2[int(n-iMs),int(j-iQs)]
                        V[ni,lj]+=rtmpV

    for n from 0<=n< V.shape[0]:
        for l from 0 <= l < V.shape[1]:        
            V[n,l]=V[n,l]/Qfak
    for n from iMs<=n <=iMf:
        #if(n==0 and cuspidal==1):
        #    continue
        nr=abs(nvec[n-iMs])
        for icusp from 0<=icusp<nc:
            ni=int(iMl*icusp+n-iMs)
            nrY2pi=nr*Y2pi
            kbes=sqrtY*besselk_dp(R,nrY2pi,pref=1)
            V[ni,ni]=V[ni,ni]-kbes

    W=dict()
    VV=mpmath_ctx.matrix(int(s),int(s),force_type=mpmath_ctx.mpf)
    for n from 0<=n< s:
        for l from 0 <= l < s:
            VV[n,l]=V[n,l]

    W['V']=VV
    W['Ms']=Ms
    W['Mf']=Mf
    W['nc']=nc
    return W




def setup_matrix_for_Maass_waveforms_np_cplx(G,double R,double Y,int M,int Q,int cuspidal=1):
    r"""
    Set up the matrix for the system of equations giving the Fourier coefficients of the Maass waveforms.

    INPUT:
    
    - ``G`` -- MySubgroup
    - ``R`` -- real (mpmath)
    - ``Y`` -- real (mpmath)
    - ``M`` -- integer
    - ``Q`` -- integer > M
    - ``cuspidal`` -- (False) logical : Assume cusp forms
    - ``sym_type``  -- (None) integer

          -  0 -- even / sine
          -  1 -- odd  / cosine

    - ``low_prec`` -- (False) : Use low (double) precision K-Bessel

    OUTPUT:
    
    - ``W`` -- dictionary
    - ``W['Mf']`` -- M start
    - ``W['nc']`` -- number of cusps
    - ``W['Ms']`` -- M stop
    - ``W['V']``  -- ((Ms-Mf+1)*num_cusps)**2 real number
    
    
    EXAMPLES::

        sage: R=9.533695261353557554344235235928770323821256395107251
        sage: Y=0.1
        sage: G=MySubgroup(Gamma0(1))
        sage: W=setup_matrix_for_Maass_waveforms_np_cplx(G,R,Y,10,20)
        sage: W.keys()
        ['Mf', 'nc', 'Ms', 'V']
        sage: print W['Ms'],W['Mf'],W['nc'];
        -10 10 1
        sage: W['V'].rows,W['V'].cols
        (21, 21)
        sage: N=set_norm_maass(1)
        sage: C=solve_system_for_Maass_waveforms(W['V'],N)
        sage: c2.real
        mpf('-1.0683335512235690597444640854010889121099661684073309')
        sage: C[0][2].real
        mpf('-1.0683335512235690597444640854010889121099661684073309')
        sage: C[0][3].real
        mpf('-0.4561973545061180270983019925431624077983243678097305')
        sage: C[0][6].real
        mpf('0.48737093979829448482994037031228389306060715149700209')
        sage: abs(C[0][6]-C[0][2]*C[0][3])
        mpf('0.00000000000002405188994758679782568437503269848800837')
        
    

    """
    cdef int l,j,icusp,jcusp,n,ni,li,iMl,iMs,iMf,iQs,iQf,iQl,s,nc
    cdef double sqrtY,Y2pi,nrY2pi,argm,argpb,twopi,two,kbes
    cdef double complex ckbes,ctmpV,iargm,iargpb,twopii,Qfak

    mpmath_ctx=mpmath.fp
    pi=mpmath_ctx.pi
    sqrtY=sqrt(Y)
    two=<double>(2)
    Y2pi=Y*pi*two
    twopi=two*pi
    if(not isinstance(G,sage.modular.maass.mysubgroup.MySubgroup)):
        GG=MySubgroup(G)
    else:
        GG=G
    nc=int(GG.ncusps())
    if(Q<M):
        Q=M+20
    Ms=-M;  Mf=M; Ml=Mf-Ms+1
    Qs=1-Q; Qf=Q
    Qfak=<double>(2*Q)
    iMl=int(Ml); iMs=int(Ms); iMf=int(Mf);
    iQf=int(Qf); iQs=int(Qs); iQl=int(Qf-Qs+1)
    # Pullpack points
    pb=pullback_pts(GG,Qs,Qf,Y,mp_ctx=mpmath_ctx)
    Xm=pb['xm']; Xpb=pb['xpb']; Ypb=pb['ypb']; Cv=pb['cvec']
    s=int(nc*iMl)
    cdef nvec=np.arange(iMs,iMf+1,dtype=DTYPE)
    cdef cnp.ndarray[CTYPE_t,ndim=2] V=np.zeros([s, s], dtype=CTYPE)
    #cdef cnp.ndarray[CTYPE_t,ndim=4] ef1=np.zeros([iQl, nc,nc,iMl], dtype=CTYPE)         
    #cdef cnp.ndarray[CTYPE_t,ndim=2] ef2=np.zeros([iMl, iQl], dtype=CTYPE)
    ef1=dict(); ef2=dict()
    for n in range(iMs,iMf+1):        
        nr=nvec[n-iMs]
        #print "nr(",n,")=",nr
        for j in range(Qs,Qf+1):
            argm=nr*Xm[j]
            ef2[n,j]=mpmath_ctx.exp(mpmath_ctx.mpc(0,-argm))
            for jcusp in range(nc):
                for icusp in range(nc):
                    if(not Xpb.has_key((icusp,jcusp,j))):
                        continue
                    #argpb=mpmath.mp(n)*Xpb[icusp,jcusp,j]
                    argpb=nr*Xpb[icusp,jcusp,j]
                    iargpb=mpmath.mp.mpc(0,argpb)
                    #if(n==1):
                    #    print "Xpb[",icusp,jcusp,j,"]=",Xpb[icusp,jcusp,j]
                    #    print "argpb =",argpb,abs(argpb-Xpb[icusp,jcusp,j])
                    #    print "iargpb=",iargpb,abs(abs(iargpb)-abs(argpb))
                    ef1[j,icusp,jcusp,n]=mpmath_ctx.exp(iargpb)

    for l from  iMs <= l <= iMf:
        if(l==0 and cuspidal==1):
            continue
        lr=nvec[l-iMs]*twopi
        for j from iQs <= j <= iQf:
            for jcusp from int(0) <= jcusp < int(nc):
                lj=iMl*jcusp+l-iMs
                for icusp from int(0) <= icusp < int(nc):
                    if(not Ypb.has_key((icusp,jcusp,j))):
                        continue
                    ypb=Ypb[icusp,jcusp,j]
                    if(ypb==0):
                        continue
                    besarg=abs(lr)*ypb
                    kbes=mpmath_ctx.sqrt(ypb)*besselk_dp(R,besarg,pref=1)
                    ckbes=kbes*ef1[j,icusp,jcusp,l]
                    for n from iMs <= n <= iMf:
                        if(n==0 and cuspidal==1):
                            continue
                        ni=iMl*icusp+n-iMs
                        ctmpV=ckbes*ef2[n,j]
                        V[ni,lj]=V[ni,lj]+ctmpV
##  if(ni==16 and lj==16):
##                             print "besarg(",j,")=",besarg
##                             print "ef2(",j,")=",ef2[n,j]
##                             print "ef1(",j,")=",ef1[j,icusp,jcusp,l]
##                             print "kbes(",j,")=",kbes
##                             print "ckbes(",j,")=",ckbes
##                             print "tmpV(",j,")=",ctmpV
##                             print "V[16,16](",j,")=",V[ni,lj]

#    print "V[16,16]=",V[16,16]
    for n from 0<=n< V.shape[0]:
        for l from 0 <= l < V.shape[1]:        
            V[n,l]=V[n,l]/Qfak

#    print "V[16,16]=",V[16,16]
#    print "sqrtY=",sqrtY
    for n from iMs<=n <=iMf:
        if(n==0 and cuspidal==1):
            continue
        nr=abs(nvec[n-iMs])
        #for icusp in range(nc):
        for icusp from 0<=icusp<nc:
            ni=iMl*icusp+n-iMs
            nrY2pi=nr*Y2pi
            ckbes=sqrtY*besselk_dp(R,nrY2pi,pref=1)
            #if(ni==16):
            #    print "arg=",nrY2pi
            #    print"ckbes(",ni,")=",besselk_dp(R,nrY2pi,pref=1)
            V[ni,ni]=V[ni,ni]-ckbes
    W=dict()
    #print "V[16,16]=",V[16,16]
    VV=mpmath_ctx.matrix(int(s),int(s),force_type=mpmath_ctx.mpc)
    for n from 0<=n< s:
        for l from 0 <= l < s:
            VV[n,l]=V[n,l]
    W['V']=VV
    W['Ms']=Ms
    W['Mf']=Mf
    W['nc']=nc
    return W

@cython.boundscheck(True)

def phase_2(G,R,Cin,int NA,int NB,ndig=10,M0=None,Yin=None, cuspidal=True,sym_type=None,low_prec=False):
    r"""
    Computes more coefficients from the given initial set.

      INPUT:
    
        - ``G`` -- MySubgroup
        - ``R`` -- real (mpmath)
        - ``Cin`` -- complex vector (mpmath)
        - ``NA`` -- integer 
        - ``NB`` -- integer 
        - ``M0`` -- integer (default None) if set we only use this many coeffficients.
        - ``ndig`` -- integer (default 10) : want new coefficients with this many digits precision.
        - ``cuspidal`` -- (False) logical : Assume cusp forms
        - ``sym_type``  -- (None) integer

              -  0 -- even / sine
              -  1 -- odd  / cosine

        - ``low_prec`` -- (False) : Use low (double) precision K-Bessel

    OUTPUT:

    - ``Cout`` -- dictionary of Fourier coefficients

    EXAMPLE::


        sage: R=mpmath.mpf(9.53369526135355755434423523592877032382125639510725198237579046413534)
        sage: IR=mpmath.mpc(0,R)
        sage: M=MaassWaveForms(Gamma0(1))
        sage: C=Maassform_coeffs(M,R)
        sage: D=phase_2(SL2Z,R,C,2,10) 

    

    """
    mpmath_ctx=mpmath.mp
    pi=mpmath.mp.pi()
    mp0=mpmath_ctx.mpf(0)
    mp1=mpmath_ctx.mpf(1)
    mp2=mpmath_ctx.mpf(2)
    eps=mpmath_ctx.power(mpmath_ctx.mpf(10),mpmath_ctx.mpf(-ndig))
    twopi=mp2*pi
    if(not isinstance(G,sage.modular.maass.mysubgroup.MySubgroup)):
        GG=MySubgroup(G)
    else:
        GG=G
    nc=int(GG.ncusps())
    IR=mpmath_ctx.mpc(0,R)
    if(M0<>None):
        M00=M0
        Mf=min ( M00, max(Cin[0].keys()))
    else:
        M00=infinity
        Mf= max(Cin[0].keys())
    if(sym_type<>None):
        Ms=1
    else:
        Ms=max ( -M00, min(Cin[0].keys()))

    Ml=Mf-Ms+1
    lvec=list()
    for l from Ms <= l <= Mf:
        lvec.append(mpmath_ctx.mpf(l))
    # starting value of Y
    if(Yin == None):
        Y0=twopi/mpmath.mp.mpf(NA)
        if(Y0>0.5):
            Y0=mp1/mp2
    else:
        Y0=Yin
    #
    NAa=NA
    ylp=0; Qp=0
    Cout=Cin # Start by assigning
    for yn in range(100):
        try:
            Yv=[Y0,Y0*mpmath.mpf(0.995)]

            Q=max(get_M_for_maass(R,Yv[1],eps)+5,Mf+10)+Qp
            Qs=1-Q; Qf=Q; Qfak=mpmath_ctx.mpf(2*Q)
            Xpb=dict(); Ypb=dict(); Cv=dict()
            pb=dict()
            pb[0]=pullback_pts(GG,Qs,Qf,Yv[0],mp_ctx=mpmath_ctx)
            pb[1]=pullback_pts(GG,Qs,Qf,Yv[1],mp_ctx=mpmath_ctx)
            Cnew=dict()
            print "Y[0]=",Yv,
            print "Ms,Mf=",Ms,Mf,"Qs,Qf=",Qs,Qf
            print "NA,NB=",NAa,NB
            kbdict=_setup_kbessel_dict(Ms,Mf,Qs,Qf,nc,IR,pb,lvec,mpmath_ctx)
            for n from NAa <= n <= NB:
                nr=mpmath.mp.mpf(n)
                for icusp from int(0) <= icusp < int(nc):
                    cvec=_get_tmpCs_twoY(n,icusp,Yv,IR,Ms,Mf,Qs,Qf,Qfak,nc,pb,Cin,kbdict,lvec,mpmath_ctx)
                    #                cvec=_get_tmpCs_twoY(n,icusp,Yv,IR,Ms,Mf,Qs,Qf,Qfak,nc,lvec,pb,Cin,kbdict,mpmath_ctx)
                    # print "c1,c2=",cvec
                    diff=abs(cvec[0]-cvec[1])
                    print "err[",n,"]=",diff
                    if(diff<eps):
                        Cnew[icusp,n]=cvec[1]
                        # # improve the coefficients we use
                        if(Cin.has_key((icusp,n))):
                            Cin[icusp,n]=cvec[1]
                    else:
                        NAa=n; ylp=ylp+1
                        if(ylp>2):
                            Qp=Qp+10; ylp=0
                        else:
                            Y0=Y0*mpmath.mpf(0.99)
                            #Qp=0

                        raise StopIteration()
                    
        except StopIteration:
            continue
        return Cnew


def _setup_kbessel_dict(Ms,Mf,Qs,Qf,nc,IR,pb,lvec=None,mpmath_ctx=mpmath.mp):
    r"""
    Setup a dictionary of K-Bessel values.

    INPUT:

    -``Ms`` -- integer
    -``Mf`` -- integer
    -``Qs`` -- integer
    -``Qs`` -- integer
    -``nc`` -- integer
    -``IR`` -- complex (purely imaginary) 
    -``pb`` -- dictionary 
    -``lvec`` -- (optional) vector of mpmath real integers in range(Ms,Mf+1)
    -``mpmath_ctx`` -- mpmath context

    OUTPUT:
    -``kbdict`` dictionary of dimensions
          Ms:Mf,Qs:Qf,0:nc-1,0:nc-1,0:1

    EXAMPLES::


    Note that this function is not in the public name space.
    Therefore we need a workaround to call it from the command line.

        sage: Yv=[mpmath.mpf(0.5),mpmath.mpf(0.45)]
        sage: IR=mpmath.mpc(0,9.53369526135355755434423523592877032382125639510725198237579046413534)
        sage: M=20; Ms=-M; Mf=M; Q=30; Qs=1-Q; Qf=Q
        sage: pb=dict()
        sage: pb[0]=pullback_pts(Gamma0(1),Qs,Qf,Yv[0])
        sage: pb[1]=pullback_pts(Gamma0(1),Qs,Qf,Yv[1])
        sage: _temp = __import__('maass_forms_alg',globals(), locals(), ['_setup_kbessel_dict'],-1)
        sage: kbdict=_temp.__dict__['_setup_kbessel_dict'](Ms,Mf,Qs,Qf,1,IR,pb)    
    """
    kbvec=dict()
    twopi=mpmath_ctx.mpf(2)*mpmath_ctx.pi()
    for l from Ms <= l <= Mf:
        if(lvec<>None):
            lr=lvec[l-Ms]
        else:
            lr=mpmath_ctx.mpf(l)
        for j from Qs <= j <= Qf:
            for icusp from int(0) <= icusp < int(nc):
                for jcusp from int(0) <= jcusp < int(nc):
                    for yj in range(2):
                        ypb=pb[yj]['ypb'][icusp,jcusp,j]
                        ypbtwopi=ypb*twopi
                        besarg=abs(lr)*ypbtwopi
                        kbes=mpmath_ctx.sqrt(ypb)*(mpmath_ctx.besselk(IR,besarg).real)
                        kbvec[l,j,icusp,jcusp,yj]=kbes
    return kbvec


def _get_tmpCs_twoY(n,ic,Yv,IR,Ms,Mf,Qs,Qf,Qfak,nc,pb,Cin,kbdict,lvec=None,mpmath_ctx=mpmath.mp):
    r"""
    Compute two coeffcients given by the two Y-values.

    INPUT:

    -``n`` -- integer
    -``ic`` -- integer
    -``Yv`` -- list of two reals Y0,Y1
    -``IR`` -- complex  IR
    -``Ms`` -- integer
    -``Mf`` -- integer
    -``Qs`` -- integer
    -``Qf`` -- integer
    -``nc`` -- integer
    -``lvec`` -- list of mpmath reals in range(Ms,Mf)
    -``pb`` -- list of dictonaries of values related to the pullback
    -``Cin`` -- list of reals : Fourier coefficients
    -``kbdict-- dictionary of K-Besssel values
    -``mpmath_ctx`` -- mpmath context

    OUTPUT:

    --``cvec`` -- list of complex numbers : Fourier coefficients  c_1[n,ic] and c_2[n,ic]


    EXAMPLES::


    Note that this function is not in the public name space.
    Therefore we need a workaround to call it from the command line.
        sage: Yv=[mpmath.mpf(0.2),mpmath.mpf(0.17)]
        sage: R=mpmath.mpf(9.53369526135355755434423523592877032382125639510725198237579046413534)
        sage: IR=mpmath.mpc(0,R)
        sage: M=MaassWaveForms(Gamma0(1))
        sage: C=Maassform_coeffs(M,R)
        sage: M=max(C[0].keys()); Ms=-M; Mf=M
        sage: Q=max(get_M_for_maass(R,Yv[1],1E-12),M+10); Qs=1-Q; Qf=Q; nc=1; ic=0
        sage: Qfak=mpmath.mpf(2*Q)
        sage: pb=dict()
        sage: pb[0]=pullback_pts(Gamma0(1),Qs,Qf,Yv[0])
        sage: pb[1]=pullback_pts(Gamma0(1),Qs,Qf,Yv[1])
        sage: _temp = __import__('maass_forms_alg',globals(), locals(), ['_setup_kbessel_dict'],-1)
        sage: kbdict=_temp.__dict__['_setup_kbessel_dict'](Ms,Mf,Qs,Qf,nc,IR,pb)  
        sage: n=10
        sage: c=_temp.__dict__['_get_tmpCs_twoY'](n,ic,Yv,IR,Ms,Mf,Qs,Qf,Qfak,nc,pb,C,kbdict)
        sage: C[0][n].real
        mpf('0.31053524297612709')
        sage: c[0]
        mpc(real='0.31053524291041912', imag='9.265227235203014e-17')
        sage: c[0]
        mpc(real='0.31053524291042567', imag='-3.8363741205534537e-17')
        sage:     sage: abs(c[0]-c[1])
        mpf('6.5516259713787926e-15')l

    
    """
    ctmp=dict(); tmpV=dict()
    nr=mpmath_ctx.mpf(n)
    twopi=mpmath_ctx.mpf(2)*mpmath_ctx.pi()
    #print "keys=",Cin[0].keys()
    #print "Yv=",Yv
    for yj in range(2):
        Y=Yv[yj]
        summa=mpmath_ctx.mpf(0)
        for jc from int(0) <= jc < int(nc):
            for l from int(Ms) <= l <= int(Mf):
                if(Cin[jc][l]==0):
                    continue
                if(lvec<>None):
                    lr=lvec[l-Ms]
                else:
                    lr=mpmath_ctx.mpf(l)
                Vnlij=mp0
                for j from int(Qs) <= j <= int(Qf):
                    if(not pb[yj]['ypb'].has_key((ic,jc,j))):
                        continue
                    ypb=pb[yj]['ypb'][ic,jc,j]
                    if(ypb==0):
                        continue
                    kbes=kbdict[l,j,ic,jc,yj]
                    xpb=pb[yj]['xpb'][ic,jc,j]
                    arg=lr*xpb-nr*pb[yj]['xm'][j]
                    tmp=mpmath.mp.exp(mpmath.mpc(0,arg))
                    Vnlij=Vnlij+kbes*tmp
                summa=summa+Vnlij*Cin[jc][l]
        tmpV[yj]=summa/Qfak
        besarg=nr*twopi*Y
        besY=mpmath_ctx.sqrt(Y)*(mpmath_ctx.besselk(IR,besarg).real)
        ctmp[yj]=tmpV[yj]/besY
        #print "summa[",yj,j,l,"]=",summa
        #print "tmpV[",yj,j,l,"]=",tmpV[yj]
        #print "C[",n,yj,"]=",ctmp[yj]
    return ctmp



def smallest_inf_norm(V):
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



def get_M_for_maass(R,Y,eps):
    r""" Computes the  truncation point for Maass waveform with parameter R.

    INPUT:
    
     - ``R`` -- spectral parameter, real
     - ``Y`` -- height, real > 0 
     - ``eps`` -- desired error
     

    OUTPUT:

     - ``M`` -- point for truncation

    EXAMPLES::ubuntuusers.de


         sage: get_M_for_maass(9.533,0.5,1E-16)
         12
         sage: get_M_for_maass(9.533,0.5,1E-25)
         18


    """
    ## Use low precision
    dold=mpmath.mp.dps
    mpmath.mp.dps=int(mpmath.ceil(abs(mpmath.log10(eps))))+5
    twopi=2*mpmath.pi()
    twopiy=twopi*mpmath.mpf(Y)
    # an extra two for the accumulation of errors
    minv=(mpmath.mpf(12)*mpmath.power(R,0.3333)+R)/twopiy
    minm=mpmath.ceil(minv)
    try:
        for m in range(minm,10000):
            erest=err_est_Maasswf(Y,m)
            if(erest<eps):
                raise StopIteration()
    except StopIteration:
        mpmath.mp.dps=dold
        return m
    mpmath.mp.dps=dold
    raise Exception," No good value for truncation was found!"


def err_est_Maasswf(Y,M):
    r"""
    Estimate the truncated series: $|\Sum_{n\ge M} a(n) \sqrt{Y} K_{iR}(2\pi nY) e(nx)|$

    CAVEATS:
    - we assume the Ramanujan bound for the coefficients $a(n)$, i.e. $|a(n)|\le2$.
    - error estimate holds for 2piMY >> R

    INPUT:

    - ``Y`` -- real > 0
    - ``M``  -- integer >0

    OUTPUT:

    - ``r``  -- error estimate

    EXAMPLES::


        sage: err_est_Maasswf(0.5,10)
        mpf('3.1837954148201557e-15')
        sage: err_est_Maasswf(0.85,10)
        mpf('5.3032651127544969e-25')

    """
    arg=mpmath.mp.pi()*mpmath.mpf(Y)
    r=mpmath.mpf(1)/arg.sqrt()
    arg=arg*mpmath.mpf(2*M)
    r=r*mpmath.gammainc(mpmath.mpf(0.5),arg)
    return r


def pullback_to_psl2z_mp(x,y,word=False):
    r"""  Pullback to the fundamental domain of SL2Z
          MPMATH version including option for solving the word problem

    INPUT:

        - ``x``  -- double
        - ``x``  -- double

    OUTPUT:

    - ``[u,v,[a,b,c,d]`` --
        u,v -- double
        [a,b,c,d] 4-tuple of integers giving the matrix A
                  such that A(x+iy)=u+iv is inside the fundamental domain of SL2Z

    EXAMPLES::


        sage: pullback_to_psl2z_mp(mpmath.mpf(0.1),mpmath.mpf(0.1))
        [mpf('6.9388939039072274e-16'), mpf('4.9999999999999991'), [5, -1, 1, 0]]
        sage: [x,y,A,w]=pullback_to_psl2z_mp(mpmath.mpf(0.1),mpmath.mpf(0.1),True)
        sage: SL2Z(A)
        [ 5 -1]
        [ 1  0]
        sage: w
        [['T', 5], ['S', 1]]
        sage: S,T=SL2Z.gens()
        sage: eval(w[0][0])**w[0][1]*eval(w[1][0])**w[1][1]
        [ 5 -1]
        [ 1  0]
        sage: T^5*S
        [ 5 -1]
        [ 1  0]

        
    """
    try:
        x.ae; y.ae
    except AttributeError:
        raise Exception,"Need x,y in mpmath format!"
    l=[]
    minus_one=mpmath.mpf(-1)
    cdef double xr,yr
    xr=<double> x
    yr=<double> y
    half=<double> 0.5
    cdef int imax=10000
    cdef int i
    cdef int a,b,c,d,aa,bb
    a=1; b=0; c=0; d=1
    cdef double dx
    cdef int numT
    for i from 0<=i<=imax:
        if(fabs(xr)>0.5):
            dx=max(xr-half,-xr-half)
            numT=ceil(dx)
            if(xr>0):
                numT=-numT
            xr=xr+<double>numT
            a=a+numT*c
            b=b+numT*d
            if(word==True):
                l.append(['T',numT])
        else:
            # We might have to flip
            absval=xr*xr+yr*yr
            if(absval<1.0-1E-15):
                xr=minus_one*xr/absval
                yr=yr/absval
                aa=a
                bb=b
                a=-c
                b=-d
                c=aa
                d=bb
                if(word==True):
                    l.append(['S',1])
            else:
                break
    # Apply the matrix to x+iy
    ar=mpmath.mpf(a);    br=mpmath.mpf(b);
    cr=mpmath.mpf(c);    dr=mpmath.mpf(d);
    
    den=(cr*x+dr)**2+(cr*y)**2
    yy=y/den
    xx=(ar*cr*(x*x+y*y)+x*(ar*dr+br*cr)+br*dr)/den
    if(word==True):
        l.reverse() # gives a word for A
        return [xx,yy,[a,b,c,d],l]
    else:
        return [xx,yy,[a,b,c,d]]

cdef tuple pullback_to_psl2z_dble(double x,double y):
    r""" Pullback to the fundamental domain of SL2Z
         Fast (typed) Cython version
    INPUT:

        - ``x``  -- double
        - ``x``  -- double

    OUTPUT:

    - ``[u,v,[a,b,c,d]`` --
        u,v -- double
        [a,b,c,d] 4-tuple of integers giving the matrix A
                  such that A(x+iy)=u+iv is inside the fundamental domain of SL2Z

    EXAMPLES::

        sage: pullback_to_psl2z_dble(0.1,0.1)
        [6.9388939039072274e-16, 4.9999999999999991, [5, -1, 1, 0]]


    """
    cdef int a,b,c,d
    [a,b,c,d]=pullback_to_psl2z_mat(x,y)
    cdef double ar=<double>a
    cdef double br=<double>b
    cdef double cr=<double>c
    cdef double dr=<double>d
    cdef double den=(cr*x+dr)**2+(cr*y)**2
    cdef double yout=y/den
    cdef double xout=(ar*cr*(x*x+y*y)+x*(ar*dr+br*cr)+br*dr)/den
    cdef tuple t=(a,b,c,d)
    cdef tuple tt=(xout,yout,t)
    return tt
    #return [xout,yout,[a,b,c,d]]
    
    
cdef tuple pullback_to_psl2z_mat(double x,double y):
    r""" Pullback to the fundamental domain of SL2Z
         Fast (typed) Cython version
    INPUT:

        - ``x``  -- double
        - ``x``  -- double

    OUTPUT:

    - ``[a,b,c,d]`` -- 4-tuple of integers giving the matrix A
                      such that A(x+iy) is inside the fundamental domain of SL2Z

    EXAMPLES::

        sage: pullback_to_psl2z_mat(0.1,0.1)
        [5, -1, 1, 0]

    
    """
    cdef int imax=10000
    cdef int done=0
    cdef int a,b,c,d,aa,bb,i
    a=1; b=0; c=0; d=1
    cdef double half=<double> 0.5
    cdef double minus_one=<double> -1.0
    cdef int numT
    cdef double absval,dx
    for i from 0<=i<=imax:
        if(fabs(x)>half):
            dx=max(x-half,-x-half)
            numT=ceil(dx)
            if(x>0):
                numT=-numT
            x=x+<double>numT
            a=a+numT*c
            b=b+numT*d
        else:
            # We might have to flip
            absval=x*x+y*y
            if(absval<1.0-1E-15):
                x=minus_one*x/absval
                y=y/absval
                aa=a
                bb=b
                a=-c
                b=-d
                c=aa
                d=bb
            else:
                break
    cdef tuple t=(a,b,c,d)
    return t




def mppr(x):
    r"""
    Print fewer digits!

    """
    
    # We don't want to print 500 digits
    # So I truncate at two digits with a stupid simple rounding
    dpold=mpmath.mp.dps
    mpmath.mp.dps=5
    s=str(x)
    mpmath.mp.dps=dpold
    (s1,dot,s2)=s.partition(".")
    (s3,es,s4)=s2.partition("e")
    if(len(s4)==0):  # No "e" format
        return s[0:15]
    try:
        ma=s3[0]
        if(len(s3)>1):
            if(s3[2]<5):
                ma=ma+s3[1]
            elif(s3[2]>=5):
                ma=ma+str(Integer(s3[1])+1)
        ss=s1+"."+ma+"e"+s4
    except:
        return s
    return ss
