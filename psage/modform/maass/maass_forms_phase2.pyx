# cython: profile=True
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
r"""
Algorithms for phase 2 for Maass waveforms
"""
include 'common_defs.pxd'


from sage.rings.complex_mpc import MPComplexField
from sage.rings.real_mpfr import RealField
from sage.all import CC
from maass_forms import Maasswaveform,maass_logger

import cython
cdef extern from "math.h":
    double fabs(double)
    double fmax(double,double)
    int ceil(double)
    double sqrt(double)
    double sin(double)
    double cos(double)
    double power(double,double)
    double M_PI

cdef extern from "complex.h":
    ### "Hack" suggested by Robert Bradshaw in 2009.
    ### TODO: Is there now a better way for pure c double complex?
    ctypedef double cdouble "double complex"
    cdef double creal(double complex)
    cdef double cimag(double complex)
    cdef double complex _Complex_I
    cdef double carg(double complex)
    cdef double cabs(double complex)
    cdef double complex cexp(double complex)
    cdef double complex csqrt(double complex)
    cdef double complex cpow(double complex,double complex)


cdef int gcd( int a, int b ):
    cdef int c
    while ( a != 0 ):
        c = a; a = b%a;  b = c
    return b

cdef double complex cexpi(double x):
    return cexp(x*_I)
#ctypedef void (*foo_t)(int)
ctypedef complex (*complex_function_type)(complex)


cdef double complex _I = _Complex_I
cdef double complex _I2 = _Complex_I + _Complex_I

cdef extern from "mpfr.h":
    int mpfr_mul_d (mpfr_t, mpfr_t, double, mp_rnd_t)

from lpkbessel import besselk_dp
from pullback_algorithms cimport pullback_pts_cplx_dp,pullback_pts_real_dp,pullback_pts_cplx_dp_sym
from lpkbessel cimport besselk_dp_c
from maass_forms_alg cimport set_Mv_Qv_symm,get_M_for_maass_dp_c
from maass_forms_alg import get_Y_from_M,get_M_and_Y
#from maass_forms import get_primitive_p,Hecke_eigenfunction_from_coeffs
from maass_forms import dict_depth


cpdef phase2(F,n,verbose=0,retf=0,fnr=-1,n_step=50,do_test=1,method='2c'):
    r""" Computes more coefficients for the maass form F """
    from sage.all import log_b,floor
    nd = floor(abs(log_b(F.test(format='float'),10)))
    C = phase_2_cplx_dp_sym(F.space(),F.eigenvalue(),2,n,M0in=int(F._M0),ndig=int(nd-1),dim=int(F._dim),cuspstart=int(0),cuspstop=int(F.group().ncusps()),fnr=int(fnr),Cin=F._coeffs,Yin=float(F._Y),verbose=verbose,retf=int(retf),n_step=int(n_step),do_test=int(do_test),method=method)
    for r in F._coeffs.keys():
        if not C.has_key(r):
            C[r]={}
            for ci in F._coeffs[r].keys():
                if not C[r].has_key(ci):
                    C[r][ci]={}
                for n in F._coeffs[r][ci].keys():
                    C[r][ci][n] = F._coeffs[r][ci][n]
        return C


@cython.boundscheck(False)
@cython.cdivision(True)
cpdef phase_2_cplx_dp_sym(S,double R,int NA,int NB,int M0in=0,int ndig=10,int dim=1,int cuspstart=0,int cuspstop=1,int fnr=-1,dict Cin={},double Yin=0,int verbose=0,int retf=0,int n_step=50,int do_test=1,method='2c',int ynmax=1000):
    #,dict c_evs={}):
    r"""
    Computes more coefficients from the given initial set.

    INPUT:

    - ``S`` -- Space of Maas waveforms
    - ``R`` -- real
    - ``Cin`` -- complex vector
    - ``NA`` -- integer
    - ``NB`` -- integer
    - ``M0`` -- integer (default None) if set we only use this many coeffficients.
    - ``ndig`` -- integer (default 10) : want new coefficients with this many digits precision.
    - ``ndig`` -- number of digits of precision
    - ``dim`` -- assumed dimension
    - ``cuspstart`` -- Compute coefficients between cusp cuspstart and cuspstop
    - ``cuspstop`` --
    - ``fnr`` -- use function nr. fnr (in case of dimension>1)
    - ``method`` -- '2c' or '2Y' for using two cusps (if possible) or two Y's for error test.

    OUTPUT:

    - ``Cout`` -- dictionary of Fourier coefficients

    EXAMPLE::


    sage: R=RR(9.53369526135355755434423523592877032382125639510725198237579046413534)
    sage: IR=CC(0,R)
    sage: M=MaassWaveForms(Gamma0(1))
    sage: C=Maassform_coeffs(M,R)
    sage: D=phase_2(SL2Z,R,C,2,10)



    """
    cdef double eps,pi,twopi,Y
    cdef int nc,sym_type,M00,Mf,Ms,Ml
    cdef int Ql
    cdef int **Mv=NULL
    cdef int **Qv=NULL
    cdef int *symmetric_cusps=NULL
    cdef double **Xm=NULL
    cdef double ****Xpb=NULL
    cdef double ****Ypb=NULL
    cdef double complex ****Cvec=NULL
        #cdef double *Xm2=NULL,***Xpb2=NULL,***Ypb2=NULL
        #cdef double complex ***Cvec2=NULL
    cdef double complex *cusp_evs=NULL
    cdef double complex ***Cold=NULL
    cdef double * Qfak=NULL
    cdef int p,n,i,N1=0,NAa=NA
    cdef int ylp,Qp,Q #,verbose
    cdef dict NN,XS
    cdef double complex ***Ctmp
    cdef double besprec,diff,diffmax,Y02
    cdef double **nr
    cdef double **kbes
    cdef double *sqrtY
    cdef double *Yv
    cdef double *Y2pi
    cdef int yi,M0
    cdef int numy=1
    G=S.group(); nc=int(G.ncusps()); sym_type=S._sym_type
    if nc==1 or S.cusp_symmetries()[1][0]<>1:
        method='TwoY'
        if method=='TwoY':
            numy=2
    Yv=<double*>sage_malloc(sizeof(double)*numy)
    Y2pi=<double*>sage_malloc(sizeof(double)*numy)
    #cdef double nr0,nr1kbes,kbes0,kbes1
    #verbose=S._verbose
    pi=M_PI
    twopi=pi*2.0
    eps=10.0**-ndig
    cdef int fstart,fstop

    if fnr<0 or fnr>dim:
        fstart=0; fstop=dim
    else:
        fstart=fnr; fstop=fnr+1
    maass_logger.debug("method={0}\n eps ={1}, NA,NB={2},{3}".format(method,eps,NA,NB))
    cdef list c_evs=S.cusp_evs()
    if not len(c_evs)==nc:
        c_evs = [0 for i in range(len(c_evs),nc)]
    cusp_evs=<double complex*>sage_malloc(sizeof(double complex)*nc)
    for i in range(nc):
        cusp_evs[i]=<double complex> CC(c_evs[i])
    Ctmp = <double complex***>sage_malloc(sizeof(double complex**)*numy)
    for yi in range(numy):
        Ctmp[yi] = <double complex**>sage_malloc(sizeof(double complex*)*(fstop-fstart))
        for i in range(fstop-fstart):
            Ctmp[yi][i] = <double complex*>sage_malloc(sizeof(double complex)*nc*2)
            for j in range(nc*2):
                Ctmp[yi][i][j]=0
    Qfak=<double*>sage_malloc(sizeof(double)*nc)
    if Qfak==NULL: raise MemoryError
    Mv=<int**>sage_malloc(sizeof(int*)*nc)
    if Mv==NULL: raise MemoryError
    Qv=<int**>sage_malloc(sizeof(int*)*nc)
    if Qv==NULL: raise MemoryError
    symmetric_cusps=<int*>sage_malloc(sizeof(int)*nc)
    if symmetric_cusps==NULL: raise MemoryError
    for i in range(nc):
        Mv[i]=<int*>sage_malloc(sizeof(int)*3)
        Qv[i]=<int*>sage_malloc(sizeof(int)*3)
    cdef int* cusp_offsets
    cusp_offsets=<int*>sage_malloc(sizeof(int)*nc)
    if dim<=0:
        if hasattr(S,"_dim"):
            dim=S._dim
        else:
            dim=1
    NN=S.set_norm(dim)
    maass_logger.debug("R,Yin,M0in,eps={0}".format((R, Yin, M0in, eps)))
    cdef int old_v=S._verbose

    S._verbose=0
    ### If we don't have any coefficients as input we have to compute some. 
    if Cin=={}:
        M0,Y = get_M_and_Y(R,Yin,M0in,eps,verbose-2)
        maass_logger.debug("Y,M0={0}".format((Y,M0)))
        maass_logger.debug("fstart,fstop={0}".format((fstart,fstop)))
        #F = #S.get_element(R,dim=dim,Mset=M0)
        F = Maasswaveform(S,R,dim=dim,Mset=M0,compute=True)
        if verbose>1:
            maass_logger.debug("dim={0} \n M0={1}".format(dim,M0))            
            if dim >1 :
                for j in range(dim):
                     maass_logger.debug("F[{0}].test()={1}".format(j,F[j].test()))
            else:
                 maass_logger.debug("F.test()={0}".format(F.test()))
        if dim>1:
            if fnr>=0 and fnr<dim:
                XS = F[fnr]._coeffs
            else:
                XS=dict()
                for i in range(fstop-fstart):
                    if isinstance(F,dict):
                        if not F.has_key(i):
                            continue
                    elif isinstance(F,list):
                        if i >= len(F):
                            continue
                    XS[i] = F[i]._coeffs[0]
        else:
            XS = F._coeffs
        if verbose>1:
            maass_logger.debug("XS={0}".format(XS))
        #sig_on()
        #XS=get_coeff_fast_cplx_dp_sym(S,R,Y,M0,Q,NN,cusp_ev=c_evs)
        #sig_off()
    else:
        M0 = M0in
        XS=Cin
    Q = M0+10
    ## Make sure that we have a real eigenvalue.
    S._verbose=old_v
    set_Mv_Qv_symm(S,Mv,Qv,Qfak,symmetric_cusps,cusp_evs,cusp_offsets,&N1,&Ml,&Ql,M0,Q,verbose)
    maass_logger.debug("Get Q={0} for Y0={1}".format(Q,Y0))
    for i in range(nc):
        maass_logger.debug("Qv[{0}]={1}".format(i,(Qv[i][0],Qv[i][1],Qv[i][2])))
        maass_logger.debug("Qfak[{0}]={1}".format(i,(Qfak[i])))
    if verbose>3:
        maass_logger.debug("XS={0}".format(XS))
        maass_logger.debug("depth={0}".format(dict_depth(XS)))
    #print "XSK=",XS.keys()
    if do_test:
        F = Maasswaveform(S,R,C=XS)
        #test = S.test_Hecke_relation(XS[0][0])
        t1=  F.test(method='Hecke',format='float')
        t2 =  F.test(method='pcoeff',format='float')
        test = min(t1,t2)
        maass_logger.debug("Hecke test={0}".format(t1))
        maass_logger.debug("P-coeff test={0}".format(t2))
        maass_logger.debug("eps={0}".format(eps))
    else:
        test=0

    if abs(test)>eps:
        maass_logger.warning("Need to improve this eigenvalue, or decrease the required precision!")
        return {}

    cdef dict Cnew=dict()
    for j in XS.keys():
        Cnew[j]={}
        for i in XS[j].keys():
            Cnew[j][i]={}
    cdef int cuspa,cuspb
    if retf or nc==1:
        cuspa=0; cuspb=nc; cuspbb=nc
        #cuspa=0; cuspb=2
    else:
        cuspa=cuspstart; cuspb=cuspstop; cuspbb=max(cuspstop,2)

    ## First add the coefficients we already have
    for j in range(fstop-fstart):
        if not XS.has_key(j+fstart): continue
        for i in range(cuspa,cuspb):
            if not XS[j+fstart].has_key(i): continue
            if not isinstance(XS[j+fstart][i],dict):
                for n in XS[j+fstart][0].keys():
                    Cnew[j][i][n]=XS[j+fstart][i]*XS[j+fstart][0][n]
            else:
                for n in XS[j+fstart][i].keys():
                    Cnew[j][i][n]=XS[j+fstart][i][n]
    #print "CK=",Cnew[2].keys()
    maass_logger.debug("after!")
    Cold=<double complex***>sage_malloc(sizeof(double complex**)*(fstop-fstart))
    for j in range(fstop-fstart):
        Cold[j]=<double complex**>sage_malloc(sizeof(double complex*)*nc)
        for i in range(nc):
            if cusp_evs[i]==0 or i==0:
                Cold[j][i]=<double complex*>sage_malloc(sizeof(double complex)*Mv[i][2])
                for n in range(Mv[i][2]):
                    Cold[j][i][n]=<double complex>CC(XS[j][i][n+Mv[i][0]])
            else:
                Cold[j][i]=<double complex*>sage_malloc(sizeof(double complex)*Mv[i][2])
                for n in range(Mv[i][2]):
                    Cold[j][i][n]=cusp_evs[i]*XS[j][0][n+Mv[i][0]]
    # using these we want to compute more
    cdef int nn0,nn1
    maass_logger.debug("got initial set of coefficients!")
    maass_logger.debug("fstart={0}".format(fstart))
    maass_logger.debug("stop={0}".format(fstop))
    if verbose>1:
        for j in range(fstop-fstart):
            for n from Mv[0][0]<=n<Mv[0][1]:
                maass_logger.debug("C[{0}][0][{1}]={2}".format(j,n,XS[j][0][n]))
                maass_logger.debug("Cold[{0}][0][{1}]={2}".format(j,n,XS[j][0][n]))                
            for i in range(nc):
                if cusp_evs[i]==0:
                    maass_logger.debug("C[{0}][{1}][1]={2}".format(j,i,XS[j][i][1]))                    
                else:
                    maass_logger.debug("Cusp_ev[{0}]={1}".format(i,cusp_evs[i]))
            for i in range(nc):
                if cusp_evs[i]==0:
                    nn0 = max(-10,Mv[i][0])
                    nn1 = min(10,Mv[i][1])
                    for n in range(nn0,nn1):
                        if n>Mv[i][2]:# or n<0:
                            #print "Cold[{0}][{1}][{2}]=?".format(j,i,n)
                            continue
                        maass_logger.debug("Cold[{0}][{1}][{2}]={3}".format(j,i,n,Cold[j][i][n-Mv[i][0]]))
                    # starting value of Y
    cdef double Y0
    if Yin <=0:
        Y0=twopi/<double>NA
        if Y0>0.5:
            Y0=0.5
    else:
        Y0=Yin
    Y0 = min(Y0,S.group().minimal_height())
    Yv[0]=Y0*0.7
    if numy>1:
        #        Yv[1]=0.85*Y0
        Yv[1]=0.995*Yv[0]
    ylp=0; Qp=0
    Ml=2*M0+1
    sig_on()
    Q=max(get_M_for_maass_dp_c(R,Y0,eps)+5,Ml+10)+Qp
    set_Mv_Qv_symm(S,Mv,Qv,Qfak,symmetric_cusps,cusp_evs,cusp_offsets,&N1,&Ml,&Ql,M0,Q,verbose)
    maass_logger.debug("Get Q={0} for Y0={1}".format(Q,Y0))
    for i in range(nc):
        maass_logger.debug("Qv[{0}]={1}".format(i,(Qv[i][0],Qv[i][1],Qv[i][2])))
        maass_logger.debug("Qfak[{0}]={1}".format(i,(Qfak[i])))
    sig_off()
    #Ql=2*Q
    Xm=<double**>sage_malloc(sizeof(double*)*numy)
    if Xm==NULL: raise MemoryError
    for yi in range(numy):
        Xm[yi]=<double*>sage_malloc(sizeof(double)*Q)
        if Xm[yi] is NULL: raise MemoryError
    maass_logger.debug("Allocated Xm of len {0} x {1}".format(numy,Q))
    Xpb = <double****> sage_malloc( sizeof(double***) * numy )
    if Xpb==NULL: raise MemoryError
    Ypb = <double****> sage_malloc( sizeof(double*** ) * numy )
    if Ypb==NULL: raise MemoryError
    Cvec = <double complex****>sage_malloc(sizeof(double complex***) * numy )
    if Cvec==NULL: raise MemoryError
    for yi in range(numy):
        Xpb[yi] = <double***> sage_malloc( sizeof(double** ) * nc )
        if Xpb[yi]==NULL: raise MemoryError
        Ypb[yi] = <double***> sage_malloc( sizeof(double** ) * nc )
        if Ypb[yi]==NULL: raise MemoryError
        Cvec[yi] = <double complex***>sage_malloc(sizeof(double complex**) * nc )
        if Cvec[yi]==NULL: raise MemoryError
        for i in range(nc):
            Xpb[yi][i]=NULL; Ypb[yi][i]=NULL; Cvec[yi][i]=NULL
            Xpb[yi][i] = <double**>sage_malloc(sizeof(double*) * nc )
            if Xpb[yi][i]==NULL: raise MemoryError
            Ypb[yi][i] = <double**>sage_malloc(sizeof(double*) * nc )
            if Ypb[yi][i]==NULL: raise MemoryError
            Cvec[yi][i] = <double complex**>sage_malloc(sizeof(double*) * nc )
            if Cvec[yi][i]==NULL: raise MemoryError
            for j in range(nc):
                Xpb[yi][i][j]=NULL; Ypb[yi][i][j]=NULL; Cvec[yi][i][j]=NULL
                Xpb[yi][i][j] = <double*>sage_malloc(sizeof(double) * Ql )
                if Xpb[yi][i][j]==NULL: raise MemoryError
                Ypb[yi][i][j] = <double*>sage_malloc(sizeof(double) * Ql )
                if Ypb[yi][i][j]==NULL: raise MemoryError
                Cvec[yi][i][j] = <double complex*>sage_malloc(sizeof(double complex) * Ql )
                if Cvec[yi][i][j]==NULL: raise MemoryError
                for n in range(Ql):
                    Xpb[yi][i][j][n]=<double>0
                    Ypb[yi][i][j][n]=<double>0
                    Cvec[yi][i][j][n]=0

    maass_logger.debug("in phase 2 with eps={0}".format(eps))
    maass_logger.debug("range: NA,NB={0}".format((NA,NB)))
    maass_logger.debug("N1={0}".format(N1))
    maass_logger.debug("Ml,Ql={0},{1}".format(Ml,Ql))
    cdef double *alphas=NULL
    alphas=<double*>sage_malloc(sizeof(double)*nc)
    for i in range(nc):
        alphas[i]=<double>S.alpha(i)[0]
    cdef double complex **** V=NULL
    #cdef int n_step=50
    cdef int n_a,n_b
    if n_step > (NB - NA +1):
        n_step = NB-NA+1
    assert n_step > 0
    n_a = NA  ## Problem: - Mv[0][0]
    n_b = n_a + n_step

    maass_logger.debug("Compute coefficients at cusps from {0} to {1}".format(cuspstart,cuspstop))

    cdef int sym_fak
    V=<double complex****>sage_malloc(sizeof(double complex***)*(numy))
    for yi in range(numy):
        V[yi]=<double complex***>sage_malloc(sizeof(double complex**)*2*nc) #(cuspbb-cuspa))
        if V[yi]==NULL: raise MemoryError
        for i in range(2*nc): #cuspbb-cuspa):
            V[yi][i]=<double complex**>sage_malloc(sizeof(double complex*)*(n_step))
            if V[yi][i]==NULL: raise MemoryError
            for l in range(n_step):
                V[yi][i][l]=<double complex*>sage_malloc(sizeof(double complex)*(N1))
                if V[yi][i][l]==NULL: raise MemoryError
                for n in range(N1):
                    V[yi][i][l][n]=0
    for yi in range(numy):
        #sig_on()
        pullback_pts_cplx_dp_sym(S,Qv,Yv[yi],Xm[yi],Xpb[yi],Ypb[yi],Cvec[yi])
        #sig_off()
        maass_logger.debug("computing first V{0}".format(yi))
        #sig_on()
        compute_V_cplx_dp_sym_for_phase2_sym(V[yi],N1,Xm[yi],Xpb[yi],Ypb[yi],
                                         Cvec[yi],
                                         cusp_evs,alphas,Mv,Qv,Qfak,
                                         symmetric_cusps,
                                         R,Yv[yi],nc,n_a,n_b,
                                         cuspa,cuspbb,1,numy,
                                         verbose-1,0)
        #sig_off()
        maass_logger.debug("after computing first V{0}".format(yi))
    besprec=1.0E-14
    cdef double besmin = 0
    cdef int ncnt=0,redov=0,fi
    nr = <double**>sage_malloc(sizeof(double*)*numy)
    for yi in range(numy):
        nr[yi] = <double*>sage_malloc(sizeof(double)*nc*2)
    kbes = <double**>sage_malloc(sizeof(double)*numy)
    for yi in range(numy):
        kbes[yi] = <double*>sage_malloc(sizeof(double)*nc*2)
    sqrtY=<double*>sage_malloc(sizeof(double)*numy)
    for yn in range(ynmax):
        try:
            try:
                sig_on()
                Q=max(get_M_for_maass_dp_c(R,Y0,eps)+5,M0+10)+Qp
                sig_off()
            except ArithmeticError:
                return {}
            for yi in range(numy):
                sqrtY[yi]=sqrt(Yv[yi])
                Y2pi[yi]=Yv[yi]*twopi
            if verbose>0:
                maass_logger.debug("Q(Y)={0}".format(Q))
                maass_logger.debug("Y[0]={0}".format(Yv[0]))
                maass_logger.debug("Y2pi={0}".format(Y2pi[0]))
                if numy>1:
                    maass_logger.debug("Y[1]={0}".format(Yv[1]))
                    maass_logger.debug("Y2pi[1]={0}".format(Y2pi[1]))
                maass_logger.debug ("NA,NB={0}".format((NAa,NB)))
                maass_logger.debug("cuspa,cuspb,cuspbb={0},{1},{2}".format(cuspa,cuspb,cuspbb))
            #kbdict=_setup_kbessel_dict(Ms,Mf,Qs,Qf,nc,IR,pb,lvec,mp_ctx)
            for n in range(NAa,NB+1):
                ## First check if the current Y is good with respect to the Bessel-factor
                maass_logger.debug("n={0} ncnt={1}".format(n,ncnt))
                decrease_y=0; redov=0
                for yi in range(numy):
                    for jcusp in range(cuspa,cuspbb):
                        nr[yi][jcusp]=fabs(<double>n+alphas[jcusp])*Y2pi[yi]
                        besselk_dp_c(&kbes[yi][jcusp],R,nr[yi][jcusp],besprec,1)
                        kbes[yi][jcusp]=sqrtY[yi]*kbes[yi][jcusp]
                        if Mv[jcusp][0]<0:
                            if alphas[jcusp]==0.0:
                                kbes[yi][jcusp+nc]=kbes[yi][jcusp]
                            else:
                                nr[yi][jcusp+nc]=fabs(<double>-n+alphas[jcusp])*Y2pi[yi]
                                besselk_dp_c(&kbes[yi][jcusp+nc],R,nr[yi][jcusp+nc],besprec,1)
                                kbes[yi][jcusp+nc]=sqrtY[yi]*kbes[yi][jcusp+nc]
                    # Compute the V matrix elements
                besmin=1
                for jcusp in range(cuspa,cuspbb):
                    if fabs(kbes[yi][0])<besmin:
                        besmin = fabs(kbes[yi][0])
                if besmin<eps:
                    if verbose>0:
                        for jcusp in range(cuspa,cuspbb):
                            maass_logger.debug("K_iR(2piY(n+a(0)))={0}".format(kbes[yi][jcusp]))
                    decrease_y=1
                elif ncnt<n_step:
                    for yi in range(numy):
                        for fi in range(fstop-fstart):
                            for jcusp in range(2*nc): #cuspa,cuspbb):
                                Ctmp[yi][fi][jcusp]=0.0
                    for yi in range(numy):
                        phase2_coefficient_sum_cplx_dp_sym(Ctmp[yi],V[yi],Cold,Mv,nc,cusp_evs,
                                                           cusp_offsets,ncnt,cuspa,cuspbb,fstart,fstop,numy,verbose+1)
                        if verbose>0:
                            for jcusp in range(2*nc):
                                maass_logger.debug("Ctmp[{0}][{1}][{2}]={3}".format(yi,0,jcusp,Ctmp[yi][0][jcusp]))
                    for yi in range(numy):
                        for fi in range(fstop-fstart):
                            for jcusp in range(cuspa,cuspbb):
                                if verbose>1:
                                    maass_logger.debug("Ctmp[{0}][{1}][{2}]={3}".format(yi,fi,jcusp,Ctmp[yi][fi][jcusp]))
                                Ctmp[yi][fi][jcusp]=Ctmp[yi][fi][jcusp]/kbes[yi][jcusp]

                                if verbose>1:
                                    maass_logger.debug("sqrt({0})K_i{1}({2})".format(sqrtY[yi],R,nr[yi][jcusp]))
                                    maass_logger.debug("K_iR(2piY(n+a(0)))={0}".format(kbes[yi][jcusp]))

                                maass_logger.debug("C{0}tmp{1}/K_iR[{2}]={3}".format(fi,jcusp,R,Ctmp[yi][fi][jcusp]))
                                if Mv[jcusp][0]<0:
                                    Ctmp[yi][fi][jcusp+nc]=Ctmp[yi][fi][jcusp+nc]/kbes[yi][jcusp+nc]
                    # We check the error for all functions
                    diff=1; diffmax=0
                    for fi in range(fstop-fstart):
                        if method=='TwoY':
                            diff=cabs(Ctmp[1][fi][0]-Ctmp[0][fi][0])
                        else:
                            if cusp_evs[1]<>0:
                                diff=cabs(Ctmp[0][fi][1]-cusp_evs[1]*Ctmp[0][fi][0])
                            if diff>eps:
                                diff=fabs(cabs(Ctmp[0][fi][1])-cabs(Ctmp[0][fi][0]))

                            maass_logger.debug("err[diff][{0}]={1}".format(n,diff))
                            if fabs(cabs(Ctmp[0][fi][1])+cabs(Ctmp[0][fi][0]))< 2*eps and n<M0-10:
                                ## This is to prevent Ctmp[1] and Ctmp[2]
                                ## to be simultaneously close to zero by accident
                                if diff<eps:
                                    if abs(Ctmp[0][fi][0]-XS[fi][0][n])>diff*1.0E6:
                                        diff=max(diff,abs(Ctmp[0][fi][0]-XS[fi][0][n]))

                                        maass_logger.debug("err[Cold][{0},{1}]={2}".format(fi,n,diff))
                        if diff>diffmax:
                            diffmax=diff
                    if verbose>0:
                        maass_logger.debug("diffmax={0}".format(diff))
                        maass_logger.debug("eps={0}".format(eps))
                    if diffmax<eps:
                        for fi in range(fstop-fstart):
                            for jcusp in range(cuspa,cuspbb):
                                Cnew[fi][jcusp][n]=Ctmp[0][fi][jcusp]
                                if Mv[jcusp][0]<0:
                                    Cnew[fi][jcusp][-n]=Ctmp[0][fi][jcusp+nc]

                            # # improve the coefficients we use
                            if n<=M0:
                                for jcusp in range(cuspa,cuspbb):
                                    maass_logger.debug("Pre: Cold[{0}][{1}][{2}]={3}".format(fi,jcusp,n,Cold[fi][jcusp][n-Mv[jcusp][0]]))
                                    Cold[fi][jcusp][n-Mv[jcusp][0]]=Cnew[fi][jcusp][n]
                                    maass_logger.debug("St: Cold[{0}][{1}][{2}]={3}".format(fi,jcusp,n,Cold[fi][jcusp][n-Mv[jcusp][0]]))
                                    if Mv[jcusp][0]<0:
                                        maass_logger.debug("Pre: Cold[{0}][{1}][{2}]={3}".format(fi,jcusp,-n,Cold[fi][jcusp][-n-Mv[jcusp][0]]))
                                        Cold[fi][jcusp][-n-Mv[jcusp][0]]=Cnew[fi][jcusp][-n]
                                        maass_logger.debug("Set: Cold[{0}][{1}][{2}]={3}".format(fi,jcusp,-n,Cold[fi][jcusp][-n-Mv[jcusp][0]]))
                                    #Cold[fi][1][n]=Cnew[fi][1][n]
                            #raise ArithmeticError
                            #if retf:
                            #    for jcusp from cuspa<=jcusp<cuspb:
                            #        XS[jcusp][n]=Cnew[jcusp][n]
                    else:
                        decrease_y=1
                    ncnt+=1
                else:
                    redov=1
                if decrease_y==1:

                    maass_logger.debug("decreasing Y!")
                    NAa=n; ylp=ylp+1
                    if ylp>4:  # We increase Q more than apriori needed in every second step
                        Qp=Qp+10; ylp=0
                    else:
                        maass_logger.debug("0.9Y0={0}, R={1}, n+n_step={2}, eps={3}".format(Y0*0.9,R,n+n_step,eps))
                        Y0 = get_good_Y_for_n(Y0*0.9,R,n+n_step,eps)
                    try:
                        maass_logger.debug("Get Q={0} for Y0={1}".format(Q,Y0))
                        sig_on()
                        Q=max(get_M_for_maass_dp_c(R,Y0,eps)+5,M0+10)+Qp
                        sig_off()
                    except ArithmeticError:
                        return {}
                    Yv[0]=Y0
                    if numy>1:
                        Yv[1]=0.995*Y0
                    #Q=max(get_M_for_maass_dp_c(R,Y0,eps)+5,Ml+10)+Qp
                    set_Mv_Qv_symm(S,Mv,Qv,Qfak,symmetric_cusps,cusp_evs,cusp_offsets,&N1,&Ml,&Ql,M0,Q,verbose-2)
                    # If we change Y we also need to recompute the pullback
                    maass_logger.debug("Here: M0,Q={0},{1}".format(M0,Q))
                    if Xm<>NULL:
                        for yi in range(numy):
                            if Xm[yi]<>NULL:
                                sage_free(Xm[yi])
                        sage_free(Xm)
                    Xm=<double**>sage_malloc(sizeof(double*)*numy)
                    if Xm==NULL: raise MemoryError
                    #if verbose>3:
                    maass_logger.debug("After deallocating Xm!")
                    for yi in range(numy):
                        Xm[yi]=NULL
                        Xm[yi]=<double*>sage_malloc(sizeof(double)*Q)                        
                        if Xm[yi] is NULL: raise MemoryError
                        maass_logger.debug("Allocated new Xm of len {0} x {1}".format(numy,Q))
                        for i in range(nc):
                            if Xpb[yi][i]<>NULL:
                                for j in range(nc):
                                    if Xpb[yi][i][j]<>NULL:
                                        sage_free(Xpb[yi][i][j])
                            if Ypb[yi][i]<>NULL:
                                for j in range(nc):
                                    if Ypb[yi][i][j]<>NULL:
                                        sage_free(Ypb[yi][i][j])
                            if Cvec[yi][i]<>NULL:
                                for j in range(nc):
                                    if Cvec[yi][i][j]<>NULL:
                                        sage_free(Cvec[yi][i][j])
                        if verbose>2:
                            maass_logger.debug("After Xm! {0}".format(yi))
                        for i in range(nc):
                            for j in range(nc):
                                if verbose>3:
                                    maass_logger.debug("Before {0}:{1}:{2}".format(yi,i,j))
                                Xpb[yi][i][j]=NULL; Ypb[yi][i][j]=NULL; Cvec[yi][i][j]=NULL
                                Xpb[yi][i][j] = <double*>sage_malloc(sizeof(double) * Ql )
                                if not Xpb[yi][i][j]: raise MemoryError
                                Ypb[yi][i][j] = <double*>sage_malloc(sizeof(double) * Ql )
                                if not Xpb[yi][i][j]: raise MemoryError
                                Cvec[yi][i][j] = <double complex*>sage_malloc(sizeof(double complex) * Ql )
                                if not Cvec[yi][i][j]: raise MemoryError
                                if verbose>3:
                                    maass_logger.debug("Before assigning! {0}:{1}".format(yi,j))
                                for l in range(Ql):
                                    Xpb[yi][i][j][l]=<double>0
                                    Ypb[yi][i][j][l]=<double>0
                                    Cvec[yi][i][j][l]=0
                                if verbose>3:
                                    maass_logger.debug("After assigning! {0}:{1}".format(yi,j))

                    #if verbose>1:
                    maass_logger.debug("Allocated Xpb Ypb!")
                    #sig_on()
                    # We need to recompute everything
                    for yi in range(numy):
                        if verbose>1:
                            maass_logger.debug("pullback new {0}".format(yi))
                        pullback_pts_cplx_dp_sym(S,Qv,Yv[yi],Xm[yi],Xpb[yi],Ypb[yi],Cvec[yi])
                    redov=1
                    maass_logger.debug("after pullback!")
                else:
                    if verbose>1:
                        maass_logger.debug("Continuing! redov={0}".format(redov))
                if redov==1:
                    NAa=n
                    #if verbose>1:
                    maass_logger.debug("Recompute V!")
                    if n>M0 and cuspb>1:  ## We don't need to compute expansions for all cusps anymore
                        cuspa=cuspstart; cuspb=cuspstop; cuspbb=max(cuspb,2)
                    for yi in range(numy):
                        for i in range(2*nc): #cuspbb-cuspa):
                            for l in range(n_step):
                                for j in range(N1):
                                    V[yi][i][l][j]=0
                        #sig_on()
                        compute_V_cplx_dp_sym_for_phase2_sym(V[yi],N1,Xm[yi],Xpb[yi],Ypb[yi],Cvec[yi],
                                                             cusp_evs,alphas,Mv,Qv,Qfak,
                                                             symmetric_cusps,
                                                             R,Yv[yi],nc,n,n+n_step,cuspa,cuspbb,1,                                                             
#Problem:                                                         R,Yv[yi],nc,n-Mv[0][0],n-Mv[0][0]+n_step,cuspa,cuspbb,1,
                                                         numy,verbose-1,0)
                        #sig_off()
                    redov=0
                    ncnt=0
                    maass_logger.debug("Recomputed V!")
                    raise StopIteration()
            raise StopIteration()
        except StopIteration:

            maass_logger.debug("Stop iteration: n={0}".format(n))
            if n>=NB:
                break #raise StopIteration()
            else:
                continue
        except ArithmeticError as arithm:
            if verbose>0:
                s = str(arithm)
                s += "\n Caught arithmetic Error for R={0},n={1},Y={2}".format(R,n,Y)
                maass_logger.debug(s)
                #print s
            return {}
    ## Deallocating stuff
    #print "too free!"
    if Yv<>NULL:
        sage_free(Yv)
    if Y2pi<>NULL:
        sage_free(Y2pi)
    if Ctmp<>NULL:
        for yi in range(numy):
            if Ctmp[yi]<>NULL:
                for j in range(fstop-fstart):
                    if Ctmp[yi][j]<>NULL:
                        sage_free(Ctmp[yi][j])
                sage_free(Ctmp[yi])
        sage_free(Ctmp)
    if Qfak<>NULL:
        sage_free(Qfak)
    if Qv<>NULL:
        for i in range(nc):
            if Qv[i]<>NULL:
                sage_free(Qv[i])
        sage_free(Qv)
    if cusp_offsets<>NULL:
        sage_free(cusp_offsets)
    if Cold<>NULL:
        for j in range(fstop-fstart):
            if Cold[j]<>NULL:
                for i in range(nc):
                    if cusp_evs[i]<>0 and i<>0:
                        continue
                    if Cold[j][i]<>NULL:
                        sage_free(Cold[j][i])
                sage_free(Cold[j])
        sage_free(Cold)
    if Xpb<>NULL:
        for yi in range(numy):
            if Xpb[yi]<>NULL:
                for i in range(nc):
                    if Xpb[yi][i]<>NULL:
                        for j in range(nc):
                            if Xpb[yi][i][j]<>NULL:
                                sage_free(Xpb[yi][i][j])
                        sage_free(Xpb[yi][i])
                sage_free(Xpb[yi])
        sage_free(Xpb)
    if Ypb<>NULL:
        for yi in range(numy):
            if Ypb[yi]<>NULL:
                for i in range(nc):
                    if Ypb[yi][i]<>NULL:
                        for j in range(nc):
                            if Ypb[yi][i][j]<>NULL:
                                sage_free(Ypb[yi][i][j])
                        sage_free(Ypb[yi][i])
                sage_free(Ypb[yi])
        sage_free(Ypb)
    if Cvec<>NULL:
        for yi in range(numy):
            if Cvec[yi]<>NULL:
                for i in range(nc):
                    if Cvec[yi][i]<>NULL:
                        for j in range(nc):
                            if Cvec[yi][i][j]<>NULL:
                                sage_free(Cvec[yi][i][j])
                        sage_free(Cvec[yi][i])
                sage_free(Cvec[yi])
        sage_free(Cvec)
    if Xm<>NULL:
        for yi in range(numy):
            if Xm[yi]<>NULL:
                sage_free(Xm[yi])
        sage_free(Xm)
    if nr<>NULL:
        for yi in range(numy):
            if nr[yi]<>NULL:
                sage_free(nr[yi])
        sage_free(nr)
    if nr<>NULL:
        for yi in range(numy):
            if kbes[yi]<>NULL:
                sage_free(kbes[yi])
        sage_free(kbes)
    if V<>NULL:
        for yi in range(numy):
            if V[yi]<>NULL:
                for i in range(2*nc): #cuspbb-cuspa):
                    if V[yi][i]<>NULL:
                        for j in range(n_step):
                            if V[yi][i][j]<>NULL:
                                sage_free(V[yi][i][j])
                        sage_free(V[yi][i])
                sage_free(V[yi])
        sage_free(V)
    if cusp_evs<>NULL:
        sage_free(cusp_evs)
    if Mv<>NULL:
        for i in range(nc):
            if Mv[i]<>NULL:
                sage_free(Mv[i])
        sage_free(Mv)

    sage_free(sqrtY)
    sage_free(alphas)
    sage_free(symmetric_cusps)
    return Cnew



@cython.cdivision(True)
cdef compute_V_cplx_dp_sym_for_phase2_sym(double complex ***V,
                                      int N1,
                                      double *Xm,
                                      double ***Xpb,
                                      double ***Ypb,
                                      double complex ***Cvec,
                                      double complex *cusp_evs,
                                      double *alphas,
                                      int **Mv,int **Qv,double *Qfak,
                                      int *symmetric_cusps,
                                      double R,double Y,
                                      int nc,
                                      int NA,
                                      int NB,
                                      int cuspA, int cuspB,
                                      int cuspidal,
                                      int numy,
                                      int verbose,
                                      int is_trivial=0):


    r"""
    Set up the matrix for the system of equations giving the Fourier coefficients of the Maass waveforms.
    INPUT:

    - ``R``   -- double (eigenvalue)
    - ``Y``   -- double (the height of the sampling horocycle)
    - ``Ms,Mf``  -- int (The coefficients we want to compute are C(n), Ms<=n<=Mf )
    - ``Qs,Qf``  -- int (The sampling points are X_m, Qs<=m<=Qf)
    - ``alphas`` -- [nc] double array (the shifts of the Fourier expansion at each cusp)
    - ``V``   -- [(Mf-Ms)*nc]^2 double complex matrix (allocated)
    - ``Xm``  -- [Qf-Qs] double array (allocated)
    - ``Xpb`` -- nc*nc*[Qf-Qs] double array (allocated)
    - ``Ypb`` -- nc*nc*[Qf-Qs] double array (allocated)
    - ``Cvec`` -- nc*nc*[Qf-Qs] double complex array (allocated)
    - `` cuspidal`` -- int (set to 1 if we compute cuspidal functions, otherwise zero)
    - `` sym_type`` -- int (set to 0/1 if we compute even/odd functions, otherwise -1)
    - ``verbose`` -- int (verbosity of output)

    """
    cdef int l,j,icusp,jcusp,n,ni,lj,Ml,Ql,s,Qs,Qf,Mf,Ms
    cdef double pi,sqrtY,Y2pi,nrY2pi,argm,argpb,twopi,two,kbes,besarg,lr,nr
    cdef double complex ckbes,ctmpV,iargm,twopii,ctmp,ctmpV1
        #cdef double *Qfak=NULL

    if not cuspidal in [0,1]:
        raise ValueError," parameter cuspidal must be 0 or 1"
    if Y<=0:
        raise ValueError," need Y>0! Got Y={0}".format(Y)
    pi=M_PI #<double>RealField(53).pi() #3.141592653589793238
    sqrtY=sqrt(Y)
    two=<double>(2)
    Y2pi=Y*pi*two
    twopi=two*pi
    if verbose>=0:
        maass_logger.debug("in compute Vnl with: R,Y={0},{1}".format(R,Y))
        maass_logger.debug("NA,NB,cuspA,cuspB,verbose={0}".format((NA,NB,cuspA,cuspB,verbose)))
        Ml=0; Ql=0
    for i in range(nc):
        if Mv[i][2]>Ml:
            Ml=Mv[i][2]
        if Qv[i][2]>Ql:
            Ql=Qv[i][2]
        if verbose>0:
            maass_logger.debug("Qv[{0}]={1}".format(i,(Qv[i][0],Qv[i][1],Qv[i][2])))
            maass_logger.debug("Mv[{0}]={1}".format(i,(Mv[i][0],Mv[i][1],Mv[i][2])))
    if verbose>2:
        maass_logger.debug("N1={0}".format(N1))
        maass_logger.debug("Ql={0}".format(Ql))
    ## This is the effective offset at the
    cdef int* cusp_offsets=NULL
    cusp_offsets=<int*>sage_malloc(sizeof(int)*nc)
    if cusp_offsets==NULL: raise MemoryError
    for jcusp from 0 <= jcusp < nc:
        cusp_offsets[jcusp]=0
        for icusp from 0 <= icusp < jcusp:
            if icusp==0 or cusp_evs[icusp]==0:
                cusp_offsets[jcusp]+=Mv[icusp][2]
        if verbose>1:
            maass_logger.debug("cusp_offsets[{0}]={1}".format(jcusp,cusp_offsets[jcusp]))
    cdef int nc_sym=0
    for jcusp from 0 <= jcusp < nc:
        if verbose>2:
            maass_logger.debug("cusp_evs[{0}]={1}".format(jcusp,cusp_evs[jcusp]))
        if jcusp==0 or cusp_evs[jcusp]<>0:
            nc_sym+=1
    cdef double **nvec=NULL
    nvec = <double**>sage_malloc(sizeof(double*)*nc)
    if not nvec: raise MemoryError
    for icusp from 0<=icusp<nc:
        nvec[icusp] = <double*>sage_malloc(sizeof(double)*Ml)
    cdef double complex ****ef1=NULL
    cdef double complex ***ef2=NULL
    ef2 = <double complex***>sage_malloc(sizeof(double complex**)*(cuspB-cuspA))
    if ef2==NULL: raise MemoryError
    cdef int n1
    cdef int iicusp

    maass_logger.debug("here0")
    for icusp in range(cuspB-cuspA):
        iicusp=icusp+cuspA
        if cusp_evs[iicusp]<>0.0 and iicusp>0 and (numy==2 or iicusp>1):
            continue
        if verbose>=0:
            maass_logger.debug("icusp={0}".format(icusp))
        if Mv[icusp][0]>=0:
            if verbose>=0:
                maass_logger.debug("test1")
            ef2[icusp]=NULL
            ef2[icusp] = <double complex**>sage_malloc(sizeof(double complex*)*(NB-NA))
            if ef2[icusp]==NULL: raise MemoryError
            for n in range(NB-NA):
                #maass_logger.debug("n={0}".format(n))
                ef2[icusp][n]=NULL
                ef2[icusp][n] = <double complex*>sage_malloc(sizeof(double complex)*Qv[icusp][2])
                if ef2[icusp][n]==NULL: raise MemoryError
        else:
            if verbose>=0:
                maass_logger.debug("test2")
            ef2[icusp]=NULL
            ef2[icusp] = <double complex**>sage_malloc(sizeof(double complex*)*2*(NB-NA))
            if ef2[icusp]==NULL: raise MemoryError
            for n in range(NB-NA):
                #maass_logger.debug("n={0}".format(n))
                n1 = n + NB-NA
                ef2[icusp][n]=NULL;  ef2[icusp][n1]=NULL
                ef2[icusp][n] = <double complex*>sage_malloc(sizeof(double complex)*Qv[icusp][2])
                if ef2[icusp][n]==NULL: raise MemoryError
                ef2[icusp][n1] = <double complex*>sage_malloc(sizeof(double complex)*Qv[icusp][2])
                if ef2[icusp][n1]==NULL: raise MemoryError
    if verbose>=0:
        maass_logger.debug("here1")
    ef1 = <double complex****>sage_malloc(sizeof(double complex***)*nc)
    if ef1==NULL: raise MemoryError
    for icusp in range(nc):
        ef1[icusp]=NULL
        ef1[icusp] = <double complex***>sage_malloc(sizeof(double complex**)*nc)
        if ef1[icusp]==NULL: raise MemoryError
        for jcusp in range(nc):
            ef1[icusp][jcusp]=NULL
            ef1[icusp][jcusp] = <double complex**>sage_malloc(sizeof(double complex*)*Mv[jcusp][2])
            if ef1[icusp][jcusp]==NULL: raise MemoryError
            for n in range(Mv[jcusp][2]):
                ef1[icusp][jcusp][n]=NULL
                ef1[icusp][jcusp][n] = <double complex*>sage_malloc(sizeof(double complex)*Qv[jcusp][2])
                if ef1[icusp][jcusp][n]==NULL: raise MemoryError
    for jcusp in range(nc):
        
        for n in range(Mv[jcusp][2]):
            nvec[jcusp][n]=<double>(n+Mv[jcusp][0])+alphas[jcusp]
    cdef double argpb1,xmm
    if verbose>=0:
        maass_logger.debug("here")
        for jcusp in range(nc):
            maass_logger.debug("symmetric_cusp[{0}]={1}".format(jcusp,symmetric_cusps[jcusp]))
    for jcusp in range(cuspB-cuspA):
        iicusp=jcusp+cuspA
        maass_logger.debug("iicusp={0}".format(iicusp))                    
        if cusp_evs[iicusp]<>0.0 and iicusp>0 and (numy==2 or iicusp>1):
            continue
        for n in range(NB-NA):
            #maass_logger.debug("n={0}".format(n))            
            nr = <double>(n+NA)+alphas[jcusp]
            if symmetric_cusps[jcusp] not in [0,1]:
                for j in range(Qv[iicusp][2]):
                    if j < Qv[0][1]:
                        #argm=nr*Xm[j]
                        xmm = Xm[j]
                        #xmm = -Xm[Qv[0][1]-1-j]
                    else:
                        #argm=-nr*Xm[2*Qv[0][1]-1-j]
                        xmm  = -Xm[2*Qv[0][1]-1-j]
                        #xmm  = Xm[j-Qv[0][1]]
                    argm = nr*xmm
                    ef2[jcusp][n][j]=cexpi(-argm)
                    if n==0 and jcusp==0 and (j==0 or j==Qv[iicusp][2]):
                        maass_logger.debug("Xm[{0}]={1}".format(j,xmm))
                        maass_logger.debug("nr={0}".format(nr))
                        maass_logger.debug("argm={0}".format(argm))                                            
                        maass_logger.debug("ef2[{0}][{1}[{2}]={3}".format(jcusp,n,j,ef2[jcusp][n][j]))                        
            else:
                for j in range(Qv[iicusp][2]):
                    xmm = Xm[j]
                    if n==0 and jcusp==0 and (j==0 or j==Qv[iicusp][2]):
                        maass_logger.debug("Xm[{0}]={1}".format(j,xmm))
                    argm = nr*xmm
                    if symmetric_cusps[jcusp]==0:
                        ef2[jcusp][n][j]=cos(argm)
                    elif symmetric_cusps[jcusp]==1:
                        ef2[jcusp][n][j]=_I*sin(-argm)

        if Mv[jcusp][0]<0: # need both c(n) and c(-n)
            for n in range(NB-NA):
                ## We don't simply use symmetry because alphas[jcusp] can be non-zero
                nr = <double>(-n-NA)+alphas[jcusp]
                n1 = n+NB-NA
                if symmetric_cusps[jcusp] not in [0,1]:
                    for j in range(Qv[iicusp][2]):
                        if j < Qv[0][1]:
                            xmm = Xm[j] #Qv[0][1]-1-j]
                        else:
                            xmm  = -Xm[2*Qv[0][1]-1-j]
                        argm = nr*xmm                            
                        ef2[jcusp][n1][j]=cexpi(-argm)
                else:
                    for j in range(Qv[iicusp][2]):
                        xmm  = Xm[j] 
                        argm = nr*xmm
                        if symmetric_cusps[jcusp]==0:
                            ef2[jcusp][n1][j]=cos(argm)
                        elif symmetric_cusps[jcusp]==1:
                            ef2[jcusp][n1][j]=_I*sin(-argm)



    for jcusp from 0 <= jcusp < nc:
        for icusp from 0<=icusp < cuspB-cuspA:
            iicusp=icusp+cuspA
            for n from 0<=n<Mv[jcusp][2]:
                for j from 0<=j<Qv[jcusp][2]: #in range(Qs,Qf+1):
                    if Ypb[iicusp][jcusp][j]==0: #not Xpb.has_key((icusp,jcusp,j):
                        continue
                    argpb=nvec[jcusp][n]*Xpb[iicusp][jcusp][j]
                    if symmetric_cusps[jcusp]==0:
                        ef1[icusp][jcusp][n][j]=cos(argpb)
                    elif symmetric_cusps[jcusp]==1:
                        ef1[icusp][jcusp][n][j]=_I*sin(argpb)
                    else:
                        ef1[icusp][jcusp][n][j]=cexpi(argpb)
                    ctmp = Cvec[iicusp][jcusp][j]
                    ef1[icusp][jcusp][n][j]=ef1[icusp][jcusp][n][j]*ctmp

    cdef double besarg_old=0.0
    cdef double y,kbes_old=1.0
    cdef double ***kbesvec=NULL
    kbesvec=<double***>sage_malloc(sizeof(double**)*nc)
    if kbesvec==NULL:
        raise MemoryError
    for jcusp from 0<=jcusp<nc:
        #print "allocating kbesvec[",jcusp,"] of size:",Mv[jcusp][2]
        kbesvec[jcusp]=<double**>sage_malloc(sizeof(double*)*Ml)
        if kbesvec[jcusp]==NULL:
            raise MemoryError
        for l from 0<=l<Ml:
            kbesvec[jcusp][l]=<double*>sage_malloc(sizeof(double)*Ql) #Qv[jcusp][2])
            if kbesvec[jcusp][l]==NULL:
                raise MemoryError

    cdef double tmpr
    cdef double besprec
    besprec=1.0E-14

    for jcusp from 0<=jcusp<nc:
        for icusp from 0 <= icusp < cuspB-cuspA:
            iicusp=icusp+cuspA
            for j from 0<=j<Qv[jcusp][2]:
                if Ypb[iicusp][jcusp][j]==0:
                    continue
                for l from 0<=l<Mv[jcusp][2]:
                    lr=nvec[jcusp][l]*twopi
                    #Mf = Mv[icusp][1]
                    besarg=fabs(lr)*Ypb[iicusp][jcusp][j]
                    if lr<>0.0:
                        besselk_dp_c(&tmpr,R,besarg,besprec,1)
                        kbesvec[icusp][l][j]=sqrt(Ypb[iicusp][jcusp][j])*tmpr
                        #if verbose>1 and j==1 and l+Mv[jcusp][0]==1:
                        if verbose>8 and icusp==0 and l==0: #and l+Mv[jcusp][0]<=2 and j>Qv[jcusp][2]-1:
                            #                        if verbose>0 and iicusp==2 and l+Mv[jcusp][1]<=1 and j==Qv[jcusp][2]:
                            maass_logger.debug("lr({0})={1}".format(l,lr))
                            maass_logger.debug("Ypb[{0}][{1}][{2}]={3}".format(iicusp,jcusp,j,Ypb[iicusp][jcusp][j]))
                            #maass_logger.debug("kbes1={0}".format(tmpr))
                            maass_logger.debug("besarg={0}".format(besarg))                            
                            maass_logger.debug("kbes2[{0}][{1}[{2}]={3}".format(icusp,l,j,kbesvec[icusp][l][j]))
                    else:
                        kbesvec[icusp][l][j]=<double>1.0
                    if verbose>8 and icusp==0 and abs(lr)<10: #and l+Mv[jcusp][0]<=2 and j>Qv[jcusp][2]-1:
                        maass_logger.debug("lr={0}".format(lr))                        
                        maass_logger.debug("Ypb[{0}][{1}][{2}]={3}".format(iicusp,jcusp,j,Ypb[iicusp][jcusp][j]))
                        maass_logger.debug("besarg={0}".format(besarg))                            
                        maass_logger.debug("kbes1={0}".format(tmpr))                        
                        
                        maass_logger.debug("kbes2[{0}][{1}][{2}]={3}".format(icusp,l,j,kbesvec[icusp][l][j]))
   
    cdef double complex cuspev
    for jcusp in range(nc):
        for l in range(Mv[jcusp][2]):
            lr=nvec[jcusp][l]*twopi
            lj=cusp_offsets[jcusp]+l
            if jcusp>0 and cusp_evs[jcusp]<>0:
                lj0=l; cuspev=cusp_evs[jcusp]
            else:
                lj0=lj;  cuspev=1.0
            if lr==0.0 and cuspidal==1:
                continue
            for j in range(Qv[jcusp][2]):
                for icusp in range(cuspB-cuspA):
                    iicusp=icusp+cuspA
                    ## If we only have one Y we need the second cusp to estimate the error
                    if cusp_evs[iicusp]<>0.0 and iicusp>0 and (numy==2 or iicusp>1):
                    #if cusp_evs[iicusp]<>0 and iicusp>0:
                        continue
                    #if verbose>0 and jcusp==2 and icusp==2 and l<=2:
                    #    print "Ypb[",iicusp,"][",jcusp,"][",j,"]=",Ypb[iicusp][jcusp][j]
                    if Ypb[iicusp][jcusp][j]==0:
                        continue
                    ckbes=kbesvec[icusp][l][j]*ef1[icusp][jcusp][l][j]
                    if ckbes==0.0:
                        continue
                    for n in range(NB-NA):
                        #if NA+n+Mv[iicusp][0]==0 and cuspidal==1:
                        #    continue
                        ctmpV=ckbes*ef2[icusp][n][j]
                        if lj0>N1:
                            raise ArithmeticError,"Index outside!"
                        V[icusp][n][lj0]+=ctmpV*cuspev
                        if Mv[icusp][0]<0:
                            n1 = n+NB-NA
                            ctmpV1=ckbes*ef2[icusp][n1][j]
                            V[icusp+nc][n][lj0]+=ctmpV1*cuspev
#                        if verbose>0 and jcusp==2 and icusp==2 and lj0+Mv[jcusp][0]<=2 and j>Qv[jcusp][2]-1:
                        if verbose>0 and n==0 and (l == 14 or l==16): #icusp==0 and n==0 and lj0==0:
                            #maass_logger.debug("kbes[{0}][{1}][{2}]={3}".format(icusp,l,j,kbesvec[icusp][l][j]))
                            maass_logger.debug("ef2[{0}][{1}][{2}]={3}".format(icusp,n,j,ef2[icusp][n][j]))
                            #maass_logger.debug("V[{0}][{1}][{2}]={3}".format(icusp,n,lj0,V[icusp][n][lj0]))
    for n in range(NB-NA):
        for l in range(Mv[jcusp][2]):
            maass_logger.debug("V[{0}][{1}][{2}]={3}".format(icusp,n,l,V[icusp][n][l]))
#    raise ArithmeticError
    if verbose>0:
        maass_logger.debug("V[0][0][0]={0}".format(V[0][0][0]))
        for jcusp from 0 <= jcusp < nc:
            maass_logger.debug("Qfak[{0}]={1}".format(jcusp,Qfak[jcusp]))

    for icusp in range(cuspB-cuspA):
        for n in range(NB-NA):
            #for 0<=l<N1:
            #V[icusp][n][l]=V[icusp][n][l]/Qfak[icusp+cuspA][0]
            for jcusp in range(nc):
                if jcusp>0 and cusp_evs[jcusp]<>0:
                    continue
                for l in range(Mv[jcusp][2]):
                    lj=cusp_offsets[jcusp]+l
                    if lj>N1: # some extra insurance...
                        maass_logger.warning("ni=icusp-cuspA={0}".format(n))
                        maass_logger.warning("lj={0}+{1}={2}".format(cusp_offsets[jcusp],l,lj))
                        raise ArithmeticError,"Index outside!"
                    V[icusp][n][lj]=V[icusp][n][lj]/Qfak[jcusp]
                    if Mv[icusp][0]<0:
                        V[icusp+nc][n][lj]=V[icusp+nc][n][lj]/Qfak[jcusp]
    if verbose>0:
        maass_logger.debug("V[0][0][0]={0}".format(V[0][0][0]))
    if kbesvec<>NULL:
        for icusp in range(nc):
            if kbesvec[icusp]<>NULL:
                for l in range(Ml):
                    if kbesvec[icusp][l]<>NULL:
                        sage_free(kbesvec[icusp][l])
                sage_free(kbesvec[icusp])
        sage_free(kbesvec)
    if verbose>0:
        maass_logger.debug("deal kbbes1")
    if ef1<>NULL:
        for jcusp in range(nc):
            if ef1[jcusp]<>NULL:
                for icusp in range(nc):
                    if ef1[jcusp][icusp]<>NULL:
                        for n in range(Mv[icusp][2]):
                            if ef1[jcusp][icusp][n]<>NULL:
                                sage_free(ef1[jcusp][icusp][n])
                        sage_free(ef1[jcusp][icusp])
                sage_free(ef1[jcusp])
        sage_free(ef1)
    #if verbose>0:
    #    print "NB-NA=",NB-NA
    if ef2<>NULL:
        for icusp in range(cuspB-cuspA):
            iicusp=icusp+cuspA
            if cusp_evs[iicusp]<>0.0 and iicusp>0 and (numy==2 or iicusp>1):
                continue

            if ef2[icusp]<>NULL:
                #if verbose>0:
                #    print "icusp=",icusp
                for n in range(NB-NA):
                    #if verbose>0:
                    #    print "n=",n
                    if ef2[icusp][n]<>NULL:
                        sage_free(ef2[icusp][n])
                    if Mv[icusp][0]<0:
                        n1 = n + NB - NA
                        if ef2[icusp][n1]<>NULL:
                            sage_free(ef2[icusp][n1])
                sage_free(ef2[icusp])
        sage_free(ef2)

    if nvec<>NULL:
        for icusp in range(nc):
            if nvec[icusp]<>NULL:
                sage_free(nvec[icusp])
        sage_free(nvec)
    if cusp_offsets<>NULL:
        sage_free(cusp_offsets)



cdef phase2_coefficient_sum_cplx_dp_sym(double complex **Cnew, double complex ***V,
                                        double complex ***Cold,int **Mv,int nc,
                                        double complex *cusp_evs,int* cusp_offsets,int n,
                                        int acusp=0,int bcusp=2,int fstart=0,int fstop=1,
                                        int numy=1,int verbose=0): #symmetric_cusps,
    cdef int icusp,jcusp,l,li,fi
    cdef double complex cc,vv
    if verbose>2:
        maass_logger.debug("in phase2 sum of coefficients ncnt={0}, acusp={1}, bcusp={2}".format(n,acusp,bcusp))
        for icusp in range(nc):
            maass_logger.debug("cusp_evs[{0}]={1}".format(icusp,cusp_evs[icusp]))
        #for icusp from 0<=icusp<nc:
        #    print "Mv[",icusp,"]=",Mv[icusp][0],Mv[icusp][1],Mv[icusp][2]
    for fi in range(fstop-fstart):
        for jcusp in range(nc):
            Cnew[fi][jcusp]=0.0
            if cusp_evs[jcusp]<>0.0 and jcusp>0 and (numy==2 or jcusp>1):
                continue
            for icusp in range(nc):
                if cusp_evs[icusp]<>0.0 and icusp>0:
                    continue
                for l in range(Mv[icusp][2]):
                    li=cusp_offsets[icusp]+l
                    vv = V[jcusp][n][li]
                    cc = Cold[fi][icusp][l]
                    #]if verbose>1:
                    maass_logger.debug("Cold[{0}][{1}][{2}]={3}".format(fi,icusp,l,cc))
                    maass_logger.debug("V[{0}][{1}][{2}]\t={3}".format(jcusp+nc,n,li,vv))                
                    Cnew[fi][jcusp] = Cnew[fi][jcusp] + vv*cc
                    maass_logger.debug("Cnew[{0}][{1}]\t={2}".format(fi,icusp,Cnew[fi][jcusp]))
            maass_logger.debug("")
                    
            if Mv[jcusp][0]<0:
                Cnew[fi][jcusp+nc]=0.0
                for icusp in range(nc):
                    if cusp_evs[icusp]<>0.0 and icusp>0:
                        continue
                    for l in range(Mv[icusp][2]):
                        li=cusp_offsets[icusp]+l
                        vv = V[jcusp+nc][n][li]
                        cc = Cold[fi][icusp][l]
                        maass_logger.debug("Cold[{0}][{1}][{2}]={3}".format(fi,icusp,l,cc))
                        maass_logger.debug("V[{0}][{1}][{2}]\t={3}".format(jcusp+nc,n,li,vv))
                        
                        Cnew[fi][jcusp+nc] = Cnew[fi][jcusp+nc] + vv*cc
                        maass_logger.debug("Cnew[{0}][{1}]\t={2}".format(fi,icusp,Cnew[fi][jcusp+nc]))
                        #if verbose>1 and vv==0:
                        #    maass_logger.debug("V[{0}][{1}][{2}]=0".format(jcusp+nc,n,li))
                        #    maass_logger.debug("Cnew[{0}][{1}]+={2}*{3}={4}".format(fi,jcusp+nc,vv,cc,Cnew[fi][jcusp+nc]))

    if verbose>1:
        maass_logger.debug("end of coefficient sum!")


    #     cdef phase2_coefficient_sum_cplx_dp_sym(double complex **Cnew, double complex ***V,
    #                                     double complex ***Cold,int **Mv,int nc,
    #                                     double complex *cusp_evs,int* cusp_offsets,int n,
    #                                     int acusp=0,int bcusp=2,int fstart=0,int fstop=1,
    #                                     int verbose=0): #symmetric_cusps,
    # cdef int icusp,jcusp,l,li,fi
    # if verbose>0:
    #     print "in phase2 sum of coefficients ncnt=",n
    #     #for icusp from 0<=icusp<nc:
    #     #    print "Mv[",icusp,"]=",Mv[icusp][0],Mv[icusp][1],Mv[icusp][2]
    # for fi from 0<=fi<fstop-fstart:
    #     for icusp from acusp<=icusp<bcusp:
    #         Cnew[fi][icusp]=0.0
    # for icusp from 0<= icusp <nc:
    #     if cusp_evs[icusp]<>0.0 and icusp>0:
    #         continue
    #     if verbose>0:
    #         print "icusp=",icusp,"offset=",cusp_offsets[icusp]
    #     for l from 0 <= l <Mv[icusp][2]:
    #         li=cusp_offsets[icusp]+l
    #         for fi from 0<=fi<fstop-fstart:
    #             for jcusp from acusp<=jcusp<bcusp:
    #                 Cnew[fi][jcusp] = Cnew[fi][jcusp] + V[jcusp][n][li] * Cold[fi][icusp][l]
    #             #Cnew[1] = Cnew[1] + V[1][n][li] * Cold[icusp][l]
    #             if verbose>3:
    #                 for jcusp from acusp<=jcusp<bcusp:
    #                     print "V{0}[{1}]*C[{2},{3}]={4}*{5}".format(jcusp,li,icusp,l+Mv[icusp][0],V[jcusp][n][li],Cold[fi][icusp][l])
    #                 #print "V1[",li,"]*C[",icusp,l+Mv[icusp][0],"]=",V[1][n][li],"*",Cold[icusp][l]
    #         #     print "Cnew0[",icusp,l,"]=",Cnew[0]
    # if verbose>0:
    #     print "end of coefficient sum!"

@cython.cdivision(True)
cdef compute_V_cplx_dp_sym_for_phase2(double complex ***V,
                                      int N1,
                                      double *Xm,
                                      double ***Xpb,
                                      double ***Ypb,
                                      double complex ***Cvec,
                                      double complex *cusp_evs,
                                      double *alphas,
                                      int **Mv,int **Qv,double *Qfak,
                                      int *symmetric_cusps,
                                      double R,double Y,
                                      int nc,
                                      int NA,
                                      int NB,
                                      int cuspA, int cuspB,
                                      int cuspidal,
                                      int numy,
                                      int verbose,
                                      int is_trivial=0):


    r"""
    Set up the matrix for the system of equations giving the Fourier coefficients of the Maass waveforms.
    INPUT:

    - ``R``   -- double (eigenvalue)
    - ``Y``   -- double (the height of the sampling horocycle)
    - ``Ms,Mf``  -- int (The coefficients we want to compute are C(n), Ms<=n<=Mf )
    - ``Qs,Qf``  -- int (The sampling points are X_m, Qs<=m<=Qf)
    - ``alphas`` -- [nc] double array (the shifts of the Fourier expansion at each cusp)
    - ``V``   -- [(Mf-Ms)*nc]^2 double complex matrix (allocated)
    - ``Xm``  -- [Qf-Qs] double array (allocated)
    - ``Xpb`` -- nc*nc*[Qf-Qs] double array (allocated)
    - ``Ypb`` -- nc*nc*[Qf-Qs] double array (allocated)
    - ``Cvec`` -- nc*nc*[Qf-Qs] double complex array (allocated)
    - `` cuspidal`` -- int (set to 1 if we compute cuspidal functions, otherwise zero)
    - `` sym_type`` -- int (set to 0/1 if we compute even/odd functions, otherwise -1)
    - ``verbose`` -- int (verbosity of output)

    """
    cdef int l,j,icusp,jcusp,n,ni,lj,Ml,Ql,s,Qs,Qf,Mf,Ms
    cdef double pi,Y2pi,nrY2pi,argm,argpb,twopi,two,kbes,besarg,lr,nr
    cdef double complex ckbes,ctmpV,iargm,twopii,ctmp,ctmpV1
        #cdef double *Qfak=NULL

    if not cuspidal in [0,1]:
        raise ValueError," parameter cuspidal must be 0 or 1"
    if Y<=0:
        raise ValueError," need Y>0! Got Y={0}".format(Y)
    pi=M_PI #<double>RealField(53).pi() #3.141592653589793238
    two=<double>(2)
    Y2pi=Y*pi*two
    twopi=two*pi
    if verbose>=0:
        maass_logger.debug("in compute Vnl with: R,Y={0},{1}".format(R,Y))
        maass_logger.debug("NA,NB,cuspA,cuspB,verbose={0}".format((NA,NB,cuspA,cuspB,verbose)))
        Ml=0; Ql=0
    for i in range(nc):
        if Mv[i][2]>Ml:
            Ml=Mv[i][2]
        if Qv[i][2]>Ql:
            Ql=Qv[i][2]
        if verbose>0:
            maass_logger.debug("Qv[{0}]={1}".format(i,(Qv[i][0],Qv[i][1],Qv[i][2])))
            maass_logger.debug("Mv[{0}]={1}".format(i,(Mv[i][0],Mv[i][1],Mv[i][2])))
    if verbose>2:
        maass_logger.debug("N1={0}".format(N1))
        maass_logger.debug("Ql={0}".format(Ql))
    ## This is the effective offset at the
    cdef int* cusp_offsets=NULL
    cusp_offsets=<int*>sage_malloc(sizeof(int)*nc)
    if cusp_offsets==NULL: raise MemoryError
    for jcusp from 0 <= jcusp < nc:
        cusp_offsets[jcusp]=0
        for icusp from 0 <= icusp < jcusp:
            if icusp==0 or cusp_evs[icusp]==0:
                cusp_offsets[jcusp]+=Mv[icusp][2]
        if verbose>1:
            maass_logger.debug("cusp_offsets[{0}]={1}".format(jcusp,cusp_offsets[jcusp]))
    cdef int nc_sym=0
    for jcusp from 0 <= jcusp < nc:
        if verbose>2:
            maass_logger.debug("cusp_evs[{0}]={1}".format(jcusp,cusp_evs[jcusp]))
        if jcusp==0 or cusp_evs[jcusp]<>0:
            nc_sym+=1
    cdef double **nvec=NULL
    nvec = <double**>sage_malloc(sizeof(double*)*nc)
    if not nvec: raise MemoryError
    for icusp from 0<=icusp<nc:
        nvec[icusp] = <double*>sage_malloc(sizeof(double)*Ml)
    cdef double complex ****ef1=NULL
    cdef double complex ***ef2=NULL
    ef2 = <double complex***>sage_malloc(sizeof(double complex**)*(cuspB-cuspA))
    if ef2==NULL: raise MemoryError
    cdef int n1
    cdef int iicusp

    maass_logger.debug("here0")
    for icusp in range(cuspB-cuspA):
        iicusp=icusp+cuspA
        if cusp_evs[iicusp]<>0.0 and iicusp>0 and (numy==2 or iicusp>1):
            continue
        if verbose>=0:
            maass_logger.debug("icusp={0}".format(icusp))
        if Mv[icusp][0]>=0:
            if verbose>=0:
                maass_logger.debug("test1")
            ef2[icusp]=NULL
            ef2[icusp] = <double complex**>sage_malloc(sizeof(double complex*)*(NB-NA))
            if ef2[icusp]==NULL: raise MemoryError
            for n in range(NB-NA):
                #maass_logger.debug("n={0}".format(n))
                ef2[icusp][n]=NULL
                ef2[icusp][n] = <double complex*>sage_malloc(sizeof(double complex)*Qv[icusp][2])
                if ef2[icusp][n]==NULL: raise MemoryError
        else:
            if verbose>=0:
                maass_logger.debug("test2")
            ef2[icusp]=NULL
            ef2[icusp] = <double complex**>sage_malloc(sizeof(double complex*)*2*(NB-NA))
            if ef2[icusp]==NULL: raise MemoryError
            for n in range(NB-NA):
                #maass_logger.debug("n={0}".format(n))
                n1 = n + NB-NA
                ef2[icusp][n]=NULL;  ef2[icusp][n1]=NULL
                ef2[icusp][n] = <double complex*>sage_malloc(sizeof(double complex)*Qv[icusp][2])
                if ef2[icusp][n]==NULL: raise MemoryError
                ef2[icusp][n1] = <double complex*>sage_malloc(sizeof(double complex)*Qv[icusp][2])
                if ef2[icusp][n1]==NULL: raise MemoryError
    if verbose>=0:
        maass_logger.debug("here1")
    ef1 = <double complex****>sage_malloc(sizeof(double complex***)*nc)
    if ef1==NULL: raise MemoryError
    for icusp in range(nc):
        ef1[icusp]=NULL
        ef1[icusp] = <double complex***>sage_malloc(sizeof(double complex**)*nc)
        if ef1[icusp]==NULL: raise MemoryError
        for jcusp in range(nc):
            ef1[icusp][jcusp]=NULL
            ef1[icusp][jcusp] = <double complex**>sage_malloc(sizeof(double complex*)*Mv[jcusp][2])
            if ef1[icusp][jcusp]==NULL: raise MemoryError
            for n in range(Mv[jcusp][2]):
                ef1[icusp][jcusp][n]=NULL
                ef1[icusp][jcusp][n] = <double complex*>sage_malloc(sizeof(double complex)*Qv[jcusp][2])
                if ef1[icusp][jcusp][n]==NULL: raise MemoryError
    for jcusp in range(nc):
        for n in range(Mv[jcusp][2]):
            nvec[jcusp][n]=<double>(n+Mv[jcusp][0])+alphas[jcusp]
    cdef double argpb1
    if verbose>=0:
        maass_logger.debug("here")
    for jcusp in range(cuspB-cuspA):
        iicusp=jcusp+cuspA
        if cusp_evs[iicusp]<>0.0 and iicusp>0 and (numy==2 or iicusp>1):
            continue
        for n in range(NB-NA):
            nr = <double>(n+NA)+alphas[jcusp]
            argm=nr*Xm[j]
            for j in range(Qv[iicusp][2]):
                argm=nr*Xm[j]
                if symmetric_cusps[jcusp]==0:
                    ef2[jcusp][n][j]=cos(argm)
                elif symmetric_cusps[jcusp]==1:
                    ef2[jcusp][n][j]=_I*sin(-argm)
                else:
                    ef2[jcusp][n][j]=cexpi(-argm)                    
        if Mv[jcusp][0]<0:
            for n in range(NB-NA):
                ## We don't simply use symmetry because alphas[jcusp] can be non-zero
                nr = <double>(-n-NA)+alphas[jcusp]
                n1 = n+NB-NA
                argm=nr*Xm[j]
                for j in range(Qv[iicusp][2]):                        
                    argm=nr*Xm[j]
                    if symmetric_cusps[jcusp]==0:
                        ef2[jcusp][n1][j]=cos(argm)
                    elif symmetric_cusps[jcusp]==1:
                        ef2[jcusp][n1][j]=_I*sin(-argm)
                    else:
                        ef2[jcusp][n1][j]=cexpi(-argm)                            
                        


    for jcusp from 0 <= jcusp < nc:
        for icusp from 0<=icusp < cuspB-cuspA:
            iicusp=icusp+cuspA
            for n from 0<=n<Mv[jcusp][2]:
                for j from 0<=j<Qv[jcusp][2]: #in range(Qs,Qf+1):
                    if Ypb[iicusp][jcusp][j]==0: #not Xpb.has_key((icusp,jcusp,j):
                        continue
                    argpb=nvec[jcusp][n]*Xpb[iicusp][jcusp][j]
                    if symmetric_cusps[jcusp]==0:
                        ef1[icusp][jcusp][n][j]=cos(argpb)
                    elif symmetric_cusps[jcusp]==1:
                        ef1[icusp][jcusp][n][j]=_I*sin(argpb)
                    else:
                        ef1[icusp][jcusp][n][j]=cexpi(argpb)
                    ctmp = Cvec[iicusp][jcusp][j]
                    ef1[icusp][jcusp][n][j]=ef1[icusp][jcusp][n][j]*ctmp

    cdef double besarg_old=0.0
    cdef double y,kbes_old=1.0
    cdef double ***kbesvec=NULL
    kbesvec=<double***>sage_malloc(sizeof(double**)*nc)
    if kbesvec==NULL:
        raise MemoryError
    for jcusp from 0<=jcusp<nc:
        #print "allocating kbesvec[",jcusp,"] of size:",Mv[jcusp][2]
        kbesvec[jcusp]=<double**>sage_malloc(sizeof(double*)*Ml)
        if kbesvec[jcusp]==NULL:
            raise MemoryError
        for l from 0<=l<Ml:
            kbesvec[jcusp][l]=<double*>sage_malloc(sizeof(double)*Ql) #Qv[jcusp][2])
            if kbesvec[jcusp][l]==NULL:
                raise MemoryError

    cdef double tmpr
    cdef double besprec
    besprec=1.0E-14

    for jcusp from 0<=jcusp<nc:
        for icusp from 0 <= icusp < cuspB-cuspA:
            iicusp=icusp+cuspA
            for j from 0<=j<Qv[jcusp][2]:
                if Ypb[iicusp][jcusp][j]==0:
                    continue
                for l from 0<=l<Mv[jcusp][2]:
                    lr=nvec[jcusp][l]*twopi
                    #Mf = Mv[icusp][1]
                    besarg=fabs(lr)*Ypb[iicusp][jcusp][j]
                    if lr<>0.0:
                        besselk_dp_c(&tmpr,R,besarg,besprec,1)
                        kbesvec[icusp][l][j]=sqrt(Ypb[iicusp][jcusp][j])*tmpr
                        #if verbose>1 and j==1 and l+Mv[jcusp][0]==1:
                        if verbose>0 and jcusp==2 and icusp==2 and l+Mv[jcusp][0]<=2 and j>Qv[jcusp][2]-1:
                            #                        if verbose>0 and iicusp==2 and l+Mv[jcusp][1]<=1 and j==Qv[jcusp][2]:
                            maass_logger.debug("lr({0})={1}".format(l,lr))
                            maass_logger.debug("Ypb[{0}][{1}][{2}]={3}".format(iicusp,jcusp,j,Ypb[iicusp][jcusp][j]))
                            maass_logger.debug("kbes1={0}".format(tmpr))
                            maass_logger.debug("kbes2[{0}][{1}[{2}]={3}".format(icusp,l,j,kbesvec[icusp][l][j]))
                    else:
                        kbesvec[icusp][l][j]=<double>1.0
    cdef double complex cuspev
    for jcusp in range(nc):
        for l in range(Mv[jcusp][2]):
            lr=nvec[jcusp][l]*twopi
            lj=cusp_offsets[jcusp]+l
            if jcusp>0 and cusp_evs[jcusp]<>0:
                lj0=l; cuspev=cusp_evs[jcusp]
            else:
                lj0=lj;  cuspev=1.0
            if lr==0.0 and cuspidal==1:
                continue
            for j in range(Qv[jcusp][2]):
                for icusp in range(cuspB-cuspA):
                    iicusp=icusp+cuspA
                    ## If we only have one Y we need the second cusp to estimate the error
                    if cusp_evs[iicusp]<>0.0 and iicusp>0 and (numy==2 or iicusp>1):
                    #if cusp_evs[iicusp]<>0 and iicusp>0:
                        continue
                    #if verbose>0 and jcusp==2 and icusp==2 and l<=2:
                    #    print "Ypb[",iicusp,"][",jcusp,"][",j,"]=",Ypb[iicusp][jcusp][j]
                    if Ypb[iicusp][jcusp][j]==0:
                        continue
                    ckbes=kbesvec[icusp][l][j]*ef1[icusp][jcusp][l][j]
                    for n in range(NB-NA):
                        #if NA+n+Mv[iicusp][0]==0 and cuspidal==1:
                        #    continue
                        ctmpV=ckbes*ef2[icusp][n][j]
                        if lj0>N1:
                            raise ArithmeticError,"Index outside!"
                        V[icusp][n][lj0]+=ctmpV*cuspev
                        if Mv[icusp][0]<0:
                            n1 = n+NB-NA
                            ctmpV1=ckbes*ef2[icusp][n1][j]
                            V[icusp+nc][n][lj0]+=ctmpV1*cuspev
#                        if verbose>0 and jcusp==2 and icusp==2 and lj0+Mv[jcusp][0]<=2 and j>Qv[jcusp][2]-1:
                        if verbose>0 and icusp==4 and n==0 and lj0==0:
                            maass_logger.debug("kbes[{0}][{1}][{2}]={3}".format(icusp,l,j,kbesvec[icusp][l][j]))
                            maass_logger.debug("V[{0}][{1}][{2}]={3}".format(icusp,n,lj0,V[icusp][n][lj0]))

    if verbose>0:
        maass_logger.debug("V[0][0][0]={0}".format(V[0][0][0]))
        for jcusp from 0 <= jcusp < nc:
            maass_logger.debug("Qfak[{0}]={1}".format(jcusp,Qfak[jcusp]))

    for icusp in range(cuspB-cuspA):
        for n in range(NB-NA):
            #for 0<=l<N1:
            #V[icusp][n][l]=V[icusp][n][l]/Qfak[icusp+cuspA][0]
            for jcusp in range(nc):
                if jcusp>0 and cusp_evs[jcusp]<>0:
                    continue
                for l in range(Mv[jcusp][2]):
                    lj=cusp_offsets[jcusp]+l
                    if lj>N1: # some extra insurance...
                        maass_logger.warning("ni=icusp-cuspA={0}".format(n))
                        maass_logger.warning("lj={0}+{1}={2}".format(cusp_offsets[jcusp],l,lj))
                        raise ArithmeticError,"Index outside!"
                    V[icusp][n][lj]=V[icusp][n][lj]/Qfak[jcusp]
                    if Mv[icusp][0]<0:
                        V[icusp+nc][n][lj]=V[icusp+nc][n][lj]/Qfak[jcusp]
    if verbose>0:
        maass_logger.debug("V[0][0][0]={0}".format(V[0][0][0]))
    if kbesvec<>NULL:
        for icusp in range(nc):
            if kbesvec[icusp]<>NULL:
                for l in range(Ml):
                    if kbesvec[icusp][l]<>NULL:
                        sage_free(kbesvec[icusp][l])
                sage_free(kbesvec[icusp])
        sage_free(kbesvec)
    if verbose>0:
        maass_logger.debug("deal kbbes1")
    if ef1<>NULL:
        for jcusp in range(nc):
            if ef1[jcusp]<>NULL:
                for icusp in range(nc):
                    if ef1[jcusp][icusp]<>NULL:
                        for n in range(Mv[icusp][2]):
                            if ef1[jcusp][icusp][n]<>NULL:
                                sage_free(ef1[jcusp][icusp][n])
                        sage_free(ef1[jcusp][icusp])
                sage_free(ef1[jcusp])
        sage_free(ef1)
    #if verbose>0:
    #    print "NB-NA=",NB-NA
    if ef2<>NULL:
        for icusp in range(cuspB-cuspA):
            iicusp=icusp+cuspA
            if cusp_evs[iicusp]<>0.0 and iicusp>0 and (numy==2 or iicusp>1):
                continue

            if ef2[icusp]<>NULL:
                #if verbose>0:
                #    print "icusp=",icusp
                for n in range(NB-NA):
                    #if verbose>0:
                    #    print "n=",n
                    if ef2[icusp][n]<>NULL:
                        sage_free(ef2[icusp][n])
                    if Mv[icusp][0]<0:
                        n1 = n + NB - NA
                        if ef2[icusp][n1]<>NULL:
                            sage_free(ef2[icusp][n1])
                sage_free(ef2[icusp])
        sage_free(ef2)

    if nvec<>NULL:
        for icusp in range(nc):
            if nvec[icusp]<>NULL:
                sage_free(nvec[icusp])
        sage_free(nvec)
    if cusp_offsets<>NULL:
        sage_free(cusp_offsets)

        
cpdef get_good_Y_for_n(double Y0,double R,int n,double eps):
    cdef int i
    cdef double arg,kbes,Y
    Y=Y0
    cdef double besprec
    besprec=1.0E-14
    for i in range(1000):
        arg=6.28*Y*n
        besselk_dp_c(&kbes,R,arg,besprec,1)
        if fabs(kbes*sqrt(Y))<eps*0.1:
            Y=Y*0.99
        else:
            break
    return Y




############## Real algorithms

@cython.boundscheck(False)
@cython.cdivision(True)
cpdef phase_2_real_dp_sym(S,double R,int NA,int NB,int M0in=0,int ndig=10,int dim=1,int cuspstart=0,int cuspstop=1,int fnr=-1,dict Cin={},double Yin=0,int verbose=0,int retf=0,int n_step=50,int do_test=1,method='2c',really_run=0):
    #,dict c_evs={}):
    r"""
    Computes more coefficients from the given initial set.

    INPUT:

    - ``S`` -- Space of Maas waveforms
    - ``R`` -- real
    - ``Cin`` -- vector
    - ``NA`` -- integer
    - ``NB`` -- integer
    - ``M0`` -- integer (default None) if set we only use this many coeffficients.
    - ``ndig`` -- integer (default 10) : want new coefficients with this many digits precision.
    - ``ndig`` -- number of digits of precision
    - ``dim`` -- assumed dimension
    - ``cuspstart`` -- Compute coefficients between cusp cuspstart and cuspstop
    - ``cuspstop`` --
    - ``fnr`` -- use function nr. fnr (in case of dimension>1)
    - ``method`` -- '2c' or '2Y' for using two cusps (if possible) or two Y's for error test.

    OUTPUT:

    - ``Cout`` -- dictionary of Fourier coefficients

    EXAMPLE::


    sage: R=RR(9.53369526135355755434423523592877032382125639510725198237579046413534)
    sage: IR=CC(0,R)
    sage: M=MaassWaveForms(Gamma0(1))
    sage: C=Maassform_coeffs(M,R)
    sage: D=phase_2(SL2Z,R,C,2,10)



    """
    if not really_run:
        raise NotImplementedError," This real method is not tested properly! Force use on your own risk!"
    cdef double eps,pi,twopi,Y
    cdef int nc,sym_type,M00,Mf,Ms,Ml,Ql
    cdef int **Mv=NULL
    cdef int **Qv=NULL
    cdef int *symmetric_cusps=NULL
    cdef double **Xm=NULL
    cdef double ****Xpb=NULL
    cdef double ****Ypb=NULL
    cdef double ****Cvec=NULL
        #cdef double *Xm2=NULL,***Xpb2=NULL,***Ypb2=NULL
        #cdef double ***Cvec2=NULL
    cdef double *cusp_evs=NULL
    cdef double complex *cusp_evs_cplx=NULL
    cdef double ***Cold=NULL
    cdef double * Qfak=NULL
    cdef int p,n,i,N1=0,NAa=NA
    cdef int ylp,Qp,Q #,verbose
    cdef dict NN,XS
    cdef double ***Ctmp
    cdef double besprec,diff,diffmax,Y02
    cdef double **nr
    cdef double **kbes
    cdef double *sqrtY
    cdef double *Yv
    cdef double *Y2pi
    cdef int yi,M0
    cdef int numy=1
    G=S.group(); nc=int(G.ncusps()); sym_type=S._sym_type
    if nc==1 or S.cusp_symmetries()[1][0]<>1:
        method='TwoY'
        if method=='TwoY':
            numy=2
        Yv=<double*>sage_malloc(sizeof(double)*numy)
    Y2pi=<double*>sage_malloc(sizeof(double)*numy)
    #cdef double nr0,nr1kbes,kbes0,kbes1
    #verbose=S._verbose
    pi=M_PI
    twopi=pi*2.0
    eps=10.0**-ndig
    cdef int fstart,fstop

    if fnr<0 or fnr>dim:
        fstart=0; fstop=dim
    else:
        fstart=fnr; fstop=fnr+1
    if verbose>0:
        print "method=",method
        print "eps=",eps
    cdef list c_evs=S.cusp_evs()
    if not len(c_evs)==nc:
        c_evs = [0 for i in range(len(c_evs),nc)]
    cusp_evs=<double*>sage_malloc(sizeof(double)*nc)
    for i in range(nc):
        cusp_evs[i]=<double> CC(c_evs[i])
    cusp_evs_cplx=<double complex*>sage_malloc(sizeof(double)*nc)
    for i in range(nc):
        cusp_evs_cplx[i]=<double complex> CC(c_evs[i])
    Ctmp = <double***>sage_malloc(sizeof(double**)*numy)
    for yi in range(numy):
        Ctmp[yi] = <double**>sage_malloc(sizeof(double*)*(fstop-fstart))
        for i in range(fstop-fstart):
            Ctmp[yi][i] = <double*>sage_malloc(sizeof(double)*nc*2)
            for j in range(nc*2):
                Ctmp[yi][i][j]=0
    Qfak=<double*>sage_malloc(sizeof(double)*nc)
    if Qfak==NULL: raise MemoryError
    Mv=<int**>sage_malloc(sizeof(int*)*nc)
    if Mv==NULL: raise MemoryError
    Qv=<int**>sage_malloc(sizeof(int*)*nc)
    if Qv==NULL: raise MemoryError
    symmetric_cusps=<int*>sage_malloc(sizeof(int)*nc)
    if symmetric_cusps==NULL: raise MemoryError
    for i in range(nc):
        Mv[i]=<int*>sage_malloc(sizeof(int)*3)
        Qv[i]=<int*>sage_malloc(sizeof(int)*3)
    cdef int* cusp_offsets
    cusp_offsets=<int*>sage_malloc(sizeof(int)*nc)
    if dim<=0:
        if hasattr(S,"_dim"):
            dim=S._dim
        else:
            dim=1
    NN=S.set_norm(dim)
    if verbose>0:
        print "R,Yin,M0in,eps=",R, Yin, M0in, eps
    M0,Y = get_M_and_Y(R,Yin,M0in,eps,verbose-2)
    cdef int old_v=S._verbose
    if verbose>0:
        print "Y,M0=",Y,M0
        print "fstart,fstop=",fstart,fstop
    Q=M0+10
    S._verbose=0
    if Cin <> {}:
        XS = Cin
    if Cin=={} or M0in <M0:
        F = S.get_element(R,dim=dim,Mset=M0)
        if verbose>1:
            print "dim=",dim
            if dim >1 :
                for j in range(dim):
                    print "F[{0}].test()={1}".format(j,F[j].test())
            else:
                print "F.test()={0}".format(F.test())
        if dim>1:
            if fnr>=0 and fnr<dim:
                XS = F[fnr]._coeffs
            else:
                XS=dict()
                for i in range(fstop-fstart):
                    if isinstance(F,dict):
                        if not F.has_key(i):
                            continue
                    elif isinstance(F,list):
                        if i >= len(F):
                            continue
                    XS[i] = F[i]._coeffs[0]
        else:
            XS = F._coeffs
        #sig_on()
        #XS=get_coeff_fast_cplx_dp_sym(S,R,Y,M0,Q,NN,cusp_ev=c_evs)
        #sig_off()
    else:
        XS=Cin
    ## Make sure that we have a real eigenvalue.
    S._verbose=old_v
    #set_Mv_Qv_symm_real(S,Mv,Qv,Qfak,symmetric_cusps,cusp_evs,cusp_offsets,&N1,&Ml,&Ql,M0,Q,verbose)
    set_Mv_Qv_symm(S,Mv,Qv,Qfak,symmetric_cusps,cusp_evs_cplx,cusp_offsets,&N1,&Ml,&Ql,M0,Q,verbose)
    #if dim>1:
    #    p=S.get_primitive_p()
    #    XS = S.Hecke_eigenfunction_from_coeffs(XS,p,fnr=fnr)
    #else:
    #    XS = {0:XS}
    if verbose>3:
        print "XS=",XS
        print "depth=",dict_depth(XS)
    #print "XSK=",XS.keys()
    if do_test:
        F = Maasswaveform(S,R,C=XS)
        #test = S.test_Hecke_relation(XS[0][0])
        t1=  F.test(method='Hecke',format='float')
        t2 =  F.test(method='pcoeff',format='float')
        test = min(t1,t2)
        if verbose>0:
            print "Hecke test=",t1
            print "P-coeff test=",t2
            print "eps=",eps
    else:
        test=0

    if abs(test)>eps:
        print "Need to improve this eigenvalue, or decrease the required precision!"
        return {}

    cdef dict Cnew=dict()
    for j in XS.keys():
        Cnew[j]={}
        for i in XS[j].keys():
            Cnew[j][i]={}
    cdef int cuspa,cuspb
    if retf or nc==1:
        cuspa=0; cuspb=nc; cuspbb=nc
        #cuspa=0; cuspb=2
    else:
        cuspa=cuspstart; cuspb=cuspstop; cuspbb=max(cuspstop,2)

    ## First add the coefficients we already have
    for j in range(fstop-fstart):
        if not XS.has_key(j+fstart): continue
        for i in range(cuspa,cuspb):
            if not XS[j+fstart].has_key(i): continue
            if not isinstance(XS[j+fstart][i],dict):
                for n in XS[j+fstart][0].keys():
                    Cnew[j][i][n]=XS[j+fstart][i]*XS[j+fstart][0][n]
            else:
                for n in XS[j+fstart][i].keys():
                    Cnew[j][i][n]=XS[j+fstart][i][n]
    #print "CK=",Cnew[2].keys()
    if verbose>0:
        print "after!"
    Cold=<double***>sage_malloc(sizeof(double**)*(fstop-fstart))
    for j in range(fstop-fstart):
        Cold[j]=<double**>sage_malloc(sizeof(double*)*nc)
        for i in range(nc):
            if cusp_evs[i]==0 or i==0:
                Cold[j][i]=<double*>sage_malloc(sizeof(double)*Mv[i][2])
                for n in range(Mv[i][2]):
                    Cold[j][i][n]=<double>CC(XS[j][i][n+Mv[i][0]])
            else:
                Cold[j][i]=<double*>sage_malloc(sizeof(double)*Mv[i][2])
                for n in range(Mv[i][2]):
                    Cold[j][i][n]=cusp_evs[i]*XS[j][0][n+Mv[i][0]]
    # using these we want to compute more
    cdef int nn0,nn1
    if verbose>0:
        print "got initial set of coefficients!"
        print "fstart=",fstart
        print "fstop=",fstop
        if verbose>1:
            for j in range(fstop-fstart):
                for n from Mv[0][0]<=n<Mv[0][1]:
                    print "C[{0}][0][{1}]={2}".format(j,n,XS[j][0][n])
                for i in range(nc):
                    if cusp_evs[i]==0:
                        print "C[{0}][{1}][1]={2}".format(j,i,XS[j][i][1])
                    else:
                        print "Cusp_ev[",i,"]=",cusp_evs[i]
                #for n from 0<=n<max(10,Mv[0][1]):
                #    print "Cold[0][",n+Mv[0][0],"]=",Cold[j][0][n]
                for i in range(nc):
                    if cusp_evs[i]==0:
                        nn0 = max(-10,Mv[i][0])
                        nn1 = min(10,Mv[i][1])
                        for n in range(nn0,nn1):
                            if n>Mv[i][2] or n<0:
                                #print "Cold[{0}][{1}][{2}]=?".format(j,i,n)
                                continue
                            print "Cold[{0}][{1}][{2}]={3}".format(j,i,n,Cold[j][i][n-Mv[i][0]])
                            #print "Cold[",i,"][1]=",Cold[j][i][1]

                    # starting value of Y
    cdef double Y0
    if Yin <=0:
        Y0=twopi/<double>NA
        if Y0>0.5:
            Y0=0.5
    else:
        Y0=Yin
    Y0 = min(Y0,S.group().minimal_height())
    Yv[0]=Y0*0.7
    if numy>1:
        #        Yv[1]=0.85*Y0
        Yv[1]=0.995*Y0
    ylp=0; Qp=0
    Ml=2*M0+1
    sig_on()
    Q=max(get_M_for_maass_dp_c(R,Y0,eps)+5,Ml+10)+Qp
    sig_off()
    Ql=2*Q
    Xm=<double**>sage_malloc(sizeof(double*)*numy)
    if Xm==NULL: raise MemoryError
    for yi in range(numy):
        Xm[yi]=<double*>sage_malloc(sizeof(double)*Ql)
    Xpb = <double****> sage_malloc( sizeof(double***) * numy )
    if Xpb==NULL: raise MemoryError
    Ypb = <double****> sage_malloc( sizeof(double*** ) * numy )
    if Ypb==NULL: raise MemoryError
    Cvec = <double****>sage_malloc(sizeof(double***) * numy )
    if Cvec==NULL: raise MemoryError
    for yi in range(numy):
        Xpb[yi] = <double***> sage_malloc( sizeof(double** ) * nc )
        if Xpb[yi]==NULL: raise MemoryError
        Ypb[yi] = <double***> sage_malloc( sizeof(double** ) * nc )
        if Ypb[yi]==NULL: raise MemoryError
        Cvec[yi] = <double***>sage_malloc(sizeof(double**) * nc )
        if Cvec[yi]==NULL: raise MemoryError
        for i in range(nc):
            Xpb[yi][i]=NULL; Ypb[yi][i]=NULL; Cvec[yi][i]=NULL
            Xpb[yi][i] = <double**>sage_malloc(sizeof(double*) * nc )
            if Xpb[yi][i]==NULL: raise MemoryError
            Ypb[yi][i] = <double**>sage_malloc(sizeof(double*) * nc )
            if Ypb[yi][i]==NULL: raise MemoryError
            Cvec[yi][i] = <double**>sage_malloc(sizeof(double*) * nc )
            if Cvec[yi][i]==NULL: raise MemoryError
            for j in range(nc):
                Xpb[yi][i][j]=NULL; Ypb[yi][i][j]=NULL; Cvec[yi][i][j]=NULL
                Xpb[yi][i][j] = <double*>sage_malloc(sizeof(double) * Ql )
                if Xpb[yi][i][j]==NULL: raise MemoryError
                Ypb[yi][i][j] = <double*>sage_malloc(sizeof(double) * Ql )
                if Ypb[yi][i][j]==NULL: raise MemoryError
                Cvec[yi][i][j] = <double*>sage_malloc(sizeof(double) * Ql )
                if Cvec[yi][i][j]==NULL: raise MemoryError
                for n in range(Ql):
                    Xpb[yi][i][j][n]=<double>0
                    Ypb[yi][i][j][n]=<double>0
                    Cvec[yi][i][j][n]=0
    set_Mv_Qv_symm(S,Mv,Qv,Qfak,symmetric_cusps,cusp_evs_cplx,cusp_offsets,&N1,&Ml,&Ql,M0,Q,verbose)
    if verbose>0:
        print "in phase 2 (real) with eps=",eps
        print "range: NA,NB=",NA,NB
        print "N1=",N1
        print "Ml,Ql=",Ml,Ql
    cdef double *alphas=NULL
    alphas=<double*>sage_malloc(sizeof(double)*nc)
    for i in range(nc):
        alphas[i]=<double>S.alpha(i)[0]
    cdef double **** V=NULL
    #cdef int n_step=50
    cdef int n_a,n_b
    n_a = NA - Mv[0][0]
    n_b = n_a + n_step
    if verbose>0:
        print "Compute coefficients at cusps from {0} to {1}".format(cuspstart,cuspstop)

    cdef int sym_fak
    V=<double****>sage_malloc(sizeof(double***)*(numy))
    for yi in range(numy):
        V[yi]=<double***>sage_malloc(sizeof(double**)*2*nc) #(cuspbb-cuspa))
        if V[yi]==NULL: raise MemoryError
        for i in range(2*nc): #cuspbb-cuspa):
            V[yi][i]=<double**>sage_malloc(sizeof(double*)*(n_step))
            if V[yi][i]==NULL: raise MemoryError
            for l in range(n_step):
                V[yi][i][l]=<double*>sage_malloc(sizeof(double)*(N1))
                if V[yi][i][l]==NULL: raise MemoryError
                for n in range(N1):
                    V[yi][i][l][n]=0
    for yi in range(numy):
        #sig_on()
        pullback_pts_real_dp(S,1-Q,Q,Yv[yi],Xm[yi],Xpb[yi],Ypb[yi],Cvec[yi])
        #sig_off()
        if verbose>0:
            print "computing first V{0}".format(yi)
        #sig_on()
        compute_V_real_dp_sym_for_phase2(V[yi],N1,Xm[yi],Xpb[yi],Ypb[yi],
                                         Cvec[yi],
                                         cusp_evs,alphas,Mv,Qv,Qfak,
                                         symmetric_cusps,
                                         R,Yv[yi],nc,n_a,n_b,
                                         cuspa,cuspbb,1,numy,
                                         verbose-1,0)
        #sig_off()
        if verbose>0:
            print "after computing first V{0}".format(yi)
    besprec=1.0E-14
    cdef double besmin = 0
    cdef int ncnt=0,redov=0,fi
    nr = <double**>sage_malloc(sizeof(double*)*numy)
    for yi in range(numy):
        nr[yi] = <double*>sage_malloc(sizeof(double)*nc*2)
    kbes = <double**>sage_malloc(sizeof(double)*numy)
    for yi in range(numy):
        kbes[yi] = <double*>sage_malloc(sizeof(double)*nc*2)
    sqrtY=<double*>sage_malloc(sizeof(double)*numy)
    for yn in range(1000):
        try:
            try:
                sig_on()
                Q=max(get_M_for_maass_dp_c(R,Y0,eps)+5,M0+10)+Qp
                sig_off()
            except ArithmeticError:
                return {}
            for yi in range(numy):
                sqrtY[yi]=sqrt(Yv[yi])
                Y2pi[yi]=Yv[yi]*twopi
            if verbose>0:
                print "Q(Y)=",Q
                print "Y[0]=",Yv[0]
                print "Y2pi=",Y2pi[0]
                if numy>1:
                    print "Y[1]=",Yv[1]
                    print "Y2pi[1]=",Y2pi[1]
                print "NA,NB=",NAa,NB
                print "cuspa,cuspb,cuspbb=",cuspa,cuspb,cuspbb
            #kbdict=_setup_kbessel_dict(Ms,Mf,Qs,Qf,nc,IR,pb,lvec,mp_ctx)
            for n in range(NAa,NB+1):
                ## First check if the current Y is good with respect to the Bessel-factor
                if verbose>0:
                    print "n=",n," ncnt=",ncnt
                decrease_y=0; redov=0
                for yi in range(numy):
                    for jcusp in range(cuspa,cuspbb):
                        nr[yi][jcusp]=fabs(<double>n+alphas[jcusp])*Y2pi[yi]
                        besselk_dp_c(&kbes[yi][jcusp],R,nr[yi][jcusp],besprec,1)
                        kbes[yi][jcusp]=sqrtY[yi]*kbes[yi][jcusp]
                        if Mv[jcusp][0]<0:
                            if alphas[jcusp]==0.0:
                                kbes[yi][jcusp+nc]=kbes[yi][jcusp]
                            else:
                                nr[yi][jcusp+nc]=fabs(<double>-n+alphas[jcusp])*Y2pi[yi]
                                besselk_dp_c(&kbes[yi][jcusp+nc],R,nr[yi][jcusp+nc],besprec,1)
                                kbes[yi][jcusp+nc]=sqrtY[yi]*kbes[yi][jcusp+nc]
                    # Compute the V matrix elements
                besmin=1
                for jcusp in range(cuspa,cuspbb):
                    if fabs(kbes[yi][0])<besmin:
                        besmin = fabs(kbes[yi][0])
                if besmin<eps:
                    if verbose>0:
                        for jcusp in range(cuspa,cuspbb):
                            print "K_iR(2piY(n+a(0)))=",kbes[yi][jcusp]
                    decrease_y=1
                elif ncnt<n_step:
                    for yi in range(numy):
                        for fi in range(fstop-fstart):
                            for jcusp in range(2*nc): #cuspa,cuspbb):
                                Ctmp[yi][fi][jcusp]=0.0
                    for yi in range(numy):
                        phase2_coefficient_sum_real_dp_sym(Ctmp[yi],V[yi],Cold,Mv,nc,cusp_evs,
                                                           cusp_offsets,ncnt,cuspa,cuspbb,fstart,fstop,numy,verbose+1)
                        if verbose>1:
                            for jcusp in range(2*nc):
                                print "CTMP[",jcusp,"]=",Ctmp[yi][0][jcusp]
                    for yi in range(numy):
                        for fi in range(fstop-fstart):
                            for jcusp in range(cuspa,cuspbb):
                                if verbose>1:
                                    print "Ctmp[{0}][{1}][{2}]={3}".format(yi,fi,jcusp,Ctmp[yi][fi][jcusp])
                                Ctmp[yi][fi][jcusp]=Ctmp[yi][fi][jcusp]/kbes[yi][jcusp]

                                if verbose>1:
                                    print "sqrt({0})K_i{1}({2})".format(sqrtY[yi],R,nr[yi][jcusp])
                                    print "K_iR(2piY(n+a(0)))=",kbes[yi][jcusp]
                                if verbose>0:
                                    print "C{0}tmp{1}/K_iR[{2}]={3}".format(fi,jcusp,R,Ctmp[yi][fi][jcusp])
                                if Mv[jcusp][0]<0:
                                    Ctmp[yi][fi][jcusp+nc]=Ctmp[yi][fi][jcusp+nc]/kbes[yi][jcusp+nc]
                    # We check the error for all functions
                    diff=1; diffmax=0
                    for fi in range(fstop-fstart):
                        if method=='TwoY':
                            diff=cabs(Ctmp[1][fi][0]-Ctmp[0][fi][0])
                        else:
                            if cusp_evs[1]<>0:
                                diff=cabs(Ctmp[0][fi][1]-cusp_evs[1]*Ctmp[0][fi][0])
                            if diff>eps:
                                diff=fabs(cabs(Ctmp[0][fi][1])-cabs(Ctmp[0][fi][0]))
                            if verbose>0:
                                print "err[diff][",n,"]=",diff
                            if fabs(cabs(Ctmp[0][fi][1])+cabs(Ctmp[0][fi][0]))< 2*eps and n<M0-10:
                                ## This is to prevent Ctmp[1] and Ctmp[2]
                                ## to be simultaneously close to zero by accident
                                if diff<eps:
                                    if abs(Ctmp[0][fi][0]-XS[fi][0][n])>diff*1.0E6:
                                        diff=max(diff,abs(Ctmp[0][fi][0]-XS[fi][0][n]))
                                        if verbose>0:
                                            print "err[Cold][",fi,n,"]=",diff
                        if diff>diffmax:
                            diffmax=diff
                    if verbose>0:
                        print "diffmax=",diff
                        print "eps=",eps
                    if diffmax<eps:
                        for fi in range(fstop-fstart):
                            for jcusp in range(cuspa,cuspbb):
                                Cnew[fi][jcusp][n]=Ctmp[0][fi][jcusp]
                                if Mv[jcusp][0]<0:
                                    Cnew[fi][jcusp][-n]=Ctmp[0][fi][jcusp+nc]

                            # # improve the coefficients we use
                            if n<=M0:
                                for jcusp in range(cuspa,cuspbb):
                                    if verbose>0:
                                        print "Pre: Cold[{0}][{1}][{2}]={3}".format(fi,jcusp,n,Cold[fi][jcusp][n-Mv[jcusp][0]])
                                    Cold[fi][jcusp][n-Mv[jcusp][0]]=Cnew[fi][jcusp][n]
                                    if verbose>0:
                                        print "Set: Cold[{0}][{1}][{2}]={3}".format(fi,jcusp,n,Cold[fi][jcusp][n-Mv[jcusp][0]])
                                    if Mv[jcusp][0]<0:
                                        if verbose>0:
                                            print "Pre: Cold[{0}][{1}][{2}]={3}".format(fi,jcusp,-n,Cold[fi][jcusp][-n-Mv[jcusp][0]])
                                        Cold[fi][jcusp][-n-Mv[jcusp][0]]=Cnew[fi][jcusp][-n]
                                        if verbose>0:
                                            print "Set: Cold[{0}][{1}][{2}]={3}".format(fi,jcusp,-n,Cold[fi][jcusp][-n-Mv[jcusp][0]])
                                    #Cold[fi][1][n]=Cnew[fi][1][n]
                            #raise ArithmeticError
                            #if retf:
                            #    for jcusp from cuspa<=jcusp<cuspb:
                            #        XS[jcusp][n]=Cnew[jcusp][n]
                    else:
                        decrease_y=1
                    ncnt+=1
                else:
                    redov=1
                if decrease_y==1:
                    if verbose>0:
                        print "decreasing Y!"
                    NAa=n; ylp=ylp+1
                    if ylp>4:  # We increase Q more than apriori needed in every second step
                        Qp=Qp+10; ylp=0
                    else:
                        Y0 = get_good_Y_for_n(Y0*0.9,R,n+n_step,eps)
                    try:
                        sig_on()
                        Q=max(get_M_for_maass_dp_c(R,Y0,eps)+5,M0+10)+Qp
                        sig_off()
                    except ArithmeticError:
                        return {}
                    Yv[0]=Y0
                    if numy>1:
                        Yv[1]=0.995*Y0
                    #Q=max(get_M_for_maass_dp_c(R,Y0,eps)+5,Ml+10)+Qp
                    set_Mv_Qv_symm(S,Mv,Qv,Qfak,symmetric_cusps,cusp_evs_cplx,cusp_offsets,&N1,&Ml,&Ql,M0,Q,verbose-2)
                    # If we change Y we also need to recompute the pullback
                    if verbose>1:
                        print "Here: M0,Q=",M0,Q
                    if Xm<>NULL:
                        for yi in range(numy):
                            if Xm[yi]<>NULL:
                                sage_free(Xm[yi])
                        sage_free(Xm)
                    Xm=<double**>sage_malloc(sizeof(double*)*numy)
                    if not Xm: raise MemoryError
                    if verbose>3:
                        print "After Xm!1"
                    for yi in range(numy):
                        Xm[yi]=NULL
                        Xm[yi]=<double*>sage_malloc(sizeof(double)*Ql)
                        if not Xm[yi]: raise MemoryError
                        for i in range(nc):
                            if Xpb[yi][i]<>NULL:
                                for j in range(nc):
                                    if Xpb[yi][i][j]<>NULL:
                                        sage_free(Xpb[yi][i][j])
                            if Ypb[yi][i]<>NULL:
                                for j in range(nc):
                                    if Ypb[yi][i][j]<>NULL:
                                        sage_free(Ypb[yi][i][j])
                            if Cvec[yi][i]<>NULL:
                                for j in range(nc):
                                    if Cvec[yi][i][j]<>NULL:
                                        sage_free(Cvec[yi][i][j])
                        if verbose>3:
                            print "After Xm! {0}".format(yi)
                        for i in range(nc):
                            for j in range(nc):
                                if verbose>3:
                                    print "Before {0}:{1}:{2}".format(yi,i,j)
                                Xpb[yi][i][j]=NULL; Ypb[yi][i][j]=NULL; Cvec[yi][i][j]=NULL
                                Xpb[yi][i][j] = <double*>sage_malloc(sizeof(double) * Ql )
                                if not Xpb[yi][i][j]: raise MemoryError
                                Ypb[yi][i][j] = <double*>sage_malloc(sizeof(double) * Ql )
                                if not Xpb[yi][i][j]: raise MemoryError
                                Cvec[yi][i][j] = <double*>sage_malloc(sizeof(double) * Ql )
                                if not Cvec[yi][i][j]: raise MemoryError
                                if verbose>3:
                                    print "Before assigning! {0}:{1}".format(yi,j)
                                for l in range(Ql):
                                    Xpb[yi][i][j][l]=<double>0
                                    Ypb[yi][i][j][l]=<double>0
                                    Cvec[yi][i][j][l]=0
                                if verbose>3:
                                    print "After assigning! {0}:{1}".format(yi,j)

                    if verbose>1:
                        print "Allocated Xpb Ypb!"
                    #sig_on()
                    # We need to recompute everything
                    for yi in range(numy):
                        if verbose>1:
                            print "pullback new {0}".format(yi)
                        pullback_pts_real_dp(S,1-Q,Q,Yv[yi],Xm[yi],Xpb[yi],Ypb[yi],Cvec[yi])
                    redov=1
                else:
                    if verbose>1:
                        print "Continuing! redov=",redov
                if redov==1:
                    NAa=n
                    if verbose>1:
                        print "Recompute V!"
                    if n>M0 and cuspb>1:  ## We don't need to compute expansions for all cusps anymore
                        cuspa=cuspstart; cuspb=cuspstop; cuspbb=max(cuspb,2)
                    for yi in range(numy):
                        for i in range(2*nc): #cuspbb-cuspa):
                            for l in range(n_step):
                                for j in range(N1):
                                    V[yi][i][l][j]=0
                        #sig_on()
                        compute_V_real_dp_sym_for_phase2(V[yi],N1,Xm[yi],Xpb[yi],Ypb[yi],Cvec[yi],
                                                         cusp_evs,alphas,Mv,Qv,Qfak,
                                                         symmetric_cusps,
                                                         R,Yv[yi],nc,n-Mv[0][0],n-Mv[0][0]+n_step,cuspa,cuspbb,1,
                                                         numy,verbose-1,0)
                        #sig_off()
                    redov=0
                    ncnt=0
                    raise StopIteration()
            raise StopIteration()
        except StopIteration:
            if verbose>0:
                print "Stop iteration: n=",n
            if n>=NB:
                break #raise StopIteration()
            else:
                continue
        except ArithmeticError as arithm:
            if verbose>0:
                s = str(arithm)
                s += "\n Caught arithmetic Error for R={0},n={1},Y={2}".format(R,n,Y)
                print s
            return {}
    ## Deallocating stuff
    #print "too free!"
    if Yv<>NULL:
        sage_free(Yv)
    if Y2pi<>NULL:
        sage_free(Y2pi)
    if Ctmp<>NULL:
        for yi in range(numy):
            if Ctmp[yi]<>NULL:
                for j in range(fstop-fstart):
                    if Ctmp[yi][j]<>NULL:
                        sage_free(Ctmp[yi][j])
                sage_free(Ctmp[yi])
        sage_free(Ctmp)
    if Qfak<>NULL:
        sage_free(Qfak)
    if Qv<>NULL:
        for i in range(nc):
            if Qv[i]<>NULL:
                sage_free(Qv[i])
        sage_free(Qv)
    if cusp_offsets<>NULL:
        sage_free(cusp_offsets)
    if Cold<>NULL:
        for j in range(fstop-fstart):
            if Cold[j]<>NULL:
                for i in range(nc):
                    if cusp_evs[i]<>0 and i<>0:
                        continue
                    if Cold[j][i]<>NULL:
                        sage_free(Cold[j][i])
                sage_free(Cold[j])
        sage_free(Cold)
    if Xpb<>NULL:
        for yi in range(numy):
            if Xpb[yi]<>NULL:
                for i in range(nc):
                    if Xpb[yi][i]<>NULL:
                        for j in range(nc):
                            if Xpb[yi][i][j]<>NULL:
                                sage_free(Xpb[yi][i][j])
                        sage_free(Xpb[yi][i])
                sage_free(Xpb[yi])
        sage_free(Xpb)
    if Ypb<>NULL:
        for yi in range(numy):
            if Ypb[yi]<>NULL:
                for i in range(nc):
                    if Ypb[yi][i]<>NULL:
                        for j in range(nc):
                            if Ypb[yi][i][j]<>NULL:
                                sage_free(Ypb[yi][i][j])
                        sage_free(Ypb[yi][i])
                sage_free(Ypb[yi])
        sage_free(Ypb)
    if Cvec<>NULL:
        for yi in range(numy):
            if Cvec[yi]<>NULL:
                for i in range(nc):
                    if Cvec[yi][i]<>NULL:
                        for j in range(nc):
                            if Cvec[yi][i][j]<>NULL:
                                sage_free(Cvec[yi][i][j])
                        sage_free(Cvec[yi][i])
                sage_free(Cvec[yi])
        sage_free(Cvec)
    if Xm<>NULL:
        for yi in range(numy):
            if Xm[yi]<>NULL:
                sage_free(Xm[yi])
        sage_free(Xm)
    if nr<>NULL:
        for yi in range(numy):
            if nr[yi]<>NULL:
                sage_free(nr[yi])
        sage_free(nr)
    if nr<>NULL:
        for yi in range(numy):
            if kbes[yi]<>NULL:
                sage_free(kbes[yi])
        sage_free(kbes)
    if V<>NULL:
        for yi in range(numy):
            if V[yi]<>NULL:
                for i in range(2*nc): #cuspbb-cuspa):
                    if V[yi][i]<>NULL:
                        for j in range(n_step):
                            if V[yi][i][j]<>NULL:
                                sage_free(V[yi][i][j])
                        sage_free(V[yi][i])
                sage_free(V[yi])
        sage_free(V)
    if cusp_evs<>NULL:
        sage_free(cusp_evs)
    if cusp_evs_cplx<>NULL:
        sage_free(cusp_evs_cplx)
    if Mv<>NULL:
        for i in range(nc):
            if Mv[i]<>NULL:
                sage_free(Mv[i])
        sage_free(Mv)

    sage_free(sqrtY)
    sage_free(alphas)
    sage_free(symmetric_cusps)
    return Cnew



@cython.cdivision(True)
cdef compute_V_real_dp_sym_for_phase2(double ***V,
                                      int N1,
                                      double *Xm,
                                      double ***Xpb,
                                      double ***Ypb,
                                      double ***Cvec,
                                      double *cusp_evs,
                                      double *alphas,
                                      int **Mv,int **Qv,double *Qfak,
                                      int *symmetric_cusps,
                                      double R,double Y,
                                      int nc,
                                      int NA,
                                      int NB,
                                      int cuspA, int cuspB,
                                      int cuspidal,
                                      int numy,
                                      int verbose,
                                      int is_trivial=0):


    r"""
    Set up the matrix for the system of equations giving the Fourier coefficients of the Maass waveforms.
    INPUT:

    - ``R``   -- double (eigenvalue)
    - ``Y``   -- double (the height of the sampling horocycle)
    - ``Ms,Mf``  -- int (The coefficients we want to compute are C(n), Ms<=n<=Mf )
    - ``Qs,Qf``  -- int (The sampling points are X_m, Qs<=m<=Qf)
    - ``alphas`` -- [nc] double array (the shifts of the Fourier expansion at each cusp)
    - ``V``   -- [(Mf-Ms)*nc]^2 double matrix (allocated)
    - ``Xm``  -- [Qf-Qs] double array (allocated)
    - ``Xpb`` -- nc*nc*[Qf-Qs] double array (allocated)
    - ``Ypb`` -- nc*nc*[Qf-Qs] double array (allocated)
    - ``Cvec`` -- nc*nc*[Qf-Qs] double array (allocated)
    - `` cuspidal`` -- int (set to 1 if we compute cuspidal functions, otherwise zero)
    - `` sym_type`` -- int (set to 0/1 if we compute even/odd functions, otherwise -1)
    - ``verbose`` -- int (verbosity of output)

    """
    cdef int l,j,icusp,jcusp,n,ni,lj,Ml,Ql,s,Qs,Qf,Mf,Ms
    cdef double pi,sqrtY,Y2pi,nrY2pi,argm,argpb,twopi,two,kbes,besarg,lr,nr
    cdef double ckbes,ctmpV,iargm,twopii,ctmp,ctmpV1
        #cdef double *Qfak=NULL

    if not cuspidal in [0,1]:
        raise ValueError," parameter cuspidal must be 0 or 1"
    if Y<=0:
        raise ValueError," need Y>0! Got Y={0}".format(Y)
    pi=M_PI #<double>RealField(53).pi() #3.141592653589793238
    sqrtY=sqrt(Y)
    two=<double>(2)
    Y2pi=Y*pi*two
    twopi=two*pi
    if verbose>=0:
        print "in compute Vnl with: R,Y",R,Y
        print "NA,NB,cuspA,cuspB,verbose=",NA,NB,cuspA,cuspB,verbose
        Ml=0; Ql=0
    for i from 0<=i<nc:
        if Mv[i][2]>Ml:
            Ml=Mv[i][2]
        if Qv[i][2]>Ql:
            Ql=Qv[i][2]
        if verbose>0:
            print "Qv[",i,"]=(",Qv[i][0],",",Qv[i][1],",",Qv[i][2],")"
            print "Mv[",i,"]=(",Mv[i][0],",",Mv[i][1],",",Mv[i][2],")"
    if verbose>2:
        print "N1=",N1
        print "Ql=",Ql
    ## This is the effective offset at the
    cdef int* cusp_offsets=NULL
    cusp_offsets=<int*>sage_malloc(sizeof(int)*nc)
    if cusp_offsets==NULL: raise MemoryError
    for jcusp from 0 <= jcusp < nc:
        cusp_offsets[jcusp]=0
        for icusp from 0 <= icusp < jcusp:
            if icusp==0 or cusp_evs[icusp]==0:
                cusp_offsets[jcusp]+=Mv[icusp][2]
        if verbose>1:
            print "cusp_offsets[",jcusp,"]=",cusp_offsets[jcusp]
    cdef int nc_sym=0
    for jcusp from 0 <= jcusp < nc:
        if verbose>2:
            print "cusp_evs[",jcusp,"]=",cusp_evs[jcusp]
        if jcusp==0 or cusp_evs[jcusp]<>0:
            nc_sym+=1
    cdef double **nvec=NULL
    nvec = <double**>sage_malloc(sizeof(double*)*nc)
    if not nvec: raise MemoryError
    for icusp from 0<=icusp<nc:
        nvec[icusp] = <double*>sage_malloc(sizeof(double)*Ml)
    cdef double ****ef1=NULL
    cdef double ***ef2=NULL
    ef2 = <double ***>sage_malloc(sizeof(double**)*(cuspB-cuspA))
    if ef2==NULL: raise MemoryError
    cdef int n1
    cdef int iicusp

    for icusp in range(cuspB-cuspA):
        iicusp=icusp+cuspA
        if cusp_evs[iicusp]<>0.0 and iicusp>0 and (numy==2 or iicusp>1):
            continue
        if Mv[icusp][0]>=0:
            ef2[icusp]=NULL
            ef2[icusp] = <double**>sage_malloc(sizeof(double*)*(NB-NA))
            if ef2[icusp]==NULL: raise MemoryError
            for n in range(NB-NA):
                ef2[icusp][n]=NULL
                ef2[icusp][n] = <double*>sage_malloc(sizeof(double)*Qv[icusp][2])
                if ef2[icusp][n]==NULL: raise MemoryError
        else:
            ef2[icusp]=NULL
            ef2[icusp] = <double**>sage_malloc(sizeof(double*)*2*(NB-NA))
            if ef2[icusp]==NULL: raise MemoryError
            for n in range(NB-NA):
                n1 = n + NB-NA
                ef2[icusp][n]=NULL;  ef2[icusp][n1]=NULL
                ef2[icusp][n] = <double*>sage_malloc(sizeof(double)*Qv[icusp][2])
                if ef2[icusp][n]==NULL: raise MemoryError
                ef2[icusp][n1] = <double*>sage_malloc(sizeof(double)*Qv[icusp][2])
                if ef2[icusp][n1]==NULL: raise MemoryError
    ef1 = <double****>sage_malloc(sizeof(double***)*nc)
    if ef1==NULL: raise MemoryError
    for icusp in range(nc):
        ef1[icusp]=NULL
        ef1[icusp] = <double***>sage_malloc(sizeof(double**)*nc)
        if ef1[icusp]==NULL: raise MemoryError
        for jcusp in range(nc):
            ef1[icusp][jcusp]=NULL
            ef1[icusp][jcusp] = <double**>sage_malloc(sizeof(double*)*Mv[jcusp][2])
            if ef1[icusp][jcusp]==NULL: raise MemoryError
            for n in range(Mv[jcusp][2]):
                ef1[icusp][jcusp][n]=NULL
                ef1[icusp][jcusp][n] = <double *>sage_malloc(sizeof(double )*Qv[jcusp][2])
                if ef1[icusp][jcusp][n]==NULL: raise MemoryError
    for jcusp in range(nc):
        for n in range(Mv[jcusp][2]):
            nvec[jcusp][n]=<double>(n+Mv[jcusp][0])+alphas[jcusp]
    cdef double argpb1

    for jcusp in range(cuspB-cuspA):
        iicusp=jcusp+cuspA
        if cusp_evs[iicusp]<>0.0 and iicusp>0 and (numy==2 or iicusp>1):
            continue
        for n in range(NB-NA):
            nr = <double>(n+NA)+alphas[jcusp]
            for j in range(Qv[jcusp+cuspA][2]):
                argm=nr*Xm[j]
                if symmetric_cusps[jcusp]==0:
                    ef2[jcusp][n][j]=cos(argm)
                elif symmetric_cusps[jcusp]==1:
                    ef2[jcusp][n][j]=sin(argm)
                else:
                    raise ValueError,"Need symmetry to be 0 or 1 here!"
        if Mv[jcusp][0]<0:
            for n in range(NB-NA):
                ## We don't simply use symmetry because alphas[jcusp] can be non-zero
                nr = <double>(-n-NA)+alphas[jcusp]
                n1 = n+NB-NA
                for j in range(Qv[jcusp+cuspA][2]):
                    argm=nr*Xm[j]
                    if symmetric_cusps[jcusp]==0:
                        ef2[jcusp][n1][j]=cos(argm)
                    elif symmetric_cusps[jcusp]==1:
                        ef2[jcusp][n1][j]=sin(argm)
                    else:
                        raise ValueError,"Need symmetry to be 0 or 1 here!"


    for jcusp from 0 <= jcusp < nc:
        for icusp from 0<=icusp < cuspB-cuspA:
            iicusp=icusp+cuspA
            for n from 0<=n<Mv[jcusp][2]:
                for j from 0<=j<Qv[jcusp][2]: #in range(Qs,Qf+1):
                    if Ypb[iicusp][jcusp][j]==0: #not Xpb.has_key((icusp,jcusp,j):
                        continue
                    argpb=nvec[jcusp][n]*Xpb[iicusp][jcusp][j]
                    if symmetric_cusps[jcusp]==0:
                        ef1[icusp][jcusp][n][j]=cos(argpb)
                    elif symmetric_cusps[jcusp]==1:
                        ef1[icusp][jcusp][n][j]=sin(argpb)
                    else:
                        raise ValueError,"Need symmetry to be 0 or 1 here!"
                    ctmp = Cvec[iicusp][jcusp][j]
                    ef1[icusp][jcusp][n][j]=ef1[icusp][jcusp][n][j]*ctmp

    cdef double besarg_old=0.0
    cdef double y,kbes_old=1.0
    cdef double ***kbesvec=NULL
    kbesvec=<double***>sage_malloc(sizeof(double**)*nc)
    if kbesvec==NULL:
        raise MemoryError
    for jcusp from 0<=jcusp<nc:
        #print "allocating kbesvec[",jcusp,"] of size:",Mv[jcusp][2]
        kbesvec[jcusp]=<double**>sage_malloc(sizeof(double*)*Ml)
        if kbesvec[jcusp]==NULL:
            raise MemoryError
        for l from 0<=l<Ml:
            kbesvec[jcusp][l]=<double*>sage_malloc(sizeof(double)*Ql) #Qv[jcusp][2])
            if kbesvec[jcusp][l]==NULL:
                raise MemoryError

    cdef double tmpr
    cdef double besprec
    besprec=1.0E-14

    for jcusp from 0<=jcusp<nc:
        for icusp from 0 <= icusp < cuspB-cuspA:
            iicusp=icusp+cuspA
            for j from 0<=j<Qv[jcusp][2]:
                if Ypb[iicusp][jcusp][j]==0:
                    continue
                for l from 0<=l<Mv[jcusp][2]:
                    lr=nvec[jcusp][l]*twopi
                    #Mf = Mv[icusp][1]
                    besarg=fabs(lr)*Ypb[iicusp][jcusp][j]
                    if lr<>0.0:
                        besselk_dp_c(&tmpr,R,besarg,besprec,1)
                        kbesvec[icusp][l][j]=sqrt(Ypb[iicusp][jcusp][j])*tmpr
                        #if verbose>1 and j==1 and l+Mv[jcusp][0]==1:
                        if verbose>0 and jcusp==2 and icusp==2 and l+Mv[jcusp][0]<=2 and j>Qv[jcusp][2]-1:
                        #if verbose>0 and iicusp==2 and l+Mv[jcusp][1]<=1 and j==Qv[jcusp][2]:
                            print "lr(",l,")=",lr
                            print "Ypb[",iicusp,jcusp,j,"=",Ypb[iicusp][jcusp][j]
                            print "kbes1=",tmpr
                            print "kbes2[",icusp,l,j,"]=",kbesvec[icusp][l][j]
                    else:
                        kbesvec[icusp][l][j]=<double>1.0
    cdef double  cuspev
    for jcusp in range(nc):
        for l in range(Mv[jcusp][2]):
            lr=nvec[jcusp][l]*twopi
            lj=cusp_offsets[jcusp]+l
            if jcusp>0 and cusp_evs[jcusp]<>0:
                lj0=l; cuspev=cusp_evs[jcusp]
            else:
                lj0=lj;  cuspev=1.0
            if lr==0.0 and cuspidal==1:
                continue
            for j in range(Qv[jcusp][2]):
                for icusp in range(cuspB-cuspA):
                    iicusp=icusp+cuspA
                    ## If we only have one Y we need the second cusp to estimate the error
                    if cusp_evs[iicusp]<>0.0 and iicusp>0 and (numy==2 or iicusp>1):
                    #if cusp_evs[iicusp]<>0 and iicusp>0:
                        continue
                    #if verbose>0 and jcusp==2 and icusp==2 and l<=2:
                    #    print "Ypb[",iicusp,"][",jcusp,"][",j,"]=",Ypb[iicusp][jcusp][j]
                    if Ypb[iicusp][jcusp][j]==0:
                        continue
                    ckbes=kbesvec[icusp][l][j]*ef1[icusp][jcusp][l][j]
                    for n in range(NB-NA):
                        #if NA+n+Mv[iicusp][0]==0 and cuspidal==1:
                        #    continue
                        ctmpV=ckbes*ef2[icusp][n][j]
                        if lj0>N1:
                            raise ArithmeticError,"Index outside!"
                        V[icusp][n][lj0]+=ctmpV*cuspev
                        if Mv[icusp][0]<0:
                            n1 = n+NB-NA
                            ctmpV1=ckbes*ef2[icusp][n1][j]
                            V[icusp+nc][n][lj0]+=ctmpV1*cuspev
                        if verbose>0 and jcusp==2 and icusp==2 and lj0+Mv[jcusp][0]<=2 and j>Qv[jcusp][2]-1:
                            print "kbes[",icusp,l,j,"]=",kbesvec[icusp][l][j]
                            print "V[",icusp,"][",n,"][",lj0,"]=",V[icusp][n][lj0]

    if verbose>0:
        print "V0[",0,0,"]=",V[0][0][0]
        for jcusp from 0 <= jcusp < nc:
            print "Qfak[",jcusp,"]=",Qfak[jcusp]

    for icusp in range(cuspB-cuspA):
        for n in range(NB-NA):
            #for 0<=l<N1:
            #V[icusp][n][l]=V[icusp][n][l]/Qfak[icusp+cuspA][0]
            for jcusp in range(nc):
                if jcusp>0 and cusp_evs[jcusp]<>0:
                    continue
                for l in range(Mv[jcusp][2]):
                    lj=cusp_offsets[jcusp]+l
                    if lj>N1: # some extra insurance...
                        print "ni=icusp-cuspA=",n
                        print "lj=",cusp_offsets[jcusp],"+",l,"=",lj
                        raise ArithmeticError,"Index outside!"
                    V[icusp][n][lj]=V[icusp][n][lj]/Qfak[jcusp]
                    if Mv[icusp][0]<0:
                        V[icusp+nc][n][lj]=V[icusp+nc][n][lj]/Qfak[jcusp]
    if verbose>0:
        print "V0[",0,0,"]=",V[0][0][0]
    if kbesvec<>NULL:
        for icusp in range(nc):
            if kbesvec[icusp]<>NULL:
                for l in range(Ml):
                    if kbesvec[icusp][l]<>NULL:
                        sage_free(kbesvec[icusp][l])
                sage_free(kbesvec[icusp])
        sage_free(kbesvec)
    if verbose>0:
        print "deal kbbes1"
    if ef1<>NULL:
        for jcusp in range(nc):
            if ef1[jcusp]<>NULL:
                for icusp in range(nc):
                    if ef1[jcusp][icusp]<>NULL:
                        for n in range(Mv[icusp][2]):
                            if ef1[jcusp][icusp][n]<>NULL:
                                sage_free(ef1[jcusp][icusp][n])
                        sage_free(ef1[jcusp][icusp])
                sage_free(ef1[jcusp])
        sage_free(ef1)
    #if verbose>0:
    #    print "NB-NA=",NB-NA
    if ef2<>NULL:
        for icusp in range(cuspB-cuspA):
            iicusp=icusp+cuspA
            if cusp_evs[iicusp]<>0.0 and iicusp>0 and (numy==2 or iicusp>1):
                continue

            if ef2[icusp]<>NULL:
                #if verbose>0:
                #    print "icusp=",icusp
                for n in range(NB-NA):
                    #if verbose>0:
                    #    print "n=",n
                    if ef2[icusp][n]<>NULL:
                        sage_free(ef2[icusp][n])
                    if Mv[icusp][0]<0:
                        n1 = n + NB - NA
                        if ef2[icusp][n1]<>NULL:
                            sage_free(ef2[icusp][n1])
                sage_free(ef2[icusp])
        sage_free(ef2)

    if nvec<>NULL:
        for icusp in range(nc):
            if nvec[icusp]<>NULL:
                sage_free(nvec[icusp])
        sage_free(nvec)
    if cusp_offsets<>NULL:
        sage_free(cusp_offsets)




cdef phase2_coefficient_sum_real_dp_sym(double **Cnew, double ***V,
                                        double ***Cold,int **Mv,int nc,
                                        double *cusp_evs,int* cusp_offsets,int n,
                                        int acusp=0,int bcusp=2,int fstart=0,int fstop=1,
                                        int numy=1,int verbose=0): #symmetric_cusps,
    cdef int icusp,jcusp,l,li,fi
    cdef double cc,vv
    if verbose>2:
        print "in phase2 sum of coefficients ncnt={0}, acusp={1}, bcusp={2}".format(n,acusp,bcusp)
        for icusp in range(nc):
            print "cusp_evs[",icusp,"]=",cusp_evs[icusp]
            #for icusp from 0<=icusp<nc:
            #    print "Mv[",icusp,"]=",Mv[icusp][0],Mv[icusp][1],Mv[icusp][2]
    for fi in range(fstop-fstart):
        for jcusp in range(nc):
            Cnew[fi][jcusp]=0.0
            if cusp_evs[jcusp]<>0.0 and jcusp>0 and (numy==2 or jcusp>1):
                continue
            for icusp in range(nc):
                if cusp_evs[icusp]<>0.0 and icusp>0:
                    continue
                for l in range(Mv[icusp][2]):
                    li=cusp_offsets[icusp]+l
                    vv = V[jcusp][n][li]
                    cc = Cold[fi][icusp][l]
                    Cnew[fi][jcusp] = Cnew[fi][jcusp] + vv*cc
            if Mv[jcusp][0]<0:
                Cnew[fi][jcusp+nc]=0.0
                for icusp in range(nc):
                    if cusp_evs[icusp]<>0.0 and icusp>0:
                        continue
                    for l in range(Mv[icusp][2]):
                        li=cusp_offsets[icusp]+l
                        vv = V[jcusp+nc][n][li]
                        cc = Cold[fi][icusp][l]
                        Cnew[fi][jcusp+nc] = Cnew[fi][jcusp+nc] + vv*cc





# cpdef get_good_Q_for_Y(double Y0,double R,eps):
#     cdef int i
#     cdef double arg,kbes,Y
#     Y=Y0

#     for n in range(10000):
#         arg=6.28*Y*n  ## ~ 2pinY
#         besselk_dp_c(&kbes,R,arg,besprec,1)
#         if abs(kbes*sqrt(Y))<eps:
#             Y=Y*0.99
#         else:
#             break
#     return Y
# @cython.boundscheck(False)
# @cython.cdivision(True)
# cpdef phase_2_cplx_dp_sym_old(S,double R,int NA,int NB,int M0=0,int ndig=10,int dim=1,int cuspstart=0,int cuspstop=1,int fnr=-1,dict Cin={},double Yin=0,int verbose=0,int retf=0,int n_step=50,int do_test=1,method='2c'):
#     #,dict c_evs={}):
#     r"""
#     Computes more coefficients from the given initial set.

#       INPUT:

#         - ``S`` -- Space of Maas waveforms
#         - ``R`` -- real
#         - ``Cin`` -- complex vector
#         - ``NA`` -- integer
#         - ``NB`` -- integer
#         - ``M0`` -- integer (default None) if set we only use this many coeffficients.
#         - ``ndig`` -- integer (default 10) : want new coefficients with this many digits precision.
#         - ``ndig`` -- number of digits of precision
#         - ``dim`` -- assumed dimension
#         - ``cuspstart`` -- Compute coefficients between cusp cuspstart and cuspstop
#         - ``cuspstop`` --
#         - ``fnr`` -- use function nr. fnr (in case of dimension>1)
#         - ``method`` -- '2c' or '2Y' for using two cusps (if possible) or two Y's for error test.

#     OUTPUT:

#     - ``Cout`` -- dictionary of Fourier coefficients

#     EXAMPLE::


#         sage: R=RR(9.53369526135355755434423523592877032382125639510725198237579046413534)
#         sage: IR=CC(0,R)
#         sage: M=MaassWaveForms(Gamma0(1))
#         sage: C=Maassform_coeffs(M,R)
#         sage: D=phase_2(SL2Z,R,C,2,10)



#     """
#     cdef double eps,pi,twopi,Y
#     cdef int nc,sym_type,M00,Mf,Ms,Ml,Ql
#     cdef int **Mv=NULL,**Qv=NULL,*symmetric_cusps=NULL
#     cdef double *Xm=NULL,***Xpb=NULL,***Ypb=NULL
#     cdef double complex ***Cvec=NULL
#     cdef double *Xm2=NULL,***Xpb2=NULL,***Ypb2=NULL
#     cdef double complex ***Cvec2=NULL
#     cdef double complex *cusp_evs=NULL, ***Cold=NULL
#     cdef double * Qfak=NULL
#     cdef int p,n,i,N1=0,NAa=NA
#     cdef int ylp,Qp,Q #,verbose
#     cdef dict NN,XS
#     cdef double complex **Ctmp
#     #cdef double nr0,nr1kbes,kbes0,kbes1
#     #verbose=S._verbose
#     pi=M_PI
#     twopi=pi*2.0
#     eps=10.0**-ndig
#     cdef int fstart,fstop
#     if fnr<0 or fnr>dim:
#         fstart=0; fstop=dim
#     else:
#         fstart=fnr; fstop=fnr+1
#     if nc==1:
#         method='TwoY'
#     cdef double Y2

#     G=S.group(); nc=int(G.ncusps()); sym_type=S._sym_type
#     cdef dict c_evs=S.cusp_evs_dict()
#     cusp_evs=<double complex*>sage_malloc(sizeof(double complex)*nc)
#     for i in range(nc):
#         cusp_evs[i]=<double complex> CC(c_evs.get(i,0))
#     Ctmp = <double complex**>sage_malloc(sizeof(double complex*)*fstop-fstart)
#     for i in range(fstop-fstart):
#         Ctmp[i] = <double complex*>sage_malloc(sizeof(double complex)*nc)
#         for j in range(nc):
#             Ctmp[i][j]=<double complex>0.0
#     Qfak=<double*>sage_malloc(sizeof(double)*nc)
#     if Qfak==NULL: raise MemoryError
#     Mv=<int**>sage_malloc(sizeof(int*)*nc)
#     if Mv==NULL: raise MemoryError
#     Qv=<int**>sage_malloc(sizeof(int*)*nc)
#     if Qv==NULL: raise MemoryError
#     symmetric_cusps=<int*>sage_malloc(sizeof(int)*nc)
#     if symmetric_cusps==NULL: raise MemoryError
#     for i in range(nc):
#         Mv[i]=<int*>sage_malloc(sizeof(int)*3)
#         Qv[i]=<int*>sage_malloc(sizeof(int)*3)
#     cdef int* cusp_offsets
#     cusp_offsets=<int*>sage_malloc(sizeof(int)*nc)
#     ### The starting coefficients
#     #dim = 1 #F._dim
#     if dim<=0:
#         if hasattr(S,"_dim"):
#             dim=S._dim
#         else:
#             dim=1
#     NN=S.set_norm(dim)
#     Y = 0.9*S._group.minimal_height()
#     if M0<=0:
#         sig_on()
#         M0 = get_M_for_maass_dp_c(R,Y,eps)
#         sig_off()

#     cdef int old_v=S._verbose
#     if verbose>0:
#         print "Y,M0=",Y,M0
#     Q=M0+10
#     S._verbose=0
#     if Cin=={}:
#         sig_on()
#         XS=get_coeff_fast_cplx_dp_sym(S,R,Y,M0,Q,NN,cusp_ev=c_evs)
#         sig_off()
#     else:
#         XS=Cin
#     ## Make sure that we have a real eigenvalue.
#     S._verbose=old_v
#     set_Mv_Qv_symm(S,Mv,Qv,Qfak,symmetric_cusps,cusp_evs,cusp_offsets,&N1,&Ml,&Ql,M0,Q,verbose)
#     if dim>1:
#         p=S.get_primitive_p()
#         XS = S.Hecke_eigenfunction_from_coeffs(XS,p,fnr=fnr)
#     else:
#         XS = {0:XS}
#     #print "XSK=",XS.keys()
#     if do_test:
#         test = S.test_Hecke_relation(XS[0])
#         if verbose>0:
#             print "Hecke test=",test
#     else:
#         test=0

#     if abs(test)>sqrt(eps):
#         return {}

#     cdef dict Cnew=dict()
#     for j in XS.keys():
#         Cnew[j]={}
#         for i in XS[j].keys():
#             Cnew[j][i]={}
#     cdef int cuspa,cuspb
#     if retf or nc==1:
#         cuspa=0; cuspb=nc; cuspbb=nc
#         #cuspa=0; cuspb=2
#     else:
#         cuspa=cuspstart; cuspb=cuspstop; cuspbb=max(cuspstop,2)

#     ## First add the coefficients we already have
#     for j in range(fstop-fstart):
#         if not XS.has_key(j+fstart): continue
#         for i in range(cuspa,cuspb):
#             if not XS[j+fstart].has_key(i): continue
#             if not isinstance(XS[j][i],dict):
#                 for n in XS[0].keys():
#                     Cnew[j][i][n]=XS[j+fstart][i]*XS[j+fstart][0][n]
#             else:
#                 for n in XS[j+fstart][i].keys():
#                     Cnew[j][i][n]=XS[j+fstart][i][n]
#     #print "CK=",Cnew[2].keys()
#     if verbose>0:
#         print "after!"
#     Cold=<double complex***>sage_malloc(sizeof(double complex**)*(fstop-fstart))
#     for j in range(fstop-fstart):
#         Cold[j]=<double complex**>sage_malloc(sizeof(double complex*)*nc)
#         for i in range(nc):
#             if cusp_evs[i]==0 or i==0:
#                 Cold[j][i]=<double complex*>sage_malloc(sizeof(double complex)*Mv[i][2])
#                 for n in range(Mv[i][2]):
#                     Cold[j][i][n]=<double complex>CC(XS[j][i][n+Mv[i][0]])
#             else:
#                 Cold[j][i]=<double complex*>sage_malloc(sizeof(double complex)*Mv[i][2])
#                 for n in range(Mv[i][2]):
#                     Cold[j][i][n]=cusp_evs[i]*XS[j][0][n+Mv[i][0]]
#     # using these we want to compute more
#     cdef int nn0,nn1
#     if verbose>0:
#         print "got initial set of coefficients!"
#         print "fstart=",fstart
#         print "fstop=",fstop
#         if verbose>1:
#             for j in range(fstop-fstart):
#                 for n from Mv[0][0]<=n<Mv[0][1]:
#                     print "C[{0}][0][{1}]={2}".format(j,n,XS[j][0][n])
#                 for i in range(nc):
#                     if cusp_evs[i]==0:
#                         print "C[{0}][{1}][1]={2}".format(j,i,XS[j][i][1])
#                     else:
#                         print "Cusp_ev[",i,"]=",cusp_evs[i]
#                 #for n from 0<=n<max(10,Mv[0][1]):
#                 #    print "Cold[0][",n+Mv[0][0],"]=",Cold[j][0][n]
#                 for i in range(nc):
#                     if cusp_evs[i]==0:
#                         nn0 = max(-10,Mv[i][0])
#                         nn1 = min(10,Mv[i][1])
#                         for n in range(nn0,nn1):
#                             if n>Mv[i][2] or n<0:
#                                 print "Cold[{0}][{1}][{2}]=?".format(j,i,n+Mv[i][0])
#                                 continue
#                             print "Cold[{0}][{1}][{2}]={3}".format(j,i,n+Mv[i][0],Cold[j][i][n])
#                             #print "Cold[",i,"][1]=",Cold[j][i][1]

#                     # starting value of Y
#     cdef double Y0
#     if Yin <=0:
#         Y0=twopi/<double>NA
#         if Y0>0.5:
#             Y0=0.5
#     else:
#         Y0=Yin
#     Y0 = min(Y0,S.group().minimal_height())
#     ylp=0; Qp=0
#     Ml=2*M0+1
#     sig_on()
#     Q=max(get_M_for_maass_dp_c(R,Y0,eps)+5,Ml+10)+Qp
#     sig_off()
#     Ql=2*Q
#     Xm=<double*>sage_malloc(sizeof(double)*Ql)
#     if Xm==NULL: raise MemoryError
#     Xpb = <double***> sage_malloc( sizeof(double** ) * nc )
#     if Xpb==NULL: raise MemoryError
#     Ypb = <double***> sage_malloc( sizeof(double** ) * nc )
#     if Ypb==NULL: raise MemoryError
#     for i in range(nc):
#         Xpb[i]=NULL; Ypb[i]=NULL
#         Xpb[i] = <double**>sage_malloc(sizeof(double*) * nc )
#         Ypb[i] = <double**>sage_malloc(sizeof(double*) * nc )
#         if Ypb[i]==NULL or Xpb[i]==NULL:
#             raise MemoryError
#         for j in range(nc):
#             Xpb[i][j]=NULL; Ypb[i][j]=NULL
#             Xpb[i][j] = <double*>sage_malloc(sizeof(double) * Ql )
#             Ypb[i][j] = <double*>sage_malloc(sizeof(double) * Ql )
#             if Ypb[i][j]==NULL or Xpb[i][j]==NULL:
#                 raise MemoryError
#             for n in range(Ql):
#                 Xpb[i][j][n]=<double>0
#                 Ypb[i][j][n]=<double>0
#     Cvec = <double complex***>sage_malloc(sizeof(double complex**) * nc )
#     if Cvec==NULL: raise MemoryError
#     for i in range(nc):
#         Cvec[i] = <double complex**>sage_malloc(sizeof(double complex*) * nc )
#         if Cvec[i]==NULL: raise MemoryError
#         for j in range(nc):
#             Cvec[i][j] = <double complex*>sage_malloc(sizeof(double complex) * Ql )
#             if Cvec[i][j]==NULL: raise MemoryError
#             for n in range(Ql):
#                 Cvec[i][j][n]=<double complex>0
#     set_Mv_Qv_symm(S,Mv,Qv,Qfak,symmetric_cusps,cusp_evs,cusp_offsets,&N1,&Ml,&Ql,M0,Q,verbose)
#     sig_on()
#     pullback_pts_cplx_dp(S,1-Q,Q,Y0,Xm,Xpb,Ypb,Cvec)
#     sig_off()
#     if method=='TwoY':
#         Xm2=<double*>sage_malloc(sizeof(double)*Ql)
#         if Xm2==NULL: raise MemoryError
#         Xpb2 = <double***> sage_malloc( sizeof(double** ) * nc )
#         if Xpb2==NULL: raise MemoryError
#         Ypb2 = <double***> sage_malloc( sizeof(double** ) * nc )
#         if Ypb2==NULL: raise MemoryError
#         for i in range(nc):
#             Xpb2[i]=NULL; Ypb[i]=NULL
#             Xpb2[i] = <double**>sage_malloc(sizeof(double*) * nc )
#             Ypb2[i] = <double**>sage_malloc(sizeof(double*) * nc )
#             if Ypb2[i]==NULL or Xpb2[i]==NULL:
#                 raise MemoryError
#             for j in range(nc):
#                 Xpb2[i][j]=NULL; Ypb[i][j]=NULL
#                 Xpb2[i][j] = <double*>sage_malloc(sizeof(double) * Ql )
#                 Ypb2[i][j] = <double*>sage_malloc(sizeof(double) * Ql )
#                 if Ypb2[i][j]==NULL or Xpb2[i][j]==NULL:
#                     raise MemoryError
#                 for n in range(Ql):
#                     Xpb2[i][j][n]=<double>0
#                     Ypb2[i][j][n]=<double>0
#         Cvec2 = <double complex***>sage_malloc(sizeof(double complex**) * nc )
#         if Cvec2==NULL: raise MemoryError
#         for i in range(nc):
#             Cvec2[i] = <double complex**>sage_malloc(sizeof(double complex*) * nc )
#             if Cvec2[i]==NULL: raise MemoryError
#             for j in range(nc):
#                 Cvec2[i][j] = <double complex*>sage_malloc(sizeof(double complex) * Ql )
#                 if Cvec2[i][j]==NULL: raise MemoryError
#             for n in range(Ql):
#                 Cvec2[i][j][n]=<double complex>0
#         Y2 = 0.98*Y0
#         sig_on()
#         pullback_pts_cplx_dp(S,1-Q,Q,Y2,Xm2,Xpb2,Ypb2,Cvec2)
#         sig_off()

#     if verbose>0:
#         print "in phase 2 with eps=",eps
#         print "range: NA,NB=",NA,NB
#         print "N1=",N1
#         print "Ml,Ql=",Ml,Ql

#     if verbose>0:
#         print "After pullback"

#     cdef double *alphas=NULL
#     alphas=<double*>sage_malloc(sizeof(double)*nc)
#     for i in range(nc):
#         alphas[i]=<double>S.alpha(i)[0]
#     # We compute coefficients at the first two cusps only
#     # and then compare
#     cdef double complex *** V=NULL
#     #cdef int n_step=50
#     cdef int n_a,n_b
#     n_a = NA - Mv[0][0]
#     n_b = n_a + n_step
#     if verbose>0:
#         print "Compute coefficients at cusps from {0} to {1}".format(cuspstart,cuspstop)
#     V=<double complex***>sage_malloc(sizeof(double complex**)*(cuspbb-cuspa))
#     if V==NULL: raise MemoryError
#     for i in range(cuspbb-cuspa):
#         V[i]=<double complex**>sage_malloc(sizeof(double complex*)*(n_step))
#         if V[i]==NULL: raise MemoryError
#         for l in range(n_step):
#             V[i][l]=<double complex*>sage_malloc(sizeof(double complex)*(N1))
#             if V[i][l]==NULL: raise MemoryError
#             for n in range(N1):
#                 V[i][l][n]=<double complex>0.0

#     if verbose>0:
#         print "computing first V"
#     sig_on()
#     compute_V_cplx_dp_sym_for_phase2(V,N1,Xm,Xpb,Ypb,Cvec,
#                                      cusp_evs,alphas,Mv,Qv,Qfak,
#                                      symmetric_cusps,
#                                      R,Y0,nc,n_a,n_b,cuspa,cuspbb,1,
#                                      verbose-1,0)
#     sig_off()
#     if verbose>0:
#         print "after computing first V"
#     cdef double complex ***V2
#     if method=='TwoY':
#         V2=<double complex***>sage_malloc(sizeof(double complex**)*(cuspbb-cuspa))
#         if V2==NULL: raise MemoryError
#         for i in range(cuspbb-cuspa):
#             V2[i]=<double complex**>sage_malloc(sizeof(double complex*)*(n_step))
#             if V2[i]==NULL: raise MemoryError
#             for l in range(n_step):
#                 V2[i][l]=<double complex*>sage_malloc(sizeof(double complex)*(N1))
#                 if V2[i][l]==NULL: raise MemoryError
#                 for n in range(N1):
#                     V2[i][l][n]=<double complex>0.0
#         if verbose>0:
#             print "before computing first V2"
#         sig_on()
#         compute_V_cplx_dp_sym_for_phase2(V2,N1,Xm2,Xpb2,Ypb2,Cvec2,
#                                      cusp_evs,alphas,Mv,Qv,Qfak,
#                                      symmetric_cusps,
#                                      R,Y2,nc,n_a,n_b,cuspa,cuspbb,1,
#                                      verbose-1,0)
#         sig_off()
#         if verbose>0:
#             print "after computing first V2"


#     cdef double besprec,Y2pi,diff,diffmax,Y2pi2,Y02
#     besprec=1.0E-14
#     cdef int ncnt=0,redov=0,fi
#     cdef double *nr,*kbes
#     nr = <double*>sage_malloc(sizeof(double)*nc)
#     kbes = <double*>sage_malloc(sizeof(double)*nc)
#     for yn in range(1000):
#         try:
#             #Y01=Y0;
#             Y02=0.995*Y0
#             try:
#                 sig_on()
#                 Q=max(get_M_for_maass_dp_c(R,Y0,eps)+5,M0+10)+Qp
#                 sig_off()
#             except ArithmeticError:
#                 return {}
#             sqrtY=sqrt(Y0)
#             Y2pi=Y0*twopi
#             if verbose>0:
#                 print "Y[0]=",Y0
#                 print "Q(Y)=",Q
#                 print "Y2pi=",Y2pi
#                 print "NA,NB=",NAa,NB
#                 print "cuspa,cuspb,cuspbb=",cuspa,cuspb,cuspbb
#             #kbdict=_setup_kbessel_dict(Ms,Mf,Qs,Qf,nc,IR,pb,lvec,mp_ctx)
#             for n in range(NAa,NB+1):
#                 ## First check if the current Y is good with respect to the Bessel-factor
#                 if verbose>0:
#                     print "n=",n," ncnt=",ncnt
#                 decrease_y=0; redov=0
#                 for jcusp in range(cuspa,cuspbb):
#                     nr[jcusp]=fabs(<double>n+alphas[jcusp])*Y2pi
#                     besselk_dp_c(&kbes[jcusp],R,nr[jcusp],besprec,1)
#                     kbes[jcusp]=sqrtY*kbes[jcusp]
#                     # Compute the V matrix elements
#                 if fabs(kbes[0])<eps or fabs(kbes[1])<eps:
#                     if verbose>0:
#                         for jcusp in range(cuspa,cuspbb):
#                             print "K_iR(2piY(n+a(0)))=",kbes[jcusp]
#                     decrease_y=1
#                 elif ncnt<n_step:
#                     phase2_coefficient_sum_cplx_dp_sym(Ctmp,V,Cold,Mv,nc,cusp_evs,
#                                                        cusp_offsets,ncnt,cuspa,cuspbb,fstart,fstop,verbose)
#                     for fi in range(fstop-fstart):
#                         for jcusp in range(cuspa,cuspbb):
#                             if verbose>1:
#                                 print "Ctmp{0}={1}".format(jcusp,Ctmp[fi][jcusp])
#                             Ctmp[fi][jcusp]=Ctmp[fi][jcusp]/kbes[jcusp]
#                             if verbose>1:
#                                 print "sqrt({0})K_i{1}({2})".format(sqrtY,R,nr[jcusp])
#                                 print "K_iR(2piY(n+a(0)))=",kbes[jcusp]
#                             if verbose>0:
#                                 print "C{0}tmp{1}/K_iR[{2}]={3}".format(fi,jcusp,R,Ctmp[fi][jcusp])
#                     # We check the error for all functions
#                     diff=1; diffmax=0
#                     for fi in range(fstop-fstart):
#                         if cusp_evs[1]<>0:
#                             diff=cabs(Ctmp[fi][1]-cusp_evs[1]*Ctmp[fi][0])
#                         if diff>eps:
#                             diff=fabs(cabs(Ctmp[fi][1])-cabs(Ctmp[fi][0]))
#                         if verbose>0:
#                             print "err[diff][",n,"]=",diff
#                         if fabs(cabs(Ctmp[fi][1])+cabs(Ctmp[fi][0]))< 2*eps and n<M0-10:
#                             ## This is to prevent Ctmp[1] and Ctmp[2]
#                             ## to be simultaneously close to zero by accident
#                             if diff<eps:
#                                 if abs(Ctmp[fi][0]-XS[fi][0][n])>diff*1.0E6:
#                                     diff=max(diff,abs(Ctmp[fi][0]-XS[fi][0][n]))
#                                     if verbose>0:
#                                         print "err[Cold][",fi,n,"]=",diff
#                         if diff>diffmax:
#                             diffmax=diff
#                     if verbose>0:
#                         print "diffmax=",diff
#                     if diffmax<eps:
#                         if verbose>0:
#                             print "Setting coeficient ",n
#                         for fi in range(fstop-fstart):
#                             for jcusp in range(cuspa,cuspbb):
#                                 Cnew[fi][jcusp][n]=Ctmp[fi][jcusp]
#                             # # improve the coefficients we use
#                             if n<=M0:
#                                 Cold[fi][0][n]=Cnew[fi][0][n]
#                                 Cold[fi][1][n]=Cnew[fi][1][n]
#                             #if retf:
#                             #    for jcusp from cuspa<=jcusp<cuspb:
#                             #        XS[jcusp][n]=Cnew[jcusp][n]
#                     else:
#                         decrease_y=1
#                     ncnt+=1
#                 else:
#                     redov=1
#                 if decrease_y==1:
#                     if verbose>0:
#                         print "decreasing Y!"
#                     NAa=n; ylp=ylp+1
#                     if ylp>4:  # We increase Q more than apriori needed in every second step
#                         Qp=Qp+10; ylp=0
#                     else:
#                         Y0 = get_good_Y_for_n(Y0*0.9,R,n+n_step,eps)

#                     try:
#                         sig_on()
#                         Q=max(get_M_for_maass_dp_c(R,Y0,eps)+5,M0+10)+Qp
#                         sig_off()
#                     except ArithmeticError:
#                         return {}
#                     #Q=max(get_M_for_maass_dp_c(R,Y0,eps)+5,Ml+10)+Qp
#                     set_Mv_Qv_symm(S,Mv,Qv,Qfak,symmetric_cusps,cusp_evs,cusp_offsets,&N1,&Ml,&Ql,M0,Q,verbose-2)
#                     # If we change Y we also need to recompute the pullback
#                     if Xm<>NULL:
#                         sage_free(Xm)
#                     Xm=<double*>sage_malloc(sizeof(double)*Ql)
#                     #if verbose>0:
#                     #   print "Allocated Xm!"
#                     for i in range(nc):
#                         if Xpb[i]<>NULL:
#                             for j in range(nc):
#                                 if Xpb[i][j]<>NULL:
#                                     sage_free(Xpb[i][j])
#                                 if Ypb[i][j]<>NULL:
#                                     sage_free(Ypb[i][j])
#                                 if Cvec[i][j]<>NULL:
#                                     sage_free(Cvec[i][j])
#                     for i in range(nc):
#                         for j in range(nc):
#                             Xpb[i][j] = <double*>sage_malloc(sizeof(double) * Ql )
#                             Ypb[i][j] = <double*>sage_malloc(sizeof(double) * Ql )
#                             Cvec[i][j] = <double complex*>sage_malloc(sizeof(double complex) * Ql )
#                             for l in range(Ql):
#                                 Xpb[i][j][l]=<double>0
#                                 Ypb[i][j][l]=<double>0
#                                 Cvec[i][j][l]=<double complex>0

#                     #if verbose>0:
#                     #    print "Allocated Xpb Ypb!"
#                     #sig_on()
#                     # We need to recompute everything
#                     pullback_pts_cplx_dp(S,1-Q,Q,Y0,Xm,Xpb,Ypb,Cvec)
#                     redov=1
#                 if redov==1:
#                     NAa=n
#                     if verbose>1:
#                         print "Recompute V!"
#                     if n>M0 and cuspb>1:  ## We don't need to compute expansions for all cusps anymore
#                         cuspa=cuspstart; cuspb=cuspstop; cuspbb=max(cuspb,2)
#                     for i in range(cuspbb-cuspa):
#                         for l in range(n_step):
#                             for j in range(N1):
#                                 V[i][l][j]=<double complex>0
#                     sig_on()
#                     compute_V_cplx_dp_sym_for_phase2(V,N1,Xm,Xpb,Ypb,Cvec,
#                                                      cusp_evs,alphas,Mv,Qv,Qfak,
#                                                      symmetric_cusps,
#                                                      R,Y0,nc,n-Mv[0][0],n-Mv[0][0]+n_step,cuspa,cuspbb,1,
#                                                      verbose-1,0)
#                     redov=0
#                     ncnt=0
#                     sig_off()
#                     raise StopIteration()
#             raise StopIteration()
#         except StopIteration:
#             #print "Stop iteration: n=",n
#             if n>NB:
#                 break
#             else:
#                 continue
#         except ArithmeticError as arithm:
#             if verbose>0:
#                 s = str(arithm)
#                 s += "\n Cought arithmetic Error for R,n,Y=",R,n,Y
#                 print s
#             return {}
#     ## Deallocating stuff
#     #print "too free!"
#     sage_free(alphas)
#     sage_free(nr)
#     sage_free(kbes)
#     sage_free(Xm)
#     sage_free(symmetric_cusps)
#     sage_free(cusp_evs)
#     if Ctmp<>NULL:
#         for j in range(fstop-fstart):
#             if Ctmp[j]<>NULL:
#                 sage_free(Ctmp[j])
#         sage_free(Ctmp)

#     #print "freed 1!"
#     if Qfak<>NULL:
#         sage_free(Qfak)
#     if Mv<>NULL:
#         for i in range(nc):
#             if Mv[i]<>NULL:
#                 sage_free(Mv[i])
#         sage_free(Mv)
#     if Qv<>NULL:
#         for i in range(nc):
#             if Qv[i]<>NULL:
#                 sage_free(Qv[i])
#         sage_free(Qv)
#     if Cold<>NULL:
#         for j in range(fstop-fstart):
#             if Cold[j]<>NULL:
#                 for i in range(nc):
#                     if Cold[j][i]<>NULL:
#                         sage_free(Cold[j][i])
#                 sage_free(Cold[j])
#         sage_free(Cold)
#     if Xpb<>NULL:
#         for i in range(nc):
#             if Xpb[i]<>NULL:
#                 for j in range(nc):
#                     if Xpb[i][j]<>NULL:
#                         sage_free(Xpb[i][j])
#                 sage_free(Xpb[i])
#         sage_free(Xpb)
#     if Ypb<>NULL:
#         for i in range(nc):
#             if Ypb[i]<>NULL:
#                 for j in range(nc):
#                     if Ypb[i][j]<>NULL:
#                         sage_free(Ypb[i][j])
#                 sage_free(Ypb[i])
#         sage_free(Ypb)
#     if Cvec<>NULL:
#         for i in range(nc):
#             if Cvec[i]<>NULL:
#                 for j in range(nc):
#                     if Cvec[i][j]<>NULL:
#                         sage_free(Cvec[i][j])
#                 sage_free(Cvec[i])
#         sage_free(Cvec)
#     if V<>NULL:
#         for i in range(cuspbb-cuspa):
#             if V[i]<>NULL:
#                 for j in range(n_step):
#                     if V[i][j]<>NULL:
#                         sage_free(V[i][j])
#                 sage_free(V[i])
#         sage_free(V)

#     return Cnew
