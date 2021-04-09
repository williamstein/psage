# cython: profile=True
r"""
Pullback algorithms optimized  for various settings.



"""
from cysignals.memory cimport sig_free,sig_malloc
from cysignals.signals cimport sig_on,sig_off

from psage.rings.mp_cimports cimport *
import logging

log = logging.getLogger(__name__)

from sage.rings.complex_mpc import MPComplexField
from sage.rings.complex_field import ComplexField
from sage.rings.real_mpfr import RealField

from sage.modules.all import vector
from sage.matrix.all import matrix
from sage.rings.integer import Integer
from sage.rings.real_mpfr import RR
import sage.structure.element
from sage.all import copy
from psage.modform.arithgroup.mysubgroup import is_Hecke_triangle_group,MySubgroup

import cython
cdef extern from "math.h":
    double fabs(double)
    double fmax(double,double)
    int ceil(double) 
    double sqrt(double)
    double pow(double,double) 


cdef extern from "complex.h":
    ### "Hack" suggested by Robert Bradshaw in 2009.
    ### TODO: Is there now a better way for pure c double complex?
    ctypedef double cdouble "double complex"
    cdef double creal(double complex)
    cdef double cimag(double complex)
    cdef double complex _Complex_I
    cdef double carg(double complex)
    cdef double cabs(double complex)
    cdef double complex cpow(double complex,double complex)
    cdef double complex cexp(double complex)
#    cdef complex CMPLXF(float,float)
#    cdef double complex CMPLX(double,double)
cdef complex CMPLXF(float x,float y):
    return x+_Complex_I*y

cdef double complex CMPLX(double x,double y):
    return x+_Complex_I*y

    
from sage.modular.arithgroup.congroup_sl2z import SL2Z

from psage.modform.arithgroup.mysubgroups_alg import apply_sl2z_map,pullback_general_group_dp,pullback_general_group
from psage.modform.arithgroup.mysubgroups_alg import normalize_point_to_cusp_mpfr,pullback_to_Gamma0N_mpfr,apply_sl2z_map_mpfr,normalize_point_to_cusp_dp,apply_sl2z_map_dp,normalize_point_to_cusp_mpmath
from psage.modform.arithgroup.mysubgroups_alg cimport _apply_sl2z_map_dp,_apply_sl2z_map_mpfr,_pullback_to_Gamma0N_dp,pullback_to_hecke_triangle_mat_c_mpfr,pullback_to_Gamma0N_mpfr_c,normalize_point_to_cusp_mpfr_c,_normalize_point_to_cusp_dp,_normalize_point_to_cusp_real_dp,SL2Z_elt,_normalize_point_to_cusp_mpfr,closest_vertex_dp_c

from sage.all import CC,save
from psage.modules.vector_real_mpfr_dense cimport Vector_real_mpfr_dense 


import mpmath


cpdef pullback_pts_dp(S,int Qs,int Qf,double Y,double weight=0,holo=False):
    cdef double* Xm_t=NULL
    cdef double*** Xpb_t=NULL
    cdef double*** Ypb_t=NULL
    cdef double complex*** Cvec_t=NULL
    cdef int Ql,i,j,n,nc
    nc= S.group().ncusps()
    Ql=Qf-Qs+1
    Xm_t=<double*>sig_malloc(sizeof(double)*Ql)
    if Xm_t==NULL: raise MemoryError
    Xpb_t = <double***> sig_malloc( sizeof(double** ) * nc )
    if Xpb_t==NULL: raise MemoryError
    Ypb_t = <double***> sig_malloc( sizeof(double** ) * nc )
    if Ypb_t==NULL: raise MemoryError
    for i in range(nc):
        Xpb_t[i] = <double**>sig_malloc(sizeof(double*) * nc )
        Ypb_t[i] = <double**>sig_malloc(sizeof(double*) * nc )
        if Ypb_t[i]==NULL or Xpb_t[i]==NULL:
            raise MemoryError
        for j in range(nc):
            Xpb_t[i][j] = <double*>sig_malloc(sizeof(double) * Ql )
            Ypb_t[i][j] = <double*>sig_malloc(sizeof(double) * Ql )
            if Ypb_t[i][j]==NULL or Xpb_t[i][j]==NULL:
                raise MemoryError
            for n in range(Ql):
                Xpb_t[i][j][n]=<double>0
                Ypb_t[i][j][n]=<double>0
    Cvec_t = <double complex***>sig_malloc(sizeof(double complex**) * nc )
    if Cvec_t==NULL: raise MemoryError
    for i in range(nc):
        Cvec_t[i] = <double complex**>sig_malloc(sizeof(double complex*) * nc )
        if Cvec_t[i]==NULL:
            raise MemoryError
        for j in range(nc):
            Cvec_t[i][j] = <double complex*>sig_malloc(sizeof(double complex) * Ql )
            if Cvec_t[i][j]==NULL:
                raise MemoryError
            for n in range(Ql):
                Cvec_t[i][j][n]=0 
    cdef int res
    res = pullback_pts_cplx_dp(S,Qs,Qf,Y,Xm_t,Xpb_t,Ypb_t,Cvec_t)
    if res==1:
        raise ArithmeticError,"Y is too large!"
    cdef dict Xm,Xpb,Ypb,Cvec,pb
    Xm=dict(); Xpb=dict() 
    Ypb=dict();Cvec=dict()
    for n in range(Ql):
        Xm[n]=Xm_t[n]
    for i in range(nc):
        Cvec[i]=dict();Xpb[i]=dict();Ypb[i]=dict()
        for j  in range(nc):
            Cvec[i][j]=dict();Xpb[i][j]=dict();Ypb[i][j]=dict()
            for n in range(Ql):
                Cvec[i][j][n]=Cvec_t[i][j][n]
                Xpb[i][j][n]=Xpb_t[i][j][n]
                Ypb[i][j][n]=Ypb_t[i][j][n]
    pb=dict()
    pb['xm']=Xm; pb['xpb']=Xpb; pb['ypb']=Ypb; pb['cvec']=Cvec
    if Ypb_t<>NULL:
        for i in range(nc):
            if Ypb_t[i]<>NULL:
                for j  in range(nc):
                    if Ypb_t[i][j]<>NULL:
                        sig_free(Ypb_t[i][j])
                sig_free(Ypb_t[i])
        sig_free(Ypb_t)
    if Xpb_t<>NULL:
        for i in range(nc):
            if Xpb_t[i]<>NULL:
                for j  in range(nc):
                    if Xpb_t[i][j]<>NULL:
                        sig_free(Xpb_t[i][j])
                sig_free(Xpb_t[i])
        sig_free(Xpb_t)
    if Cvec_t<>NULL:
        for i in range(nc):
            if Cvec_t[i]<>NULL:
                for j  in range(nc):
                    if Cvec_t[i][j]<>NULL:
                        sig_free(Cvec_t[i][j])
                sig_free(Cvec_t[i])
        sig_free(Cvec_t)
    return pb
@cython.cdivision(True)
#cdef int pullback_pts_cplx_dp(S,int Qs,int Qf,double Y,double *Xm,double *** Xpb,double*** Ypb,double complex ***Cvec):
cdef int pullback_pts_cplx_dp_sym(S,int **Qv,double Y,double *Xm,double *** Xpb,double*** Ypb,double complex ***Cvec):
    r""" Computes a whole array of pullbacked points using double precision

    INPUT:

    - ``S``  -- AutomorphicFormSpace
    - ``Qs`` -- integer
    - ``Qf`` -- integer
    - ``Y``  -- real >0
    - ``weight`` -- (optional) real
    - ``holo``   -- (False) logical

    
    OUTPUT:

    
    - ``pb`` -- dictionary with entries:
       - 'xm'   -- real[Qf-Qs+1]  : x_m=2*pi*(1-2*m)/2(Qf-Qs+1)
       - 'xpb'  -- real[0:nc,0:nc,Qf-Qs+1]  : real part of pullback of x_m+iY
       - 'ypb'  -- real[0:nc,0:nc,Qf-Qs+1]  : imag. part of pullback of x_m+iY
       - 'cvec' -- complex[0:nc,0:nc,Qf-Qs+1] : mult. factor

    RETURNS::


    - 0 -- if everything is ok
    - 1 -- if Y is too large

    EXAMPLES::


        sage: G=MySubgroup(Gamma0(2))
        sage: [Xm,Xpb,Ypb,Cv]=pullback_pts(G,1,3,0.1)

    Note that we need to have $Y<Y_{0}$ so that all pullbacked points
    are higher up than Y.

        sage: pullback_pts(G,1,10,0.44)
        Traceback (most recent call last):
        ...
        ArithmeticError: Need smaller value of Y

    
    """
#    sig_on()
    G = S.group()
    multiplier = S.multiplier()
    #character = S.character
    holo = S.is_holomorphic()
    cdef double weight
    weight = S.weight()
    if weight<>0 or not multiplier.is_trivial():
        non_trivial=True
    else:
        non_trivial=False
    trivial_mult = multiplier.is_trivial()
    cdef double pi=3.1415926535897932384626433
    cdef double twopi,twoQl
    cdef double complex twopii,mp_i,m1,m2,m3,tmp,v
    cdef double ep0=1E-10
    cdef double tmpabs,tmparg,tmparg1,tmparg2,tmparg3
    twopi=2.0*pi
    twopii = twopi* _Complex_I
    #=complex(0,twopi)
    if(holo and weight==0):
        raise Exception,"Must support a weight for holomorphic forms!"
    cdef double fourQ
    cdef int wi,wj
    cdef int A[4]
    cdef int pba,pbb,pbc,pbd
    cdef list U
    #    AA=[1,0,0,1]
#    Tj=[1,0,0,1]
    cdef double swi,swj,ar,br,cr,dr
    cdef int ci,cj
    cdef int verbose=S._verbose-1
    cdef double x,y,x1,y1,x2,y2,x3,y3,xx,yy
    cdef int** normalizers=NULL
    cdef int** cusp_maps=NULL
    cdef int** vertex_maps=NULL
    cdef int* vertex_cusp=NULL
    cdef double *vertex_widths=NULL
    cdef int a,b,c,d,i,j,Ql,vi,vj
    cdef double wi_d,wj_d
    cdef ComplexNumber ctmp
    cdef int nc=G._ncusps
    cdef int nv=G._nvertices
    cdef int use_int=1
    cdef int is_Gamma0,nreps
    is_Gamma0=<int>G._is_Gamma0
    cdef int*** reps=NULL
    cdef int N
    cdef double xmm
    if verbose>0:
        print "In pullback_pts_cplx_dp_sym!"
    if is_Gamma0==1:
        nreps=G._index
        N=G._level
        reps= <int ***> sig_malloc(sizeof(int**) * nreps)
        for j in range(nreps):
            reps[j]=NULL
            reps[j]=<int **> sig_malloc(sizeof(int*) * 2)
            if reps[j]==NULL:
                raise MemoryError
            reps[j][0]=NULL;reps[j][1]=NULL
            reps[j][0]=<int *> sig_malloc(sizeof(int) * 2)
            if reps[j][0]==NULL:
                raise MemoryError
            reps[j][1]=<int *> sig_malloc(sizeof(int) * 2)
            if reps[j][1]==NULL:
                raise MemoryError
            reps[j][0][0]=G.coset_reps()[j][0]
            reps[j][0][1]=G.coset_reps()[j][1]
            reps[j][1][0]=G.coset_reps()[j][2]
            reps[j][1][1]=G.coset_reps()[j][3]
#            reps[j][0][0]=G._coset_reps_v0[j][0]
#            reps[j][0][1]=G._coset_reps_v0[j][1]
#            reps[j][1][0]=G._coset_reps_v0[j][2]
#            reps[j][1][1]=G._coset_reps_v0[j][3]

            
    if verbose>0:
        print "nc=",nc
    if not isinstance(G._cusp_data[0]['width'],(int,Integer)):
        use_int=0
    normalizers=<int**>sig_malloc(sizeof(int*)*nc)
    if normalizers==NULL: raise MemoryError
    for i in range(nc):
        normalizers[i]=<int*>sig_malloc(sizeof(int)*4)
        a,b,c,d=G._cusp_data[i]['normalizer']
        normalizers[i][0]=a
        normalizers[i][1]=b
        normalizers[i][2]=c
        normalizers[i][3]=d
        if verbose>2:
            print "Normalizer[",i,"]=",a,b,c,d
    vertex_maps=<int**>sig_malloc(sizeof(int*)*nv)
    if vertex_maps==NULL: raise MemoryError
    cusp_maps=<int**>sig_malloc(sizeof(int*)*nv)
    if cusp_maps==NULL: raise MemoryError
    vertex_widths=<double*>sig_malloc(sizeof(double)*nv)
    if vertex_widths==NULL: raise MemoryError
    vertex_cusp=<int*>sig_malloc(sizeof(int)*nv)    
    if vertex_cusp==NULL: raise MemoryError
    for i in range(nv):
        cusp_maps[i]=<int*>sig_malloc(sizeof(int)*4)
        cusp_maps[i][0]=<int>G._cusp_maps[i][0]
        cusp_maps[i][1]=<int>G._cusp_maps[i][1]
        cusp_maps[i][2]=<int>G._cusp_maps[i][2]
        cusp_maps[i][3]=<int>G._cusp_maps[i][3]
        vertex_maps[i]=<int*>sig_malloc(sizeof(int)*4)
        vertex_maps[i][0]=<int>G._vertex_maps[i].a()
        vertex_maps[i][1]=<int>G._vertex_maps[i].b()
        vertex_maps[i][2]=<int>G._vertex_maps[i].c()
        vertex_maps[i][3]=<int>G._vertex_maps[i].d()
        vertex_widths[i]=<double>G._vertex_widths[i]
        vertex_cusp[i]=<int>G._vertex_data[i]['cusp']
    cdef double* widths
    widths=<double*>sig_malloc(sizeof(double)*nc)
    for i in range(nc):
        if verbose>3:
            print "width[",i,"]=",G._cusp_data[i]['width'],type(G._cusp_data[i]['width'])
        widths[i]=<double>G._cusp_data[i]['width']
    cdef int dir_char=0
    cdef double complex *charvec=NULL
    cdef int modulus
    #if verbose>0:
    #    print "Here1.5"
    if not multiplier.is_trivial():
        dir_char=1
        modulus = multiplier._character.modulus()
        charvec=<double complex*>sig_malloc(sizeof(double complex)*modulus)
        for i in range(modulus):
            mc = multiplier._character(i)
            if hasattr(mc,"complex_embedding"):
                charvec[i]=<double complex> multiplier._character(i).complex_embedding()
            else:
                charvec[i]=<double complex>CC(multiplier._character(i))
            if verbose>1:
                print "chi(",i,")=",charvec[i]
    #if verbose>0:
    #    print "Here2"
    Ql = 0
#    cdef int some_nonsymmetric_cusp = 0
    for ci in range(nc):
        if Qv[ci][2]>Ql:
            Ql = Qv[ci][2]
#                Ql=Qf-Qs+1
#    if Qs<0:
#        fourQ=<double>2*(Qf-Qs+1)
#    else:
#        fourQ=<double>(4*Qf)
    cdef int Q = Qv[0][1]
    fourQ = <double> (4*Q)
    for j in range(Qv[0][1]): #from Qs<= j <= Qf:
        Xm[j]=<double>(2*(j+1)-1)/fourQ        
    if verbose>2:
        print "use_int=",use_int
        print "Q=",Q
        print "fourQ=",fourQ
    res={}
                
    for ci in range(nc): #in G._cusps:
        # if ci==0:
        #     verbose=3
        # else:
        #     verbose=0
        #ci = G._cusps[i]
        if verbose>2:
            print "------------------------------------------"
            print "ci=",ci
        #cii=G._cusps.index(ci)
        if use_int==1:
            wi = <int>widths[ci] #G._cusp_data[ci]['width'] #(ci)
            swi=sqrt(<double>wi)
        else:
            wi_d=widths[ci] #<double>G._cusp_data[ci]['width'] #(ci)
            swi=sqrt(wi_d)
        for j in range(Qv[ci][2]): #Ql, 1-Q,Q):
            if j < Q:
                xmm = Xm[j]
            else:
                xmm = -Xm[2*Q-1-j]
            if verbose>1:
                print "Xm[{0}]={1}".format(j,xmm)
            x = xmm
            y=Y
            a=normalizers[ci][0]; b=normalizers[ci][1]
            c=normalizers[ci][2]; d=normalizers[ci][3]
            if verbose>2:
                print "normalizer[",ci,"]=",a,b,c,d
                print "width = ",wi
                print "x0,y0=",x,y
            if use_int==1:
                _normalize_point_to_cusp_dp(&x,&y,a,b,c,d,wi)
            else:
                _normalize_point_to_cusp_real_dp(&x,&y,a,b,c,d,wi_d)
            #[x,y]   = normalize_point_to_cusp_dp(G,ci,Xm[j],Y)
            if verbose>2:
                print "N(z0)=",x,y
            if is_Gamma0 == 1:
                #x1,y1,pba,pbb,pbc,pbd =  G.pullback(x,y,ret_mat=0)
                _pullback_to_Gamma0N_dp(reps,nreps,N,&x,&y,
                                        &pba,&pbb,&pbc,&pbd,0)
                x1=x; y1=y
            else:
#                x1,y1,pba,pbb,pbc,pbd =  G.pullback(x,y,ret_mat=0)
                x1,y1,pba,pbb,pbc,pbd =  pullback_general_group_dp(G,x,y,ret_mat=1,verbose=verbose)
                if verbose>2:
                    print "Pbz=",x1,y1
            ## Want to replace this with a cpdef'd function
            #cj,vj= G.closest_cusp(x1,y1,vertex=1)
            vj = closest_vertex_dp_c(nv,vertex_maps,vertex_widths,&x1,&y1)
            cj = vertex_cusp[vj]
            if pba<=0 and pbb<=0 and pbc<=0 and pbd<=0:
                pba=-pba; pbb=-pbb; pbc=-pbc; pbd=-pbd
            if verbose>2:
                print "x0,y0=",x,y
                print "Q=",Q
                print "j=",j+Qv[ci][0]
                print "Xm[{0}],Y={1},{2}".format(j,xmm,Y)
                print "x,y=",x,y
                print "Pbmap=",pba,pbb,pbc,pbd
                print "x1,y1=",x1,y1
                print "cj,vj=",cj,vj
            #cjj=G._cusps.index(cj)
            if use_int==1:
                wj=<int>widths[cj] #G._cusp_data[cj]['width']
                swj=sqrt(<double>wj) 
            else:
                wj_d=widths[cj] #G._cusp_data[cj]['width']
                swj=sqrt(wj_d)
            a = cusp_maps[vj][0]   #a,b,c,d=G._vertex_data[cj]['cusp_map'] #[v]
            b = cusp_maps[vj][1]
            c = cusp_maps[vj][2]
            d = cusp_maps[vj][3]
            if verbose>2:
                if cj==1:
                    print ">>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<"
                print "cj=",cj
                print "U1(",vj,")=",a,b,c,d
                print "x00,y00=",x1,y1
            x2=x1; y2=y1
            if a<>1 or b<>0 or c<>0 or d<>1:
                #if verbose>2:
                #    print "here1 type=",type(a)
                #[x2,y2]=_apply_sl2z_map_dp(x1,y1,U[0],U[1],U[2],U[3])
                _apply_sl2z_map_dp(&x2,&y2,a,b,c,d)
            if verbose>2:
                print "x01,y01=",x2,y2
            #    print "cj=",cj
            #x2=x1; y2=y1
            x3=x2; y3=y2
            a=normalizers[cj][3]; b=-normalizers[cj][1]
            c=-normalizers[cj][2]; d=normalizers[cj][0]
            if verbose>2:
                print "Normalizer: a,b,c,d=",a,b,c,d
            #_normalize_point_to_cusp_dp(&x3,&y3,a,b,c,d,wj,1)
            if use_int==1:
                _normalize_point_to_cusp_dp(&x3,&y3,a,b,c,d,wj,1)
            else:
                _normalize_point_to_cusp_real_dp(&x3,&y3,a,b,c,d,wj_d,1)
            if verbose>2:
                print "width=",wj
                print "x02,y02=",x3,y3
            #[x3,y3] = normalize_point_to_cusp_dp(G,cj,x2,y2,inv=1)
            Xpb[ci][cj][j]=x3*twopi
            # Recall that Ypb must be greater than Y otherwise we need a better Y
            if y3>=Y+ep0:
                Ypb[ci][cj][j]=y3
            else:
                if verbose>0:
                    print "ci,cj=",ci,cj
                    print "Tj=",pba,pbb,pbc,pbd
                    print "wj=",wj
                    print "wj_d=",wj_d
                    print "x,y=",x,y
                    print "x1,y1=",x1,y1
                    print "x2,y2=",x2,y2
                    print "x3,y3=",x3,y3
                    print "Normalizer=",a,b,c,d
                    print "ep0=",ep0
                return 1
                #raise ArithmeticError,"Need smaller value of Y. Got:{0}".format(Y) 
            if verbose>2:
                print "ci,cj=",ci,cj
                print "Tj=",pba,pbb,pbc,pbd
                print "Xm=",xmm
                print "x,y=",x,"\n",y
                print "x1,y1=",x1,"\n",y1
                print "x2,y2=",x2,"\n",y2
                print "x3,y3=",x3,"\n",y3
                print "Xpb=",Xpb[ci][cj][j]
                print "Ypb=",Ypb[ci][cj][j]
            # We also get the multiplier if we need it
            if non_trivial:
                if weight<>0:
                    if verbose>2:
                        print "here holo=",holo," weight=",weight
                        print "A=",A[0],A[1],A[2],A[3]
                    cr = <double>(normalizers[ci][2])*swi
                    dr = <double>(normalizers[ci][3])/swi
                    if verbose>2:
                        print "cr,dr=",cr,dr
                    j_fak_dp_c(cr,dr,xmm,Y,weight,holo,&tmpabs,&tmparg1)
                    m1=cexp(-_Complex_I*tmparg1)
                    if verbose>2:
                        print "-tmparg1=",-tmparg1
                        print "m1=",m1
                    cr = <double>(normalizers[cj][2])*swj
                    dr = <double>(normalizers[cj][3])/swj
                    j_fak_dp_c(cr,dr,x3,y3,weight,holo,&tmpabs,&tmparg2)
                    m2=cexp(_Complex_I*tmparg2)
                    if verbose>2:
                        print "tmparg2=",tmparg2
                    #l=mat_mul_lit(cusp_maps[vj],Tj)
                    #A[0]=l[0]; A[1]=l[1]; A[2]=l[2]; A[3]=l[3]
                    _mat_mul_list(cusp_maps[vj][0],cusp_maps[vj][1],cusp_maps[vj][2],cusp_maps[vj][3],pba,pbb,pbc,pbd,A)
                    #A=A**-1
                    cr=<double>(-A[2]); dr=<double>(A[0])
                    j_fak_dp_c(cr,dr,x2,y2,weight,holo,&tmpabs,&tmparg3)
                    m3=cexp(_Complex_I*tmparg3)
                    tmp=cexp(_Complex_I*(tmparg2-tmparg1+tmparg3))
                    if verbose>2:
                        print "tmparg3=",tmparg3
                        print "A[2],A[3]=",A[2],A[3]
                        print "A[3]=",A[3],"% modulus:",vi
                        print "m1,m2,m3=",m1,m2,m3
                else:
                    _mat_mul_list(cusp_maps[vj][0],cusp_maps[vj][1],cusp_maps[vj][2],cusp_maps[vj][3],pba,pbb,pbc,pbd,A)
                    tmp=CMPLX(1.0,0.0) #<double complex>1.0
                if dir_char:
                    vi=A[3] % modulus
                    if vi<0:
                        vi=vi+modulus
                    v=charvec[vi]
                    if verbose>2:
                        print "A[3]=",A[3],"% modulus:",vi
                        print "m1,m2,m3=",m1,m2,m3
                        print "tmp=",tmp
                        print "mult=",v 

                    #v = <double complex> ctmp.real()+_Complex_I*ctmp.imag() #CC(multiplier(AA)) #character(AA[1,1])a
                    tmp = tmp/v
                elif not trivial_mult:
                    if verbose>1:
                        print "mult=",multiplier([A[0],A[1],A[2],A[3]])
                    ctmp=multiplier([A[0],A[1],A[2],A[3]]).complex_embedding()
                    v = <double complex> ctmp.real()+_Complex_I*ctmp.imag() #CC(multiplier(AA)) #character(AA[1,1])a
                    tmp = tmp/v
                Cvec[ci][cj][j]=tmp # mp_ctx.mpc(0,tmp)
            else:
                Cvec[ci][cj][j]=CMPLX(1.0,0.0) #<double complex>1.0   #mp_ctx.mpc(0,tmp)
    #for j from 0<=j<Ql:
    for j in range(Qv[0][1]): 
        Xm[j]=Xm[j]*twopi
    if charvec<>NULL:
        sig_free(charvec)
    if normalizers<>NULL:
        for i in range(nc):
            if normalizers[i]<>NULL:
                sig_free(normalizers[i])
        sig_free(normalizers)
    if cusp_maps<>NULL:
        for i in range(nv):
            if cusp_maps[i]<>NULL:
                sig_free(cusp_maps[i])
        sig_free(cusp_maps)
    if vertex_maps<>NULL:
        for i in range(nv):
            if vertex_maps[i]<>NULL:
                sig_free(vertex_maps[i])
        sig_free(vertex_maps)
    if vertex_widths<>NULL:
        sig_free(vertex_widths)
    if vertex_cusp<>NULL:
        sig_free(vertex_cusp)
    if widths<>NULL:
        sig_free(widths)
    if is_Gamma0==1:
        if reps<>NULL:
            for j in range(nreps):
                if reps[j]<>NULL:
                    if reps[j][0]<>NULL:
                        sig_free(reps[j][0])
                    if reps[j][1]<>NULL:
                        sig_free(reps[j][1])
                    sig_free(reps[j])
            sig_free(reps)
    return 0
#    sig_off()
@cython.cdivision(True)
cdef int pullback_pts_cplx_dp(S,int Qs,int Qf,double Y,double *Xm,double *** Xpb,double*** Ypb,double complex ***Cvec):
    r""" Computes a whole array of pullbacked points using double precision

    INPUT:

    - ``S``  -- AutomorphicFormSpace
    - ``Qs`` -- integer
    - ``Qf`` -- integer
    - ``Y``  -- real >0
    - ``weight`` -- (optional) real
    - ``holo``   -- (False) logical

    
    OUTPUT:

    
    - ``pb`` -- dictionary with entries:
       - 'xm'   -- real[Qf-Qs+1]  : x_m=2*pi*(1-2*m)/2(Qf-Qs+1)
       - 'xpb'  -- real[0:nc,0:nc,Qf-Qs+1]  : real part of pullback of x_m+iY
       - 'ypb'  -- real[0:nc,0:nc,Qf-Qs+1]  : imag. part of pullback of x_m+iY
       - 'cvec' -- complex[0:nc,0:nc,Qf-Qs+1] : mult. factor

    RETURNS::


    - 0 -- if everything is ok
    - 1 -- if Y is too large

    EXAMPLES::


        sage: G=MySubgroup(Gamma0(2))
        sage: [Xm,Xpb,Ypb,Cv]=pullback_pts(G,1,3,0.1)

    Note that we need to have $Y<Y_{0}$ so that all pullbacked points
    are higher up than Y.

        sage: pullback_pts(G,1,10,0.44)
        Traceback (most recent call last):
        ...
        ArithmeticError: Need smaller value of Y

    
    """
#    sig_on()
    G = S.group()
    multiplier = S.multiplier()
    #character = S.character
    holo = S.is_holomorphic()
    cdef double weight
    weight = S.weight()
    if weight<>0 or not multiplier.is_trivial():
        non_trivial=True
    else:
        non_trivial=False
    trivial_mult = multiplier.is_trivial()
    cdef double pi=3.1415926535897932384626433
    cdef double twopi,twoQl
    cdef double complex twopii,mp_i,m1,m2,m3,tmp,v
    cdef double ep0=1E-10
    cdef double tmpabs,tmparg,tmparg1,tmparg2,tmparg3
    twopi=2.0*pi
    twopii = twopi* _Complex_I
    #=complex(0,twopi)
    if(holo and weight==0):
        raise Exception,"Must support a weight for holomorphic forms!"
    cdef double fourQ
    cdef int wi,wj
    cdef int A[4]
    cdef int pba,pbb,pbc,pbd
    cdef list U
    #    AA=[1,0,0,1]
#    Tj=[1,0,0,1]
    cdef double swi,swj,ar,br,cr,dr
    cdef int ci,cj
    cdef int verbose=S._verbose
    cdef double x,y,x1,y1,x2,y2,x3,y3,xx,yy
    cdef int** normalizers=NULL
    cdef int** cusp_maps=NULL
    cdef int** vertex_maps=NULL
    cdef int* vertex_cusp=NULL
    cdef double *vertex_widths=NULL
    cdef int a,b,c,d,i,j,Ql,vi,vj
    cdef double wi_d,wj_d
    cdef ComplexNumber ctmp
    cdef int nc=G._ncusps
    cdef int nv=G._nvertices
    cdef int use_int=1
    cdef int is_Gamma0,nreps
    is_Gamma0=<int>G._is_Gamma0
    cdef int*** reps=NULL
    cdef int N
    if verbose>1:
        print "In pullback_pts_cplx_dp!"
    if is_Gamma0==1:
        nreps=G._index
        N=G._level
        reps= <int ***> sig_malloc(sizeof(int**) * nreps)
        for j in range(nreps):
            reps[j]=NULL
            reps[j]=<int **> sig_malloc(sizeof(int*) * 2)
            if reps[j]==NULL:
                raise MemoryError
            reps[j][0]=NULL;reps[j][1]=NULL
            reps[j][0]=<int *> sig_malloc(sizeof(int) * 2)
            if reps[j][0]==NULL:
                raise MemoryError
            reps[j][1]=<int *> sig_malloc(sizeof(int) * 2)
            if reps[j][1]==NULL:
                raise MemoryError
            reps[j][0][0]=G.coset_reps()[j][0]
            reps[j][0][1]=G.coset_reps()[j][1]
            reps[j][1][0]=G.coset_reps()[j][2]
            reps[j][1][1]=G.coset_reps()[j][3]
#            reps[j][0][0]=G._coset_reps_v0[j][0]
#            reps[j][0][1]=G._coset_reps_v0[j][1]
#            reps[j][1][0]=G._coset_reps_v0[j][2]
#            reps[j][1][1]=G._coset_reps_v0[j][3]

            
    if verbose>1:
        print "nc=",nc
    if not isinstance(G._cusp_data[0]['width'],(int,Integer)):
        use_int=0
    normalizers=<int**>sig_malloc(sizeof(int*)*nc)
    if normalizers==NULL: raise MemoryError
    for i in range(nc):
        normalizers[i]=<int*>sig_malloc(sizeof(int)*4)
        a,b,c,d=G._cusp_data[i]['normalizer']
        normalizers[i][0]=a
        normalizers[i][1]=b
        normalizers[i][2]=c
        normalizers[i][3]=d
        if verbose>2:
            print "Normalizer[",i,"]=",a,b,c,d
    vertex_maps=<int**>sig_malloc(sizeof(int*)*nv)
    if vertex_maps==NULL: raise MemoryError
    cusp_maps=<int**>sig_malloc(sizeof(int*)*nv)
    if cusp_maps==NULL: raise MemoryError
    vertex_widths=<double*>sig_malloc(sizeof(double)*nv)
    if vertex_widths==NULL: raise MemoryError
    vertex_cusp=<int*>sig_malloc(sizeof(int)*nv)    
    if vertex_cusp==NULL: raise MemoryError
    for i in range(nv):
        cusp_maps[i]=<int*>sig_malloc(sizeof(int)*4)
        cusp_maps[i][0]=<int>G._cusp_maps[i][0]
        cusp_maps[i][1]=<int>G._cusp_maps[i][1]
        cusp_maps[i][2]=<int>G._cusp_maps[i][2]
        cusp_maps[i][3]=<int>G._cusp_maps[i][3]
        vertex_maps[i]=<int*>sig_malloc(sizeof(int)*4)
        vertex_maps[i][0]=<int>G._vertex_maps[i].a()
        vertex_maps[i][1]=<int>G._vertex_maps[i].b()
        vertex_maps[i][2]=<int>G._vertex_maps[i].c()
        vertex_maps[i][3]=<int>G._vertex_maps[i].d()
        vertex_widths[i]=<double>G._vertex_widths[i]
        vertex_cusp[i]=<int>G._vertex_data[i]['cusp']
    cdef double* widths
    widths=<double*>sig_malloc(sizeof(double)*nc)
    for i in range(nc):
        if verbose>3:
            print "width[",i,"]=",G._cusp_data[i]['width'],type(G._cusp_data[i]['width'])
        widths[i]=<double>G._cusp_data[i]['width']
    cdef int dir_char=0
    cdef double complex *charvec=NULL
    cdef int modulus
    if verbose>0:
        print "Here1.5"
    if not multiplier.is_trivial():
        dir_char=1
        modulus = multiplier._character.modulus()
        charvec=<double complex*>sig_malloc(sizeof(double complex)*modulus)
        for i in range(modulus):
            mc = multiplier._character(i)
            if hasattr(mc,"complex_embedding"):
                charvec[i]=<double complex> multiplier._character(i).complex_embedding()
            else:
                charvec[i]=<double complex>CC(multiplier._character(i))
            if verbose>1:
                print "chi(",i,")=",charvec[i]
    #if verbose>0:
    #    print "Here2"
    Ql=Qf-Qs+1
    if Qs<0:
        fourQ=<double>2*(Qf-Qs+1)
    else:
        fourQ=<double>(4*Qf)
    for j in range(Qs,Qf+1):
        Xm[j-Qs]=<double>(2*j-1)/fourQ        
    if verbose>0:
        print "use_int=",use_int
    res={}
                
    for ci in range(nc): #in G._cusps:
        # if ci==0:
        #     verbose=3
        # else:
        #     verbose=0
        #ci = G._cusps[i]
        if verbose>2:
            print "------------------------------------------"
            print "ci=",ci
        #cii=G._cusps.index(ci)
        if use_int==1:
            wi = <int>widths[ci] #G._cusp_data[ci]['width'] #(ci)
            swi=sqrt(<double>wi)
        else:
            wi_d=widths[ci] #<double>G._cusp_data[ci]['width'] #(ci)
            swi=sqrt(wi_d)
        for j in range(Ql): #Ql, 1-Q,Q):
            if verbose>0:
                print "j=",j
            x=Xm[j]; y=Y
            
            a=normalizers[ci][0]; b=normalizers[ci][1]
            c=normalizers[ci][2]; d=normalizers[ci][3]
            if verbose>2:
                print "normalizer[",ci,"]=",a,b,c,d
                print "width = ",wi
                print "x0,y0=",x,y
            if use_int==1:
                _normalize_point_to_cusp_dp(&x,&y,a,b,c,d,wi)
            else:
                _normalize_point_to_cusp_real_dp(&x,&y,a,b,c,d,wi_d)
            #[x,y]   = normalize_point_to_cusp_dp(G,ci,Xm[j],Y)
            if verbose>2:
                print "N(z0)=",x,y
            if is_Gamma0 == 1:
                #x1,y1,pba,pbb,pbc,pbd =  G.pullback(x,y,ret_mat=0)
                _pullback_to_Gamma0N_dp(reps,nreps,N,&x,&y,
                                        &pba,&pbb,&pbc,&pbd,0)
                x1=x; y1=y
            else:
#                x1,y1,pba,pbb,pbc,pbd =  G.pullback(x,y,ret_mat=0)
                x1,y1,pba,pbb,pbc,pbd =  pullback_general_group_dp(G,x,y,ret_mat=1)
                if verbose>2:
                    print "Pbz=",x1,y1
            ## Want to replace this with a cpdef'd function
            #cj,vj= G.closest_cusp(x1,y1,vertex=1)
            if verbose>2:
                print "x1,y1=",x1,y1
            vj = closest_vertex_dp_c(nv,vertex_maps,vertex_widths,&x1,&y1)
            if verbose>2:
                print "vj=",vj
            cj = vertex_cusp[vj]
            if pba<=0 and pbb<=0 and pbc<=0 and pbd<=0:
                pba=-pba; pbb=-pbb; pbc=-pbc; pbd=-pbd
            if verbose>2:
                print "x0,y0=",x,y
                print "j=",j+Qs
                print "Xm,Y=",Xm[j],Y
                print "x,y=",x,y
                print "Pbmap=",pba,pbb,pbc,pbd
                print "x1,y1=",x1,y1
                print "cj,vj=",cj,vj
            #cjj=G._cusps.index(cj)
            if use_int==1:
                wj=<int>widths[cj] #G._cusp_data[cj]['width']
                swj=sqrt(<double>wj) 
            else:
                wj_d=widths[cj] #G._cusp_data[cj]['width']
                swj=sqrt(wj_d)
            a = cusp_maps[vj][0]   #a,b,c,d=G._vertex_data[cj]['cusp_map'] #[v]
            b = cusp_maps[vj][1]
            c = cusp_maps[vj][2]
            d = cusp_maps[vj][3]
            if verbose>2:
                if cj==1:
                    print ">>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<"
                print "cj=",cj
                print "U1(",vj,")=",a,b,c,d
                print "x00,y00=",x1,y1
            x2=x1; y2=y1
            if a<>1 or b<>0 or c<>0 or d<>1:
                #if verbose>2:
                #    print "here1 type=",type(a)
                #[x2,y2]=_apply_sl2z_map_dp(x1,y1,U[0],U[1],U[2],U[3])
                _apply_sl2z_map_dp(&x2,&y2,a,b,c,d)
            if verbose>2:
                print "x01,y01=",x2,y2
            #    print "cj=",cj
            #x2=x1; y2=y1
            x3=x2; y3=y2
            a=normalizers[cj][3]; b=-normalizers[cj][1]
            c=-normalizers[cj][2]; d=normalizers[cj][0]
            if verbose>2:
                print "Normalizer: a,b,c,d=",a,b,c,d
            #_normalize_point_to_cusp_dp(&x3,&y3,a,b,c,d,wj,1)
            if use_int==1:
                _normalize_point_to_cusp_dp(&x3,&y3,a,b,c,d,wj,1)
            else:
                _normalize_point_to_cusp_real_dp(&x3,&y3,a,b,c,d,wj_d,1)
            if verbose>2:
                print "width=",wj
                print "x02,y02=",x3,y3
                print "ci,cj,j=",ci,cj,j
            #[x3,y3] = normalize_point_to_cusp_dp(G,cj,x2,y2,inv=1)
            Xpb[ci][cj][j]=x3*twopi
            if verbose>2:
                print "Xpb[{0}][{1}][{2}]={3}".format(ci,cj,j,Xpb[ci][cj][j])
            # Recall that Ypb must be greater than Y otherwise we need a better Y
            if y3>=Y+ep0:
                Ypb[ci][cj][j]=y3
                if verbose>2:
                    print "Ypb[{0}][{1}][{2}]={3}".format(ci,cj,j,Ypb[ci][cj][j])                
            else:
                if verbose>0:
                    print "ci,cj=",ci,cj
                    print "Tj=",pba,pbb,pbc,pbd
                    print "wj=",wj
                    print "wj_d=",wj_d
                    print "x,y=",x,y
                    print "x1,y1=",x1,y1
                    print "x2,y2=",x2,y2
                    print "x3,y3=",x3,y3
                    print "Normalizer=",a,b,c,d
                    print "ep0=",ep0
                return 1
                #raise ArithmeticError,"Need smaller value of Y. Got:{0}".format(Y) 
            if verbose>2:
                print "ci,cj=",ci,cj
                print "Tj=",pba,pbb,pbc,pbd
                print "Xm=",Xm[j]
                print "x,y=",x,"\n",y
                print "x1,y1=",x1,"\n",y1
                print "x2,y2=",x2,"\n",y2
                print "x3,y3=",x3,"\n",y3
                print "Xpb=",Xpb[ci][cj][j]
                print "Ypb=",Ypb[ci][cj][j]
            # We also get the multiplier if we need it
            if non_trivial:
                if weight<>0:
                    if verbose>2:
                        print "here holo=",holo," weight=",weight
                        print "A=",A[0],A[1],A[2],A[3]
                    cr = <double>(normalizers[ci][2])*swi
                    dr = <double>(normalizers[ci][3])/swi
                    if verbose>2:
                        print "cr,dr=",cr,dr
                    j_fak_dp_c(cr,dr,Xm[j],Y,weight,holo,&tmpabs,&tmparg1)
                    m1=cexp(-_Complex_I*tmparg1)
                    if verbose>2:
                        print "-tmparg1=",-tmparg1
                        print "m1=",m1
                    cr = <double>(normalizers[cj][2])*swj
                    dr = <double>(normalizers[cj][3])/swj
                    j_fak_dp_c(cr,dr,x3,y3,weight,holo,&tmpabs,&tmparg2)
                    m2=cexp(_Complex_I*tmparg2)
                    if verbose>2:
                        print "tmparg2=",tmparg2
                    #l=mat_mul_lit(cusp_maps[vj],Tj)
                    #A[0]=l[0]; A[1]=l[1]; A[2]=l[2]; A[3]=l[3]
                    _mat_mul_list(cusp_maps[vj][0],cusp_maps[vj][1],cusp_maps[vj][2],cusp_maps[vj][3],pba,pbb,pbc,pbd,A)
                    #A=A**-1
                    cr=<double>(-A[2]); dr=<double>(A[0])
                    j_fak_dp_c(cr,dr,x2,y2,weight,holo,&tmpabs,&tmparg3)
                    m3=cexp(_Complex_I*tmparg3)
                    tmp=cexp(_Complex_I*(tmparg2-tmparg1+tmparg3))
                    if verbose>2:
                        print "tmparg3=",tmparg3
                        print "A[2],A[3]=",A[2],A[3]
                        print "A[3]=",A[3],"% modulus:",vi
                        print "m1,m2,m3=",m1,m2,m3
                else:
                    _mat_mul_list(cusp_maps[vj][0],cusp_maps[vj][1],cusp_maps[vj][2],cusp_maps[vj][3],pba,pbb,pbc,pbd,A)
                    tmp=CMPLX(1.0,0.0) #<double complex>1.0
                if dir_char:
                    vi=A[3] % modulus
                    if vi<0:
                        vi=vi+modulus
                    v=charvec[vi]
                    if verbose>2:
                        print "A[3]=",A[3],"% modulus:",vi
                        print "m1,m2,m3=",m1,m2,m3
                        print "tmp=",tmp
                        print "mult=",v 

                    #v = <double complex> ctmp.real()+_Complex_I*ctmp.imag() #CC(multiplier(AA)) #character(AA[1,1])a
                    tmp = tmp/v
                elif not trivial_mult:
                    if verbose>1:
                        print "mult=",multiplier([A[0],A[1],A[2],A[3]])
                    ctmp=multiplier([A[0],A[1],A[2],A[3]]).complex_embedding()
                    v = <double complex> ctmp.real()+_Complex_I*ctmp.imag() #CC(multiplier(AA)) #character(AA[1,1])a
                    tmp = tmp/v
                Cvec[ci][cj][j]=tmp # mp_ctx.mpc(0,tmp)
            else:
                Cvec[ci][cj][j]=CMPLX(1.0,0.0) #<double complex>1.0   #mp_ctx.mpc(0,tmp)
    for j in range(Ql):
        Xm[j]=Xm[j]*twopi
    if charvec<>NULL:
        sig_free(charvec)
    if normalizers<>NULL:
        for i in range(nc):
            if normalizers[i]<>NULL:
                sig_free(normalizers[i])
        sig_free(normalizers)
    if cusp_maps<>NULL:
        for i in range(nv):
            if cusp_maps[i]<>NULL:
                sig_free(cusp_maps[i])
        sig_free(cusp_maps)
    if vertex_maps<>NULL:
        for i in range(nv):
            if vertex_maps[i]<>NULL:
                sig_free(vertex_maps[i])
        sig_free(vertex_maps)
    if vertex_widths<>NULL:
        sig_free(vertex_widths)
    if vertex_cusp<>NULL:
        sig_free(vertex_cusp)
    if widths<>NULL:
        sig_free(widths)
    if is_Gamma0==1:
        if reps<>NULL:
            for j in range(nreps):
                if reps[j]<>NULL:
                    if reps[j][0]<>NULL:
                        sig_free(reps[j][0])
                    if reps[j][1]<>NULL:
                        sig_free(reps[j][1])
                    sig_free(reps[j])
            sig_free(reps)
    return 0
#    sig_off()


@cython.cdivision(True)
cdef int pullback_pts_real_dp(S,int Qs,int Qf,double Y,double *Xm,double *** Xpb,double*** Ypb,double ***Cvec):
    r"""
    Computes a whole array of pullbacked points using double precision
    The multiplier/Character is supposed to be real in this case.
    Wa also suppose weight = 0 
    NOTE: For docs. See pullback_pts_dp

    
    """
    sig_on()
    G = S.group()
    multiplier = S.multiplier()
    #character = S.character
    cdef double weight
    weight = S.weight()
    holo = S.is_holomorphic()
    if weight<>0 or not S._use_real:
        raise ValueError,"This (real) algorithm should not be used in this situation!"
    cdef int non_trivial=0
    if not multiplier.is_trivial():
        non_trivial=1
    cdef double pi=3.1415926535897932384626433
    cdef double twopi=6.28318530717958647692528676656
    cdef double twoQl
    #twopi=2.0*pi
    #=complex(0,twopi)
    if(holo and weight==0):
        raise Exception,"Must support a weight for holomorphic forms!"
    cdef double fourQ
    cdef int wi,wj
    cdef int A[4]
    cdef int pba,pbb,pbc,pbd
    cdef list U
#    AA=[1,0,0,1]
#    Tj=[1,0,0,1]
    cdef double swi,swj,ar,br,cr,dr,xx,yy
    cdef double tmp,v
    cdef int ci,cj
    cdef int verbose=S._verbose
    cdef double x,y,x1,y1,x2,y2,x3,y3
    cdef int** normalizers=NULL
    cdef int** cusp_maps=NULL
    cdef int a,b,c,d,i,j,Ql,vi,vj
    cdef double wi_d,wj_d
    cdef double ep0=1E-10
    #cdef RComplexNumber ctmp
    cdef int nc=G._ncusps
    cdef int nv=G._nvertices
    cdef int use_int=1
    cdef int nreps,N,is_Gamma0
    cdef int*** reps=NULL
    is_Gamma0=G._is_Gamma0
    if is_Gamma0==1:
        nreps=G._index
        N=G._level
        reps= <int ***> sig_malloc(sizeof(int**) * nreps)
        for j in range(nreps):
            reps[j]=NULL
            reps[j]=<int **> sig_malloc(sizeof(int*) * 2)
            if reps[j]==NULL:
                raise MemoryError
            reps[j][0]=NULL;reps[j][1]=NULL
            reps[j][0]=<int *> sig_malloc(sizeof(int) * 2)
            if reps[j][0]==NULL:
                raise MemoryError
            reps[j][1]=<int *> sig_malloc(sizeof(int) * 2)
            if reps[j][1]==NULL:
                raise MemoryError
            reps[j][0][0]=G._coset_reps_v0[j][0]
            reps[j][0][1]=G._coset_reps_v0[j][1]
            reps[j][1][0]=G._coset_reps_v0[j][2]
            reps[j][1][1]=G._coset_reps_v0[j][3]
    #if verbose>0:
    #    print "Here1"
    if not isinstance(G._cusp_data[0]['width'],(int,Integer)):
        use_int=0
    normalizers=<int**>sig_malloc(sizeof(int*)*nc)
    if normalizers==NULL: raise MemoryError
    for i in range(nc):
        normalizers[i]=<int*>sig_malloc(sizeof(int)*4)
        [a,b,c,d]=G._cusp_data[i]['normalizer']
        normalizers[i][0]=a
        normalizers[i][1]=b
        normalizers[i][2]=c
        normalizers[i][3]=d
    cusp_maps=<int**>sig_malloc(sizeof(int*)*nv)
    #if verbose>0:
    #    print "Here2 nc=",nc,type(nc)
    if cusp_maps==NULL: raise MemoryError
    for i in range(nv):
        cusp_maps[i]=<int*>sig_malloc(sizeof(int)*4)
        cusp_maps[i][0]=G._cusp_maps[i][0]
        cusp_maps[i][1]=G._cusp_maps[i][1]
        cusp_maps[i][2]=G._cusp_maps[i][2]
        cusp_maps[i][3]=G._cusp_maps[i][3]
    cdef double* widths
    widths=<double*>sig_malloc(sizeof(double)*nc)
    for i in range(nc):
        if verbose>0:
            print "width[",i,"]=",G._cusp_data[i]['width'],type(G._cusp_data[i]['width'])
        widths[i]=<double>G._cusp_data[i]['width']
    cdef int dir_char=0
    cdef double *charvec=NULL
    cdef int modulus
    #if verbose>0:
    #    print "Here1.5"
    if non_trivial:
        modulus = multiplier._character.modulus()
        charvec=<double*>sig_malloc(sizeof(double)*modulus)
        for i in range(modulus):
            charvec[i]=<double> RR(multiplier._character(i))
            if verbose>1:
                print "chi(",i,")=",charvec[i]
    #if verbose>0:
    #    print "Here2"
    Ql=Qf-Qs+1
    if Qs<0:
        fourQ=<double>2*(Qf-Qs+1)
    else:
        fourQ=<double>(4*Qf)
    for j in range(Qs,Qf+1):
        Xm[j-Qs]=<double>(2*j-1)/fourQ        
    if verbose>2:
        print "use_int=",use_int
    for ci in range(nc): #in G._cusps:
        #ci = G._cusps[i]
        if verbose>2:
            print "ci=",ci
        #cii=G._cusps.index(ci)
        if use_int==1:
            wi = <int>widths[ci] #G._cusp_data[ci]['width'] #(ci)
            swi=sqrt(<double>wi)
        else:
            wi_d=widths[ci] #<double>G._cusp_data[ci]['width'] #(ci)
            swi=sqrt(wi_d)
        for j in range(Ql): #1-Q,Q):
            x=Xm[j]; y=Y
            a=normalizers[ci][0]
            b=normalizers[ci][1]
            c=normalizers[ci][2]
            d=normalizers[ci][3]
            if use_int==1:
                _normalize_point_to_cusp_dp(&x,&y,a,b,c,d,wi)
            else:
                _normalize_point_to_cusp_real_dp(&x,&y,a,b,c,d,wi_d)
            #[x,y]   = normalize_point_to_cusp_dp(G,ci,Xm[j],Y)
            if verbose>2:
                print "j=",j
                print "Xm,Y=",Xm[j],Y
                print "x,y=",x,y
                print "a,b,c,d,wi=",a,b,c,d,wi
            if is_Gamma0==1:
                #x1,y1,pba,pbb,pbc,pbd =  G.pullback(x,y,ret_mat=0)
                _pullback_to_Gamma0N_dp(reps,nreps,N,&x,&y,
                                        &pba,&pbb,&pbc,&pbd,0)
                x1=x; y1=y
            else:
                x1,y1,pba,pbb,pbc,pbd =  G.pullback(x,y,ret_mat=0)
            if verbose>2:
                print "x1,y1=",x1,y1
            cj,vj= G.closest_cusp(x1,y1,vertex=1)
            #vj = closest_vertex(vertex_maps,vertex_widths,nv,x1,y1)
            #scj = vertex_cusp[vj]
            if verbose>2:
                print "cj,vj=",cj,vj
            #cjj=G._cusps.index(cj)
            if use_int==1:
                wj=<int>widths[cj] #G._cusp_data[cj]['width']
                swj=sqrt(<double>wj) 
            else:
                wj_d=widths[cj] #G._cusp_data[cj]['width']
                swj=sqrt(wj_d)
            a = cusp_maps[vj][0]   #a,b,c,d=G._vertex_data[cj]['cusp_map'] #[v]
            b = cusp_maps[vj][1]
            c = cusp_maps[vj][2]
            d = cusp_maps[vj][3]
            x2=x1; y2=y1
            if verbose>3:
                print "cusp_map[",vj,"]=",a,b,c,d
            if a<>1 or b<>0 or c<>0 or d<>1:
                _apply_sl2z_map_dp(&x2,&y2,a,b,c,d)
            x3=x2; y3=y2
            a=normalizers[cj][3]; b=-normalizers[cj][1]
            c=-normalizers[cj][2]; d=normalizers[cj][0]
            if use_int==1:
                _normalize_point_to_cusp_dp(&x3,&y3,a,b,c,d,wj,1)
            else:
                _normalize_point_to_cusp_real_dp(&x3,&y3,a,b,c,d,wj_d,1)
            #[x3,y3] = normalize_point_to_cusp_dp(G,cj,x2,y2,inv=1)
            Xpb[ci][cj][j]=x3*twopi
            # Recall that Ypb must be greater than Y otherwise we need a better Y
            if y3>=Y+ep0:
                Ypb[ci][cj][j]=y3
            else:
                if verbose>0:
                    print "ci,cj=",ci,cj
                    print "Tj=",pba,pbb,pbc,pbd
                    print "wj=",wj
                    print "wj_d=",wj_d
                    print "x,y=",x,y
                    print "x1,y1=",x1,y1
                    print "x2,y2=",x2,y2
                    print "x3,y3=",x3,y3
                    print "Normalizer=",a,b,c,d
                return 1
            #raise ArithmeticError,"Need smaller value of Y. Got: y3={0} and Y-ep={1}".format(y3,Y-ep0) 
            if verbose>2:
                print "ci,cj=",ci,cj
                print "Tj=",pba,pbb,pbc,pbd
                print "Xm=",Xm[j]
                print "x,y=",x,"\n",y
                print "x1,y1=",x1,"\n",y1
                print "x2,y2=",x2,"\n",y2
                print "x3,y3=",x3,"\n",y3
                print "Xpb=",Xpb[ci][cj][j]
                print "Ypb=",Ypb[ci][cj][j]
            # We also get the multiplier if we need it
            Cvec[ci][cj][j]=<double>1.0   
            if non_trivial:
                _mat_mul_list(cusp_maps[vj][0],cusp_maps[vj][1],cusp_maps[vj][2],cusp_maps[vj][3],pba,pbb,pbc,pbd,A)
                vi=A[3] % modulus
                if vi<0:
                    vi=vi+modulus
                v=charvec[vi]
                if verbose>1:
                    print "A[3]=",A[3],"% modulus:",vi
                    print "mult=",v #multiplier([A[0],A[1],A[2],A[3]])
                Cvec[ci][cj][j]=Cvec[ci][cj][j]/v
                

    for j in range(Ql):
        Xm[j]=Xm[j]*twopi
    if non_trivial:
        if charvec<>NULL:
            sig_free(charvec)
    if normalizers<>NULL:
        sig_free(normalizers)
    if cusp_maps<>NULL:
        sig_free(cusp_maps)
    if widths<>NULL:
        sig_free(widths)
    if is_Gamma0==1:
        if reps<>NULL:
            for j in range(nreps):
                if reps[j]<>NULL:
                    if reps[j][0]<>NULL:
                        sig_free(reps[j][0])
                    if reps[j][1]<>NULL:
                        sig_free(reps[j][1])
                    sig_free(reps[j])
            sig_free(reps)
    sig_off()
    return 0


def pullback_pts_fp(S,Qs,Qf,Y,weight=0,holo=False):
    r""" Computes a whole array of pullbacked points using floating point precision

    INPUT:

    - ``S``  -- AutomorphicFormSpace
    - ``Qs`` -- integer
    - ``Qf`` -- integer
    - ``Y``  -- real >0
    - ``weight`` -- (optional) real
    - ``holo``   -- (False) logical

    
    OUTPUT:

    - ``pb`` -- dictionary with entries:
       - 'xm'   -- real[Qf-Qs+1]  : x_m=2*pi*(1-2*m)/2(Qf-Qs+1)
       - 'xpb'  -- real[0:nc,0:nc,Qf-Qs+1]  : real part of pullback of x_m+iY
       - 'ypb'  -- real[0:nc,0:nc,Qf-Qs+1]  : imag. part of pullback of x_m+iY
       - 'cvec' -- complex[0:nc,0:nc,Qf-Qs+1] : mult. factor

    EXAMPLES::


        sage: G=MySubgroup(Gamma0(2))
        sage: [Xm,Xpb,Ypb,Cv]=pullback_pts(G,1,3,0.1)

    Note that we need to have $Y<Y_{0}$ so that all pullbacked points
    are higher up than Y.

        sage: pullback_pts(G,1,10,0.44)
        Traceback (most recent call last):
        ...
        ArithmeticError: Need smaller value of Y

    
    """
    G = S.group()
    multiplier = S.multiplier()
    #character = S.character
    holo = S.is_holomorphic()
    if weight<>0 or not multiplier.is_trivial():
        non_trivial=True
    else:
        non_trivial=False

    trivial_mult = multiplier.is_trivial()
    cdef float pi=3.1415926535897932384626433
    cdef float twopi,twoQl
    cdef complex twopii,mp_i
    cdef float ep0=1E-10
    twopi=2.0*pi
    twopii = twopi* _Complex_I
    #=complex(0,twopi)
    Xm=dict()  
    Xpb=dict() 
    Ypb=dict() 
    Cvec=dict()
    if(holo and weight==0):
        raise Exception,"Must support a weight for holomorphic forms!"
    twoQl=<double>(2*(Qf+1-Qs))
    if(Qs<0):
        for j in range(Qs,Qf+1):
            Xm[j]=<double>(2*j-1)/twoQl        
    else:
        for j in range(Qs,Qf+1):
            Xm[j]=<double>(2*j-1)/twoQl        

    cdef int cii,wi
    cdef float swi
    cdef complex tmp
    cdef tuple ci
    cdef int verbose=S._verbose
    for ci in G._cusps:
        cii=G._cusps.index(ci)
        wi = G._cusp_data[ci]['width'] #(ci)
        swi=sqrt(<float>wi)
        for j in range(Qs,Qf+1):
            [x,y]   = normalize_point_to_cusp_dp(G,ci,Xm[j],Y)
            if verbose>1:
                print "x,y=",x,y
            [x1,y1,Tj] =  G.pullback(x,y)
            if verbose>1:
                print "x1,y1=",x1,y1
            #vi = G.closest_vertex(x1,y1)
            cj= G.closest_cusp(x1,y1) #G._vertex_data[G._vertices[vi]]['cusp'] #cusp_representative[v]
            #v = G.closest_vertex(x1,y1)
            #cj= G._vertex_data[v]['cusp'] #cusp_representative[v]
            cjj=G._cusps.index(cj)
            swj=sqrt(<float>(G._cusp_data[cj]['width']))
            U = G._vertex_data[v]['cusp_map'] #[v]
            if U<>SL2Z_elt(1,0,0,1):
                [x2,y2] = apply_sl2z_map_mpfr(x1,y1,U[0],U[1],U[2],U[3])
            else:
                x2=x1; y2=y1;
            [x3,y3] = normalize_point_to_cusp_dp(G,cj,x2,y2,inv=1)
            Xpb[cii,cjj,j]=x3*twopi
            # Recall that Ypb must be greater than Y otherwise we need a better Y
            if y3>=Y-ep0:
                Ypb[cii,cjj,j]=y3
            else:
                if verbose>0:
                    print "ci,cj=",ci,cj
                    print "Tj=",Tj
                    print "x,y=",x,y
                    print "x1,y1=",x1,y1
                    print "x2,y2=",x2,y2
                    print "x3,y3=",x3,y3
                raise ArithmeticError,"Need smaller value of Y. Got:{0}".format(Y) 
            if verbose>1:
                print "ci,cj=",ci,cj
                print "Xm=",Xm[j]
                print "x,y=",x,"\n",y
                print "x1,y1=",x1,"\n",y1
                print "x2,y2=",x2,"\n",y2
                print "x3,y3=",x3,"\n",y3
                print "Xpb=",Xpb[cii,cjj,j]
                print "Ypb=",Ypb[cii,cjj,j]
            # We also get the multiplier if we need it
            if non_trivial:
                if weight<>0:
                    c = <float>(G._cusp_data[ci]['normalizer'][0])*swi
                    d = <float>(G._cusp_data[ci]['normalizer'][1])/swi
                    m1=j_fak_dp(c,d,Xm[j],Y,-weight,holo)
                    c = <float>(G._cusp_data[cj]['normalizer'][2])*swj
                    d = <float>(G._cusp_data[cj]['normalizer'][3])/swj
                    m2=j_fak_dp(c,d,x3,y3,weight,holo)
                    A=(U*Tj)
                    #A=A**-1
                    c=<float>(-A[1,0]); d=<float>(A[0,0])
                    m3=j_fak_dp(c,d,x2,y2,weight,holo)
                    tmp=m1*m2*m3
                else:
                    A=(U*Tj)
                    tmp=CMPLXF(1.0,0.0) #<complex>1.0
                if not trivial_mult:
                    AA=A**-1
                    if verbose>1:
                        print "mult=",multiplier(AA),type(multiplier(AA))
                    mA = multiplier(AA)
                    if hasattr(mA,"complex_embedding"):
                        ctmp = mA.complex_embedding()
                    else:
                        ctmp = CC(mA)
                    #try:
                    #    ctmp=multiplier(AA).complex_embedding()
                    #except:
                    #    ctmp=CC(AA)
                    v = <complex> ctmp.real()+_Complex_I*ctmp.imag() #CC(multiplier(AA)) #character(AA[1,1])a
                    tmp = tmp*v
                Cvec[cii,cjj,j]=tmp # mp_ctx.mpc(0,tmp)
            else:
                Cvec[cii,cjj,j]=CMPLXF(1.0,0.0) #<complex>1.0   #mp_ctx.mpc(0,tmp)

    for j in range(Qs,Qf+1): #1-Q,Q):
        Xm[j]=Xm[j]*twopi
    pb=dict()
    pb['xm']=Xm; pb['xpb']=Xpb; pb['ypb']=Ypb; pb['cvec']=Cvec
    return pb



cpdef pullback_pts_mpc(S,int Qs,int Qf,RealNumber Y,deb=False):
    r""" Computes a whole array of pullbacked points-
         using MPFR/MPC types.
    INPUT:

    - ``S``  -- AutomorphicFormSpace
    - ``Qs`` -- integer
    - ``Qf`` -- integer
    - ``Y``  -- real >0
    - ``deb``    -- (False) logical

    
    OUTPUT:

    - ``pb`` -- dictionary with entries:
       - 'xm'   -- real[Qf-Qs+1]  : x_m=2*pi*(1-2*m)/2(Qf-Qs+1)
       - 'xpb'  -- real[0:nc,0:nc,Qf-Qs+1]  : real part of pullback of x_m+iY
       - 'ypb'  -- real[0:nc,0:nc,Qf-Qs+1]  : imag. part of pullback of x_m+iY
       - 'cvec' -- complex[0:nc,0:nc,Qf-Qs+1] : mult. factor

    EXAMPLES::


        sage: G=MySubgroup(Gamma0(2))
        sage: [Xm,Xpb,Ypb,Cv]=pullback_pts(G,1,3,0.1)   

    Note that we need to have $Y<Y_{0}$ so that all pullbacked points
    are higher up than Y.

        sage: pullback_pts(G,1,10,0.44)
        Traceback (most recent call last):
        ...
        ArithmeticError: Need smaller value of Y

    
    """
    # WE let the input determine the precision
    #print "pullback mpc!"
    cdef int prec=Y.parent().prec()
    CF=ComplexField(prec)
    RF=RealField(prec)
    cdef RealNumber x0,y0,x1,y1,x2,y2,weight
    x0=RealNumber(RF,1)
    y0=RealNumber(RF,1)
    x1=RealNumber(RF,1)
    x2=RealNumber(RF,1)
    y1=RealNumber(RF,1)
    y2=RealNumber(RF,1)
    cdef int a,b,c,d
    cdef RealNumber ar,br,cr,dr
    #weight = RF(wt)
    G = S.group()
    multiplier = S.multiplier()
    #character = S.character
    holo = S.is_holomorphic()
    cdef int ui
    ui = S._unitary_action
    #if holo:
    #    ui=0
    #else:
    #    ui=1
    weight = RF(S.weight())
    if weight<>0 or not multiplier.is_trivial():
        non_trivial=True
    else:
        non_trivial=False
    trivial_mult = multiplier.is_trivial()
    twopi=RF(2)*RF.pi() 
    twopii=CF(0,twopi)
    mp_i=CF(0,1)
    Xm=dict()  #vector(RF,2*Q)
    Xpb=dict()  #vector(RF,2*Q)
    Ypb=dict()  #vector(RF,2*Q)
    Cvec=dict()
    if holo and weight==0 and not S._weak:
        raise Exception,"Must support a non-zero weight for holomorphic forms!"
    verbose = S._verbose
    if verbose>3:
        fp0=open("xm.txt","w")
        fp1=open("xpb.txt","w")
        fp2=open("ypb.txt","w")
        fp3=open("cv.txt","w")
    #twoQl=RF(2*(Qf+1-Qs))
    if Qs<0 :
        Qfak = RF(2*(Qf-Qs+1))
        for j in range(Qs,Qf+1): #1-Q,Q):
            Xm[j]=RF(2*j-1)/Qfak
            if verbose>3:
                s=str(Xm[j])+"\n"
                fp0.write(s)
    else:
        Qfak = RF(4*(Qf-Qs+1))
        for j in range(Qs,Qf+1):
            Xm[j]=RF(2*j-1)/Qfak
            if verbose>3:
                s=str(Xm[j])+"\n"
                fp0.write(s)
    cdef tuple ci,cj
    for ci in G._cusps:
        cii=G._cusps.index(ci)
        if verbose>1:
            print "cusp =",ci
        swi=RF(G._cusp_data[ci]['width']).sqrt()
        for j in range(Qs,Qf+1):
            #if verbose > 0:
            #    print "Y before=",Y
            if cii==0:
                x0=Xm[j]; y0=Y
            else:
                x0,y0   = normalize_point_to_cusp_mpfr(G,ci[0],ci[1],Xm[j],Y,inv=0)
            #if verbose > 0:
            #    print "Y after=",Y
            #[x1,y1,Tj] =  pullback_to_G(S,x,y,mp_ctx=mp_ctx)
             
            x1,y1,a,b,c,d =  pullback_to_Gamma0N_mpfr(G,x0,y0)
            #if verbose > 0:
            #    print "Y after2=",Y
            Tj=SL2Z([a,b,c,d])
            #vi = G.closest_vertex(x1,y1)
            cj= G.closest_cusp(x1,y1) #G._vertex_data[G._vertices[vi]]['cusp'] #cusp_representative[v]
            #v = G.closest_vertex(x1,y1)
            #cj= G._vertex_data[v]['cusp'] 
            cjj=G._cusps.index(cj)
            swj=RF(G._cusp_data[cj]['width']).sqrt()
            #swj=mp_ctx.sqrt(mp_ctx.mpf(G._cusp_width[cj][0]))
            U = G._vertex_data[v]['cusp_map']
            if U<>SL2Z_elt(1,0,0,1):
                #mpfr_set(mpx2,x1.value,rnd_re)
                #mpfr_set(mpy2,y1.value,rnd_re)
                x2=x1; y2=y1
                _apply_sl2z_map_mpfr(x2.value,y2.value,U[0,0],U[0,1],U[1,0],U[1,1])
                #[x2,y2] = apply_sl2_map(x1,y1,U)
            else:
                x2=x1; y2=y1;
            #[x3,y3] = normalize_point_to_cusp(G,x2,y2,cj,inv=True,mp_ctx=mp_ctx)
            if cjj<>0:
                [x3,y3] = normalize_point_to_cusp_mpfr(G,cj[0],cj[1],x2,y2,inv=1)
            else:
                x3=x2; y3=y2
            #Xpb[cii,cjj,j]=mp_ctx.mpc(0,x3*twopi)
            Xpb[cii,cjj,j]=x3*twopi
            # Recall that Ypb must be greater than Y otherwise we need a better Y
            if(y3>Y):
                Ypb[cii,cjj,j]=y3
            else:
                if verbose > 0:
                    print "Y=",Y
                    print "ci,cj=",ci,cj
                    print "Xm=",Xm[j]
                    print "x,y=",x0,"\n",y0
                    print "x1,y1=",x1,"\n",y1
                    print "x2,y2=",x2,"\n",y2
                    print "x3,y3=",x3,"\n",y3
                    print "Xpb=",Xpb[cii,cjj,j]
                    #print "Ypb=",Ypb[cii,cjj,j]
                raise ArithmeticError,"Need smaller value of Y" 

            # We also get the multiplier if we need it
            if(non_trivial):
                if weight<>0:
                    cr = RF(G._cusp_data[ci]['normalizer'][2])*swi
                    dr = RF(G._cusp_data[ci]['normalizer'][3])/swi
                    #print "test:",cr,dr,Xm[j],Y,-weight
                    #print type(c),type(d)
                    m1=j_fak_mpc(cr,dr,Xm[j],Y,-weight,ui)
                    cr = RF(G._cusp_data[cj]['normalizer'][2])*swj
                    dr = RF(G._cusp_data[cj]['normalizer'][3])/swj
                    #print "test:",cr,dr,x3,y3,-weight
                    m2=j_fak_mpc(cr,dr,x3,y3,weight,ui)
                    A=(U*Tj)
                    #A=A**-1
                    cr=RF(-A[1,0]); dr=RF(A[0,0])
                    m3=j_fak_mpc(cr,dr,x2,y2,weight,ui)
                    tmp=m1*m2*m3
                    #print "m1,m2,m3=",m1,m2,m3
                else:
                    A=(U*Tj)
                    tmp=CF(1)
                if not trivial_mult :
                    #if multiplier<>None and not multiplier.is_trivial():
                    AA=A**-1
                    #v = S.character(AA[1,1]).complex_embedding(CF.prec())
                    v = multiplier(AA)
                    #if v<>multiplier(AA):
                    #    v1=v; v2=multiplier(AA)
                    #    raise ArithmeticError," {0} <> {1}".format(v1,v2)
                    #if hasattr(v,'complex_embedding'):
                    #    v=v.complex_embedding(CF.prec())
                    #else:
                    v=CF(v)
                    #print "AA=",AA
                    #print "v=",v
                    tmp = tmp*v
                Cvec[cii,cjj,j]=tmp # mp_ctx.mpc(0,tmp)
            else:
                Cvec[cii,cjj,j]=CF(1) #CF(0,tmp)
            if verbose>3:
                fp1.write(str(x3)+"\n")
                fp2.write(str(y3)+"\n")
                fp3.write(str(tmp)+"\n")

        #        Xmi[j]=mp_ctx.mpc(0,Xm[j]*twopi)
    for j in range(Qs,Qf+1):
        Xm[j]=Xm[j]*twopi
    if verbose>3:
        fp1.close()
        fp2.close()
        fp3.close()
    pb=dict()
    pb['xm']=Xm; pb['xpb']=Xpb; pb['ypb']=Ypb; pb['cvec']=Cvec
    return pb


cpdef pullback_pts_mpc_new(S,int Qs,int Qf,RealNumber Y,deb=False):
    r""" Computes a whole array of pullbacked points-
         using MPFR/MPC types.
    INPUT:

    - ``S``  -- AutomorphicFormSpace
    - ``Qs`` -- integer
    - ``Qf`` -- integer
    - ``Y``  -- real >0
    - ``deb``    -- (False) logical

    
    OUTPUT:

    - ``pb`` -- dictionary with entries:
       - 'xm'   -- real[Qf-Qs+1]  : x_m=2*pi*(1-2*m)/2(Qf-Qs+1)
       - 'xpb'  -- real[0:nc,0:nc,Qf-Qs+1]  : real part of pullback of x_m+iY
       - 'ypb'  -- real[0:nc,0:nc,Qf-Qs+1]  : imag. part of pullback of x_m+iY
       - 'cvec' -- complex[0:nc,0:nc,Qf-Qs+1] : mult. factor

    EXAMPLES::


        sage: G=MySubgroup(Gamma0(2))
        sage: [Xm,Xpb,Ypb,Cv]=pullback_pts(G,1,3,0.1)   

    Note that we need to have $Y<Y_{0}$ so that all pullbacked points
    are higher up than Y.

        sage: pullback_pts(G,1,10,0.44)
        Traceback (most recent call last):
        ...
        ArithmeticError: Need smaller value of Y

    
    """
    # WE let the input determine the precision
    #print "pullback mpc!"
    cdef int prec=Y.parent().prec()
    CF=MPComplexField(prec)
    RF=RealField(prec)
    cdef RealNumber x1,y1,x2,y2,x3,y3,weight
    x1=RealNumber(RF,1)
    x2=RealNumber(RF,1)
    x3=RealNumber(RF,1)
    y1=RealNumber(RF,1)
    y2=RealNumber(RF,1)
    y3=RealNumber(RF,1)
    cdef int a,b,c,d,v0,v1,itmp,i,j,ui
    cdef mpfr_t swj,swi,x0,y0,YY
    cdef Vector_real_mpfr_dense Xm
    cdef RealNumber ar,br,cr,dr,twopi,xm,ep0
    mpfr_init2(swj,prec)
    mpfr_init2(swi,prec)
    mpfr_init2(x0,prec)
    mpfr_init2(y0,prec)
    mpfr_init2(YY,prec)
    mpfr_set(YY,Y.value,rnd_re)
    #mpfr_init2(x1,prec)
    #mpfr_init2(y1,prec)
    ui = S._unitary_action
    weight = RF(S.weight())
    G = S.group()
    multiplier = S.multiplier()
    if weight<>0 or not multiplier.is_trivial():
        non_trivial=True
    else:
        non_trivial=False
    trivial_mult = multiplier.is_trivial()
    twopi=RF(2)*RF.pi() 
    twopii=CF(0,twopi)
    mp_i=CF(0,1)
    cdef RealNumber one
    one=RF(1)
    ep0=RF(1E-10)
    #Xm=vector(RF,Qf-Qs+1)
    VS = vector(RF,Qf-Qs+1).parent()
    Xm = Vector_real_mpfr_dense(VS,0)
    Xpb=dict() #vector(RF,2*Q)
    Ypb=dict()  #vector(RF,2*Q)
    Cvec=dict()
    #cdef mpc_t ***Cvec
    
    if S._holomorphic and not S._weak and weight==0:
        raise Exception,"Must support a non-zero weight for non-weak holomorphic forms!"
    verbose = S._verbose
    if verbose>3:
        fp0=open("xm.txt","w")
        fp1=open("xpb.txt","w")
        fp2=open("ypb.txt","w")
        fp3=open("cv.txt","w")
    #twoQl=RF(2*(Qf+1-Qs))
    if Qs<0:
        Qfak = 2*(Qf-Qs+1)
        for j in range(Qs,Qf+1):
            mpfr_set_si(Xm._entries[j-Qs],2*j-1,rnd_re)
            mpfr_div_si(Xm._entries[j-Qs],Xm._entries[j-Qs],Qfak,rnd_re)
    else:
        Qfak = 4*(Qf-Qs+1)
        for j in range(Qs,Qf+1):
            mpfr_set_si(Xm._entries[j-Qs],2*j-1,rnd_re)
            mpfr_div_si(Xm._entries[j-Qs],Xm._entries[j-Qs],Qfak,rnd_re)
    cdef int ci,cj
    cdef int cia,cib,cja,cjb,ciindex,cjindex
    cdef int nc=G._ncusps
    cdef int nv=G._nvertices
    #nc = len(G._cusps)
    #cdef RealNumber swi
    cdef int is_Gamma0,nreps
    is_Gamma0=<int>G._is_Gamma0
    cdef int*** reps=NULL
    cdef int N
    if is_Gamma0==1:
        nreps=G._index
        N=G._level
        reps= <int ***> sig_malloc(sizeof(int**) * nreps)
        for j in range(nreps):
            reps[j]=NULL
            reps[j]=<int **> sig_malloc(sizeof(int*) * 2)
            if reps[j]==NULL:
                raise MemoryError
            reps[j][0]=NULL;reps[j][1]=NULL
            reps[j][0]=<int *> sig_malloc(sizeof(int) * 2)
            if reps[j][0]==NULL:
                raise MemoryError
            reps[j][1]=<int *> sig_malloc(sizeof(int) * 2)
            if reps[j][1]==NULL:
                raise MemoryError
            reps[j][0][0]=G._coset_reps_v0[j][0]
            reps[j][0][1]=G._coset_reps_v0[j][1]
            reps[j][1][0]=G._coset_reps_v0[j][2]
            reps[j][1][1]=G._coset_reps_v0[j][3]
    
    #if verbose>2:
    #    print "nc=",nc
    if not isinstance(G._cusp_data[0]['width'],(int,Integer)):
        use_int=0
    cdef int** normalizers=NULL
    cdef int** cusp_maps=NULL
    normalizers=<int**>sig_malloc(sizeof(int*)*nc)
    if normalizers==NULL: raise MemoryError
    for i in range(nc):
        normalizers[i]=<int*>sig_malloc(sizeof(int)*4)
        a,b,c,d=G._cusp_data[i]['normalizer']
        normalizers[i][0]=a
        normalizers[i][1]=b
        normalizers[i][2]=c
        normalizers[i][3]=d
        if verbose>2:
            print "Normalizer[",i,"]=",a,b,c,d
    vertex_maps=<int**>sig_malloc(sizeof(int*)*nv)
    if vertex_maps==NULL: raise MemoryError
    cusp_maps=<int**>sig_malloc(sizeof(int*)*nv)
    if cusp_maps==NULL: raise MemoryError
    vertex_widths=<double*>sig_malloc(sizeof(double)*nv)
    if vertex_widths==NULL: raise MemoryError
    vertex_cusp=<int*>sig_malloc(sizeof(int)*nv)    
    if vertex_cusp==NULL: raise MemoryError
    if cusp_maps==NULL: raise MemoryError
    for i in range(nv):
        cusp_maps[i]=<int*>sig_malloc(sizeof(int)*4)
        cusp_maps[i][0]=G._cusp_maps[i][0]
        cusp_maps[i][1]=G._cusp_maps[i][1]
        cusp_maps[i][2]=G._cusp_maps[i][2]
        cusp_maps[i][3]=G._cusp_maps[i][3]
        vertex_maps[i]=<int*>sig_malloc(sizeof(int)*4)
        vertex_maps[i][0]=<int>G._vertex_maps[i].a()
        vertex_maps[i][1]=<int>G._vertex_maps[i].b()
        vertex_maps[i][2]=<int>G._vertex_maps[i].c()
        vertex_maps[i][3]=<int>G._vertex_maps[i].d()
        vertex_widths[i]=<double>G._vertex_widths[i]
        vertex_cusp[i]=<int>G._vertex_data[i]['cusp']
    cdef double* widths
    widths=<double*>sig_malloc(sizeof(double)*nc)
    for i in range(nc):
        if verbose>3:
            print "width[",i,"]=",G._cusp_data[i]['width'],type(G._cusp_data[i]['width'])
        widths[i]=<double>G._cusp_data[i]['width']
    cdef int wi,wj,vi,vj
    cdef SL2Z_elt A,Tj
    cdef double dp_x,dp_y
    for ci in range(nc):
        wi = <int>widths[ci] #G._cusp_data[ci]['width'] #(ci)
        if verbose>1:
            print "cusp =",ci
        if wi<>1:
            mpfr_set_ui(swi,wi,rnd_re)
            mpfr_sqrt(swi,swi,rnd_re)
        else:
            mpfr_set_ui(swi,1,rnd_re)
        for j in range(Qs,Qf+1):
            #if ciindex==0:
            mpfr_set(x0,Xm._entries[j-Qs],rnd_re)
            mpfr_set(y0,Y.value,rnd_re)
            if ci<>0:
                a=normalizers[ci][0]; b=normalizers[ci][1]
                c=normalizers[ci][2]; d=normalizers[ci][3]
                _normalize_point_to_cusp_mpfr(x0,y0,a,b,c,d,wi)
            pullback_to_Gamma0N_mpfr_c(G,x1.value,y1.value,x0,y0,&a,&b,&c,&d)
            Tj=SL2Z_elt(a,b,c,d)
            #dp_x=mpfr_get_d(x1,rnd_re)
            dp_x = float(x1); dp_y = float(y1)
            #dp_y=mpfr_get_d(y1,rnd_re)
            vj = closest_vertex_dp_c(nv,vertex_maps,vertex_widths,&dp_x,&dp_y)
            cj = vertex_cusp[vj]
            #cj,vj= G.closest_cusp(x1,y1,vertex=1)
            #cja,cjb= G.closest_cusp(x1,y1) #G._vertex_data[G._vertices[vi]]['cusp'] 
            #cjindex=G._cusps.index((cja,cjb))
            wj=<int>widths[cj] #G._cusp_data[cj]['width']
            if wj<>1:
                mpfr_set_ui(swj,wj,rnd_re)
                mpfr_sqrt(swj,swj,rnd_re)
            else:
                mpfr_set_ui(swj,1,rnd_re)
                #swj=sqrt(<double>wj)
            #wj = G._cusp_data[(cja,cjb)]['width']
            a = cusp_maps[vj][0]   #a,b,c,d=G._vertex_data[cj]['cusp_map'] #[v]
            b = cusp_maps[vj][1]
            c = cusp_maps[vj][2]
            d = cusp_maps[vj][3]  #U = G._vertex_data[vj]['cusp_map']
            x2=x1; y2=y1
            if a<>1 or b<>0 or c<>0 or d<>1: #U<>SL2Z_elt(1,0,0,1):
                _apply_sl2z_map_mpfr(x2.value,y2.value,a,b,c,d) #U[0,0],U[0,1],U[1,0],U[1,1])
            if cj>0:
                x3=x2; y3=y2
                a=normalizers[cj][3]; b=-normalizers[cj][1]
                c=-normalizers[cj][2]; d=normalizers[cj][0]
                if verbose>2:
                    print "Normalizer: a,b,c,d=",a,b,c,d
                _normalize_point_to_cusp_mpfr(x3.value,y3.value,a,b,c,d,wj,1)
                #[x3,y3] = normalize_point_to_cusp_mpfr(G,cja,cjb,x2,y2,inv=1)
            else:
                x3=x2; y3=y2
            Xpb[ci,cj,j]=x3*twopi
            # Recall that Ypb must be greater than Y otherwise we need a better Y
            if y3>=Y-ep0:
                Ypb[ci,cj,j]=y3*one ## Ridiculous but necessary due to bug in cython...
            else:
                if verbose > 0:
                    print "ci,cj=",ci,cj
                    print "Xm=",Xm[j-Qs]
                    #print "x,y=",x0,"\n",y0
                    print "x1,y1=",x1,"\n",y1
                    print "x2,y2=",x2,"\n",y2
                    print "x3,y3=",x3,"\n",y3
                    print "Xpb=",Xpb[ci,cj]
                    print "Ypb=",Ypb[ci,cj,j]
                raise ArithmeticError,"Need smaller value of Y" 
            if verbose > 2:
                print "ci,cj=",ci,cj
                print "Xm=",Xm[j-Qs]
                print "x,y=",mpfr_get_d(x0,rnd_re),"\n",mpfr_get_d(y0,rnd_re)
                print "x1,y1=",x1,"\n",y1
                print "x2,y2=",x2,"\n",y2
                print "x3,y3=",x3,"\n",y3
                print "Xpb=",Xpb[ci,cj,j]
                print "Ypb=",Ypb[ci,cj,j]
            # We also get the multiplier if we need it
            if(non_trivial):
                if weight<>0:
                    c = normalizers[ci][2]
                    d = normalizers[ci][3]
                    cr = RF(c); dr = RF(d)
                    mpfr_mul_si(cr.value,swi,c,rnd_re)
                    mpfr_set_si(dr.value,d,rnd_re)
                    mpfr_div(dr.value,dr.value,swi,rnd_re)
                    #cr = RF(G._cusp_data[ci]['normalizer'][1,0])*swi
                    #dr = RF(G._cusp_data[ci]['normalizer'][1,1])/swi
                    #print "test:",cr,dr,Xm[j],Y,-weight
                    #print type(c),type(d)
                    m1=j_fak_mpc(cr,dr,Xm[j-Qs],Y,-weight,ui)
                    c = normalizers[cj][2]
                    d = normalizers[cj][3]
                    cr = RF(c); dr = RF(d)
                    mpfr_mul_si(cr.value,swj,c,rnd_re)
                    mpfr_set_si(dr.value,d,rnd_re)
                    mpfr_div(dr.value,dr.value,swj,rnd_re)
                    #cr = RF(G._cusp_data[cj]['normalizer'][1,0])*swj
                    #dr = RF(G._cusp_data[cj]['normalizer'][1,1])/swj
                    #print "test:",cr,dr,x3,y3,-weight
                    m2=j_fak_mpc(cr,dr,x3,y3,weight,ui)
                    #aa = U[0,0]*a+U[0,1]*c
                    A = G._vertex_data[vj]['cusp_map']*Tj
                    #A=U.matrix()*Tj

                    #A=A**-1
                    cr=RF(-A.c()); dr=RF(A.a())
                    m3=j_fak_mpc(cr,dr,x2,y2,weight,ui)
                    tmp=m1*m2*m3
                    #print "m1,m2,m3=",m1,m2,m3
                else:
                    #print "U:",U.matrix(),type(U.matrix())
                    #print "Tj:",Tj,type(Tj)
                    #A=U.matrix()*Tj
                    if not trivial_mult :
                        A = G._vertex_data[vj]['cusp_map']*Tj
                    #print "A:",A,type(A)
                    tmp=CF(1)
                if not trivial_mult :
                    #if multiplier<>None and not multiplier.is_trivial():
                    AA=A**-1
                    #v = S.character(AA[1,1]).complex_embedding(CF.prec())
                    v = CF(multiplier(AA).complex_embedding(CF.prec()))
                    #v=CF(v)
                    tmp = tmp*v
                Cvec[ci,cj,j]=tmp # mp_ctx.mpc(0,tmp)
                #print "C[",cii,cjj,j,"]=",tmp
            else:
                Cvec[ci,cj,j]=CF(1) #CF(0,tmp)
            #if verbose>3:
            #    fp1.write(str(x3)+"\n")
            #    fp2.write(str(y3)+"\n")
            #    fp3.write(str(tmp)+"\n")

        #        Xmi[j]=mp_ctx.mpc(0,Xm[j]*twopi)
    for j in range(Qs,Qf+1): 
        mpfr_mul(Xm._entries[j-Qs],Xm._entries[j-Qs],twopi.value,rnd_re)
    if verbose>3:
        fp1.close()
        fp2.close()
        fp3.close()
    pb=dict()
    if normalizers<>NULL:
        for i in range(nc):
            if  normalizers[i]<>NULL:
                sig_free(normalizers[i])
    sig_free(normalizers)            
    pb['xm']=Xm; pb['xpb']=Xpb; pb['ypb']=Ypb; pb['cvec']=Cvec
    mpfr_clear(swi); mpfr_clear(swj)
    mpfr_clear(x0);  mpfr_clear(y0)
    mpfr_clear(YY)
    return pb

cpdef pullback_pts_mp(S,int Qs,int Qf,RealNumber Y,int holo=0):
    r"""
    Pullback the horocycle at height Y.
    If Qs=1 and Qf =Q then we use symmetry, otherwise not.
    """
    cdef mpfr_t* Xm_t
    cdef mpfr_t*** Xpb_t=NULL
    cdef mpfr_t*** Ypb_t=NULL
    cdef mpc_t*** Cvec_t=NULL
    cdef mpfr_t**** RCvec_t=NULL
    cdef int*** CSvec_t=NULL
    cdef int Ql,i,j,k,n,nc,prec
    cdef RealField_class RF
    cdef RealNumber weight,tmp
#    cdef RealNumber rtmp
    cdef list l
    nc= S.group().ncusps()
    Ql=Qf-Qs+1
    prec=Y.prec()
    if S._verbose>1:
        print "Yin=",Y
        print "prec=",prec
    RF = RealField(prec)
    CF = MPComplexField(prec)
#    rtmp=RF(0)
    weight = RF(S._weight)
    G=S._group
 
    tmp=RF(0)
    Xm_t = <mpfr_t*> sig_malloc( sizeof(mpfr_t) * Ql )
    for n in range(Ql):
        mpfr_init2(Xm_t[n],prec)
    Xpb_t = <mpfr_t***> sig_malloc( sizeof(mpfr_t** ) * nc )
    if Xpb_t==NULL: raise MemoryError
    Ypb_t = <mpfr_t***> sig_malloc( sizeof(mpfr_t** ) * nc )
    if Ypb_t==NULL: raise MemoryError
    for i in range(nc):
        Xpb_t[i] = <mpfr_t**>sig_malloc(sizeof(mpfr_t*) * nc )
        Ypb_t[i] = <mpfr_t**>sig_malloc(sizeof(mpfr_t*) * nc )
        if Ypb_t[i]==NULL or Xpb_t[i]==NULL:
            raise MemoryError
        for j in range(nc):
            Xpb_t[i][j] = <mpfr_t*>sig_malloc(sizeof(mpfr_t) * Ql )
            Ypb_t[i][j] = <mpfr_t*>sig_malloc(sizeof(mpfr_t) * Ql )
            if Ypb_t[i][j]==NULL or Xpb_t[i][j]==NULL:
                raise MemoryError
            for n in range(Ql):
                mpfr_init2(Xpb_t[i][j][n],prec) 
                mpfr_init2(Ypb_t[i][j][n],prec) 
                mpfr_set_si(Xpb_t[i][j][n],0,rnd_re)
                mpfr_set_si(Ypb_t[i][j][n],0,rnd_re)
                #Ypb_t[i][j][n]=<double>0
    if Qs<0 or is_Hecke_triangle_group(G):
        Cvec_t = <mpc_t***>sig_malloc(sizeof(mpc_t**) * nc )
        if Cvec_t==NULL: raise MemoryError
        for i in range(nc):
            Cvec_t[i] = <mpc_t**>sig_malloc(sizeof(mpc_t*) * nc )
            if Cvec_t[i]==NULL:
                raise MemoryError
            for j  in range(nc):
                Cvec_t[i][j] = <mpc_t*>sig_malloc(sizeof(mpc_t) * Ql )
                if Cvec_t[i][j]==NULL:
                    raise MemoryError
                for n in range(Ql):
                    mpc_init2(Cvec_t[i][j][n],prec)
                    mpc_set_si(Cvec_t[i][j][n],0,rnd)
    else:
        print "HERE0"        
        CSvec_t = <int***>sig_malloc(sizeof(int**) * nc )
        if CSvec_t==NULL: raise MemoryError
        for i in range(nc):
            CSvec_t[i] = <int**>sig_malloc(sizeof(int*) * nc )
            if CSvec_t[i]==NULL:
                raise MemoryError
            for j  in range(nc):
                CSvec_t[i][j] = <int*>sig_malloc(sizeof(int) * Ql )
                if CSvec_t[i][j]==NULL:
                    raise MemoryError
        print "HERE"
        RCvec_t = <mpfr_t****>sig_malloc(sizeof( Xpb_t) * nc )
        if RCvec_t==NULL: raise MemoryError
        for i in range(nc):
            RCvec_t[i] = <mpfr_t***>sig_malloc(sizeof(mpfr_t**) * nc )
            if RCvec_t[i]==NULL:
                raise MemoryError
            for j  in range(nc):
                RCvec_t[i][j] = <mpfr_t**>sig_malloc(sizeof(mpfr_t*) * Ql )
                if RCvec_t[i][j]==NULL:
                    raise MemoryError
                for n in range(Ql):
                    RCvec_t[i][j][n] = <mpfr_t*>sig_malloc(sizeof(mpfr_t) * 3 )
                    mpfr_init2(RCvec_t[i][j][n][0],prec)
                    mpfr_init2(RCvec_t[i][j][n][1],prec)
                    mpfr_init2(RCvec_t[i][j][n][2],prec)
                    
                    mpfr_set_si(RCvec_t[i][j][n][0],1,rnd_re)
                    mpfr_set_si(RCvec_t[i][j][n][1],0,rnd_re)
                    mpfr_set_si(RCvec_t[i][j][n][2],0,rnd_re)
                #Cvec_t[i][j][n]=<double complex>0
    print "HERE2"                
    cdef int is_hecke = 0
    cdef int res = 0
    if is_Hecke_triangle_group(G):
        print "HERE3"                
    
        res = pullback_pts_hecke_triangle_mpc_new_c(S,Qs,Qf,Y,Xm_t,Xpb_t,Ypb_t,Cvec_t)
        print "HERE4"                
    
    else:
        if Qs<0:
            res = pullback_pts_mpc_new_c(S,Qs,Qf,Y.value,Xm_t,Xpb_t,Ypb_t,Cvec_t)
        else:
            res = pullback_pts_mpc_new_c_sym(S,Qs,Qf,Y.value,Xm_t,Xpb_t,Ypb_t,RCvec_t,CSvec_t)
    if res==1:
        raise ArithmeticError,"Need smaller Y!"
    cdef dict Xm,Xpb,Ypb,Cvec,pb,RCvec,CSvec
    Xm=dict(); Xpb=dict() 
    Ypb=dict();Cvec=dict()
    RCvec=dict(); CSvec={}
    cdef MPComplexNumber ctmp
    cdef RealNumber rtmp,rtmp1,rtmp2,one
    one=RF(1)
    rtmp = RF(0);rtmp1 = RF(1); rtmp2 = RF(0)
    ctmp = CF(0)
    for n in range(Ql):
        mpfr_set(tmp.value,Xm_t[n],rnd_re)
        Xm[n]=tmp
    for i in range(nc):
        Xpb[i]=dict();Ypb[i]=dict()
        for j in range(nc):
            Xpb[i][j]=dict();Ypb[i][j]=dict()
            for n in range(Ql):
                mpfr_set(rtmp.value,Xpb_t[i][j][n],rnd_re)
                Xpb[i][j][n]=rtmp*one
                mpfr_set(rtmp.value,Ypb_t[i][j][n],rnd_re)
                Ypb[i][j][n]=rtmp*one
    print "here 5"
    pb=dict()
    pb['xm']=Xm; pb['xpb']=Xpb; pb['ypb']=Ypb 
    if Qs<0:
        for i in range(nc):
            Cvec[i]=dict()
            for j in range(nc):
                Cvec[i][j]=dict()
                for n in range(Ql):
                    mpc_set(ctmp.value,Cvec_t[i][j][n],rnd)
                    ## multiplying by one is a simple way to assuer a new element is created
                    ## thereby avoiding an issue where all elements point to the same value...
                    Cvec[i][j][n]=CF(ctmp)*one 
                    #print "RCvec0=",ctmp,"*",one

        pb['cvec']=Cvec
    else:
        for i in range(nc):
            RCvec[i]=dict()
            CSvec[i]=dict()
            for j in range(nc):
                RCvec[i][j]=dict()
                CSvec[i][j]=dict()
                for n in range(Ql):
                    #RCvec[i][j][n]=[]
                    #print "RCvec0=",RCvec[i][j][n]
                    rtmp1=RF(1);rtmp2=RF(0)
                    mpfr_set(rtmp1.value,RCvec_t[i][j][n][0],rnd_re)
                    mpfr_set(rtmp2.value,RCvec_t[i][j][n][1],rnd_re)
                    #CSvec[i][j][n]=CSvec_t[i][j][n]
                    #k = mpfr_get_si(RCvec_t[i][j][n][2],rnd_re)
                    #if S._verbose>1:
                    #    print "abs(cvec)=",rtmp1
                    #    print "arg(cvec)=",rtmp2                    
                    #l = [rtmp1,rtmp2,k]
                    #if S._verbose>1:
                    #    print "Rcvec[{0}][{1}][{2}]={3}".format(i,j,n,l)
                    RCvec[i][j][n]=[rtmp1,rtmp2]
                    #if S._verbose>1:
                    #    print "Rcvec[{0}][{1}][{2}]={3}".format(i,j,n,RCvec[i][j][n])
        pb['cvec']=RCvec
    #if S._verbose>1:
    #    print "Rcvec[{0}][{1}][{2}]={3}".format(1,1,8,RCvec[1][1][8])
    if Ypb_t<>NULL:
        for i in range(nc):
            if Ypb_t[i]<>NULL:
                for j in range(nc):
                    if Ypb_t[i][j]<>NULL:
                        sig_free(Ypb_t[i][j])
                sig_free(Ypb_t[i])
        sig_free(Ypb_t)
    if Xpb_t<>NULL:
        for i in range(nc):
            if Xpb_t[i]<>NULL:
                for j in range(nc):
                    if Xpb_t[i][j]<>NULL:
                        sig_free(Xpb_t[i][j])
                sig_free(Xpb_t[i])
        sig_free(Xpb_t)
    if Cvec_t<>NULL:
        for i in range(nc):
            if Cvec_t[i]<>NULL:
                for j in range(nc):
                    if Cvec_t[i][j]<>NULL:
                        sig_free(Cvec_t[i][j])
                sig_free(Cvec_t[i])
        sig_free(Cvec_t)
    if RCvec_t<>NULL:
        for i in range(nc):
            if RCvec_t[i]<>NULL:
                for j in range(nc):
                    if RCvec_t[i][j]<>NULL:
                        for n in range(Ql):
                            if RCvec_t[i][j][n]<>NULL:
                                mpfr_clear(RCvec_t[i][j][n][0])
                                mpfr_clear(RCvec_t[i][j][n][1])
                                mpfr_clear(RCvec_t[i][j][n][2])
                                sig_free(RCvec_t[i][j][n])
                        sig_free(RCvec_t[i][j])
                sig_free(RCvec_t[i])
        sig_free(RCvec_t)
    if CSvec_t<>NULL:
        for i in range(nc):
            if CSvec_t[i]<>NULL:
                for j in range(nc):
                    if CSvec_t[i][j]<>NULL:
                        sig_free(CSvec_t[i][j])
                sig_free(CSvec_t[i])
        sig_free(CSvec_t)
    return pb

@cython.cdivision(True)
cdef int pullback_pts_mpc_new_c(S,int Qs,int Qf,mpfr_t Y, mpfr_t* Xm,mpfr_t*** Xpb, mpfr_t*** Ypb,mpc_t ***Cvec):

    r""" Computes a whole array of pullbacked points-
         using MPFR/MPC types.
    INPUT:

    - ``S``  -- AutomorphicFormSpace
    - ``Qs`` -- integer
    - ``Qf`` -- integer
    - ``Y``  -- real >0
    - ``deb``    -- (False) logical

    
    OUTPUT:

    - ``pb`` -- dictionary with entries:
       - 'xm'   -- real[Qf-Qs+1]  : x_m=2*pi*(1-2*m)/2(Qf-Qs+1)
       - 'xpb'  -- real[0:nc,0:nc,Qf-Qs+1]  : real part of pullback of x_m+iY
       - 'ypb'  -- real[0:nc,0:nc,Qf-Qs+1]  : imag. part of pullback of x_m+iY
       - 'cvec' -- complex[0:nc,0:nc,Qf-Qs+1] : mult. factor

    EXAMPLES::


        sage: G=MySubgroup(Gamma0(2))
        sage: [Xm,Xpb,Ypb,Cv]=pullback_pts(G,1,3,0.1)   

    Note that we need to have $Y<Y_{0}$ so that all pullbacked points
    are higher up than Y.

        sage: pullback_pts(G,1,10,0.44)
        Traceback (most recent call last):
        ...
        ArithmeticError: Need smaller value of Y

     Returns 1 if Y is too large
    
    """
    # WE let the input determine the precision
    #print "pullback mpc!"
    cdef int prec=mpfr_get_prec(Y)
    CF=MPComplexField(prec)
    RF=RealField(prec)
    cdef RealNumber x1,y1,x2,y2,x3,y3,one,weight
    x1=RF(1)
    x2=RF(1)
    x3=RF(1)
    y1=RF(1)
    y2=RF(1)
    y3=RF(1)
    if S._verbose>1:
        print "Y in pb =",mpfr_get_d(Y,rnd_re)
        print "prec=",prec
    cdef int a,b,c,d,v0,v1,itmp,i,j,ui,Ql
    cdef mpfr_t swj,swi,x0,y0,YY
    cdef RealNumber ar,br,cr,dr,twopi,xm,ep0,tmp,tmp1
    cdef MPComplexNumber ctmp,m1,m2,m3
    m1=CF(1); m2=CF(1); m3=CF(1)
    ar=RF(1); br=RF(0); cr=RF(0); dr=RF(1)
    tmp=RF(0); tmp1=RF(0)
    mpfr_init2(swi,prec)
    mpfr_init2(swj,prec)
    mpfr_init2(x0,prec)
    mpfr_init2(y0,prec)
    mpfr_init2(YY,prec)
    mpfr_set(YY,Y,rnd_re)
    #mpfr_init2(x1,prec)
    #mpfr_init2(y1,prec)
    Ql = Qf - Qs + 1
    ui = S._unitary_action
    weight = RF(S.weight())
    G = S.group()
    multiplier = S.multiplier()
    if weight<>0 or not multiplier.is_trivial():
        non_trivial=True
    else:
        non_trivial=False
    trivial_mult = multiplier.is_trivial()
    twopi=RF(2)*RF.pi() 
    twopii=CF(0,twopi)
    mp_i=CF(0,1)
    ctmp= CF(0)
    one=RF(1)
    ep0=RF(1E-10)
    #Xm=vector(RF,Qf-Qs+1)
    #VS = vector(RF,Qf-Qs+1).parent()
    #Xm = Vector_real_mpfr_dense(VS,0)
    #cdef mpc_t ***Cvec
    
    if S._holomorphic and not S._weak and weight==0:
        raise Exception,"Must support a non-zero weight for non-weak holomorphic forms!"
    verbose = S._verbose
    if Qs<0:
        Qfak = 2*(Qf-Qs+1)
        for j in range(Qs,Qf+1):
            mpfr_set_si(Xm[j-Qs],2*j-1,rnd_re)
            mpfr_div_si(Xm[j-Qs],Xm[j-Qs],Qfak,rnd_re)
    else:
        Qfak = 4*(Qf-Qs+1)
        for j in range(Qs,Qf+1):
            mpfr_set_si(Xm[j-Qs],2*j-1,rnd_re)
            mpfr_div_si(Xm[j-Qs],Xm[j-Qs],Qfak,rnd_re)
    cdef int ci,cj
    cdef int cia,cib,cja,cjb,ciindex,cjindex
    cdef int nc=G._ncusps
    cdef int nv=G._nvertices
    #nc = len(G._cusps)
    #cdef RealNumber swi
    cdef int is_Gamma0,nreps
    is_Gamma0=<int>G._is_Gamma0
    cdef int*** reps=NULL
    cdef int N
    if is_Gamma0==1:
        nreps=G._index
        N=G._level
        reps= <int ***> sig_malloc(sizeof(int**) * nreps)
        for j in range(nreps):
            reps[j]=NULL
            reps[j]=<int **> sig_malloc(sizeof(int*) * 2)
            if reps[j]==NULL:
                raise MemoryError
            reps[j][0]=NULL;reps[j][1]=NULL
            reps[j][0]=<int *> sig_malloc(sizeof(int) * 2)
            if reps[j][0]==NULL:
                raise MemoryError
            reps[j][1]=<int *> sig_malloc(sizeof(int) * 2)
            if reps[j][1]==NULL:
                raise MemoryError
            reps[j][0][0]=G._coset_reps_v0[j][0]
            reps[j][0][1]=G._coset_reps_v0[j][1]
            reps[j][1][0]=G._coset_reps_v0[j][2]
            reps[j][1][1]=G._coset_reps_v0[j][3]
    
    #if verbose>2:
    #    print "nc=",nc
    if not isinstance(G._cusp_data[0]['width'],(int,Integer)):
        use_int=0
    cdef int** normalizers=NULL
    cdef int** cusp_maps=NULL
    normalizers=<int**>sig_malloc(sizeof(int*)*nc)
    if normalizers==NULL: raise MemoryError
    for i in range(nc):
        normalizers[i]=<int*>sig_malloc(sizeof(int)*4)
        a,b,c,d=G._cusp_data[i]['normalizer']
        normalizers[i][0]=a
        normalizers[i][1]=b
        normalizers[i][2]=c
        normalizers[i][3]=d
        if verbose>2:
            print "Normalizer[",i,"]=",a,b,c,d
    vertex_maps=<int**>sig_malloc(sizeof(int*)*nv)
    if vertex_maps==NULL: raise MemoryError
    cusp_maps=<int**>sig_malloc(sizeof(int*)*nv)
    if cusp_maps==NULL: raise MemoryError
    vertex_widths=<double*>sig_malloc(sizeof(double)*nv)
    if vertex_widths==NULL: raise MemoryError
    vertex_cusp=<int*>sig_malloc(sizeof(int)*nv)    
    if vertex_cusp==NULL: raise MemoryError
    if cusp_maps==NULL: raise MemoryError
    for i in range(nv):
        cusp_maps[i]=<int*>sig_malloc(sizeof(int)*4)
        cusp_maps[i][0]=G._cusp_maps[i][0]
        cusp_maps[i][1]=G._cusp_maps[i][1]
        cusp_maps[i][2]=G._cusp_maps[i][2]
        cusp_maps[i][3]=G._cusp_maps[i][3]
        vertex_maps[i]=<int*>sig_malloc(sizeof(int)*4)
        vertex_maps[i][0]=<int>G._vertex_maps[i].a()
        vertex_maps[i][1]=<int>G._vertex_maps[i].b()
        vertex_maps[i][2]=<int>G._vertex_maps[i].c()
        vertex_maps[i][3]=<int>G._vertex_maps[i].d()
        vertex_widths[i]=<double>G._vertex_widths[i]
        vertex_cusp[i]=<int>G._vertex_data[i]['cusp']
    cdef double* widths
    widths=<double*>sig_malloc(sizeof(double)*nc)
    for i in range(nc):
        if verbose>3:
            print "width[",i,"]=",G._cusp_data[i]['width'],type(G._cusp_data[i]['width'])
        widths[i]=<double>G._cusp_data[i]['width']
    cdef int wi,wj,vi,vj
    cdef SL2Z_elt A,Tj
    cdef double dp_x,dp_y
    for ci in range(nc):
        wi = <int>widths[ci] #G._cusp_data[ci]['width'] #(ci)
        if verbose>1:
            print "cusp =",ci
        if wi<>1:
            mpfr_set_ui(swi,wi,rnd_re)
            mpfr_sqrt(swi,swi,rnd_re)
        else:
            mpfr_set_ui(swi,1,rnd_re)
        for j in range(Ql): #Qs,Qf+1): #from Qs<= j <= Qf:
            #if ciindex==0:
            mpfr_set(x0,Xm[j],rnd_re)
            mpfr_set(y0,Y,rnd_re)
            if ci<>0:
                a=normalizers[ci][0]; b=normalizers[ci][1]
                c=normalizers[ci][2]; d=normalizers[ci][3]
                _normalize_point_to_cusp_mpfr(x0,y0,a,b,c,d,wi)
            pullback_to_Gamma0N_mpfr_c(G,x1.value,y1.value,x0,y0,&a,&b,&c,&d)
            Tj=SL2Z_elt(a,b,c,d)
            #if verbose > 2:
            #    print "pull back of {0},{1} is {2},{3}".format(Xm[j],mpfr_get_d(Y,rnd_re),x1,y1)
            #    print "map:",a,b,c,d
            #dp_x=mpfr_get_Ad(x1,rnd_re)
            dp_x = float(x1); dp_y = float(y1)
            #dp_y=mpfr_get_d(y1,rnd_re)
            vj = closest_vertex_dp_c(nv,vertex_maps,vertex_widths,&dp_x,&dp_y)
            cj = vertex_cusp[vj]
            #cj,vj= G.closest_cusp(x1,y1,vertex=1)
            #cja,cjb= G.closest_cusp(x1,y1) #G._vertex_data[G._vertices[vi]]['cusp'] 
            #cjindex=G._cusps.index((cja,cjb))
            wj=<int>widths[cj] #G._cusp_data[cj]['width']
            if wj<>1:
                mpfr_set_ui(swj,wj,rnd_re)
                mpfr_sqrt(swj,swj,rnd_re)
            else:
                mpfr_set_ui(swj,1,rnd_re)
                #swj=sqrt(<double>wj)
            #wj = G._cusp_data[(cja,cjb)]['width']
            a = cusp_maps[vj][0]   #a,b,c,d=G._vertex_data[cj]['cusp_map'] #[v]
            b = cusp_maps[vj][1]
            c = cusp_maps[vj][2]
            d = cusp_maps[vj][3]  #U = G._vertex_data[vj]['cusp_map']
            mpfr_set(x2.value,x1.value,rnd_re)
            mpfr_set(y2.value,y1.value,rnd_re)
            #x2=x1; y2=y1
            if a<>1 or b<>0 or c<>0 or d<>1: #U<>SL2Z_elt(1,0,0,1):
                _apply_sl2z_map_mpfr(x2.value,y2.value,a,b,c,d) #U[0,0],U[0,1],U[1,0],U[1,1])
            if cj>0:
                mpfr_set(x3.value,x2.value,rnd_re)
                mpfr_set(y3.value,y2.value,rnd_re)
                #x3=x2; y3=y2
                a=normalizers[cj][3]; b=-normalizers[cj][1]
                c=-normalizers[cj][2]; d=normalizers[cj][0]
                if verbose>2:
                    print "Normalizer: a,b,c,d=",a,b,c,d
                _normalize_point_to_cusp_mpfr(x3.value,y3.value,a,b,c,d,wj,1)
            else:
                mpfr_set(x3.value,x2.value,rnd_re)
                mpfr_set(y3.value,y2.value,rnd_re)
                #x3=x2; y3=y2
            #x3 = x3*twopi
            mpfr_set(Xpb[ci][cj][j],x3.value,rnd_re)
            mpfr_mul(Xpb[ci][cj][j],Xpb[ci][cj][j],twopi.value,rnd_re)
            # Recall that Ypb must be greater than Y otherwise we need a better Y
            tmp=RF(y3-ep0)
            if mpfr_cmp(tmp.value,Y)>=0:
                mpfr_set(Ypb[ci][cj][j],y3.value,rnd_re)
            else:
                if verbose > 0:
                    print "ci,cj=",ci,cj
                    print "Xm=",mpfr_get_d(Xm[j],rnd_re)
                    #print "x,y=",x0,"\n",y0
                    print "x1,y1=",x1,"\n",y1
                    print "x2,y2=",x2,"\n",y2
                    print "x3,y3=",x3,"\n",y3
                    print "Xpb=",mpfr_get_d(Xpb[ci][cj][j],rnd_re)
                    print "Ypb=",mpfr_get_d(Ypb[ci][cj][j],rnd_re)
                return 1 #raise ArithmeticError,"Need smaller value of Y" 
            if verbose > 2:
                print "ci,cj=",ci,cj
                print "Xm=",mpfr_get_d(Xm[j],rnd_re)
                print "x,y=",mpfr_get_d(x0,rnd_re),"\n",mpfr_get_d(y0,rnd_re)
                print "x1,y1=",x1,"\n",y1
                print "x2,y2=",x2,"\n",y2
                print "x3,y3=",x3,"\n",y3
                print "Xpb=",mpfr_get_d(Xpb[ci][cj][j],rnd_re)
                print "Ypb=",mpfr_get_d(Ypb[ci][cj][j],rnd_re)
                print "non_trivial=",non_trivial
            # We also get the multiplier if we need it
            if(non_trivial):
                if weight<>0:
                    c = normalizers[ci][2]
                    d = normalizers[ci][3]
                    #if c==0 and d==1:
                    #    #m1=CF(1)
                    #    mpc_set_ui(m1.value,1,rnd)
                    #else:
                    # cr = RF(c) #; dr = RF(d)
                    mpfr_mul_si(cr.value,swi,c,rnd_re)
                    mpfr_set_si(dr.value,d,rnd_re)
                    mpfr_div(dr.value,dr.value,swi,rnd_re)
                    mpfr_set(tmp.value,Xm[j],rnd_re)
                    mpfr_set(tmp1.value,Y,rnd_re)
                    m1=j_fak_mpc(cr,dr,tmp,tmp1,-weight,ui)
                    #m1=j_fak_old(cr,dr,Xm[j],Y,-weight,ui)
                    c = normalizers[cj][2]
                    d = normalizers[cj][3]
                    if verbose>2:
                        print "c=",c
                        print "d=",d
                        print "weight=",weight
                        #print "Xm[",j,"]=",x3 #Xm[j]
                        #print "Y=",y3
                    #if c==0 and d==1:
                    #    mpc_set_ui(m2.value,1,rnd)
                    #    #m2 = CF(1)
                    #else:
                        #cr = RF(c); dr = RF(d)
                    mpfr_mul_si(cr.value,swj,c,rnd_re)
                    mpfr_set_si(dr.value,d,rnd_re)
                    mpfr_div(dr.value,dr.value,swj,rnd_re)
                    if verbose>2:
                        print "cr=",cr
                        print "dr=",dr
                        print "x3,y3=",x3,y3
                        print "weight=",weight
                        print "ui=",ui
                    m2=j_fak_mpc(cr,dr,x3,y3,weight,ui)
                    if verbose>2:
                        print "m2=",m2
                    #m2=j_fak_old(cr,dr,x3,y3,weight,ui)
                    # aa = U[0,0]*a+U[0,1]*c
                    A = G._vertex_data[vj]['cusp_map']*Tj
                    #A=U.matrix()*Tj
                    #A=A**-1
                    c = -A.c(); d= A.a()
                    #cr=RF(-A.c()); dr=RF(A.a())
                    #if c==0 and d==1:
                    #    mpc_set_si(m3.value,1,rnd)
                    #    #m3 = CF(1)
                    #else:
                    cr=RF(c); dr=RF(d)
                    m3=j_fak_mpc(cr,dr,x2,y2,weight,ui)
                    #m3=j_fak_old(cr,dr,x2,y2,weight,ui)
                    ctmp=m1*m2*m3
                    if verbose>2:
                        print "m1=",m1
                        print "m2=",m2
                        print "m3=",m3
                else:
                    #A=U.matrix()*Tj
                    if not trivial_mult :
                        A = G._vertex_data[vj]['cusp_map']*Tj
                    #print "A:",A,type(A)
                    ctmp=CF(1)
                if not trivial_mult :
                    #if multiplier<>None and not multiplier.is_trivial():
                    AA=A**-1
                    #v = S.character(AA[1,1]).complex_embedding(CF.prec())
                    mA = multiplier(AA)
                    if hasattr(mA,"complex_embedding"):
                        v = CF(mA.complex_embedding(CF.prec()))
                    else:
                        v = CF(mA)
                    #v=CF(v)
                    #if verbose>2:
                    #print "v=",v
                    #print "ctmp=",ctmp
                    ctmp = ctmp*v
                mpc_set(Cvec[ci][cj][j],ctmp.value,rnd)
            else:
                mpc_set_ui(Cvec[ci][cj][j],1,rnd)
        #        Xmi[j]=mp_ctx.mpc(0,Xm[j]*twopi)
    for j in range(Ql):
        mpfr_mul(Xm[j],Xm[j],twopi.value,rnd_re)
    if normalizers<>NULL:
        for i in range(nc):
            if  normalizers[i]<>NULL:
                sig_free(normalizers[i])
        sig_free(normalizers)            
    mpfr_clear(swi); mpfr_clear(swj)
    mpfr_clear(x0);  mpfr_clear(y0)
    mpfr_clear(YY)
    return 0

@cython.cdivision(True)
cdef int pullback_pts_mpc_new_c_sym(S,int Qs,int Qf,mpfr_t Y, mpfr_t* Xm,mpfr_t*** Xpb, mpfr_t*** Ypb,mpfr_t ****RCvec,int*** CSvec):

    r""" Computes a whole array of pullbacked points-
         using MPFR/MPC types.
    INPUT:

    - ``S``  -- AutomorphicFormSpace
    - ``Qs`` -- integer
    - ``Qf`` -- integer
    - ``Y``  -- real >0
    - ``deb``    -- (False) logical

    
    OUTPUT:

    - ``pb`` -- dictionary with entries:
       - 'xm'   -- real[Qf-Qs+1]  : x_m=2*pi*(1-2*m)/2(Qf-Qs+1)
       - 'xpb'  -- real[0:nc,0:nc,Qf-Qs+1]  : real part of pullback of x_m+iY
       - 'ypb'  -- real[0:nc,0:nc,Qf-Qs+1]  : imag. part of pullback of x_m+iY
       - 'cvec' -- Real[0:nc,0:nc,Qf-Qs+1,0:2] : arg, abs and
       - 'sign' -- integer [0:nc,0:nc,Qf-Qs+1] : determines if we use cos or sine

    EXAMPLES::


        sage: G=MySubgroup(Gamma0(2))
        sage: [Xm,Xpb,Ypb,Cv]=pullback_pts(G,1,3,0.1)   

    Note that we need to have $Y<Y_{0}$ so that all pullbacked points
    are higher up than Y.

        sage: pullback_pts(G,1,10,0.44)
        Traceback (most recent call last):
        ...
        ArithmeticError: Need smaller value of Y

    
    """
    # WE let the input determine the precision
    #print "pullback mpc!"
    cdef int prec=mpfr_get_prec(Y)
    CF=MPComplexField(prec)
    RF=RealField(prec)
    cdef RealNumber x1,y1,x2,y2,x3,y3,one,weight,rtmp
    x1=RF(1)
    x2=RF(1)
    x3=RF(1)
    y1=RF(1)
    y2=RF(1)
    y3=RF(1)
    if S._verbose>1:
        print "Y in pb =",mpfr_get_d(Y,rnd_re)
        print "prec=",prec
    cdef int a,b,c,d,v0,v1,itmp,i,j,ui,Ql,k,ch_m1
    cdef mpfr_t swj,swi,x0,y0,YY,tmpab1,tmpar1,tmpab2,tmpar2,tmpab3,tmpar3,tmpar,tmpab,weight_t
    cdef RealNumber ar,br,cr,dr,twopi,xm,ep0,tmp
    cdef MPComplexNumber ctmp,m1,m2,m3
    
    m1=CF(1); m2=CF(1); m3=CF(1)
    ar=RF(1); br=RF(0); cr=RF(0); dr=RF(1)
    rtmp = RF(0); tmp=RF(0)
    mpfr_init2(swi,prec)
    mpfr_init2(swj,prec)
    mpfr_init2(x0,prec)
    mpfr_init2(y0,prec)
    mpfr_init2(YY,prec)
    mpfr_set(YY,Y,rnd_re)
    mpfr_init2(tmpab,prec); mpfr_init2(tmpar,prec)
    mpfr_init2(tmpab1,prec); mpfr_init2(tmpar1,prec)
    mpfr_init2(tmpab2,prec); mpfr_init2(tmpar2,prec)
    mpfr_init2(tmpab3,prec); mpfr_init2(tmpar3,prec)
    mpfr_init2(weight_t,prec)
    #mpfr_init2(x1,prec)
    #mpfr_init2(y1,prec)
    Ql = Qf - Qs + 1
    ui = S._unitary_action
    weight = RF(S.weight())
    k = int(weight)
    if abs(k-weight)>1e-10:
        raise ValueError,"Only use this routine for integer weight! (for now)"
    mpfr_set(weight_t,weight.value,rnd_re)
    G = S.group()
    multiplier = S.multiplier()
    if weight<>0 or not multiplier.is_trivial():
        non_trivial=True
    else:
        non_trivial=False
    if multiplier.character()(-1)==1:
        ch_m1=0
    else:
        ch_m1=1
    trivial_mult = multiplier.is_trivial()
    cdef mpfr_t mppi
    mpfr_init2(mppi,prec)
    twopi=RF.pi()
    mpfr_set(mppi,twopi.value,rnd_re)
    twopi=RF(2)*twopi
    twopii=CF(0,twopi)
    mp_i=CF(0,1)
    ctmp= CF(0)
    one=RF(1)
    ep0=RF(1E-10)
    #Xm=vector(RF,Qf-Qs+1)
    #VS = vector(RF,Qf-Qs+1).parent()
    #Xm = Vector_real_mpfr_dense(VS,0)
    #cdef mpc_t ***Cvec
    
    if S._holomorphic and not S._weak and weight==0:
        raise Exception,"Must support a non-zero weight for non-weak holomorphic forms!"
    verbose = S._verbose
    if Qs<0:
        Qfak = 2*(Qf-Qs+1)
        for j in range(Qs,Qf+1):
            mpfr_set_si(Xm[j-Qs],2*j-1,rnd_re)
            mpfr_div_si(Xm[j-Qs],Xm[j-Qs],Qfak,rnd_re)
    else:
        Qfak = 4*(Qf-Qs+1)
        for j in range(Qs,Qf+1):
            mpfr_set_si(Xm[j-Qs],2*(j-Qf)-1,rnd_re)
            mpfr_div_si(Xm[j-Qs],Xm[j-Qs],Qfak,rnd_re)
    cdef int ci,cj
    cdef int cia,cib,cja,cjb,ciindex,cjindex
    cdef int nc=G._ncusps
    cdef int nv=G._nvertices
    #nc = len(G._cusps)
    #cdef RealNumber swi
    cdef int is_Gamma0,nreps
    is_Gamma0=<int>G._is_Gamma0
    cdef int*** reps=NULL
    cdef int N
    if is_Gamma0==1:
        nreps=G._index
        N=G._level
        reps= <int ***> sig_malloc(sizeof(int**) * nreps)
        for j in range(nreps):
            reps[j]=NULL
            reps[j]=<int **> sig_malloc(sizeof(int*) * 2)
            if reps[j]==NULL:
                raise MemoryError
            reps[j][0]=NULL;reps[j][1]=NULL
            reps[j][0]=<int *> sig_malloc(sizeof(int) * 2)
            if reps[j][0]==NULL:
                raise MemoryError
            reps[j][1]=<int *> sig_malloc(sizeof(int) * 2)
            if reps[j][1]==NULL:
                raise MemoryError
            reps[j][0][0]=G._coset_reps_v0[j][0]
            reps[j][0][1]=G._coset_reps_v0[j][1]
            reps[j][1][0]=G._coset_reps_v0[j][2]
            reps[j][1][1]=G._coset_reps_v0[j][3]
    
    #if verbose>2:
    #    print "nc=",nc
    if not isinstance(G._cusp_data[0]['width'],(int,Integer)):
        use_int=0
    cdef int** normalizers=NULL
    cdef int** cusp_maps=NULL
    normalizers=<int**>sig_malloc(sizeof(int*)*nc)
    if normalizers==NULL: raise MemoryError
    for i in range(nc):
        normalizers[i]=<int*>sig_malloc(sizeof(int)*4)
        a,b,c,d=G._cusp_data[i]['normalizer']
        normalizers[i][0]=a
        normalizers[i][1]=b
        normalizers[i][2]=c
        normalizers[i][3]=d
        if i>0 and d<>0:
            raise ArithmeticError,"Can not use this symmetrized version of the algorithm!"
        if verbose>2:
            print "Normalizer[",i,"]=",a,b,c,d
    vertex_maps=<int**>sig_malloc(sizeof(int*)*nv)
    if vertex_maps==NULL: raise MemoryError
    cusp_maps=<int**>sig_malloc(sizeof(int*)*nv)
    if cusp_maps==NULL: raise MemoryError
    vertex_widths=<double*>sig_malloc(sizeof(double)*nv)
    if vertex_widths==NULL: raise MemoryError
    vertex_cusp=<int*>sig_malloc(sizeof(int)*nv)    
    if vertex_cusp==NULL: raise MemoryError
    if cusp_maps==NULL: raise MemoryError
    cdef int delta
    delta = 0
    for i in range(nv):
        cusp_maps[i]=<int*>sig_malloc(sizeof(int)*4)
        cusp_maps[i][0]=G._cusp_maps[i][0]
        cusp_maps[i][1]=G._cusp_maps[i][1]
        cusp_maps[i][2]=G._cusp_maps[i][2]
        cusp_maps[i][3]=G._cusp_maps[i][3]
        vertex_maps[i]=<int*>sig_malloc(sizeof(int)*4)
        vertex_maps[i][0]=<int>G._vertex_maps[i].a()
        vertex_maps[i][1]=<int>G._vertex_maps[i].b()
        vertex_maps[i][2]=<int>G._vertex_maps[i].c()
        vertex_maps[i][3]=<int>G._vertex_maps[i].d()
        vertex_widths[i]=<double>G._vertex_widths[i]
        vertex_cusp[i]=<int>G._vertex_data[i]['cusp']
    cdef double* widths
    widths=<double*>sig_malloc(sizeof(double)*nc)
    for i in range(nc):
        if verbose>3:
            print "width[",i,"]=",G._cusp_data[i]['width'],type(G._cusp_data[i]['width'])
        widths[i]=<double>G._cusp_data[i]['width']
    cdef int wi,wj,vi,vj
    cdef SL2Z_elt A,Tj
    cdef double dp_x,dp_y
    #print "Y={0}".format(Y)
    for ci in range(nc):
        wi = <int>widths[ci] #G._cusp_data[ci]['width'] #(ci)
        if verbose>1:
            print "cusp =",ci
        if wi<>1:
            mpfr_set_ui(swi,wi,rnd_re)
            mpfr_sqrt(swi,swi,rnd_re)
        else:
            mpfr_set_ui(swi,1,rnd_re)
        for j in range(Ql): #Qs,Qf+1): #from Qs<= j <= Qf:
            #if ciindex==0:
            mpfr_set(x0,Xm[j],rnd_re)
            mpfr_set(y0,YY,rnd_re)
            mpfr_set(y1.value,YY,rnd_re)
            #print "Xm[j]={0}, Y={1}".format(Xm[j],y1)
            if ci<>0:
                a=normalizers[ci][0]; b=normalizers[ci][1]
                c=normalizers[ci][2]; d=normalizers[ci][3]
                _normalize_point_to_cusp_mpfr(x0,y0,a,b,c,d,wi)
            if verbose > 2:
                mpfr_set(x1.value,x0,rnd_re)
                mpfr_set(y1.value,y0,rnd_re)
                print "x1={0},y1={1}".format(x1,y1)
            pullback_to_Gamma0N_mpfr_c(G,x1.value,y1.value,x0,y0,&a,&b,&c,&d)
            Tj=SL2Z_elt(a,b,c,d)
            if verbose > 2:
                print "pull back of {0},{1} is {2},{3}".format(mpfr_get_d(Xm[j],rnd_re),mpfr_get_d(Y,rnd_re),x1,y1)
                print "map:",a,b,c,d
            #dp_x=mpfr_get_Ad(x1,rnd_re)
            dp_x = float(x1); dp_y = float(y1)
            #dp_y=mpfr_get_d(y1,rnd_re)
            vj = closest_vertex_dp_c(nv,vertex_maps,vertex_widths,&dp_x,&dp_y)
            cj = vertex_cusp[vj]
            #cj,vj= G.closest_cusp(x1,y1,vertex=1)
            #cja,cjb= G.closest_cusp(x1,y1) #G._vertex_data[G._vertices[vi]]['cusp'] 
            #cjindex=G._cusps.index((cja,cjb))
            wj=<int>widths[cj] #G._cusp_data[cj]['width']
            if wj<>1:
                mpfr_set_ui(swj,wj,rnd_re)
                mpfr_sqrt(swj,swj,rnd_re)
            else:
                mpfr_set_ui(swj,1,rnd_re)
                #swj=sqrt(<double>wj)
            #wj = G._cusp_data[(cja,cjb)]['width']
            a = cusp_maps[vj][0]   #a,b,c,d=G._vertex_data[cj]['cusp_map'] #[v]
            b = cusp_maps[vj][1]
            c = cusp_maps[vj][2]
            d = cusp_maps[vj][3]  #U = G._vertex_data[vj]['cusp_map']
            mpfr_set(x2.value,x1.value,rnd_re)
            mpfr_set(y2.value,y1.value,rnd_re)
            #x2=x1; y2=y1
            if a<>1 or b<>0 or c<>0 or d<>1: #U<>SL2Z_elt(1,0,0,1):
                _apply_sl2z_map_mpfr(x2.value,y2.value,a,b,c,d) #U[0,0],U[0,1],U[1,0],U[1,1])
            if cj>0:
                mpfr_set(x3.value,x2.value,rnd_re)
                mpfr_set(y3.value,y2.value,rnd_re)
                #x3=x2; y3=y2
                a=normalizers[cj][3]; b=-normalizers[cj][1]
                c=-normalizers[cj][2]; d=normalizers[cj][0]
                if verbose>2:
                    print "Normalizer: a,b,c,d=",a,b,c,d
                _normalize_point_to_cusp_mpfr(x3.value,y3.value,a,b,c,d,wj,1)
            else:
                mpfr_set(x3.value,x2.value,rnd_re)
                mpfr_set(y3.value,y2.value,rnd_re)
                #x3=x2; y3=y2
            #x3 = x3*twopi
            mpfr_set(Xpb[ci][cj][j],x3.value,rnd_re)
            mpfr_mul(Xpb[ci][cj][j],Xpb[ci][cj][j],twopi.value,rnd_re)
            # Recall that Ypb must be greater than Y otherwise we need a better Y
            tmp = RF(y3-ep0)
            if mpfr_cmp(tmp.value,Y)>=0:
                mpfr_set(Ypb[ci][cj][j],y3.value,rnd_re)
            else:
                if verbose > 0:
                    print "ci,cj=",ci,cj
                    print "Xm=",mpfr_get_d(Xm[j],rnd_re)
                    #print "x,y=",x0,"\n",y0
                    print "x1,y1=",x1,"\n",y1
                    print "x2,y2=",x2,"\n",y2
                    print "x3,y3=",x3,"\n",y3
                    print "Xpb=",mpfr_get_d(Xpb[ci][cj][j],rnd_re)
                    print "Ypb=",mpfr_get_d(Ypb[ci][cj][j],rnd_re)
                #raise ArithmeticError,"Need smaller value of Y"
                return 1
            if verbose > 2:
                print "ci,cj=",ci,cj
                print "Xm=",mpfr_get_d(Xm[j],rnd_re)
                print "x,y=",mpfr_get_d(x0,rnd_re),"\n",mpfr_get_d(y0,rnd_re)
                print "x1,y1=",x1,"\n",y1
                print "x2,y2=",x2,"\n",y2
                print "x3,y3=",x3,"\n",y3
                print "Xpb=",mpfr_get_d(Xpb[ci][cj][j],rnd_re)
                print "Ypb=",mpfr_get_d(Ypb[ci][cj][j],rnd_re)
            # We also get the multiplier if we need it
            if(non_trivial):
                delta = ch_m1
                if verbose>2:
                    print "delta0=",delta

                if weight<>0:
                    c = normalizers[ci][2]
                    d = normalizers[ci][3]
                    mpfr_mul_si(cr.value,swi,c,rnd_re)
                    mpfr_set_si(dr.value,d,rnd_re)
                    mpfr_div(dr.value,dr.value,swi,rnd_re)
                    #m1=j_fak_mpc(cr,dr,Xm[j],Y,-weight,ui)
                    j_fak_mpc_c(cr.value,dr.value,Xm[j],YY,weight_t,ui,tmpab1,tmpar1)
                    if c<>0:  ## 
                        delta -= k
                    if verbose>2:
                        print "delta1=",delta
                    c = normalizers[cj][2]
                    d = normalizers[cj][3]
                    if c<>0:  ## 
                        delta += k
                    #else:
                    #    delta += 0
                    if verbose>2:
                        print "delta2=",delta
                    if verbose>2:
                        print "c=",c
                        print "d=",d
                        print "weight=",weight
                    mpfr_mul_si(cr.value,swj,c,rnd_re)
                    mpfr_set_si(dr.value,d,rnd_re)
                    mpfr_div(dr.value,dr.value,swj,rnd_re)
                    if verbose>2:
                        print "cr=",cr
                        print "dr=",dr
                        print "x3,y3=",x3,y3
                        print "weight=",weight
                        print "ui=",ui
                    #m2=j_fak_mpc(cr,dr,x3,y3,weight,ui)
                    j_fak_mpc_c(cr.value,dr.value,x3.value,y3.value,weight_t,ui,tmpab2,tmpar2)
                    if verbose>2:
                        print "m2=",m2
                    A = G._vertex_data[vj]['cusp_map']*Tj
                    c = -A.c(); d= A.a()
                    #if c<>0:  ## 
                    delta += k
                    if verbose>2:
                        print "delta3=",delta
                    #else:
                    #    delta += 0
                    cr=RF(c); dr=RF(d)
                    j_fak_mpc_c(cr.value,dr.value,x2.value,y2.value,weight_t,ui,tmpab3,tmpar3)
                    #m3=j_fak_mpc(cr,dr,x2,y2,weight,ui)
                    #ctmp=m1*m2*m3
                    mpfr_sub(tmpar,tmpar2,tmpar1,rnd_re)
                    mpfr_add(tmpar,tmpar,tmpar3,rnd_re)
                    mpfr_div(tmpab,tmpab2,tmpab1,rnd_re)
                    mpfr_mul(tmpab,tmpab,tmpab3,rnd_re)                    
                    #if verbose>2:
                    #    #mpfr_set(m1.value,)
                    #    print "m1=",m1
                    #    print "m2=",m2
                    #    print "m3=",m3
                else:
                    #A=U.matrix()*Tj
                    if not trivial_mult :
                        A = G._vertex_data[vj]['cusp_map']*Tj
                    #print "A:",A,type(A)
                    #ctmp=CF(1)
                    mpfr_set_si(tmpar,0,rnd_re)
                    mpfr_set_si(tmpab,1,rnd_re)
                if not trivial_mult :
                    #if multiplier<>None and not multiplier.is_trivial():
                    AA=A**-1
                    #v = S.character(AA[1,1]).complex_embedding(CF.prec())
                    mA = multiplier(AA)
                    if hasattr(mA,"complex_embedding"):
                        v = CF(mA.complex_embedding(CF.prec()))
                        rtmp = v.argument()
                        mpfr_add(tmpar,tmpar,rtmp.value,rnd_re)
                        rtmp = v.abs()
                        mpfr_mul(tmpab,tmpab,rtmp.value,rnd_re)
                    else: ## v is 1 or -1
                        if mA==-1:
                            mpfr_add(tmpar,tmpar,mppi,rnd_re)
                        #v = CF(mA)
                    #v=CF(v)
                    if verbose>2:
                        print "icusp,jcusp,j=",ci,cj,":",j
                        print "Tj=",Tj
                        print "mA=",mA
                        print "ctmp=",ctmp
                        print "delta=",delta
                        mpfr_set(rtmp.value,tmpab,rnd_re)
                        print "|v|=",rtmp
                        mpfr_set(rtmp.value,tmpar,rnd_re)
                        print "arg(v)=",rtmp
                        print "Y=",mpfr_get_d(Y,rnd_re)
                    #ctmp = ctmp*v
                delta = delta % 2
                mpfr_set(RCvec[ci][cj][j][0],tmpab,rnd_re)
                mpfr_set(RCvec[ci][cj][j][1],tmpar,rnd_re)
                #mpfr_set_si(RCvec[ci][cj][j][2],delta,rnd_re)
                CSvec[ci][cj][j]=delta % 2
            else:
                mpfr_set_ui(RCvec[ci][cj][j][0],1,rnd_re)
                mpfr_set_si(RCvec[ci][cj][j][1],0,rnd_re)
                CSvec[ci][cj][j]=0
                #mpfr_set_si(RCvec[ci][cj][j][2],0,rnd_re)                

    for j in range(Ql):
        mpfr_mul(Xm[j],Xm[j],twopi.value,rnd_re)
    if normalizers<>NULL:
        for i in range(nc):
            if  normalizers[i]<>NULL:
                sig_free(normalizers[i])
        sig_free(normalizers)
    mpfr_clear(weight_t)
    mpfr_clear(tmpar);    mpfr_clear(tmpab)
    mpfr_clear(tmpar1);    mpfr_clear(tmpab1)
    mpfr_clear(tmpar2);    mpfr_clear(tmpab2)
    mpfr_clear(tmpar3);    mpfr_clear(tmpab3)
    mpfr_clear(mppi)
    mpfr_clear(swi); mpfr_clear(swj)
    mpfr_clear(x0);  mpfr_clear(y0)
    mpfr_clear(YY)
    return 0

@cython.cdivision(True)
cdef int pullback_pts_hecke_triangle_mpc_new_c(S,int Qs,int Qf,RealNumber Y, mpfr_t* Xm,mpfr_t*** Xpb, mpfr_t*** Ypb,mpc_t ***Cvec):

    r""" Computes a whole array of pullbacked points-
         using MPFR/MPC types. For Hecke triangle groups

    INPUT:

    - ``S``  -- AutomorphicFormSpace
    - ``Qs`` -- integer
    - ``Qf`` -- integer
    - ``Y``  -- real >0
    - ``deb``    -- (False) logical

    
    OUTPUT:

    - ``pb`` -- dictionary with entries:
       - 'xm'   -- real[Qf-Qs+1]  : x_m=2*pi*(1-2*m)/2(Qf-Qs+1)
       - 'xpb'  -- real[0:nc,0:nc,Qf-Qs+1]  : real part of pullback of x_m+iY
       - 'ypb'  -- real[0:nc,0:nc,Qf-Qs+1]  : imag. part of pullback of x_m+iY
       - 'cvec' -- complex[0:nc,0:nc,Qf-Qs+1] : mult. factor

    EXAMPLES::


        sage: G=MySubgroup(Gamma0(2))
        sage: [Xm,Xpb,Ypb,Cv]=pullback_pts(G,1,3,0.1)   

    Note that we need to have $Y<Y_{0}$ so that all pullbacked points
    are higher up than Y.

        sage: pullback_pts(G,1,10,0.44)
        Traceback (most recent call last):
        ...
        ArithmeticError: Need smaller value of Y

    
    """
    # We let the input determine the precision
    #print "pullback mpc!"
    cdef int prec=Y.parent().prec()
    CF=MPComplexField(prec)
    RF=RealField(prec)
    cdef RealNumber tmpx,tmpy,weight,tmpc,tmpd
    cdef MPComplexNumber tmpz
    tmpz=CF(0)
    tmpx=RF(0); tmpy=RF(0); tmpc=RF(0); tmpd=RF(0)
    cdef mpfr_t a,b,c,d,x0,y0,lambdaq,Yep,twopi
    mpfr_init2(twopi,prec);mpfr_init2(Yep,prec)
    mpfr_init2(a,prec);    mpfr_init2(b,prec)
    mpfr_init2(c,prec);    mpfr_init2(d,prec)
    mpfr_init2(x0,prec);   mpfr_init2(y0,prec)
    mpfr_init2(lambdaq,prec)
    cdef int ui,Ql
    cdef RealNumber ep0
    Ql = Qf - Qs + 1
    ui = S._unitary_action
    weight = RF(S.weight())
    G = S.group()
    #if weight<>0 or not multiplier.is_trivial():
    #    non_trivial=True
    #else:
    #    non_trivial=False
    ep0=RF(1E-10)
    if S._holomorphic and not S._weak and weight==0:
        raise Exception,"Must support a non-zero weight for non-weak holomorphic forms!"
    verbose = S._verbose
    print "now 1"
    if Qs<0:
        Qfak = 2*(Qf-Qs+1)
        for j in range(Qs,Qf+1):
            mpfr_set_si(Xm[j-Qs],2*j-1,rnd_re)
            mpfr_div_si(Xm[j-Qs],Xm[j-Qs],Qfak,rnd_re)
    else:
        Qfak = 4*(Qf-Qs+1)
        for j in range(Qs,Qf+1):
            mpfr_set_si(Xm[j-Qs],2*j-1,rnd_re)
            mpfr_div_si(Xm[j-Qs],Xm[j-Qs],Qfak,rnd_re)
    print "now 2"            
    cdef int nc=G._ncusps
    cdef int nv=G._nvertices
    if nc>1 or nv>1:
        raise ValueError,"Need to be called with a Hecke triangle group!"
    cdef int q
    try: 
        q = <int>G._q
    except AttributeError:
        log.critical("Routine called with non-Hecke-triangle group!")
        q = 3 ## the modular group
    mpfr_const_pi(lambdaq,rnd_re)
    mpfr_div_ui(lambdaq,lambdaq,q,rnd_re)
    mpfr_cos(lambdaq,lambdaq,rnd_re)
    mpfr_mul_ui(lambdaq,lambdaq,2,rnd_re)
    mpfr_set(tmpx.value,lambdaq,rnd_re)
    # Recall that the width is lambda
    mpfr_sub(Yep,Y.value,ep0.value,rnd_re)
    mpfr_const_pi(twopi,rnd_re)
    mpfr_mul_ui(twopi,twopi,2,rnd_re)
    if verbose>0:
        print "lambda=",tmpx
    print "now 3"
    for j in range(Ql): #Qs,Qf+1): #from Qs<= j <= Qf:
        print j
        mpfr_set(x0,Xm[j],rnd_re)
        mpfr_set(y0,Y.value,rnd_re)
        mpfr_mul(x0,x0,lambdaq,rnd_re)
        mpfr_mul(y0,y0,lambdaq,rnd_re)
        ## If we want to normalize points we need z-> z /lambda
        pullback_to_hecke_triangle_mat_c_mpfr(x0,y0,lambdaq,a,b,c,d)
        # If we want to normalize points we need z0-> z0 /lambda
        # Recall that Ypb must be greater than Y otherwise we need a better Y
        if mpfr_cmp(y0,Yep)<=0:
            return 1
            #if verbose > 0:
        #        print "Xm,Y=",Xm[j],Y
        #        mpfr_set(tmpx.value,x0,rnd_re)
        #        mpfr_set(tmpy.value,y0,rnd_re)
        #        print "xpb,ypb=",tmpx,"\n",tmpy
        #raise ArithmeticError,"Need smaller value of Y" 
                # We also get the multiplier if we need it
        mpfr_div(x0,x0,lambdaq,rnd_re)
        mpfr_div(y0,y0,lambdaq,rnd_re)
        mpfr_set(Ypb[0][0][j],y0,rnd_re)
        mpfr_set(Xpb[0][0][j],x0,rnd_re)
        mpfr_mul(Xpb[0][0][j],Xpb[0][0][j],twopi,rnd_re)
        print "in here"
        #if verbose > 2:
        #    print "Xm,Y=",Xm[j],Y
        #    mpfr_set(tmpx.value,x0,rnd_re)
        #    mpfr_set(tmpy.value,y0,rnd_re)
        #    print "xpb,ypb=",tmpx,"\n",tmpy
        ## We didn't implement non-trivial mul
        ## but we will soon allow for non-zero weight...
        if weight<>0:
            print "weight !=0"
            mpfr_set(tmpx.value,x0,rnd_re)
            mpfr_set(tmpy.value,y0,rnd_re)
            mpfr_set(tmpc.value,c,rnd_re)
            mpfr_set(tmpd.value,d,rnd_re)
            print "j_fak+"
            tmpz=j_fak_mpc(tmpc,tmpd,tmpx,tmpy,weight,ui)
            mpc_set(Cvec[0][0][j],tmpz.value,rnd)
        else:
            mpc_set_ui(Cvec[0][0][j],1,rnd)
        print "end"
    print "now 4"
    for j in range(Ql):
        mpfr_mul(Xm[j],Xm[j],twopi,rnd_re)
    mpfr_clear(x0);  mpfr_clear(y0)
    mpfr_clear(Yep); mpfr_clear(lambdaq);mpfr_clear(twopi)
    mpfr_clear(a);mpfr_clear(b);mpfr_clear(c);mpfr_clear(d)
    return 0
# @cython.cdivision(True)
# cdef pullback_pts_mpc_new2(S,int Qs,int Qf,RealNumber Y,RealNumber weight, Vector_real_mpfr_dense Xm,mpfr_t*** Xpb, mpfr_t*** Ypb,mpc_t ***Cvec):
#     r""" Computes a whole array of pullbacked points-
#          using MPFR/MPC types.
#     INPUT:

#     - ``S``  -- AutomorphicFormSpace
#     - ``Qs`` -- integer
#     - ``Qf`` -- integer
#     - ``Y``  -- real >0
#     - ``deb``    -- (False) logical

    
#     OUTPUT:

#     - ``pb`` -- dictionary with entries:
#        - 'xm'   -- real[Qf-Qs+1]  : x_m=2*pi*(1-2*m)/2(Qf-Qs+1)
#        - 'xpb'  -- real[0:nc,0:nc,Qf-Qs+1]  : real part of pullback of x_m+iY
#        - 'ypb'  -- real[0:nc,0:nc,Qf-Qs+1]  : imag. part of pullback of x_m+iY
#        - 'cvec' -- complex[0:nc,0:nc,Qf-Qs+1] : mult. factor

#     EXAMPLES::


#         sage: G=MySubgroup(Gamma0(2))
#         sage: [Xm,Xpb,Ypb,Cv]=pullback_pts(G,1,3,0.1)   

#     Note that we need to have $Y<Y_{0}$ so that all pullbacked points
#     are higher up than Y.

#         sage: pullback_pts(G,1,10,0.44)
#         Traceback (most recent call last):
#         ...
#         ArithmeticError: Need smaller value of Y

    
#     """
#     # WE let the input determine the precision
#     #print "pullback mpc!"
#     cdef int prec=Y.parent().prec()
#     CF=MPComplexField(prec)
#     RF=RealField(prec)
#     cdef MPComplexNumber tmpz,tmp
#     cdef RealNumber x1,y1,x2,y2,x3,y3,weight,tmpx,tmpy
#     x1=RealNumber(RF,1)
#     x2=RealNumber(RF,1)
#     x3=RealNumber(RF,1)
#     y1=RealNumber(RF,1)
#     y2=RealNumber(RF,1)
#     y3=RealNumber(RF,1)
#     cdef int a,b,c,d,v0,v1,itmp,i,j,ui
#     cdef mpfr_t swj,swi,x0,y0,YY
#     #cdef Vector_real_mpfr_dense Xm
#     cdef RealNumber ar,br,cr,dr,twopi,xm
#     mpfr_init2(swi,prec)
#     mpfr_init2(swj,prec)
#     mpfr_init2(x0,prec)
#     mpfr_init2(y0,prec)
#     mpfr_init2(YY,prec)
#     mpfr_set(YY,Y.value,rnd_re)
#     #mpfr_init2(x1,prec)
#     #mpfr_init2(y1,prec)
#     ui = S._unitary_action
#     weight = RF(S.weight())
#     G = S.group()
#     multiplier = S.multiplier()
#     if weight<>0 or not multiplier.is_trivial():
#         non_trivial=True
#     else:
#         non_trivial=False
#     trivial_mult = multiplier.is_trivial()
#     tmpx=RF(0)
#     tmpy=RF(0)
#     twopi=RF(2)*RF.pi() 
#     twopii=CF(0,twopi)
#     mp_i=CF(0,1)
#     cdef RealNumber one
#     one=RF(1)
#     tmpz=CF(0)
#     tmp=CF(0)
#     if S._holomorphic and not S._weak and weight==0:
#         raise Exception,"Must support a non-zero weight for non-weak holomorphic forms!"
#     verbose = S._verbose
#     if verbose>3:
#         fp0=open("xm.txt","w")
#         fp1=open("xpb.txt","w")
#         fp2=open("ypb.txt","w")
#         fp3=open("cv.txt","w")
#     #twoQl=RF(2*(Qf+1-Qs))
#     if Qs<0:
#         Qfak = 2*(Qf-Qs+1)
#         for j in range(Qs,Qf+1):
#             mpfr_set_si(Xm._entries[j-Qs],2*j-1,rnd_re)
#             mpfr_div_si(Xm._entries[j-Qs],Xm._entries[j-Qs],Qfak,rnd_re)
#     else:
#         Qfak = 4*(Qf-Qs+1)
#         for j in range(Qs,Qf+1):
#             mpfr_set_si(Xm._entries[j-Qs],2*j-1,rnd_re)
#             mpfr_div_si(Xm._entries[j-Qs],Xm._entries[j-Qs],Qfak,rnd_re)
#     cdef tuple ci,cj
#     cdef int cia,cib,cja,cjb,ciindex,cjindex,nc
#     nc = len(G._cusps)
#     #cdef RealNumber swi
#     cdef int** normalizers=NULL
#     normalizers = <int **>sig_malloc (sizeof(int*)*nc)
#     if normalizers==NULL:
#         raise MemoryError
#     for i from 0 <= i < nc:
#         normalizers[i]=NULL
#         normalizers[i]=<int *> sig_malloc (sizeof(int)*4)
#         if normalizers[i]==NULL:
#             raise MemoryError
#         normalizers[i][0]=<int>G._cusp_data[G._cusps[i]]['normalizer'][0,0]
#         normalizers[i][1]=<int>G._cusp_data[G._cusps[i]]['normalizer'][0,1]
#         normalizers[i][2]=<int>G._cusp_data[G._cusps[i]]['normalizer'][1,0]
#         normalizers[i][3]=<int>G._cusp_data[G._cusps[i]]['normalizer'][1,1]
#         #print "normalizer ",i,":", normalizers[i][0],normalizers[i][1],normalizers[i][2],normalizers[i][3]
#     cdef int wi,wj
    
#     MS = sage.matrix.matrix_integer_2x2.MatrixSpace_ZZ_2x2()
#     cdef Matrix_integer_2x2 Tj,A
#     for cia,cib in G._cusps:
#         ciindex=G._cusps.index((cia,cib))
#         wi = G._cusp_data[(cia,cib)]['width']
#         if verbose>1:
#             print "cusp =",cia,cib
#         if wi<>1:
#             mpfr_set_ui(swi,wi,rnd_re)
#             mpfr_sqrt(swi,swi,rnd_re)
#         else:
#             mpfr_set_ui(swi,1,rnd_re)
#         for j from Qs<= j <= Qf:
#             #if ciindex==0:
#             mpfr_set(x0,Xm._entries[j-Qs],rnd_re)
#             mpfr_set(y0,Y.value,rnd_re)
#             if ciindex<>0:
#                 normalize_point_to_cusp_mpfr_c(x0,y0,normalizers[ciindex],wi,x0,y0,0)
#             pullback_to_Gamma0N_mpfr_c(G,x1.value,y1.value,x0,y0,&a,&b,&c,&d)
#             Tj=MS([a,b,c,d])
#             #vi = G.closest_vertex(x1,y1)
#             cja,cjb= G.closest_cusp(x1,y1)  #G._vertex_data[G._vertices[vi]]['cusp'] 
#             cjindex=G._cusps.index((cja,cjb))
#             wj = G._cusp_data[(cja,cjb)]['width']
#             U = G._vertex_data[v]['cusp_map']
#             if U<>SL2Z_elt(1,0,0,1):
#                 x2=x1; y2=y1
#                 _apply_sl2z_map_mpfr(x2.value,y2.value,U[0,0],U[0,1],U[1,0],U[1,1])
#             else:
#                 x2=x1; y2=y1;
#             #[x3,y3] = normalize_point_to_cusp(G,x2,y2,cj,inv=True,mp_ctx=mp_ctx)
#             if cjindex<>0:
#                 normalize_point_to_cusp_mpfr_c(x3.value,y3.value,normalizers[cjindex],wj,x2.value,y2.value,1)
#                 #[x3,y3] = normalize_point_to_cusp_mpfr(G,cja,cjb,x2,y2,inv=1)
#             else:
#                 x3=x2
#                 y3=y2
#              # Recall that Ypb must be greater than Y otherwise we need a better Y
#             if y3>=Y:
#                 #Ypb[ciindex,cjindex,j]=y3*one ## Ridiculous but necessary due to bug in cython...
#                 mpfr_set(Xpb[ciindex][cjindex][j-Qs],x3.value,rnd_re)
#                 mpfr_mul(Xpb[ciindex][cjindex][j-Qs],Xpb[ciindex][cjindex][j-Qs],twopi.value,rnd_re)
#                 mpfr_set(Ypb[ciindex][cjindex][j-Qs],y3.value,rnd_re)
#                 #mpfr_mul(Zpb[ciindex][cjindex[j][1],Zpb[ciindex][cjindex[j][0],twopi.value,rnd_re)
                                      
#             else:
#                 if verbose > 0:
#                     print "ci,cj=",cia,cib,":",cja,cjb
#                     print "Xm=",Xm[j-Qs]
#                     #print "x,y=",x0,"\n",y0
#                     print "x1,y1=",x1,"\n",y1
#                     print "x2,y2=",x2,"\n",y2
#                     print "x3,y3=",x3,"\n",y3
#                     mpfr_set(tmpx.value,Xpb[ciindex][cjindex][j-Qs],rnd_re)
#                     mpfr_set(tmpy.value,Ypb[ciindex][cjindex][j-Qs],rnd_re)
#                     print "Xpb=",tmpx
#                     print "Ypb=",tmpy
#                 raise ArithmeticError,"Need smaller value of Y" 
#             if verbose > 2:
#                 print ciindex,cjindex,j
#                 print "ci,cj=",cia,cib,":",cja,cjb
#                 print "Xm=",Xm[j-Qs]
#                 print "x,y=",mpfr_get_d(x0,rnd_re),"\n",mpfr_get_d(y0,rnd_re)
#                 print "x1,y1=",x1,"\n",y1
#                 print "x2,y2=",x2,"\n",y2
#                 print "x3,y3=",x3,"\n",y3
#                 mpfr_set(tmpx.value,Xpb[ciindex][cjindex][j-Qs],rnd_re)
#                 mpfr_set(tmpy.value,Ypb[ciindex][cjindex][j-Qs],rnd_re)
#                 print "Xpb=",tmpx
#                 print "Ypb=",tmpy
#             # We also get the multiplier if we need it
#             if(non_trivial):
#                 if weight<>0:
#                     c = G._cusp_data[(cia,cib)]['normalizer'][1,0]
#                     d = G._cusp_data[(cia,cib)]['normalizer'][1,0]
#                     mpfr_mul_si(cr.value,swi,c,rnd_re)
#                     mpfr_set_si(dr.value,d,rnd_re)
#                     mpfr_div(dr.value,dr.value,swi,rnd_re)
#                     #cr = RF(G._cusp_data[ci]['normalizer'][1,0])*swi
#                     #dr = RF(G._cusp_data[ci]['normalizer'][1,1])/swi
#                     #print "test:",cr,dr,Xm[j],Y,-weight
#                     #print type(c),type(d)
#                     m1=j_fak_mpc(cr,dr,Xm[j-Qs],Y,-weight,ui)
#                     c = G._cusp_data[cjindex]['normalizer'][1,0]
#                     d = G._cusp_data[cjindex]['normalizer'][1,1]
#                     mpfr_mul_si(cr.value,swj,c,rnd_re)
#                     mpfr_set_si(dr.value,d,rnd_re)
#                     mpfr_div(dr.value,dr.value,swj,rnd_re)
#                     #cr = RF(G._cusp_data[cj]['normalizer'][1,0])*swj
#                     #dr = RF(G._cusp_data[cj]['normalizer'][1,1])/swj
#                     #print "test:",cr,dr,x3,y3,-weight
#                     m2=j_fak_mpc(cr,dr,x3,y3,weight,ui)
#                     #aa = U[0,0]*a+U[0,1]*c
#                     A=U.matrix()*Tj
#                     #A=A**-1
#                     cr=RF(-A[1,0]); dr=RF(A[0,0])
#                     m3=j_fak_mpc(cr,dr,x2,y2,weight,ui)
#                     tmp=m1*m2*m3
#                     #print "m1,m2,m3=",m1,m2,m3
#                 else:
#                     #print "U:",U.matrix(),type(U.matrix())
#                     #print "Tj:",Tj,type(Tj)
#                     A=U.matrix()*Tj
#                     #print "A:",A,type(A)
#                     tmp=CF(1)
#                 if not trivial_mult :
#                     #if multiplier<>None and not multiplier.is_trivial():
#                     AA=A**-1
#                     #v = S.character(AA[1,1]).complex_embedding(CF.prec())
#                     v = CF(multiplier(AA).complex_embedding(CF.prec()))
#                     #v=CF(v)
#                     tmpz = CF(tmp)*v
#                 #Cvec[ciindex,cjindex,j]=tmp # mp_ctx.mpc(0,tmp)
#                 mpc_set(Cvec[ciindex][cjindex][j-Qs],tmpz.value,rnd)
#                 #print "C[",cii,cjj,j,"]=",tmp
#             else:
#                 mpc_set_ui(Cvec[ciindex][cjindex][j-Qs],1,rnd)
#             if verbose>3:
#                 fp1.write(str(x3)+"\n")
#                 fp2.write(str(y3)+"\n")
#                 fp3.write(str(tmp)+"\n")

#         #        Xmi[j]=mp_ctx.mpc(0,Xm[j]*twopi)
#     for j in range(Qs,Qf+1): 
#         mpfr_mul(Xm._entries[j-Qs],Xm._entries[j-Qs],twopi.value,rnd_re)
#     if verbose>3:
#         fp1.close()
#         fp2.close()
#         fp3.close()
#     pb=dict()
#     if normalizers<>NULL:
#         for i from 0 <= i < nc:
#             if  normalizers[i]<>NULL:
#                 sig_free(normalizers[i])
#     sig_free(normalizers)            
#     #pb['xm']=Xm; pb['xpb']=Xpb; pb['ypb']=Ypb; pb['cvec']=Cvec
#     mpfr_clear(swi); mpfr_clear(swj)
#     mpfr_clear(x0);  mpfr_clear(y0)
#     mpfr_clear(YY)
#     #return pb


cpdef j_fak_old(RealNumber c,RealNumber d,RealNumber x,RealNumber y,RealNumber k,int unitary=1):
    CF=MPComplexField(x.prec())
    return (c*CF(x,y)+d)**k


cpdef j_fak_dp(double c,double d,double x,double y,double k,int unitary=1):
    cdef double tmpabs,tmparg
    j_fak_dp_c(c,d,x,y,k,unitary,&tmpabs,&tmparg)
    if unitary==0:
        return tmpabs,tmparg
    else:
        return tmparg
    
cdef void j_fak_dp_c(double c,double d,double x,double y,double k,int unitary,double *tmpabs,double *tmparg):
    r"""
    Computes the argument (with principal branch)  of the
    automorphy factor j_A(z;k) defined by:
    For A=(a b ; c d ) we have j_A(z)=(cz+d)^k=|cz+d|^k*exp(i*k*Arg(cz+d))

    INPUT:

    - ``c`` -- integer or real
    - ``d`` -- integer or real
    - ``x`` -- real 
    - ``y`` -- real  > 0
    - ``k`` -- real 
    - ``unitary`` -- (default True) logical
               = True  => use the unitary version (i.e. for Maass forms)
               = False => use the non-unitary version (i.e. for holomorphic forms)


    OUTPUT:

    - ``t``     -- real: t = k*Arg(cz+f) (if unitary==True)
    - ``[r,t]`` -- [real,real]: r=|cz+d|^k, t=k*Arg(cz+f) (if unitary==False)

    EXAMPLES::


    """
    cdef double complex tmp
    if  (c==0.0 and d==1.0) or (k==0.0):
        tmpabs[0]=1.0
        tmparg[0]=0.0
    else:
        tmp = c*x+d + _Complex_I * c*y
        tmparg[0]=carg(tmp)
        tmparg[0]=tmparg[0]*k
        if unitary==0:
            tmpabs[0]=cabs(tmp)
            tmpabs[0]=pow(tmpabs[0],k)


cpdef j_fak_mpc(RealNumber c,RealNumber d,RealNumber x,RealNumber y,RealNumber k,int unitary=1):
    cdef mpfr_t tmpabs,tmparg
    cdef mpfr_t xx,yy
    cdef RealNumber tmpab,tmpar
    cdef MPComplexNumber res
    cdef int prec = x.prec()
    mpfr_init2(tmpabs,prec)
    mpfr_init2(tmparg,prec)
    mpfr_init2(xx,prec)
    mpfr_init2(yy,prec)
    mpfr_set(xx,x.value,rnd_re)
    mpfr_set(yy,y.value,rnd_re)
    j_fak_mpc_c(c.value,d.value,xx,yy,k.value,unitary,tmpabs,tmparg)
    #mpfr_set(tmpab.value,tmpabs,rnd_re)
    #mpfr_set(tmpar.value,tmparg,rnd_re)
    mpfr_clear(xx)
    mpfr_clear(yy)
    CF = MPComplexField(prec)
    if unitary==0:
        res = CF(0)
        mpc_set_fr(res.value,tmparg,rnd)
        mpfr_swap(mpc_realref(res.value),mpc_imagref(res.value))
        mpc_exp(res.value,res.value,rnd)
        mpc_mul_fr(res.value,res.value,tmpabs,rnd)
        mpfr_clear(tmpabs)
        mpfr_clear(tmparg)
        return res
    elif unitary==1:
        tmpar=x.parent()(1)
        mpfr_set(tmpar.value,tmparg,rnd_re)
        mpfr_clear(tmpabs)
        mpfr_clear(tmparg)
        return CF(tmpar,0)
    elif unitary==2:
        tmpar=x.parent()(1); tmpab=x.parent()(1)
        mpfr_set(tmpar.value,tmparg,rnd_re)
        mpfr_set(tmpab.value,tmpabs,rnd_re)
        mpfr_clear(tmpabs)
        mpfr_clear(tmparg)
        return CF(tmpab,tmpar)  ## We want to return the same type always... 
    #res = tmpab*CF(0,tmpar).exp()
    #return res
    # return tmpab,tmpar       
    #else:
    #    return tmpar
 

    
        
cdef void j_fak_mpc_c(mpfr_t c,mpfr_t d,mpfr_t x,mpfr_t y,mpfr_t k,int unitary,mpfr_t tmpabs,mpfr_t tmparg):
    r"""
    Computes the argument (with principal branch)  of the
    automorphy factor j_A(z;k) defined by:
    For A=(a b ; c d ) we have j_A(z)=(cz+d)^k=|cz+d|^k*exp(i*k*Arg(cz+d))

    INPUT:

    - ``c`` -- integer or real
    - ``d`` -- integer or real
    - ``x`` -- real 
    - ``y`` -- real > 0
    - ``k`` -- real 
    - ``unitary`` -- (default True) logical
               = True  => use the unitary version (i.e. for Maass forms)
               = False => use the non-unitary version (i.e. for holomorphic forms)

    OUTPUT:

    - ``t``     -- real: t = k*Arg(cz+f) (if unitary==True)
    - ``[r,t]`` -- [real,real]: r=|cz+d|^k, t=k*Arg(cz+f) (if unitary==False)

    EXAMPLES::


    """
    cdef int prec = mpfr_get_prec(x)
    cdef mpc_t tmp
    cdef mpfr_t xx,yy
    mpfr_init2(xx,prec)
    mpfr_init2(yy,prec)
    mpc_init2(tmp,prec)
    if (mpfr_zero_p(c)<>0 and mpfr_cmp_ui(d,1)==0) or (mpfr_zero_p(k)<>0):
            mpfr_set_ui(tmparg,0,rnd_re)
            mpfr_set_ui(tmpabs,1,rnd_re)
    if (mpfr_zero_p(c)<>0 and mpfr_cmp_ui(d,-1)==0):
        mpfr_const_pi(tmparg,rnd_re)
        mpfr_mul(tmparg,tmparg,k,rnd_re)
        mpfr_set_ui(tmpabs,1,rnd_re)
    else:
        mpfr_mul(xx,x,c,rnd_re)
        mpfr_add(xx,xx,d,rnd_re)
        mpfr_mul(yy,y,c,rnd_re)
        mpc_set_fr_fr(tmp,xx,yy,rnd)
        mpc_arg(tmparg,tmp,rnd_re)
        mpfr_mul(tmparg,tmparg,k,rnd_re)
        if unitary==0:
            mpc_abs(tmpabs,tmp,rnd_re)
            mpfr_pow(tmpabs,tmpabs,k,rnd_re)
    mpc_clear(tmp); mpfr_clear(xx); mpfr_clear(yy)

cpdef j_fak_int_mpfr(int c,int d,RealNumber x,RealNumber y,RealNumber k,int unitary=0):
    cdef mpfr_t tmpabs,tmparg
    cdef mpfr_t xx,yy
    cdef mpc_t cparg
    cdef RealNumber tmpab,tmpar
    cdef MPComplexNumber res
    cdef int prec = x.prec()
    CF = MPComplexField(prec)
    mpfr_init2(tmpabs,prec)
    mpfr_init2(tmparg,prec)
    mpfr_init2(xx,prec)
    mpfr_init2(yy,prec)
    mpc_init2(cparg,prec)
    mpfr_set(xx,x.value,rnd_re)
    mpfr_set(yy,y.value,rnd_re)
    j_fak_int_mpfr_c(c,d,xx,yy,k.value,unitary,tmpabs,tmparg)
    #tmpab=x.parent()(1)
    #tmpar=x.parent()(1)
    # arg -> i*arg
    mpc_set_fr(cparg,tmparg,rnd)
    mpfr_swap(mpc_realref(cparg),mpc_imagref(cparg))
    res = CF(0)
    mpc_exp(cparg,cparg,rnd)
    mpc_set(res.value,cparg,rnd)
    #mpfr_set(tmpab.value,tmpabs,rnd_re)
    #mpfr_set(tmpar.value,tmparg,rnd_re)
    mpfr_clear(xx)
    mpfr_clear(yy)
    mpfr_clear(tmparg)
    mpc_clear(cparg)
    if unitary==0:
        mpc_mul_fr(res.value,res.value,tmpabs,rnd)
        mpfr_clear(tmpabs)
        #res = tmpab*CF(0,tmpar).exp()
        return res
    # return tmpab,tmpar
    else:
        #res = CF(0,tmpar).exp()
        return res



cdef void j_fak_int_mpfr_c(int c,int d,mpfr_t x,mpfr_t y,mpfr_t k,int unitary,mpfr_t tmpabs,mpfr_t tmparg):
    r"""
    Computes the argument (with principal branch)  of the
    automorphy factor j_A(z;k) defined by:
    For A=(a b ; c d ) we have j_A(z)=(cz+d)^k=|cz+d|^k*exp(i*k*Arg(cz+d))

    INPUT:

    - ``c`` -- integer 
    - ``d`` -- integer
    - ``x`` -- real 
    - ``y`` -- real > 0
    - ``k`` -- real 
    - ``unitary`` -- (default True) logical
               = True  => use the unitary version (i.e. for Maass forms)
               = False => use the non-unitary version (i.e. for holomorphic forms)

    OUTPUT:

    - ``t``     -- real: t = k*Arg(cz+f) (if unitary==True)
    - ``[r,t]`` -- [real,real]: r=|cz+d|^k, t=k*Arg(cz+f) (if unitary==False)

    EXAMPLES::


    """
    cdef int prec = mpfr_get_prec(x)
    cdef mpc_t tmp
    mpc_init2(tmp,prec)
    if (c==0 and d==1) or (mpfr_zero_p(k)<>0):
        mpfr_set_ui(tmparg,0,rnd_re)
        mpfr_set_ui(tmpabs,1,rnd_re)
    if (c==0 and d==-1):
        mpfr_const_pi(tmparg,rnd_re)
        mpfr_mul(tmparg,tmparg,k,rnd_re)
        mpfr_set_ui(tmpabs,1,rnd_re)
    else:
        mpfr_mul_si(x,x,c,rnd_re)
        mpfr_add_si(x,x,d,rnd_re)
        mpfr_mul_si(y,y,c,rnd_re)
        mpc_set_fr_fr(tmp,x,y,rnd)
        mpc_arg(tmparg,tmp,rnd_re)
        mpfr_mul(tmparg,tmparg,k,rnd_re)
        if unitary==0:
            mpc_abs(tmpabs,tmp,rnd_re)
            mpfr_pow(tmpabs,tmpabs,k,rnd_re)
    mpc_clear(tmp)



    #tmpabs=cabs(tmp)**k
    #return (tmpabs,tmparg)
    #else:
    #    mpfr_set_ui(tmpabs,1,mpfr_rnd_t)
    #    return tmparg



def j_fak_mpmath(c,d,x,y,k,unitary=False,mp_ctx=mpmath.mp):
    tmp=mp_ctx.mpc(c*x+d,c*y)
    tmparg=mp_ctx.arg(tmp)
    tmparg=tmparg*mp_ctx.mpf(k)
    if(not unitary):
        tmpabs=mp_ctx.power(mp_ctx.absmax(tmp),mp_ctx.mpf(k))
        return (tmpabs,tmparg)
    else:
        return tmparg



def pullback_pts_mpmath(H,Qs,Qf,Y):
    r""" Computes a whole array of pullbacked points-
    INPUT:

    - ``H``  -- Space of Harmonic Weak Maass forms
    - ``Qs`` -- integer
    - ``Qf`` -- integer
    - ``Y``  -- real (mpmath) 
    
    OUTPUT:
    
    - ``Xm``   -- real[Qf-Qs+1]  : x_m=2*pi*(1-2*m)/2(Qf-Qs+1)
    - ``Xpb``  -- real[0:nc,0:nc,Qf-Qs+1]  : real part of pullback of x_m+iY
    - ``Ypb``  -- real[0:nc,0:nc,Qf-Qs+1]  : imag. part of pullback of x_m+iY
    - ``Cvec`` -- complex[0:nc,0:nc,Qf-Qs+1] : mult. factor

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
    import mpmath
    GG=H._group
    holo=H._holomorphic
    mpmath_ctx=H._mp_ctx
    weight=mpmath_ctx.mpf(H._weight)
    RF = RealField(mpmath_ctx.prec)
    #print "mpmath_ctx=",mpmath_ctx
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
    if holo and weight==None and not H._weak:
        raise Exception,"Must support a weight for holomorphic forms!"
    elif holo:
        try:
            weight.ae
        except AttributeError:
            raise TypeError,"Need weight in mpmath format for pullback!"
    #Qs=1-Q; Qf=Q
    if(H._verbose>1):
        fp0=open("xm.txt","w")
        fp1=open("zpb.txt","w")
        #fp2=open("ypb.txt","w")
        #fp3=open("cv.txt","w")
    twoQl=mpmath_ctx.mpf(2*(Qf+1-Qs))
    #print "Qs,Qf=",Qs,Qf
    if(Qs<0):
        for j in range(Qs,Qf+1): #1-Q,Q):
            Xm[j]=mpmath_ctx.mpf(2*j-1)/twoQl        
            #if(H._verbose>1):
            #    s=str(Xm[j])+"\n"
            #    #fp0.write(s)
    else:
        for j in range(Qs,Qf+1): #1-Q,Q):
            Xm[j]=mpmath_ctx.mpf(2*j-1)/twoQl        
            #print "Xm[",j,"]=",Xm[j]
            #if(H._verbose>1):
            #    s=str(Xm[j])+"\n"
            #    #fp0.write(s)
            
    for ci in GG._cusps:
        cii=GG._cusps.index(ci)
        if(H._verbose>1):
            print "cusp :=",ci
        swi=mpmath_ctx.sqrt(mpmath_ctx.mpf(GG._cusp_data[ci]['width']))
        for j in range(Qs,Qf+1): #1-Q,Q):
            #print "here:",Xm[j],Y,ci,mpmath_ctx
            [x,y]   = normalize_point_to_cusp_mpmath(GG,Xm[j],Y,ci,mp_ctx=mpmath_ctx)
            #[x1,y1,Tj] =  pullback_to_G(GG,x,y,mp_ctx=mpmath_ctx)
            [x1,y1,Tj] =  GG.pullback(RF(x),RF(y),prec=mpmath_ctx.prec)
            x1=mpmath.mp.mpf(x1)
            y1=mpmath.mp.mpf(y1)
            #v0,v1 = GG.closest_vertex(x1,y1)
            cj=  GG.closest_cusp(x1,y1)  #GG._vertex_data[(v0,v1)]['cusp']
            cjj=GG._cusps.index(cj)
            swj=mpmath_ctx.sqrt(mpmath_ctx.mpf(GG._cusp_data[cj]['width']))
            U = GG._vertex_data[v]['cusp_map']
            if U<>SL2Z_elt(1,0,0,1):
                [x2,y2] = apply_sl2z_map(x1,y1,U,mp_ctx=mpmath_ctx)
            else:
                x2=x1; y2=y1;
            [x3,y3] = normalize_point_to_cusp_mpmath(GG,x2,y2,cj,inv=True,mp_ctx=mpmath_ctx)
            #Xpb[cii,cjj,j]=mpmath_ctx.mpc(0,x3*twopi)
            Xpb[cii,cjj,j]=x3*twopi
            # Recall that Ypb must be greater than Y otherwise we need a better Y
            if(y3>Y):
                Ypb[cii,cjj,j]=y3
            else:
                raise ArithmeticError,"Need smaller value of Y than %s" % Y 
            if(H._verbose>1 and j==-59):
                print "ci,cj=",ci,cj
                print "ci,cj=",cii,cjj
                print "X[",j,"]=",Xm[j]
                print "x,y=",x,"\n",y
                print "x1,y1=",x1,"\n",y1
                print "x2,y2=",x2,"\n",y2
                print "x3,y3=",x3,"\n",y3
                print "Xpb=",Xpb[cii,cjj,j]/twopi
                print "Ypb=",Ypb[cii,cjj,j]
            # We also get the multiplier if we need it
            if(True):
                c = mpmath_ctx.mpf(GG._cusp_data[ci]['normalizer'][1,0])*swi
                d = mpmath_ctx.mpf(GG._cusp_data[ci]['normalizer'][1,1])/swi
                m1=j_fak_mpmath(c,d,Xm[j],Y,-weight,unitary=False,mp_ctx=mpmath_ctx)
                c = mpmath_ctx.mpf(GG._cusp_data[cj]['normalizer'][1,0])*swj
                d = mpmath_ctx.mpf(GG._cusp_data[cj]['normalizer'][1,1])/swj
                m2=j_fak_mpmath(c,d,x3,y3,weight,unitary=False,mp_ctx=mpmath_ctx)
                A=(U*Tj)
                #A=A**-1
                c=mpmath_ctx.mpf(-A[1,0]); d=mpmath_ctx.mpf(A[0,0])
                m3=j_fak_mpmath(c,d,x2,y2,weight,unitary=False,mp_ctx=mpmath_ctx)
                tmp=m1[0]*m2[0]*m3[0]
                tmparg=m1[1]+m2[1]+m3[1]
                if(H._verbose>1 and j==-59):
                    print "m1=",m1
                    print "m2=",m2
                    print "m3=",m3
                    print "tmp=",tmp
                    print "tmparg=",tmparg
                #print "tmp=",tmparg
                #if(abs(tmparg)>mpmath.mpf(0)):
                tmp=tmp*mpmath_ctx.exp(mpmath_ctx.mpc(0,tmparg))
                if(H._rep):
                    AA=A**-1
                    v=H.rep(AA)
                    tmp=tmp*mpmath_ctx.mpc(v)
                    if(H._verbose>1 and j==-59):
                        print "v=",AA[1,0],AA[1,1],v
                        print "ch=",tmp
                Cvec[cii,cjj,j]=tmp # mpmath_ctx.mpc(0,tmp)
            if(H._verbose>1 and cii==0):
                #fp0.write(str(Xm[j])+" "+str(Y)+"\n")
                fp0.write(str(x1)+" "+str(y1)+"\n")
                if(cii==0):
                    fp1.write(str(x3)+" "+str(y3)+"\n")
                #fp3.write(str(tmp)+"\n")

    for j in range(Qs,Qf+1): #1-Q,Q):
        #        Xmi[j]=mpmath_ctx.mpc(0,Xm[j]*twopi)
        Xm[j]=Xm[j]*twopi
    if(H._verbose>1):
        fp0.close()
        fp1.close()
        #fp3.close()
    pb=dict()
    pb['xm']=Xm
    pb['xpb']=Xpb
    pb['ypb']=Ypb
    pb['cvec']=Cvec
    return pb
# return [Xm,Xpb,Ypb,Cvec]

cpdef mat_mul_list(list l1,list l2):
    cdef int a,b,c,d,x,y,z,w,res[4]
    a=int(l1[0]); b=int(l1[1]);c=int(l1[2]); d=int(l1[3])
    x=int(l2[0]); y=int(l2[1]);z=int(l2[2]); w=int(l2[3])
    _mat_mul_list(a,b,c,d,x,y,z,w,res)
    return [res[0],res[1],res[2],res[3]]

cdef  _mat_mul_list(int a,int b,int c,int d,int x,int y,int z,int w,int res[4]):
    res[0]=a*x+b*z
    res[1]=a*y+b*w
    res[2]=c*x+d*z
    res[3]=c*y+d*w
