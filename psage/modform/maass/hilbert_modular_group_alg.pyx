# cython: profile=False
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
Cython algorithms for Hilbert modular groups 
Used by routines in development.

"""
# include 'sage/ext/stdsage.pxi'
# include "sage/ext/cdefs.pxi"
# include 'sage/ext/interrupt.pxi'
# include "sage/ext/gmp.pxi"
# #include "sage/rings/mpc.pxi"


#from sage.all import real,NFCusp,copy,RR,CC,RealNumber,ComplexNumber,real,imag,vector,Matrix
#from sage.all import ComplexField,RealField,Infinity,Rational,Integer,QQ
from sage.rings.infinity import Infinity
from sage.matrix.all import Matrix
from sage.modules.all import vector
from sage.rings.real_mpfr import is_RealNumber
from sage.rings.complex_number import is_ComplexNumber
from sage.groups.perm_gps.permgroup_element import is_PermutationGroupElement
from sage.rings.number_field.number_field_element import is_NumberFieldElement
from sage.modular.cusps_nf import NFCusp
from sage.rings.real_mpfr import RR
from sage.libs.mpfr cimport *
cdef mpc_rnd_t rnd
cdef mpfr_rnd_t rnd_re
rnd = MPC_RNDNN
rnd_re = GMP_RNDN
from sage.functions.other import real,imag
from sage.rings.rational import Rational
from sage.rings.rational_field import QQ
from sage.rings.integer import Integer

from sage.rings.complex_field import ComplexField
from sage.rings.complex_mpc cimport MPComplexNumber
from sage.rings.complex_mpc import MPComplexField
from sage.rings.real_mpfr cimport RealNumber,RealField_class
from sage.rings.real_mpfr import RealField
#ctypedef Hn 
#include "gen.pxd"
#try:
#from gen cimport gen
#except:
from sage.libs.pari.gen cimport gen
import cython
#include "hilbert_modular_group_algs.pxd"
cdef extern from "math.h":
    double fabs(double)
    double fmax(double,double)
    int ceil(double)
    int floor(double) 
    double sqrt(double)
    double sin(double)
    double cos(double)
    double log(double)
    double power(double,double)

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

cdef double complex _I = _Complex_I


    

cpdef get_closest_cusp(Hn z,G,int denom_max=3,int verbose=0):
    r""" Locate the closest cusp of G to z.
    """
    cdef int degree,ns,nr,i,j,nsigmamax
    degree=G._degree
    #basis = G._OK.basis()
    cdef double d,ny,cK,delta,delta_min,delta0,x
    cdef double GdeltaK=G.deltaK()

    cdef double *xv, *yv
    cdef int inrhomax,inrhomin,break_rho,break_sigma
    cdef double nrhomax,nrhomin,tmp,tmp1,tmp2,tmpH,tmpL
    cdef double *nrhomax_loc=NULL, *nrhomin_loc=NULL,*s2y2=NULL,*rhoemb=NULL,*semb=NULL
    cdef double *c1=NULL,*c2=NULL,*c3=NULL
    cdef double *nsigmamax_loc
    cdef int np,do_cont
    nsigmamax_loc = <double *>sage_malloc(degree*sizeof(double))
    nrhomax_loc = <double *>sage_malloc(degree*sizeof(double))
    nrhomin_loc = <double *>sage_malloc(degree*sizeof(double))
    s2y2 = <double *>sage_malloc(degree*sizeof(double))
    rhoemb = <double *>sage_malloc(degree*sizeof(double))
    semb = <double *>sage_malloc(degree*sizeof(double))
    xv = <double *>sage_malloc(degree*sizeof(double))
    yv = <double *>sage_malloc(degree*sizeof(double))
    c1 = <double *>sage_malloc(degree*sizeof(double))
    c2 = <double *>sage_malloc(degree*sizeof(double))
    c3 = <double *>sage_malloc(degree*sizeof(double))
    nrhomax = 1.0;  nrhomin = 1.0
    
    if xv ==NULL or yv==NULL or semb==NULL or rhoemb==NULL or s2y2==NULL or nsigmamax_loc==NULL or nrhomax_loc==NULL or nrhomin_loc==NULL or c1==NULL or c2==NULL or c3==NULL:
        raise MemoryError
    for i in range(degree):
        xv[i]=z.x(i); yv[i]=z.y(i)
    #basis = [G._K(1)]; basis.extend(G._K.gens())
    cdef list power_basis
    power_basis = G._K.power_basis()
    if verbose>0:
        print "basis=",power_basis
    cdef int basislen = len(power_basis)
    cdef double **basis_embeddings = NULL
    basis_embeddings = <double**>sage_malloc(basislen*sizeof(double*))
    if basis_embeddings==NULL:
        raise MemoryError    
    for i in range(basislen):
        basis_embeddings[i] = <double*>sage_malloc(degree*sizeof(double*))
        for j in range(degree):
            xi =  power_basis[i].complex_embeddings()[j].real()
            #cbase.append(xi.real())        
            basis_embeddings[i][j]=float(xi)
            if verbose>0:
                print "basis_embeddings[{0}][{1}]={2}".format(i,j,basis_embeddings[i][j])
    ## Get initial estimate by comparing with 
    ny = z.imag_norm()
    cK = G.cK()
    if verbose>0:
        print "cK={0}".format(cK)
        print "N(y)=",ny
    delta0 = 2.0**(3-G._prec)
    cdef int *rho_v, *sigma_v
    rho_v = <int*>sage_malloc(degree*sizeof(int))
    sigma_v = <int*>sage_malloc(degree*sizeof(int))
    cdef double *** cusp_reps=NULL
    cdef int nc = G.ncusps()
    cusp_reps = <double***>sage_malloc(nc*sizeof(double**))
    for i in range(nc):
        cusp_reps[i] = <double**>sage_malloc(2*sizeof(double*))
        cusp_reps[i][0] = <double*>sage_malloc(degree*sizeof(double))
        cusp_reps[i][1] = <double*>sage_malloc(degree*sizeof(double))
        if i==0:
            for j in range(degree):
                cusp_reps[i][0][j]=1.0
                cusp_reps[i][1][j]=0.0
        else:
            for j in range(degree):
                cusp_reps[i][0][j]=G.cusps()[i].numerator().complex_embeddings()[j]
                cusp_reps[i][1][j]=G.cusps()[i].denominator().complex_embeddings()[j]

    # containse both the basis matrix and its inverse
    cdef double*** integer_lattice=NULL
    integer_lattice =  <double***>sage_malloc(2*sizeof(double**))
    if integer_lattice == NULL: raise MemoryError
    integer_lattice[0] =  <double**>sage_malloc(degree*sizeof(double*))
    integer_lattice[1] =  <double**>sage_malloc(degree*sizeof(double*))
    B0 = G.translation_module(0,inv=0)
    B1 = G.translation_module(0,inv=1)
    BA = G.translation_module(0,ret_type='alg')
    for i in range(degree):
        integer_lattice[0][i] =  <double*>sage_malloc(degree*sizeof(double))
        integer_lattice[1][i] =  <double*>sage_malloc(degree*sizeof(double))
        for j in range(degree):
            integer_lattice[0][i][j]=float(B0[i][j])
            integer_lattice[1][i][j]=float(B1[i][j])
    ## We know that if Norm(y)>1 then oo is the closest cusp
    ny = yv[0]
    for i in range(1,degree):
        ny=ny*yv[i]
    if ny>1.0:
        if verbose>0:
            print "oo is the closest cusp since Ny={0}>1".format(ny)
        c = NFCusp(G._K,G._K(1),G._K(0))
        delta_min = ny**-0.5
        return c,delta_min
    get_initial_cusp_distance_c(xv,yv,ny,degree,rho_v,sigma_v,&d,nc,cusp_reps,integer_lattice,denom_max,delta0,verbose)
    cdef int h = G._K.class_number()
    #cdef NumberFieldElement_quadratic rho_min,sigma_min
    rho_min = G._K(0) # NumberFieldElement_quadratic(G._K,0)
    sigma_min = G._K(0) #NumberFieldElement_quadratic(G._K,0)
    rho_min  = BA[0]*rho_v[0] 
    sigma_min= BA[0]*sigma_v[0]
    if verbose>0:
        for i in range(degree):
            print "rho_v[{0}]={1}".format(i,rho_v[i])
            print "sigma_v[{0}]={1}".format(i,sigma_v[i])
    for i in range(1,degree):
        if rho_v[i]<>0:
            rho_min+=BA[i]*rho_v[i]
        if sigma_v[i]<>0:
            sigma_min+=BA[i]*sigma_v[i]
    if verbose>0:
        print "Initial values:"
        print "rho_min=",rho_min
        print "sigma_min=",sigma_min
        print "d=",d
    np = 0
    nsigmamax = ceil(cK**degree*d/ny**0.5)
    #nsigmamax_loc = [] #; yv=copy(z.y); xv=copy(z.x)
    #tmp = cK*ny**(-1/(2*n))
    tmp = float(RR(cK)*RR(d)**(RR(1)/RR(degree)))
    delta_min = d #ny**-0.5
    for i in range(degree):
        nsigmamax_loc[i]=tmp/yv[i]
    # and the bounds for rho (TODO: add ref to notation in paper?)
    for i in range(degree):
        c1[i]=cK**2*yv[i]*d**(RR(2)/RR(degree))
        c2[i]=yv[i]**2
        c3[i]=z.x(i)
    if verbose>0:
        print "sigmamax=",nsigmamax
        for i in range(degree):
            print "sigmamax_loc[{0}]={1}".format(i,nsigmamax_loc[i])
    F = G._K.pari_bnf()
    ## Get positive units
    cdef list pos_norm_units
    pos_norm_units=[G._K(1).list()]
    for u in G._K.units():
        if u.norm()>0:
            pos_norm_units.append(u)
    
    # Gens for K as a vector space over R
    ## Make sure that we don't have a prevsious field in cache
    check_cached_field(F)
    for ns in range(-nsigmamax,nsigmamax):
        if ns==0: continue
        break_sigma=0
        list_of_sigmas =  elements_of_norm(F,ns,degree,basis_embeddings)
        if verbose>1:
            print "ns=",ns
        if verbose>2:
            print "sigmas of norm ns=",list_of_sigmas
        #for sigma in list_of_sigmas:
        for sigmat in list_of_sigmas:            
            for i in range(degree):
                semb[i] = float(sigmat[1][i])
            if verbose>1:
                sigma = G._K(0)
                for i in range(degree):
                    sigma+=G._K(sigmat[0][i])*power_basis[i]
                print "semb="
                for i in range(degree):
                    print "semb[{0}]={1}".format(i,semb[i])
                print "sigmat[0]=",sigmat[0]
                print "basis=",power_basis
                print "semb0=",G._K(sigma).complex_embeddings()
                print "sigma=",sigma
            do_cont=-1
            for i in range(degree):
                if fabs(semb[i])>nsigmamax_loc[i]:
                    do_cont=i
                    break
            if do_cont>=0:
                continue
            elif verbose>1:                    
                print "sigma[{0}]={1} OK comp to :{2}!".format(do_cont,semb[do_cont],nsigmamax_loc[i])            
            ## List the rhos
            nrhomax = 1.0;  nrhomin = 1.0
            for i in range(degree):
                nrhomax_loc[i]=0
                nrhomin_loc[i]=0
            #nrhomax_loc=[]; nrhomin_loc=[]
                
            for i in range(degree):
                if c1[i]-c2[i]*semb[i]<0:
                    s="c2:{0}*{1}={2}".format(i,c2[i],semb[i],c2[i]*semb[i])
                    if verbose>0:
                        print "c1[{0}]={1}".format(i,c1[i])
                        print s
                            
                    raise ArithmeticError,s
                tmp1 = sqrt(c1[i]-c2[i]*semb[i])
                tmp2 = semb[i]*xv[i]
                tmpL = tmp2-tmp1
                tmpH = tmp2+tmp1
                nrhomax_loc[i]=tmpH; nrhomax *= tmpH
                nrhomin_loc[i]=tmpL; nrhomin *= tmpL
            if nrhomax<nrhomin:
                tmp = nrhomin
                nrhomin = nrhomax
                nrhomax = tmp
            if nrhomax<0:
                inrhomax=floor(nrhomax)
            else:
                inrhomax=ceil(nrhomax)
            if nrhomin<0:
                inrhomin=floor(nrhomin)
            else:
                inrhomin=ceil(nrhomin)
            if verbose>1: # and ns==-1:
                print "nrhomin({0})={1}".format(sigma,nrhomin)
                print "nrhomax({0})={1}".format(sigma,nrhomax)
                for i in range(degree):
                    print "nrhomin_loc({0})".format(nrhomin_loc[i])
                    print "nrhomax_loc({0})".format(nrhomax_loc[i])
            for i in range(degree):
                s2y2[i]=float(semb[i]**2*yv[i]**2)
            for nr in range(inrhomin,inrhomax):
                if nr==0:
                    continue
                break_rho=0
                if verbose>2:
                    print "nr=",nr
                list_of_rhos =  elements_of_norm(F,nr,degree,basis_embeddings)
                for rhot in list_of_rhos:
                    for i in range(degree):
                        rhoemb[i]=RR(rhot[1][i])
                    do_cont = -1
                    for i in range(degree):
                        if rhoemb[i]-nrhomax_loc[i]>0.0:
                            if nr==5 and ns==-1 and verbose>1:
                                print "rhoemb too large!"

                                print "rhoemb=",rhoemb[i],type(rhoemb[i])
                                print "nrhomax_loc=",nrhomax_loc[i],type(nrhomax_loc[i])
                                print "diff={0}".format(rhoemb[i]-nrhomax_loc[i])
                            do_cont = i
                            break
                    if do_cont > -1:
                        if verbose>1:
                            print "rho[{0}]={1} too large comp to :{2}!".format(do_cont,rhoemb[do_cont],nrhomax_loc[i])
                        continue
                    for i in range(degree):
                        if rhoemb[i]<nrhomin_loc[i]:
                            do_cont = i

                            if nr==5 and ns==-1 and verbose>0:
                                print "rhoemb too small!"
                                #print "v=",v
                                print "rhoemb=",rhoemb[i],type(rhoemb[i])
                                #print "emb/K=",G._K(rho).complex_embeddings()
                                print "nrhomin_loc=",nrhomin_loc[i],type(nrhomin_loc[i])
                            break
                    if do_cont > -1:
                        if verbose>1:
                            print "rho[{0}]={1} too small comp to :{2}!".format(do_cont,rhoemb[do_cont],nrhomin_loc[i])
                        continue
                    delta = 1                        
                    for i in range(degree):
                        delta*=(-semb[i]*xv[i]+rhoemb[i])**2+s2y2[i]
                    delta=sqrt(delta/ny)
                    np+=1
                    if verbose>0:
                        print "delta=",delta
                    if delta < delta_min-delta0:
                        sigma_min = G._K(0)
                        rho_min = G._K(0)
                        for i in range(degree):
                            sigma_min+=G._K(sigmat[0][i])*power_basis[i]
                            rho_min+=G._K(rhot[0][i])*power_basis[i]
                        delta_min = delta
                        if verbose>0:                                
                            print "Got smaller delta at s={0} r={1}, delta={2}".format(sigma_min,rho_min,delta_min)
                        if delta<GdeltaK:
                            ## WE don't have to search anymore...
                            if verbose>0:
                                print "Tested {0} pairs!".format(np)
                                print "We have distance smaller than delta_K={0}".format(G.deltaK())
                            #c = NFCusp(G._K,rho_min,sigma_min)
                            break_rho=1
                            break
                            #return c,delta_min
                if break_rho==1:
                    break_sigma=1
                    break
            if break_sigma==1:
                break
        if break_sigma==1:
            break        
    #return sigma_min,rho_min,delta_min
    if verbose>0:
        print "Tested {0} pairs!".format(np)
        print "rho_min=",rho_min
        print "sigma_min=",sigma_min
    # Reduce:
    if sigma_min<>0:
        cc = rho_min / sigma_min
        rho_min = cc.numerator()
        sigma_min = cc.denominator()
    ci = G._K.ideal(rho_min,sigma_min)
    if ci<>1:
        g = ci.gens_reduced()[0]
        rho_min = rho_min/g
        sigma_min = sigma_min/g
    c = NFCusp(G._K,G._K(rho_min),G._K(sigma_min))
    sage_free(semb)
    sage_free(rhoemb)
    sage_free(nrhomax_loc)
    sage_free(nrhomin_loc)
    sage_free(s2y2)
    sage_free(rho_v); sage_free(sigma_v)
    sage_free(c1);sage_free(c2);sage_free(c3)
    if cusp_reps<>NULL:
        for j in range(nc):
            if cusp_reps[j]<>NULL:
                if cusp_reps[j][0]<>NULL:
                    sage_free(cusp_reps[j][0])
                if cusp_reps[j][1]<>NULL:
                    sage_free(cusp_reps[j][1])
                sage_free(cusp_reps[j])
        sage_free(cusp_reps)
    if integer_lattice<>NULL:
        if integer_lattice[0]<>NULL:
            for i in range(degree):
                if integer_lattice[0][i]<>NULL:
                    sage_free(integer_lattice[0][i])
            sage_free(integer_lattice[0])
        if integer_lattice[1]<>NULL:
            for i in range(degree):
                if integer_lattice[1][i]<>NULL:
                    sage_free(integer_lattice[1][i])
            sage_free(integer_lattice[1])
        sage_free(integer_lattice)
    ## We do not want to reduce c at the end.
    #for cc in G.cusps():
    #    if c.is_Gamma0_equivalent(cc,G._K.ideal(1)):
    return c,delta_min
    #raise ArithmeticError,"Could not get reduced cusp!"

elements_of_F_with_norm = {}
units_of_F_of_positive_norm = []
current_cached_field=0

cpdef check_cached_field(F):
    r"""
    Call this routine before in any call to elements_of_norm
    if you want to ensure that the field is correct.
    """
    global current_cached_field,elements_of_F_with_norm,units_of_F_of_positive_norm

    if current_cached_field<>0:
        if F<>current_cached_field:
            current_cached_field = F
            elements_of_F_with_norm={}
            units_of_F_of_positive_norm=[]
    else:
        current_cached_field = F

cpdef list_elements_of_norm(K,int norm,int check=0,int verbose=0):
    cdef gen F
    F = K.pari_bnf()
    cdef int degree = K.degree()
    cdef list power_basis
    power_basis = K.power_basis()
    if verbose>0:
        print "basis=",power_basis
    cdef int basislen = len(power_basis)
    cdef double **basis_embeddings = NULL
    basis_embeddings = <double**>sage_malloc(basislen*sizeof(double*))
    if basis_embeddings==NULL:
        raise MemoryError    
    for i in range(basislen):
        basis_embeddings[i] = <double*>sage_malloc(degree*sizeof(double*))
        for j in range(degree):
            xi =  power_basis[i].complex_embeddings()[j].real()
            #cbase.append(xi.real())        
            basis_embeddings[i][j]=float(xi)
    cdef list res
    res = elements_of_norm(F,norm,degree,basis_embeddings,check)
    for i in range(basislen):
        sage_free(basis_embeddings[i])
    sage_free(basis_embeddings)
    return res
    
## Cached version ##
## I think this is better than 
cdef list elements_of_norm(gen F,int n,int degree,double ** basis,int check=0):
    r"""
    Return all elements of norm n, modulo units square
    """
    global elements_of_F_with_norm,units_of_F_of_positive_norm
    cdef int i,j
    cdef int numv
    cdef list v,emb
    cdef double x
    cdef gen a,elts
    if check==1:
        check_cached_field(F)
    if units_of_F_of_positive_norm ==[]:
        units_of_F_of_positive_norm=[F.bnfunit()[0]**0]
        #i = 0
        #print "F=",F
        #print "unit=",F.bnfunit()
        for a in F.bnfunit():
            if a.norm()>0:
                units_of_F_of_positive_norm.append(a)
    print "units=",units_of_F_of_positive_norm
    ## How to check totally positivity in pari?
    if not elements_of_F_with_norm.has_key(n):
        if n==0:
            elts = F.nfbasistoalg(0).list()
        else:            
            elts = F.bnfisintnorm(n)
        elements_of_F_with_norm[n]=[]
        for a in elts:
            for ep in units_of_F_of_positive_norm:
                #print "n=",n
                #print "a=",a
                #print "ep=",ep
                v = F.nfbasistoalg_lift(ep*a).list()
                numv = len(v)
                if numv<degree:
                    for i in range(numv,degree):
                        v.append(0)
                emb = []
                for i in range(degree):
                    x = 0
                    for j in range(numv):
                        print "v[{0}]={1}".format(j,v[j])
                        print "basis[{0}][{1}]={2}".format(j,i,basis[j][i])
                        x+=float(v[j])*basis[j][i]
                    emb.append(float(x))

                elements_of_F_with_norm[n].append((v,emb))
    return elements_of_F_with_norm[n]

cpdef get_initial_cusp_distance(x,y,G,denom_max=3,verbose=0):
    r"""
    Compute initial estimates for the cusp distance for the cusps equivalent to the representatives

    """
    cdef int degree = G._K.degree()
    cdef double *xv=NULL,*yv=NULL
    cdef double d=1.0,ny=1.0
    cdef int i,i0,i1,j0,j1
    cdef int *rho_min=NULL,*sigma_min=NULL
    cdef list rho_v=[],sigma_v=[]
    xv = <double*>sage_malloc(degree*sizeof(double))
    yv = <double*>sage_malloc(degree*sizeof(double))
    rho_min = <int*>sage_malloc(degree*sizeof(int))
    sigma_min = <int*>sage_malloc(degree*sizeof(int))    
    if xv==NULL or yv==NULL or rho_min==NULL or sigma_min==NULL:
        raise MemoryError
    cdef double eps0=2.0**(3-G._prec)
    for i in range(degree):
        xv[i]=x[i]; yv[i]=y[i]
        ny = ny*y[i]

    cdef double *** cusp_reps=NULL
    cdef int nc = G.ncusps()
    cusp_reps = <double***>sage_malloc(nc*sizeof(double**))
    for i in range(nc):
        cusp_reps[i] = <double**>sage_malloc(2*sizeof(double*))
        cusp_reps[i][0] = <double*>sage_malloc(degree*sizeof(double))
        cusp_reps[i][1] = <double*>sage_malloc(degree*sizeof(double))
        if i==0:
            for j in range(degree):
                cusp_reps[i][0][j]=1.0
                cusp_reps[i][1][j]=0.0
        else:
            for j in range(degree):
                cusp_reps[i][0][j]=G.cusps()[i].numerator().complex_embeddings()[j]
                cusp_reps[i][1][j]=G.cusps()[i].denominator().complex_embeddings()[j]
    
    # containse both the basis matrix and its inverse
    cdef double*** integer_lattice=NULL
    integer_lattice =  <double***>sage_malloc(2*sizeof(double**))
    if integer_lattice == NULL: raise MemoryError
    integer_lattice[0] =  <double**>sage_malloc(degree*sizeof(double*))
    integer_lattice[1] =  <double**>sage_malloc(degree*sizeof(double*))
    B0 = G.translation_module(0,inv=0)
    B1 = G.translation_module(0,inv=1)
    BA = G.translation_module(0,ret_type='alg')
    for i in range(degree):
        integer_lattice[0][i] =  <double*>sage_malloc(degree*sizeof(double))
        integer_lattice[1][i] =  <double*>sage_malloc(degree*sizeof(double))
        for j in range(degree):
            integer_lattice[0][i][j]=float(B0[i][j])
            integer_lattice[1][i][j]=float(B1[i][j])
       
        
    get_initial_cusp_distance_c(xv,yv,ny,degree,rho_min,sigma_min,&d,nc,cusp_reps,integer_lattice,denom_max,eps0,verbose)
    sigma_res=G._K(0);    rho_res=G._K(0)
    for i in range(degree):
        sigma_res+=BA[i]*sigma_min[i]
        rho_res+=BA[i]*rho_min[i]
        if verbose>0:
            print "rho_min[{0}]={1}".format(i,rho_min[i])
            print "sigma_min[{0}]={1}".format(i,sigma_min[i])
            print "BA[{0}]={1}".format(i,BA[i])
        #sigma_v.append(sigma_min[i])
        #rho_v.append(rho_min[i])
    sage_free(rho_min); sage_free(sigma_min)
    sage_free(xv); sage_free(yv)
    if cusp_reps<>NULL:
        for j in range(nc):
            if cusp_reps[j]<>NULL:
                if cusp_reps[j][0]<>NULL:
                    sage_free(cusp_reps[j][0])
                if cusp_reps[j][1]<>NULL:
                    sage_free(cusp_reps[j][1])
                sage_free(cusp_reps[j])
        sage_free(cusp_reps)

    if integer_lattice<>NULL:
        if integer_lattice[0]<>NULL:
            for i in range(degree):
                if integer_lattice[0][i]<>NULL:
                    sage_free(integer_lattice[0][i])
            sage_free(integer_lattice[0])
        if integer_lattice[1]<>NULL:
            for i in range(degree):
                if integer_lattice[1][i]<>NULL:
                    sage_free(integer_lattice[1][i])
            sage_free(integer_lattice[1])
        sage_free(integer_lattice)                           
        
    return rho_res,sigma_res,d

cdef get_initial_cusp_distance_c(double* x,double *y,double ny,int degree,int *rho_min,int *sigma_min,double* d,int nc,double*** cusp_reps,double*** integer_lattice, int denom_max=3,double eps0=1e-12,int verbose=0):
    cdef double *rho=NULL, *sigma=NULL,*xv=NULL,*yv=NULL
    cdef double d0,d1,d2
    cdef int i,j

    rho = <double*>sage_malloc(degree*sizeof(double))
    sigma = <double*>sage_malloc(degree*sizeof(double))
    xv = <double*>sage_malloc(degree*sizeof(double))
    yv = <double*>sage_malloc(degree*sizeof(double))
    if xv==NULL or yv==NULL or rho==NULL or sigma==NULL:
        raise MemoryError
    ## Check cusp at infinity
    ## Check cusp at infinity
    if verbose>0:
        for i in range(degree):
            print "x[{0}]={1}".format(i,x[i])
            print "y[{0}]={1}".format(i,y[i])
            print "ny={0}".format(ny)

    ## First check the cusp representatives
    cdef double dmin = 0.0
    for j in range(nc):    
        for i in range(degree):
            rho[i]=cusp_reps[j][0][i]
            sigma[i]=cusp_reps[j][1][i]
        d0 = delta_cusp_c(x,y,rho,sigma,ny,degree)
        if verbose>0:
            print "dist(z,c{0})={1}".format(j,d0)
        if j==0:
            dmin = d0
            for i in range(degree):
                rho_min[i]=int(rho[i])
                sigma_min[i]=int(sigma[i])
                if verbose>0:
                    print "0rho_min[{0}]={1}".format(i,rho_min[i])
                    print "0sigma_min[{0}]={1}".format(i,sigma_min[i])
        else:
            if d0<dmin-eps0:
                dmin = d0
                for i in range(degree):
                    rho_min[i]=int(rho[i])
                    sigma_min[i]=int(sigma[i])
                    if verbose>0:
                        print "1rho_min[{0}]={1}".format(i,rho_min[i])
                        print "1sigma_min[{0}]={1}".format(i,sigma_min[i])

    ## Check cusp at zero
    for i in range(degree):
        rho[i]=0.0; sigma[i]=1.0
    d1 = delta_cusp_c(x,y,rho,sigma,ny,degree)
    if verbose>0:
        print "dist(z,0)=",d1
    if d1<dmin-eps0:
        dmin = d1
        for i in range(degree):
            rho_min[i]=int(rho[i])
            sigma_min[i]=int(sigma[i])
    ## Check cusp at a point 'approximating' x
    ## TODO: Do this in a systematic way. Rational approximation?
    cdef int jmin=0
    cdef int ii,ci
    cdef double dist,dist_min=1
    cdef double* xcoord
    xcoord = <double*>sage_malloc(degree*sizeof(double))
    for i in range(degree):
        xcoord[i]=0.0
        for j in range(degree):
            xcoord[i]+=x[j]*integer_lattice[1][i][j]
        if verbose>0:
            print "xcoord[{0}]={1}".format(i,xcoord[i])
    ## Then we try to go through points in P_1(K) with rational denominators and which are close to the given point
    ## These have coordinates of the form y_i = p/q with |y_i - x_i|<=1/2
    ## i.e. for fixed q:   q*x_i-q/2 <= p <= q*x+q/2

    cdef int p,q
    cdef double qhalf
    ## Specialize to quadratic number fields
    for q in range(1,denom_max):
        pmax = []; pmin=[]
        qhalf = float(q)*0.5
        for i in range(degree):
            pmax.append(xcoord[i]*float(q)+qhalf)
            pmin.append(xcoord[i]*float(q)-qhalf)
        if verbose>0:
            print "pmin=",pmin
            print "pmax=",pmax
        pcoeffs = get_vectors_integer_in_range(degree,pmin,pmax,verbose)
        if verbose>1 and len(pcoeffs)<100:
            print "pcoeffs=",pcoeffs
        for v in pcoeffs:
            for i in range(degree):
                rho[i]=0
                for j in range(degree):
                    rho[i]+=v[j]*integer_lattice[0][i][j]
                sigma[i]=<double>q
            d2 = delta_cusp_c(x,y,rho,sigma,ny,degree)
            if verbose>0:
                rhol=[];vv=[]
                for i in range(degree):
                    rhol.append(rho[i]/float(q))
                    vv.append(QQ(v[i])/QQ(q))
                print "dist(z,{0}:{1})={2}\t {3}\t{4}".format(rhol,q,d2,v,vv)
            if d2<dmin-eps0:
                dmin = d2
                for i in range(degree):
                    rho_min[i]=v[i];
                    if i>0:
                        sigma_min[i]=0 
                    else:
                        sigma_min[0]=q #sigma[i]
            # for ci in range(1,nc):
            #     for i in range(degree):
            #         rho[i]=cusp_reps[ci][0][i]+<double>j
            #         sigma[i]=cusp_reps[ci][1][i]+<double>ii
            #     d2 = delta_cusp_c(x,y,rho,sigma,ny,degree)
            #     if verbose>0:
            #         print "dist(z,c{0}+{1})={2}".format(ci,float(j)/float(ii),d2)

    sage_free(xv)
    sage_free(yv)
    sage_free(rho)
    sage_free(sigma)
    if verbose>0:
        print "dmin=",dmin
        for i in range(degree):
            print "rho_min[{0}]={1}".format(i,rho_min[i])
            print "sigma_min[{0}]={1}".format(i,sigma_min[i])
    d[0] = dmin
    #    print "d1=",d1
    #    print "d2=",d2
    #for i in range(degree):
    #    rho_min[i]=0; sigma_min[i]=0
    #if d2 < d1 and d2 < d0:
    #    rho_min[0]=j; sigma_min[0]=1; d[0]=d2
    #    #return j,0,1,0,d2    
    #elif d1 < d0:
    #    rho_min[0]=0; sigma_min[0]=1; d[0]=d1
    #    #return 0,0,1,0,d1
    #else:
    #    rho_min[0]=1; sigma_min[0]=0; d[0]=d0
    #return 1,0,0,0,d0
    
cpdef delta_cusp(Hn z,ca,cb,int degree,int verbose=0):
    cdef double *x,*y,*rho,*sigma
    cdef int i
    x=NULL;y=NULL;rho=NULL;sigma=NULL
    cdef double delta,ny=1.0
    x = <double*>sage_malloc(degree*sizeof(double))
    y = <double*>sage_malloc(degree*sizeof(double))
    rho = <double*>sage_malloc(degree*sizeof(double))
    sigma = <double*>sage_malloc(degree*sizeof(double))
    if x==NULL or y==NULL or rho==NULL or sigma==NULL:
        raise MemoryError
    for i in range(degree):
        x[i]=z.x(i)
        y[i]=z.y(i)
        rho[i]=float(ca.complex_embeddings()[i])
        sigma[i]=float(cb.complex_embeddings()[i])
        ny = ny*y[i]
    delta = delta_cusp_c(x,y,rho,sigma,ny,degree)
    sage_free(x)
    sage_free(y)
    return delta

@cython.cdivision(True)
cdef double delta_cusp_c(double *x, double *y,double* rho, double* sigma,double ny, int degree):
    r"""
    Compute distance to cusp given by rho/sigma in a number field.
    INPUT:
    
    - `sigma` --  vector of embeddinga of sigma
    - `rho` --  vector of embeddinga of rho
    - `z` --  vector of pairs [x_i,y_i] in the upper half-plane
    """
    cdef double delta = 1.0
    cdef int i
    for i in range(degree):
        delta*=(-sigma[i]*x[i]+rho[i])**2+sigma[i]**2*y[i]**2
    delta=delta/ny
    return sqrt(delta)


cpdef in_SL2OK(M):
    if M.det()<>1:
        return False
    if not M[0,0].is_integral():
        return False
    if not M[0,1].is_integral():
        return False
    if not M[1,0].is_integral():
        return False
    if not M[1,1].is_integral():
        return False
    return True


# cpdef get_bounded_integers_quadratic(z,basis,bounds):
#     r"""
#     Get list of integers in a number field with bounded embeddings            #v = sigma.list()
            #semb = []
            # if len(v)<>G._degree:
            #     for i in range(len(v),degree):
            #         v.append(0)
            # for i in range(degree):
            #     x = 0
            #     for j in range(basislen):
            #         #for b in basis_embeddings:
            #         x+=float(v[j])*basis_embeddings[j][i]
            #     semb[i]=float(x)
            #print "basis=",basis

#     TODO: Figure out how to do tha
#     """
#     ## Get estimates.        #self._elliptic_points[l]=ell_pts

#     cdef int itmax = 10000
#     cdef int i,j
#     cdef list res
#     cdef double b00,b01,b10,b11
#     ### Find the max range of i and j we need to check
    
#     for i in range(itmax):
#         for j in range(itmax):
#             for k in range(2):
#                 a0 = i*basis[0][0]+j*basis[0][1]
#                 if a0 < bounds[0][0] or a0 > bounds[0][1]:
#                     continue
#                 a1 = i*basis[1][0]+j*basis[1][1]
#                 if a1 < bounds[1][0] or a0 > bounds[1][1]:
#                     continue
#                 res.append((i,j))

cdef class  Hn_dble(object):
    r"""
    Class of elements in products of complex upper  half-planes
    with additional ring structure given by:
    - component-wise multiplication and division
    - multiplication and division by elements of a number field of the same degree as the dimension.

    TODO: Inherit from Ring or something? (for speed probably NO!)

    """
    def __cinit__(self,x,y=[],degree=0,prec=53,verbose=0):
        r"""
        Init self from a list or an algebraic point.
        Currently we only work with double (53 bits) precision.
        """
        self._verbose = verbose
        u = []
        if isinstance(x,list) and isinstance(y,list) and y<>[] and x<>[]:
            self._degree = len(x)
            assert len(y)==self._degree
        elif not isinstance(x,list):
            y=[]
            if hasattr(x,'complex_embeddings'):
                for a,b in  x.complex_embeddings(prec):
                    if b>0:
                        u.append(a)
                        y.append(b)
                self._prec=prec
                degree=len(u)
            elif hasattr(x,'_is_Hn'):
                u = x.x(); y = x.y()                
                self._prec=x._prec
            else:
                raise TypeError,"Can not construct point in H^n from %s!" % x

        else:
            ## WE now have to see if we are called with one of the three possibilities:
            ## v=[xlist,ylist]=[[x[0],x[1],...,],[y[0],y[1],...]]
            y = []
            if isinstance(x[0],list):
                if verbose>0:
                    print "x[0] is list"
                if(len(x))<>2:
                    raise ValueError
                u = x[0]; y = x[1]
                degree = len(x)
            ## v=zlist=[z[0],z[1],...,]
            ## Fast cases first
            elif isinstance(x[0],complex):
                if verbose>0:
                    print "x[0] is complex"
                degree = len(x)
                for i in range(degree):
                    u.append(x[i].real)
                    y.append(x[i].imag)
            elif is_ComplexNumber(x[0]):
                if verbose>0:
                    print "x[0] is ComplexNumber"
                degree = len(x)
                for i in range(degree):
                    u.append(x[i].real())
                    y.append(x[i].imag())
            elif isinstance(x[0],float) or is_RealNumber(x[0]):
                degree = len(x)
                for i in range(degree):
                    u.append(x[i])
                    y.append(0)
            elif isinstance(x,list):
                degree = len(x)
                try: 
                    for i in range(degree):
                        u.append(real(x[i]))
                        y.append(imag(x[i]))
                    self._degree = len(x)
                except:
                    raise TypeError,"Can not construct point in H^n from %s!" % x
            else:
                raise TypeError,"Can not construct point in H^n from %s!" % x
        if u<>[]:
            self._xlist = u
        else:
            self._xlist = x
        self._ylist = y
        self._prec = 53  ## Nothing else implemented for now
        if verbose>0:
            print "Set x=",self._xlist
            print "Set y=",self._ylist
            print "set degree=",degree
        self._degree = degree
        assert len(self._xlist)==self._degree and len(self._ylist)==self._degree
        ## norms and indicators
        self._imag_norm = 0.0
        self._real_norm = 0.0
        self._norm = 0.0
        self._imag_norm_set = 0
        self._real_norm_set = 0
        self._norm_set = 0
        self.c_new(self._xlist,self._ylist)
        if verbose>0:
            print "c_new successful!"
        

    def __init__(self,x,y=[],degree=0,prec=53,verbose=0):
        pass
        
    cdef c_new(self,list x,list y):
        self._x = NULL; self._y=NULL
        self._x = <double*>sage_malloc(sizeof(double)*self._degree)
        if self._x==NULL:
            raise MemoryError
        self._y = <double*>sage_malloc(sizeof(double)*self._degree)
        if self._y==NULL:
            raise MemoryError
        cdef int i
        for i in range(self._degree):
            self._x[i] = <double>x[i]
            self._y[i] = <double>y[i]
            if y[i]<0:
                raise ValueError,"Not in H^n^*! y[{0}]={1}".format(i,y[i])
        if self._verbose>0:
            print "allocated x and y!"


    def __richcmp__(self, right, int op):
        res=1
        if op <>2 and op <>3:
            raise NotImplementedError,"Ordering of points in H^n is not implemented!"
        if type(self)<>type(right):
            res=0
        elif right.degree() <> self.degree():
            res=0
        else:
            for i in range(self.degree()):
                if self.x(i)<>right.x(i) or self.y(i)<>right.y(i):
                    res = 0
                    break
            
        if op==3: # != # op=2 is ==
            if res==0:
                return True
            else:
                return False
        if res==0:
            return False
        else:
            return True

    def __copy__(self):
        x = self.real_list()
        y = self.imag_list()
        return Hn_dble(x,y,self._degree,self._prec,self._verbose)
        
    def __repr__(self):
        return str(list(self))

    def __getitem__(self,i):        
        if isinstance(i,(int,Integer)) and i>=0 and i<self._degree:
            return self.z(i)
        else:
            raise IndexError
        
    def __dealloc__(self):
        self._free()
            
    cdef _free(self):
        if self._x<>NULL:
            sage_free(self._x)
        if self._y<>NULL:
            sage_free(self._y)
        self._degree = 0        

    def __add__(self,other):
        r"""
        Add two points in the upper half-plane produce another point
        """
        if not is_Hn_dble(other):
            other = Hn_dble(other)
        x = self.real_list(); y = self.imag_list()
        assert self.degree()==other.degree()
        for i in range(self.degree()):
            x[i]+=other.x(i)
            y[i]+=other.y(i)
        return Hn_dble(x,y)

    def __sub__(self,other):
        r"""
        Substract two points in the upper half-plane may sometimes produce another point.
        """
        if not is_Hn_dble(other):
            other = Hn_dble(other)
        x = self.real_list(); y = self.imag_list()
        assert self.degree()==other.degree()
        for i in range(self.degree()):
            x[i]-=other.x(i)
            y[i]-=other.y(i)
        return Hn_dble(x,y)
                
    def __imag__(self):
        return self.imag()

    def __real__(self):
        return self.real()
    
    def __len__(self):
        return self.degree()


    def addto(self,other):
        #x = self.real_list(); y = self.imag_list()
        if hasattr(other,'complex_embeddings'):
            c = other.complex_embeddings(self._prec)
            assert len(c)==self.degree()
            for i in range(self.degree()):
                self._x[i]+=c[i].real()
                self._y[i]+=c[i].imag()
        elif hasattr(other,'_xlist'):
            assert self.degree()==other.degree()
            for i in range(self.degree()):
                self._x[i]+=other.x(i)
                self._y[i]+=other.y(i)
        elif is_ComplexNumber(other):
            for i in range(self.degree()):
                self._x[i]+=other
        else:
            raise NotImplementedError,"Can not add Hn and %s of type %s!" %(other,type(other))
        for i in range(self.degree()):
            self._xlist[i]=self._x[i]
            self._ylist[i]=self._y[i]
            
        #return Hn(res)

    def prec(self):
        return self._prec

    cpdef imag_norm(self):
        cdef int i
        if self._imag_norm_set==0:
            self._imag_norm=1.0
            for i in range(self._degree):
                self._imag_norm = self._imag_norm*self._y[i]
            self._imag_norm_set=1
        return self._imag_norm


    cpdef real_norm(self):        
        cdef int i
        if self._real_norm_set==0:
            self._real_norm=1.0
            for i in range(self._degree):
                self._real_norm = self._real_norm*self._x[i]
            self._real_norm_set=1
        return self._real_norm

    cpdef imag_list(self):
        return self._ylist


    cpdef real_list(self):
        return self._xlist

    cpdef norm(self):
        cdef int i
        if self._norm_set==0:
            self._norm = <double complex>1.0
            for i in range(self._degree):
                self._norm = self._norm*(self._x[i]+_I*self._y[i])
            self._norm_set=1
        return self._norm 

    cpdef vector_norm(self):
        r"""
        Return the Euclidean norm of self as a vector in C^n
        """
        cdef double res = 1.0
        cdef int i
        for i in range(self._degree):
            res = res+(cabs(self._x[i]+_I*self._y[i]))**2
        return sqrt(res)

    cpdef trace(self):
        cdef double complex res
        res = 1.0
        cdef int i
        for i in range(self._degree):
            res = res+self._x[i]+_I*self._y[i]
        return res

    cpdef real_trace(self):
        cdef double res = 1.0
        cdef int i
        for i in range(self._degree):
            res = res+self._x[i]
        return res

    cpdef imag_trace(self):
        cdef double res = 1.0
        cdef int i
        for i in range(self._degree):
            res = res+self._y[i]
        return res
    

    cpdef degree(self):
        return self._degree
    
    
    cpdef x(self,int i):
        if i>=0 and i<self._degree:
            return self._x[i]
        else:
            raise IndexError

    cpdef y(self,int i):
        if i>=0 and i<self._degree:
            return self._y[i]
        else:
            raise IndexError

    cpdef z(self,int i):
        if i>=0 and i<self._degree:
            return self._x[i]+_I*self._y[i]
        else:
            raise IndexError
        
    cpdef imag(self):        
        return Hn_dble(self._ylist)

    cpdef real(self):
        return Hn_dble(self._xlist)

    def __list__(self):
        res = []
        for i in range(self.degree()):
            res.append(self.z(i))
        return res


    
    def permute(self,s):
        r"""
        The permutation group acts by permuting the components.
        """
        assert is_PermutationGroupElement(s)
        znew = [0 for i in range(self._degree)]
        for i in range(self._degree):
            si = s.list()[i]-1
            znew[si]=self.z(i)
        return Hn_dble(znew)




    
    def acton_by(self,A):
        r"""
        act on self by the matrix A in GL^+(2,K)
        """
        ## Make sure that we act on upper half-planes
        assert A.det().is_totally_positive() 
        res=[]
        prec = self._prec
        if isinstance(A[0,0],Rational):
            a = [A[0,0] for i in range(self._degree)]
            b = [A[0,1] for i in range(self._degree)]
            c = [A[1,0] for i in range(self._degree)]
            d = [A[1,1] for i in range(self._degree)]
        elif is_NumberFieldElement(A[0,0]):
            a = A[0,0].complex_embeddings(prec)
            b = A[0,1].complex_embeddings(prec)
            c = A[1,0].complex_embeddings(prec)
            d = A[1,1].complex_embeddings(prec)
        for i in range(self._degree):
            den = self.z(i)*c[i]+d[i]
            if den<>0:
                w = (self.z(i)*a[i]+b[i])/den
            else:
                w = Infinity
            res.append(w)
        return Hn(res)





    
    def parent(self):
        return None

        
        
    def is_in_upper_half_plane(self):
        r"""
        Returns true if all components are in the upper half-plane.
        """
        for y in self.y:
            if y<0:
                return False
        return True
    
    #def log(self):
    #    res = map(log,list(self))
    #    return Hn(res)
        
    def hyp_dist(self,w,dtype=1):
        r"""
        If self and w are in the upper half-plane then
        return d(self.z(i),z.z(i))
        dtype = 0 => dist = dist1+dist2
        dtype = 1 => dist = max(dist1,dist2)
        TODO: FASTER VERSION...
        """
        if dtype == 1:
            maxd = 0
        if hasattr(w,'_xlist') and w._degree == self._degree:
            if w.z == self.z:
                return 0
            distances=[]
            for i in range(self.degre()):
                ab1 =abs(self.z(i)-w.z(i))
                ab2 =abs(self.z(i)-MPComplexField(self._prec)(w.x(i),-w.y(i)))
                l = log((ab1+ab2)/(ab2-ab1))
                distances.append(l)

            if dtype==0:
                return sum(distances)
            else:
                return max(distances)

    def reflect(self):
        r"""
        z -> -conjugate(z)
        """
        for i in range(self.degree()):
            self._x[i]=-self._x[i]
            self._xlist[i]=-self._xlist[i]

    cpdef diff_abs_norm(self,z):
        r"""
        Return the euclidean norm of the difference of self and z as vectors in C^n: ||self-z||
        """
        cdef double norm = 0.0
        if hasattr(z,'complex_embeddings'):
            c = z.complex_embeddings(self._prec)
            for i in range(self.degree()):
                norm += abs(self.z(i)-c[i])**2
        elif is_Hn_dble(z):
            for i in range(self.degree()):
                norm += abs(self.z(i)-z.z(i))**2
        elif not isinstance(z,list):  ## Assume z is a complex number...
            try:
                for i in range(self.degree()):
                    norm += abs(self.z(i)-z)**2
            except:
                raise NotImplementedError,"Can not substract point in Hn and %s of type %s!" %(z,type(z))
        else:
            raise NotImplementedError,"Can not substract point in Hn and %s of type %s!" %(z,type(z))
        return sqrt(norm)




    ### Below this point some algorithms might be inapproprate leftovers... 
    ### 
    cpdef addto_re(self,x):
        r"""
        Add the real x to self
        """
        if hasattr(x,'complex_embeddings'):
            c = x.complex_embeddings(self._prec)
            for i in range(self.degree()):
                self._x[i]+=c[i]
                self._xlist[i]+=c[i]
        elif is_Hn_dble(x):
            for i in range(self.degree()):
                self._x[i]+=x._x[i]
                self._xlist[i]+=x._xlist[i]
        elif is_RealNumber(x):
            for i in range(self.degree()):
                self._x[i]+=x
                self._xlist[i]+=x
        else:
            raise NotImplementedError,"Can not add Hn and %s of type %s!" %(x,type(x))
        
    def __div__(self,other):
        if hasattr(other,'_is_Hn') and self._degree == other._degree:            
            res = [self.z(i)/other.z(i) for i in xrange(self._degree) ]
        elif hasattr(other,'complex_embeddings') and other.parent().degree()==self._degree:
            w = other.complex_embeddings(self._prec)
            res = [self.z(i)/w[i] for i in range(self.degree()) ]
        elif isinstance(other,list) and len(other)==self.degree():
                res = [self.z(i)/other[i] for i in range(self.degree()) ]
        else:
            raise NotImplementedError,"Can not divide Hn by %s of type %s!" %(other,type(other))
        return Hn_dble(res)


    def __rmul__(self,other):
        if isinstance(other,type(self)):
            assert self._degree == other._degree
            res = [self.z(i)*other.z(i) for i in xrange(self._degree) ]
        elif hasattr(other,'complex_embeddings'):
            w = other.complex_embeddings(self._prec)
            if len(w) == self._degree:
                res = [w[i]*self.z(i) for i in xrange(self._degree) ]
        elif isinstance(other,list) and len(other)==self._degree:
                res = [self.z(i)*other[i] for i in xrange(self._degree) ]
        else:
            res = [self.z(i)*other for i in xrange(self._degree) ]
        return Hn_dble(res)

    def __lmul__(self,other):
        if isinstance(other,type(self)):
            assert self.degree() == other.degree()
            res = [self.z(i)*other.z(i) for i in range(self.degree()) ]
        elif hasattr(other,'complex_embeddings'):
            w = other.complex_embeddings(self._prec)
            assert len(w) == self.degree()
            res = [w[i]*self.z(i) for i in range(self.degree()) ]
        elif isinstance(other,list) and len(other)==self.degree():
            res = [self.z(i)*other[i] for i in range(self.degree()) ]
        else:
            res = [self.z(i)*other for i in range(self.degree()) ]
        return Hn_dble(res)


    def __mul__(self,other):
        if isinstance(other,type(self)):
            assert self.degree() == other._degree
            res = [self.z(i)*other.z(i) for i in range(self.degree()) ]
        elif hasattr(other,'complex_embeddings'):
            w = other.complex_embeddings(self.prec())
            assert len(w) == self.degree()
            res = [w[i]*self.z(i) for i in range(self.degree()) ]
        elif isinstance(other,list) and len(other)==self.degree():
            res = [self.z(i)*other[i] for i in range(self.degree()) ]
        else:
            res = [self.z(i)*other for i in range(self.degree()) ]
        return Hn_dble(res)



cdef class  Hn(object):
    r"""
    Class of elements in products of complex upper  half-planes
    with additional ring structure given by:
    - component-wise multiplication and division
    - multiplication and division by elements of a number field of the same degree as the dimension.

    TODO: Inherit from Ring or something? (for speed probably NO!)

    """
    def __cinit__(self,x,y=[],degree=0,prec=53,verbose=0):
        r"""
        Init self from a list or an algebraic point.
        """
        self._verbose = verbose
        u = []
        if isinstance(x,list) and isinstance(y,list) and y<>[] and x<>[]:
            self._degree = len(x)
            assert len(y)==self._degree
        elif not isinstance(x,list):
            y=[]
            if hasattr(x,'complex_embeddings'):
                for a,b in  x.complex_embeddings(prec):
                    if b>0:
                        u.append(a)
                        y.append(b)
                self._prec=prec
                degree=len(u)
            elif hasattr(x,'_is_Hn'):
                u = x.x(); y = x.y()                
                self._prec=x._prec
            else:
                raise TypeError,"Can not construct point in H^n from %s!" % x

        else:
            ## WE now have to see if we are called with one of the three possibilities:
            ## v=[xlist,ylist]=[[x[0],x[1],...,],[y[0],y[1],...]]
            y = []
            if isinstance(x[0],list):
                if verbose>0:
                    print "x[0] is list"
                if(len(x))<>2:
                    raise ValueError
                u = x[0]; y = x[1]
                degree = len(x)
            ## v=zlist=[z[0],z[1],...,]
            ## Fast cases first
            elif isinstance(x[0],complex):
                if verbose>0:
                    print "x[0] is complex"
                degree = len(x)
                for i in range(degree):
                    u.append(x[i].real)
                    y.append(x[i].imag)
            elif hasattr(x[0],"imag"):
                if verbose>0:
                    print "x[0] is ComplexNumber"
                degree = len(x)
                for i in range(degree):
                    u.append(x[i].real())
                    y.append(x[i].imag())
            elif isinstance(x[0],float) or is_RealNumber(x[0]):
                degree = len(x)
                for i in range(degree):
                    u.append(x[i])
                    y.append(0)
            elif isinstance(x,list):
                degree = len(x)
                try: 
                    for i in range(degree):
                        u.append(real(x[i]))
                        y.append(imag(x[i]))
                    self._degree = len(x)
                except:
                    raise TypeError,"Can not construct point in H^n from %s! of type %s " % (x,type(x[0]))
            else:
                raise TypeError,"Can not construct point in H^n from %s!" % x
        if u<>[]:
            self._xlist = u
        else:
            self._xlist = x
        self._ylist = y
        self._prec = 53  ## Nothing else implemented for now
        if verbose>0:
            print "Set x=",self._xlist
            print "Set y=",self._ylist
            print "set degree=",degree
        self._degree = degree
        assert len(self._xlist)==self._degree and len(self._ylist)==self._degree
        ## norms and indicators
        mpfr_set_ui(self._imag_norm,0,rnd_re)
        mpfr_set_ui(self._real_norm,0,rnd_re)
        mpc_set_ui(self._norm,0,rnd)
        self._imag_norm_set = 0
        self._real_norm_set = 0
        self._norm_set = 0
        self.c_new(self._xlist,self._ylist)
        if verbose>0:
            print "c_new successful!"
        

    def __init__(self,x,y=[],degree=0,prec=53,verbose=0):
        pass
        
    cdef c_new(self,list x,list y):
        cdef RealNumber tmpx
        RF = RealField(self._prec)
        self._x = NULL; self._y=NULL
        self._x = <mpfr_t*>sage_malloc(sizeof(mpfr_t)*self._degree)
        if self._x==NULL:
            raise MemoryError
        self._y = <mpfr_t*>sage_malloc(sizeof(mpfr_t)*self._degree)
        if self._y==NULL:
            raise MemoryError
        self._z = <mpc_t*>sage_malloc(sizeof(mpc_t)*self._degree)
        if self._z==NULL:
            raise MemoryError
        cdef int i
        for i in range(self._degree):
            mpfr_init2(self._x[i],self._prec)
            tmpx = RF(x[i])
            mpfr_set(self._x[i],tmpx.value,rnd_re)
            mpfr_init2(self._y[i],self._prec)
            tmpx = RF(y[i])
            mpfr_set(self._y[i],tmpx.value,rnd_re)            
            #self._y[i] = <double>y[i]
            if tmpx<0:
                raise ValueError,"Not in H^n^*! y[{0}]={1}".format(i,tmpx)
            mpc_init2(self._z[i],self._prec)
            mpc_set_fr_fr(self._z[i],self._x[i],self._y[i],rnd)
        if self._verbose>0:
            print "allocated x and y!"

    def __dealloc__(self):
        self._free()
            
    cdef _free(self):
        cdef int i
        if self._x<>NULL:
            for i in range(self._degree):
                mpfr_clear(self._x[i])
            sage_free(self._x)
        if self._y<>NULL:
            for i in range(self._degree):
                mpfr_clear(self._y[i])
            sage_free(self._y)
        if self._z<>NULL:
            for i in range(self._degree):
                mpc_clear(self._z[i])
            sage_free(self._z)
        self._degree = 0        
        if self._imag_norm_set==1:
            mpfr_clear(self._imag_norm)
        if self._real_norm_set==1:
            mpfr_clear(self._real_norm)
        if self._norm_set==1:
            mpc_clear(self._norm)
            
            
        
    def __richcmp__(self, right, int op):
        res=1
        if op <>2 and op <>3:
            raise NotImplementedError,"Ordering of points in H^n is not implemented!"
        if type(self)<>type(right):
            res=0
        elif right.degree() <> self.degree():
            res=0
        else:
            for i in range(self.degree()):
                if self.x(i)<>right.x(i) or self.y(i)<>right.y(i):
                    res = 0
                    break
            
        if op==3: # != # op=2 is ==
            if res==0:
                return True
            else:
                return False
        if res==0:
            return False
        else:
            return True

    def __copy__(self):
        x = self.real_list()
        y = self.imag_list()
        return Hn(x,y,self._degree,self._prec,self._verbose)
        
    def __repr__(self):
        return str(list(self))

    def __getitem__(self,i):        
        if isinstance(i,(int,Integer)) and i>=0 and i<self._degree:
            return self.z(i)
        else:
            raise IndexError
        
  
    def __add__(self,other):
        r"""
        Add two points in the upper half-plane produce another point
        """
        if not is_Hn(other):
            other = Hn(other)
        x = self.real_list(); y = self.imag_list()
        assert self.degree()==other.degree()
        for i in range(self.degree()):
            x[i]+=other.x(i)
            y[i]+=other.y(i)
        return Hn(x,y)

    def __sub__(self,other):
        r"""
        Substract two points in the upper half-plane may sometimes produce another point.
        """
        if not is_Hn(other):
            other = Hn(other)
        x = self.real_list(); y = self.imag_list()
        assert self.degree()==other.degree()
        for i in range(self.degree()):
            x[i]-=other.x(i)
            y[i]-=other.y(i)
        return Hn(x,y)
                
    def __imag__(self):
        return self.imag()

    def __real__(self):
        return self.real()
    
    def __len__(self):
        return self.degree()


    def addto(self,other):
        #x = self.real_list(); y = self.imag_list()
        cdef MPComplexNumber tmpc
        cdef RealNumber tmpx
        RF = RealNumber(self._prec)
        CF = MPComplexField(self._prec)
        tmpx = RF(0)
        tmpc = CF(0)
        if hasattr(other,'complex_embeddings'):
            c = other.complex_embeddings(self._prec)
            assert len(c)==self.degree()
            for i in range(self.degree()):
                tmpx = RF(c[i].real())
                mpfr_add(self._x[i],self._x[i],tmpx.value,rnd_re)
                tmpx = RF(c[i].imag())
                mpfr_add(self._y[i],self._y[i],tmpx.value,rnd_re)
                tmpc = RF(c[i])
                mpc_add(self._z[i],self._z[i],tmpc.value,rnd)
                #self._x[i]+=c[i].real()
                #self._y[i]+=c[i].imag()
        elif hasattr(other,'_xlist'):
            assert self.degree()==other.degree()
            for i in range(self.degree()):
                tmpx = RF(other.x(i))
                mpfr_add(self._x[i],self._x[i],tmpx.value,rnd_re)
                tmpx = RF(other.y(i))
                mpfr_add(self._y[i],self._y[i],tmpx.value,rnd_re)
                tmpc = RF(other.z(i))
                mpc_add(self._z[i],self._z[i],tmpc.value,rnd)
                #self._y[i]+=other.y(i)
        elif is_ComplexNumber(other):
            for i in range(self.degree()):
                tmpc = CF(other)
                mpfr_add(self._x[i],self._x[i],tmpc.value.re,rnd_re)
                mpfr_add(self._y[i],self._y[i],tmpc.value.im,rnd_re)
                mpc_add(self._z[i],self._z[i],tmpc.value,rnd)
        else:
            raise NotImplementedError,"Can not add Hn and %s of type %s!" %(other,type(other))
        for i in range(self.degree()):
            mpfr_set(tmpx.value,self._x[i],rnd_re)
            self._xlist[i]=tmpx
            mpfr_set(tmpx.value,self._y[i],rnd_re)
            self._ylist[i]=tmpx
            
        #return Hn(res)

    def prec(self):
        return self._prec

    cpdef imag_norm(self):
        cdef int i
        cdef RealNumber res
        res = RealField(self._prec)(0)
        if self._imag_norm_set==0:
            mpfr_init2(self._imag_norm,self._prec)
            mpfr_set_ui(self._imag_norm,1,rnd_re)
            for i in range(self._degree):
                mpfr_mul(self._imag_norm,self._imag_norm,self._y[i],rnd_re)
                #self._imag_norm = self._imag_norm*self._y[i]
            self._imag_norm_set=1
        mpfr_set(res.value,self._imag_norm,rnd_re)
        return res


    cpdef real_norm(self):        
        cdef int i
        cdef RealNumber res
        if self._real_norm_set==0:
            mpfr_init2(self._real_norm,self._prec)
            mpfr_set_ui(self._real_norm,1,rnd_re)
            for i in range(self._degree):
                mpfr_mul(self._real_norm,self._real_norm,self._x[i],rnd_re)
            self._real_norm_set=1
        res = RealField(self._prec)(0)
        mpfr_set(res.value,self._real_norm,rnd_re)
        return res

    cpdef imag_list(self):
        return self._ylist


    cpdef real_list(self):
        return self._xlist

    cpdef norm(self):
        cdef int i
        cdef MPComplexNumber ctmp 
        if self._norm_set==0:
            mpc_init2(self._norm,self._prec)
            mpc_set_si(self._norm,1,rnd)
            for i in range(self._degree):
                mpc_mul(self._norm,self._norm,self._z[i],rnd)
            self._norm_set=1
        ctmp = MPComplexField(self._prec)
        mpc_set(ctmp.value,self._norm,rnd)
        return ctmp

    cpdef vector_norm(self):
        r"""
        Return the Euclidean norm of self as a vector in C^n
        """
        cdef int i
        cdef mpfr_t res_t,tmpx
        cdef RealNumber res
        mpfr_init2(res_t,self._prec)
        mpfr_init2(tmpx,self._prec)
        mpfr_set_si(res_t,0,rnd_re)
        for i in range(self._degree):
            mpc_norm(tmpx,self._z[i],rnd_re)
            mpfr_add(res_t,res_t,tmpx,rnd_re)
        mpfr_sqrt(res_t,res_t,rnd_re)
        res = RealField(self._prec)(0)
        mpfr_set(res.value,res_t,rnd_re)
        mpfr_clear(res_t)
        mpfr_clear(tmpx)
        return res

    cpdef trace(self):
        cdef MPComplexNumber res
        cdef mpc_t res_t
        mpc_init2(res_t,self._prec)
        mpc_set_si(res_t,0,rnd)
        cdef int i
        for i in range(self._degree):
            mpc_add(res_t,res_t,self._z[i],rnd)
            #res = res+self._x[i]+_I*self._y[i]
        res = MPComplexField(self._prec)
        mpc_set(res.value,res_t,rnd)
        return res

    cpdef real_trace(self):
        cdef mpfr_t res_t
        cdef RealNumber res
        mpfr_init2(res_t,self._prec)
        mpfr_set_si(res_t,0,rnd_re)
        cdef int i
        for i in range(self._degree):
            mpfr_add(res_t,res_t,self._x[i],rnd_re)
            #res = res+self._x[i]
        res = RealField(self._prec)(0)
        mpfr_set(res.value,res_t,rnd_re)
        mpfr_clear(res_t)
        return res

    cpdef imag_trace(self):
        cdef mpfr_t res_t
        cdef RealNumber res
        mpfr_init2(res_t,self._prec)
        mpfr_set_si(res_t,0,rnd_re)
        cdef int i
        for i in range(self._degree):
            mpfr_add(res_t,res_t,self._y[i],rnd_re)
            #res = res+self._x[i]
        res = RealField(self._prec)(0)
        mpfr_set(res.value,res_t,rnd_re)
        mpfr_clear(res_t)
        return res

   

    cpdef degree(self):
        return self._degree
    
    
    cpdef x(self,int i):
        cdef RealNumber x
        x = RealField(self._prec)(0)
        if i>=0 and i<self._degree:
            mpfr_set(x.value,self._x[i],rnd_re)
            return x
        else:
            raise IndexError

    cpdef y(self,int i):
        cdef RealNumber y
        y = RealField(self._prec)(0)
        if i>=0 and i<self._degree:
            mpfr_set(y.value,self._y[i],rnd_re)
            return y
        else:
            raise IndexError

    cpdef z(self,int i):
        cdef MPComplexNumber z
        z = MPComplexField(self._prec)(0)
        if i>=0 and i<self._degree:
            mpc_set(z.value,self._z[i],rnd)
            return z
        else:
            raise IndexError
        
    cpdef imag(self):        
        return Hn(self._ylist)

    cpdef real(self):
        return Hn(self._xlist)

    def __list__(self):
        res = []
        for i in range(self.degree()):
            res.append(self.z(i))
        return res


    
    def permute(self,s):
        r"""
        The permutation group acts by permuting the components.
        """
        assert is_PermutationGroupElement(s)
        znew = [0 for i in range(self._degree)]
        for i in range(self._degree):
            si = s.list()[i]-1
            znew[si]=self.z(i)
        return Hn(znew)




    
    def acton_by(self,A):
        r"""
        act on self by the matrix A in GL^+(2,K)
        """
        ## Make sure that we act on upper half-planes
        assert A.det().is_totally_positive() 
        res=[]
        prec = self._prec
        if isinstance(A[0,0],Rational):
            a = [A[0,0] for i in range(self._degree)]
            b = [A[0,1] for i in range(self._degree)]
            c = [A[1,0] for i in range(self._degree)]
            d = [A[1,1] for i in range(self._degree)]
        elif is_NumberFieldElement(A[0,0]):
            a = A[0,0].complex_embeddings(prec)
            b = A[0,1].complex_embeddings(prec)
            c = A[1,0].complex_embeddings(prec)
            d = A[1,1].complex_embeddings(prec)
        for i in range(self._degree):
            den = self.z(i)*c[i]+d[i]
            if den<>0:
                w = (self.z(i)*a[i]+b[i])/den
            else:
                w = Infinity
            res.append(w)
        return Hn(res)





    
    def parent(self):
        return None

        
        
    def is_in_upper_half_plane(self):
        r"""
        Returns true if all components are in the upper half-plane.
        """
        for y in self._ylist:
            if y<0:
                return False
        return True
    
    #def log(self):
    #    res = map(log,list(self))
    #    return Hn(res)
        
    def hyp_dist(self,w,dtype=1):
        r"""
        If self and w are in the upper half-plane then
        return d(self.z(i),z.z(i))
        dtype = 0 => dist = dist1+dist2
        dtype = 1 => dist = max(dist1,dist2)
        TODO: FASTER VERSION...
        """
        if dtype == 1:
            maxd = 0
        if hasattr(w,'_xlist') and w._degree == self._degree:
            if w.z == self.z:
                return 0
            distances=[]
            for i in range(self.degre()):
                ab1 =abs(self.z(i)-w.z(i))
                ab2 =abs(self.z(i)-ComplexField(self._prec)(w.x(i),-w.y(i)))
                l = log((ab1+ab2)/(ab2-ab1))
                distances.append(l)

            if dtype==0:
                return sum(distances)
            else:
                return max(distances)

    def reflect(self):
        r"""
        z -> -conjugate(z)
        """
        for i in range(self.degree()):
            mpfr_neg(self._x[i],self._x[i],rnd_re) #=-self._x[i]
            self._xlist[i]=-self._xlist[i]

    cpdef diff_abs_norm(self,z):
        r"""
        Return the euclidean norm of the difference of self and z as vectors in C^n: ||self-z||
        """
        cdef RealNumber norm
        norm = RealField(self._prec)(0)
        if hasattr(z,'complex_embeddings'):
            c = z.complex_embeddings(self._prec)
            for i in range(self.degree()):
                norm += abs(self.z(i)-c[i])**2
        elif is_Hn(z):
            for i in range(self.degree()):
                norm += abs(self.z(i)-z.z(i))**2
        elif not isinstance(z,list):  ## Assume z is a complex number...
            try:
                for i in range(self.degree()):
                    norm += abs(self.z(i)-z)**2
            except:
                raise NotImplementedError,"Can not substract point in Hn and %s of type %s!" %(z,type(z))
        else:
            raise NotImplementedError,"Can not substract point in Hn and %s of type %s!" %(z,type(z))
        return sqrt(norm)




    ### Below this point some algorithms might be inapproprate leftovers... 
    ### 
    cpdef addto_re(self,x):
        r"""
        Add the real x to self
        """
        cdef RealNumber tmpx
        tmpx = RealField(self._prec)(0)
        #cdef Hn tmphn
        if hasattr(x,'complex_embeddings'):
            c = x.complex_embeddings(self._prec)
            for i in range(self.degree()):
                tmpx = RealField(self._prec)(c[i].real())
                mpfr_add(self._x[i],self._x[i],tmpx.value,rnd_re)
                self._xlist[i]+=c[i].real()
        elif is_Hn(x):
            for i in range(self.degree()):
                tmpx = RealField(self._prec)(x.x(i))
                #self._x[i]+=x._x[i]
                mpfr_add(self._x[i],self._x[i],tmpx.value,rnd_re)
                self._xlist[i]+=tmpx
        elif is_RealNumber(x):
            for i in range(self.degree()):
                tmpx = RealField(self._prec)(x)
                mpfr_add(self._x[i],self._x[i],tmpx.value,rnd_re)
                #self._x[i]+=x
                self._xlist[i]+=x
        else:
            raise NotImplementedError,"Can not add Hn and %s of type %s!" %(x,type(x))
        
    def __div__(self,other):
        if hasattr(other,'_is_Hn') and self._degree == other._degree:            
            res = [self.z(i)/other.z(i) for i in xrange(self._degree) ]
        elif hasattr(other,'complex_embeddings') and other.parent().degree()==self._degree:
            w = other.complex_embeddings(self._prec)
            res = [self.z(i)/w[i] for i in range(self.degree()) ]
        elif isinstance(other,list) and len(other)==self.degree():
                res = [self.z(i)/other[i] for i in range(self.degree()) ]
        else:
            raise NotImplementedError,"Can not divide Hn by %s of type %s!" %(other,type(other))
        return Hn(res)


    def __rmul__(self,other):
        if isinstance(other,type(self)):
            assert self._degree == other._degree
            res = [self.z(i)*other.z(i) for i in xrange(self._degree) ]
        elif hasattr(other,'complex_embeddings'):
            w = other.complex_embeddings(self._prec)
            if len(w) == self._degree:
                res = [w[i]*self.z(i) for i in xrange(self._degree) ]
        elif isinstance(other,list) and len(other)==self._degree:
                res = [self.z(i)*other[i] for i in xrange(self._degree) ]
        else:
            res = [self.z(i)*other for i in xrange(self._degree) ]
        return Hn(res)

    def __lmul__(self,other):
        if isinstance(other,type(self)):
            assert self.degree() == other.degree()
            res = [self.z(i)*other.z(i) for i in range(self.degree()) ]
        elif hasattr(other,'complex_embeddings'):
            w = other.complex_embeddings(self._prec)
            assert len(w) == self.degree()
            res = [w[i]*self.z(i) for i in range(self.degree()) ]
        elif isinstance(other,list) and len(other)==self.degree():
            res = [self.z(i)*other[i] for i in range(self.degree()) ]
        else:
            res = [self.z(i)*other for i in range(self.degree()) ]
        return Hn(res)


    def __mul__(self,other):
        if isinstance(other,type(self)):
            assert self.degree() == other._degree
            res = [self.z(i)*other.z(i) for i in range(self.degree()) ]
        elif hasattr(other,'complex_embeddings'):
            w = other.complex_embeddings(self.prec())
            assert len(w) == self.degree()
            res = [w[i]*self.z(i) for i in range(self.degree()) ]
        elif isinstance(other,list) and len(other)==self.degree():
            res = [self.z(i)*other[i] for i in range(self.degree()) ]
        else:
            res = [self.z(i)*other for i in range(self.degree()) ]
        return Hn(res)
 
        

cpdef is_Hn(z):
    if hasattr(z,"real_list"):
        #if not isinstance(z._imag_norm,float):
        return 1
    return 0

cpdef is_Hn_dble(z):
    if hasattr(z,"real_list"):
        #if isinstance(z._imag_norm,float):
        return 1
    return 0




cpdef cusp_coordinates(G,cuspi,Hn z,int prec=53,int verbose=0):
    r"""
    Gives the coordinate of the CP z withe respcet to the cusp, as in e.g.
    Siegel "Advanced analytic Number theory", p. 171 (147)

    INPUT:
     - z -- CP
     - cuspi  -- integer giving a cusp of self

    """
    if prec==None:
        if is_Hn(z):
            prec = z.prec()
        else:
            prec = G._prec
    if not cuspi in range(G._ncusps): #G.cusps():
        raise ValueError,"Need to call with a cusp representative of the group G! Got:{0}".format(cuspi)
    RF=RealField(prec)
    cdef Hn zz
    if not is_Hn(z):           
        zz = Hn(z)
    else:
        zz = z
    if not cuspi ==0: ## Cheaper than .is_infinity():         
        A = G.cusp_normalizing_map(cuspi,inverse=1)
        zz = z.acton_by(A)
        if verbose>0:
            print "A(zz)=",zz
    cdef double tmp,Y0,ny
    cdef int degree = zz.degree()
    ny = zz.imag_norm()    
    Y0 = sqrt(ny)**-1
    if degree==2:
        #tmp = Y0
        Y = [log(zz.y(0)*Y0)/RF(2)/G._Lambda[0,0]]
    else:
        tmp = ny**(-1.0/<double>degree)  
        Yrhs = vector(RF,degree-1)
        Ymat = Matrix(RF,degree-1,degree-1)
        for i in range(degree-1):
            for j in range(degree-1):
                Ymat[i,j]=G._Lambda[i,j]
            if verbose>0:
                print "ystar[{0}]={1}".format(i,zz.y(i))
            Yrhs[i] = log(zz.y(i)*tmp)/RF(2)
        Y = Ymat.solve_right(Yrhs)

    Xmat = Matrix(RF,degree,degree)
    Xrhs = vector(RF,degree)
    if verbose>2 and degree<>2:
            print "Ymat=",Ymat
            print "Yrhs=",Yrhs
    Xmat = G.translation_module(G.cusps()[cuspi])
    for i in range(degree):
        Xrhs[i] = zz.x(i)
    X = Xmat.solve_right(Xrhs)
    if verbose>2:
        print "Xmat=",Xmat
        print "Xrhs=",Xrhs
    if verbose>0:
        print "Y=",Y
        print "X=",X
    return (Y0,Y,X)



# cpdef cusp_coordinates_arb_cusp(G,cusp,Hn z,int prec=53,int verbose=0):
#     r"""
#     Gives the coordinate of the CP z withe respcet to the cusp, as in e.g.
#     Siegel "Advanced analytic Number theory", p. 171 (147)
#     Note: For arbitrary cusp
#     INPUT:
#      - z -- CP
#      - cusp  -- integer giving a cusp of self

#     """
#     if prec==None:
#         if is_Hn(z):
#             prec = z.prec()
#         else:
#             prec = G._prec
#     if not cusp in range(G._ncusps): #G.cusps():
#         raise ValueError,"Need to call with a cusp representative of the group G!"
#     RF=RealField(prec)
#     cdef Hn zz
#     if not is_Hn(z):           
#         zz = Hn(z)
#     else:
#         zz = z
#     if not cuspi ==0: ## Cheaper than .is_infinity():         
#         A = G.cusp_normalizing_map(cuspi,inverse=1)
#         zz = z.acton_by(A)
#         if verbose>0:
#             print "A(zz)=",zz
#     cdef double tmp,Y0,ny
#     cdef int degree = zz.degree()
#     ny = zz.imag_norm()    
#     Y0 = sqrt(ny)**-1
#     if degree==2:
#         #tmp = Y0
#         Y = [log(zz.y(0)*Y0)/RF(2)/G._Lambda[0,0]]
#     else:
#         tmp = ny**(-1.0/<double>degree)  
#         Yrhs = vector(RF,degree-1)        
#         Ymat = matrix(RF,degree-1)
#         for i in range(degree-1):
#             for j in range(degree-1):
#                 Ymat[i,j]=G._Lambda[i,j]
#             if verbose>0:
#                 print "ystar[{0}]={1}".format(i,zz.y(i))
#             Yrhs[i] = log(zz.y(i)*tmp)/RF(2)
#         Y = Ymat.solve_right(Yrhs)

#     Xmat = matrix(RF,degree)
#     Xrhs = vector(RF,degree)
#     if verbose>2:
#         print "Ymat=",Ymat
#         print "Yrhs=",Yrhs
#     Xmat = G.translation_module(cuspi)
#     for i in range(degree):
#         Xrhs[i] = zz.x(i)
#     X = Xmat.solve_right(Xrhs)
#     if verbose>2:
#         print "Xmat=",Xmat
#         print "Xrhs=",Xrhs
#     if verbose>0:
#         print "Y=",Y
#         print "X=",X
#     return (Y0,Y,X)
    
cpdef get_nearest_integers(list Y):
    r"""
    Compute coordinates k1,...,kn so that -1/2 < Y[j]-kj < 1/2
    
    """
    cdef int i,n,k
    cdef list kv
    n = len(Y)
    kv = []
    for i in range(n):
        if Y[i]>0.5:
            k = -ceil(Y[i]-0.5)
        elif Y[i]<-0.5:
            k = ceil(-0.5-Y[i])         
        else:
            k=0
        kv.append(k)
    return kv

from sage.combinat.combination import IntegerVectors
from sage.combinat.integer_list import IntegerListsLex
cpdef get_vectors_integer_in_range(len,liml,limu,verbose=0): #numv,vectors,verbose=0):
    r"""
    Gives a list of vectors of length len with components integers betweek liml and limu
    """
    cdef list l=[]
    ## Do the easy cases first
    if len==1:
        for i in range(floor(liml[0]),ceil(limu[0])+1):
            l.append([i])
        return l
    elif len==2:
        for i1 in range(floor(liml[0]),ceil(limu[0])+1):
            for i2 in range(floor(liml[1]),ceil(limu[1])+1):
                l.append([i1,i2])
        return l
    ## Brute-force like method. Not efficient
    maxlim = 0; minlim=liml[0]
    for i in range(len):
        if limu[i]>maxlim: maxlim = limu[i]
        if liml[i]<minlim: minlim = liml[i]
    if verbose>0:
        print "maxlim={0}".format(maxlim)
        print "minlim={0}".format(minlim)
    maxlim1 = abs(maxlim)
    ## Produce lists of 0/1 corresponding to signs... 
    signs = []
    for i in range(len):
        signs.extend(IntegerListsLex(i,len,min_part=0,max_part=1).list())
    nsigns = len(signs)
    res=[]
    for sumv in range(maxlim1+1):
        ## Recall that integer vectors are all positive so we have to test
        ## combinations with negative elements as well... 
        IV = IntegerVectors(sumv,len)
        for v in IV:
            for i in range(nsigns):
                vv = []
                is_in = 1
                for j in range(len):
                    if signs[i][j]==1:
                        vv[j]=-v[j]
                    else:
                        vv[j]=v[j]
                    ## check if vv satisfies the conditions
                    if vv[j]<liml[j] or vv[j]>limu[j]:
                        is_in=0
                        break
                if is_in==1:
                    res.append(vv)
    return res
                    
                
    
