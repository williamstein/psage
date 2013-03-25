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

#include "math_inc.pxi"
#from sage.libs.ntl.ntl_ZZX cimport ntl_ZZX
#from sage.libs.ntl.ntl_ZZ cimport ntl_ZZ
#from sage.libs.ntl.ntl_ZZX cimport *
#from sage.libs.ntl.ntl_ZZ cimport ZZ_c 
#from sage.libs.ntl.ntl_ZZ_decl cimport  *
#include "sage/rings/mpc.pxi"
#include "sage/libs/ntl/decl.pxi"
#from sage.libs.ntl.ntl_ZZX cimport ntl_ZZX
#from sage.libs.ntl.ntl_ZZ cimport ntl_ZZ 
#from sage.libs.ntl.ntl_ZZX_decl cimport ZZX_c


from sage.rings.number_field.number_field_element import is_NumberFieldElement
#from sage.rings.number_field.number_field_element cimport NumberFieldElement
from sage.matrix.matrix2 import Matrix

from sage.all import real,NFCusp,copy,RR,CC,RealNumber,ComplexNumber,real,imag,vector,RealField,Infinity,ComplexField,QQ

from psage.modform.hilbert.hn_class cimport Hn
from psage.modform.hilbert.hn_class import is_Hn

#from sage.libs.pari.gen import pari
from sage.libs.pari.gen cimport gen
import cython

#from psage.rings.double_prec_math cimport *
#include "../../rings/double_prec_math.pxi"
from libc.math cimport sqrt,fabs,fmax,ceil,floor,sin,cos,log

cdef test_nf(NumberFieldElement x):
    return x

cpdef cusp_coordinates(G,cuspi,z,int verbose=0):
    r"""
    Gives the coordinate of the CP z withe respcet to the cusp, as in e.g.
    Siegel "Advanced analytic Number theory", p. 171 (147)

    INPUT:
     - z -- CP
     - cuspi  -- integer giving a cusp of self

    """
    print "here!"
    if cuspi<0 or cuspi >= G.ncusps(): #G.cusps():
        raise ValueError,"Need to call with a cusp representative of the group G! Got:{0}".format(cuspi)
    print "here1!"
    cdef Hn zz
    print "here2!"
    if not hasattr(z,"real_norm"):
        zz = Hn(z)
    else:
        print "here3!"
        zz = copy(z)
        print "here4!"
    if verbose>0:
        print "z=",zz
    if not cuspi ==0: ## Cheaper than .is_infinity():         
        A = G.cusp_normalizing_map(cuspi,inverse=1)
        zz = z.acton_by(A)
        if verbose>0:
            print "A(zz)=",zz
    cdef double tmp,Y0,ny
    cdef int degree = zz.degree()
    ny = zz.imag_norm()    
    Y0 = ny**(-1.0/<double>degree)
    #Y0 = sqrt(ny)**-1
    RF=RealField(53)
    if degree==2:
        #tmp = Y0
        Y = [log(zz.y(0)*Y0)/RF(2)/G._Lambda[0,0]]
    else:
        Yrhs = vector(RF,degree-1)        
        Ymat = Matrix(RF,degree-1)
        for i in range(degree-1):
            for j in range(degree-1):
                Ymat[i,j]=G._Lambda[i,j]
            if verbose>0:
                print "ystar[{0}]={1}".format(i,zz.y(i))
            Yrhs[i] = log(zz.y(i)*Y0)/RF(2)
        Y = Ymat.solve_right(Yrhs)

    Xmat = Matrix(RF,degree)
    Xrhs = vector(RF,degree)
    if verbose>2:
        print "Ymat=",Ymat
        print "Yrhs=",Yrhs
    Xmat = G.translation_module(cuspi)
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


cpdef get_closest_cusp(Hn z,G,int denom_max=3,int verbose=0):
    r""" Locate the closest cusp of G to z.

    """
    cdef int degree,ns,nr,i,j,nsigmamax
    degree=G._degree
    #basis = G._OK.basis()
    cdef double d,ny,cK,delta,delta_min,delta0,x
    cdef double GdeltaK=G.deltaK()
    cdef gen F
    F = G._K.pari_bnf()
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

    cdef list power_basis
    power_basis = G._K.power_basis()
    if verbose>0:
        print "basis=",power_basis
    cdef int basislen = len(power_basis)
    cdef double **basis_embeddings = NULL
    basis_embeddings = <double**>sage_malloc(basislen*sizeof(double*))
    cdef list basis_list,tmp_list
    if basis_embeddings==NULL:
        raise MemoryError    
    basis_list=[]
    for i in range(basislen):
        basis_embeddings[i] = <double*>sage_malloc(degree*sizeof(double*))
        tmp_list=[]
        for j in range(degree):
            xi =  power_basis[i].complex_embeddings()[j].real()
            #cbase.append(xi.real())        
            basis_embeddings[i][j]=float(xi)
            tmp_list.append(float(xi))
            if verbose>0:
                print "basis_embeddings[{0}][{1}]={2}".format(i,j,basis_embeddings[i][j])
        basis_list.append(tmp_list)
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
    #cdef double *** cusp_reps_coord=NULL
    cdef double disc
    disc = <double> float(G.number_field().discriminant())
    cdef int nc = G.ncusps()
    cusp_reps = <double***>sage_malloc(nc*sizeof(double**))
    #cusp_reps_coord = <double***>sage_malloc(nc*sizeof(double**))
    B0 = G.translation_module(0)
    BDA = G.translation_module(0,ret_type='alg')

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
#    get_initial_cusp_distance_c(F,xv,yv,ny,degree,rho_v,sigma_v,&d,nc,cusp_reps,cusp_reps_coord,integer_lattice,disc,denom_max,delta0,verbose)
    get_initial_cusp_distance_c(xv,yv,ny,degree,rho_v,sigma_v,&d,nc,cusp_reps,integer_lattice,disc,denom_max,delta0,verbose)
    cdef int h = G._K.class_number()
    #cdef NumberFieldElement_quadratic rho_min,sigma_min
    rho_min = G._K(0) # NumberFieldElement_quadratic(G._K,0)
    sigma_min = G._K(0) #NumberFieldElement_quadratic(G._K,0)
    rho_min  = BA[0]*rho_v[0] 
    sigma_min= BDA[0]*sigma_v[0]
    if verbose>0:
        for i in range(degree):
            print "After get_initial_cusp_dist:"
            print "rho_v[{0}]={1}".format(i,rho_v[i])
            print "sigma_v[{0}]={1}".format(i,sigma_v[i])
    
    for i in range(1,degree):
        if rho_v[i]<>0:
            rho_min+=BA[i]*rho_v[i]
        if sigma_v[i]<>0:
            sigma_min+=BDA[i]*sigma_v[i]

                
    if verbose>0:
        print "Initial values:"
        print "rho_min=",rho_min
        print "sigma_min=",sigma_min
        print "d=",d
    #verbose=0
    np = 0
    nsigmamax = <int>ceil(cK**degree*d/ny**0.5)
    #nsigmamax_loc = [] #; yv=copy(z.y); xv=copy(z.x)
    #tmp = cK*ny**(-1/(2*n))
    #if ida<>None:
    #    cK = cK*float(ida.norm())
    tmp = float(RR(cK)*RR(d)**(RR(1)/RR(degree)))
    delta_min = d #ny**-0.5
    for i in range(degree):
        nsigmamax_loc[i]=tmp/sqrt(yv[i])
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
    # Gens for K as a vector space over R
    ## Make sure that we don't have a prevsious field in cache
    check_cached_field(F)
    nsigmamax = nsigmamax
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
            if verbose>2:
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
                for i in range(degree):
                    print "norm_max_loc[",i,"]=",nsigmamax_loc[i]
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
                    s="c1[{0}]={1} \n".format(i,c1[i])
                    s+="c2[{0}]:{1}*{2}={3}".format(i,c2[i],semb[i],c2[i]*semb[i])
                    if verbose>0:
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
                inrhomax=<int>floor(nrhomax)
            else:
                inrhomax=<int>ceil(nrhomax)
            if nrhomin<0:
                inrhomin=<int>floor(nrhomin)
            else:
                inrhomin=<int>ceil(nrhomin)
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
                        rhoemb[i]=<double>float(rhot[1][i])
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
                        if verbose>2:
                            print "rho[{0}]={1} too large comp to :{2}!".format(do_cont,rhoemb[do_cont],nrhomax_loc[i])
                        continue
                    for i in range(degree):
                        if rhoemb[i]<nrhomin_loc[i]:
                            do_cont = i
                            if verbose>2:
                                print "rhoemb too small!"
                            if nr==5 and ns==-1 and verbose>1:
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
                    if verbose>1:
                        print "delta=",delta
                    if delta < delta_min-delta0:
                        sigma_min = G._K(0)
                        rho_min = G._K(0)
                        for i in range(degree):
                            sigma_min+=G._K(sigmat[0][i])*power_basis[i]
                            rho_min+=G._K(rhot[0][i])*power_basis[i]
                        ### We only want 'relatively prime' elements
                        ctest = rho_min/sigma_min
                        c_num = G._K(ctest.numerator())
                        c_den = G._K(ctest.denominator())
                        if verbose>0:
                            print "Possible min at c=",ctest
                            print "numerator=",c_num,type(c_num)
                            print "denominator=",c_den,type(c_den)
                        delta_test = delta_cusp(z,c_num,c_den,1.0,degree)
                        if delta_test < delta_min-delta0:
                            delta_min = delta_test
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
        print "delta_min=",delta_min
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
    ### Check to see if we have a cusp representative
    for i in range(degree):
        if c==G.cusps()[i]:
            c = G.cusps()[i]
            break        
        
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
    ## Reduce c at the end.    
    return c,delta_min

elements_of_F_with_norm = {}
current_cached_field=0

def check_cached_field(F):
    r"""
    Call this routine before in any call to elements_of_norm
    if you want to ensure that the field is correct.
    """
    global current_cached_field
    if current_cached_field<>0:
        if F<>current_cached_field:
            current_cached_field = F
            elements_of_norm={}
            elements_of_norm_in_ideal={}
    else:
        current_cached_field = F
## Cached version ##
## I think this is better than 
cdef list elements_of_norm(gen F,int n,int degree,double ** basis,int check=0):
    global elements_of_F_with_norm
    cdef int i,j
    cdef int numv
    cdef list v,emb
    cdef double x
    cdef gen a,elts
    if check==1:
        check_cached_field(F)
    if not elements_of_F_with_norm.has_key(n):
        elts = F.bnfisintnorm(n)
        elements_of_F_with_norm[n]=[]

        for a in elts:
            v = a.list()
            numv = len(v)
            if numv<degree:
                for i in range(numv,degree):
                    v.append(0)
            emb = []
            for i in range(degree):
                x = 0
                for j in range(numv):
                    x+=float(v[j])*basis[j][i]
                emb.append(float(x))

            elements_of_F_with_norm[n].append((v,emb))
    return elements_of_F_with_norm[n]

elements_of_F_with_norm_in_ideal = {}
current_cached_field=0


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
    cdef double disc
    cdef gen F
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
    disc = <double> float(G.number_field().discriminant())
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
        
    get_initial_cusp_distance_c(xv,yv,ny,degree,rho_min,sigma_min,&d,nc,cusp_reps,integer_lattice,disc,denom_max,eps0,verbose)
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


cdef get_initial_cusp_distance_c(double* x,double *y,double ny,int degree,int *rho_min,int *sigma_min,double* d,int nc,double*** cusp_reps,double*** integer_lattice, double disc, int denom_max=3,double eps0=1e-12,int verbose=0):
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
    #vr = pari.vector(degree)
    #vs = pari.vector(degree)
    cdef gen idp
    for j in range(nc):

        for i in range(degree):
            rho[i]=cusp_reps[j][0][i]
            sigma[i]=cusp_reps[j][1][i]
        norm = 1.0 #get_ideal_norm(F,vr,vs)
        d0 = delta_cusp_c(x,y,rho,sigma,ny,norm,degree)
        if verbose>0:
            print "dist(z,c{0})={1}".format(j,d0)
        if j==0:
            dmin = d0
            for i in range(degree):
                rho_min[i]=<int>round(rho[i])
                sigma_min[i]=<int>round(sigma[i])
                if verbose>0:
                    print "0rho_min[{0}]={1}".format(i,rho_min[i])
                    print "0sigma_min[{0}]={1}".format(i,sigma_min[i])
        else:
            if d0<dmin-eps0:
                dmin = d0
                for i in range(degree):
                    rho_min[i]=<int>round(rho[i])
                    sigma_min[i]=<int>round(sigma[i])
                    if verbose>0:
                        print "1rho_min[{0}]={1}".format(i,rho_min[i])
                        print "1sigma_min[{0}]={1}".format(i,sigma_min[i])
    ## Check cusp at zero
    for i in range(degree):
        rho[i]=0.0; sigma[i]=1.0
        
    d1 = delta_cusp_c(x,y,rho,sigma,ny,1.0,degree)
    if verbose>0:
        print "dist(z,0)=",d1
    if d1<dmin-eps0:
        dmin = d1
        for i in range(degree):
            rho_min[i]=<int>round(rho[i])
            sigma_min[i]=<int>round(sigma[i])
    ## Check cusp at a point 'approximating' x
    ## TODO: Do this in a systematic way. Rational approximation?
    cdef int jmin=0
    cdef int ii,ci
    cdef double dist,dist_min=1
    cdef double* xcoord=NULL
    xcoord = <double*>sage_malloc(degree*sizeof(double))
    if xcoord == NULL:
        raise MemoryError
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
            norm = 1.0 #get_ideal_norm(F,vr,vs)
            d2 = delta_cusp_c(x,y,rho,sigma,ny,norm,degree)
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

    sage_free(xv)
    sage_free(yv)
    sage_free(rho)
    sage_free(sigma)
    sage_free(xcoord)
    if verbose>0:
        print "dmin=",dmin
        for i in range(degree):
            print "rho_min[{0}]={1}".format(i,rho_min[i])
            print "sigma_min[{0}]={1}".format(i,sigma_min[i])
    d[0] = dmin
 
cpdef delta_cusp(Hn z,ca,cb,double norm, int degree,int verbose=0):
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
    
    delta = delta_cusp_c(x,y,rho,sigma,ny,norm,degree,verbose)
    sage_free(x)
    sage_free(y)
    sage_free(rho)
    sage_free(sigma)
    return delta

@cython.cdivision(True)
cdef double delta_cusp_c(double *x, double *y,double* rho, double* sigma,double ny,
                         double norm,int degree,int verbose=0):
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
    return sqrt(delta/ny/norm) 


from sage.combinat.combination import IntegerVectors
from sage.combinat.integer_list import IntegerListsLex
cpdef get_vectors_integer_in_range(len,liml,limu,verbose=0): #numv,vectors,verbose=0):
    r"""
    Gives a list of vectors of length len with components integers betweek liml and limu
    """
    cdef list l=[]
    ## Do the easy cases first
    if len==1:
        for i in range(<int>floor(liml[0]),<int>ceil(limu[0])+1):
            l.append([i])
        return l
    elif len==2:
        for i1 in range(<int>floor(liml[0]),<int>ceil(limu[0])+1):
            for i2 in range(<int>floor(liml[1]),<int>ceil(limu[1])+1):
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

###
### Helper functions
### 

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

cpdef get_nearest_integers(list Y,int sgn=1):
    r"""
    Return a list of nearest integers to the entries of the list Y.
    I.e.   kv=[k1,...,kn] where -1/2 < Y[j]-kj < 1/2
    INPUT:

    - 'Y' list
    - 'sgn' -- integer
    
    OUTPUT:

    - 'M' -- list

    EXAMPLE:

    
    
    """
    cdef int i,n,k
    cdef list res=[]
    cdef double *Yv=NULL
    cdef int* NI=NULL
    n = len(Y)
    Yv = <double*> sage_malloc(sizeof(double)*n)
    NI = <int*> sage_malloc(sizeof(int)*n)
    if Yv==NULL: raise MemoryError
    for i in range(n):
        Yv[i]=<double>Y[i]
    get_nearest_integers_dp(n,Yv,NI)
    for i in range(n):
        res.append(NI[i])
    sage_free(Yv)
    sage_free(NI)
    return res


cdef get_nearest_integers_dp(int len,double *Y,int* NI):
    r"""
    Return a list of nearest integers to the entries of the list Y.
    I.e.   kv=[k1,...,kn] where -1/2 < Y[j]-kj < 1/2

    INPUT:

    - 'Y' list

    OUTPUTL

    - 'M' -- list

    EXAMPLE:

    
    
    """
    if NI==NULL: raise MemoryError,"Need to allocate NI!"
    cdef int i
    for i in range(len):
        if Y[i]>0.5:
            NI[i] = <int>ceil(Y[i]-0.5)
        elif Y[i]<-0.5:
            NI[i] = <int>(-ceil(-0.5-Y[i]))
        else:
            NI[i]=0



