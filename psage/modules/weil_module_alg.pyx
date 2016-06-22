# cython: profile=False
# -*- coding: utf-8 -*-
#*****************************************************************************
#  Copyright (C) 2010 Fredrik StrÃ¶mberg <stroemberg@mathematik.tu-darmstadt.de>,
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
#****************************************************************************

r"""
Algorithms and classes for permutations representing subgroups of the modular group, as implemented in 'MySubgroup'.

CLASSES:


AUTHOR:

 - Stephan Ehlen <stephan.j.ehlen@gmail.com>
 - Fredrik Stroemberg <fredrik314@gmail.com>


"""

#include "cysignals/signals.pxi"
from sage.modules.vector_integer_dense import *
from sage.misc.functional import is_even
from sage.rings.arith import kronecker,odd_part,gcd,valuation,is_prime
from sage.rings.integer import Integer
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from psage.matrix.matrix_complex_dense import Matrix_complex_dense
from psage.matrix.matrix_complex_dense cimport Matrix_complex_dense
from sage.all import MatrixSpace,is_odd
from sage.rings.complex_number cimport ComplexNumber
from sage.rings.complex_field import ComplexField
from sage.rings.complex_mpc cimport MPComplexNumber
from sage.rings.complex_mpc import MPComplexField
from sage.rings.real_mpfr cimport RealNumber

cpdef cython_el_index(c, gen_orders):
    cdef long ii, jj = 0
    cdef md = 1
    ii=0
    md=1
    for jj in range(0,len(gen_orders)):
        m=gen_orders[jj]
        #print jj,n,c[jj]
        ii=ii+(c[jj]%int(m))*md
        md=md*m
    return ii

cpdef cython_elt(long ii,gen_orders):
    elt = list()
    cdef long md = 1
    cdef long jj = 0
    cdef long c = 0
    for jj in range(0,len(gen_orders)):
        md = gen_orders[jj]
        c = ii%md
        elt.append(c)
        ii = ii-c
        ii = ii/md
    return elt

cpdef cython_neg_index(long ii, gen_orders):
    cdef long jj=0
    elt = cython_elt(ii,gen_orders)
    for jj in range(0,len(gen_orders)):
        elt[jj]=-elt[jj]
    return cython_el_index(elt,gen_orders)

cpdef cython_neg_indices(long n, gen_orders):
    cdef long jj=0
    l=list()
    for jj in range(0,n):
        l.append(cython_neg_index(jj,gen_orders))
    return l

cpdef cython_xis(int a, int b, int c, int d, W):
    cdef int dc=1
    cdef int argl=0
    JD=W._QM.jordan_decomposition()
    cdef long absD = W._n
    cdef long p = 0
    oddity=W._inv['total oddity']
    oddities=W._inv['oddity']
    pexcesses=W._inv['p-excess']
    sign=W._inv['signature']
    z8=W._z8
    gammafactor=1
    xis=dict()
    for comp in JD:
        p=comp[1 ][0 ]
        xis[p]=1 
    if(a*d <>0 ):
        if(is_even(sign)):
            argl=- 2 *sign
            xis[0 ]=z8**argl
        else:
            dc=kronecker(-a,c)
            argl=-2 *sign
            if(is_even(c)):
                argl=argl+(a+1 )*(odd_part(c)+1 )
            xis[0 ]=z8**argl*dc
    else:
        argl=-sign
        xis[0 ]=z8**argl
        return xis
    if(xis.keys().count(2 )>0 ):
        if(is_odd(c)):
            argl=(c*oddity) % 8 
        else:
            argl=(-(1 +a)*oddity) % 8 
        xis[2]=z8**argl
    for comp in JD:
        [p,n,r,ep]=comp[1 ][:4 ]
        t = None if (len(comp[1 ])==4 ) else comp[1 ][4 ]
        q=p**n
        qc=gcd(q,c)
        qqc=Integer(q/qc)
        ccq=Integer(c/qc)
        nq=valuation(qqc,p)
        gammaf=1 
        dc=1 
        #if(c % q == 0): # This only contributes a trivial factor
        #    continue
        if(p==2 ):
            if(is_even(c)):
                if(c % q<>0 ):
                    odt=self._get_oddity(p,nq,r,ep,t)
                    argl=-a*ccq*odt
                    gammaf=z8**argl
                    dc=kronecker(-a*ccq,qqc**r)
                dc=dc*kronecker(-a,q**r)
                xis[ 2 ]=xis[ 2 ]*gammaf*dc
            else:
                dc=kronecker(c,q**r)
                xis[ 2 ]=xis[ 2 ]*dc
        else:
            if(c%q <>0 ):
                exc=W._get_pexcess(p,nq,r,ep)
                argl=-exc
                gammaf=z8**argl
                dc=kronecker(ccq,qqc**r)
            dc=dc*kronecker(-d,qc**r)
            xis[p]=xis[p]*gammaf*dc
        return xis

cpdef cython_get_gammaN_conjcl_reps(long N):
#    matrix=matrix22()
    if not is_prime(N) and N>2:
        raise NotImplementedError()
    #n=self._QM.order()
    cdef int p=N
    cdef int eps = 0
    F=FiniteField(p)
    four = F(4)
    if p == 3:
        eps = 2
    for ii in range(0,(p-1)/2+1):
        print ii, kronecker(ii,p)
        if ii != 0 and kronecker(ii,p)==-1:
            eps=ii
            print eps, kronecker(ii,p)
            break
    cl=list()
    cl.append([[1,0,0,1],1])
    cl.append([[-1,0,0,-1],1])
    cl.append([[1,1,0,1],(p-1)*(p+1)/2])
    cl.append([[-1,1,0,-1],(p-1)*(p+1)/2])
    cl.append([[1,eps,0,1],(p-1)*(p+1)/2])
    cl.append([[-1,eps,0,-1],(p-1)*(p+1)/2])
    done=list()
    if kronecker(-1,p)==-1:
        cl.append([[0,-1,1,0],p*(p-1)])
    for a in F:
        if a != 0 and a != 1 and a != -1:
            m=[a.lift(),0,0,(a**-1).lift()]
            c=p*(p+1)
            mp=[m,c]
            if not (a**-1) in done:
                cl.append(mp)
                done.append(a)
        x = a**2-four
        if a!= 0 and kronecker(x,p)==-1:
            cl.append([[0,-1,1,a.lift()],p*(p-1)])
    return cl


### Numerical action algorithms

cpdef weil_rep_matrix_mpc(W,a,b,c,d,filter=None,prec=53,verbose=0):
    r"""
    Compute the matrix for the Weil representation
    """
    # Need a WeilModuleElement to compute the matrix
    e=W.basis()[0]
    if filter<>None:
        return action_of_SL2Z_formula_mpc(e,a,b,c,d,filter,prec=prec,verbose=verbose)
    else:
        return action_of_SL2Z_formula_mpc(e,a,b,c,d,prec=prec,verbose=verbose)


## Action by special (simple) elements of SL(2,Z)   
cpdef action_of_T_mpc(W,int b=1,int sign=1,filter=None,prec=53,verbose=0):
    r""" Action by the generator sign*T^pow=[[a,b],[0,d]]
    where a=d=sign
    """
    cdef Matrix_complex_dense r
    cdef int n,ii,level
    n = W._n; level = W._level
    CF = MPComplexField(prec)
    MS = MatrixSpace(CF,n)
    r = Matrix_complex_dense(MS,0)
    cdef MPComplexNumber zl,si
    cdef ComplexNumber z
    cdef RealNumber fac
    z = W._zl.complex_embedding(prec)
    zl = CF(z.real(),z.imag())
    #r = matrix(W._K,W._n)
    if sign==-1:
        z = (W._QM.sigma_invariant()**2).complex_embedding(prec)
        si = CF(z.real(),z.imag())
        #*sigma(Z,T)=1
        #si=W._QM.sigma_invariant()**2 
    else:
        si=CF(1)
    for ii in range(0,n):
        if filter<>None and filter[ii,ii]<>1:
            continue
        if sign==1:
            r[ii,ii] = zl**(b*int(level*W.Q(ii)))
        else:
            r[W._minus_element(ii),ii] = si*zl**(b*int(level*W.Q(ii)))
    fac = CF.base_ring(1)
    return [r,fac]

cpdef action_of_S_mpc(W,filter=None,int sign=1,int mult_by_fact=0,prec=53,verbose=0):
    r"""
    Action by the generator S=[[0,-1],[1,0]]
    """
    cdef Matrix_complex_dense r
    CF = MPComplexField(prec)
    cdef RealNumber fac
    cdef int n,b,ii,level
    n = W._n; level = W._level
    MS = MatrixSpace(CF,n)
    r = Matrix_complex_dense(MS,0)
    cdef MPComplexNumber zl,si
    cdef ComplexNumber z
    z = W._zl.complex_embedding(prec)
    zl = CF(z.real(),z.imag())
    #r = matrix(W._K,W._n)
    if sign==-1:
        z = (W._QM.sigma_invariant()**3).complex_embedding(prec)
        si = CF(z.real(),z.imag())
        if is_odd(W.parent().signature()):
            si= -si  # Here c>0  sigma_Z_A(c,d)
        #si = CFW._K(W._QM.sigma_invariant()**3)
    else:
        z = W._QM.sigma_invariant().complex_embedding(prec)
        si = CF(z.real(),z.imag())
        #si = W._K(W._QM.sigma_invariant())
    for ii in range(0 ,W._n):
        for jj in range(0 ,W._n):
            if filter<>None and filter[ii,jj]<>1:
                continue
            arg = -sign*int(level*W.Bi(ii,jj))
            #arg = -W._level*W._QM.B(W._L[ii],W._L[jj])
            r[ii,jj] = si*zl**arg
    #r = r*
    fac = CF.base_ring()(W._sqn)**-1
    return [r,fac]

cpdef action_of_STn_mpc(W,int pow=1,int sign=1,filter=None,prec=53,verbose=0):
    r""" Action by  S*T^pow
    NOTE: we do not divide by |D|
    """
    ## Have to find a basefield that also contains the sigma invariant
    if pow==0:
        return action_of_S_mpc(W,filter,sign)
    cdef Matrix_complex_dense r
    cdef RealNumber fac
    cdef int n,b,ii,jj,level,argl
    n = W._n; level = W._level
    CF = MPComplexField(prec)
    MS = MatrixSpace(CF,n)
    r = Matrix_complex_dense(MS,0)
    cdef MPComplexNumber zl,si
    cdef ComplexNumber z
    z = W._zl.complex_embedding(prec)
    zl = CF(z.real(),z.imag())
    #r  = matrix(W._K,W._n)
    if sign==-1:
        z = (W._QM.sigma_invariant()**3).complex_embedding(prec)
        si = CF(z.real(),z.imag())
        if is_odd(W.parent().signature()):
            si = -si   #  Here c>0*sigma_Z_A(c,d)
        #si = CFW._K(W._QM.sigma_invariant()**3)
    else:
        z = W._QM.sigma_invariant().complex_embedding(prec)
        si = CF(z.real(),z.imag())
        #si=W._QM.sigma_invariant()
    for ii in range(n):
        # for x in W._L:
        # for y in W._L:
        for jj in range(n):
            argl=int(level*(pow*W.Q(jj)-sign*W.Bi(ii,jj)))
            #ii = W._L.index(x); jj= W._L.index(j)
            if filter<>None and filter[ii,jj]<>1:
                continue
            r[ii,jj] = si*zl**argl
    fac = CF.base_ring()(W._sqn)**-1
    return [r,fac]

cpdef action_of_Z_mpc(W,filter=None,prec=53,verbose=0):
    r""" Action by  Z=-Id
    NOTE: we do not divide by |D|
    """
    ## Have to find a basefield that also contains the sigma invariant
    cdef Matrix_complex_dense r
    CF = MPComplexField(prec)
    cdef RealNumber fac
    cdef int n,ii,level
    n = W._n; level = W._level
    MS = MatrixSpace(CF,n)
    r = Matrix_complex_dense(MS,0)
    cdef MPComplexNumber si
    cdef ComplexNumber z
    #r = matrix(W._K,W._n)
    z = (W._QM.sigma_invariant()**2).complex_embedding(prec)
    si = CF(z.real(),z.imag())
    for ii in range(n):
        if filter<>None and filter[ii,ii]<>1:
            continue
        jj=W._W._neg_index(ii)
        r[ii,jj]=si #CF(1)

    #r = r*W._QM.sigma_invariant()**-2 
    fac = CF.base_ring()(1)
    return [r,fac]

cpdef action_of_Id_mpc(W,filter=None,prec=53,verbose=0):
    r""" Action by  Z=-Id
    NOTE: we do not divide by |D|
    """
    ## Have to find a basefield that also contains the sigma invariant
    #r = matrix(W._K,W._n)
    cdef Matrix_complex_dense r
    CF = MPComplexField(prec)
    cdef RealNumber fac
    cdef int ii
    n = W._n
    MS = MatrixSpace(CF,n)
    r = Matrix_complex_dense(MS,0)
    for ii in range(n):
        if filter<>None and filter[ii,ii]<>1:
            continue
        r[ii,ii]=1 
    #r = r*W._QM.sigma_invariant()**2 
    fac = CF.base_ring()(1)
    return [r,fac]


# Action by Gamma_0(N) through formula
cpdef action_of_Gamma0_mpc(W,a,b,c,d,filter=None,prec=53,verbose=0):
    r"""
    Action by A in Gamma_0(l) 
    where l is the level of the FQM
    INPUT:
       A in SL2Z with A[1,0] == 0 mod l
       act ='r' or 'l' : do we act from left or right'
    filter = |D|x|D| integer matrix with entries 0 or 1
                     where 1 means that we compute this entry
                     of the matrix rho_{Q}(A) 
    """

    cdef Matrix_complex_dense r
    cdef int n,sign,ii,jj,level,argl
    n = W._n; level = W._level
    CF = MPComplexField(prec)
    MS = MatrixSpace(CF,n)
    r = Matrix_complex_dense(MS,0)
    cdef RealNumber fac
    cdef MPComplexNumber zl,si,CI
    cdef ComplexNumber z
    z = W._zl.complex_embedding(prec)
    zl = CF(z.real(),z.imag())
    CI = CF(0,1)
    if c % W._level <>0:
        raise ValueError, "Must be called with Gamma0(l) matrix! not A=" %([a,b,c,d])
    #r = matrix(W._K,W._n)
    for ii in range(n):
        for jj in range(n):
            if(W._L[ii]==d*W._L[jj] and (filter==None or filter[ii,jj]==1 )):
                argl=b*d*int(level*W.Q(W._L[jj]))
                r[ii,jj]=zl**argl
    # Compute the character 
    signature = inv['signature']
    if W._level % 4  == 0:
        test = (signature + kronecker(-1 ,W._n)) % 4
        if is_even(test):
            if test==0:
                power=1 
            elif test==2:
                power=-1 
            if d % 4  == 1:
                chi = 1 
            else:
                chi=CI**power
            chi=chi*kronecker(c,d)
        else:
            if(test==3):
                chi= kronecker(-1 ,d)
            else:
                chi=1 
        chi = chi*kronecker(d,n*2**signature)
    else:
        chi = kronecker(d,n*2**signature)
    r=r*chi
    fac = CF.base_ring()(1)
    return [r,fac]

# Now we want the general action

cpdef action_of_SL2Z_formula_mpc(W,int a,int b,int c,int d,filter=None,int prec=53,verbose=0):
    r"""
    The Action of A in SL2(Z) given by a matrix rho_{Q}(A)
    as given by the formula
    filter = |D|x|D| integer matrix with entries 0 or 1
                     where 1 means that we compute this entry
                     of the matrix rho_{Q}(A) 
    """
    
    cdef Matrix_complex_dense r
    cdef int n,sign,ii,jj,level,arg,ngamma_c,gi,nbm,nb,na
    cdef MPComplexNumber zl,si,zxi
    cdef ComplexNumber z
    cdef RealNumber fac
    sign=1
    if a*d-b*c<>1:
        raise ValueError,"Need matrix in SL(2,Z)!"
    if c==0:
        if b==0:
            if a<0:
                return action_of_Z_mpc(W,filter,prec)
            else:
                return action_of_Id_mpc(W,filter,prec)
        if a<1:
            sign=-1
        else:
            sign=1
        return action_of_T_mpc(W,b,sign,filter,prec)
    cdef int sgnc
    if abs(c)==1 and a==0:
        sgnc=1
        if c<0: sgnc=-1
        return action_of_STn_mpc(W,pow=d*sgnc,sign=sgnc,filter=filter,prec=prec)        
    # These are all known easy cases
    # recall we assumed the formula
    if c<0 or (c==0 and d<0): # change to -A
        a=-a; b=-b; c=-c; d=-d
        #A=SL2Z(matrix(ZZ,2,2,[a,b,c,d]))
        sign=-1 
    else:
        sign=1
    cdef dict xis
    xis=W._get_xis(a,b,c,d)
    cdef ComplexNumber xi
    xi=ComplexField(prec)(1) 
    for q in xis.keys():
        if hasattr(xis[q],"complex_embedding"):
            xi=xi*xis[q].complex_embedding(prec)
        else:
            xi=xi*xis[q]
    norms_c=W._get_all_norm_alpha_cs(c)
    #norms_c_old=W._get_all_norm_alpha_cs_old(c)
    if verbose>0:
        print "xi=",xi
        print "norms=",norms_c
    #print "11"

    n = W._n; level = W._level
    CF = MPComplexField(prec)
    MS = MatrixSpace(CF,n)
    r = Matrix_complex_dense(MS,0)
    z = W._zl.complex_embedding(prec)
    zl = CF(z.real(),z.imag())
    zxi = CF(xi.real(),xi.imag())

    #r = matrix(W._K,W._n)
    if sign==-1:
        z = (W._QM.sigma_invariant()**2).complex_embedding(prec)
        si = CF(z.real(),z.imag())
        if is_odd(W.parent().signature()):
            si = si*sigma_Z_A(c,d)
    else:
        si=CF(1)
    if verbose>0:
        print "si=",si
    for na in range(n):
        for nb in range(n):
            if filter <> None and filter[na,nb]==0:
                continue
            if sign==-1:
                #print type(nb)
                nbm=W._minus_element(nb) #-alpha
            else:
                nbm=nb
            gi=W.lin_comb(na,-d,nbm)
            try:
                ngamma_c=int(norms_c[gi]*level)
            except KeyError:
                #print alpha," not in D^c*"                    
                continue
            arg=a*ngamma_c+b*int(level*W.Bi(gi,nbm))-b*d*int(level*W.Q(nbm))
            if verbose>0 and na==1 and nb==0:
                print "ngamma_c[{0}]={1}".format(gi,ngamma_c)
                print "arg=",               a,"*",ngamma_c,"+",b,"*",W.Bi(gi,nbm),"-",b,"*",d,"*",W.Q(nbm)
                print "arg=",arg
            r[na,nb]=si*xi*zl**(arg)
            if verbose>0 and na==1 and nb==0:
                print "r[",na,nb,"]=",r[na,nb]
    fac = CF.base_ring()(W._get_lenDc(c))
    fac = fac / CF.base_ring()(n)
    fac = fac.sqrt()
    #fac = (fac/n)).sqrt()]
    #print "12"
    return [r,fac]

cpdef sigma_Z_A(c,d):
    r"""
    Return sigma(Z,A) where Z=(-1 0 // 0 -1) and A=(a b // c d)
    """
    if c>0:
        return -1
    if c==0 and d<0:
        return -1
    return 1
