# cython: profile = True, boundscheck = False, cdivision = False, cdivision_warnings=False, wraparound = False
# -*- coding: utf-8 -*-
#*****************************************************************************
#  Copyright (C) 2013 Stephan Ehlen, Nils Skoruppa
#
#  Distributed under the terms of the GNU General Public License (GPLv2)
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
 Implements an algorithm to compute the invariants
 of the group ring attached to a finite quadratic module
 under the Weil representation.

 AUTHORS::

  - Nils Skoruppa <nils.skoruppa@gmail.com>
  - Stephan Ehlen <stephan.j.ehlen@gmail.com>

 EXAMPLES::

  sage: A = FiniteQuadraticModule('5^2')
  sage: %time X=cython_invariants(A)
  QQbar
  0.036002: table
  [5, 5]
  Gram matrix conversion: 0.0
  0.000000: +- reps
  ni = 5
  n = 13
  0.000000: sorting
  0.028001: init of H
  0.616040: kernel
  0.072004: span
  CPU times: user 1.16 s, sys: 0.04 s, total: 1.20 s
  Wall time: 1.25 s
  sage: X
  ([(0, 0, 2), (7, 0, 1), (8, 0, 1), (11, 0, 1), (14, 0, 1)],
  Vector space of degree 5 and dimension 2 over Algebraic Field
  Basis matrix:
  [ 1.000000000000000?                   0  1.000000000000000?  1.000000000000000?             0.?e-16]
  [                  0  1.000000000000000? -1.000000000000000? -1.000000000000000?  1.000000000000000?])
"""

#from psage.modules.finite_quadratic_module import *

from cysignals.memory cimport sig_free,sig_malloc
from cysignals.signals cimport sig_on,sig_off

from sage.modules.free_module import span
from sage.matrix.constructor import Matrix
from sage.rings.qqbar import QQbar
from sage.all import copy, exp, Integer, pi, I, walltime, CyclotomicField, ZZ, QQ, is_prime_power, \
    kronecker, vector, CC, GF, next_prime, lcm, sqrt, cached_function, MatrixSpace#, sig_malloc, sig_free, ZZ
from sage.rings.number_field.number_field import NumberField_cyclotomic
from cython.parallel import prange
from sage.matrix.matrix_modn_dense_double cimport Matrix_modn_dense_double
from sage.matrix.matrix_modn_dense_float cimport Matrix_modn_dense_float
import numpy as np


cdef long _el_index(long * c, long *ed, int r) nogil:
    r"""
    Returns the canonical index of the element specified by the list c
    in a finite abelian group with elementary divisors specified by `ed`.

    NOTES::
        Note that the index does not only depend on the divisors but also
        on the ordering of these. As long as an ordering is used that is guaranteed to be stable,
        everything will work.
    """
    cdef long ii, jj = 0
    cdef long md = 1
    cdef long m = 0
    ii=0
    md=1
    for jj in range(0,r):
        m=ed[jj]
        #print jj,n,c[jj]
        ii = ii + (c[jj]%m)*md
        md=md*m
    return ii

cdef long* _elt(long ii, long *ed, int r) nogil:
    r"""
    Returns the element corresponding to the canonical index `ii`
    in a finite abelian group with elementary divisors specified by `ed`.
    """
    cdef long* eltl = NULL
    eltl = <long*> sig_malloc(sizeof(long)*r)
    cdef long md = 1
    cdef long jj = 0
    cdef long c = 0
    for jj in range(0,r):
        md = ed[jj]
        c = ii%md
        eltl[jj] = c
        ii = ii - c
        ii = ii/md
    return eltl

cdef long _neg_index(long ii, long *ed, int r) nogil:
    r"""
    Returns the negative of the element corresponding to the canonical index `ii`
    in a finite abelian group with elementary divisors specified by `ed`.
    """
    cdef long jj=0
    cdef long* eltl = _elt(ii, ed, r)
    for jj in range(0,r):
        eltl[jj]=-eltl[jj]
    return _el_index(eltl, ed, r)

cpdef norm_cmp(x, y):
    if x[1] < y[1]:
        return int(-1)
    else:
        return int(1)
    
cdef long B(long i, long j, long **JJ, long *ed, int r) nogil:
    r"""
    Returns the bilinear form of the elements with canonical index `i` and `j`
    with respect to the Gram matrix JJ and elementary divisors `ed`.
    """
    #sig_on()
    cdef long* ll = _elt(j, ed, r)
    cdef long* kk = _elt(i, ed, r)
    cdef long ii, jj = 0
    cdef long res = 0
    for ii in range(r):
        for jj in range(r):
            res = res + ll[ii]*kk[jj]*JJ[ii][jj]
    #sig_off()
    return res

cdef long* Bl(long i, long **JJ, long *ed, int r) nogil:
    cdef long* ll = _elt(i, ed, r)
    cdef long ii, jj = 0
    cdef long* kk = NULL
    kk = <long*> sig_malloc(sizeof(long)*r)
    for ii in range(r):
        kk[ii] = 0
        for jj in range(r):
            kk[ii] = kk[ii] + ll[jj]*JJ[jj][ii]
    return kk

cdef long BBl(long j, long *Bi, long *ed, int r) nogil:
    cdef long ii, res = 0
    cdef long* ll = _elt(j, ed, r)
    for ii in range(r):
        res = res + Bi[ii]*ll[ii]
    return res

cpdef cython_invariants_dim(FQM, use_reduction = True, proof = True, debug=1):
    if debug > 0:
        print "Computing invariants dimension of {}".format(FQM)
    if FQM.signature() % 2 != 0:
        return 0
    dim = 0
    if not is_prime_power(FQM.level()):
        for p,n in FQM.level().factor():
            N = None
            for j in range(1,n+1):
                J = FQM.jordan_decomposition()
                C = J.constituent(p**j)[0]
                if N == None and C.level() != 1:
                    N = C
                elif C.level() != 1:
                    N = N + C
            if dim == 0:
                dim = cython_invariants_dim(N, use_reduction, proof, debug)
            else:
                dim = dim*cython_invariants_dim(N, use_reduction ,proof, debug)
            if dim == 0:
                return 0
    else:
        Sp = cython_invariants(FQM, use_reduction, proof, debug)[1]
        dim = dim + Sp.dimension()
    return dim

cpdef cython_invariants_matrices(FQM, K = QQbar, proof = True, debug=0, return_H = False):
    cdef long i, j = 0
    cdef long p = 0
    cdef int l = long(FQM.level())
    cdef long n = long(FQM.order())
    cdef long ii,jj = 0
    if debug > 1: print "l = ", l
    
    try:
        s = FQM.sigma_invariant()
        s2 = Integer( s**2)
    except:
        return span( [], K)

    if debug > 0: t = walltime()
    cdef long[:] table = np.ndarray(l, dtype=long)
    cdef list table0 = list()

    cdef long q = K.characteristic()
    if q > 0:
        if debug > 0: print 'positive characteristic: ', q
        if 1 != q % l:
            raise ValueError( '%d must be = 1 modulo %d' %(q, l))
        if not q % 4 == 1:
                raise ValueError( '%d must be = 1 modulo 4.' %(q))
        if not is_prime_power(l):
            raise NotImplementedError('This function can only be called with p-modules.')

        pr = K.primitive_element()
        # now we choose I, sqrt(|FQM|) compatible with the choice of a primitive element
        I = pr**((q-1)/4)
        if not s2 == 1:
            if FQM.signature() == 2:
                s = -I
            else:
                s = I
        z = pr**((q-1)/l)
        A = Integer(FQM.order())
        #print A
        if A.is_square():
            w = K(sqrt(A))
        else:
            AA = FQM.order()/FQM.order().squarefree_part()
            AA = sqrt(AA)
            AA = K(Integer(AA))
            P = A.prime_factors()[0]
            if P > 2:
                PP = kronecker(-1,P)*P
                zz = z**(l/P)
                w = sum([kronecker(PP,ii)*zz**ii for ii in range(P)])
                eps = 1 if PP > 0 else -I
                w = eps*w
                w = w*AA
            else:
                zz = pr**((q-1)/8)
                w = (zz+zz**(-1))
                w = w*AA
        if debug > 0:
            print "w = {}".format(w)
        for ii in range(l):
            zt = z**ii
            if debug > 0:
                print "zt = {}**{} = {}".format(z, ii, zt)
            table[ii] = long(K(s)*K(zt + s2 * zt**-1)/K(w))
        if proof and q > 0:
            if debug > 0: tt = walltime()
            if debug > 0: print "proof"
            L = QQbar
            zl = L.zeta(l)
            if s.parent() != ZZ:
                z8 = L.zeta(8)
                sl = z8**(-FQM.signature())
            else:
                sl = s
            # print sl
            wl = L(FQM.order()).sqrt()
            table0 = [sl*(zl**p)/wl for p in range(l)]
            if debug > 0: print table0
            if debug > 0: print "table0: {0}".format(walltime(tt))
    else:
        try:
            w = K(FQM.order()).sqrt()
        except:
            raise RuntimeError("K = {0} does not contain a square-root of |FQM| = {1}".format(K,FQM.order()))
    if debug > 0: print q,w

    if 0 == q:
        if isinstance(K,NumberField_cyclotomic):
            if debug > 0: print 'cyclotomic'
            o = K.gen().multiplicative_order()
            K = CyclotomicField(o, embedding=CC(QQbar.zeta(o)))
            z = K.gen()
            if not Integer(l).divides(o):
                raise ValueError("K has to contain {0}th root of unity.".format(l))
            z = z**(Integer(o)/Integer(l))
            # ensure we have the correct sqrt of FQM.order()
            w = K(FQM.order()).sqrt()
            if w.complex_embedding().real().sign() < 0:
                w = -w
            if 1 == s2: 
                table = [s*((z**p) + (z**p).conjugate())/w for p in range(l)]
            else:
                table = [s*((z**p) - (z**p).conjugate())/w for p in range(l)]
        else:
            if K == QQbar:
                if debug > 0: print 'QQbar'
                z = K.zeta(l)
                if debug > 1: print "z=",z
                if s.parent() != ZZ:
                    z8 = K.zeta(8)
                    s = z8**(-FQM.signature())
                    if debug > 0: print "s={0}".format(s)
            else:
                z = K(exp(2*pi*I/l))
            if 1 == s2: 
                table = [2*s*(z**p).real()/w for p in range(l)]
            else:
                table = [s*(z**p-z**(-p))/w for p in range(l)]
    if debug > 0: print len(table), [table[i] for i in range(l)]
    if debug > 0: print '%f: init, table'%(walltime(t))

    if debug > 0: t = walltime()
    cdef long* ed = NULL
    fed = FQM.elementary_divisors()
    cdef int r = len(fed)
    ed = <long*> sig_malloc(sizeof(long) * r)
    if ed is NULL:
        raise MemoryError('Cannot allocate memory.')
    for i,d in enumerate(FQM.elementary_divisors()):
        ed[i] = long(d)
    if debug > 1: print fed

    J = FQM.__dict__['_FiniteQuadraticModule_ambient__J']
    #sig_on()
    t = walltime()
    cdef long** JJ = NULL
    JJ = <long**> sig_malloc(sizeof(long*) * r)
    if JJ is NULL:
        raise MemoryError('Cannot allocate memory.')
    for i in range(r):
        JJ[i] = NULL
        JJ[i] = <long*> sig_malloc(sizeof(long)*r)
        if JJ[i] == NULL:
            raise MemoryError('Cannot allocate memory.')
        for j in range(r):
            JJ[i][j] =  long((2*l*J[i,j]))
    #sig_off()
    if debug > 0: print 'Gram matrix conversion: {0}'.format(walltime(t))

    cdef list Ml = list()
    cdef long ni = 0
    cdef long kk = 0
    cdef int f = 0
    cdef list skip_list = list()
    if debug > 0: t = walltime()
    
    for i in range(n):
        #x = cython_elt(i,ed)
        if debug > 1: print "B=", B(i,i,JJ,ed,r)
        j = int(B(i,i,JJ,ed,r)/2) % l
        if debug > 1: print "j=",j
        kk = _neg_index(i,ed,r)
        #print i, kk, j
        f = 1
        if s2 == 1 or kk != i:
            if not i in skip_list:
                skip_list.append(i)
                skip_list.append(kk)
                if i != kk:
                    f = 2
                #print skip_list
                Ml.append((i,j,f))
                if j == 0:
                    ni = ni + 1
                    #print 'ni: ', ni
    if debug > 0: print '%f: +- reps'%(walltime(t)); t = walltime()
    if debug > 0: print 'ni = %d'%(ni)
    n = len(Ml)

    Ml.sort(norm_cmp)
    cdef long *Mli = NULL
    Mli = <long*> sig_malloc(sizeof(long)*n)
    cdef long[:] Mlf = np.ndarray(n, dtype=int)
    for ii, xx in enumerate(Ml):
        Mli[ii] = xx[0]
        Mlf[ii] = xx[2]
    
    if debug > 0: print 'n = %d'%(n)
    if debug > 0: print '%f: sorting'%(walltime(t)); t = walltime()

    #Now compute and cache bilinear forms
    #we use that it is symmetric and that we only need
    #it for isotropic elements paired with any other element
    cdef long** BB = NULL
    cdef long* Bli = NULL

    BB = <long**> sig_malloc(sizeof(long*)*n)
    if BB is NULL:
        raise MemoryError('Cannot allocate memory.')
    for i in prange(ni, nogil=True):
        BB[i] = NULL
        BB[i] = <long*> sig_malloc(sizeof(long)*(n-i))
        #if BB[i] == NULL:
        #    raise MemoryError('Cannot allocate memory.')
        Bli = Bl(Mli[i], JJ, ed, r)
        #sig_on()
        #if debug > 1:
        #    print [Bli[ii] for ii in range(r)]
        #sig_off()
        for j in range(n-i):
            #BB[i][j] = B(Ml[i+j][0],Ml[i][0],JJ,ed,r) % l
            BB[i][j] = (l-BBl(Mli[i+j],Bli,ed,r)) % l

    if not JJ is NULL:
        for i in range(r):
            if not JJ[i] is NULL:
                sig_free(JJ[i])
        sig_free(JJ)
     
    if debug > 0:
        print 'bilinear form computations: {0}'.format(walltime(t))
        t = walltime()

    if debug > 0: t = walltime()
    #A = np.ndarray((n,ni), dtype=long)
    #cdef long[:,:] HH = A
    cdef H = Matrix(K,n,ni)
    for j in range(ni):
        #print j, f
        for i in range(j, n):
            p = BB[j][i-j]
            H[i,j] = table[p]*Mlf[j]
            if i==j:
                H[j,j] += 2
            elif i<ni and i>j:
                H[j,i] = table[p]*Mlf[i]
            #if debug > 1: print "i={0}, j={1}, H[i,j] = {2}, p = {3}".format(i,j,HH[i,j],p)
    if debug > 0: print '%f: init of H'%(walltime(t)); t = walltime()
    #if debug > 0: print A
    #if debug > 1: print A.str()
    #print H.str()
    #H = Matrix(A)
    #H = H.change_ring(K)
    #if debug > 0: print '%f: conversion to H'%(walltime(t)); t = walltime(t)
    
    if return_H: return Ml, ni, H

    U = H.matrix_from_rows(range(ni,n))
    V = H.matrix_from_rows(range(ni))
    
    if proof and q > 0:
        t = walltime()
        M = Matrix(L,n,ni)
        eps = 1 if s2 == 1 else -1
        i = j = 0
        for j in range(ni):
            for i in range(j,n):
                p = BB[j][i-j]
                if debug > 1:
                    print "i={}, j={}, p={}".format(i, j, p)
                M[i,j] = table0[p]
                if Ml[j][2] == 2:
                    M[i,j] += eps*table0[l-p]
                if i<ni and i!=j:
                   M[j,i] = M[i,j]
                   if Ml[i][2] == 2:
                       M[j,i] += eps*p
        if debug > 1: print table0, M
        R = (Ml, ni, U, V, M)
        if debug > 0: print "Matrix M: {0}".format(walltime(t))
    else:
        R = (Ml, ni, U,V)

    if not BB is NULL:
        for i in range(ni):
            if not BB[i] is NULL:
                sig_free(BB[i])
        sig_free(BB)

    return R

cpdef reconstruction(x):
    if x in [0,1,-1]:
        return x.lift_centered()
    else:
        return 1/QQ((~x).lift_centered())

cpdef cython_invariants(FQM, use_reduction = True, proof = False, checks = False, debug=0, K = None):
    if use_reduction and K == None:
        found = False
        p = FQM.level()*10
        while not found:
            if p % lcm(4,FQM.level()) == 1:
                found = True
            else:
                p = next_prime(p)
        #print p
        K = GF(p)
    else:
        if K == None:
            K = CyclotomicField(lcm(8,FQM.level()))
    #print K
    I = cython_invariants_matrices(FQM, K, proof = checks or proof, debug=debug)
    if type(I)==list or type(I) == tuple:
        if not proof or K.characteristic() == 0:
            Ml, ni, U, V = I
        else:
            Ml, ni, U, V, M = I
    else:
        return I
    
    t = walltime()
    X = U.right_kernel()
    if debug > 0:
        print '%f: kernel'%(walltime(t))
        print "dimension kernel: {0}".format(X.dimension())

    t = walltime()
    Sp = span([V*x for x in X.basis()], K)
    if debug > 0: print '%f: span'%(walltime(t))

    if debug > 2:
        return U,V,X
    if (proof or checks) and K.characteristic() > 0:
        N = M.matrix_from_rows(range(ni))
    V1 = None
    if proof and K.characteristic() > 0 and Sp.dimension() > 0:
        if debug > 0: tt = walltime()
        l = FQM.level()
        NN = M.matrix_from_rows(range(ni,M.dimensions()[0]))
        V1 = []
        for v in Sp.basis():
            if debug > 0: print v
            vv = vector([reconstruction(x) for x in v])
            a = N*vv-vv
            b = NN*vv
            if not (a == 0 and b == 0):
                #trying another lift
                vv = vector([x.lift_centered() for x in v])
                a = N*vv-vv
                b = NN*vv
            if debug >1: print vv, a, b
            if not (a == 0 and b == 0):
                raise RuntimeError("Could not show that mod {0} invariant {1} lifts.".format(K.characteristic(),vv))
            else:
                V1.append(vv)
        V1 = span(V1)
        if debug > 0: print "proof: {0}".format(walltime(tt))
    if checks and K.characteristic() > 0:
        # Also check the multiplicity of the eigenvalue 1:
        f = N.characteristic_polynomial()
        val = f.valuation(f.parent().gen()-1)
        if val != V1.dimension():
            raise RuntimeError("The multiplicity of the eigenvalue 1 of N is {} but the dimension of the eigenspace we computed is {}".format(val, V1.dimension()))
    if V1 is not None:
        return Ml[:ni], V1
    else:
        return Ml[:ni], Sp

cpdef invariants(FQM, use_reduction = True, proof = False, checks=False, debug = 0):
    #print 'use_reduction = ', use_reduction
    I = cython_invariants(FQM, use_reduction, proof=proof, checks=checks, debug=debug)
    if type(I) == list or type(I) == tuple:
        Ml, Sp = I
    else:
        return I
    Mll = list()

    cdef long* ed = NULL
    fed = FQM.elementary_divisors()
    cdef int r = len(fed)
    ed = <long*> sig_malloc(sizeof(long) * r)
    for i,d in enumerate(FQM.elementary_divisors()):
        ed[i] = long(d)

    cdef long* vv = NULL
    for v in Ml:
        vv = _elt(v[0],ed,r)
        vvl = [vv[i] for i in range(r)]
        Mll.append(FQM(vvl))
        
    return Mll, Sp
            
            
            
        
