# cython: profile = False
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

include "stdsage.pxi"

from sage.modules.free_module import span
from sage.matrix.constructor import Matrix
from sage.rings.qqbar import QQbar
from sage.all import exp, Integer, pi, I, cputime, CyclotomicField, ZZ, is_prime_power #, sage_malloc, sage_free, ZZ
from sage.rings.number_field.number_field import NumberField_cyclotomic


cdef long cython_el_index(list c, list ed):
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
    cdef int r = len(ed)
    for jj in range(0,r):
        m=ed[jj]
        #print jj,n,c[jj]
        ii = ii + (c[jj]%m)*md
        md=md*m
    return ii

cdef list cython_elt(long ii,ed):
    r"""
    Returns the element corresponding to the canonical index `ii`
    in a finite abelian group with elementary divisors specified by `ed`.
    """
    elt = list()
    cdef long md = 1
    cdef long jj = 0
    cdef int r = len(ed)
    cdef int c = 0
    for jj in range(0,r):
        md = ed[jj]
        c = ii%md
        elt.append(c)
        ii = ii - c
        ii = ii/md
    return elt

cdef long cython_neg_index(long ii, ed):
    r"""
    Returns the negative of the element corresponding to the canonical index `ii`
    in a finite abelian group with elementary divisors specified by `ed`.
    """
    cdef long jj=0
    cdef list elt = cython_elt(ii,ed)
    cdef int r = len(ed)
    for jj in range(0,r):
        elt[jj]=-elt[jj]
    return cython_el_index(elt,ed)

cpdef norm_cmp(x, y):
    if x[1] < y[1]:
        return int(-1)
    else:
        return int(1)

cdef int B(int i, int j, int **JJ, list ed):
    r"""
    Returns the bilinear form of the elements with canonical index `i` and `j`
    with respect to the Gram matrix JJ and elementary divisors `ed`.
    """
    cdef int r = len(ed)
    cdef list ll = cython_elt(j, ed)
    cdef list kk = cython_elt(i, ed)
    cdef int ii, jj = 0
    cdef int res = 0
    for ii in range(r):
        for jj in range(r):
            res = res + ll[ii]*kk[jj]*JJ[ii][jj]
    return res

cpdef cython_invariants_dim(FQM, K = QQbar, debug=0):
    if FQM.signature() != 0:
        return 0
    dim = 0
    if not is_prime_power(FQM.level()):
        for p,n in FQM.level().factor():
            N = None
            for j in xrange(1,n+1):
                J = FQM.jordan_decomposition()
                C = J.constituent(p**j)[0]
                if N == None and C.level() != 1:
                    N = C
                elif C.level() != 1:
                    N = N + C
            if dim == 0:
                dim = cython_invariants_dim(N,K,debug)
            else:
                dim = dim*cython_invariants_dim(N,K,debug)
            if dim == 0:
                return 0
    else:
        Sp = cython_invariants(FQM, K, debug)[1]
        dim = dim + Sp.dimension()
    return dim

cpdef cython_invariants_matrices(FQM, K = QQbar, debug=0):
    cdef int i, j = 0
    cdef int l = int(FQM.level())
    if debug > 1: print "l = ", l
    
    try:
        s = FQM.sigma_invariant()
        s2 = Integer( s**2)
    except:
        return span( [], K)

    q = K.characteristic()
    w = K(FQM.order()).sqrt()

    t = cputime()

    cdef list table = list()
    if 0 == q:
        if isinstance(K,NumberField_cyclotomic):
            if debug > 0: print 'cyclotomic'
            z = K.gen()
            if 1 == s2: 
                table = [2*s*((z**p) + (z**p).conjugate())/w for p in range(l)]
            else:
                table = [2*s*((z**p) - (z**p).conjugate())/w for p in range(l)]
        else:
            if K == QQbar:
                if debug > 0: print 'QQbar'
                z = K.zeta(l)
                if debug > 1: print "z=",z
                if s.parent() != ZZ:
                    z8 = s.parent().gen()
                    pw = z8.coordinates_in_terms_of_powers()(s)
                    z8 = K.zeta(8)
                    s = sum([pw[i]*z8**i for i in range(4)])
            else:
                z = K( exp(2*pi*I/l))
            if 1 == s2: 
                table = [2*s*(z**p).real()/w for p in range(l)]
            else:
                table = [2*s*(z**p).imag()/w for p in range(l)]
    if q > 0:
        if debug > 0: print 'positive characteristic: ', q
        if 1 != q % l:
            raise ValueError( '%d: must be = 1 modulo %d' %(q, l))
        pr = K.primitive_element()
        z = pr**((q-1)/l)
        table = [z**p for p in range(l)]
        for i in range(l):
            zt = table[i]
            table[i] = K(zt + s2 * zt**-1)/K(w)
    if debug > 1: print len(table), table

    if debug > 0: print '%f: table'%(cputime(t))

    t = cputime()
    ed = list()
    for i in FQM.elementary_divisors():
        ed.append(Integer(i))
    if debug > 1: print ed

    cdef int r = len(ed)
    J = FQM.__dict__['_FiniteQuadraticModule_ambient__J']
    t = cputime()
    cdef int** JJ = NULL
    JJ = <int**> sage_malloc(sizeof(int*) * r)
    if JJ == NULL:
        raise MemoryError('Cannot allocate memory.')
    for i in range(r):
        JJ[i] = NULL
        JJ[i] = <int*> sage_malloc(sizeof(int)*r)
        if JJ[i] == NULL:
            raise MemoryError('Cannot allocate memory.')
        for j in range(r):
            JJ[i][j] =  int(2*l*J[i,j])
    if debug > 0: print 'Gram matrix conversion: {0}'.format(cputime(t))

    t = cputime()
    Ml = list()
    cdef int ni = 0
    cdef long kk = 0
    cdef int f = 0
    cdef list skip_list = list()
    for i in range(FQM.order()):
        #x = cython_elt(i,ed)
        if debug > 1: print "B=", B(i,i,JJ,ed)
        j = int(B(i,i,JJ,ed)/2) % l
        if debug > 1: print "j=",j
        kk = cython_neg_index(i,ed)
        #print i, kk, j
        f = 1
        if s2 == 1 or kk != i:
            if not i in skip_list:
                skip_list.append(i)
                skip_list.append(kk)
                if i == kk:
                    f = 2
                #print skip_list
                Ml.append((i,j,f))
                if j == 0:
                    ni = ni + 1
                    #print 'ni: ', ni
    if debug > 0: print '%f: +- reps'%(cputime(t))
    if debug > 0: print 'ni = %d'%(ni)
    cdef int n = len(Ml)
    
    t = cputime()
    Ml.sort(norm_cmp)
    n = len(Ml)
    if debug > 0: print 'n = %d'%(n)
    if debug > 0: print '%f: sorting'%(cputime(t))

    t = cputime()
    cdef int ii,jj = 0
    H = Matrix(K, n, ni)
    for j in range(ni):
        H[j,j] = Ml[j][2]
        #print j, f
        for i in range(n):
            p = -B(Ml[i][0],Ml[j][0], JJ, ed) % l
            H[i,j] += table[p]
    if debug > 0: print '%f: init of H'%(cputime(t))
    if debug > 1: print H.str()
    #print H.str()

    U = H.matrix_from_rows(range(ni,n))
    V = H.matrix_from_rows(range(ni))

    if not JJ is NULL:
        for i in range(r):
            if not JJ[i] is NULL:
                sage_free(JJ[i])
        sage_free(JJ)
    return Ml, ni, U,V

cpdef cython_invariants(FQM, K = QQbar, debug=0):    
    I = cython_invariants_matrices(FQM, K, debug)
    if type(I)==list or type(I) == tuple:
        Ml, ni, U,V = I
    else:
        return I
    
    t = cputime()
    X = U.right_kernel()
    if debug > 0: print '%f: kernel'%(cputime(t))

    t = cputime()
    Sp = span([V*x for x in X.basis()], K)
    if debug > 0: print '%f: span'%(cputime(t))

    if debug > 1:
        return U,V,X
    return Ml[:ni], Sp

cpdef invariants(FQM, K = QQbar, debug = 0):
    I = cython_invariants(FQM, K, debug)
    if type(I) == list or type(I) == tuple:
        Ml, Sp = I
    else:
        return I
    Mll = list()
    ed = list()
    for i in FQM.elementary_divisors():
        ed.append(Integer(i))
    for v in Ml:
        vv = cython_elt(v[0],ed)
        Mll.append(FQM(vv))
    return Mll, Sp
            
            
            
        
