#################################################################################
#
# (c) Copyright 2011 R. Andrew Ohana
#
#  This file is part of PSAGE
#
#  PSAGE is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  PSAGE is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################

"""
Fast Cython code to choose a canonical minimal model.
"""

from sage.ext.stdsage cimport PY_NEW
from sage.libs.gmp.all cimport *

from sage.schemes.elliptic_curves.constructor import EllipticCurve
from sage.rings.integer cimport Integer
from libc.stdint cimport int64_t, uint64_t

cdef void f_map64(int64_t *o):
    cdef bint s
    cdef uint64_t x[2],y[2],g[3],f[3],h[3],t
    o[2] = 1
    o[3] = 0
    s = (o[0] < 0L)
    # to simplify boolean checks
    o[0] *= 1-2*s
    o[1] *= 1-2*s
    x[1] = o[0] >> 32
    x[0] = o[0]&0xFFFFFFFF
    # to avoid multiplying by negative ints
    y[1] = (1-2*(o[1]<0L))*o[1] >> 32
    y[0] = ((1-2*(o[1]<0L))*o[1])&0xFFFFFFFF
    # f = o[0]*o[1]
    f[0] = x[0]*y[0]
    t = x[0]*y[1] + x[1]*y[0] + (f[0]>>32)
    f[1] = x[1]*y[1] + (t>>32)
    f[0] = (f[0]&0xFFFFFFFF) + (t<<32)
    # fix sign if needed
    if o[1] < 0L:
        f[0] = -f[0]
        f[1] = -f[1]-(f[0]>0)
        f[2] = -1
    else:
        f[2] = 0
    # g = o[0]*o[0]
    g[0] = x[0]*x[0]
    t = ((x[0]*x[1])<<1)+(g[0]>>32)
    g[1] = x[1]*x[1] + (t>>32)
    g[0] = (g[0]&0xFFFFFFFF) + (t<<32)
    # g *= 36
    g[2] = 9*(g[1]>>4)>>58
    g[1] = 36*g[1] + (9*(g[0]>>4)>>58)
    g[0] *= 36
    # h = o[1]*o[1]
    h[0] = y[0]*y[0]
    t = ((y[0]*y[1])<<1)+(h[0]>>32)
    h[1] = y[1]*y[1] + (t>>32)
    h[0] = (h[0]&0xFFFFFFFF) + (t<<32)
    # h *= 180
    h[2] = 45*(h[1]>>6)>>56
    h[1] = 180*h[1] + (45*(h[0]>>6)>>56)
    h[0] *= 180
    # g += h
    g[2] += h[2]
    g[1] += h[1]
    g[0] += h[0]
    g[1] += g[0] < h[0]
    g[2] += g[1] < h[1] or (g[1] == h[1] and g[0] < h[0])
    # h = 161*f
    h[2] = 161*f[2] + (161*(f[1]>>8)>>56)
    h[1] = 161*f[1] + (161*(f[0]>>8)>>56)
    h[0] = 161*f[0]
    # f *= 2
    f[2] = (f[2]<<1) + (f[1]>>63)
    f[1] = (f[1]<<1) + (f[0]>>63)
    f[0] <<= 1
    # f = 2*o[0]*o[1]
    if o[1] > 0L:
        # g = -g
        g[0] = -g[0]
        g[1] = -g[1]-(g[0]>0)
        g[2] = -g[2]-(g[1]>0 or g[0]>0)
        # g += h
        g[2] += h[2]
        g[1] += h[1]
        g[0] += h[0]
        g[1] += g[0] < h[0]
        g[2] += g[1] < h[1] or (g[1] == h[1] and g[0] < h[0])
        # g = 161*o[0]*o[1]-36*o[0]*o[0]-180*o[1]*o[1]
        while g[2] < 0x8000000000000000L:
            t = 161*o[0]-360*o[1]
            o[1] = 161*o[1]-72*o[0]
            o[0] = t
            # f = 644*g-f
            h[2] = 644*g[2] + (161*(g[1]>>8)>>54)
            h[1] = 644*g[1] + (161*(g[0]>>8)>>54)
            h[0] = 644*g[0]
            f[2] = h[2]-f[2]
            f[1] = h[1]-f[1]
            f[0] = h[0]-f[0]
            f[1] -= f[0] > h[0]
            f[2] -= f[1] > h[1] or (f[1] == h[1] and f[0] > h[0])
            # g = 161*f-g
            h[2] = 161*f[2] + (161*(f[1]>>8)>>56)
            h[1] = 161*f[1] + (161*(f[0]>>8)>>56)
            h[0] = 161*f[0]
            g[2] = h[2]-g[2]
            g[1] = h[1]-g[1]
            g[0] = h[0]-g[0]
            g[1] -= g[0] > h[0]
            g[2] -= g[1] > h[1] or (g[1] == h[1] and g[0] > h[0])
            o[3] += o[2]
            o[2] = o[3]-o[2]
    elif o[1] < 0L:
        # g += h
        g[2] += h[2]
        g[1] += h[1]
        g[0] += h[0]
        g[1] += g[0] < h[0]
        g[2] += g[1] < h[1] or (g[1] == h[1] and g[0] < h[0])
        # g = 161*o[0]*o[1]+36*o[0]*o[0]+180*o[1]*o[1]
        while g[2] >= 0x8000000000000000L:
            t = 161*o[0]+360*o[1]
            o[1] = 72*o[0]+161*o[1]
            o[0] = t
            # f = 644*g-f
            h[2] = 644*g[2] + (161*(g[1]>>8)>>54)
            h[1] = 644*g[1] + (161*(g[0]>>8)>>54)
            h[0] = 644*g[0]
            f[2] = h[2]-f[2]
            f[1] = h[1]-f[1]
            f[0] = h[0]-f[0]
            f[1] -= f[0] > h[0]
            f[2] -= f[1] > h[1] or (f[1] == h[1] and f[0] > h[0])
            # g = 161*f-g
            h[2] = 161*f[2] + (161*(f[1]>>8)>>56)
            h[1] = 161*f[1] + (161*(f[0]>>8)>>56)
            h[0] = 161*f[0]
            g[2] = h[2]-g[2]
            g[1] = h[1]-g[1]
            g[0] = h[0]-g[0]
            g[1] -= g[0] > h[0]
            g[2] -= g[1] > h[1] or (g[1] == h[1] and g[0] > h[0])
            o[2] = o[3]-o[2]
            o[3] -= o[2]
    o[0] = o[0]-o[1]>>1
    if not o[0]%13 and not o[1]%8 and o[0]/13 == -(o[1]/8):
        t = o[0]/13
        o[0] = 5*t
        o[1] = 8*t
        o[2] = o[3]-o[2]
        o[3] -= o[2]
    elif not o[0]%29 and not o[1]%18 and o[0]/29 == -(o[1]/18):
        t = o[1]/18
        o[0] = 11*t
        o[1] = 18*t
        o[2] = o[3]-o[2]
        o[3] -= o[2]
    o[0] *= 1-2*s
    o[1] *= 1-2*s

cdef void f_mapC(mpz_t rop1, mpz_t rop2, mpz_t rop3, mpz_t rop4, mpz_t op1, mpz_t op2):
    """
    Let `G` be the function
        `G(\delta) = 2*a*a + 2*a*b + 3*b*b,`
    where `\delta = a+b*\gamma`, then the output are the unique values such that
        `G(rop1+rop2*\gamma) <= G(v^12 * (rop1+rop2*\gamma))`
    for any unit `v`, and
        `op1+op2*\gamma = (rop3+rop4*\gamma)^12 * (rop1+rop2*\gamma).`
    The implementation transforms to a new space, then iteratively finds these
    values by storing either
        `G(\gamma^12*\delta) - G(\delta)`
    or
        `G(\gamma^-12*\delta) - G(\delta)`
    and terminating when this difference is positive, and transforming back.
    """
    cdef bint s
    cdef int64_t o[4]
    cdef mpz_t f,g,t
    mpz_add(rop1,op1,op1)
    mpz_add(rop1,rop1,op2)
    if not mpz_sgn(rop1):
        mpz_set(rop1, op1)
        mpz_set(rop2, op2)
        mpz_set_si(rop3, 1)
        mpz_set_si(rop4, 0)
        return
    mpz_set(rop2,op2)
    if mpz_sizeinbase(rop1, 2) < (sizeof(o[0])<<3):
        if mpz_sizeinbase(rop2, 2) < (sizeof(o[1])<<3):
            mpz_export(&o[0], NULL, -1, sizeof(o[0]), -1, 0, rop1)
            mpz_export(&o[1], NULL, -1, sizeof(o[1]), -1, 0, rop2)
            o[0] *= mpz_sgn(rop1)
            o[1] *= mpz_sgn(rop2)
            f_map64(o)
            mpz_set_si(rop1,o[0])
            mpz_set_si(rop2,o[1])
            mpz_set_si(rop3,o[2])
            mpz_set_si(rop4,o[3])
            return
    mpz_set_si(rop3,1)
    mpz_set_si(rop4,0)
    s = mpz_sgn(rop1) < 0
    # to simplify boolean checks
    if s:
        mpz_neg(rop1,rop1)
        mpz_neg(rop2,rop2)
    mpz_init(f); mpz_init(g); mpz_init(t)
    
    # g = 36*rop1*rop1+180*rop2*rop2
    mpz_mul(f,rop2,rop2)
    mpz_mul(g,rop1,rop1)
    mpz_mul_si(g,g,36)
    mpz_addmul_ui(g,f,180)

    mpz_mul(f,rop1,rop2)
    if mpz_sgn(rop2) > 0:
        # g = 36*rop1*rop1-161*rop1*rop2+180*rop2*rop2
        mpz_submul_ui(g,f,161)
        mpz_mul_si(f,f,2)
        while mpz_sgn(g) < 0:
            # rop1+rop2*g /= g^12
            mpz_mul_si(t,rop1,161)
            mpz_submul_ui(t,rop2,360)
            mpz_mul_si(rop2,rop2,161)
            mpz_submul_ui(rop2,rop1,72)
            mpz_set(rop1,t)
            
            # fix f and g for new rop1+rop2*g
            mpz_addmul_ui(f,g,644)
            mpz_neg(f,f)
            mpz_addmul_ui(g,f,161)
            mpz_neg(g,g)
            
            # rop3+rop4*g *= g
            mpz_add(rop4,rop3,rop4)
            mpz_sub(rop3,rop4,rop3)
    elif mpz_sgn(rop2) < 0:
        # g = 36*rop1*rop1+161*rop1*rop2+180*rop2*rop2
        mpz_addmul_ui(g,f,161)
        mpz_mul_si(f,f,2)
        while mpz_sgn(g) < 0:
            # rop1+rop2*g *= g^12
            mpz_mul_si(t,rop1,161)
            mpz_addmul_ui(t,rop2,360)
            mpz_mul_si(rop2,rop2,161)
            mpz_addmul_ui(rop2,rop1,72)
            mpz_set(rop1,t)
            
            # fix f and g for new rop1+rop2*g
            mpz_submul_ui(f,g,644)
            mpz_neg(f,f)
            mpz_submul_ui(g,f,161)
            mpz_neg(g,g)
            
            # rop3+rop4*g /= g
            mpz_sub(rop3,rop4,rop3)
            mpz_sub(rop4,rop4,rop3)
    # fix sign if needed
    mpz_clear(f)
    mpz_sub(rop1,rop1,rop2)
    mpz_divexact_ui(rop1,rop1,2ul)
    if not mpz_mod_ui(g, rop1, 13ul) and not mpz_mod_ui(t, rop2, 8ul):
        mpz_divexact_ui(t, rop1, 13ul)
        mpz_divexact_ui(g, rop2, 11ul)
        mpz_neg(g,g)
        if not mpz_cmp(t,g):
            mpz_mul_ui(rop1, t, 5ul)
            mpz_mul_ui(rop2, t, 8ul)
            mpz_sub(rop3,rop4,rop3)
            mpz_sub(rop4,rop4,rop3)
    elif not mpz_mod_ui(g, rop1, 29ul) and not mpz_mod_ui(t, rop2, 18ul):
        mpz_divexact_ui(g, rop1, 29ul)
        mpz_divexact_ui(t, rop2, 18ul)
        mpz_neg(t,t)
        if not mpz_cmp(t,g):
            mpz_mul_ui(rop1, g, 11ul)
            mpz_mul_ui(rop2, g, 18ul)
            mpz_sub(rop3,rop4,rop3)
            mpz_sub(rop4,rop4,rop3)
    mpz_clear(g); mpz_clear(t)
    if s:
        mpz_neg(rop1,rop1)
        mpz_neg(rop2,rop2)

cpdef f_map(Integer a, Integer b):
    """
    Computes 'cannonical' representitives for the equivalence relation
        `\lambda / \mu \in U_F^{12},`
    on the ring of integers for `F=Q(\sqrt{5})`.
    
    INPUT:
    
    - ``a+b*\gamma`` -- an element of the ring of integers of `F`
    
    OUTPUT:
    
    A 'cannonical' element `\delta*u^{-12}` of the maximal order of `F` and
    the unit `u`
    
    EXAMPLES::

        sage: from psage.ellcurve.minmodel.sqrt5 import f_map
        sage: K.<a> = NumberField(x^2-x-1)
        sage: E = EllipticCurve(K,[a,a+1,a+2,a+3,a+4])
        sage: olddelta = E.discriminant()
        sage: olddeltapair = (Integer(olddelta[0]),Integer(olddelta[1]))
        sage: temp = f_map(olddeltapair[0],olddeltapair[1])
        sage: temp[0] == olddeltapair[0] and temp[1] == olddeltapair[1]
        True
        sage: F = E.change_weierstrass_model([a^5,0,0,0])
        sage: delta = F.discriminant()
        sage: deltapair = (Integer(delta[0]),Integer(delta[1]))
        sage: temp = f_map(deltapair[0],deltapair[1])
        sage: temp[0] == deltapair[0] and temp[1] == deltapair[1]
        False
        sage: temp[0] == olddeltapair[0] and temp[1] == olddeltapair[1]
        True
        sage: F.change_weierstrass_model([temp[2]+a*temp[3],0,0,0]) == E
        True
    
    """
    cdef Integer c,d,e,f
    c = <Integer>PY_NEW(Integer)
    d = <Integer>PY_NEW(Integer)
    e = <Integer>PY_NEW(Integer)
    f = <Integer>PY_NEW(Integer)
    f_mapC(c.value, d.value, e.value, f.value, a.value, b.value)
    return c,d,e,f

def canonical_model(E):
    """
    Given a global minimal model E over Q(sqrt5) returns a canonical global
    minimal model.

    EXAMPLES::

      sage: from psage.ellcurve.minmodel.sqrt5 import canonical_model
      sage: K.<a> = NumberField(x**2-x-1)
      sage: E = EllipticCurve(K,[a,a-1,a,-1,-a+1])
      sage: F = canonical_model(E)
      sage: E == F
      False
      sage: E.is_isomorphic(F)
      True
      sage: F == canonical_model(F)
      True
    """
    try:
        K = E.base_field()
        if K.disc() != 5:
            raise TypeError, "E must be over Q(sqrt(5))"
        D1, D2 = E.discriminant()
        a1,a2,a3,a4,a6 = E.a_invariants()
    except AttributeError:
        raise TypeError, "E must be an elliptic curve"
    cdef Integer u1, u2
    u1 = Integer(D1)
    u2 = Integer(D2)
    cdef mpz_t t1,t2
    mpz_init(t1); mpz_init(t2)
    f_mapC(t1, u1.value, t2, u2.value, u1.value, u2.value)
    mpz_add(u1.value, t2, u2.value)
    mpz_mul(t1, u1.value, t2)
    mpz_submul(t1, u2.value, u2.value)
    if mpz_sgn(t1) > 0:
        mpz_neg(u2.value, u2.value)
    else:
        mpz_neg(u1.value, u1.value)
    mpz_clear(t1); mpz_clear(t2)
    u = K(u1)+K.gen(0)*K(u2)
    a1,a2,a3,a4,a6 = u*a1,u*u*a2,u*u*u*a3,a4*u**4,a6*u**6
    P2 = K.ideal(2)
    P3 = K.ideal(3)
    a1p = a1.mod(P2)
    s = (a1p - a1)/K(2)
    sa1 = s*a1
    a2p = (a2 - sa1 - s**2).mod(P3)
    r = (a2p - a2 + sa1 + s**2)/K(3)
    ra1p = r*a1p
    a3p = (a3+ra1p).mod(P2)
    t = r*s + (a3p - a3 - ra1p)/K(2)
    a4p = a4 - s*a3 + 2*r*a2 - (t+r*s)*a1 + 3*r**2 - 2*s*t
    a6p = a6 + r*a4 + r**2*a2 + r**3 - t*a3 - t**2 - r*t*a1
    return  EllipticCurve(K, [a1p, a2p, a3p, a4p, a6p])

cpdef canonical_model_c_invariants(Integer a, Integer b, Integer c, Integer d):
    """
    Given the pair (c4,c6) of a global minimal model, canonical_model returns
    a 'cannonical' representation of the curve in terms of (c4,c6).
    
    INPUT:
    
    - ``c4=a+b*\gamma`` -- an element of the ring of integers of `F=Q(\sqrt{5})`
    - ``c6=c+d*\gamma`` -- an element of the ring of integers of `F`
    
    OUTPUT:
    
    A 'cannonical' model of a curve in terms of (c4,c6)
    
    EXAMPLES::
    
        sage: K.<a> = NumberField(x^2-x-1)
        sage: E = EllipticCurve(K,[a,a+1,a+2,a+3,a+4])
        sage: F = E.change_weierstrass_model(a,a+1,a+2,a+3)
        sage: E.c_invariants() == F.c_invariants()
        False
        sage: Ecs = [Integer(i) for i in E.c4()]
        sage: Ecs += [Integer(i) for i in E.c6()]
        sage: Fcs = [Integer(i) for i in F.c4()]
        sage: Fcs += [Integer(i) for i in F.c6()]
        sage: from psage.ellcurve.minmodel.sqrt5 import canonical_model_c_invariants
        sage: canonical_model_c_invariants(Ecs[0],Ecs[1],Ecs[2],Ecs[3])
        (-118, -45, -3001, 116)
        sage: canonical_model_c_invariants(Fcs[0],Fcs[1],Fcs[2],Fcs[3])
        (-118, -45, -3001, 116)
    

    """
    cdef Integer x,y,z,w
    cdef mpz_t p,q
    x = <Integer>PY_NEW(Integer)
    y = <Integer>PY_NEW(Integer)
    z = <Integer>PY_NEW(Integer)
    w = <Integer>PY_NEW(Integer)
    mpz_init(p)
    mpz_init(q)
    
    # p+q*g = (a+b*g)^3
    mpz_mul(x.value,a.value,a.value)
    mpz_mul(y.value,b.value,b.value)
    mpz_mul(p,y.value,a.value)
    mpz_mul_si(p,p,3)
    mpz_mul(q,x.value,b.value)
    mpz_mul_si(q,q,3)
    mpz_add(q,q,p)
    mpz_addmul(p,x.value,a.value)
    mpz_mul(y.value,y.value,b.value)
    mpz_add(p,p,y.value)
    mpz_addmul_ui(q,y.value,2)
    
    
    # compute discriminant 
    
    # p+q*g -= (c+d*g)^2
    mpz_submul(p,c.value,c.value)
    mpz_mul(x.value,c.value,d.value)
    mpz_submul_ui(q,x.value,2)
    mpz_mul(x.value,d.value,d.value)
    mpz_sub(p,p,x.value)
    mpz_sub(q,q,x.value)
    
    # p+q*g /= 1728
    mpz_divexact_ui(p,p,1728)
    mpz_divexact_ui(q,q,1728)
    
    
    # find the needed transformation, and x+y*g to corresponding unit
    f_mapC(p,q,x.value,y.value,p,q)
    
    
    # compute the transform
    
    # p+q*g = (x+y*g)^-2
    mpz_add(x.value,x.value,y.value)
    mpz_mul(p,y.value,y.value)
    mpz_mul(q,x.value,y.value)
    mpz_mul_si(q,q,2)
    mpz_sub(q,q,p)
    mpz_addmul(p,x.value,x.value)
    mpz_neg(q,q)
    
    # x+y*g = (p+q*g)^3
    mpz_mul(z.value,p,p)
    mpz_mul(w.value,q,q)
    mpz_mul(x.value,w.value,p)
    mpz_mul_si(x.value,x.value,3)
    mpz_mul(y.value,z.value,q)
    mpz_mul_si(y.value,y.value,3)
    mpz_add(y.value,y.value,x.value)
    mpz_addmul(x.value,z.value,p)
    mpz_mul(w.value,w.value,q)
    mpz_add(x.value,x.value,w.value)
    mpz_addmul_ui(y.value,w.value,2)
    
    # z+w*g = (x+y*g)*(c+d*g)
    mpz_mul(z.value,y.value,d.value)
    mpz_mul(w.value,x.value,d.value)
    mpz_addmul(w.value,y.value,c.value)
    mpz_add(w.value,w.value,z.value)
    mpz_addmul(z.value,x.value,c.value)
    
    # p+q*g *= p+q*g
    mpz_mul(x.value,p,q)
    mpz_mul(q,q,q)
    mpz_mul(p,p,p)
    mpz_add(p,p,q)
    mpz_addmul_ui(q,x.value,2)
    
    # x+y*g = (p+q*g)*(a+b*g)
    mpz_mul(x.value,q,b.value)
    mpz_mul(y.value,p,b.value)
    mpz_addmul(y.value,q,a.value)
    mpz_add(y.value,y.value,x.value)
    mpz_addmul(x.value,p,a.value)
    
    mpz_clear(p)
    mpz_clear(q)

    return x,y,z,w
