include "wrapper0.pyx"
include 'stdsage.pxi'

cdef class GroupModp:
    pass

cdef class FrobPoly:
    cdef long a
    cdef unsigned long p
    def __init__(self, a, p):
        self.a = a
        self.p = p

    def __reduce__(self):
        return FrobPoly_new, (self.a, self.p)

    def __repr__(self):
        return "Frobenius in characteristic %s with charpoly %s"%(self.p, self.characteristic_polynomial())

    def prime(self):
        return int(self.p)

    def trace(self):
        return -self.a

    def characteristic_polynomial(self, var='x'):
        from sage.rings.all import ZZ
        return ZZ[var]([self.p, self.a, 1])

    charpoly = characteristic_polynomial

cpdef FrobPoly_new(long a, unsigned long p):
    cdef FrobPoly L = PY_NEW(FrobPoly)
    L.a = a; L.p = p
    return L


cdef int callback_Lpolys(smalljac_Qcurve_t c, unsigned long p, int good,  long a[],  int n, void *arg):
    cdef int i
    if good:
        (<SmallJac>arg).tmp[int(p)] = FrobPoly_new(a[0], p)
    else:
        (<SmallJac>arg).tmp[int(p)] = None
    return 1

cdef int callback_Lpolys_single(smalljac_Qcurve_t c, unsigned long p, int good,  long a[],  int n, void *arg):
    if good:
        (<SmallJac>arg).tmp = FrobPoly_new(a[0], p)
    else:
        (<SmallJac>arg).tmp = None
    return 1
    
cdef int callback_ap_dict(smalljac_Qcurve_t c, unsigned long p, int good,  long a[],  int n, void *arg):
    if good:
        (<SmallJac>arg).tmp[int(p)] = -a[0]
    else:
        (<SmallJac>arg).tmp[int(p)] = None
    return 1

cdef int callback_ap_single(smalljac_Qcurve_t c, unsigned long p, int good,  long a[],  int n, void *arg):
    if good:
        (<SmallJac>arg).tmp = -a[0]
    else:
        (<SmallJac>arg).tmp = None
    return 1

