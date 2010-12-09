cdef extern from "smalljac.h":
    ctypedef void *smalljac_Qcurve_t
    int smalljac_Lpoly (long* a, char *curve, unsigned long q, unsigned long flags)
    void smalljac_Qcurve_clear (smalljac_Qcurve_t c)
    smalljac_Qcurve_t smalljac_Qcurve_init (char *str, int *err)
    char *smalljac_Qcurve_str (smalljac_Qcurve_t c)    
    int smalljac_Qcurve_genus (smalljac_Qcurve_t c)
    long smalljac_Lpolys (smalljac_Qcurve_t c,
                          unsigned long start,
                          unsigned long end,				
                          unsigned long flags,
                          int (*callback)(smalljac_Qcurve_t c, unsigned long p,	int good,  long a[],  int n, void *arg),
                          void *arg)

cdef class SmallJac:
    cdef smalljac_Qcurve_t c
    cdef object tmp
    def __cinit__(self):
        self.c = <smalljac_Qcurve_t>0
        
    def __init__(self, v):
        curve = str(v)
        cdef int err
        self.c = smalljac_Qcurve_init(curve, &err)
        if err:
            # todo -- see smalljac.h for how to parse this
            raise RuntimeError, "Error code %s"%err

    def __repr__(self):
        return "Smalljac curve defined by %s"%smalljac_Qcurve_str(self.c)

    def __dealloc__(self):
        if self.c:
            smalljac_Qcurve_clear(self.c)

    cpdef int genus(self):
        """
        EXAMPLES::

            sage: import psage.libs.smalljac.wrapper as smalljac
            sage: smalljac.SmallJac('x^3+x+1').genus()
            1
            sage: smalljac.SmallJac('x^5+x+1').genus()
            2
            sage: smalljac.SmallJac('x^7+x+1').genus()
            3
            sage: smalljac.SmallJac('x^9+x+1').genus()
            4
            sage: smalljac.SmallJac('x^11+x+1').genus()
            Traceback (most recent call last):
            ...
            RuntimeError: Error code -4
        """
        return smalljac_Qcurve_genus(self.c)

    def ap_dict(self, unsigned long a, unsigned long b=0):
        """
        Return a dictionary of key:value pairs, with keys the good
        primes <= maxp, and values the traces of Frobenius at the good
        primes.

        EXAMPLES::

            sage: import psage.libs.smalljac.wrapper as smalljac
            sage: C = smalljac.SmallJac(x^3 - 13392*x - 1080432)
            sage: C.ap_dict(37)
            {5L: 1, 7L: -2, 13L: 4, 17L: -2, 37L: 3, 19L: 0, 23L: -1, 29L: 0, 31L: 7}
            sage: C.ap_dict(19, 37)
            {31L: 7, 37L: 3, 19L: 0, 29L: 0, 23L: -1}
            sage: C = smalljac.SmallJac("x^5 - 1")
            sage: C.ap_dict(37)
            Traceback (most recent call last):
            ...
            ValueError: ap_dict only makes sense for genus 1 curves
        """
        if b == 0:
            a, b = b, a
        if self.genus() != 1:
            raise ValueError, "ap_dict only makes sense for genus 1 curves"
        self.tmp = {}
        smalljac_Lpolys(self.c, a, b, 0, callback_ap_dict, <void*>self)
        return self.tmp

    def Lpolys(self, unsigned long a, unsigned long b=0, flags=None):
        """
        EXAMPLES::

        
        """
        if b == 0:
            a, b = b, a
        cdef unsigned long _flags = 0
        self.tmp = {}
        smalljac_Lpolys(self.c, a, b, _flags, callback_Lpolys, <void*>self)
        return self.tmp

    def groups(self, unsigned long start, unsigned long end, flags=None):
        raise NotImplementedError

cdef int callback_Lpolys(smalljac_Qcurve_t c, unsigned long p, int good,  long a[],  int n, void *arg):
    cdef int i
    (<SmallJac>arg).tmp[int(p)] = [a[i] for i in range(n)]
    return 1
    
cdef int callback_ap_dict(smalljac_Qcurve_t c, unsigned long p, int good,  long a[],  int n, void *arg):
    if good:
        (<SmallJac>arg).tmp[int(p)] = -a[0]
    return 1

def Lpoly(f, unsigned long q):
    """
    EXAMPLES::

    
    """
    curve = str(f)
    cdef int i, n
    cdef long a[10]
    n = smalljac_Lpoly(a, curve, q, 0)
    if n == 0:
        raise ValueError, "bad reduction"
    if n < 0:
        raise RuntimeError
    return [a[i] for i in range(n)]

    
 
