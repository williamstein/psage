from defs cimport *

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

    def ap(self, unsigned long p, unsigned long b=0):
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
            # Compute for a single p only
            smalljac_Lpolys(self.c, p, p, SMALLJAC_A1_ONLY, callback_ap_single, <void*>self)
            return self.tmp
        # Compute for a range of primes, starting with p.
        self.tmp = {}
        smalljac_Lpolys(self.c, p, b,  SMALLJAC_A1_ONLY, callback_ap_dict, <void*>self)
        return self.tmp

    def frob(self, unsigned long p, unsigned long b=0):
        """
        EXAMPLES::
        """
        if b == 0:
            # Compute for a single p only
            smalljac_Lpolys(self.c, p, p, 0, callback_Lpolys_single, <void*>self)
            return self.tmp
        self.tmp = {}
        smalljac_Lpolys(self.c, p, b, 0, callback_Lpolys, <void*>self)
        return self.tmp

    def group(self, unsigned long p, unsigned long b=0):
        raise NotImplementedError



