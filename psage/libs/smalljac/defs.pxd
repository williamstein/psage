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
    cdef int SMALLJAC_A1_ONLY

    int smalljac_Lpoly (long a[], char *curve, unsigned long q, unsigned long flags)
