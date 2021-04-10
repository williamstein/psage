#################################################################################
#
# (c) Copyright 2011 William Stein
#
#  This file is part of PSAGE.
#
#  PSAGE is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
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

from cysignals.signals cimport sig_on,sig_off

from sage.finance.time_series cimport TimeSeries
from sage.stats.intlist cimport IntList

cdef extern:
    double log(double)
    double sqrt(double)
    double cos(double)
    double sin(double)

cdef class J:
    cdef long N_max
    cdef TimeSeries a, s
    cdef IntList pv
    
    def __init__(self, int N_max=10**5):
        from sage.all import prime_powers, factor
        self.N_max   = N_max
        PP = prime_powers(N_max+1)[1:]

        cdef double logp
        cdef int i, n  = len(PP)
        self.a   = TimeSeries(n)
        self.s   = TimeSeries(n)
        self.pv  = IntList(n)
        i = 0
        for pv in PP:
            F = factor(pv)
            p, v = F[0]
            self.pv._values[i] = pv
            logp = log(p)
            self.a._values[i] = logp/sqrt(pv)
            self.s._values[i] = v*logp
            i += 1

    def __repr__(self):
        return "The function J(t,N) of the Mazur-Stein game with N_Max=%s"%self.N_max

    cpdef double F(self, double t, int N):
        cdef double ans = 0
        cdef int i = 0
        if N > self.N_max:
            raise ValueError("J not defined for N > %s"%self.N_max)
        while 1:
            if self.pv._values[i] >= N:
                return ans
            ans += self.a._values[i] * cos(self.s._values[i]*t)
            i += 1
            
    cpdef double e(self, double t, int N):
        return sqrt(N)/(.25 + t*t) * \
               ( cos(t*log(N)) + 2*t*sin(t*log(N)) )

    cpdef double H(self, double t, int N):
        cdef double ans = 0, F = 0
        cdef int i = 0, n
        if N > self.N_max:
            raise ValueError("J not defined for N > %s"%self.N_max)
        n = 1
        sig_on()
        while 1:
            if self.pv._values[i] > N:
                # time to return.  But first we have to add a few more
                # values to ans up to N (where F is constant).
                while n <= N:
                    ans += F + self.e(t,n)
                    #print (n, F, self.e(t,n))
                    n += 1
                # Finally return after dividing by N.
                sig_off()                
                return ans / N
            # At this point, F is equal to self.F(t, self.pv._values[i]),
            # and now we add to our running total a range of values of

            #    F(t, n) + e(t,n)
            while n <= self.pv._values[i]:
                ans += F + self.e(t,n)
                #print (n, F, self.e(t,n))                
                n += 1

            F += self.a._values[i] * cos(self.s._values[i]*t)

            # move to next prime power
            i += 1

    cpdef double J(self, double t, int N):
        return self.H(t,N) / log(N)

    def __call__(self, double t, int N):
        return self.J(t,N)

    cpdef double H2(self, double t, int N):
        cdef int n
        return sum([self.F(t,n) + self.e(t,n) for n in range(1,N+1)]) / N

    cpdef double J2(self, double t, int N):
        return self.H2(t, N) / log(N)

    ######################################################
    # Cesaro of Cesaro...
    ######################################################    
    def Jc(self, double t, int N):
        cdef double ans = 0, F = 0, ans2 = 0
        cdef int i = 0, n
        if N > self.N_max:
            raise ValueError("J not defined for N > %s"%self.N_max)
        n = 1
        sig_on()
        while 1:
            if self.pv._values[i] > N:
                # time to return.  But first we have to add a few more
                # values to ans up to N (where F is constant).
                while n <= N:
                    ans += F + self.e(t,n)
                    if n >= 2:
                        ans2 += ans / (n*log(n))                    
                    #print (n, F, self.e(t,n))
                    n += 1
                # Finally return after dividing by N.
                sig_off()                
                return ans2 / N
            # At this point, F is equal to self.F(t, self.pv._values[i]),
            # and now we add to our running total a range of values of

            #    F(t, n) + e(t,n)
            while n <= self.pv._values[i]:
                ans += F + self.e(t,n)
                if n >= 2:
                    ans2 += ans / (n*log(n))
                #print (n, F, self.e(t,n))                
                n += 1

            F += self.a._values[i] * cos(self.s._values[i]*t)

            # move to next prime power
            i += 1

def cesaro_sum(TimeSeries x, TimeSeries y):
    """
    Given the graph of a function defined by (x[i], y[i]), compute its
    Cesaro smoothing.  This function returns a new series giving the y
    values of the Cesaro smoothing evaluated at the same values of x.
    """
    assert len(x) == len(y)
    cdef Py_ssize_t n
    cdef TimeSeries yc = TimeSeries(len(y))
    yc._values[0] = 0
    for n in range(1, len(x)):
        yc._values[n] = (x._values[n] - x._values[n-1])*y._values[n] + yc._values[n-1]
        
    for n in range(1,len(x)):
        yc._values[n] /= (x._values[n] - x._values[0])
    return yc

def average_value_function(v):
    """Given a list v = [...,(x,y),...] of points with the x's
    increasing that defines a piecewise linear function y=f(x),
    compute and return the average value function (1/z)*int_{v[0][0]}^z f(x)dx,
    expressed again as a list with the same x's. 
    """
    raise NotImplementedError

def plot_with_averages(v):
    """Given a list v = [...,(x,y),...] of points with the x's increasing, plot the function
    through these points, along with the graph in red of the average value of this function
    """
    from sage.all import line
    a = average_value_function(v)
    return line(v) + line(a, color='red')
    

###############################################################################
# J with a character
###############################################################################

cdef class J_chi(J):
    cdef bint _e
    def __init__(self, chi, int N_max=10**5, e=False):
        self._e = e
        # This is exactly the init from J, but the
        # definition of self.a._values[i] is changed below.
        from sage.all import prime_powers, factor
        self.N_max   = N_max
        PP = prime_powers(N_max+1)[1:]

        cdef double logp
        cdef int i, n  = len(PP)
        self.a   = TimeSeries(n)
        self.s   = TimeSeries(n)
        self.pv  = IntList(n)
        i = 0
        for pv in PP:
            F = factor(pv)
            p, v = F[0]
            self.pv._values[i] = pv
            logp = log(p)
            self.a._values[i] = logp/sqrt(pv) * chi(pv)
            self.s._values[i] = v*logp
            i += 1

    cpdef double e(self, double t, int N):
        if self._e:
            return sqrt(N)/(.25 + t*t) * \
                   ( cos(t*log(N)) + 2*t*sin(t*log(N)) )
        else:
            return 0

