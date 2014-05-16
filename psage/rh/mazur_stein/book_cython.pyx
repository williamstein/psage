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

include "stdsage.pxi"
include "interrupt.pxi"

def mult_parities(int bound, bint verbose=False):
    """
    Return list v of length bound such that v[n] is 0 if n is
    multiplicative even, and v[n] is 1 if n is multiplicatively odd.
    Also v[0] is None.

    This goes up to bound=`10^7` in about 2 seconds, but uses a lot of
    memory beyond this size.  On a machine with sufficient `10^8`
    takes around 20 seconds. If you have about 20-30GB of RAM, you can
    go to `10^9` in about 6 minutes; this is far enough to witness the
    first crossover point.
    """
    from sage.all import log, floor, prime_range
    cdef int i, j, k, n, p, last_len, cur_ptr, last_parity, cur_parity
    cdef long long m

    cdef int* v = <int*> sage_malloc(sizeof(int)*bound)
    if not v: raise MemoryError

    cdef int* last = <int*> sage_malloc(sizeof(int)*bound)
    if not last: raise MemoryError
    
    cdef int* cur = <int*> sage_malloc(sizeof(int)*bound)
    if not cur: raise MemoryError    

    for i in range(bound):
        v[i] = -1

    v[1] = 0
    
    P = prime_range(bound)
    cdef int* primes = <int*> sage_malloc(sizeof(int)*len(P))
    if not primes: raise MemoryError
    for i in range(len(P)):
        primes[i] = P[i]
        v[primes[i]] = 1
        last[i] = primes[i]

    cdef int len_P = len(P)
    last_len = len_P
    last_parity = 1
    cdef int loops = floor(log(bound,2))+1
    for k in range(loops):
        cur_ptr = 0
        cur_parity = (last_parity+1)%2
        if verbose:
            print "loop %s (of %s);  last = %s"%(k,loops, last_len)
        sig_on()
        for n in range(last_len):
            for j in range(len_P):
                m = (<long long> last[n]) * (<long long> primes[j])
                if m >= bound:
                    break
                if v[m] == -1:
                    v[m] = cur_parity
                    cur[cur_ptr] = m
                    cur_ptr+=1
        sig_off()
        last_parity = cur_parity
        last_len = cur_ptr
        sig_on()
        for n in range(last_len):
            last[n] = cur[n]
        sig_off()

    ans = [v[i] for i in range(bound)]
    sage_free(v)
    sage_free(last)
    sage_free(cur)
    sage_free(primes)
    ans[0] = None
    return ans
