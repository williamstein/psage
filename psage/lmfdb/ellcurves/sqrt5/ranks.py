from __future__ import print_function
from __future__ import absolute_import
#################################################################################
#
# (c) Copyright 2011 William Stein
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

from sage.all import EllipticCurve

from psage.modform.hilbert.sqrt5.sqrt5 import F

def r_an(ainvs, eps=1e-4):
    E = EllipticCurve(F, ainvs)
    L = E.lseries()
    D = L.dokchitser(30)
    f = D.taylor_series(1)
    r = 0
    while abs(f[r]) < eps:
        r += 1
    if D.eps == 1:
        assert r%2 == 0
    else:
        assert r%2 == 1
    return r
    
def r_alg(a_invs):
    return EllipticCurve(F, a_invs).rank()    

from sage.all import fork
@fork
def rank(a_invs):
    try:
        return r_alg(a_invs)
    except ValueError:
        return r_an(a_invs)
    
def compute_conjectural_ranks(level_norms, address='localhost:29000'):
    """
    For each level norm in the input list, compute the
    *conjectural* rank of the elliptic curve, and put
    that data in a field "r?" in the database.
    This uses Simon 2-descent if it works, and otherwise
    uses Dokchitser.  This should not be considered
    super-reliable!
    """
    from sage.all import sage_eval

    from . import util
    C = util.ellcurves_sqrt5(address)
    for N in level_norms:
        for E in C.find({'Nlevel':int(N), 'r?':{'$exists':False}, 'weq':{'$exists':True}}):
            print(E)
            weq = sage_eval(E['weq'], {'a':F.gen()})
            E['r?'] = int(rank(weq))
            print(E)
            C.update({'_id':E['_id']}, E, safe=True)
            
            

    
