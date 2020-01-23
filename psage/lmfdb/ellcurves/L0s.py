from __future__ import print_function
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

from builtins import str
from builtins import range
import sage.parallel.ncpus
from psage.lmfdb.auth import userpass

def populate_db(address, level_min, level_max, num_zeros=100,
                ncpus=sage.parallel.ncpus.ncpus()):
    """
    Compute and insert into the MongoDB database with the given
    address the imaginary parts of the first num_zeros zeros of the
    L-series of each optimal elliptic curve in the given level range.

    Only curves with L0s not yet set are affected by this function.
    The key on the database is "L0s".
    """
    user, password = userpass()

    import math, random
    from sage.all import parallel, EllipticCurve

    from sage.libs.lcalc.lcalc_Lfunction import Lfunction_from_elliptic_curve
    
    level_min = int(level_min); level_max = int(level_max)
    s = int(math.ceil((level_max - level_min)/float(ncpus)))
    blocks = [(level_min+i*s, min(level_max,level_min+(i+1)*s)) for i in range(ncpus)]
    
    @parallel(ncpus)
    def f(l_min, l_max):
        from pymongo import Connection
        C = Connection(address).research
        C.authenticate(user, password)
        C = C.ellcurves
        for v in C.find({'level':{'$gte':level_min, '$lt':level_max},
                         'number':1,
                         'L0s':{'$exists':False}}):
            L = Lfunction_from_elliptic_curve(EllipticCurve(eval(v['weq'])), 10**5)
            z = L.find_zeros_via_N(num_zeros)
            L0s =  dict([(str(i),float(z[i])) for i in range(len(z))])
            C.update({'_id':v['_id']}, {'$set':{'L0s':L0s}})

    for ans in f(blocks):
        print(ans)


"""
EXAMPLE QUERIES:

from pymongo import Connection
db = Connection(port=int(29000)).research
e = db.ellcurves
v = e.find({'level':{'$lt':100r}, 'L0s':{'$exists':True}}, )
v.count()
7

This counts the number of optimal curves for which the 0-th zero (the
first one) is >1 and the number for which it is < 1.

sage: v = e.find({'level':{'$lt':100r}, 'L0s.0':{'$gt':int(1)}}, )
sage: v.count()
75
sage: v = e.find({'level':{'$lt':100r}, 'L0s.0':{'$lt':int(1)}}, )
sage: v.count()
17

"""

