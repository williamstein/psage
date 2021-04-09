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

def populate_db(address, level_min, level_max, pmax=100,
                ncpus=sage.parallel.ncpus.ncpus()):
    """
    Compute and insert into the MongoDB database with the given
    address the Fourier coefficients a_p for p up to pmax for the optimal
    elliptic curves of the given range of levels (top level not
    included), using the given number of threads.

    Only curves with ap not yet set are affected by this function.
    """
    user, password = userpass()
    import math, random
    from sage.all import prime_range, parallel, pari

    level_min = int(level_min); level_max = int(level_max)
    P = prime_range(pmax)
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
                         'ap':{'$exists':False}}):
            E = pari('ellinit(%s,1)'%v['weq'])
            ap = dict([(str(p),int(E.ellap(p))) for p in P])
            C.update({'_id':v['_id']}, {'$set':{'ap':ap}})

    for ans in f(blocks):
        print(ans)


"""
EXAMPLE QUERIES:

from pymongo import Connection
db = Connection(port=int(29000)).research
e = db.ellcurves
v = e.find({'level':{'$lt':100r}, 'ap.2':-2r}, )
sage: v.count()
7

sage: v = e.find({'level':{'$lt':100r}, 'ap.2':-2r}, ['level', 'weq']) 
sage: v.next()
{u'weq': u'[0,-1,1,-10,-20]', u'_id': ObjectId('4c9258841e8b55611895b170'), u'level': 11}
sage: v.next()
{u'weq': u'[0,0,1,-1,0]', u'_id': ObjectId('4c9258841e8b55611895b1bc'), u'level': 37}

sage: v = e.find({'level':{'$lt':100r}, 'ap.2':{'$mod':[int(2),int(0)]}}, ['level', 'weq', 'ap.2', 'ap.3']) 
sage: v.next()
{u'ap': {u'3': -1, u'2': -2}, u'weq': u'[0,-1,1,-10,-20]', u'_id': ObjectId('4c9258841e8b55611895b170'), u'level': 11}
"""

