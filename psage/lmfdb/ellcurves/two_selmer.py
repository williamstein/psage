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


from sage.all import mwrank_EllipticCurve

def selmer2(a_invariants, max_time=None):
    if max_time is None:
        E = mwrank_EllipticCurve(a_invariants)
        return int(E.selmer_rank())
    else:
        try:
            from sage.all import fork  # not in official sage.
        except ImportError:
            raise ImportError, "You need to apply the patch from http://trac.sagemath.org/sage_trac/ticket/9631"
        @fork(timeout=max_time)
        def f():
            E = mwrank_EllipticCurve(a_invariants)
            return int(E.selmer_rank())
        return f()

import sage.parallel.ncpus

def populate_db(address, level_min, level_max, 
                ncpus=sage.parallel.ncpus.ncpus(),
                max_time=None):
    """
    Compute and insert into the database the 2-selmer ranks of all the
    curves in a give range of levels, for which 2-selmer ranks aren't
    already known.
    """
    import math, random
    from sage.all import prime_range, parallel, pari

    level_min = int(level_min); level_max = int(level_max)
    s = int(math.ceil((level_max - level_min)/float(ncpus)))
    blocks = [(level_min+i*s, min(level_max,level_min+(i+1)*s)) for i in range(ncpus)]
    
    @parallel(ncpus)
    def f(l_min, l_max):
        from pymongo import Connection
        C = Connection(address).research.ellcurves
        for v in C.find({'level':{'$gte':level_min, '$lt':level_max},
                         'sel2':{'$exists':False}}):
            sel2 = selmer2(eval(v['weq']), max_time)
            C.update({'_id':v['_id']}, {'$set':{'sel2':sel2}})

    for ans in f(blocks):
        print ans
    
"""
EXAMPLE QUERIES:

from pymongo import Connection
db = Connection(port=int(29000)).research
e = db.ellcurves
v = e.find({'level':{'$lt':100r}, 'sel2':{'$exists':True}})
sage: v.count()
7
"""
