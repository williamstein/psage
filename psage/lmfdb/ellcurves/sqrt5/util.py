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

def ellcurves_sqrt5(address='localhost:29000', username=None, password=None):
    from sage.databases.cremona import cremona_letter_code
    from psage.modform.hilbert.sqrt5.sqrt5 import F
    from .aplists import labeled_primes_of_bounded_norm
    labels, primes = labeled_primes_of_bounded_norm(F, 100)

    from pymongo import Connection
    C = Connection(address).research
    if username is None or password is None:
        from psage.lmfdb.auth import userpass
        username, password = userpass()
        
    if not C.authenticate(username, password):
        raise RuntimeError("failed to authenticate")
    
    return C.ellcurves_sqrt5

def find_isogeneous_curves(ellcurves_sqrt5, E):
    """
    INPUT:
        - ellcurves_sqrt5 -- MongoDB collection
        - E -- an elliptic curve over Q(sqrt(5))
    OUTPUT:
        - cursor iterating over entries in the collection that have
          the same good a_p, for p of norm up to 100.
    """
    from .aplists import aplist
    w = aplist(E, 100)
    v = dict([('ap.%s'%p, a) for p, a in w.items()])
    from psage.modform.hilbert.sqrt5.tables import canonical_gen    
    v['level'] = str(canonical_gen(E.conductor())).replace(' ','')
    return ellcurves_sqrt5.find(v)

    
       
    
