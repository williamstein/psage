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

def str_to_apdict(s, labels):
    return dict([(labels[a], int(b)) for a,b in enumerate(s.split()) if b != '?'])

def labeled_primes_of_bounded_norm(F, B):
    """
    Return a list [prime][letter code] of strings, corresponding to
    the primes of F of bounded norm of F, ordered by residue
    characteristic (not norm).
    """
    from sage.databases.cremona import cremona_letter_code
    from psage.modform.hilbert.sqrt5.sqrt5 import primes_of_bounded_norm
    labels = []
    last_c = 1
    number = 0
    primes = primes_of_bounded_norm(F, B)
    for p in primes:
        c = p.smallest_integer() # residue characteristic
        if c != last_c:
            last_c = c
            number = 0
        else:
            number += 1
        labels.append('%s%s'%(c,cremona_letter_code(number)))
    return labels, primes


def import_table(address, aplists_txt_filename, max_level=None):
    """
    Import a text table of eigenvalues, using upsert to avoid
    replication of data.
    """
    from psage.lmfdb.auth import userpass
    user, password = userpass()

    from sage.databases.cremona import cremona_letter_code
    from psage.modform.hilbert.sqrt5.sqrt5 import F
    labels, primes = labeled_primes_of_bounded_norm(F, 100)

    from pymongo import Connection
    C = Connection(address).research
    if not C.authenticate(user, password):
        raise RuntimeError("failed to authenticate")
    e = C.ellcurves_sqrt5
    for X in open(aplists_txt_filename).read().splitlines():
        if X.startswith('#'):
            continue
        Nlevel, level, iso_class, ap = X.split('\t')
        ap = str_to_apdict(ap, labels)        
        Nlevel = int(Nlevel)
        iso_class = cremona_letter_code(int(iso_class))
        v = {'level':level, 'iso_class':iso_class,
             'number':1, 'Nlevel':Nlevel, 'ap':ap}
        if max_level and Nlevel > max_level: break
        print(v)
        spec = dict(v)
        del spec['ap']
        e.update(spec, v, upsert=True, safe=True)
    return e
        



def aplist(E, B=100):
    """
    Compute aplist for an elliptic curve E over Q(sqrt(5)), as a
    string->number dictionary.

    INPUT:
        - E -- an elliptic curve
        - B -- a positive integer (default: 100)

    OUTPUT:
        - dictionary mapping strings (labeled primes) to Python ints,
          with keys the primes of P with norm up to B such that the
          norm of the conductor is coprime to the characteristic of P.
    """
    from psage.modform.hilbert.sqrt5.tables import canonical_gen    
    v = {}
    from psage.modform.hilbert.sqrt5.sqrt5 import F
    labels, primes = labeled_primes_of_bounded_norm(F, B)

    from sage.all import ZZ
    N = E.conductor()
    try:
        N = ZZ(N.norm())
    except:
        N = ZZ(N)
    
    for i in range(len(primes)):
        p = primes[i]
        k = p.residue_field()
        if N.gcd(k.cardinality()) == 1:
            v[labels[i]] = int(k.cardinality() + 1 - E.change_ring(k).cardinality())
    return v


