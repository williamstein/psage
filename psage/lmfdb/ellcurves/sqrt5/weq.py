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


def import_table(address, table_filename, max_level=None):
    """
    Import a text table of weq's, using upsert to avoid
    replication of data.  Format is like this:

    ::
    
      31      5*a-2       0   -3 2 -2 2 4 -4 4 -4 -2 -2 ? ? -6 -6 12 -4 6 -2 -8 0 0 16 10 -6              [1,a+1,a,a,0]
      31      5*a-3       0   -3 2 -2 2 -4 4 -4 4 -2 -2 ? ? -6 -6 -4 12 -2 6 0 -8 16 0 -6 10              [1,-a-1,a,0,0]
      36      6           0   ? ? -4 10 2 2 0 0 0 0 -8 -8 2 2 -10 -10 2 2 12 12 0 0 10 10                 [a,a-1,a,-1,-a+1]
    """
    from psage.modform.hilbert.sqrt5.sqrt5 import F
    from sage.databases.cremona import cremona_letter_code
    from aplists import labeled_primes_of_bounded_norm, str_to_apdict
    primes = labeled_primes_of_bounded_norm(F, 100)

    from psage.lmfdb.auth import userpass
    user, password = userpass()

    from pymongo import Connection
    C = Connection(address).research
    if not C.authenticate(user, password):
        raise RuntimeError, "failed to authenticate"
    e = C.ellcurves_sqrt5


    for X in open(table_filename).read().splitlines():
        if X.startswith('#'):
            continue
        z = X.split()
        Nlevel = z[0]; level = z[1]; iso_class = z[2]; weq = z[-1]
        ap = ' '.join(z[3:-1])
        ap = str_to_apdict(ap, primes)
        Nlevel = int(Nlevel)
        iso_class = cremona_letter_code(int(iso_class))
        v = {'level':level, 'iso_class':iso_class,
             'number':1, 'Nlevel':Nlevel, 'ap':ap,
             'weq':weq}
        if max_level and Nlevel > max_level: break
        print v
        spec = dict(v)
        del spec['weq']
        e.update(spec, v, upsert=True, safe=True)
 
        

