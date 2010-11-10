"""

 (c) Copyright 2009-2010 Salman Baig and Chris Hall

 This file is part of ELLFF

 ELLFF is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 ELLFF is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""

#######################
# ELLFF DATABASE CODE #
#######################
import sage.databases.db
from sage.rings.all import PolynomialRing, GF

class jCurveDatabase(sage.databases.db.Database):
    def __init__(self, read_only=True):
        """
        Initialize the database.
        
        INPUT:
        
        
        -  ``read_only`` - bool (default: True), if True, then
           the database is read_only and changes cannot be committed to
           disk.
        """
        sage.databases.db.Database.__init__(self, name='jcurve_euler_tables', read_only=read_only)

    # may not need __getitem__
    """
    def __getitem__(self, key):
        If key=q is an integer, return all data about FF_q(t) in the database.
        If key=(q,n) is a pair of integers, return corresponding euler table,
        if it is in the database.
        
        INPUT:
        
        
        -  ``key`` - int or list of two ints
        
        
        OUTPUT: dict (if key is an int) or list (if key is a list)

        if isinstance(key, list) and len(key) > 1:
            return sage.databases.db.Database.__getitem__(self, key[0])[key[1]]

        return sage.databases.db.Database.__getitem__(self, key)
    """

    def __repr__(self):
        """
        String representation of this database. OUTPUT: str
        """
        return 'Database of euler tables for the versal j-curve'

    # TODO: which q are in db, which n are in db given q

    def _update_(self, q, n, table, force=False, verbose=False):
        r"""
        Upate the database for self over $\mathbb{F}_{q^n}$, forcing
        overwrite if so desired. If the table already exists and is
        forced to be overwritten, the two tables are compared for
        equality. If they are not equal, the old table is replaced
        with the new one.

        INPUT:

            - q -- an integer the size of F_q
            - n -- the degree of the extension of F_q
            - force -- boolean that forces overwrite
        
        """
        if self.read_only:
            raise RuntimeError, 'The database must not be read_only.'
        if self.has_key(q):
            if self[q].has_key(n) and force:
                if verbose:
                    print 'Already have this table; forcing overwrite'
                if not self[q][n] == table:
                    print 'Tables mismatch; replacing preexisting table with new given one'
                self[q][n] = table
                self.changed(q)
                # self[q][n] = self[q][n] # so database knows that self[q][n] changed
            else:
                self[q][n] = table
                self.changed(q)
                # self[q][n] = self[q][n] # so database knows that self[q][n] changed
        else:
            self[q] = {}
            self[q][n] = table
            self.changed(q)
            # self[q][n] = self[q][n] # so database knows that self[q][n] changed
        self.commit()

_jdb = None
def jCurveEulerTables(read_only=False):
    r"""
    Create the database of euler factors for the versal j curve.
    """
    
    global _jdb    
    if _jdb != None:
        return _jdb
    if _jdb == None:
        _jdb = jCurveDatabase(read_only)
        return _jdb

class LocalEulerDatabase(sage.databases.db.Database):
    def __init__(self, read_only=True):
        """
        Initialize the database.
        
        INPUT:
        
        
        -  ``read_only`` - bool (default: True), if True, then
           the database is read_only and changes cannot be committed to
           disk.
        """
        sage.databases.db.Database.__init__(self, name='local_euler_tables', read_only=read_only)

    def __repr__(self):
        """
        String representation of this database. OUTPUT: str
        """
        return 'Database of euler tables for user curves'

    # TODO: which q are in db, which n are in db given q, which ainvs are in db

    def _update_(self, ainvs, q, n, table, force=False, verbose=False):
        r"""
        Upate the database for self over $\mathbb{F}_{q^n}$, forcing
        overwrite if so desired. If the table already exists and is
        forced to be overwritten, the two tables are compared for
        equality. If they are not equal, the old table is replaced
        with the new one.

        INPUT:

            - q -- an integer the size of F_q
            - n -- the degree of the extension of F_q
            - force -- boolean that forces overwrite
        
        """
        if self.read_only:
            raise RuntimeError, 'The database must not be read_only.'
        if self.has_key(ainvs):
            if self[ainvs].has_key(q):
                if self[ainvs][q].has_key(n):
                    if verbose:
                        print 'Already have this table;',
                    if force:
                        if verbose:
                            print 'forcing overwrite'
                        if not self[ainvs][q][n] == table:
                            print 'Tables mismath; replacing preexisting table with new given one'
                        self[ainvs][q][n] = table
                        self.changed(ainvs)
                else:
                    self[ainvs][q][n] = table
                    self.changed(ainvs)
            else:
                self[ainvs][q] = {}
                self[ainvs][q][n] = table
                self.changed(ainvs)
        else:
            self[ainvs] = {}
            self[ainvs][q] = {}
            self[ainvs][q][n] = table
            self.changed(ainvs)
        self.commit()

_ldb = None
def LocalEulerTables(read_only=False):
    r"""
    Create the database of euler factors for the `user` curve
    """
    
    global _ldb    
    if _ldb != None:
        return _ldb
    if _ldb == None:
        _ldb = LocalEulerDatabase(read_only)
        return _ldb                       

def _save_euler_table(self, n, verbose=False):
    r"""
    Save the euler table for self over the degree n extension of
    $\mathbb{F}_q$ to disk. This is currently implemented with
    sage.database.db, which uses ZODB. If self is the versal j-curve,
    it stores the table in the database
    
    SAGE_ROOT/data/jcurve_euler_tables .

    Otherwise, the tables are stored in the `user` table

    SAGE_ROOT/data/local_euler_tables .

    It currently doesn't check if the table already is stored; it
    merely writes over it in that case. This should eventually be
    implemented using MongoDB.

    INPUT:
    
        - n -- the degree of the extension of F_q

    EXAMPLES::

        sage: import psage
        sage: K.<t> = psage.FunctionField(GF(11))
        sage: E = psage.ellff_EllipticCurve(K,[0,0,0,-27*t/(t-1728),54*t/(t-1728)])
        sage: E._build_euler_table(1)
        sage: E._euler_table(1)
        [0, 0, 4, -6, 3, 5, 1, -2, 4, -2, 3, 1]
        sage: E._build_euler_table(2)
        sage: E._build_euler_table(3)
        sage: E._euler_table(1)
        [0, 0, 4, -6, 3, 5, 1, -2, 4, -2, 3, 1]
        sage: E._save_euler_table(1)
        sage: E._save_euler_table(2)
        sage: E._save_euler_table(3)
        
    """
        
    import os
    SAGE_ROOT = os.environ['SAGE_ROOT']
    
    K = self.K
    R = self.R
    t = K.gens()[0]
    p = self.p
    d = self.d
    q = self.q
    R2 = PolynomialRing(GF(q), 's')
    a1n = R2(0)
    a1d = R2(1)
    a2n = R2(0)
    a2d = R2(1)
    a3n = R2(0)
    a3d = R2(1)
    a4n = R2(self.a4.numerator().coeffs())
    a4d = R2(self.a4.denominator().coeffs())
    a6n = R2(self.a6.numerator().coeffs())
    a6d = R2(self.a6.denominator().coeffs())

    ainvs = [0, 0, 0, self.a4, self.a6]
    ainvs_pairs = ((a1n, a1d), (a2n, a2d), (a3n, a3d), (a4n, a4d), (a6n, a6d))

    # recognize if self is j-curve and use special repository
    if ainvs == [0,0,0,-27*t*(t-1728)**3,54*t*(t-1728)**5]:
        if verbose:
            print 'j-curve recognized; saving euler table to database'
        if not os.path.exists(SAGE_ROOT + '/data/jcurve_euler_tables/jcurve_euler_tables'):
            print 'Database does not exist; starting a new one'
            if not os.path.exists(SAGE_ROOT + '/data/jcurve_euler_tables'):
                os.makedirs(SAGE_ROOT + '/data/jcurve_euler_tables/')
            filedb = open(SAGE_ROOT + '/data/jcurve_euler_tables/jcurve_euler_tables', "wb")
            filedb.close()
                
        euler_db = jCurveEulerTables(read_only = False)
        euler_db._update_(q, n, self._euler_table(n))
        if verbose:
            print euler_db.as_dict()
        euler_db.commit()

            
        # work with user's repository of euler tables
    else:
        if not os.path.exists(SAGE_ROOT + '/data/local_euler_tables/local_euler_tables'):
            print 'Database does not exist; creating a new one'
            if not os.path.exists(SAGE_ROOT + '/data/local_euler_tables'):
                os.makedirs(SAGE_ROOT + '/data/jcurve_euler_tables/')
            filedb = open(SAGE_ROOT + '/data/local_euler_tables/local_euler_tables', "wb")
            filedb.close()
                
        local_euler_db = LocalEulerTables(read_only = False)
        local_euler_db._update_(ainvs_pairs, q, n, self._euler_table(n))
        if verbose:
            print local_euler_db.as_dict()
        local_euler_db.commit()

def _load_euler_table(self, n, force=False, verbose=False):
    r"""
    Load the euler table for self over the degree n extension of
    $\mathbb{F}_q$ to disk. If self is the versal j-curve, the table
    is pulled from
    
    SAGE_ROOT/data/jcurve_euler_tables .
    
    Otherwise, the table is pulled from the `user` table

    SAGE_ROOT/data/local_euler_tables .

    This should eventually be implemented using MongoDB.

    It currently doesn't check if the key exist. If the key doesn't
    exist, a RuntimeError is raised by sage.database.db. This
    RuntimeError should be sufficient, so key checking may not be
    necessary.

    INPUT:

        - n -- the degree of the extension of F_q
        - force -- boolean that overwrites self's euler table with
                   one from database

    EXAMPLES::

        sage: import psage
        sage: K.<t> = psage.FunctionField(GF(11))
        sage: E = psage.ellff_EllipticCurve(K,[0,0,0,-27*t/(t-1728),54*t/(t-1728)])
        sage: E._euler_table(1)
        Traceback (most recent call last):
            ...
        RuntimeError: table is empty
        sage: E._load_euler_table(1)
        sage: E._euler_table(1)
        [0, 0, 4, -6, 3, 5, 1, -2, 4, -2, 3, 1]
            
    """
        
    import os
    SAGE_ROOT = os.environ['SAGE_ROOT']
        
    K = self.K
    R = self.R
    t = K.gens()[0]
    p = self.p
    d = self.d
    q = self.q
    R2 = PolynomialRing(GF(q), 's')
    s = R2.gens()[0]
    a1n = R2(0)
    a1d = R2(1)
    a2n = R2(0)
    a2d = R2(1)
    a3n = R2(0)
    a3d = R2(1)
    a4n = R2(self.a4.numerator().coeffs())
    a4d = R2(self.a4.denominator().coeffs())
    a6n = R2(self.a6.numerator().coeffs())
    a6d = R2(self.a6.denominator().coeffs())

    ainvs = [0, 0, 0, self.a4, self.a6]
    ainvs_pairs = ((a1n, a1d), (a2n, a2d), (a3n, a3d), (a4n, a4d), (a6n, a6d))

    # recognize if self is j-curve and use special repository
    if ainvs == [0,0,0,-27*t*(t-1728)**3,54*t*(t-1728)**5]:
        if verbose:
            print 'j-curve recognized; saving euler table to database'
        if not os.path.exists(SAGE_ROOT + '/data/jcurve_euler_tables/jcurve_euler_tables'):
            print 'Database does not exist; cannot load from it'
        else:
            euler_db = jCurveEulerTables()
            # check that keys exist?
            self._set_euler_table(n, euler_db[q][n], force)
            
    # work with user's repository of euler tables
    else:
        if not os.path.exists(SAGE_ROOT + '/data/local_euler_tables/local_euler_tables'):
            print 'Database does not exist; cannot load from it'
        else:                
            local_euler_db = LocalEulerTables()
            # check that keys exist?
            self._set_euler_table(n, local_euler_db[ainvs_pairs][q][n], force)
