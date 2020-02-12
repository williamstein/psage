#################################################################################
#
# (c) Copyright 2010 William Stein
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


"""
This module implements a simple key:value store using SQLite3 and
cPickle, and other useful tools built on top of it.
"""

from future import standard_library
standard_library.install_aliases()
from builtins import next
from builtins import str
from builtins import range
from builtins import object
import sqlite3, zlib
try:
    import cPickle as pickle
except:
    import _pickle as pickle
# A key:value store

class SQLiteKeyValueStore(object):
    def __init__(self, file, compress=False):
        """
        Create or open the SQLite3-based key:value database stored in the given file.

        INPUTS:
            - file -- string; the name of a file.
            - compress -- bool (default: False); if True, by default compress all
              pickled values using zlib

        You do not have to be consistent with the compress option.   The database will still
        work if you switch back and forth between compress=True and compress=False.
        """
        self._db = sqlite3.connect(file)
        self._cursor = self._db.cursor()
        self._file = file
        self._compress = compress
        try:
            next(self._cursor.execute("select * from sqlite_master"))
        except StopIteration:
            # This exception will occur if the database is brand new (has no tables yet)
            try: 
                self._cursor.execute("CREATE TABLE cache (key BLOB, value BLOB, compressed INTEGER, UNIQUE(key))")
                self._cursor.execute("CREATE INDEX cache_idx ON cache(key)")
                self._db.commit()
            except sqlite3.OperationalError: 
                pass  # failure could happen if another process maybe created 
                      # and initialized the database at the same time.  That's fine.
            
    def __del__(self):
        """Called when the database is freed to close the connection."""        
        self._db.close()        
    
    def __repr__(self):
        """String representation of the database."""
        return "SQLite3-based key:value database stored in '%s'"%self._file
    
    def has_key(self, key):    
        """Returns True if database has the given key."""
        return self._cursor.execute("SELECT count(*) FROM cache WHERE key=?", (self._dumps(key),)).next()[0] > 0
            
    def __getitem__(self, key):
        """Return item in the database with given key, or raise KeyError."""        
        s = self._cursor.execute("SELECT value,compressed FROM cache WHERE key=?", (self._dumps(key),))
        try:
            v = next(s)
            return self._loads(str(v[0]), bool(v[1]))
        except StopIteration:
            raise KeyError(str(key))
        
    def __setitem__(self, key, value):
        """Sets an item in the database.  Call commit to make this permanent."""
        self._cursor.execute("INSERT OR REPLACE INTO cache VALUES(?, ?, ?)", (
            self._dumps(key), self._dumps(value, self._compress), self._compress))
        
    def __delitem__(self, key):
        """Removes an item from the database.  Call commit to make this permanent."""
        self._cursor.execute("DELETE FROM cache WHERE key=?", (self._dumps(key),))
        
    def _dumps(self, x, compress=False):
        """Converts a Python object to a binary string that can be stored in the database."""
        s = pickle.dumps(x,2)
        if compress:
            s = zlib.compress(s)
        return sqlite3.Binary(s)
        
    def _loads(self, x, compress=False):
        """Used internally to turn a pickled object in the database into a Python object."""
        if compress:
            x = zlib.decompress(x)
        return pickle.loads(x)
    
    def keys(self):
        """Return list of keys in the database."""
        return [self._loads(str(x[0])) for x in self._cursor.execute( "SELECT key FROM cache")]
        
    def commit(self):
        """Write assignments made to the database to disk."""
        self._db.commit()

        
def test_sqlite_keyval_1():
    """A straightforward test."""
    import tempfile
    file = tempfile.mktemp()
    try:
        for compress in [False, True]:
            db = SQLiteKeyValueStore(file, compress)

            db[2] = 3
            db[10] = {1:5, '17a':[2,5]}
            assert list(db.keys()) == [2,10]
            assert db[10] == {1:5, '17a':[2,5]}
            assert db[2] == 3
            db.commit()
            db[5] = 18  # does not get committed

            db = SQLiteKeyValueStore(file, not compress)
            assert list(db.keys()) == [2,10]
            assert db[10] == {1:5, '17a':[2,5]}
            assert db[2] == 3

            assert 2 in db
            assert 3 not in db
            del db
            import os; os.unlink(file)

    finally:
        if os.path.exists(file):
            import os; os.unlink(file)


# A SQLite cached function decorator

class sqlite_cached_function(object):
    """
    Use this like so::

         @sqlite_cached_function('/tmp/foo.sqlite', compress=True)
         def f(n,k=5):
             return n+k

    Then whenever you call f, the values are cached in the sqlite
    database /tmp/foo.sqlite.  This will persist across different
    sessions, of course.     Moreover, f.db is the underlying
    SQLiteKeyValueStore and f.keys() is a list of all keys computed
    so far (normalized by ArgumentFixer). 
    """
    def __init__(self, file, compress=False):
        self.db = SQLiteKeyValueStore(file, compress=compress)
        
    def __call__(self, f):
        """Return decorated version of f."""
        from sage.misc.function_mangling import ArgumentFixer
        A = ArgumentFixer(f)
        def g(*args, **kwds):
            k = A.fix_to_named(*args, **kwds)
            try:
                return self.db[k]
            except KeyError: pass    
            x = self.db[k] = f(*args, **kwds)  
            self.db.commit()
            return x            
        def keys():
            return list(self.db.keys())
        g.keys = keys    
        g.db = self.db
        return g

def test_sqlite_cached_function_1():
    try:
        import tempfile
        file = tempfile.mktemp()
        @sqlite_cached_function(file)
        def f(a, b=10):
            return a + b
        assert f(2) == 12
        assert f(2,4) == 6
        assert f(2) == 12
        assert f(2,4) == 6
    finally:
        import os; os.unlink(file)

def test_sqlite_cached_function_2():
    try:
        from sage.all import sleep, walltime
        import tempfile
        file = tempfile.mktemp()

        @sqlite_cached_function(file, compress=True)
        def f(a, b=10):
            sleep(1)
            return a + b
        f(2)
        f(2,b=4)

        t = walltime()
        assert f(2) == 12
        assert f(b=4,a=2) == 6
        assert walltime() - t < 1, "should be fast!"

        # Make new cached function, which will now use the disk cache first.
        @sqlite_cached_function(file, compress=True)
        def f(a, b=10):
            sleep(1)

        t = walltime()
        assert f(2) == 12
        assert f(b=4,a=2) == 6
        assert walltime() - t < 1, "should be fast!"

    finally:
        import os; os.unlink(file)

def test_sqlite_cached_function_3():
    import tempfile
    file = tempfile.mktemp()

    try:
        from sage.all import parallel, sleep

        # This "nasty" test causes 10 processes to be spawned all at once,
        # and simultaneously try to initialize and write to the database,
        # repeatedly.  This tests that we're dealing with concurrency robustly.

        @parallel(10)
        def f(a, b=10):
            @sqlite_cached_function(file)
            def g(a, b):
                sleep(.5)
                return a + b
            return g(a, b)

        for X in f(list(range(1,30))):
            assert X[1] == X[0][0][0] + 10

    finally:
        import os; os.unlink(file)

