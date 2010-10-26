"""
This module defines a light wrapper around a MongoDB collection in the
MFDB database.
"""

class Collection:
    def __init__(self, collection, db):
        self.collection = collection
        self.db = db

    def backup(self, outdir=None):
        """Dump this collection to outdir.  If outdir is None,
        dumps to backup/year-month-day-hour-minute."""
        import os
        if outdir is None:
            import time
            outdir = os.path.join('backup',time.strftime('%Y%m%d-%H%M'))
        cmd = 'time mongodump -c "%s" -h %s:%s -d mfdb -o "%s"'%(
            self.collection.name, self.db.host, self.db.port, outdir)
        print cmd
        os.system(cmd)

    def find(self, *args, **kwds):
        """Perform a query on the collection.  See the help for self.collection.find."""
        return self.collection.find(*args, **kwds)

