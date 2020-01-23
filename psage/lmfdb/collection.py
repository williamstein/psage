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
This module defines a light wrapper around a MongoDB collection in the
MFDB database.
"""
from __future__ import print_function

from builtins import object
class Collection(object):
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
        print(cmd)
        os.system(cmd)

    def find(self, *args, **kwds):
        """Perform a query on the collection.  See the help for self.collection.find."""
        return self.collection.find(*args, **kwds)

