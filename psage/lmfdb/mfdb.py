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
This module implements the ModularFormsDataBase (MFDB) class.

The MFDB class represents a connection to an MongoDB modular forms
database running on a local or remote server.  It is instantianted
with a hostname and port, and just makes a MongoDB connection with
that port, then grabs a reference to the mfdb database there.

The newforms method returns a Python object that can be used to query
the database about classical GL2 newforms over QQ, and compute new
data about such newforms.

The backup method backs up the whole mfdb database.
"""

class MFDB:
    def __init__(self, host='localhost', port=29000):
        # Open conection to the MongoDB
        from pymongo import Connection
        self.connection = Connection(host, port)
        self.db = self.connection.mfdb
        self.port = port
        self.host = host
        from objectdb import ObjectDB
        self.objectdb = ObjectDB(self.db)

    def __repr__(self):
        return "Modular Forms Database\n%s"%self.connection

    def newforms(self):
        """Returns object that can be used for querying about GL2
        newforms over QQ and populating the database with them."""
        from newforms import NewformCollection
        return NewformCollection(self.db.newforms, self)

    def backup(self, outdir=None):
        """Dump the whole database to outdir.  If outdir is None,
        dumps to backup/year-month-day-hour-minute."""
        import os
        if outdir is None:
            import time
            outdir = os.path.join('backup',time.strftime('%Y%m%d-%H%M'))
        cmd = 'time mongodump -h %s:%s -d mfdb -o "%s"'%(
            self.host, self.port, outdir)
        print cmd
        os.system(cmd)

    

