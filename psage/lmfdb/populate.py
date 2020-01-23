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
This module implements an object that is used to coordinate populating
the database.

This is an abstract base class for other classes, e.g., for populating
the database with newforms.
"""
from __future__ import print_function

from builtins import object
class Populate(object):
    def __init__(self, collection):
        self.collection = collection
        
    def percent_done(self):
        return 100*float(self.count()) / self.collection.count()

    def populate_all(self, verbose=True):
        while True:
            if self.count() == self.collection.count():
                break
            d = self.percent_done()
            if verbose: print("Percent done: {0:.2f}%".format(d))
            self.populate_one(verbose=verbose)

