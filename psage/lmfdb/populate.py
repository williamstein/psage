"""
This module implements an object that is used to coordinate populating
the database.

This is an abstract base class for other classes, e.g., for populating
the database with newforms.
"""

class Populate:
    def __init__(self, collection):
        self.collection = collection
        
    def percent_done(self):
        return 100*float(self.count()) / self.collection.count()

    def populate_all(self, verbose=True):
        while True:
            if self.count() == self.collection.count():
                break
            d = self.percent_done()
            if verbose: print "Percent done: %.2f%%"%d
            self.populate_one(verbose=verbose)

