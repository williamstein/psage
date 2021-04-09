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
This module implement a class that represents the collection of
newforms.  

The code in this module defines classes both for working with the
collection of all newforms in the database, and for populating this
collection with new data about newforms.   The NewformCollection
class is instantiated by the main MFDB database object.
The classes for populating the newforms table further or in
turn instantiated by the NewformCollection object.

The newform collection has documents that contain data about newforms.
It also has a subcollection 'counts', which records the number of
newforms with given level, weight, and character.  """
from __future__ import print_function
from __future__ import absolute_import
from .populate import Populate
from past.builtins import cmp
from builtins import str
from .collection import Collection

class NewformCollection(Collection):
    """
    """
    def __repr__(self):
        return "Collection of newforms"

    def count(self):
        """Return number of newforms in the newforms collection."""
        return self.find().count()

    def spaces(self, key=0):
        """Return sorted (by either level (key=0), weight (key=1), etc.)
        list of triples
               (level, weight, character.order(), count)
        for which all count newforms in the corresponding space are known."""
        key = int(key)
        C = self.collection.counts
        # Query C for the 4-tuples, as described in the docstring above.
        Q = [(x['level'], x['weight'], x['character']['order'], x['count']) for x in
             C.find({},['level','weight','character.order','count'])]
        Q.sort(key=lambda x: x[key])
        return Q

    def normalize(self, level, weight, character):
        """
        Return normalized level, weight, character, and a MongoDB
        document representing them.  Normalized means the level and
        weight are Python ints, and the character is a Sage Dirichlet
        character (in particular, it is not None).
        """
        from .converter import to_db
        level  = to_db(level)
        weight = to_db(weight)
        if character is None:
            from sage.all import trivial_character
            character = trivial_character(level)
        e = to_db(character)
        return level, weight, character, {'level':level, 'weight':weight, 'character':e}

    def populate_newform_eigenvalue_field(self):
        """Return object that organizes populating the eigenvalue field
        property of newform documents."""
        return PopulateNewformEigenvalueField(self)

    def populate_newforms(self):
        """Return object that organizes populating the newform
        documents.  Get this object if you want to add new newforms to
        the database."""
        return PopulateNewforms(self)




class PopulateNewforms(Populate):
    def __repr__(self):
        return "Populate Newforms"

    def count(self):
        """Return number of newforms in the database (with degree field set)."""
        return self.collection.newforms.find({'degree':{'$exists':True}}).count()
    
    def populate(self, level, weight=2, character=None, verbose=True):
        nc = self.collection # newform collection
        level, weight, character, D = nc.normalize(level, weight, character)
        if verbose: print(D)
        # Check the counts subcollection
        if nc.collection.counts.find(D).count() > 0:
            # Don't bother
            if verbose: print("Skipping since counts subcollection already has an entry.")
            return
        from psage.modform.rational.newforms import degrees
        degs = degrees(level, weight, character)
        for num, d in enumerate(degs):
            deg = int(d)
            # Update the document for the given newform
            query = dict(D)
            query['number'] = num
            nc.collection.update(query, {'$set':{'degree':deg}},
                                 upsert=True, safe=True)
        D['count'] = len(degs)
        nc.collection.counts.insert(D, safe=True)

    def populate_all_characters(self, level, weight, verbose=True):
        from sage.all import DirichletGroup
        G = DirichletGroup(level)
        B = G.galois_orbits()
        B = [character[0].minimize_base_ring() for character in B]
        for character in B:
            self.populate(level, weight, character, verbose)

    def populate_quadratic_characters(self, level, weight, verbose=True):
        from sage.all import DirichletGroup, QQ
        G = DirichletGroup(level,QQ)
        B = G.galois_orbits()
        B = [character[0].minimize_base_ring() for character in B
             if character[0].order()==2]
        for character in B:
            self.populate(level, weight, character, verbose)
                

class PopulateNewformEigenvalueField(Populate):
    def __repr__(self):
        return "Populating newform Hecke eigenvalue fields"
    
    def count(self):
        """Return number of newforms with eigenvalue field computed."""
        return self.collection.collection.find({'eigenvalue_field':{'$exists':True}}).count()

    def populate_one(self, verbose=True):
        """
        Compute Hecke eigenvalues for one unknown level,weight,character.
        If all data is known, raise a ValueError.
        """
        A = self.collection.collection.find_one({'eigenvalue_field':{'$exists':False}})
        if A is None:
            raise ValueError("All Hecke eigenvalue fields are currently known.")
        from .converter import db_converter
        self.populate(A['level'], A['weight'],
                      db_converter.to_dirichlet_character(A['character']),
                      verbose=verbose)
    
    def populate(self, level, weight, character=None, verbose=True):
        nc = self.collection
        level, weight, character, D = nc.normalize(level, weight, character)
        if verbose: print(D)
        C = nc.collection.counts.find(D)
        if C.count() == 0:
            if verbose: print("Skipping -- no newforms known yet (counts=0)")
            return
        cnt = C.next()['count']
        if cnt == 0:
            # There are no newforms, so don't do any further work.
            if verbose: print("No newforms")
            return
        # Now check to see if all of the eigenvalue fields are known
        # by doing a query for all forms of the given level, weight, and
        # character for which the eigenvalue_field key is set.
        E = dict(D)
        E['eigenvalue_field'] = {'$exists':True}
        if nc.collection.find(E).count() == cnt:
            if verbose: print("All eigenvalue fields already known")
            return 
        
        from psage.modform.rational.newforms import eigenvalue_fields
        from sage.all import set_random_seed, QQ

        set_random_seed(0)
        fields = eigenvalue_fields(level, weight, character)
        
        for num, K in enumerate(fields):
            # Update the document for the given newform
            query = dict(D)
            query['number'] = num
            if K == QQ:
                f = 'x-1'
            else:
                f = str(K.defining_polynomial()).replace(' ','')
            if verbose:
                print(f)
            nc.collection.update(query, {'$set':{'eigenvalue_field':f}}, safe=True)
