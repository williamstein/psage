"""
Queries

This file contains code that carries out interesting queries
regarding the database of elliptic curves.
"""

def create_counts_table(levels, address):
    """
    QUESTION: What proportion of curves in the database have
    squarefree conductor, as a function of the conductor?

    To answer, make another table ellcurves.counts with documents:

       {'_id':N,  'c':number_of_isogeny_classes_of_curves_of_conductor_N, 'ss':True}

    where 'ss':True is set only if N is squarefree.

    Once we have this table, the rest should be relatively easy.
    """
    from psage.lmfdb.auth import userpass
    user, password = userpass()
    from pymongo import Connection
    from pymongo.connection import DuplicateKeyError
    C = Connection(address).research
    C.authenticate(user, password)
    ellcurves = C.ellcurves
    db_counts = ellcurves.counts
    
    from sage.all import is_squarefree
    
    def counts(N):
        N = int(N)
        c = ellcurves.find({'level':N, 'number':1}).count()
        doc = {'_id':N, 'c':c}
        if is_squarefree(N):
            doc['ss'] = True
        try:
            db_counts.insert(doc, safe=True)
        except DuplicateKeyError:
            print '[%s]'%N,
        else:
            print N,
        import sys; sys.stdout.flush()
        
    for N in levels:
        counts(N)


    


