from sqrt5 import *

def ideals_of_bounded_norm(B):
    return sum([v for _, v in F.ideals_of_bdd_norm(B).iteritems()], [])
    
def ideals_of_norm(v):
    z = F.ideals_of_bdd_norm(max(v))
    return sum([z[n] for n in v],[])
