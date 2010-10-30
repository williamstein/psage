# -*- coding: utf-8 -*-
r"""
Algorithms for use together with MySubgroup -- an extension to the standard implementation of subgroups of the modular group.

AUTHOR:

 - Fredrik Stroemberg

 


"""
#*****************************************************************************
#  Copyright (C) 2010 Fredrik Strömberg <stroemberg@mathematik.tu-darmstadt.de>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#****************************************************************************

include "sage/ext/interrupt.pxi"  # ctrl-c interrupt block support
include "sage/ext/stdsage.pxi"  # ctrl-c interrupt block support
include "sage/ext/cdefs.pxi"

from copy import deepcopy
from sage.combinat.permutation import Permutation_class
from sage.groups.perm_gps.permgroup_element import PermutationGroupElement

def are_transitive_permutations(E,R):
    r""" Check that E,R are transitive permutations, i.e. that <E,R>=S_N

    INPUT:

         - ``E`` -- permutation on N letters 
         - ``R`` -- permutation on N letters

             - E and R can be in any of the following formats:

                 - list [a1,a2,...,aN]
                 - member of Permutations(N)
                 - member of SymmetricGroup(N)

     OUTPUT:

     - bool  


     EXAMPLES::

         sage: E=Permutations(4)([1,2,4,3]); E.to_cycles()
         [(1,), (2,), (3, 4)]
         sage: R=Permutations(4)([2,1,3,4]); R.to_cycles()
         [(1, 2), (3,), (4,)]
         sage: are_transitive_permutations(E,R)
         False
         sage: R=Permutations(4)([2,3,1,4]); R.to_cycles()
         [(1, 2, 3), (4,)]
         sage: are_transitive_permutations(E,R)
         True
         sage: ES=SymmetricGroup(4)([1,2,4,3]);ES
         (3,4)
         sage: ER=SymmetricGroup(4)([2,3,1,4]);ER
         (1,2,3)
         sage: are_transitive_permutations(ES,RS)
         True

    """
    gotten=list()    
    t0=isinstance(E,list)
    if(t0):
        El=E; Rl=R
    else:
        t1=isinstance(E,Permutation_class) # constructed from Permutations(N)
        if(t1):
            El=list(E)
            Rl=list(R)
        else:
            t2=isinstance(E,PermutationGroupElement) # constructed from SymmetricGroup(N)
            if(t2):
                El=E.list()
                Rl=R.list()           
            else:
                raise TypeError, "Indata need to be Permuatations! Got E=%s of type=%s" %(E,type(E))
    N=len(El)
    gotten.append(Rl[0])
    gotten0=deepcopy(gotten)
    for j in range(N):
        for x in gotten0:
            y=El[x-1]
            if(gotten.count(y)==0):
                gotten.append(y)
            yy=Rl[y-1]
            if(gotten.count(yy)==0):
                gotten.append(yy)
            yyy=Rl[yy-1]
            if(gotten.count(yyy)==0):
                gotten.append(yyy)
        if(gotten0==gotten):
            if(len(gotten)<>N):
                return False
            else:
                return True
        gotten0=deepcopy(gotten)
    return False




