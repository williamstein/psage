from __future__ import print_function
#################################################################################
#
# (c) Copyright 2010 Fredrik Stroemberg
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
################################################################################

from psage.modform.maass import *


def test_vv(N=11):
    return 
    WR=WeilRepDiscriminantForm(11,dual=True)    
    M=VVHarmonicWeakMaassForms(WR,0.5,15)
    PP={(7/22,0):1}
    F=M.get_element(PP,12)
    print(F)


def test_construction_of_space(N=1):
    M = MaassWaveForms(N)
    M.get_element_in_range(9,10)

def test_group_list():
    l=get_list_of_valid_signatures(6)
    list_all_admissable_pairs(l[0],verbose=2)


def _test_permutation_iterator_1():
    r"""
    Make sure that the permutation iterator works as it should.
    Type of iterator: entire S_n
    """
    PI = MyPermutationIterator(4)
    assert len(PI.list())==24
    for x in PI:
        pass
    return

def _test_permutation_iterator_2():
    r"""
    Make sure that the permutation iterator works as it should.
    Type of iterator: set of permutations with a given set of fixed points.
    """
    PI=MyPermutationIterator(6,fixed_pts=[1,2])
    assert len(PI.list())==9
    for x in PI:
        pass
    return

def _test_permutation_iterator_3():
    r"""
    Make sure that the permutation iterator works as it should.
    Type of iterator: permutations without fixed points.
    """
    PI=MyPermutationIterator(6,num_fixed=0)
    assert len(PI.list())==265
    for x in PI:
        pass
    return


def _list_of_perms_test(N,type=1):
    r"""
    Make sure that the permutation iterator works as it should.
    Type of iterator: permutations without fixed points.
    """
    PI=MyPermutationIterator(N,num_fixed=0)
    if type==1:
        for x in PI:
            y=x
            pass
    elif type==2:
        try:
            while True: #for x in range(0,PI.max_num()):
                x = PI.current_perm()
                PI._next()
        except StopIteration:
            pass
    else:
        try:
            while True: #for x in range(0,PI.max_num()):
                x = PI.current_perm()
                PI.__next__()
        except StopIteration:
            pass        
    return

