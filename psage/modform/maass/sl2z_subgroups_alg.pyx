# cython: profile=True
# -*- coding: utf-8 -*-
#*****************************************************************************
#  Copyright (C) 2010 Fredrik Stromberg <stroemberg@mathematik.tu-darmstadt.de>,
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
r"""

Algorithms for enumerating and classifying subgroups of the modular group in terms of pairs of permutations.


AUTHOR:

 - Fredrik Stroemberg



"""

include "sage/ext/interrupt.pxi" 
include "sage/ext/stdsage.pxi"  
include "sage/ext/cdefs.pxi"

from permutation_alg cimport MyPermutation,MyPermutationIterator,CycleCombinationIterator
from permutation_alg cimport print_vec,_conjugate_perm,_are_eq_vec,transposition,_mult_perm_unsafe,are_transitive_perm_c,perm_to_cycle_c,are_conjugate_perm,get_conjugating_perm_list,get_conjugating_perm_ptr_unsafe,num_cycles_c

from psage.modform.maass.mysubgroup import MySubgroup
from psage.modform.maass.mysubgroups_alg cimport SL2Z_elt
from sage.modules.vector_integer_dense cimport Vector_integer_dense

from permutation_alg import verbosity,sort_perms
from sage.all import deepcopy,copy,ZZ,vector,subsets
from sage.groups.perm_gps.permgroup import PermutationGroup
from sage.groups.perm_gps.permgroup_named import SymmetricGroup

from sage.combinat.combinat import tuples
from sage.all import SL2Z
from sage.interfaces.all import gap
from sage.rings.integer cimport Integer
from time import clock, time

from logging import warning,error

## Namingscheme: oldX is more recent with higher number of X

cpdef list_all_admissable_pairs(sig,int get_details=1,int verbose=0,int get_one_rep=0,int congruence=-1,int do_strict=0,int check=0):
    r"""
    List all possible pairs (up to conjugacy) of admissible permutations E,R
    corresponding to groups G with signature = sig
    get_details = True/False --> compute number of conjugates of the group G
    INPUT:
    - sig --  signature
    - get_details -- logical
    - verbose -- binary, additive

        - 8 -- set verbosity of MyPermutationIterator

    - get_one_rep -- set to one if you just want one representative for this signature.
    - congruence -- Integer. 1 or 0 to find a congruence or a non-congruence subgroup.
    
    """
    cdef int mu,h,e2,e3,g
    try:
        [mu,h,e2,e3,g]=sig
    except:
        raise ValueError, "Indata not of correct format! sig=%s" %(sig)
    # For simplicity we fix the permutation of order two, E first:
    ## We should check that the signature is valid (corresponds to a group)
    if 12*g<>12+mu-6*h-3*e2-4*e3:
        print "lhs=",12*g
        print "rhs=",12+mu-6*h-3*e2-4*e3
        print "Not a valid signature!"
        return
    cdef MyPermutationIterator PRI,PONEI,PONEI2
    cdef int *Sptr, *Tptr, *cycle, *cycle_lens
    cdef int *rr2, *rs
#    cdef MyPermutation rss,pres,rss_tmp,ee,r,rr,epp,rpp,ep,rp,eppp,rppp
    cdef MyPermutation S,R,Sp,Tp,Spc,Stest,Rtest,Rp,Rpc
    cdef MyPermutation pR,p,Spp,S_canonical

    Sptr = <int*>sage_malloc(sizeof(int)*mu)
    if not Sptr: raise MemoryError
    Tptr = <int*>sage_malloc(sizeof(int)*mu)
    if not Tptr: raise MemoryError
    cycle = <int*> sage_malloc(sizeof(int)*mu)
    if not cycle: raise MemoryError
    cycle_lens = <int*> sage_malloc(sizeof(int)*mu)
    if not cycle_lens: raise MemoryError
    rr2=<int *>sage_malloc(sizeof(int)*mu)
    if not rr2: raise MemoryError
    rs=<int *>sage_malloc(sizeof(int)*mu)
    if not rs: raise MemoryError
    cdef int j
    cdef long i,ii
    cdef int is_in
    cdef int end_fc # end of cycles in S containing fixed points of R
    end_fc = e2+2*e3+1
    #sig_on()
    for j from 0<=j <=e2-1:
        Sptr[j]=j+1
    j = e2
    while j<mu-1:
        Sptr[j]=j+2
        Sptr[j+1]=j+1        
        j=j+2
    if verbose>0:
        print print_vec(mu,Sptr)
    cdef list Sl,Tl
    Sl=[];  Tl=[]
    for j from 0 <= j <mu:
        Sl.append(Sptr[j])
    spS=str(Sl)
    Sp = MyPermutation(Sl)
    ## The canonical PSL(2,Z)-representative of S
    S_canonical = MyPermutation(Sl)
    ### Note that there are only two choices of S which are inequivalent modulo permutations fixing 1
    ### S = (1)...(e2)(1+e2 2+e2)...  or
    ### S = (1 2)(3)...(1+e2)(2+e2 3+e2)...  
    ### In the first stage we are only concerned about getting PSL(2,Z) conjugacy classes so we fix S of the first form.
    cdef int a,b,c
    if verbose>0:
        print "signature=",sig
    ## Then we have to find matching order 3 permutations R wih e3 fixed points
    ## Without loss of generality, up to PSL(2,Z) conjugacy we can choose those fixed points as
    ## e2+1, e2+3,...,e2+2*e3-1 
    cdef list rfx_list = []
    if mu>=3:
        rfx_list=range(e2+1,e2+2*e3,2)
    elif mu==2:
        rfx_list = [1,2]
    else:
        rfx_list = [1]
        
    cdef int* rfx
    rfx = <int*>sage_malloc(sizeof(int)*e3)
    for i in range(e3):
        if i < len(rfx_list):
            rfx[i]=rfx_list[i]
        else:
            rfx[i]=0
    if verbose>0:
        print "fixed pts for R={0}. No. fixed e3={1}".format(rfx_list,e3)
    if len(rfx_list)<>e3:
        raise ValueError, "Did not get correct number of fixed points!"
    cdef list list_of_R,list_of_Rs
    list_of_R=[]
    list_of_Rs=[]
    if verbosity(verbose,3):
        mpi_verbose= verbose % 8
    else:
        mpi_verbose=0
    if verbose>0:
        print "verbose=",verbose
        print "mpi_verbose=",mpi_verbose
    cdef int max_fixed_by_R = 0
    if rfx_list<>[]:
        max_fixed_by_R = max(rfx_list)+2
    else:
        max_fixed_by_R = e2 + 1
    cdef int first_non_fixed_elt = 1
    if 1 in rfx_list:
        first_non_fixed_elt = 2 
    PRI = MyPermutationIterator(mu,order=3,fixed_pts=rfx_list,num_fixed=len(rfx_list),verbose=mpi_verbose)
    if verbose>0:
        print "PRI.list=",printf("%p ", PRI._list_of_perms)
    #for R in P.filter(lambda x: x.fixed_points()==rfx):
    max_num = PRI.max_num()
    if verbose>0:
        print "max. num of possible R=",max_num        
        print "original fixedpts=",PRI.fixed_pts()
    cdef int num_rcycles
    cdef int* rcycle_lens=NULL
    cdef int* Rcycles=NULL
    cdef int rcycle_beg=0
    cdef int* gotten=NULL
    cdef int* used=NULL
    cdef int x
    cdef list rc
    cdef list Rcycles_list
    cdef int t1,t2,t3
    cdef list list_of_groups=[]
    cdef list list_of_groups_tmp=[]
    cdef list list_R_tmp=[]
    cdef dict dict_of_R_modpsl
    cdef int do_cnt,t,k,l,llg
    list_of_groups=[] 
    cdef int len_list_of_R 

    rcycle_lens = <int*>sage_malloc(sizeof(int)*mu)
    if not rcycle_lens: raise MemoryError
    Rcycles = <int*>sage_malloc(sizeof(int*)*mu*mu)
    if not Rcycles: raise MemoryError
    gotten = <int *>sage_malloc(sizeof(int)*mu)
    if not gotten: raise MemoryError
    #DEB
    sig_on()
    #cdef MyPermutation ptest
    #ptest = MyPermutation([[1, 3, 4], [2,  5,7], [8, 6,9], [10, 11, 12]])
    if used==NULL:
        used = <int *>sage_malloc(sizeof(int)*mu)
    if verbose>=0:
        start = time()
    cdef list Tcyc=[]
    cdef MyPermutation pT
    cdef int h_tmp=0
    cdef long checked=0
  
    if verbose>1:
        print "First cycle check test if R({0})>{1}".format(first_non_fixed_elt,max_fixed_by_R)
    #for pR in PRI:
    ## Avoid infinite loops
    cdef Integer counter
    counter = ZZ(0)
    #    mpz_init(counter)
    pR = MyPermutation(length=mu)
    #print "max =", PRI._max_num
    t = 0
    while mpz_cmp(counter.value,PRI._max_num.value)<=0:        
        if t==1:
            break
        pR = MyPermutation(length=mu)        
        mpz_add_ui (counter.value,counter.value,1)
        #print "counter=",counter
        checked+=1
        #for i in range(mu):
        pR.set_entries(PRI._current_perm)
#       pR._entries = PRI._current_perm
        t = PRI._next()
        #print "pR=",pR
        ## If a is the first element not fixed by R and x = max of fixed elments by R
        ## then we can always assume R(a)<=x+2
        if verbose>1:
            print "Checking first cycle! "            
        if pR._entries[first_non_fixed_elt-1]>max_fixed_by_R:
            #raise ArithmeticError," Should not have gotten this far! p = {0}".format(pR)
            continue
        if verbose>1:
            print "S=",S_canonical.cycles() #print_vec(mu,Sptr)
            print "R=",pR.cycles() #print_vec(mu,<int *>PRI._current_perm)
        if verbose>1:
            print "Checking the number of cusps!"
        _mult_perm_unsafe(mu,Sptr,<int *>pR._entries,Tptr)
        #T=E*R
        #Tcyc=perm_to_cycle_c(mu,Tptr)
        h_tmp = num_cycles_c(mu,Tptr)
        if verbose>1:
            pT = MyPermutation(length=mu)
            pT.set_entries(Tptr)
            print "Tp=",Tcyc
            print "current fixedpts=",PRI.fixed_pts()
            print "number of cusps=",h_tmp
        if h_tmp<>h:
            if verbose>1:
                print "Incorrect number of cusps. remove!"
                print "next!"
            continue
        if verbose>1:
            print "Checking transitivity!"
        # If we have only one cusp we are always transitive
        if h_tmp<>1 and not are_transitive_perm_c(<int*>S_canonical._entries,<int*>pR._entries,gotten,mu,mpi_verbose):
#        if not are_transitive_perm_c(<int*>S_canonical._entries,<int*>PRI._current_perm,gotten,mu,mpi_verbose):            
            continue


        if verbose>1:
            print "current fixedpts1=",PRI.fixed_pts()
            print "rfx=",print_vec(e3,rfx)
        if get_one_rep:
            if congruence<>-1:
                G=MySubgroup(o2=S_canonical,o3=pR)
                if G.is_congruence() and congruence==1:
                    return S_canonical,pR
                elif not G.is_congruence() and congruence==0:
                    return S_canonical,pR
            else:
                return S_canonical,pR
            continue
        ## We will now further shorten the list by trying to pick representatives for R modulo permutations fixing 1
        # We will do this by choosing the lexicographically 'smallest' possible 3-cycles
        # Recall that we have chosen:
        # S = (1)...(e2)(e2+1 e2+2)(e2+3 e2+4)...(mu-1 mu)
        # R = (e2+1)(e2+3)...(e2+2*e3-1)(a b c)...
        for x in range(mu):
            used[x]=0
        for i in range(e3):
            used[rfx[i]-1]=1            
        rcycle_beg=0
        Rcycles_list=pR.cycles_as_lists() #perm_to_cycle_c(mu,<int*>pR._entries)
        do_cont = 0
#        cdef int do_strict = 1
        for rc in Rcycles_list:
            if len(rc)<>3:
                continue
            # only the 3-cycles are relevant here (not the 1-cycles)
            # cy = (a b c)
            a=<int>rc[0]; b=<int>rc[1]; c=<int>rc[2]
            if verbose>1:
                print "a,b,c=",a,b,c
            used[a-1]=1 # We have fixed a
            if verbose>1:
                print "used=",print_vec(mu,used)
            ## If b and c are equivalent with respect to conjugation by a perm. p which preserves S, i.e. if
            ##   i) S(b)=b and S(c)=c, or
            ##  ii) (b b') and (c c') are two cycles of S with c and c' not used previously and not fixed points of R
            ## then we choose b < c
            #if (b<=e2 and c<=e2) or ((b>=end_fc and used[Sptr[b-1]-1]==0) and (c>=end_fc and used[Sptr[c-1]-1]==0)):
            if equivalent_integers_mod_fixS(b,c,mu,e2,end_fc,Sptr,used)==1:
                if verbose>1:
                    print b," (b c) and ",c," are equivalent in",rc
                if b>c:
                    if verbose>1:
                        print "remove (b>c)!"
                    do_cont = 1
                    break
            ### Apply this also to a 'rotated' version of the cycle, i.e. (b c a )

            if do_strict==1:
                if equivalent_integers_mod_fixS(a,c,mu,e2,end_fc,Sptr,used)==1:
                    if verbose>1:
                        print a," (a c) and ",c," are equivalent in",rc
                    if b>c:
                        if verbose>1:
                            print "remove (b>c)!"
                        do_cont = 1
                        break
#            for j from a+1<= j < b:
            for j in range(a+1,b):
               # I want to see if there is a j, equivalent to b, smaller than b
               if equivalent_integers_mod_fixS(b,j,mu,e2,end_fc,Sptr,used)==1:
                    #if t1 or (t2 and t3):
                    if verbose>1:
                        print b," (b) and ",j," are equivalent in",rc
                    if j<b:
                        if verbose>1:
                            print "remove!"
                        do_cont = 1
                        break #raise StopIteration()
            if do_cont==1:
                if verbose>1:
                    print "breaking0!"
                break
            used[b-1]=1
            if verbose>1:
                print "used=",print_vec(mu,used)
            for j in range(a+1,c):
                # I want to see if there is a j, equivalent to c, smaller than c and which is not used
                if (c<=e2 and j<=e2 and used[j-1]==0) or ((c>e2+e3+1 and used[Sptr[c-1]-1]==0 and used[c-1]==0) and (j>e2+e3+1 and used[Sptr[j-1]-1]==0 and used[j-1]==0)):
                    if verbose>1:
                        print j," (c) and ",c," are equivalent in",rc
                    if j<c:
                        if verbose>1:
                            print "remove!"
                        do_cont = 1 #raise StopIteration()
                        break
            if do_cont==1:
                if verbose>1:
                    print "breaking1!"
                break
            used[c-1]=1
            if verbose>1:
                print "used=",print_vec(mu,used)
        if do_cont==1:
            if verbose>1:
                print "next and continue!"
            continue
        else:
            pass
        ## If we are here, R is a true candidate.
        list_of_R.append(copy(pR))
        if verbose>1:
            print "added pR=",pR
            print "Checked {0} Rs out of max. {1}".format(checked,max_num)
        if verbose>1:
            print "current list of Rs and Ts:"
            for r in list_of_R:
                print ":::::::::::::::: ",r.cycles(),";",(S_canonical*r).cycles()

    if gotten<>NULL:
        sage_free(gotten)
#    mpz_clear(counter)
    if get_one_rep:
        if congruence==1:
            if verbose>=0:
                print "We didn't find any congruence subgroups..."
            return []
        elif congruence==0:
            if verbose>=0:
                print "We didn't find any non-congruence subgroups..."
            return []
        else:
            if verbose>=0:
                print "We didn't find any subgroups..."
            return []
    if used<>NULL:
        sage_free(used)
        used=NULL
    PRI._c_dealloc()
    if Rcycles<>NULL:
        sage_free(Rcycles)
        Rcycles=NULL
    if rcycle_lens<>NULL:
        sage_free(rcycle_lens)
    # Now list_of_R contains at least one representative for each group with the correct signature
    #DEB
    sig_off()
    list_of_R.sort()
    len_list_of_R = len(list_of_R)
    if verbose>=0:
        print "Time:",time()-start
        print "Original list of R="
        for i in range(len_list_of_R):
            print list_of_R[i]
    if verbose>=0:
        print "Number of original R's:",len_list_of_R
    # Then add all conjugates mod PSL(2,Z) to get all groups before filtering.
    # For uniformity we conjugate so that if 1 is not fixed then S(1)=2
    ## Do a temporary irst filter
    list_of_R_tmp = copy(list_of_R)
    if verbose>=0:
        start = time()
    dict_of_R_modpsl={}
    #cdef MyPermutation ptmp1,ptmp2
    #ptmp1 = MyPermutation([[1, 2, 8],[3, 5, 10],[4],[6],[7, 9, 11]])
    #ptmp2 = MyPermutation([[1, 7, 2], [3, 4, 9], [5, 11, 10], [6, 8, 12]])
    cdef int* plist=NULL
    plist = <int*>sage_malloc(sizeof(int)*mu)  
    for k in range(len_list_of_R):
        Rp = copy(list_of_R[k])
        Spc = copy(S_canonical)
        if Rp not in list_of_R_tmp:
            #print "cont"            
            continue
        if verbose>1:
            print "Compare=",k,Rp
        for l in range(k+1,len(list_of_R)):                        
            Rpc= list_of_R[l]
            if Rpc not in list_of_R_tmp:
                #print "cont"
                continue
            if verbose>1:
                print "With ",l,Rpc
            #if k==0 and l == 2:
            #    t=are_conjugate_pairs_of_perms_c(S_canonical,Rp,Spc,Rpc,plist,0,0,verbose=2) #,map_from=1,map_to=1#)
            #else:
            t=are_conjugate_pairs_of_perms_c(S_canonical,Rp,Spc,Rpc,plist,0,0) #,map_from=1,map_to=1)
            if t<>0:
                if verbose>0:
                    p = MyPermutation(length=mu)
                    p.set_entries(plist)
                    print "{0} is equivalent to {1} with p={2}".format(Rp,Rpc,p)
                list_of_R_tmp.remove(Rpc)
    if verbose>=0:
        print "Tmp list of R="
        for i from 0 <= i < len(list_of_R_tmp):
            print list_of_R_tmp[i]
    if verbose>=0:
        print "Number of reduced R's:",len(list_of_R_tmp)
        print "Time for zeroth filter= ",time()-start 
    list_of_R = copy(list_of_R_tmp)
    if verbose>=0:
        start = time()
    cdef int cntt = 0
    ### Recall that the original list of R's represent PSL-inequivalent groups
    ### We now want to add back PSL-conjugates
    cdef list list_of_j,pl
    cdef int do_cnt1
    len_list_of_R = len(list_of_R)
    for k in range(len_list_of_R):
        Rp = copy(list_of_R[k])
        ## We want to check if conjugating (S,R) with (1 j) gives a new group not in the list
        if verbose>0:
            print "----------------------------------"
            print "R_{0}={1}".format(k,Rp)
        do_cnt1 = 0        
        for l in range(len(list_of_groups)):
            Rpc = copy(<MyPermutation>list_of_groups[l][1])
            Spc = copy(<MyPermutation>list_of_groups[l][0])
            #t,p=are_conjugate_pairs_of_perms(S_canonical,Rp,Spc,Rpc,ret='perm',map_from=1,map_to=1)
            t=are_conjugate_pairs_of_perms_c(S_canonical,Rp,Spc,Rpc,plist,1,1)
            p = MyPermutation(length=mu)
            p.set_entries(plist)
            
            cntt+=1
            if t==1:
                if verbose>0:
                    print "{0} and {1} are equivalent with p={2}!".format(k,l,p)
                    print "Spc,Rpc=",Spc,Rpc
                do_cnt1 = 1
                break
        if do_cnt1==0:
            if verbose>0:
                print "Append Sp,Rp ({0})={1},{2}".format(len(list_of_groups),S_canonical,Rp)
            list_of_groups.append((S_canonical,Rp))
        else:
            pass #continue
        ## first get a list of j's for which the current (S,R)^( 1 j ) are inequivalent
        #if verbose>0:
        #    print "HERE!"Q
        dict_of_R_modpsl[Rp]={}
        list_of_j = []
        for  j in range(2,mu+1):
            do_cnt = 0
            #if dict_of_R_modpsl[Rp].has_key(j):
            #    if verbose>0:
            #        print "Append Sp,Rp ({0}) from dict={1},{2}".format(len(list_of_groups),S_canonical,dict_of_R_modpsl[Rp][j][0])
            #    list_of_groups.append(dict_of_R_modpsl[Rp][j][0])
            #    #continue
            for l in range(1,j):
                t,p=are_conjugate_pairs_of_perms(S_canonical,Rp,S_canonical,Rp,ret='perm',map_from=l,map_to=j)
                if t==1:
                    do_cnt = 1
                    if verbose>0:
                        print "(S,R)^(1 {0}) =_1 (S,R)^(1 {1})".format(j,l)
                    break
            if do_cnt == 0:
                list_of_j.append(j)
        if verbose>0:
            print "list_of_j=",list_of_j
        for j in range(2,mu+1):
            Spp = copy(S_canonical)
            Rpp = copy(Rp)
            Rpp.conjugate_with_transposition(1,j)
            Spp.conjugate_with_transposition(1,j)
            do_cnt = 0
            if j not in list_of_j:
                if verbose>0:
                    
                    print "Not adding: S,R^(1 {0})={1},{2}".format(j,Spp,Rpp)
                    #t,p=are_conjugate_pairs_of_perms(Sp,Rp,Spc,Rpp,ret='perm',map_from=1,map_to=1)
                    #print "Are indeed same group!"
                #continue
            for l in range(len(list_of_groups)):
                Rpc = copy(list_of_groups[l][1])
                Spc = copy(list_of_groups[l][0])
#                if Rpc==Rpp and Spc==Spp:
#                    continue
                t,p=are_conjugate_pairs_of_perms(Spp,Rpp,Spc,Rpc,ret='perm',map_from=1,map_to=1)
                cntt+=1
                if t==1:
                    if verbose>0:
                        print "({0})^{1} = {2},{3} \n is equivalent with p={4}!".format(Rp,j,Spp,Rpp,p)
                        print "to: Spc,Rpc=",Spc,Rpc                        
                    do_cnt = 1
                    break
            if do_cnt == 1:
                continue
            if verbose>0:
                print "R,S^(1 {0})       ={1},{2}".format(j,Rpp,Spc)
            ## We try to normalize the fixed points of S:
            pl = range(1,Rp.N()+1)
            pl[0]=j; pl[j-1]=1
            p = MyPermutation(pl)
            if Spp(1)<>1 and Spp(1)<>2:
                l = Spp(1)
                p.conjugate_with_transposition(2,l)
                Rpp.conjugate_with_transposition(2,l)
                Spp.conjugate_with_transposition(2,l)
                if verbose>1:
                    print "j=",j
                    print "Rpp=",Rpp
                    print "Spp=",Spp
                    print "p^{0}={1}".format(l,p)
                    print "Rpp^{0}={1}".format(l,Rpp)
                    print "Spp^{0}={1}".format(l,Spp)
                #if verbose>0:
                #Rpp.conjugate_with_transposition(fix1[l],fix2[l])
                #Spp.conjugate_with_transposition(fix1[l],fix2[l])
            if (j not in list_of_j and not dict_of_R_modpsl[Rp].has_key(j)) or do_cnt1==1:
                raise ArithmeticError," Contradiction !! "
            if verbose>0:                
                print "Appending (S,R){0}={1},{2},p={3}".format(len(list_of_groups),Spp,Rpp,p)
            #if (Spc,Rp) not in list_of_groups:
            dict_of_R_modpsl[Rp][j]=(Spp,Rpp),p
            list_of_groups.append((Spp,Rpp))
    if verbose>=0:
        print "Time for zeroth filter= ",time()-start
        print "num of comparisons:",cntt
    list_of_groups.sort(cmp=sort_perms2)
    if verbose>=0:
        print "List of different groups:"
        if verbose>0:
            print  "Dict:"
            for R in dict_of_R_modpsl.keys():
                print "D[{0}]=".format(R)
                for j in dict_of_R_modpsl[R].keys():
                    (Spp,Rpp),p = dict_of_R_modpsl[R][j]
                    print "\t {0} \t ({1}, {2}); {3}".format(j,Spp,Rpp,p)
        for (Rp,Sp) in list_of_groups:
            print (Rp,Sp)
        print "Number of groups:",len(list_of_groups)
    cdef int nlr = len(list_of_R)
    if plist<>NULL:
        sage_free(plist)

#### CHECK
    if check==1:
        list_of_groups_tmp = copy(list_of_groups)
        if verbose>=0:
            start = time()
        for k in range(len(list_of_groups)):
            Sp = copy(list_of_groups[k][0])
            Rp = copy(list_of_groups[k][1])
            for l in range(k+1,len(list_of_groups)):                        
                Spc = copy(list_of_groups[l][0])
                Rpc = copy(list_of_groups[l][1])
                if (Spc,Rpc) not in list_of_groups_tmp:
                    continue
                t,p=are_conjugate_pairs_of_perms(Sp,Rp,Spc,Rpc,ret='perm') #,map_from=1,map_to=1)
                if t==1 and p(1)==1:                
                    raise ArithmeticError,"ERROR: {0},{1} is equivalent to {2},{3} with p={4}".format(Rp,Sp,Rpc,Spc,p)
        if verbose>=0:
            print "Number of original tmp R's:",len(list_of_groups_tmp)
            print "Time for check filter= ",time()-start 

## END CHECK
        
    indicator_list=range(1,len(list_of_R)+1)
    list_of_R_tmp=[]
    cdef dict lc_psl,lc_pgl,lc_psl_maps,lc_pgl_maps
    lc_psl=dict()      # list of conjugates
    lc_psl_maps=dict() # maps between conjugates
    if verbose>=0:
        start = time()
    ## At this point dict_of_groups contain information about PSL(2,Z) conjugates
    # conjugates = {}
    # conjugate_maps = {}
    # cdef dict tmpd1,tmpd2
    # Gmodpsl = []
    # for Rpp in dict_of_R_modpsl:
    #     Gmodpsl.append(S_canonical,Rpp)
    #     tmpd1 = {'psl':[],'pgl':[]}
    #     tmpd2 = {'psl':[],'pgl':[]}
    #     for j in dict_of_R_modpsl[Rpp].dict():
    #         Sp,Rp,p = dict_of_R_modpsl[Rpp][j][0]
    #         tmpd1['psl'].append((Sp,Rp))
    #         tmpd2['psl'].append((Sp,Rp))
    #     conjugates[(S_canonical,Rpp)]=tmpd   
    # ## Add PGL-conjugates as well.
    # Gmodpgl = copy(Gmodpsl)
    # for S,R in Gmodpsl:
    #     Rs = R.square().conjugate(S)        
    #     for Sp,Rp in Gmodpsl:
    #         t,A = are_conjugate_pairs_of_perms(S,Rs,Sp,Rp)
    #         if t==1:
    #             conjugates[(S_canonical,Rpp)]=tmpd   
    #sig_on()
    conjugates,conjugate_maps,Gmodpsl,Gmodpgl=find_conjugate_pairs(list_of_groups,mu,verbose=verbose-1) #verbose-1)
    if check==1:
        # check conjugates:
        for S,R in conjugates.keys():
            for i in range(len(conjugates[(S,R)]['psl'])):
                S1,R1 = conjugates[(S,R)]['psl'][i]
                p = conjugate_maps[(S,R)]['psl'][i]
                S2 = copy(S); R2 = copy(R)
                S2 = S2.conjugate(p)
                R2 = R2.conjugate(p)
                if S1<>S2 or R1<>R2:
                    s="Error check failed! conjugates and maps do not match for PSL! \n"
                    s+="S,R= {0},{1}\n".format(S,R)
                    s+="p={0}\n".format(p)
                    s+="S1,R1= {0},{1}\n".format(S1,R1)
                    s+="S,R^p= {0},{1}\n".format(S2,R2)
                    error(s)
            for i in range(len(conjugates[(S,R)]['pgl'])):
                S1,R1 = conjugates[(S,R)]['pgl'][i]
                p = conjugate_maps[(S,R)]['pgl'][i]
                S2 = copy(S); R2 = copy(R)
                R2 = R2.pow(2).conjugate(S)
                S2 = S2.conjugate(p)
                R2 = R2 = R2.conjugate(p)
                if S1<>S2 or R1<>R2:
                    s="Error check failed! conjugates and maps do not match for PGL!\n"
                    s+="S*,R*= {0},{1}\n".format(S,R.pow(2).conjugate(S))
                    s+="p={0}\n".format(p)
                    s+="S1,R1= {0},{1}\n".format(S1,R1)
                    s+="(S*,R*)^p= {0},{1}\n".format(S2,R2)
                    error(s)


#sig_off()
#    lc_psl,lc_psl_maps,list_of_R=filter_list_mod_psl_mod_S_new(list_of_R,mu,e2,e3,Sp,verbose)
    ## in the list of PGL-conjugates we want to only have representatives
    ## of classes modulo PSL
    if verbose>=0:    
        print "Time for second filter= ",time()-start
        print "List of conjugacy classes mod PGL(2,Z):" 
        for S,R in Gmodpgl:
            print S,R
            #if conjugates.has_key((S,R)):
            #    t = conjugates[(S,R)]['psl']
            #    if t<>[]:
            #        print t
    ## We finally want to add a list of reflected groups
    reflections={}
    for S,R in Gmodpsl:
        reflections[(S,R)]={}
        Rs = R.square().conjugate(S)
        conj_pgl = conjugates[(S,R)]['pgl']
        conj_psl = conjugates[(S,R)]['psl']
        do_cont = 0
        for S1,R1 in conj_psl:
            t,A,p = are_conjugate_pairs_of_perms(S,Rs,S1,R1)
            if t==1:
                reflections[(S,R)]={'group':(S1,R1),'map':A,'perm':p}
                do_cont = 1
                break
        if do_cont==1:
            continue
        for S1,R1 in conj_pgl:
            t,A,p = are_conjugate_pairs_of_perms(S,Rs,S1,R1)
            if t==1:
                reflections[(S,R)]={'group':(S1,R1),'map':A,'perm':p}
                do_cont = 1
                break
        if do_cont==1:
            continue
        for Stest,Rtest in Gmodpsl:
            t,A,p = are_conjugate_pairs_of_perms(S,Rs,Stest,Rtest)
            if t==1:
                reflections[(S,R)]={'group':(Stest,Rtest),'map':A,'perm':p}
                break
    lens_of_cc=[]    
    for Stest,Rtest in Gmodpsl:
        lens_of_cc.append(len(conjugates[(Stest,Rtest)]['psl']))
    d = dict()
    d['sig']=sig
    d['numg']=len(list_of_groups)
    d['num_cc_psl']=len(Gmodpsl)
    d['num_cc_pgl']=len(Gmodpgl)
    d['groups']=list_of_groups
    d['groups_mod_psl']=Gmodpsl
    d['groups_mod_pgl']=Gmodpgl
    d['conjugates']=conjugates
    d['conjugate_maps']=conjugate_maps
    d['reflections']=reflections

    
    # check
    ## for S,R in Gmodpsl:
    ##     S1,R1 = reflections[(S,R)]['group'] 
    ##     p = reflections[(S,R)]['perm'] 
    ##     if (S1,R1) <> (S,R):
    ##         d['numg_mod_pgl']+=1
    ##     if (S1,R1) not in conjugates[(S,R)]['pgl'] or (S1,R1) not in conjugates[(S,R)]['psl']:
    ##         print "S1,R1=",S1,R1
    ##         print "psl=",conjugates[(S,R)]['psl']
    ##         print "pgl=",conjugates[(S,R)]['pgl']
    ##         print "Possible error!" ## Should issue a warning with warning module....
    
    for key in d.keys():
        if isinstance(d[key],MyPermutation):
            d[key].set_rep(0)
        elif isinstance(d[key],dict):
            for key1 in d[key].keys():
                if isinstance(d[key][key1],MyPermutation):
                    d[key][key1].set_rep(0)
                elif isinstance(d[key][key1],list):
                    for v in d[key]:
                        if isinstance(v,MyPermutation):
                            v.set_rep(0)
        elif isinstance(d[key],list):
            for v in d[key]:
                if isinstance(v,MyPermutation):
                    v.set_rep(0)
    if rr2<>NULL:
        sage_free(rr2)
        rr2=NULL
    if rs<>NULL:
        sage_free(rs)
        rs=NULL
    if Sptr<>NULL:
        sage_free(Sptr)
    if Tptr<>NULL:
        sage_free(Tptr)
    if cycle<>NULL:
        sage_free(cycle); cycle=NULL
    if cycle_lens<>NULL:
        sage_free(cycle_lens); cycle_lens=NULL

    #sig_off()
    return d


cdef int equivalent_integers_mod_fixS(int a,int b,int mu,int e2,int end_fix,int* Sptr,int* used):
    r"""
    Check if a and b are equivalent under a permutation which fixes S under conjugation. I.e. Check if
    i) S(b)=b and S(c)=c, or
    ii) (b b') and (c c') are two cycles of S with c and c' not used previously and not fixed points of R
     
    """
    if used[b-1]==1:
        return 0
    if a<=e2 and b<=e2:
        return 1
    if (a>=end_fix and used[Sptr[a-1]-1]==0) and (b>=end_fix and used[Sptr[b-1]-1]==0):
        return 1
    return 0
             






cpdef tuple find_conjugate_pairs_old(list listGin, int mu, int verbose=0,int mpi_verbose=0,int do_pgl=0):
    r""" Removes duplicates in listR modulo permutations of the form (1j) which keeps S invariant.

    The purpose of this step is to reduce the list we later have to check completely.
    If do_pgl = 1 we reduce modulo PGL(2,Z)/PSL(2,Z) instead, i.e. we check R ~ ER^2E 
    """
    cdef list Gmodpsl=[],Gmodpgl=[]
    cdef dict conjugates={},conjugate_maps={}
    cdef MyPermutation R,S,R1,Sc,Rpsl,Rpgl,Scp,pp,pi
    cdef int numr = len(listGin)
    cdef int* checked=NULL
    cdef int t=0
    pp = MyPermutation(length=mu,rep=0)
    checked = <int*>sage_malloc(numr*sizeof(int))
    for i in range(numr):
        checked[i]=0
    if verbose>0:
        print "Groups in:",listGin
    for i in range(0,numr):
        S,R = listGin[i]
        conjugates[(S,R)]={'psl':[(S,R)],'pgl':[]}
        conjugate_maps[(S,R)]={'psl':[pp],'pgl':[]}
        
    for i in range(0,numr):
        S,R = listGin[i]
        if verbose>0:
            print "Test R[{0}]={1}".format(i,R)
        # checked[i] = 0 if R[i] is not (yet) found to be a conjugate
        # checked[i] = 1 if R[i] is conjugate mod PSL but not (yet) PGL\PSL
        # checked[i] = 2 if R[i] is conjugate mod PGL\PSL but not (yet) PSL
        # checked[i] = 3 if R[i] is conjugate mod PGL\PSL and PSL                           
        if checked[i]<2:
            Gmodpgl.append((S,R))
        if checked[i] in [0,2]:
            Gmodpsl.append((S,R))
        if verbose>0:
            for j in range(numr):
                print "checked[{0}]={1}".format(j,checked[j])
        for j in range(i,numr):
            if checked[j]==3:
                continue
            Sc,Rpsl = listGin[j]
            Rpgl = R.square().conjugate(Sc)
            p = are_conjugate_perm(R,Rpsl)
            if p==0:
                print " {0} and {1} are not conjugate!".format(R,Rpsl)
            if verbose>0:
                print "------------------------------------------------------------------"
                print "Comp PSL R[{0}]={1}".format(j,Rpsl)
                print "p=",p #
            Scp = S.conjugate(p)
            if verbose>0:
                print "R=",R
                print "S=",S
                print "S^p=",Scp
            pp=are_conjugate_wrt_stabiliser(Rpsl,Scp,Sc,p,&t)
            if verbose>0:
                print "t=",t
                print "j=",j
            if t==1:
                if verbose>0:                    
                    print "p=",p
                    print "pp=",pp
                pp = p*pp
                Scp = S.conjugate(pp)
                if verbose>0:                    
                    print "pp*p=",pp
                    print "S^pp=",Sc
                Rcp = R.conjugate(pp)
                if Rcp<>Rpsl or Scp<>Sc:
                    raise ArithmeticError,"Error with PSL-conjugating map!"
                if listGin[j] not in conjugates[(S,R)]['psl']:
                    conjugates[(S,R)]['psl'].append(listGin[j])
                    pp.set_rep(0)
                    conjugate_maps[(S,R)]['psl'].append(pp)
                checked[j]+=1
                if (S,R) not in conjugates[listGin[j]]['psl']:
                    conjugates[listGin[j]]['psl'].append((S,R))
                    pi = pp.inverse(); pi.set_rep(0)
                    conjugate_maps[listGin[j]]['psl'].append(pi)
            if checked[j]<2:
                # Check PGL(2,Z)
                p = are_conjugate_perm(R,Rpgl)
                if verbose>0:
                    print "------------------------------------------------------------------"
                    print "Comp PGL R[{0}]={1}".format(j,Rpgl)
                    print "p=",p
                Scp = S.conjugate(p)
                if verbose>1 or (i==2 and j==2):
                    print "R=",R
                    print "S=",S
                    print "S^p=",Scp
                    pp=are_conjugate_wrt_stabiliser(Rpgl,Scp,Sc,p,&t,0,0,1)
                else:
                    pp=are_conjugate_wrt_stabiliser(Rpgl,Scp,Sc,p,&t)
                if t==1:
                    pp = p*pp
                    if verbose>0:                    
                        print "pp=",pp
                    Scp = S.conjugate(pp)
                    Rcp= R.conjugate(pp)
                    if Rcp<>Rpgl or Scp<>Sc:
                        raise ArithmeticError,"Error with PGL-conjugating map!"
                    conjugates[(S,R)]['pgl'].append(listGin[j])
                    pp.set_rep(0)
                    conjugate_maps[(S,R)]['pgl'].append(pp)
                    checked[j]+=2
                    if i==0 and j==3:
                        print "R=",R
                        print "newlist=",listGin[j]
                        print "conjugates=",conjugates[listGin[j]]['pgl']
                    if (S,R) not in conjugates[listGin[j]]['pgl']:
                        conjugates[listGin[j]]['pgl'].append((S,R))
                        if i==0 and j==3:
                            print "Appending!"
                            print "conjugates=",conjugates
                            
                        pi = pp.inverse(); pi.set_rep(0)
                        conjugate_maps[listGin[j]]['pgl'].append(pi)
    #if R in conjugates.keys():
    #    for R1 in conjugates[R]['psl']:
    ## in the list of PGL-conjugates we want to only have representatives
    ## of classes modulo PSL. Observe that prevsiously we have basically only checked
    ## conjugation modulo PGL / PSL = <J>
    cdef tuple grp,grp1
    cdef list Gmodpgl_tmp,Gmodpsl_tmp
    ### We want to choose a representative with S in canonical form.
    Gmodpsl_tmp = copy(Gmodpsl)
    for S,R in Gmodpsl:
        if (S.num_fixed()>0 and S(1)<>1) or S(1)<>2:
            i = Gmodpsl_tmp.index((S,R))
            print "Checking PSL:",S,R
            for S1,R1 in conjugates[(S,R)]['psl']:
                if (S1.num_fixed()>0 and S1(1)==1) or S1(1)==2:
                    Gmodpsl_tmp[i]=(S1,R1)
                    print "S1,R1 is ok:",S1,R1
                else:
                    print "S1 not ok:",S1
                if S1(1)==1:  # this is the best option
                    break
            if Gmodpsl_tmp[i]==(S,R):
                print "Could not find appropriate PSL representative for {0}".format((S,R))
    Gmodpsl = Gmodpsl_tmp    
    Gmodpgl_tmp = copy(Gmodpgl)
    cdef list conj_temp
    for S,R in Gmodpgl:
        if (S.num_fixed()>0 and S(1)<>1) or S(1)<>2 or (S,R) not in Gmodpsl:
            i = Gmodpgl_tmp.index((S,R))
            print "Checking PGL:",S,R
            conj_tmp = conjugates[(S,R)]['pgl']
            conj_tmp.extend(conjugates[(S,R)]['psl'])
            for S1,R1 in conj_tmp:
                if S1(1) in [1,2] and (S1,R1) in Gmodpsl:
                    Gmodpgl_tmp[i]=(S1,R1)
                    print "S1,R1 is ok:",S1,R1
                else:
                    print "S1 not ok:",S1
                if S1(1)==1:
                    break
            if Gmodpgl_tmp[i]==(S,R):
                print "Could not find appropriate PGL representative for {0}".format((S,R))
    Gmodpgl = Gmodpgl_tmp
    
    # for grp in conjugate_maps.keys():        
    #     testlist = copy(conjugates[grp]['pgl'])
    #     print "testlist=",testlist
    #     for grp1 in testlist:
    #         if grp1 not in Gmodpsl:
    #             i = conjugates[grp]['pgl'].index(grp1)
    #             print "grp=",grp
    #             print "remove grp1=",grp1
    #             print "i=",i
    #             conjugates[grp]['pgl'].remove(grp1)
    #             conjugate_maps[grp]['pgl'].pop(i)
    #             print "new list=",conjugates[grp]
    #             print "new list=",conjugate_maps[grp]
    if checked<>NULL:
        sage_free(checked)
    return conjugates,conjugate_maps,Gmodpsl,Gmodpgl



cpdef tuple find_conjugate_pairs(list listGin, int mu, int verbose=0,int mpi_verbose=0,int do_pgl=0):
    r""" Removes duplicates in listR modulo permutations of the form (1j) which keeps S invariant.

    The purpose of this step is to reduce the list we later have to check completely.
    If do_pgl = 1 we reduce modulo PGL(2,Z)/PSL(2,Z) instead, i.e. we check R ~ ER^2E 
    """
    cdef list Gmodpsl=[],Gmodpgl=[]
    cdef dict conjugates={},conjugate_maps={}
    cdef MyPermutation R,S,R1,Sc,Rpsl,Rpgl,Scp,pp,pi
    cdef int numr = len(listGin)
    cdef int* checked=NULL
    cdef int t=0
    cdef list Gmodpsl_tmp
    cdef int i,j,k
    pp = MyPermutation(length=mu,rep=0)
    checked = <int*>sage_malloc(numr*sizeof(int))
    for i in range(numr):
        checked[i]=0
    if verbose>=0:
        print "Groups in:",len(listGin)
    #for i in range(0,numr):
    #    S,R = listGin[i]
    ## First sort the list into PSL(2,Z) conjugacy classes
#    cdef MyPermutation SS,RR
#    SS = MyPermutation('(1)(2 3)(4 5)(6 7)(8 9)(10 11)')
#    RR = MyPermutation('(1 3 6)(2)(4)(5 7 8)(9 10 11)')
    cdef MyPermutation ppp
    ppp = MyPermutation(length=mu,rep=3)
    for i in range(0,numr):
        S,R = copy(listGin[i])
        pp = MyPermutation(length=mu,rep=3)
        if verbose>0:
            print "Test R[{0}]={1}".format(i,R)
        # checked[i] = 0 if R[i] is not (yet) found to be a conjugate
        # checked[i] = 1 if R[i] is conjugate mod PSL.
        if checked[i]==1:
            continue
        Gmodpsl.append((S,R))
        conjugates[(S,R)]={'psl':[(S,R)],'pgl':[]}
        conjugate_maps[(S,R)]={'psl':[pp],'pgl':[]}
        checked[i]=1
        for j in range(i,numr):
            if checked[j]==1:
                continue
            Sc,Rpsl = copy(listGin[j])
            #print "Sc,Rpsl=",j,listGin[j]
            p = are_conjugate_perm(R,Rpsl)
            if p==0:
                print " {0} and {1} are not conjugate!".format(R,Rpsl)
            if verbose>0:
                print "------------------------------------------------------------------"
                print "Comp PSL R[{0}]={1}".format(j,Rpsl)
                print "p=",p #
            Scp = S.conjugate(p)
            if verbose>0:
                print "R=",R
                print "S=",Sc
                print "S^p=",Scp
            pp=are_conjugate_wrt_stabiliser(Rpsl,Scp,Sc,p,&t,0,0) #,verbose-1)
            if t==1:
                if verbose>0:                    
                    print "p=",p
                    print "pp=",pp
                #  pp = p o pp
                #pp = p*pp
                for k in range(mu):
                    ppp._entries[k]=pp._entries[k]
                _mult_perm_unsafe(mu,(<MyPermutation>p)._entries,ppp._entries,pp._entries)

                Scp = copy(S.conjugate(pp))
                Rcp = copy(R.conjugate(pp))
                if verbose>0:                    
                    print "pp*p=",pp
                    print "S^pp=",Scp
                    print "R^pp=",Rcp
                if Rcp<>Rpsl or Scp<>Sc:
                    raise ArithmeticError,"Error with PSL-conjugating map!"
                if listGin[j] not in conjugates[(S,R)]['psl']:
                    conjugates[(S,R)]['psl'].append(listGin[j])
                    pp.set_rep(3)
                    conjugate_maps[(S,R)]['psl'].append(pp)
                    checked[j]=1
#    print "S,R=",SS,RR
#    SS,RR = Gmodpsl[0]
#    print "Conjugates(S,R)=",conjugates.get((SS,RR),[])
#    print "Conjugate_maps(S,R)=",conjugate_maps.get((SS,RR),[])      
#    print "keys:"
#    i = 0
#    for S,R in conjugates.keys():
#        print S,R#
#        if i==0:
#            print "Values=",conjugates[(S,R)]
#        i+=1
    ## We now have a list with PSL(2,Z) representatives and their conjugacy classes.
    ## We now have to find out which conjugacy classes merge when considered modulo PGL instead
    ## First, however, we change to a representative of canonical form.
    Gmodpsl_tmp = copy(Gmodpsl)
    cdef dict tmp_dict
    cdef list tmp_list
    for S,R in Gmodpsl:
        if (S.num_fixed()>0 and S(1)<>1) or S(1)<>2: # choose another representative
            i = Gmodpsl_tmp.index((S,R))
            #print "Checking PSL:",S,R
            for S1,R1 in conjugates[(S,R)]['psl']:                
                if (S1.num_fixed()>0 and S1(1)==1) or S1(1)==2:
                    SS = copy(S1); RR=copy(R1)
                    Gmodpsl_tmp[i]=(SS,RR)
                    #print "S1,R1 is ok:",S1,R1
                #else:
                #    print "S1 not ok:",S1
                if S1(1)==1:  # this is the best option
                    break
            if Gmodpsl_tmp[i]==(S,R):
                s = "Could not find appropriate PSL representative for {0}".format((S,R))
                warning(s)
            else:
                j = conjugates[(S,R)]['psl'].index((S1,R1))
                p = conjugate_maps[(S,R)]['psl'][j]
                
                tmp_dict = conjugates.pop((S,R))
                conjugates[(S1,R1)] = tmp_dict
                tmp_list = conjugate_maps[(S,R)]['psl']
                #tmp_dict = conjugate_maps.pop((S,R))
                for k in range(len(tmp_list)):
                    pp = tmp_list[k]
                    ptmp = p.inverse()*pp
                    if S1.conjugate(ptmp) == tmp_dict['psl'][k][0] and R1.conjugate(ptmp) == tmp_dict['psl'][k][1]:
                        if not conjugate_maps.has_key((S1,R1)):
                            conjugate_maps[(S1,R1)] = {'psl': range(len(tmp_list)),'pgl': []}
                            
                        conjugate_maps[(S1,R1)]['psl'][k]=ptmp
                    else:
                        s = "Did not conjugate correctly! \n"
                        s+= "S,R={0},{1} \n".format(S,R)
                        s+="p=",p.inverse(),'*',pp,'=',ptmp
                        s+= "S1,R1={0},{1} \n".format(S1,R1)
                        s+= "S1,R1^p={0},{1} \n".format(S1.conjugate(ptmp),R1.conjugate(ptmp))
                        s+ "conjugates[S,R][{0}]={1}".format(k,tmp_dict['psl'][k])
                        warning(s)
                #conjugate_maps[(S1,R1)] = tmp_dict 
    Gmodpsl = Gmodpsl_tmp
#    print "Conjugates(S,R)=",conjugates.get((SS,RR),[])
#    print "Conjugate_maps(S,R)=",conjugate_maps.get((SS,RR),[])    
#    print "Groups mod PSL:",
#    for S,R in Gmodpsl:
#        print S,R
#    print "keys2:"
#    for S,R in conjugates.keys():
#        print S,R
    if checked<>NULL:
        sage_free(checked)
    numr = len(Gmodpsl)
    checked = <int*>sage_malloc(numr*sizeof(int))
    if checked==NULL: raise MemoryError
    for i in range(numr):
        checked[i]=0
    for i in range(numr):
        if checked[i]==1:
            continue
        S,R = Gmodpsl[i]
        Gmodpgl.append((S,R))
        Rpgl = R.square().conjugate(S)
        if verbose>0:
            print "------------------------------------------------------------------"
            print "Comp PGL: R={0}; S={1}".format(R,S)
            print "ER^2S[{0}]={1}".format(i,Rpgl)
        for j in range(i,numr):
            if checked[j]==1:
                continue
            # Check PGL(2,Z)
            S1,R1 = Gmodpsl[j]
            p = are_conjugate_perm(Rpgl,R1)
            Scp = S1.conjugate(p)
            if verbose>0:
                print "p=",p
                print "R=",R
                print "S=",S
                print "S^p=",Scp
                pp=are_conjugate_wrt_stabiliser(R1,Scp,S,p,&t,0,0)
            else:
                pp=are_conjugate_wrt_stabiliser(R1,Scp,S,p,&t)
            if t==1:
                pp = p*pp
                if verbose>0:                    
                    print "pp=",pp
                Scp = S.conjugate(pp)
                Rcp= Rpgl.conjugate(pp)
                if verbose>0:
                    print "Rpgl=",Rpgl
                    print "S=",S
                    
                    print "S^pp=",Scp
                    print "R^pp=",Rcp
                if Rcp<>R1 or Scp<>S1:
                    raise ArithmeticError,"Error with PGL-conjugating map!"
                conjugates[(S,R)]['pgl'].append((S1,R1))
                #tmp_dict = conjugates.pop((S,R))
                conjugate_maps[(S,R)]['pgl'].append(pp)
                checked[j]=1
    if verbose>=0:
        print "Groups out:",len(Gmodpsl)
    if checked<>NULL:
        sage_free(checked)
    return conjugates,conjugate_maps,Gmodpsl,Gmodpgl



cpdef are_conjugate_groups(G1,G2,ret='SL2Z',coset_rep=1,check=0,verbose=0):
    r"""
    Determine whether G1 and G2 are conjugate in PSL(2,Z) and return either a permutation or a matrix A in SLZ which performs the conjugation, i.e. A^-1 G1 A = G2.

    In the first step we find p s.t. p R1 p^(-1) = R2
    and then we find s s.t. s R2 s^-1 = R2 and s (p S1 p^(-1)) s^-1 = S2
    """
    cdef MyPermutation R1,R2,S1,S2,p,Sc,pp
    cdef int t
    if G1.signature()<>G2.signature():
        if ret=='SL2Z':
            return 0,SL2Z([1,0,0,1])
        else:
            return 0, MyPermutation(length=G1.index())
    R1 = G1.permR; R2=G2.permR
    S1 = G1.permS; S2=G2.permS    
    t,pp = are_conjugate_pairs_of_perms(S1,R1,S2,R2,ret='perm',verbose=verbose)
    A = G1.coset_reps()[pp(1)-1]
    return t,A

cpdef tuple are_conjugate_pairs_of_perms(MyPermutation S1,MyPermutation R1,MyPermutation S2,MyPermutation R2,str ret='all',int map_from=0,int map_to=0,int verbose=0):
    r"""
    Check if the pairs of permutations (S1,R1) and (S2,R2) are conjugate.
    If map_to>0 and map_frpm>0 we require any conjugating permutation p to satisfy p(map map_from)=map_to

    INPUT::

    - 'S1' -- Permutation
    - 'R1' -- Permutation
    - 'S2' -- Permutation
    - 'R2' -- Permutation    
    - 'ret'-- string (default 'all') Decides the format of the return value.
              possible values: 'perm','SL2Z','all'
    - 'map_from' -- integer (default 0) 
    - 'map_to' -- integer (default 0) 
    - 'verbose' -- integer (default 0) 
    """

    p = are_conjugate_perm(R1,R2)
    cdef MyPermutation pp0,pp1,pp,Sc,Scp
    if p==0:
        if ret == 'perm':
            return 0,MyPermutation(length=S1.N())
        elif ret=='SL2Z':
            return 0,SL2Z_elt(1,0,0,1)
        else:
            return 0, SL2Z_elt(1,0,0,1),MyPermutation(length=S1.N())
    pp = <MyPermutation?>p
    Sc = S1.conjugate(pp)
    if verbose>0:
        print "R1,R2=",R1,R2
        print "S1,S2=",S1,S2
        print "p=",pp
        print "S^p=",Sc
    cdef int j,t=0
    if map_to<>0:        
        j = pp(map_from)
        if verbose>0:
            print"Need sigma mapping {0} to {1}".format(map_from,j)
        pp=are_conjugate_wrt_stabiliser(R2,Sc,S2,pp,&t,j,map_to,verbose-1)
        if verbose>0:
            print "want total map to take {0} to {1}".format(map_from,map_to)
            print "p=",p
            print "pp=",pp
            if pp<>None:
                print "pp*p=",pp*p
                print "p*pp=",p*pp
    else:
        pp=are_conjugate_wrt_stabiliser(R2,Sc,S2,pp,&t,0,0,verbose-1)
    if t == 0:
        if ret=='perm':
            return 0, MyPermutation(length=S1.N())
        elif ret=='SL2Z':
            return 0, SL2Z_elt(1,0,0,1)        
        else:
            return 0, SL2Z_elt(1,0,0,1),MyPermutation(length=S1.N())
    if verbose>0:
        pp0 = copy(pp)
        pp1 = p*pp
    pp = p*pp
    if map_from>0:
        if pp._entries[map_from-1]<>map_to:
            raise ArithmeticError,"Could not find permutation conjugating and sending {0} to {1}".format(map_from,map_to,pp)
    if verbose>0:
        print "pp = pp*p=",pp
    Scp = S1.conjugate(pp)
    if verbose>0:
        print "S1=",S1
        print "S^pi=",Sc
        print "S^(pi*sigma)=",Scp
        Sc = Sc.conjugate(pp0)
        print "(S^pi)^sigma)=",Sc
        print "S^(sigma*pi)=",S1.conjugate(pp1)
    if Scp<>S2:
        raise ArithmeticError,"Conjugation did not work!"
    if ret=='perm':
        return 1,pp
    if pp.is_identity():
        if ret=='SL2Z':
            return 1,SL2Z([1,0,0,1])
        else:
            return 1,SL2Z([1,0,0,1]),pp
#    return 1,matrix_from_perm((S1,R1),pp,verbose)
    cdef SL2Z_elt V
    V = MySubgroup(o2=S1,o3=R1).coset_reps()[pp(1)-1]
    #A = V.SL2Z()
    if ret=='SL2Z':
        return 1,V
    else:
        return 1,V,pp


cpdef are_conjugate_pairs_of_perms_newer(MyPermutation S1,MyPermutation R1,MyPermutation S2,MyPermutation R2,str ret='all',int map_from=0,int map_to=0,int verbose=0):
    cdef int* plist
    plist = <int*>sage_malloc(sizeof(int)*S1.N())
    return  are_conjugate_pairs_of_perms_c(S1,R1,S2,R2,plist,map_from,map_to,verbose)

cdef int are_conjugate_pairs_of_perms_c(MyPermutation S1,MyPermutation R1,MyPermutation S2,MyPermutation R2,int* plist,int map_from=0,int map_to=0,int verbose=0):
    r"""

    If map_to>0 we require any permutation p to satisfy p(map map_from)=to map_to
    """

    #S1,R1 = pair1; S2,R2 = pair2
    ## First check if R1 and R2 are conjugate
    cdef int i,j,t=0
    cdef MyPermutation pp0,pp1,pp,Sc,Scp
    p = are_conjugate_perm(R1,R2)
    for i in range(S1.N()):        
        plist[i] = i+1
    if isinstance(p,int):
        return 0
    pp = <MyPermutation?>p
    Sc = S1.conjugate(pp)
    if verbose>0:
        print "R1,R2=",R1,R2
        print "S1,S2=",S1,S2
        print "p=",pp
        print "S^p=",Sc
    if map_to<>0:        
        j = pp(map_from)
        if verbose>0:
            print"Need sigma mapping {0} to {1}".format(map_from,j)
        pp=are_conjugate_wrt_stabiliser(R2,Sc,S2,pp,&t,j,map_to,verbose-1)
        if verbose>0:
            print "want total map to take {0} to {1}".format(map_from,map_to)
            print "p=",p
            print "pp=",pp
            if pp<>None:
                print "pp*p=",pp*p
                print "p*pp=",p*pp
    else:
        pp=are_conjugate_wrt_stabiliser(R2,Sc,S2,pp,&t,0,0,verbose-1)
    if t == 0:
        return 0
    if verbose>0:
        pp0 = copy(pp)
        pp1 = p*pp
    pp = p*pp
    if map_from>0:
        if pp._entries[map_from-1]<>map_to:
            raise ArithmeticError,"Could not find permutation conjugating and sending {0} to {1}".format(map_from,map_to,pp)
    if verbose>0:
        print "pp = pp*p=",pp
    Scp = S1.conjugate(pp)
    if verbose>0:
        print "S1=",S1
        print "S^pi=",Sc
        print "S^(pi*sigma)=",Scp
        Sc = Sc.conjugate(pp0)
        print "(S^pi)^sigma)=",Sc
        print "S^(sigma*pi)=",S1.conjugate(pp1)
    if Scp<>S2:
        raise ArithmeticError,"Conjugation did not work!"
    for i in range(S1._N):        
        plist[i] = pp._entries[i]
    return 1
#    return 1,matrix_from_perm((S1,R1),pp,verbose)
#   A = MySubgroup(o2=S1,o3=R1).coset_reps()[pp(1)-1].SL2Z()
#    return 1,A
    

cpdef matrix_from_perm(gens,p,verbose=0):
    r"""
    If possible express the permutation p as a word in gens[0]=S and gens[1]=R and translate into a matrix in SL(2,Z).
    If desired we also find a representative modulo the subgroup given by the pair (S,R).
    """
    S,R = gens
    cdef int mu = R.N()
    G = SymmetricGroup(mu)
    G1= MySubgroup(o2=S,o3=R)
    g = G(p.list())
    H = gap.Group([G(S.list()),G(R.list())])
    if verbose>0:
        print "gens=({0},{1})".format(S,R)
        print "G=",G
        print "mu=",mu
        print "p=",p
        print "g=",g
        print "H=",H
    if g in H:
        s = H.EpimorphismFromFreeGroup().PreImagesRepresentative(g)
    else:
        return None
    if verbose>0:
        print "s=",s
    s = str(s)
    A = str_to_slz(s)
    ## If possible we also try to simplify the conjugating matrix
    for i in range(G1.generalised_level()):
        Ti = SL2Z([1,-i,0,1])
        if A*Ti in G1:
            return SL2Z([1,i,0,1])
    return A 



cpdef str_to_slz(str):
    S,T = SL2Z.gens()
    R = T*S
    sl = str.split('*')
    res = SL2Z([1,0,0,1])
    for t in sl:
        tt = t.split('^')
        if len(tt)==1:
            k = 1
        elif len(tt)==2:
            k = tt[1]
        else:
            raise ArithmeticError
        if tt[0]=='x1':
            res = res*S**k
        else:
            res = res*R**k
    return res


cpdef are_mod1_equivalent(MyPermutation R1,MyPermutation S1, MyPermutation R2,MyPermutation S2,int verbose=0):
    cdef MyPermutation p
    cdef int* pres=NULL
    cdef int res
    cdef int N=R1._N
    assert S1._N==N and R2._N==N and S2._N==N
    pres = <int*>sage_malloc(sizeof(int)*N)
    p=MyPermutation(length=N)
    #DEB sig_on()
    res = are_mod1_equivalent_c(R1._N,S1,R1,S2,R2,pres,verbose)
    #DEB sig_off()
    p.set_entries(pres)
    if verbose>0:
        print "res=",res
        print "pres=",print_vec(N,pres)
        print "p=",p
    if pres<>NULL:
        sage_free(pres)
    return res,p


cdef int are_mod1_equivalent_c(int N,MyPermutation S1,MyPermutation R1,MyPermutation S2,MyPermutation R2,int* pres,int verbose=0,int check=0):
    r"""
    Check if S1 ~_1 S2 and R1 ~_1 R2 (i.e. equiv. mod 1)


    If check=0 we assume that the cycles are compatible (as in the case if we call this from routine above)
    """
    
    cdef int i,j,ii,cont
    cdef MyPermutationIterator PONEI
    cdef int *epp, *rpp
    cdef int fixS,fixR,res
    epp = NULL
    rpp=NULL
    for i in range(N):
        pres[i]=i+1
    #pres = MyPermutation(length=N)._entries
    res=1
    if _are_eq_vec(N,S1._entries,S2._entries) and _are_eq_vec(N,R1._entries,R2._entries):
        return res
        #return True,MyPermutation(range(1,N+1))
    # Then check if they are compatble at all:
    fixS=0; fixR=0
    if check:
        cl1=map(len,S1.cycles())
        cl2=map(len,S2.cycles())
        cl1.sort(); cl2.sort()
        if cl1<>cl2:
            res=0
        else:
            if verbose>0:
                #print "S1cyc=",S1.to_cycles()
                print "cl1=",cl1
            fixS=cl1.count(1)
            cl1=map(len,R1.cycles())
            cl2=map(len,R2.cycles())
            cl1.sort(); cl2.sort()
            if cl1<>cl2:
                res=0
            else:
                if verbose>0:
                    #print "R2cyc=",R2.to_cycles()
                    print "cl2=",cl2
                fixR=cl2.count(1)
    else:
        fixS=S1.num_fixed_c()
        fixR=R1.num_fixed_c()
    if res==0:
        return res
    res=0
    if verbose>0:
        print "testing G1:={0}:{1}".format(S1,R1)
        print "testing G2:={0}:{1}".format(S2,R2)
    epp =  <int *>sage_malloc(sizeof(int)*N)
    if not epp:  raise MemoryError
    rpp =  <int *> sage_malloc(sizeof(int)*N)
    if not rpp:  raise MemoryError
    res = 0
    if verbose>0:
        print "HERE 0!"
    #pres = MyPermutation(length=N)
    cdef list thisset,fixptsets,listRout
    cdef list fixptsS1=S1.fixed_elements()
    cdef list fixptsR1=R1.fixed_elements()
    cdef list fixptsS2=S2.fixed_elements()
    cdef list fixptsR2=R2.fixed_elements()
    if verbose>0:
        print "HERE 1!"
    cdef MyPermutation p,p0
#    p0=transposition(12,5,7)*transposition(12,6,8)*transposition(12,10,7)*transposition(12,9,8)*transposition(12,10,12)*transposition(12,9,11)*transposition(12,11,12)
    #print "pres._entries=",printf("%p ",pres._entries)
    # Check which sets are allowed:
    cdef list iteratorlist=[]
    #print "fixS1=",fixptsS1
    #print "fixS2=",fixptsS1
    #print "fixR1=",fixptsR1
    #print "fixR2=",fixptsR2
    for thisset in subsets(range(1,N+1)):
        if len(thisset)>=N-1:  # Then all pts are fixed and we have a trivial permutation only
            continue
        if 1 not in thisset:
            continue
        if verbose>1:
            print "fixset test.:",thisset
        cont=False
        for a in thisset:
            if ((a in fixptsS1) and (a not in fixptsS2)) or ((a in fixptsS2) and (a not in fixptsS1)):
                cont=True
                if verbose>2:
                    print "a"
                break
            if ((a in fixptsR1) and (a not in fixptsR2)) or ((a in fixptsR2) and (a not in fixptsR1)):
                cont=True
                if verbose>2:
                    print "b"
                break
        if cont:
            continue # fixptsets.remove(thisset)
        if verbose>1:
            print "fixset ok.:",thisset
        PONEI = MyPermutationIterator(N,fixed_pts=thisset,num_fixed=len(thisset))
        iteratorlist.append(PONEI)
    #sig_on()
    if verbose>0:
        print "HERE!"
    for PONEI in iteratorlist:      
        #if verbose>1:
        #    print "Checking perms in ",PONEI
        for p in PONEI:
            if verbose>2:
                print "p=",p
            #if p==p0:
            #    print "p0 is here!"
            _conjugate_perm(N,epp,S1._entries,p._entries) # epp=p*E*p.inverse()
            #if p==p0:
            #    print "epp=",print_vec(N,epp)
            if not _are_eq_vec(N,epp,S2._entries):
                continue
            _conjugate_perm(N,rpp,R1._entries,p._entries) # rpp=p*R*p.inverse()
            #if p==p0:
            #    print "rpp=",print_vec(N,rpp)
            if _are_eq_vec(N,rpp,R2._entries):
                res = 1
                for i in range(N):
                    pres[i]=p._entries[i]
                break
        if res==1:
            break
        #PONEI._c_dealloc()
    #print "end pres cur=",pres
    #print "pres._entries=",printf("%p ",pres._entries)
    if verbose>0:
        print "Are mod 1 equivalent:",res
        print "pres=",print_vec(N,pres)
    if epp<>NULL:
        sage_free(epp)
    if rpp<>NULL:
        sage_free(rpp)
    #sig_off()

    return res







cpdef is_consistent_signature(int ix,int nc=0,int e2=0,int e3=0,int g=0):
    #print "rhs=",12+ix-6*nc-3*e2-4*e3
    return int(12*g == 12+ix-6*nc-3*e2-4*e3)








cpdef stabiliser_of_R_StoSp(MyPermutation pR,MyPermutation pS,MyPermutation pS1,int verbose=0):
    r"""
    Return a the list of permutations which preserves pR and maps pS to pS1 under conjugation.
    """
    assert pS.order() in [1,2]
    assert pS1.order() in [1,2]
    assert pR.order() in [1,3]
    cdef list Rctype = pR.cycle_type()
    cdef list Rcycles = pR.cycles()
    cdef lR = pR.cycles_dered_as_list()
    cdef list cycles1,cycles3
    cdef int nf,nf1,d1,d3
    cdef list fixedS,fixedS1,nonfixedS,nonfixedS1,res,one_cycles,three_cycles
    d1 = Rctype.count(1)
    d3 = Rctype.count(3)
    cycles1 = []; cycles3=[]
    for c in Rcycles:
        if len(c)==1:
            cycles1.append(c)
        elif len(c)==3:
            cycles3.append(c)
    if verbose>0:
        print "cycles1=",cycles1
        print "cycles3=",cycles3
        print "lR=",lR
    nf = pS.num_fixed(); nf1 = pS1.num_fixed()
    if nf<>nf1:
        return []
    fixedS=pS.fixed_elements();  nonfixedS=pS.non_fixed_elements()
    fixedS1=pS1.fixed_elements();  nonfixedS1=pS1.non_fixed_elements()
    if verbose>0:
        print "fixedS =",fixedS
        print "non-fixedS =",nonfixedS
        print "fixedS1=",fixedS1
        print "non-fixedS1=",nonfixedS1
    ## These permuations permutes orders of the 1- and 3-cycles separately
    cdef MyPermutationIterator perms1,perms3
    cdef MyPermutation p,pp
    perms1 = MyPermutationIterator(max(d1,1))
    perms3 = MyPermutationIterator(max(d3,1))
    ## Then we also have to include permutations by the 3-cycles of R
    three_cycles_perm =[]
    for p in pR.cycles(type='perm'):
        if p.is_order_eq(3)==1:            
            three_cycles_perm.append(p)
    cdef list three_cycles_R = []
    perm_tuples = tuples([0,1,2],len(cycles3))
    for i in range(len(perm_tuples)):
        p = MyPermutation(length=pS.N()) # identity
        for j in range(len(cycles3)):
            if perm_tuples[i][j]==1:
                p = p*three_cycles_perm[j]
            elif perm_tuples[i][j]==2:
                p = p*three_cycles_perm[j].square()
            #print "cy=",cy
        three_cycles_R.append(p)
    res = [] 
    for p1 in perms1:
        if verbose>0:
            print "perm1=",p1
            print "cy1=",cycles1
        if cycles1<>[]:
            one_cycles = p1(cycles1) #permute_list(cycles[1],p1)
        else:
            one_cycles = []
        l0 = copy(one_cycles)
        for p3 in perms3:
            l = copy(l0)
            three_cycles = p3(cycles3) #permute_list(cycles[3],p3)
            l.extend(three_cycles)
            perm = MyPermutation(l)
            ll = flatten_list2d(l)
            p = get_conjugating_perm_list(ll,lR)
            if verbose>0:
                print "ll=",ll
                print "perm3=",p3
                print "perm0=",p
            ## We now act with all combinations of cycles from R
            for ptmp in three_cycles_R:
                do_cont = 0
                pp = p*ptmp
                if verbose>0:
                    print "ptmp=",ptmp
                    print "pp=",pp
                # if S=pSp^-1 and Sx=x then S(px)=pS(x)=px
                for x in fixedS:
                    if pp(x) not in fixedS1:
                        do_cont = 1
                        break
                if do_cont==1:
                    if verbose>0:
                        print "Fails at preserving fixed points!"
                    continue
                for x in nonfixedS:
                    if pp(x) not in nonfixedS1:
                        do_cont = 1
                        break
                if do_cont==1:
                    if verbose>0:
                        print "Fails at preserving non-fixed points!"
                    continue        
                pScon = pS.conjugate(pp)
                if pScon<>pS1:
                    if verbose>0:
                        print "Fails at preserving S!"
                        print "ppSpp^-1=",pScon
                        print "pS1=",pS1
                    continue
                res.append(pp)
            ## We only consider permutations which fixes S
        ## We assume that p ionly has cycles of length 1 or 3 (i.e. it has order 3)
    return res

cpdef  are_conjugate_wrt_stabiliser2(MyPermutation pR,MyPermutation pS,MyPermutation pS1,MyPermutation p_in,int map_from=0,int map_to=0,int verbose=0):
     cdef MyPermutation res
     cdef int tt
     res = are_conjugate_wrt_stabiliser(pR,pS,pS1,p_in,&tt,map_from,map_to,verbose)
     return tt,res


cdef MyPermutation are_conjugate_wrt_stabiliser(MyPermutation pR,MyPermutation pS,MyPermutation pS1,MyPermutation p_in,int* t,int map_from=0,int map_to=0,int verbose=0):
    r"""
    Return a the list of permutations which preserves pR and maps pS to pS1 under conjugation.
    To be more precise. We seek permutation p s.t.
    p^-1 pR p = pR and p^-1 pS p = pS1   and p(map_from)=map_to

    """
    cdef list Rctype,Rcycles,lR #,l0,l
    cdef list cycles1,cycles3
    cdef int nf,nf1,d1,d3,i,j,mu,do_cont,dd1,dd3,ir,nir
    cdef list fixedS,fixedS1,nonfixedS,nonfixedS1,res,one_cycles,three_cycles
    cdef MyPermutationIterator perms1,perms3
    cdef MyPermutation p,pp,ptmp,pScon,p1,p3,perm
    t[0] = 0
    cdef int * ll = NULL
    mu = <int>pR.N()
    cdef int num_one_cycles = pR.num_fixed()
    cdef int num_three_cycles = (mu - num_one_cycles)/3
    ll = <int*> sage_malloc(sizeof(int)*mu)
    if ll == NULL:
        raise MemoryError
    if pS.is_order_eq_c(2)<>1 or pS1.is_order_eq_c(2)<>1 or pR.is_order_eq_c(3)<>1:
        pp = MyPermutation(length=mu)
        return pp
    nf = pS.num_fixed_c(); nf1 = pS1.num_fixed_c()
    if nf<>nf1:
        pp = MyPermutation(length=mu)
        return pp
    pp = MyPermutation(length=mu)
    lR = pR.cycles_ordered_as_list()
    Rcycles = pR.cycles() #o_cycles() #'list')
    Rctype = pR._cycle_type
    cdef int * llR=NULL
    llR = <int*> sage_malloc(sizeof(int)*mu)
    if llR == NULL:
        raise MemoryError
    for i in range(mu):
        llR[i]=<int>lR[i]
    d1 = num_one_cycles # Rctype.count(1)
    d3 = num_three_cycles #Rctype.count(3)
    cycles1 = []; cycles3=[]
    cdef int* cycles1ptr = NULL
    cycles1ptr = <int*>sage_malloc(sizeof(int)*num_one_cycles)
    if cycles1ptr == NULL:
        raise MemoryError
    cdef int** cycles3ptr = NULL
    cycles3ptr = <int**>sage_malloc(sizeof(int*)*num_three_cycles)
    if cycles3ptr == NULL:
        raise MemoryError
    #for i in range(num_one_cycles):
    #    cycles1ptr[i]=cycles1[i][0]
    for i in range(num_three_cycles):
        cycles3ptr[i] = NULL
        cycles3ptr[i] = <int*>sage_malloc(sizeof(int)*3)
        if cycles3ptr[i] == NULL:
            raise MemoryError
    cdef int i1=0,i3=0
    for c in Rcycles:
        if len(c)==1:
            #cycles1.append(c)
            cycles1ptr[i1]=c[0]
            i1+=1
        elif len(c)==3:
            #cycles3.append(c)    
            cycles3ptr[i3][0] = c[0]
            cycles3ptr[i3][1] = c[1]
            cycles3ptr[i3][2] = c[2]
            i3+=1
    if verbose>0:
        for i in range(num_one_cycles):
            cycles1.append(cycles1ptr[i])
        for i in range(num_three_cycles):
            cycles3.append([cycles3ptr[i][0],cycles3ptr[i][1],cycles3ptr[i][2]])
        print "cycles1=",cycles1
        print "cycles3=",cycles3
        if map_from>0:
            print "Want permutation mapping {0} to {1}".format(map_from,map_to)
    dd1 = int(max(d1,1)); dd3=int(max(d3,1))
    perms1 = MyPermutationIterator(dd1)
    perms3 = MyPermutationIterator(dd3)
    cdef CycleCombinationIterator CCIR
#    cdef list three_cycles_R    
    CCIR = CycleCombinationIterator(pR)
    nir = CCIR.length()
    cdef int**  three_cycles_R_ptr = NULL
    three_cycles_R_ptr = <int**> sage_malloc(sizeof(int*)*(nir))
    if three_cycles_R_ptr == NULL:
        raise MemoryError
    for i in range(nir):
        three_cycles_R_ptr[i] = NULL
        three_cycles_R_ptr[i] = <int*> sage_malloc(sizeof(int)*mu)
        if three_cycles_R_ptr[i] == NULL:
            raise MemoryError
        CCIR.permutation_nr_c_ptr(i,three_cycles_R_ptr[i])      
#    three_cycles_R = CCIR.list()
    if verbose>0:
        #print "three_cycles=",three_cycles_R
        for i in range(nir):
            for j in range(mu):
                print "three_cycles_R[{0}][{1}]={2}".format(i,j,three_cycles_R_ptr[i][j])
#    cdef int len_l0
#    l0 = [0 for i in range(mu)]
#    l = [0 for i in range(mu)]
    perm = MyPermutation(length=mu,check=0,init=0)
    for p1 in perms1:
        if verbose>0:
            print "perm1=",p1
            print "cy1=",cycles1
        if num_one_cycles>0:
            for i in range(num_one_cycles):
                #one_cycles = p1(cycles1) #permute_list(cycles[1],p1)
                #ll[i] = one_cycles[i][0]
                ll[i]=cycles1ptr[p1._entries[i]-1]
        for p3 in perms3:            
            #l = copy(l0)
            if num_three_cycles>0:
                #three_cycles = p3(cycles3) #permute_list(cycles[3],p3)
                for i in range(num_three_cycles):
                    j = p3._entries[i]-1
                    ll[num_one_cycles+i*3+0]=cycles3ptr[j][0]
                    ll[num_one_cycles+i*3+1]=cycles3ptr[j][1]
                    ll[num_one_cycles+i*3+2]=cycles3ptr[j][2]
            perm.set_entries(ll)
            #ll = flatten_list2d(l)
            #p = get_conjugating_perm_list(l,lR)
            if verbose>1:
                for i in range(mu):
                    print "ll=",i,ll[i]
                    print "llR=",i,llR[i]
            p = get_conjugating_perm_ptr_unsafe(mu,ll,llR)
            if verbose>0:
                #print "l=",l
                print "perm3=",p3
                print "perm0=",p
            ## We now act with all combinations of cycles from R
            for ir in range(nir):
                #ptmp = <MyPermutation> three_cycles_R[ir]
#                ptmp = CCIR.permutation_nr(ir)
                do_cont = 0
#                pp = p._mult_perm(CCIR.permutation_nr(ir))
#                _mult_perm_unsafe(mu,p._entries,(<MyPermutation>three_cycles_R[ir])._entries,pp._entries)
                _mult_perm_unsafe(mu,p._entries,three_cycles_R_ptr[ir],pp._entries)

                #pp = p._mult_perm(three_cycles_R[ir])
                if verbose>0:
                    ptmp = CCIR.permutation_nr_c(ir)
                    print "ptmp=",ptmp
                    print "pp=",pp
                # if S=pSp^-1 and Sx=x then S(px)=pS(x)=px
                for i in range(mu):
                     j = pp._entries[i]
                     if pS._entries[i]==i+1:
                         if pS1._entries[j-1]<>j:
                             do_cont = 1
                             if verbose>0:
                                 print "Fails at preserving fixed points!"
                     else:
                         if pS1._entries[j-1]==j:
                             do_cont = 1   #if pp(x) not in nonfixedS1:
                             if verbose>0:
                                 print "Fails at preserving non-fixed points!"
                     if do_cont==1:
                         break
                if do_cont==1:
                    continue
                # If we are here we may have a conjugating map
                if map_from > 0 and pp._entries[map_from-1]<>map_to:
                    if verbose>0:
                        print "Fails at mapping desired points!"
                    continue
                pScon = pS._conjugate(pp)
                if verbose>0:
                    print "pS^pp=",pScon
                if pScon.eq_c(pS1)==1:
                    t[0] = 1
                    #p_out = pp
                    if verbose>0:
                        print "p_out = ",pp
                    if ll <> NULL:
                        sage_free(ll)
                    if llR <> NULL:
                        sage_free(llR)
                    if cycles1ptr<>NULL:
                        sage_free(cycles1ptr)
                    if cycles3ptr<>NULL:
                        for i in range(num_three_cycles):
                            if cycles3ptr[i]<>NULL:
                                sage_free(cycles3ptr[i])
                        sage_free(cycles3ptr)
                    return pp
            ## We only consider permutations which fixes S
        ## We assume that p ionly has cycles of length 1 or 3 (i.e. it has order 3)
        #perms[i]=MyPermutation(i)
    if verbose>0:
        print "t=",t[0]
    if ll <> NULL:
        sage_free(ll)
    if llR <> NULL:
        sage_free(llR)
    if cycles1ptr<>NULL:
        sage_free(cycles1ptr)    
    if cycles3ptr<>NULL:
        for i in range(num_three_cycles):
            if cycles3ptr[i]<>NULL:
                sage_free(cycles3ptr[i])
        sage_free(cycles3ptr)
    if three_cycles_R_ptr<>NULL:
        for i in range(nir):
            if three_cycles_R_ptr[i] <> NULL:
                sage_free(three_cycles_R_ptr[i])
        sage_free(three_cycles_R_ptr)

cpdef list flatten_list2d(list l):
    r"""
    Takes a list of list and flattens it.
    """
    cdef int i,n
    n = len(l)
    cdef list res
    res = copy(l[0])
    for i in range(1,n):
        res.extend(l[i])
    return res


def sort_perms2(a,b):
    r"""
    Sort a list of pairs of permutations.
    """
    if a<>b:
        sa=str(a[0])
        sb=str(b[0])
        ra=str(a[0])
        rb=str(b[0])
        if sa<sb:
            return -1
        elif sa==sb:
            if ra<rb:
                return -1
            else:
                return 1
        else:
            return 1
    else:
        return 0
