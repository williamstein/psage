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
Helper algorithms for enumerating and classifying subgroups of the modular group in terms of pairs of permutations.

CLASSES:


AUTHOR:

 - Fredrik Stroemberg



"""

include "sage/ext/interrupt.pxi" 
include "sage/ext/stdsage.pxi"  
include "sage/ext/cdefs.pxi"

from permutation_alg cimport MyPermutation,MyPermutationIterator
from permutation_alg cimport print_vec,_conjugate_perm,_are_eq_vec,transposition,_mult_perm,are_transitive_perm_c,perm_to_cycle_c

from psage.modform.maass.mysubgroup import MySubgroup

from sage.modules.vector_integer_dense cimport Vector_integer_dense

from permutation_alg import verbosity,print_cycles,perm_to_cycle,sort_perms
from sage.all import deepcopy,copy,ZZ,vector,subsets
from sage.groups.perm_gps.permgroup import PermutationGroup
from sage.groups.perm_gps.permgroup_named import SymmetricGroup

from sage.combinat.combinat import tuples
from sage.all import SL2Z
from sage.interfaces.all import gap

from time import clock, time

cpdef list_all_admissable_pairs(sig,int get_details=1,int verbose=0,int get_one_rep=0,int congruence=-1,int do_new=0):
    r"""
    List all possible pairs (up to conjugacy) of admissible permutations E,R
    correcponsing to groups G with signature = sig
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
    cdef MyPermutation rss,pres,rss_tmp,ee,r,rr,epp,rpp,ep,rp,eppp,rppp

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
    cdef MyPermutation Sp,Tp,Spc
    Sp = MyPermutation(Sl)
    ### Note that there are only two choices of S which are inequivalent modulo permutations fixing 1
    ### S = (1)...(e2)(1+e2 2+e2)...  or
    ### S = (1 2)(3)...(1+e2)(2+e2 3+e2)...  

    rss=MyPermutation(length=mu)
    pres=MyPermutation(length=mu)
    ee=MyPermutation(length=mu)
    r=MyPermutation(length=mu)
    rr=MyPermutation(length=mu)
    ep=MyPermutation(length=mu)
    rp=MyPermutation(length=mu)
    epp=MyPermutation(length=mu)
    rpp=MyPermutation(length=mu)
    eppp=MyPermutation(length=mu)
    rppp=MyPermutation(length=mu)
    cdef int a,b,c
    if verbose>0:
        print "signature=",sig
    ## Then we have to find matching order 3 permutations R wih e3 fixed points
    ## Without loss of generality we can choose those fixed points as either 
    ## e2+1, e2+3,...,e2+2*e3-1  or
    ## 1,e2+1, e2+3,...,e2+2*e3-1
    cdef list rfx_list
    rfx_list=list()
    rfx_list=range(e2+1,e2+2*e3,2)
    cdef int* rfx
    rfx = <int*>sage_malloc(sizeof(int)*e3)
    for i from 0<= i < e3:
        if i < len(rfx_list):
            rfx[i]=rfx_list[i]
        else:
            rfx[i]=0
    if verbose>0:
        print "fixed pts for R=",print_vec(e3,rfx)
    if len(rfx_list)<>e3:
        raise ValueError, "Did not get correct number of fixed points!"
    cdef list list_of_R,list_of_Rs
    #print "rfx=",print_vec(e3,rfx)
    #sprint "rfx_list=",rfx_list
    #cdef dict list_of_T
    list_of_R=[]
    list_of_Rs=[]
    
    #list_of_T=dict()
    if verbosity(verbose,3):
        mpi_verbose= verbose % 8
    else:
        mpi_verbose=0
    if verbose>0:
        print "mpi_verbose=",mpi_verbose
    PRI = MyPermutationIterator(mu,order=3,fixed_pts=rfx_list,verbose=mpi_verbose)
    if verbose>0:
        print "PRI.list=",printf("%p ", PRI._list_of_perms)
    #for R in P.filter(lambda x: x.fixed_points()==rfx):
    max_num = PRI.max_num()
    if verbose>0:
        print "max. num of possible R=",max_num
        print "original fixedpts=",PRI.fixed_pts()
    cdef MyPermutation pR,p
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
    rcycle_lens = <int*>sage_malloc(sizeof(int)*mu)
    if not rcycle_lens: raise MemoryError
    Rcycles = <int*>sage_malloc(sizeof(int*)*mu*mu)
    if not Rcycles: raise MemoryError
    gotten = <int *>sage_malloc(sizeof(int)*mu)
    if not gotten: raise MemoryError
    sig_on()
    for pR in PRI: #ii from 0<= ii <max_num:
        #_to_cycles2(mu,pR._entries,Rcycles,rcycle_lens,num_rcycles)
        # we might also make some apriori requirements on R
        # 1) R can not contain all fixed points of E in one cycle
        if verbose>0:
            print "S=",Sp.to_cycles() #print_vec(mu,Sptr)
            print "R=",pR.to_cycles() #print_vec(mu,<int *>PRI._current_perm)
            #print "pr=",pR
        if verbose>0:
            print "Checking transitivity!"
        #if not are_transitive_permutations(Sl,pR,mu):
        if not are_transitive_perm_c(<int*>Sp._entries,<int*>pR._entries,gotten,mu,mpi_verbose):
            #sage_free(gotten)
            continue
        #sage_free(gotten)        
        if verbose>0:
            print "Checking the number of cusps!"
        _mult_perm(mu,Sptr,<int *>pR._entries,Tptr)
        #T=E*R
        Tcyc=perm_to_cycle_c(mu,Tptr)
        if verbose>0:
            print "Tp=",Tcyc
            print "current fixedpts=",PRI.fixed_pts()
            print "number of cusps=",len(Tcyc)
        if len(Tcyc)<>h:
            if verbose>0:
                print "Incorrect number of cusps. remove!"
                print "next!"
            continue
        if verbose>0:
            #print "ursed=",print_vec(mu,used)
            print "current fixedpts1=",PRI.fixed_pts()
            print "rfx=",print_vec(e3,rfx)
        if get_one_rep:
            if congruence<>-1:
                G=MySubgroup(o2=Sp,o3=pR)
                if G.is_congruence() and congruence==1:
                    return Sp,pR
                elif not G.is_congruence() and congruence==0:
                    return Sp,pR
            else:
                return Sp,pR
            continue
        ## If we are making a complete list then 
        ## we will now try to pick representatives for R modulo permutations fixing 1
        ## by elimininating as many possibilities as possible
        # Recall that we chose:
        # S = (1)...(e2)(e2+1 e2+2)(e2+3 e2+4)...(mu-1 mu)
        # R = (e2+1)(e2+3)...(e2+2*e3-1)(a b c)...
        if used==NULL:
            used = <int *>sage_malloc(sizeof(int)*mu)
        for x from 0 <= x < mu: 
            used[x]=0
        for i from 0<= i < e3:
            used[rfx[i]-1]=1            
        rcycle_beg=0
        Rcycles_list=perm_to_cycle_c(mu,<int*>pR._entries)
        do_cont = 0
        for rc in Rcycles_list:
            if len(rc)<>3:
                continue
            # only the 3-cycles are relevant here (not the 1-cycles)
            # cy = (a b c)
            a=rc[0]; b=rc[1]; c=rc[2]
            if verbose>0:
                print "a,b,c=",a,b,c
            used[a-1]=1 # We have fixed a
            if verbose>0:
                print "used=",print_vec(mu,used)
            ## If b and c are equivalent with respect to conjugation by a perm. p which preserves S, i.e. if
            ##   i) S(b)=b and S(c)=c, or
            ##  ii) (b b') and (c c') are two cycles of S with c and c' not used previously and not fixed points of R
            ## then we choose b < c
            #if (b<=e2 and c<=e2) or ((b>=end_fc and used[Sptr[b-1]-1]==0) and (c>=end_fc and used[Sptr[c-1]-1]==0)):
            if equivalent_integers_mod_fixS(b,c,mu,e2,end_fc,Sptr,used)==1:
                if verbose>0:
                    print b," (b c) and ",c," are equivalent in",rc
                if b>c:
                    if verbose>0:
                        print "remove (b>c)!"
                    do_cont = 1
                    break
            for j from a+1<= j < b:  # I want to see if there is a j, equivalent to b, smaller than b
                #t1 = (b<=e2 and j<=e2)  # if j and b are both fixed by S
                #t2 = (b>=end_fc and used[Sptr[b-1]-1]==0 and used[b-1]==0)
                #t3 = (j>=end_fc and used[Sptr[j-1]-1]==0 and used[j-1]==0)
                
                #if (b<=e2 and j<=e2) or ((b>end_fc and used[Sptr[b-1]-1]==0 and used[b-1]==0) and (j>end_fc and used[Sptr[j-1]-1]==0 and used[j-1]==0)):
               if equivalent_integers_mod_fixS(b,j,mu,e2,end_fc,Sptr,used)==1:
                    #if t1 or (t2 and t3):
                    if verbose>0:
                        print j," (b) and ",b," are equivalent in",rc
                    if j<b:
                        if verbose>0:
                            print "remove!"
                        do_cont = 1
                        break #raise StopIteration()
            if do_cont==1:
                if verbose>0:
                    print "breaking0!"
                break
            used[b-1]=1
            if verbose>0:
                print "used=",print_vec(mu,used)
            for j from a+1 <= j < c: # I want to see if there is a j, equivalent to c, smaller than c and which is not used
                if (c<=e2 and j<=e2 and used[j-1]==0) or ((c>e2+e3+1 and used[Sptr[c-1]-1]==0 and used[c-1]==0) and (j>e2+e3+1 and used[Sptr[j-1]-1]==0 and used[j-1]==0)):
                    if verbose>0:
                        print j," (c) and ",c," are equivalent in",rc
                    if j<c:
                        if verbose>0:
                            print "remove!"
                        do_cont = 1 #raise StopIteration()
                        break
            if do_cont==1:
                if verbose>0:
                    print "breaking1!"
                break
            used[c-1]=1
            if verbose>0:
                print "used=",print_vec(mu,used)
        if do_cont==1:
            if verbose>0:
                print "next and continue!"
            continue
        else:
            pass
        ## If we are here, R is a true candidate.
        list_of_R.append(pR)
        if verbose>0:
            print "added pR=",pR
            print "current list of Rs and Ts:"
            for r in list_of_R:
                #print ":::::::::::::::: ",r.cycles(),";",(Sp*r).cycles()
                print ":::::::::::::::: ",r.to_cycles(),";",(Sp*r).to_cycles()
    if gotten<>NULL:
        sage_free(gotten)
    if get_one_rep:
        if congruence==1:
            print "We didn't find any congruence subgroups..."
            return []
        elif congruence==0:
            print "We didn't find any non-congruence subgroups..."
            return []
        else:
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
    #print "after deallPRI.list=",printf("%p ", PRI._list_of_perms)
    # Now list_of_R contains at least one representative for each group with the correct signature
    sig_off()
    if verbose>0:
        #print "Original list of R=",list_of_R
        print "Original list of R="
        for i from 0 <= i < len(list_of_R):
            print perm_to_cycle(list_of_R[i])
    list_of_R.sort(cmp=sort_perms)
    # list_of_R  will contain a list of representatives for groups
    # but might contain more than one representative for each group.
    # Therefore we now filter away the identical groups in PSL(2,Z), i.e. mod 1
    cdef list list_of_R_tmp=[]
    cdef int nlr = len(list_of_R)
    if verbose>=0:
        start = time()

    if len(list_of_R)>1:
        if do_new==1:
            list_of_R_tmp=filter_list_mod_1_mod_S_new(list_of_R,mu,e2,Sp,verbose)
        else:
            list_of_R_tmp=filter_list_mod_1_mod_S(list_of_R,mu,e2,Sp,verbose)
        list_of_R=list_of_R_tmp
    if verbose>=0:    
        print "Time for first filter= ",time()-start
        print "Preliminary list of DIFFERENT subgroups in PSL(2,Z):" #,map(lambda x:S(x),list_of_R)
        print "(Note: these may or may not be PSL or PGL conjugate)"
        for i from 0 <= i < len(list_of_R):
            print perm_to_cycle(list_of_R[i])    
    list_of_groups=copy(list_of_R)
    indicator_list=range(1,len(list_of_R)+1)
    list_of_R_tmp=[]
    cdef dict lc_psl,lc_pgl,lc_psl_maps,lc_pgl_maps
    lc_psl=dict()      # list of conjugates
    lc_psl_maps=dict() # maps between conjugates
    if verbose>=0:
        start = time()
    sig_on()
    conjugates,conjugate_maps,Rmodpsl,Rmodpgl=filter_list_mod_S(list_of_R,mu,e2,e3,Sp,verbose)
    sig_off()
#    lc_psl,lc_psl_maps,list_of_R=filter_list_mod_psl_mod_S_new(list_of_R,mu,e2,e3,Sp,verbose)
    ## in the list of PGL-conjugates we want to only have representatives
    ## of classes modulo PSL
    if verbose>=0:    
        print "Time for second filter= ",time()-start
        print "List of conjugacy classes mod PGL(2,Z):" 
        for R in Rmodpgl:
            if conjugates.has_key(R):
                t = conjugates[R]['psl']
                if t<>[]:
                    print t
    ## We finally want to add a list of reflected groups
    reflections={}
    for R in Rmodpsl:
        reflections[R]={}
        Rs = R.square().conjugate(Sp)
        conj_pgl = conjugates[R]['pgl']
        conj_psl = conjugates[R]['psl']
        if R in conj_pgl:
            t,A = are_conjugate_pairs_of_perms(Sp,Rs,Sp,R)
            reflections[R]={'group':R,'map':A}
        elif conj_pgl<>[]:
            Rtest = conjugates[R]['pgl'][0]
            if Rtest in conj_psl:  ## G^* is conjugate to G
                p = conjugate_maps[R]['psl'][0]
                t,A = are_conjugate_pairs_of_perms(Sp,Rs,Sp,R)
                if t==1:
                    reflections[R]={'group':R,'map':A}
            else:
                t,A = are_conjugate_pairs_of_perms(Sp,Rs,Sp,Rtest)
                if t==1:
                    reflections[R]={'group':Rtest,'map':A}
        if reflections[R]=={}:
            for Rtest in Rmodpsl:
                t,A = are_conjugate_pairs_of_perms(Sp,Rs,Sp,Rtest)
                if t==1:
                    reflections[R]={'group':Rtest,'map':A}
    ##  Modulo *:  E->E, R->ER^2E
#    lc_pgl,lc_pgl_maps,list_of_R=filter_list_mod_psl_mod_S_new(lc_psl.keys(),mu,e2,e3,Sp,verbose,do_pgl=1)
    d = dict()
    d['sig']=sig
    d['S']=Sp
    d['numg']=len(list_of_groups)
    d['groups']=list_of_groups
    d['groups_mod_psl']=Rmodpsl
    d['groups_mod_pgl']=Rmodpgl
    d['conjugates']=conjugates
    d['conjugate_maps']=conjugate_maps
    d['reflections']=reflections
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



cpdef list_all_admissable_pairs_new(sig,int get_details=1,int verbose=0,int get_one_rep=0,int congruence=-1,int do_new=0):
    r"""
    List all possible pairs (up to conjugacy) of admissible permutations E,R
    correcponsing to groups G with signature = sig
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
    cdef MyPermutation rss,pres,rss_tmp,ee,r,rr,epp,rpp,ep,rp,eppp,rppp

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
    cdef MyPermutation Sp,Tp,Spc
    Sp = MyPermutation(Sl)
    ### Note that there are only two choices of S which are inequivalent modulo permutations fixing 1
    ### S = (1)...(e2)(1+e2 2+e2)...  or
    ### S = (1 2)(3)...(1+e2)(2+e2 3+e2)...  
    ### In the first stage we are only concerned about getting PSL(2,Z) conjugacy classes so we fix S of the first form.
    rss=MyPermutation(length=mu)
    pres=MyPermutation(length=mu)
    ee=MyPermutation(length=mu)
    r=MyPermutation(length=mu)
    rr=MyPermutation(length=mu)
    ep=MyPermutation(length=mu)
    rp=MyPermutation(length=mu)
    epp=MyPermutation(length=mu)
    rpp=MyPermutation(length=mu)
    eppp=MyPermutation(length=mu)
    rppp=MyPermutation(length=mu)
    cdef int a,b,c
    if verbose>0:
        print "signature=",sig
    ## Then we have to find matching order 3 permutations R wih e3 fixed points
    ## Without loss of generality, up to PSL(2,Z) conjugacy we can choose those fixed points as
    ## e2+1, e2+3,...,e2+2*e3-1 
    cdef list rfx_list
    rfx_list=list()
    rfx_list=range(e2+1,e2+2*e3,2)
    cdef int* rfx
    rfx = <int*>sage_malloc(sizeof(int)*e3)
    for i from 0<= i < e3:
        if i < len(rfx_list):
            rfx[i]=rfx_list[i]
        else:
            rfx[i]=0
    if verbose>0:
        print "fixed pts for R=",print_vec(e3,rfx)
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
        print "mpi_verbose=",mpi_verbose
    PRI = MyPermutationIterator(mu,order=3,fixed_pts=rfx_list,verbose=mpi_verbose)
    if verbose>0:
        print "PRI.list=",printf("%p ", PRI._list_of_perms)
    #for R in P.filter(lambda x: x.fixed_points()==rfx):
    max_num = PRI.max_num()
    if verbose>0:
        print "max. num of possible R=",max_num
        print "original fixedpts=",PRI.fixed_pts()
    cdef MyPermutation pR,p,Rp
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
    rcycle_lens = <int*>sage_malloc(sizeof(int)*mu)
    if not rcycle_lens: raise MemoryError
    Rcycles = <int*>sage_malloc(sizeof(int*)*mu*mu)
    if not Rcycles: raise MemoryError
    gotten = <int *>sage_malloc(sizeof(int)*mu)
    if not gotten: raise MemoryError
    sig_on()
    for pR in PRI: 
        if verbose>0:
            print "S=",Sp.to_cycles() #print_vec(mu,Sptr)
            print "R=",pR.to_cycles() #print_vec(mu,<int *>PRI._current_perm)
        if verbose>0:
            print "Checking transitivity!"
        if not are_transitive_perm_c(<int*>Sp._entries,<int*>pR._entries,gotten,mu,mpi_verbose):
            continue
        if verbose>0:
            print "Checking the number of cusps!"
        _mult_perm(mu,Sptr,<int *>pR._entries,Tptr)
        #T=E*R
        Tcyc=perm_to_cycle_c(mu,Tptr)
        if verbose>0:
            print "Tp=",Tcyc
            print "current fixedpts=",PRI.fixed_pts()
            print "number of cusps=",len(Tcyc)
        if len(Tcyc)<>h:
            if verbose>0:
                print "Incorrect number of cusps. remove!"
                print "next!"
            continue
        if verbose>0:
            print "current fixedpts1=",PRI.fixed_pts()
            print "rfx=",print_vec(e3,rfx)
        if get_one_rep:
            if congruence<>-1:
                G=MySubgroup(o2=Sp,o3=pR)
                if G.is_congruence() and congruence==1:
                    return Sp,pR
                elif not G.is_congruence() and congruence==0:
                    return Sp,pR
            else:
                return Sp,pR
            continue
        ## We will now further shorten the list by trying to pick representatives for R modulo permutations fixing 1
        # We will do this by choosing the lexicographically 'smallest' possible 3-cycles
        # Recall that we have chosen:
        # S = (1)...(e2)(e2+1 e2+2)(e2+3 e2+4)...(mu-1 mu)
        # R = (e2+1)(e2+3)...(e2+2*e3-1)(a b c)...
        if used==NULL:
            used = <int *>sage_malloc(sizeof(int)*mu)
        for x from 0 <= x < mu: 
            used[x]=0
        for i from 0<= i < e3:
            used[rfx[i]-1]=1            
        rcycle_beg=0
        Rcycles_list=perm_to_cycle_c(mu,<int*>pR._entries)
        do_cont = 0
        for rc in Rcycles_list:
            if len(rc)<>3:
                continue
            # only the 3-cycles are relevant here (not the 1-cycles)
            # cy = (a b c)
            a=rc[0]; b=rc[1]; c=rc[2]
            if verbose>0:
                print "a,b,c=",a,b,c
            used[a-1]=1 # We have fixed a
            if verbose>0:
                print "used=",print_vec(mu,used)
            ## If b and c are equivalent with respect to conjugation by a perm. p which preserves S, i.e. if
            ##   i) S(b)=b and S(c)=c, or
            ##  ii) (b b') and (c c') are two cycles of S with c and c' not used previously and not fixed points of R
            ## then we choose b < c
            #if (b<=e2 and c<=e2) or ((b>=end_fc and used[Sptr[b-1]-1]==0) and (c>=end_fc and used[Sptr[c-1]-1]==0)):
            if equivalent_integers_mod_fixS(b,c,mu,e2,end_fc,Sptr,used)==1:
                if verbose>0:
                    print b," (b c) and ",c," are equivalent in",rc
                if b>c:
                    if verbose>0:
                        print "remove (b>c)!"
                    do_cont = 1
                    break
            for j from a+1<= j < b:  # I want to see if there is a j, equivalent to b, smaller than b
               if equivalent_integers_mod_fixS(b,j,mu,e2,end_fc,Sptr,used)==1:
                    #if t1 or (t2 and t3):
                    if verbose>0:
                        print j," (b) and ",b," are equivalent in",rc
                    if j<b:
                        if verbose>0:
                            print "remove!"
                        do_cont = 1
                        break #raise StopIteration()
            if do_cont==1:
                if verbose>0:
                    print "breaking0!"
                break
            used[b-1]=1
            if verbose>0:
                print "used=",print_vec(mu,used)
            for j from a+1 <= j < c: # I want to see if there is a j, equivalent to c, smaller than c and which is not used
                if (c<=e2 and j<=e2 and used[j-1]==0) or ((c>e2+e3+1 and used[Sptr[c-1]-1]==0 and used[c-1]==0) and (j>e2+e3+1 and used[Sptr[j-1]-1]==0 and used[j-1]==0)):
                    if verbose>0:
                        print j," (c) and ",c," are equivalent in",rc
                    if j<c:
                        if verbose>0:
                            print "remove!"
                        do_cont = 1 #raise StopIteration()
                        break
            if do_cont==1:
                if verbose>0:
                    print "breaking1!"
                break
            used[c-1]=1
            if verbose>0:
                print "used=",print_vec(mu,used)
        if do_cont==1:
            if verbose>0:
                print "next and continue!"
            continue
        else:
            pass
        ## If we are here, R is a true candidate.
        list_of_R.append(pR)
        if verbose>0:
            print "added pR=",pR
            print "current list of Rs and Ts:"
            for r in list_of_R:
                print ":::::::::::::::: ",r.to_cycles(),";",(Sp*r).to_cycles()
    if gotten<>NULL:
        sage_free(gotten)
    if get_one_rep:
        if congruence==1:
            print "We didn't find any congruence subgroups..."
            return []
        elif congruence==0:
            print "We didn't find any non-congruence subgroups..."
            return []
        else:
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
    #print "after deallPRI.list=",printf("%p ", PRI._list_of_perms)
    # Now list_of_R contains at least one representative for each group with the correct signature
    sig_off()
    if verbose>0:
        #print "Original list of R=",list_of_R
        print "Original list of R="
        for i from 0 <= i < len(list_of_R):
            print perm_to_cycle(list_of_R[i])
    # Then add conjugates mod PSL(2,Z) to get all groups before filtering.
    cdef list list_of_groups=[]
    for R in list_of_R:
        list_of_groups.append((Sp,R))
        for j in range(2,mu):
            Rp = copy(R); Spc = copy(Sp)
            Rp.conjugate_with_transposition(1,j)
            Spc.conjugate_with_transposition(1,j)
            if (Spc,Rp) not in list_of_groups:
                list_of_groups.append((Spc,Rp))
    list_of_groups.sort(cmp=sort_perms2)
    # list_of_R  will contain a list of representatives for groups
    # but might contain more than one representative for each group.
    # Therefore we now filter away the identical groups in PSL(2,Z), i.e. mod 1
    if verbose>=0:
        print "Preliminary list of groups:"
        for t in list_of_groups:
            print t
    cdef list list_of_groups_tmp=[]
    cdef int nlr = len(list_of_R)
    if verbose>=0:
        start = time()

    if len(list_of_groups)>1:
        list_of_groups_tmp=filter_list_of_pairs_mod_1_mod_S_new(list_of_groups,mu,verbose)
    #list_of_R_tmp=filter_list_mod_1_mod_S(list_of_R,mu,e2,Sp,verbose)
        list_of_groups=list_of_groups_tmp
    if verbose>=0:    
        print "Time for first filter= ",time()-start
        print "Preliminary list of DIFFERENT subgroups in PSL(2,Z):" #,map(lambda x:S(x),list_of_R)
        print "(Note: these may or may not be PSL or PGL conjugate)"
        for i from 0 <= i < len(list_of_groups):
            print list_of_groups[i]
        for i from 0 <= i < len(list_of_groups):            
            Sp = list_of_groups[i][0]
            R = list_of_groups[i][1]
            print "S=",Sp.to_cycles(),R.to_cycles()
    indicator_list=range(1,len(list_of_R)+1)
    list_of_R_tmp=[]
    cdef dict lc_psl,lc_pgl,lc_psl_maps,lc_pgl_maps
    lc_psl=dict()      # list of conjugates
    lc_psl_maps=dict() # maps between conjugates
    if verbose>=0:
        start = time()
    sig_on()
    conjugates,conjugate_maps,Gmodpsl,Gmodpgl=find_conjugate_pairs(list_of_groups,mu,verbose)
    sig_off()
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
            t,A = are_conjugate_pairs_of_perms(S,Rs,S1,R1)
            if t==1:
                reflections[(S,R)]={'group':(S1,R1),'map':A}
                do_cont = 1
                break
        if do_cont==1:
            continue
        for S1,R1 in conj_pgl:
            t,A = are_conjugate_pairs_of_perms(S,Rs,S1,R1)
            if t==1:
                reflections[(S,R)]={'group':(S1,R1),'map':A}
                do_cont = 1
                break
        if do_cont==1:
            continue
        for Stest,Rtest in Gmodpsl:
            t,A = are_conjugate_pairs_of_perms(S,Rs,Stest,Rtest)
            if t==1:
                reflections[(S,R)]={'group':(Stest,Rtest),'map':A}
                break
            
    ##  Modulo *:  E->E, R->ER^2E
#    lc_pgl,lc_pgl_maps,list_of_R=filter_list_mod_psl_mod_S_new(lc_psl.keys(),mu,e2,e3,Sp,verbose,do_pgl=1)
    d = dict()
    d['sig']=sig
    d['S']=Sp
    d['numg']=len(list_of_groups)
    d['groups']=list_of_groups
    d['groups_mod_psl']=Gmodpsl
    d['groups_mod_pgl']=Gmodpgl
    d['conjugates']=conjugates
    d['conjugate_maps']=conjugate_maps
    d['reflections']=reflections
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
    if a<=e2 and b<=e2 and used[b-1]==0:
        return 1
    if (a>=end_fix and used[Sptr[a-1]-1]==0) and (b>=end_fix and used[Sptr[b-1]-1]==0):
        return 1
    return 0
             




cpdef tuple filter_list_mod_S(list listRin, int mu, int e2,int e3,MyPermutation S,int verbose=0,int mpi_verbose=0,int do_pgl=0):
    r""" Removes duplicates in listR modulo permutations of the form (1j) which keeps S invariant.

    The purpose of this step is to reduce the list we later have to check completely.
    If do_pgl = 1 we reduce modulo PGL(2,Z)/PSL(2,Z) instead, i.e. we check R ~ ER^2E 
    """
    cdef list Rmodpsl=[],Rmodpgl=[]
    cdef dict conjugates={},conjugate_maps={}
    cdef MyPermutation R,R1,Sc,Rpsl,Rpgl,Scp,pp,pi
    cdef int numr = len(listRin)
    cdef int* checked=NULL
    cdef int t=0
    pp = MyPermutation(length=mu,rep=0)
    checked = <int*>sage_malloc(numr*sizeof(int))
    for i in range(numr):
        checked[i]=0
    if verbose>0:
        print "R's:",listRin
    for i in range(0,numr):
        R = listRin[i]
        conjugates[R]={'psl':[R],'pgl':[]}
        conjugate_maps[R]={'psl':[pp],'pgl':[]}
        
    for i in range(0,numr):
        #if checked[i]==1:
        #    continue
        R = listRin[i]
        if verbose>0:
            print "Test R[{0}]={1}".format(i,R)
        # checked[i] = 0 if R[i] is not (yet) found to be a conjugate
        # checked[i] = 1 if R[i] is conjugate mod PSL but not (yet) PGL\PSL
        # checked[i] = 2 if R[i] is conjugate mod PGL\PSL but not (yet) PSL
        # checked[i] = 3 if R[i] is conjugate mod PGL\PSL and PSL                           
        if checked[i]<2:
            Rmodpgl.append(R)
        if checked[i] in [0,2]:
            Rmodpsl.append(R)
        if verbose>0:
            for j in range(numr):
                print "checked[{0}]={1}".format(j,checked[j])
        for j in range(i,numr):
            if checked[j]==3:
                continue
            Rpsl = listRin[j]
            Rpgl = listRin[j].square().conjugate(S)
            p = are_conjugate_perm(R,Rpsl)
            if verbose>0:
                print "------------------------------------------------------------------"
                print "Comp PSL R[{0}]={1}".format(j,Rpsl)
                print "p=",p.to_cycles()
            Scp = S.conjugate(p)
            if verbose>0:
                print "R=",R.to_cycles()
                print "S=",S.to_cycles()
                print "S^p=",Scp.to_cycles()
            pp=are_conjugate_wrt_stabiliser(Rpsl,Scp,S,p,&t)
            if verbose>0:
                print "t=",t
                print "j=",j
            if t==1:
                if verbose>0:                    
                    print "p=",p.to_cycles()
                    print "pp=",pp.to_cycles()
                pp = p*pp
                Sc = S.conjugate(pp)
                if verbose>0:                    
                    print "pp*p=",pp.to_cycles()
                    print "S^pp=",Sc.to_cycles()
                
                if Sc==S: ## the pair is conjugate
                    if listRin[j] not in conjugates[R]['psl']:
                        conjugates[R]['psl'].append(listRin[j])
                        pp.set_rep(0)
                        conjugate_maps[R]['psl'].append(pp)
                    checked[j]+=1
                    if R not in conjugates[listRin[j]]['psl']:
                        conjugates[listRin[j]]['psl'].append(R)
                        pi = pp.inverse(); pi.set_rep(0)
                        conjugate_maps[listRin[j]]['psl'].append(pi)
            if checked[j]<2:
                # Check PGL(2,Z)
                p = are_conjugate_perm(R,Rpgl)
                if verbose>0:
                    print "------------------------------------------------------------------"
                    print "Comp PGL R[{0}]={1}".format(j,Rpgl)
                    print "p=",p.to_cycles()
                Scp = S.conjugate(p)
                if verbose>1 or (i==2 and j==2):
                    print "R=",R.to_cycles()
                    print "S=",S.to_cycles()
                    print "S^p=",Scp.to_cycles()
                    pp=are_conjugate_wrt_stabiliser(Rpgl,Scp,S,p,&t,1)
                else:
                    pp=are_conjugate_wrt_stabiliser(Rpgl,Scp,S,p,&t)
                if t==1:
                    pp = p*pp
                    if verbose>0:                    
                        print "pp=",pp.to_cycles()
                    Sc = S.conjugate(pp)
                    if Sc==S: ## the pair is conjugate
                        conjugates[R]['pgl'].append(listRin[j])
                        pp.set_rep(0)
                        conjugate_maps[R]['pgl'].append(pp)
                        checked[j]+=2
                        if i==0 and j==3:
                            print "R=",R
                            print "newlist=",listRin[j]
                            print "conjugates=",conjugates[listRin[j]]['pgl']
                        if R not in conjugates[listRin[j]]['pgl']:
                            conjugates[listRin[j]]['pgl'].append(R)
                            if i==0 and j==3:
                                print "Appending!"
                                print "conjugates=",conjugates

                            pi = pp.inverse(); pi.set_rep(0)
                            conjugate_maps[listRin[j]]['pgl'].append(pi)
    #if R in conjugates.keys():
    #    for R1 in conjugates[R]['psl']:
    ## in the list of PGL-conjugates we want to only have representatives
    ## of classes modulo PSL. Observe that prevsiously we have basically only checked
    ## conjugation modulo PGL / PSL = <J>
    for R in Rmodpgl:
        if R not in Rmodpsl:
            Rmodpgl.remove(R)
    for R in conjugate_maps.keys():
        testlist = copy(conjugates[R]['pgl'])
        for R1 in testlist:
            if R1 not in Rmodpsl:
                i = conjugates[R]['pgl'].index(R1)
                conjugates[R]['pgl'].remove(R1)
                conjugate_maps[R]['pgl'].pop(i)
                
    if checked<>NULL:
        sage_free(checked)
    return conjugates,conjugate_maps,Rmodpsl,Rmodpgl


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
                print " {0} and {1} are not conjugate!".format(R.to_cycles(),Rpsl.to_cycles())
            if verbose>0:
                print "------------------------------------------------------------------"
                print "Comp PSL R[{0}]={1}".format(j,Rpsl)
                print "p=",p #.to_cycles()
            Scp = S.conjugate(p)
            if verbose>0:
                print "R=",R.to_cycles()
                print "S=",S.to_cycles()
                print "S^p=",Scp.to_cycles()
            pp=are_conjugate_wrt_stabiliser(Rpsl,Scp,Sc,p,&t)
            if verbose>0:
                print "t=",t
                print "j=",j
            if t==1:
                if verbose>0:                    
                    print "p=",p.to_cycles()
                    print "pp=",pp.to_cycles()
                pp = p*pp
                Scp = S.conjugate(pp)
                if verbose>0:                    
                    print "pp*p=",pp.to_cycles()
                    print "S^pp=",Sc.to_cycles()
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
                    print "p=",p.to_cycles()
                Scp = S.conjugate(p)
                if verbose>1 or (i==2 and j==2):
                    print "R=",R.to_cycles()
                    print "S=",S.to_cycles()
                    print "S^p=",Scp.to_cycles()
                    pp=are_conjugate_wrt_stabiliser(Rpgl,Scp,Sc,p,&t,1)
                else:
                    pp=are_conjugate_wrt_stabiliser(Rpgl,Scp,Sc,p,&t)
                if t==1:
                    pp = p*pp
                    if verbose>0:                    
                        print "pp=",pp.to_cycles()
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
        if S.num_fixed()>0 and S(1)<>1:
            i = Gmodpsl_tmp.index((S,R))
            for S1,R1 in conjugates[(S,R)]['pgl']:
                if S1.num_fixed()>0 and S1(1)==1:
                    Gmodpsl_tmp[i]=(S1,R1)
                    break
    Gmodpsl = Gmodpsl_tmp    
    Gmodpgl_tmp = copy(Gmodpgl)
    for S,R in Gmodpgl:
        if (S.num_fixed()>0 and S(1)<>1) or (S,R) not in Gmodpsl:
            i = Gmodpgl_tmp.index((S,R))
            for S1,R1 in conjugates[(S,R)]['pgl']:
                if (S1.num_fixed()==0 or S1(1)==1) and (S1,R1) in Gmodpsl:
                    Gmodpgl_tmp[i]=(S1,R1)
                    break
    Gmodpgl = Gmodpgl_tmp
    for grp in conjugate_maps.keys():        
        testlist = copy(conjugates[grp]['pgl'])
        for grp1 in testlist:
            if grp1 not in Gmodpsl:
                i = conjugates[grp]['pgl'].index(grp1)
                conjugates[grp]['pgl'].remove(grp1)
                conjugate_maps[grp]['pgl'].pop(i)
                
    if checked<>NULL:
        sage_free(checked)
    return conjugates,conjugate_maps,Gmodpsl,Gmodpgl



cpdef are_conjugate_groups(G1,G2,ret='SL2Z',coset_rep=1,check=0,verbose=0):
    r"""
    Determine whether G1 and G2 are conjugate in PSL(2,Z) and return either a permutation or a matrix A in SLZ which performs the conjugation, i.e. A^-1 G1 A = G2. 
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

cpdef are_conjugate_pairs_of_perms(MyPermutation S1,MyPermutation R1,MyPermutation S2,MyPermutation R2,str ret='SL2Z',int fix_one=0,int verbose=0):
    r"""

    fix_one = 1 if we want a permutation which fixes 1
    """

    #S1,R1 = pair1; S2,R2 = pair2
    ## First check if R1 and R2 are conjugate
    p = are_conjugate_perm(R1,R2)
    cdef MyPermutation pp,Sc
    if p==0:
        if ret<>'SL2Z':
            return 0,MyPermutation(length=S1.N())
        else:
            return 0,SL2Z([1,0,0,1])
    pp = <MyPermutation?>p
    Sc = S1.conjugate(pp)
    if verbose>0:
        print "p=",pp.to_cycles()
        print "S^p=",Sc.to_cycles()
    cdef int j,t=0
    if fix_one==1:
        j = pp.inverse()(1)
        pp=are_conjugate_wrt_stabiliser(R2,Sc,S2,pp,&t,j,verbose)
    else:
        pp=are_conjugate_wrt_stabiliser(R2,Sc,S2,pp,&t,verbose)
    if t == 0:
        if ret=='SL2Z':
            return 0, SL2Z([1,0,0,1])
        else:
            return 0, MyPermutation(length=S1.N())
    pp = p*pp
    Sc = S1.conjugate(pp)
    if Sc<>S2:
        raise ArithmeticError,"Conjugation did not work!"
    if ret<>'SL2Z':
        return 1,pp
    if pp.is_identity():
        return 1,SL2Z([1,0,0,1])
#    return 1,matrix_from_perm((S1,R1),pp,verbose)
    A = MySubgroup(o2=S1,o3=R1).coset_reps()[pp(1)-1].SL2Z()
    return 1,A
    

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
        print "gens=({0},{1})".format(S.to_cycles(),R.to_cycles())
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




#cpdef list filter_list_mod_1_mod_E(list listRout,list listRin, int mu, int e2,MyPermutation E,int verbose=0,int mpi_verbose=0):
cpdef list filter_list_mod_1_mod_S(list listRin, int mu, int e2,MyPermutation S,int verbose=0,int mpi_verbose=0):
    r"""
    Removes duplicates in listR modulo permutations which keeps S invariant and fixes 1.
    """
    ## We iterate over possible sets of fixed points.
    cdef MyPermutation Spc,r,Rpc,p
    cdef MyPermutationIterator PONEI
    cdef int ir,irx,nrin,nrout
    cdef list thisset,fixptsets,listRout
    fixptsets = list(subsets(range(1,mu+1)))
    nrin = len(listRin)
    nrout=nrin
    listRout = []
    cdef  Vector_integer_dense Rix
    Rix =  vector(ZZ,nrin)
    #print "In: ",listRin
    Spc=MyPermutation(length=mu)
    Rpc=MyPermutation(length=mu)
    for thisset in fixptsets:
        ## Check if this set of fixed points is valid for preserving S
        if 1 not in thisset or len(thisset)==mu:
            continue
        # if p(j)=j for j=e2+1,...,mu  then p(j+1)=j+1
        cont = False
        for i in range(e2+1,mu,2):
            if i in thisset and i+1 not in thisset:
                cont=True
                break
        if cont:
            continue

        if verbose>1:
            print "Checking fixed point set:",thisset
        PONEI = MyPermutationIterator(mu,fixed_pts=thisset,verbose=mpi_verbose)
        #print "Rix=",Rix
        for ir from 0<=ir<nrin:
            if Rix[ir]<>0:
                continue
            r=listRin[ir]
            for p in PONEI:
                #Spc=Sp._conjugate(p)
                Spc._conjugate_ptr(S._entries,p._entries) 
                #if Sp._conjugate(p)<>Sp: #test = conjugate_perm(Sl,p)
                if not Spc.eq(S):
                    continue
                #rp = r._conjugate(p)
                Rpc._conjugate_ptr(r._entries,p._entries)
                if verbose>2:
                    print p.cycles_non_triv()," preserve S"
                # rp=p*r*p.inverse()
                is_in=0
                for i from 0<= i < nrin:
                    # Need only check those not already removed
                    if Rix[ir]<>0:
                        continue
                    rppp=listRin[i]
                    if rppp.eq(Rpc):
                        is_in=1
                        break
                if is_in==1:
                    if verbose>0:
                        print "G:",print_cycles(r)
                        print "is conjugate mod 1 to: nr.",i
                        print "G':",print_cycles(Rpc)
                        #print print_cycles(Sl)," and ",print_cycles(rp)," are mod-1-conjugate!"
                        print "via p=",print_cycles(p)," p*r*p^-1=",Rpc
                        print "rp is in list!"
                    #list_of_R.remove(rp)
                    nrout=nrout-1
                    Rix[i]=1
        PONEI._c_dealloc()
        if nrout==1:
            break
    listRout=[]
    #print "Out0: ",listRout
    for ir from 0<=ir<nrin:
        if Rix[ir]<>0:
            continue
        r=listRin[ir]
        listRout.append(r)
        #print "Out1: ",listRout
    #print "Out: ",listRout
    return listRout



cpdef list filter_list_mod_1_mod_S_new(list listRin, int mu, int e2,MyPermutation S,int verbose=0):
    r"""
    Removes duplicates in listR modulo permutations which keeps S invariant and fixes 1.
    """
    ## We iterate over possible sets of fixed points.
    cdef MyPermutation Spc,r,Rpc,p,R,R1
    cdef MyPermutationIterator PONEI
    cdef int ir,irx,nrin,nrout
    cdef int i,j,nineq
    cdef list thisset,fixptsets,listRout
    cdef int *checked
    cdef int res
    cdef int* pres=NULL
    nrin = len(listRin)
    listRout = []
    checked = <int*>sage_malloc(sizeof(int)*nrin)
    pres = <int*>sage_malloc(sizeof(int)*mu)
    for i in range(nrin):
        checked[i]=0
    nineq = 0
    for i in range(nrin):
        if checked[i]==1:
            continue
        R = listRin[i]
        for j in range(i+1,nrin):
            if checked[j]==1:
                continue
            R1 = listRin[j]
            #are_mod1_equivalent_c(mu,S,R,S,R1,&res,pres,verbose)
            t,pp = are_conjugate_pairs_of_perms(S,R,S,R1,fix_one=1,verbose=verbose)
            if t==1:
                if verbose>0:
                    print "conjugating map:",pp
                    print "Remove ",R1
                checked[j]=1
                continue
        listRout.append(R)
    if pres<>NULL:
        sage_free(pres)
    if checked<>NULL:
        sage_free(checked)
    return listRout


cpdef list filter_list_of_pairs_mod_1_mod_S_new(list listGin, int mu,int verbose=0):
    r"""
    Removes duplicates in listR modulo permutations which keeps S invariant and fixes 1.
    """
    ## We iterate over possible sets of fixed points.
    cdef MyPermutation Spc,r,Rpc,p,R,R1
    cdef MyPermutationIterator PONEI
    cdef int ir,irx,nrin,nrout
    cdef int i,j,nineq
    cdef list thisset,fixptsets,listRout
    cdef int *checked
    cdef int res
    cdef int* pres=NULL
    nrin = len(listGin)
    listGout = []
    checked = <int*>sage_malloc(sizeof(int)*nrin)
    pres = <int*>sage_malloc(sizeof(int)*mu)
    for i in range(nrin):
        checked[i]=0
    nineq = 0
    for i in range(nrin):
        if checked[i]==1:
            continue
        S = listGin[i][0]
        R = listGin[i][1]
        listGout.append((S,R))
        assert S.order() in [1,2]
        assert R.order() in [1,3]
        for j in range(i+1,nrin):
            if checked[j]==1:
                continue
            S1,R1 = listGin[j]
            #are_mod1_equivalent_c(mu,S,R,S,R1,&res,pres,verbose)
            t,pp = are_conjugate_pairs_of_perms(S,R,S1,R1,fix_one=1,verbose=verbose-1)
            if t==1:
                if verbose>0:
                    print "conjugating map:",pp
                    print "Remove ",listGin[j]
                checked[j]=1
                continue
    if pres<>NULL:
        sage_free(pres)
    if checked<>NULL:
        sage_free(checked)
    return listGout


cpdef are_mod1_equivalent(MyPermutation R1,MyPermutation S1, MyPermutation R2,MyPermutation S2,verbose=0):
    cdef MyPermutation p
    cdef int* pres=NULL
    cdef int res[0]
    cdef int N=R1._N
    assert S1._N==N and R2._N==N and S2._N==N
    pres = <int*>sage_malloc(sizeof(int)*N)
    p=MyPermutation(length=N)
    are_mod1_equivalent_c(R1._N,S1,R1,S2,R2,res,pres,verbose)
    p.set_entries(pres)
    if verbose>0:
        print "res=",res[0]
        print "pres=",print_vec(N,pres)
        print "p=",p.to_cycles()
    if pres<>NULL:
        sage_free(pres)
    return res[0],p


cdef are_mod1_equivalent_c(int N,MyPermutation S1,MyPermutation R1,MyPermutation S2,MyPermutation R2,int* res,int* pres,int verbose=0,check=1):
    r"""
    Check if S1 ~_1 S2 and R1 ~_1 R2 (i.e. equiv. mod 1)


    If check=0 we assume that the cycles are compatible (as in the case if we call this from routine above)
    """
    
    cdef int i,j,ii,cont
    cdef MyPermutationIterator PONEI
    cdef int *epp, *rpp
    cdef int fixS,fixR
    epp = NULL
    rpp=NULL
    for i in range(N):
        pres[i]=i
    #pres = MyPermutation(length=N)._entries
    res[0]=1
    if _are_eq_vec(N,S1._entries,S2._entries) and _are_eq_vec(N,R1._entries,R2._entries):
        return 
        #return True,MyPermutation(range(1,N+1))
    # Then check if they are compatble at all:
    fixS=0; fixR=0
    if check:
        cl1=map(len,S1.to_cycles())
        cl2=map(len,S2.to_cycles())
        cl1.sort(); cl2.sort()
        if cl1<>cl2:
            res[0]=0
        else:
            if verbose>0:
                #print "S1cyc=",S1.to_cycles()
                print "cl1=",cl1
            fixS=cl1.count(1)
            cl1=map(len,R1.to_cycles())
            cl2=map(len,R2.to_cycles())
            cl1.sort(); cl2.sort()
            if cl1<>cl2:
                res[0]=0
            else:
                if verbose>0:
                    #print "R2cyc=",R2.to_cycles()
                    print "cl2=",cl2
                fixR=cl2.count(1)
    else:
        fixS=S1.num_fixed_c()
        fixR=R1.num_fixed_c()
    if res[0]==0:
        return
    res[0]=0
    if verbose>0:
        print "testing G1:={0}:{1}".format(S1.to_cycles(),R1.to_cycles())
        print "testing G2:={0}:{1}".format(S2.to_cycles(),R2.to_cycles())
    epp =  <int *>sage_malloc(sizeof(int)*N)
    if not epp:  raise MemoryError
    rpp =  <int *> sage_malloc(sizeof(int)*N)
    if not rpp:  raise MemoryError
    res[0] = 0
    #pres = MyPermutation(length=N)
    cdef list thisset,fixptsets,listRout
    cdef list fixptsS1=S1.fixed_elements()
    cdef list fixptsR1=R1.fixed_elements()
    cdef list fixptsS2=S2.fixed_elements()
    cdef list fixptsR2=R2.fixed_elements()
    cdef MyPermutation p,p0
    p0=transposition(12,5,7)*transposition(12,6,8)*transposition(12,10,7)*transposition(12,9,8)*transposition(12,10,12)*transposition(12,9,11)*transposition(12,11,12)
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
        if verbose>2:
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
        if verbose>2:
            print "fixset ok.:",thisset
        PONEI = MyPermutationIterator(N,fixed_pts=thisset)
        iteratorlist.append(PONEI)
    sig_on()
    for PONEI in iteratorlist:      
        if verbose>1:
            print "Checking perms in ",PONEI
        for p in PONEI:
            if verbose>2:
                print "p=",p
            if p==p0:
                print "p0 is here!"
            _conjugate_perm(N,epp,S1._entries,p._entries) # epp=p*E*p.inverse()
            if p==p0:
                print "epp=",print_vec(N,epp)
            if not _are_eq_vec(N,epp,S2._entries):
                continue
            _conjugate_perm(N,rpp,R1._entries,p._entries) # rpp=p*R*p.inverse()
            if p==p0:
                print "rpp=",print_vec(N,rpp)
            if _are_eq_vec(N,rpp,R2._entries):
                res[0] = 1
                for i in range(N):
                    pres[i]=p._entries[i]
                break
        if res[0]==1:
            break
        #PONEI._c_dealloc()
    #print "end pres cur=",pres
    #print "pres._entries=",printf("%p ",pres._entries)
    if verbose>0:
        print "Are mod 1 equivalent:",res[0]
        print "pres=",print_vec(N,pres)
    if epp<>NULL:
        sage_free(epp)
    if rpp<>NULL:
        sage_free(rpp)
    #return res
    sig_off()
    
cpdef is_consistent_signature(int ix,int nc=0,int e2=0,int e3=0,int g=0):
    #print "rhs=",12+ix-6*nc-3*e2-4*e3
    return int(12*g == 12+ix-6*nc-3*e2-4*e3)


cpdef are_conjugate_perm(MyPermutation A,MyPermutation B):
   r"""
   Check if the permutation A is conjugate to B
   """
   cdef int N = A._N
   if N<>B._N:
       return 0
   cdef list ctA,ctB
   ctA = A.cycle_type()
   ctB = B.cycle_type()
   if ctA<>ctB:
       return 0
   lA = A.cycles_as_perm()
   lB = B.cycles_as_perm()
   return get_conjugating_perm_list(lA,lB)


cpdef get_conjugating_perm_list(list Al,list Bl):
    r"""
    Find permutation which conjugates Al to Bl, where Al and Bl are lists given by cycle structures.
    """
    cdef dict cperm
    cdef int i
    cperm = {}
    if len(Al)<>len(Bl):
        raise ValueError,"Need lists of the same length! Got:{0} and {1}".format(Al,Bl)
    for i in range(len(Al)):
        cperm[Al[i]-1]=Bl[i]
    return MyPermutation(cperm.values())



cpdef stabiliser_of_R_StoSp(MyPermutation pR,MyPermutation pS,MyPermutation pS1,int verbose=0):
    r"""
    Return a the list of permutations which preserves pR and maps pS to pS1 under conjugation.
    """
    assert pS.order() in [1,2]
    assert pS1.order() in [1,2]
    assert pR.order() in [1,3]
    cdef list Rctype = pR.cycle_type()
    cdef list Rcycles = pR.cycles('list')
    cdef lR = pR.cycles_as_perm()
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
        if p.order()==3:            
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
                print "perm0=",p.to_cycles()
            ## We now act with all combinations of cycles from R
            for ptmp in three_cycles_R:
                do_cont = 0
                pp = p*ptmp
                if verbose>0:
                    print "ptmp=",ptmp.to_cycles()
                    print "pp=",pp.to_cycles()
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
                        print "ppSpp^-1=",pScon.to_cycles()
                        print "pS1=",pS1.to_cycles()
                    continue
                res.append(pp)
            ## We only consider permutations which fixes S
        ## We assume that p ionly has cycles of length 1 or 3 (i.e. it has order 3)
        #perms[i]=MyPermutation(i)
        
    return res

# cdef are_conjugate_wrt_stabiliser(MyPermutation pR,MyPermutation pS,MyPermutation pS1,MyPermutation p_in,int verbose=0):
#     cdef int t
#     cdef MyPermutation pp
#     t = 0; pp=MyPermutation(length=pS.N())
#     are_conjugate_wrt_stabiliser_c(pR,pS,pS1,p_in,pp,&t,verbose)
#     return t,pp

   
#cpdef are_conjugate_wrt_stabiliser(MyPermutation pR,MyPermutation pS,MyPermutation pS1,MyPermutation p_in,MyPermutation p_out,int* t,int verbose=0):
cdef MyPermutation are_conjugate_wrt_stabiliser(MyPermutation pR,MyPermutation pS,MyPermutation pS1,MyPermutation p_in,int* t,int map_one_to=0,int verbose=0):
    r"""
    Return a the list of permutations which preserves pR and maps pS to pS1 under conjugation.
    """
    assert pS.order() in [1,2]
    assert pS1.order() in [1,2]
    assert pR.order() in [1,3]
    cdef list Rctype,Rcycles,lR,ll,l0,l
    cdef list cycles1,cycles3
    cdef int nf,nf1,d1,d3,i,j,mu,do_cont,dd1,dd3
    cdef list fixedS,fixedS1,nonfixedS,nonfixedS1,res,one_cycles,three_cycles
    cdef MyPermutationIterator perms1,perms3
    cdef MyPermutation p,pp,ptmp,pScon
    lR = pR.cycles_as_perm()
    Rctype = pR.cycle_type()
    Rcycles = pR.cycles('list')
    d1 = Rctype.count(1)
    d3 = Rctype.count(3)
    t[0] = 0
    mu = <int>pR.N()
    p_out = MyPermutation(length=mu)
    cycles1 = []; cycles3=[]
    for c in Rcycles:
        if len(c)==1:
            cycles1.append(c)
        elif len(c)==3:
            cycles3.append(c)
    if verbose>0:
        print "cycles1=",cycles1
        print "cycles3=",cycles3
    nf = pS.num_fixed(); nf1 = pS1.num_fixed()
    if nf<>nf1:
        pp = MyPermutation(length=pS.N())
        return pp
    dd1 = int(max(d1,1)); dd3=int(max(d3,1))
    perms1 = MyPermutationIterator(dd1)
    perms3 = MyPermutationIterator(dd3)
    ## Then we also have to include permutations by the 3-cycles of R
    three_cycles_perm =[]
    for p in pR.cycles(type='perm'):
        if p.order()==3:            
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
    if verbose>0:
        print "three_cycles_R=",three_cycles_R
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
            if cycles3<>[]:
                three_cycles = p3(cycles3) #permute_list(cycles[3],p3)
            else:
                three_cycles=[]
            l.extend(three_cycles)
            perm = MyPermutation(l)
            ll = flatten_list2d(l)
            p = get_conjugating_perm_list(ll,lR)
            if verbose>0:
                print "ll=",ll
                print "perm3=",p3
                print "perm0=",p.to_cycles()
            ## We now act with all combinations of cycles from R
            for ptmp in three_cycles_R:
                do_cont = 0
                pp = p._mult_perm(ptmp)
                if verbose>0:
                    print "ptmp=",ptmp.to_cycles()
                    print "pp=",pp.to_cycles()
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

                # for x in fixedS:
                #     if pp(x) not in fixedS1:
                #         do_cont=1
                #         if verbose>0:
                #             print "Fails at preserving fixed points!"                        
                #         break                
                # if do_cont==1:
                #     continue
                # for x in nonfixedS:
                #     if pp(x) not in nonfixedS1:
                #         do_cont=1
                #         if verbose>0:
                #             print "Fails at preserving non-fixed points!"
                #         break                
                if do_cont==1:
                    continue
                # If we are here we may have a conjugating map
                if map_one_to>0 and pp(1)<>map_one_to:
                    continue
                pScon = pS.conjugate(pp)
                if verbose>0:
                    print "pS^pp=",pScon.to_cycles()
                if pScon==pS1:
                    t[0] = 1
                    #p_out = pp
                    if verbose>0:
                        print "p_out = ",pp.to_cycles()
                    return pp
            ## We only consider permutations which fixes S
        ## We assume that p ionly has cycles of length 1 or 3 (i.e. it has order 3)
        #perms[i]=MyPermutation(i)       
#    return 0,0


                    
   
cpdef flatten_list2d(list l):
    r"""
    Takes a list of list and flattens it.
    """
    cdef int i,n = len(l)
    cdef list res = copy(l[0])
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
