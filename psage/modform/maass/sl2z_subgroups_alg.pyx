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

cpdef list_all_admissable_pairs(list sig,int get_details=1,int verbose=0,int get_one_rep=0,int get_cong=0):
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
    - get_cong -- set to one if you want one congruence subgroup
    
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
    ### Note that it is not a restriction to fix S, once the cycle structure is fixed.
    ### It is simply equivalent to renaming the digits. 

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
    ## Without loss of generality we may take those as
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
    for pR in PRI: #ii from 0<= ii <max_num:
        #_to_cycles2(mu,pR._entries,Rcycles,rcycle_lens,num_rcycles)
        # we might also make some apriori requirements on R
        # 1) R can not contain all fixed points of E in one cycle
        if verbose>0:
            print "S=",Sp #.cycles() #print_vec(mu,Sptr)
            print "R=",pR #s.cycles() #print_vec(mu,<int *>PRI._current_perm)
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
        #Tp= Sp*pR  #MyPermutation(length=mu)
        #_mult_perm(mu,Sp._entries,pR._entries,Tp._entries)
        #Tcyc=perm_to_cycle_c(mu,<int*>Tp._entries)
        Tcyc=perm_to_cycle_c(mu,Tptr)
        if verbose>0:
            print "Tp=",Tcyc
            print "current fixedpts=",PRI.fixed_pts()
            print "number of cusps=",len(Tcyc)
        if len(Tcyc)<>h:
            if verbose>0:
                print "remove!"
                print "next!"
            continue
        #if used<>NULL:
        #    sage_free(used)
        #    used=NULL
        if verbose>0:
            #print "used=",print_vec(mu,used)
            print "current fixedpts1=",PRI.fixed_pts()
            print "rfx=",print_vec(e3,rfx)

        if get_one_rep:
            if get_cong:
                G=MySubgroup(o2=Sp,o3=pR)
                if G.is_congruence():
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
        try:
            for rc in Rcycles_list:
                #for i from 0<=i<=num_rcycles:
                #if rcycle_lens[i]<>3:
                if len(rc)<>3:
                    #rcycle_beg+=rcycle_lens[i]
                    continue
                # only the 3-cycles are relevant here (not the 1-cycles)
                # cy = (a b c)
                if verbose>0:
                    print "rc=",rc
                a=rc[0]
                b=rc[1]
                c=rc[2]
                #a=Rcycles[rcycle_beg] #rc[0]
                #b=Rcycles[rcycle_beg+1] #rc[1]
                #c=Rcycles[rcycle_beg+2] #rc[2]
                #rcycle_beg+=rcycle_lens[i]
                if verbose>0:
                    print "a,b,c=",a,b,c
                used[a-1]=1 # We have fixed a
                if verbose>0:
                    print "used=",print_vec(mu,used)
                ## If b and c are equivalent with respect to S, i.e. if
                ##   i) S(b)=b and S(c)=c, or
                ##  ii) (b b') and (c c') are two cycles of S with c and c' not used previously and not fixed points of R
                ## iii) 
                ## then we choose b < c
                if (b<=e2 and c<=e2) or ((b>=end_fc and used[Sptr[b-1]-1]==0) and (c>=end_fc and used[Sptr[c-1]-1]==0)):
                    if verbose>0:
                        print b," (b c) and ",c," are equivalent in",rc
                    if b>c:
                        if verbose>0:
                            print "remove (b>c)!"
                        raise StopIteration()
                for j from a+1<= j < b:  # I want to see if there is a j, equivalent to b, smaller than b
                    t1 = (b<=e2 and j<=e2)
                    t2 = (b>=end_fc and used[Sptr[b-1]-1]==0 and used[b-1]==0)
                    t3 = (j>=end_fc and used[Sptr[j-1]-1]==0 and used[j-1]==0)
                    #if (b<=e2 and j<=e2) or ((b>end_fc and used[Sptr[b-1]-1]==0 and used[b-1]==0) and (j>end_fc and used[Sptr[j-1]-1]==0 and used[j-1]==0)):
                    if t1 or (t2 and t3):
                        if verbose>0:
                            print j," (b) and ",b," are equivalent in",rc
                        if j<b:
                            if verbose>0:
                                print "remove!"
                            raise StopIteration()
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
                            raise StopIteration()

                used[c-1]=1
                if verbose>0:
                    print "used=",print_vec(mu,used)
        except StopIteration:
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
        if get_cong:
            # We didn't find any congruence subgroups...
            print "We didn't find any congruence subgroups..."
            return
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
    if len(list_of_R)>1:
        list_of_R_tmp=filter_list_mod_1_mod_S(list_of_R,mu,e2,Sp,verbose)
        #print "returned:",list_of_R_tmp

        list_of_R=list_of_R_tmp
    if verbose>=0:    
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
    lc_psl,lc_psl_maps,list_of_R=filter_list_mod_psl_mod_S(list_of_R,mu,e2,e3,Sp,verbose)
    if verbose>=0:    
        print "List of conjugacy classes mod PSL(2,Z):" 
        for i from 0 <= i < len(list_of_R):
            print perm_to_cycle(list_of_R[i])    
    ##  Next part is to check outer automorphisms, i.e. we reduce the list of conjugacy classses
    ##  Modulo *:  E->E, R->ER^2E
    lc_pgl,lc_pgl_maps,list_of_R=filter_list_mod_pgl_mod_S(list_of_R,mu,e2,e3,Sp,verbose)
    d = dict()
    d['sig']=sig
    d['numg']=len(list_of_groups)
    d['groups']=dict() #list_of_groups
    d['groups_list']=list_of_groups
    d['num_psl_cc']=len(lc_psl.keys())
    d['len_psl_cc']=map(len, lc_psl.values())
    d['psl_conj']=lc_psl
    d['psl_conj_maps']=lc_psl_maps
    d['pgl_conj']=lc_pgl
    d['pgl_conj_maps']=lc_pgl_maps
    d['num_pgl_cc']=len(lc_pgl.keys())
    d['len_pgl_cc']=map(len, lc_pgl.values())
    ## The most reduced ist of groups.
    for rr in lc_pgl.keys():
        d['groups'][Sp,rr]=dict()
        d['groups'][Sp,rr]['psl_conj']=lc_psl[rr]
        d['groups'][Sp,rr]['num_psl_cc']=len(lc_psl[rr])
        d['groups'][Sp,rr]['pgl_conj']=lc_pgl[rr]
        d['groups'][Sp,rr]['num_pgl_cc']=len(lc_psl[rr])
    #res['T']=list_of_T
    #d['Gstar']=Gstar
    #d['Gstarmap']=Gstarmap
    #res=dict()
    #res[tuple(sig)]=d
    #print "res=",res
    #print "-----END---"
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


cpdef tuple filter_list_mod_pgl_mod_S(list listRin, int mu, int e2,int e3,MyPermutation S,int verbose=0,int mpi_verbose=0):
    r""" Removes duplicates in listR modulo squares and permutations which keeps S invariant.

    The purpose of this step is to reduce the list we later have to check completely.
    """
    if verbose>0:
        print "Filter mod PGL!"
    cdef MyPermutation r,r0,Rpc,p,Spc
    cdef MyPermutationIterator PONEI
    cdef int i,ir,irx,nrin,nrout
    cdef list thisset,fixptsets,listRout

    nrin = len(listRin)
    nrout=nrin
    listRout = []
    cdef  Vector_integer_dense Rix    
    Rix =  vector(ZZ,nrin)
    Rpc=MyPermutation(length=mu)
    Spc=MyPermutation(length=mu)
    Rpc=MyPermutation(length=mu)
    conjugates=dict()
    conjugate_maps=dict()
    cdef list rfx_list=range(e2+1,e2+2*e3,2)
    for ir from 0<=ir<nrin:
        conjugates[listRin[ir]]=[listRin[ir]]
        #conjugate_maps[listRin[ir]]=list()
    cdef list iterator_list=[]
    # Remove unsuitable sets of fixed points
    fixptsets = list(subsets(range(1,mu+1)))
    #print "fixptsets0=",fixptsets
    for thisset in fixptsets:
        ## Check if this set of fixed points is valid for preserving S
        if len(thisset)>=mu-1:  # trivial
            continue
        # if p(j)=j for j=e2+1,...,mu  then p(j+1)=j+1
        cont = False
        for i in range(e2+1,mu,2):
            if i in thisset and i+1 not in thisset:
                cont=True
                break
        if len(rfx_list)==1 and rfx_list[0] not in thisset:
            cont=True
        if cont:
            continue # fixptsets.remove(thisset)
        PONEI = MyPermutationIterator(mu,fixed_pts=thisset,verbose=mpi_verbose)
        iterator_list.append(PONEI)
        
    #print "iterator_list="
    #for PONEI in iterator_list:
    #    print PONEI
    for ir from 0<=ir<nrin:
        if Rix[ir]<>0:
            continue
        r0=listRin[ir]*listRin[ir]
        if verbose>0:
            print "r0=",r0.to_cycles()
        for irx from 0<=irx<nrin:
            if Rix[irx]<>0:
                continue
            r=listRin[irx]
            if nrout<=1:
                break
            for PONEI in iterator_list:
                if verbose>1:
                    print "Checking fixed point set:",PONEI.fixed_pts()
                for p in PONEI:
                    cont=False
                    for i from 0<=i<e3:
                        if p._entries[rfx_list[i]-1] not in rfx_list:
                            cont=True
                            break
                    if cont:
                        continue                    
                    Rpc._conjugate_ptr(r0._entries,p._entries)
                    Spc._conjugate_ptr(S._entries,p._entries)
                    if not Spc.eq(S):
                        continue
                    if Rpc.eq(r):
                        if verbose>0:
                            print "r0=",listRin[ir].to_cycles()
                            print "r0^2=",r0.to_cycles()
                            print "r=",r.to_cycles()
                        # and add to list of conjugates
                        Rix[irx]=1
                        if not r.eq(listRin[ir]):
                            if verbose>0:
                                print "Remove this conjugate from list of R's and insert into conjugates!"
                            conjugates[listRin[ir]].extend(conjugates.pop(r))
                        #conjugates.pop(r)
                        conjugate_maps[(listRin[ir],r)]=S*p
                        #if conjugates.has_key(r):
                        #    conjugates.pop(r)
                        nrout-=1
                        break
                if Rix[irx]<>0:
                    break
                # Otherwise check how mu
                #print "R^(1j)=",Rpc

        if nrout<=1:
            break
    for ir from 0<=ir<nrin:
        r=listRin[ir]
        if conjugates.has_key(r):
            listRout.append(r)
    if verbose>1:
        print "conjugates1:",conjugates
        print "conjugate_maps:",conjugate_maps
        print "conjugates.keys:",conjugates.keys()
    return conjugates,conjugate_maps,listRout    




cpdef tuple filter_list_mod_psl_mod_S(list listRin, int mu, int e2,int e3,MyPermutation S,int verbose=0,int mpi_verbose=0):
    r""" Removes duplicates in listR modulo permutations of the form (1j) which keeps S invariant.

    The purpose of this step is to reduce the list we later have to check completely.
    """
    cdef MyPermutation r,r0,Rpc,p,Spc
    cdef MyPermutationIterator PONEI
    cdef int i,ir,irx,nrin,nrout
    cdef list thisset,fixptsets,listRout
    fixptsets = list(subsets(range(1,mu+1)))
    nrin = len(listRin)
    nrout=nrin
    listRout = []
    cdef  Vector_integer_dense Rix    
    Rix =  vector(ZZ,nrin)
    Rpc=MyPermutation(length=mu)
    Spc=MyPermutation(length=mu)
    Rpc=MyPermutation(length=mu)
    conjugates=dict()
    conjugate_maps=dict()
    cdef list rfx_list=range(e2+1,e2+2*e3,2)
    for ir from 0<=ir<nrin:
        conjugates[listRin[ir]]=[listRin[ir]]
        #conjugate_maps[listRin[ir]]=list()
    for thisset in fixptsets:
        ## Check if this set of fixed points is valid for preserving S
        if len(thisset)==mu:
            continue
        # if p(j)=j for j=e2+1,...,mu  then p(j+1)=j+1
        cont = False
        for i in range(e2+1,mu,2):
            if i in thisset and i+1 not in thisset:
                cont=True
                break
        # We also can't move the fixed points of R (since all R's have the same fixed-point set)
        #if len(rfx_list)==1 and rfx_list[0] not in thisset:
        #    cont=True
        if cont:
            continue
        if verbose>1:
            print "Checking fixed point set:",thisset
        PONEI = MyPermutationIterator(mu,fixed_pts=thisset,verbose=mpi_verbose)
        for p in PONEI:
            cont=False
            for i from 0<=i<e3:
                if p._entries[rfx_list[i]-1] not in rfx_list:
                    cont=True
            #if cont:
            #    continue
            for ir from 0<=ir<nrin:
                if Rix[ir]<>0:
                    continue
                r0=listRin[ir]
                #print "r=",r
                Rpc._conjugate_ptr(r0._entries,p._entries)
                Spc._conjugate_ptr(S._entries,p._entries)
                if not Spc.eq(S):
                    continue
                #print "R^(1j)=",Rpc
                for irx from 0<=irx<nrin:
                    if Rix[irx]<>0:
                        continue
                    r=listRin[irx]
                    if Rpc.eq(r):
                        if r0.eq(r):
                           continue 
                        if verbose>0:
                            print "Remove this conjugate from list of R's and insert into conjugates!"
                            print "r0=",r0
                            print "r=",r
                            print "p=",p
                        # and add to list of conjugates
                        Rix[irx]=1
                        if not r0.eq(r):
                            conjugates[r0].extend(conjugates.pop(r))
                        #conjugates.pop(r)
                        conjugate_maps[(r0,r)]=p
                        #if conjugates.has_key(r):
                        #    conjugates.pop(r)
                        nrout-=1
                        break
                # Otherwise check how mu
            if nrout<=1:
                break
        if nrout<=1:
            break
    if verbose>0:
        print "conjugates0:",conjugates
    for ir from 0<=ir<nrin:
        r=listRin[ir]
        #if Rix[ir]<>0:
        if conjugates.has_key(r):
            #conjugates.pop(r)
            #continue
            listRout.append(r)
    if verbose>0:
        print "conjugates1:",conjugates
        print "conjugate_maps:",conjugate_maps
        print "conjugates.keys:",conjugates.keys()
    return conjugates,conjugate_maps,listRout



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
