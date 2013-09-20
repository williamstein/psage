# cython: profile=True
# -*- coding: utf-8 -*-
r"""
Algorithms and classes for permutations representing subgroups of the modular group, as implemented in 'MySubgroup'.

CLASSES:


AUTHOR:

 - Fredrik Stroemberg



"""

#*****************************************************************************
#  Copyright (C) 2010 Fredrik StrÃ¶mberg <stroemberg@mathematik.tu-darmstadt.de>,
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


### higher levels of debugging for development and which should not be controlled by user
DEF debug = 0


#cdef extern from "convert.h":
#    cdef void t_INT_to_ZZ( mpz_t value, long *g )
import cProfile; import profile


from sage.modules.vector_integer_dense cimport Vector_integer_dense
from sage.rings.integer cimport Integer
import sage.combinat.permutation_cython
from sage.combinat.permutation_cython cimport reset_swap,next_swap
from sage.combinat.permutation import Permutation_class
from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
from sage.all import deepcopy,copy,ZZ,vector,subsets
from sage.all import binomial #,lcm
import cython
from psage.modform.maass.mysubgroup import MySubgroup
from psage.modform.maass.mysubgroups_alg cimport SL2Z_elt

from libc.stdlib cimport malloc, free

#from Cython.Utils import cached_function, cached_method
#from sage.all import cached_method
cdef extern from "math.h":
    double fabs(double)
    double fmax(double,double)
    int ceil(double)
    double M_LN10
    double log(double)
    


cdef extern from "mpz_pylong.h":
    cdef mpz_get_pylong(mpz_t src)
    cdef mpz_get_pyintlong(mpz_t src)
    cdef int mpz_set_pylong(mpz_t dst, src) except -1

    
cdef class MyPermutation(SageObject):
    def __cinit__(self,entries=[],int length=0,int init=1,int verbose=0,int check=1,int rep=3):
        r"""
        If init=0 then we do not initialize the list.
        TODO?: faster check of consistency of input?
        """
        #cdef int i1
        #cdef long i2
        #cdef Integer i3
        #i1=2**200
        #i2=2**200
        #i3=Integer(2)**Integer(200)
        cdef int* used=NULL
        cdef int i,ok,ei
        cdef str s
        self._entries = NULL
        self._str=''
        self._rep=rep
        self._list = []
        self._cycles_ordered_as_list=[]
        self._order=-1
        self._cycles_list_is_ordered = 0
        self._verbose=verbose
        ## If we have a list of lists we are trying to initialize from cycle structure.
        #if self._cycles<>NULL:
        #    sage_free(self._cycles)
        self._cycles = NULL
        self._cycles_list = []
        if verbose > 0:
            print "THERE"
        self._cycles_permutations = []
        if verbose>0:
            print "HERE0!"
        self._cycle_type=[]
        if verbose>0:
            print "HERE1!"
        #if self._cycle_lens<>NULL:
        #    sage_free(self._cycle_lens)
        self._cycle_lens = NULL
        if verbose>0:
            print "HERE2!"

        self._num_cycles = 0
        cdef list entries_list=[]
        if verbose>0:
            print "entries={0} of type {1}".format(entries,type(entries))
        if entries==[]:
            if length>0:
                self._N = length
            else:
                raise ValueError,"Invalid Input for permutation! entries: {0} and length:{1}".format(entries,length)           
        else:
            if isinstance(entries,basestring):
                # if the string represents a list
                s = <str> entries
#               if entries.count('[')==1 and entries.count(']')==1:
                # The string can either be given as a string of cycles or a string giving a list
                if s.count('[')==1 and s.count(']')==1:
                    entries_list=eval(s)
                elif s.count('(')>=1 and s.count(')')>=1:
                    ## Assume it is a string representation of cycles
                    entries_list=self._cycles_from_str(entries)
                else:
                    # It could be given by self.export_as_string()
                    entries_list=self.list_from_string(s)
            elif hasattr(entries,list):
                entries_list=entries.list()
            elif isinstance(entries,(list,tuple)):
                if isinstance(entries[0],(int,Integer)):
                    entries_list = <list>entries
                elif isinstance(entries[0],(list,tuple)):
                    try:
                        entries_list = self._list_from_cycles(entries)
                    except:
                        pass
            if len(entries_list)==0 or (length>0 and length<>len(entries_list)):
                raise ValueError,"Invalid Input for permutation!! entries: {0}".format(entries)           
            if length>0:
                self._N=length
            else:
                self._N = len(entries_list)
        self._init = 1        
        if self._N>0:
            if check==1:
                #print "N=",self._N
                used=<int*>malloc(sizeof(int)*self._N)
                if used==NULL:
                    raise MemoryError
                for i in range(self._N):
                    used[i]=0
                ok=1
            self._entries = <int*> sage_malloc(sizeof(int)*self._N)
            if self._entries==NULL:
                raise MemoryError
            if len(entries_list)>0:
                # Check consistency: we need every element 1,2,...,N exactly once.
                if check==1 and used<>NULL:
                    for i in range(self._N):
                        ei = <int>entries_list[i]
                        self._entries[i]=ei
                        if ei>self._N or ei<1:
                            free(self._entries); self._entries=NULL                            
                            free(used)
                            raise ValueError,"Invalid Input for permutation!!!!! entries: {0}".format(entries)
                        if used[ei-1]>0:
                            free(self._entries); self._entries=NULL
                            free(used)
                            raise ValueError,"Invalid Input for permutation!!!!!! entries:{0}".format(entries)
                        else:
                            used[ei-1]=1
                else:
                    for i in range(self._N):
                        self._entries[i]=<int>entries_list[i]
                self._list = []
                for i in range(self._N):
                    self._list.append(<int>entries_list[i])
            elif init==1:
                for i in range(self._N):
                    self._entries[i]=i+1
                self._list = range(1,self._N+1)
            else:
                self._init = 0
        self._hash = 0  
        if check==1 and used<>NULL and self._N>0:
            free(used)

    def __init__(self,entries=[],int length=0,int init=0,int verbose=0,int check=1,int rep=0):
        r"""
        The possible representations are the following: if self=MyPermutation(length=7)
        rep = 0 => 1234567 (compact list)
        rep = 1 => [1, 2, 3, 4, 5, 6, 7] (verbose list)
        rep = 2 => [[1], [2], [3], [4], [5], [6], [7]] (cycles)
        """
        pass


            
    def _cycles_from_str(self,s):
        r"""
        Convert a string of the form (1 3 5)(2 7 9)(4 8 11)(6 10 12)
        to corresponding list representing this permutation.
        """
        s=s.replace(")(","),(")
        s=s.replace(" ",",")
        s=s.replace("(","[")
        s=s.replace(")","]")
        s="["+s+"]"
        cycles=eval(s) # list of lists
        self._cycles_list = cycles
        return self._list_from_cycles(cycles)

    cpdef list _list_from_cycles(self,list cycles):
        cdef int i,mu,mutmp
        cdef list cycle,res
        mu=1
        # first search for the potential N
        for cycle in cycles:
            mutmp=max(cycle)
            if mutmp>mu:
                mu=mutmp
        res=range(1,mu+1)
        if self._verbose>0:
            print "res=",res
        for cycle in cycles:
            for i in range(len(cycle)-1):
                res[cycle[i]-1]=cycle[i+1]
            if self._verbose>0:
                print "cycle[-1]-1=",cycle[-1]-1
            res[cycle[-1]-1]=cycle[0]                    
        if self._verbose>0:
            print "res=",res
        return res
        
        
    def __copy__(self):
        r"""
        Copy self.
        """
        #l = self.list()
        #res = MyPermutation(l)
        return self._copy_c()

    cpdef _copy_c(self):
        r" copy self"
        cdef MyPermutation res
        res = MyPermutation(length=self._N)
        res.set_entries(self._entries)
        return res

    #def __deepcopy__(self):
    #    res = MyPermutation(self._list)
    #    return res

    cpdef str export_as_string(self,str sep='0'):
        r"""
        Export self as string.
        """
        cdef str s=''
        cdef int base
        base=self._N+1
        cdef int i
        # We don't want to use alphanumeric characters as separators except for 0
        if sep.isalnum() and sep<>'0':
            raise ValueError,"Need a non-alphanumeric separator! Got {0}".format(sep)
        if base<=36 and (sep=='0' or sep==""):
            for i from 0<=i<self._N:
                s+=Integer(self._entries[i]).str(base)
        else: # If we have more than one symbol per char we insert zeros instead
            for i from 0<=i<self._N-1:
                s+=str(self._entries[i])+sep
            s+=str(self._entries[self._N-1])
        if sep<>"0" and sep<>"":
            s+=".{0}".format(sep)  ## If not "0" we give the separator explicitly
        return s

    cpdef list list_from_string(self,s):
        r"""
        Try to make a list out of a string  returned by the above function
        """
        base = 10
        if s.count(".")==1:
            s,sep = s.split(".")
        elif s.count("0")>0:
            sep = "0"
        else:
            sep=""  ##
            base = len(s)+1 ### The base
        if sep<>"":
            slist = s.split(sep)
        else:
            slist = s
        res = []
        for x in slist:
            res.append(Integer(x,base))
        return res
    
    cdef void c_new(self,int* entries):
        cdef int i
        if self._entries<>NULL:
            for i from 0 <= i< self._N:
                self._entries[i]=entries[i]          


    def __reduce__(self):
        return(MyPermutation,(self.list(),self._N,self._init,self._verbose))

    def _pickle(self):
        if self._verbose>0:
            print "in pickle!"
        return self._pickle_version0(), 0
    
    def _unpickle(self, data, int version=0):
        if self._verbose>0:
            print "in _unpickle!"
        if version == 0:
            self._unpickle_version0(data)
        else:
            raise RuntimeError, "unknown permutation version (=%s)"%version
        
    cdef _pickle_version0(self):
        if self._verbose>0:
            print "in _pickle v0!"
        s = str(self._N)
        s+=" "+self.export_as_string(sep='x')
        return s

    cdef _unpickle_version0(self, data):
        r"""
        Data is a pickled string
        """
        if self._verbose>0:
            print "in _unpickle_v0!"
        print "data=",data
        if not isinstance(data,basestring):
            raise ValueError,"Can not unpickle object: {0} as MyPermutation!".format(data)
        smu,sentries=data.split()
        self._N = int(smu)
        cdef int i
        svec=sentries.split("0")
        cdef list entries=[]
        for i from 0<=i<self._N:
            if svec[i]<>'':
                entries.append(int(svec[i]))

        if self._entries<>NULL:
            sage_free(self._entries)
        self._entries = <int*> sage_malloc(sizeof(int)*self._N)
        #self.set_entries(entries)
        for i from 0<=i<self._N:
            self._entries[i]=entries[i]
        self._list = entries

        


        
    def __richcmp__(self, right, int op):
        #if op <>2 and op <>3:
        #    raise NotImplementedError,"Ordering of permutations is not implemented!"
        if type(self)<>type(right) or (<MyPermutation>right).N() <> self.N():
            if op==2:
                return False
            elif op==3:
                return True
            else:
                raise ValueError,"Can not compare {0} with {1}!".format(self,right)
        if op in [2,3]:        
            t = _are_eq_vec(self.N(),(<MyPermutation>self)._entries,(<MyPermutation>right)._entries)
            if op==2 and t==1:
                return True
            if op==3 and t==0: # != # op=2 is ==
                return True
            return False
        # If we are here we need to compare
        #print "A = ",self
        #print "B=",right
        lA = self.cycles_ordered_as_list()
        lB = right.cycles_ordered_as_list()
        if op == 0:
            return lA < lB
        if op == 1:
            return lA <= lB
        if op == 4:
            return lA > lB
        if op == 5:
            return lA >= lB
        raise ValueError,"Can not compare permutation with op={1}!".format(op)
            

    cdef int set_entries(self, int *entries):
        r"""
        Unsafe setting of entries via pointer.
        """
        cdef int i
        if entries==NULL:
            raise MemoryError
        for i in range(self._N):
            self._entries[i] = entries[i]
            if self._verbose>0:
                print "setting {0} to {1}".format(i,entries[i])
        # and reset the list and string of self too
        self._list = []
        self._str=''
        self._init=1
        self._order = -1
        self._cycles_ordered_as_list=[]
        if self._cycles<>NULL:
            sage_free(self._cycles)
        self._cycles=NULL
        self._cycles_list=[]
        self._cycles_permutations = []
        self._cycle_type = []
        if self._cycle_lens<>NULL:
            sage_free(self._cycle_lens)
        self._cycle_lens = NULL
        self._num_cycles = 0
        return 0

    def parent(self):
        return "MyPermutation on {0} letters".format(self._N)
        
    def __call__(self,i):
        if isinstance(i,(int,Integer)) and 1 <= i and i<= self._N:
            return self._entries[i-1]
        elif isinstance(i,(tuple,list)) and len(i)==self._N:
            res = [0 for j in range(len(i))]
            for j in range(len(i)):        
                res[j]=i[self._entries[j]-1]
            if isinstance(i,tuple):
                return tuple(res)
            return res
        else:
            raise NotImplementedError,"Can not apply permutation to {0}".format(i)

    def __getitem__(self, i):
        """
        EXAMPLES::
        

        """
        if isinstance(i,(int,Integer)) and i>=0 and i<self._N and self._entries<>NULL:
            return self._entries[i] 
        else:
            raise NotImplementedError,"Indexing outside bound:  i={0}".format(i)

    cpdef int _get_unsafe(self,int i):
        cdef int j
        j = self.get_unsafe(i)
        return j

    cdef int get_unsafe(self,int i):
        r"""
        Unsafe get. No checking performed!
        """
        return self._entries[i]


    def __len__(self):
        return self._N

    def __list__(self):
        #print "in list!"
        return self.list()
    
    def __hash__(self):
        """
        For hashing.
        EXAMPLES::
        
            sage: d = {}
            sage: p = MyPermutation([1,2,3])
            sage: d[p] = 1
            sage: d[p]
            1

        """
        if self._hash == 0:
            self._hash = str(self).__hash__()
        return self._hash

    def __str__(self):
        if self._str=='':
            if self._rep==0:
                res = self.export_as_string(sep='') #str(self.list())
            elif self._rep==1:
                res = str(self.list())
            elif self._rep==2:
                res = str(self.cycles())
            elif self._rep==3:
                res = str(self.cycles())
                res=res.replace("], [",")(")
                res=res.replace("[[","(")
                res=res.replace("]]",")")
                res=res.replace(",","")
            else:
                raise NotImplementedError
            self._str = res
        return self._str
            
        #return repr(self)

    def  __dealloc__(self):
        self._dealloc_c()

    def _dealloc_c(self):
        if self._entries <> NULL:
            sage_free(self._entries)
            self._entries = NULL
        if self._cycles <>NULL:
            sage_free(self._cycles)
            self._cycles = NULL
        if self._cycle_lens<>NULL:
            sage_free(self._cycle_lens)
            self._cycle_lens = NULL

    
            
    def set_rep(self,rep=0):
        if rep<>self._rep:
            self._str = ''
        self._rep=rep

    def get_rep(self):
        return self._rep
        
    def __repr__(self):
        r"""
        Represent self.
        """
        return str(self)

    def __latex__(self):
        r"""
        Latex representation of self.
        """
        res = str(self.cycles())
        res=res.replace("], [",")(")
        res=res.replace("[[","(")
        res=res.replace("]]",")")
        res=res.replace(",","\,")
        res = "$ {0}$".format(res)
        return res
    
    def set_rep(self,rep):
        self._rep=rep
        self._str=''


    def order(self):
        return self.__order__()
        
    def __order__(self):
        if self._order==-1:
            c = self.cycles()
            self._order = 1
            for l in self.cycle_type():
                self._order = lcm(self._order,l)  #lcm( map(len,c))
        return self._order

    cpdef is_order_eq(self,int o):
        r"""
        Testing if order of self is o without computing the order if not necessary
        """
        return self.is_order_eq_c(o)
    
    cdef int is_order_eq_c(self,int o):
        r"""
        Testing if order of self is o without computing the order if not necessary
        """
        if self._order<>-1:
            if o % self._order == 0 :
                return 1
            else:
                return 0
        cdef int i,j,k,l,ok
        for l in range(1,o+1):
            if o % l <>0:
                continue
            ok = 1
            for i in range(1,self._N+1):
                k = i
                for j in range(l):
                    k = self._entries[k-1]
                if k<>i:
                    ok = 0
                    break
            if ok == 1:
                return 1
        return 0
    
    def is_identity(self):
        r"""
        Test if self is the identity
        """
        return self.is_order_eq(1)==1
    
    cpdef _get_order(self):
        return self._order

    def __pow__(self,k,m):
        ## TODO: More efficient...
        if m<>None:
            raise RuntimeError, "__pow__ dummy argument ignored"
        if not isinstance(k,(int,Integer)):
            raise TypeError,"Can only take integer powers! not:{0}".format(k)
        if k>=1:
            o = self._get_order() # it is not efficient to calculate the order
            if o>0 and k>o: # unless it is set at the beginning
                k = k % o
            if k%2==0:
                tmp=self.square()
                res=tmp
                #for i in range(k/2-1):
                for i in xrange(k/2-1):
                    res=res*tmp
            else:
                tmp=self.square()
                res=self
                for i in xrange((k-1)/2):
                    res=res*tmp
            return res
        elif k==-1:
            return self.inverse()
        elif k<0:
            k = -k
            o = self._get_order()
            if o>0 and k>o:
                k = k % o
            return self.inverse().__pow__(k)
        elif k==0:
            return MyPermutation(length=self.N(),init=1)

    cpdef pow(self,int k):
        cdef MyPermutation res,fac
        #cdef int* entries_res 
        #res = MyPermutation(length=self._N)
        cdef int p,m,j,o
        if k==0:
            return MyPermutation(length=self._N)
        if k==1:
            return self
        if k==-1:
            return self.inverse()
        if k<0:
            k = -k
            fac = self.inverse()
        else:
            fac = self
        o = self.order()
        k = k % o
        res = MyPermutation(length=self._N)

        for p,m in ZZ(k).factor():
            if m==1:
                if p==2:
                    res=res*fac.square()
                else:
                    for j from 0<=j<p:
                        res=res*fac
            else:
                for j from 0<=j<m:
                    res=res*fac.pow(p)
        return res

        # self._pow_c(entries_res,k)
        #res.set_entries(entries_res)
        #return self._pow


    
    # cdef _pow(self,int k, int* res):
    #     cdef int i
    #     if k==0:
    #         for i from 0<=i<=self._N:
    #             res[i]=i+1
    #         return
    #     if k==1:
    #         for i from 0<=i<=self._N:
    #             res[i]=self._entries[i]
    #         return
    #     cdef int inv,o
    #     o = self._order
    #     if k >1:
    #         inv = 0 
    #     if k>o:
    #         k = k % o
    #         inv = 1
    #     if k%2==0:
    #         tmp=self.square()
    #         res=tmp
    #         #for i in range(k/2-1):
    #         for i from 0<=i<k/2-1:
    #             res=res*tmp
    #         else:
    #             tmp=self.square()
    #             res=self
    #             for i in xrange((k-1)/2):
    #                 res=res*tmp
    #         return res
    #     elif k<0:
    #         k = -k
    #         o = self.order()
    #         if k>o:
    #             k = k % o
    #         return self.inverse().__pow__(k)
    #     elif k==0:
    #         return MyPermutation(length=self.N(),init=1)  
        
    
    def inverse(self):
        return self._inverse()
        # l=range(1,self.N()+1)
        # for i in range(self.N()):
        #     j=self._entries[i]
        #     l[j-1]=i+1
        # return MyPermutation(l)
        
    cpdef _inverse(self):
        cdef int* res_ent=NULL
        cdef MyPermutation res
        res_ent=<int*> sage_malloc(sizeof(int)*self._N)
        if res_ent==NULL:
            raise MemoryError
        res = MyPermutation(length=self._N)
        for i in range(self.N()):
            j=self._entries[i]
            res_ent[j-1]=i+1
        res.set_entries(res_ent)        
        if res_ent<>NULL:
            sage_free(res_ent)
            res_ent=NULL   
        return res
    
    cdef int num_fixed_c(self):
        r"""
        The number of elementes which are fixed by self.
        """
        cdef int i,res
        res=0
        for i from 0<=i<self._N:
            if self._entries[i]==i+1:
                res+=1
        return res

    cpdef num_fixed(self):
        r"""
        The number of elementes which are fixed by self.
        """
        return self.num_fixed_c()

    cpdef fixed_elements(self):
        r"""
        The elements fixed by self.
        """
        cdef int i
        cdef list res=[]
        for i in range(self._N):
            if self._entries[i]==i+1:
                res.append(i+1)
        return res

    cpdef non_fixed_elements(self):
        r"""
        The elements not fixed by self.
        """
        cdef int i,num
        cdef list res=[]
        cdef int* non_fixed=NULL
        non_fixed = <int*>sage_malloc(sizeof(int)*self._N)
        num = 0
        for i in range(self._N):
            if self._entries[i]<>i+1:
                non_fixed[num]=i+1
                num+=1
            else:
                non_fixed[i]=0
        res = range(num)
        for i in range(num):
            res[i]=non_fixed[i]
        if non_fixed<>NULL:
            sage_free(non_fixed)
        return res
    
    def __mul__(self,other):
        if isinstance(other,MyPermutation):
            if other.N() == self.N():
                return MyPermutation._mult_perm(self,other)            
        elif isinstance(other,list):
            if len(other)==self.N():
                return MyPermutation._mult_list(self,other)
        #elif isinstance(other,type(self)):
            #print "other=",other
        elif hasattr(other,"list"):
            if len(other.list())==self.N():
                return MyPermutation._mult_list(self,other.list())
            
        raise TypeError,"Can not multiply self with {0}".format(other)

        
    def _mult_list(self,list other):
        cdef int *entries, N
        cdef MyPermutation res
        cdef int *res_ent=NULL
        N = len(other)
        if N<>self._N:
            raise NotImplementedError,"Can not multiply permutation on {0} letters with list of length {1}".format(self._N,N)             
        entries = <int*>sage_malloc(sizeof(int)*N)
        if entries==NULL: raise MemoryError
        for i from 0 <= i < N:
            entries[i]=<int>other[i]
        res_ent = <int*>sage_malloc(sizeof(int)*N)
        if entries==NULL: raise MemoryError
        _mult_perm_unsafe(N,self._entries,entries,res_ent)
        res = MyPermutation(length=N,init=0,check=0)
        res.set_entries(res_ent)        
        if entries<>NULL:
            sage_free(entries)
            entries=NULL
        if res_ent<>NULL:
            sage_free(res_ent)
            res_ent=NULL        
        return res

    cpdef _mult_perm(self,MyPermutation other):
        cdef MyPermutation res
        cdef int *res_ent=NULL
        cdef int N = self._N
        res_ent = <int*>malloc(sizeof(int)*N)
        if res_ent==NULL: raise MemoryError
        _mult_perm_unsafe(N,self._entries,other._entries,res_ent)
        res = MyPermutation(length=N,init=0,check=0)
        res.set_entries(res_ent)
        if res_ent<>NULL:
            free(res_ent)
            res_ent=NULL
        return res

 

    def square(self):
        cdef MyPermutation res
        cdef int i
        res = MyPermutation(length=self._N,init=0)
        cdef int* ent
        ent = <int*>sage_malloc(sizeof(int)*self._N)
        for i from 0 <= i < self._N:
            #print "i=",i
            #print "j=",self._entries[i]-1
            #print "k=",self._entries[self._entries[i]-1]
            #print "res[",i,"]0=",res._entries[i]
            ent[i]=self._entries[self._entries[i]-1]
            #print "res[",i,"]1=",res._entries[i]
        res.set_entries(ent)
        sage_free(ent)
        return res

    cpdef int N(self):
        return self._N

    cdef int N_c(self):
        return self._N

    def list(self):
        cdef list res
        cdef int i
        if self._list<>[]:
            return self._list
        res = []
        #print "in self.list!"
        if self._init==0:
            ## Then we have to init it
            for i in range(self._N):
                self._entries[i]=i+1
        for i in range(self._N):
            res.append(self._entries[i])
        self._list=res
        if self._verbose>0:
            print "self._list=",self._list
            print "type(list[0])=",type(self._list[0])
        return res

    cpdef conjugate(self,other):
        r"""
        Conjugate self with other, i.e. compute other*self*other**-1
        """
        cdef int l,i
        cdef MyPermutation res
        cdef int *entries=NULL
        if isinstance(other,type(self)):
            l = (<MyPermutation>other)._N
        elif isinstance(other,list):
            l = len(other)        
        if l<>self._N:
            raise TypeError,"Conjugate with permutation of same length!"
        res = MyPermutation(length=l)
        entries = <int *>sage_malloc(sizeof(int)*self._N)
        if isinstance(other,MyPermutation):
            _conjugate_perm(self._N,entries,self._entries,<int *>(<MyPermutation>other)._entries)
        elif isinstance(other,list):        
            _conjugate_perm_list(self._N,entries,self._entries,other)        
        else:
            if entries<>NULL:
                sage_free(entries)
            raise TypeError,"Can not conjugate self with {0}!".format(other)    
        res.set_entries(entries)
        if entries<>NULL:
            sage_free(entries)
            entries=NULL
        ## Reset cycles info
        
        return res

    cdef MyPermutation _conjugate(self,MyPermutation other):
        r"""
        COnjugate self by other: res = other*self*other**-1
        """
        cdef int* entries=NULL
        cdef MyPermutation res
        entries = <int *>sage_malloc(sizeof(int)*self._N)
        _conjugate_perm(self._N,entries,self._entries,other._entries)  
        res = MyPermutation(length=self._N)
        if entries<>NULL:
            res.set_entries(entries)        
            sage_free(entries)
            entries=NULL        
        return res

    cdef void _conjugate_ptr(self,int *entries,int *other):
        r"""
        Assume that self, entries ansd other are all allocated of the same length. 
        Reset the entries of self with the conjugate of other.
        """
        _conjugate_perm(self._N,self._entries,entries,other)  
        self._list=[]
        self._str=''
        self._hash=0
        

    cdef MyPermutation _conjugate_list(self,list other):
        r"""
        Conjugate self by other: res = other*self*other**-1
        Note: unsafe! do not use unless you are me! 
        """
        cdef int* entries=NULL
        cdef MyPermutation res
        entries = <int *>sage_malloc(sizeof(int)*self._N)
        _conjugate_perm_list(self._N,entries,self._entries,other)
        #_conjugate_perm(self._N,entries,self._entries,other._entries)  
        res = MyPermutation(length=self._N)
        if entries<>NULL:
            res.set_entries(entries)        
            sage_free(entries)
            entries=NULL        
        return res

    cdef MyPermutation _conjugate_vec(self,int * other):
        r"""
        Conjugate self by other: res = other*self*other**-1
        Note: unsafe! do not use unless you are me! 
        """
        cdef int* entries=NULL
        cdef MyPermutation res
        entries = <int *>sage_malloc(sizeof(int)*self._N)
        #_conjugate_perm_list(self._N,entries,self._entries,other)
        _conjugate_perm(self._N,entries,self._entries,other)  
        res = MyPermutation(length=self._N)
        if entries<>NULL:
            res.set_entries(entries)        
            sage_free(entries)
            entries=NULL        
        return res

    cdef void _conj_w_transp(self,int a, int b,int verbose=0):
        r"""
        Conjugate self by the transposition (a b)
        """
        cdef int i,ia,ib,j,sa,sb
        sa = self._entries[a-1]
        sb = self._entries[b-1]
        if verbose>0:
            print "Swap: {0} <-> {1}".format(a,b)
            print "a=",a, " s(a)=",sa
            print "b=",b, " s(b)=",sb
        if (a==sa and b==sb) or (sa==b and sb==a):
            pass
        elif a==sa:
            if verbose>0:
                print "a=sa: setting S[{0}] to {1}".format(b-1,b)
                print "      setting S[{0}] to {1}".format(a-1,sb)
            self._entries[b-1]=b   ## Now b is fixed 
            self._entries[a-1]=sb
            for i in range(self._N): 
                if self._entries[i]==b and i<>b-1:
                    self._entries[i]=a
                    break
        elif b==sb:
            self._entries[a-1]=a
            self._entries[b-1]=sa
            for i in range(self._N): 
                if self._entries[i]==a and i<>a-1:
                    self._entries[i]=b
                    break
        else:
            self._entries[a-1] = sb
            self._entries[b-1] = sa
            ia=-1; ib=-1
            for i in range(self._N): #from 1<= i <= self._N:
                if ia < 0 and self._entries[i]==a:
                        ia = i
                if ib<0 and self._entries[i]==b:
                        ib = i
                if ia>=0 and ib>=0:
                    break
            self._entries[ia]=b
            self._entries[ib]=a
        # for i in range(1,self._N+1): #from 1<= i <= self._N:
        #     #print "i=",i
        #     if i==a:
        #         #print "i=a"
        #         if ib==a:
        #             self._entries[i-1]=b
        #         elif ib==b:
        #             self._entries[i-1]=a
        #         else:
        #             self._entries[i-1]=ib
        #     elif i==b:
        #         #print "i=b"
        #         if ia==a:
        #             self._entries[i-1]=b
        #         elif ia==b:
        #             self._entries[i-1]=a
        #         else:
        #             self._entries[i-1]=ia
        #     else:
        #         #print "i not in {a,b}"
        #         if self._entries[i-1]==a:
        #             #print "have s(i)=a"
        #             self._entries[i-1]=b
        #         elif self._entries[i-1]==b:
        #             #print "have s(i)=b"
        #             self._entries[i-1]=a

            #print "s(i)=",self._entries[i-1]
        # Conjugate the list only if necessary
        self._list = []
        self._str=''
        self._hash=0
        cdef int reorder
        ## We have to recompute any cycle involving a and b
        ## For simplicity we just remove the cycles and recompute them if necessary
        ## A start at a more efficient way is below...
        if self._num_cycles>0 and self._cycles<>NULL:
            self._num_cycles=0
            sage_free(self._cycles)
            self._cycles = NULL
            # ia = -1; ib=-1
            # for i in range(self._N):
            #     if ia<0 and self._cycles[i]==a:
            #         ia = i
            #     if ib<0 and self._cycles[i]==b:
            #         ib = i
            #     if ib>=0 and ia>=0:
            #         break        
            # # now have the indices of a and b
            # self._cycles[ia]=b
            # self._cycles[ib]=a
            # possibly have to reorder the cycles...
            # reorder = 0
            # for i in range(self._num_cycles):
            #     for j in range(self._cycle_lens[i]):

        self._cycles_list = []
        self._cycles_permutations = []
        self._cycles_ordered_as_list = []
        #print "self=",self

    cpdef conjugate_with_transposition(self,int a,int b,int verbose=0):
        r"""
        Conjugate in place
        """
        #self._list = []
        self._conj_w_transp(a,b,verbose)


    cpdef is_conjugate_to(self,MyPermutation other,int ret_perm=0):
        r"""
        Check if self is conjugate to other.
        """
        cdef int N = other._N
        cdef MyPermutation p
        if ret_perm==1:
            p = MyPermutation(length=self._N)
        if N<>other._N:
            if ret_perm==1:
                return 0,p
            else:
                return 0
        cdef list ctA,ctB
        ctA = self.cycle_type()
        ctB = other.cycle_type()
        if ctA<>ctB:
            if ret_perm==1:
                return 0,p
            else:
                return 0
        if ret_perm==0:
            return 1
        lA = self.cycles_ordered_as_list()
        lB = other.cycles_ordered_as_list()
        return 1,get_conjugating_perm_list(lA,lB)

    
    cdef int eq_c(self,MyPermutation other):
        if self._N<>other._N:
            return 0
        cdef int i
        for i from 0 <= i< self._N:
            if self._entries[i]<>other._entries[i]:
                return 0
        return 1


    cpdef int eq(self,MyPermutation other):
        if self._N<>other._N:
            return 0
        cdef int i
        for i from 0 <= i< self._N:
            if self._entries[i]<>other._entries[i]:
                return 0
        return 1

#    cpdef cycles(self, int ordered = 1,int type=1): #'list'):
    cpdef cycles_as_permutations(self): #'list'):
        r"""

        Gives the cycles of self as list of permutations.
        Type = 0 'perm' or 1 'list'. If the type is list we can also choose ordered=1 to return a list with
        cycles given by increasing length.
        """
        cdef int i
        cdef list l,pl
        cdef MyPermutation p0
        if self._cycles_permutations == []:
            l = self.cycles_as_lists()
            for pl in l:
                p0 = MyPermutation(length=self._N)
                if len(pl)==1:
                    self._cycles_permutations.append(p0)
                else:
                    # Have a non-trivial cycle
                    for i in range(len(pl)-1):
                        p0._entries[pl[i]-1]=pl[i+1]
                    p0._entries[pl[-1]-1]=pl[0]
                    self._cycles_permutations.append(p0)
        return self._cycles_permutations

    cpdef cycles_as_lists(self):
        r"""
        Returns the cycles of self as a list of lists.
        """
        cdef list tmp
        cdef int i,j,bd
        #print "in cycles_as_list!"
        if self._cycles_list == []:
            self._cycles_list_is_ordered=0
            if self._num_cycles == 0 or self._cycles==NULL or self._cycle_lens==NULL:
                self.set_cycle_info()
            tmp = []
            for j in range(self._cycle_lens[0]):
                tmp.append(self._cycles[j])
            self._cycles_list=[tmp]
            bd = 0
            for i in range(1,self._num_cycles):
                tmp = []
                bd += self._cycle_lens[i-1]
                for j in range(self._cycle_lens[i]):
                    tmp.append(self._cycles[j+bd])
                self._cycles_list.append(tmp)
            #print "cycles_list = ",self._cycles_list
        return self._cycles_list
    
    def cycles(self,type='list',order=0,verbose=0):
        r"""
        Returns the cycles of self as a list of lists or permutations.
        order = 0 => ordered lexicographically
        order = 1 => ordered according to length first
        """
        if verbose>0:
            print "type=",type
            print "order=",order
            print "is_ordered=",self._cycles_list_is_ordered
            print "cycles_as_list=",self._cycles_list
        if type=='list':
            if self._cycles_list == []:
                self.cycles_as_lists()
            if order == self._cycles_list_is_ordered:
                return self._cycles_list
            if order==1:
                self._cycles_list.sort(key = lambda x: len(x))
                self._cycles_list_is_ordered=1
            elif order==0:
                self._cycles_list=[]
                self.cycles_as_lists()                
            return self._cycles_list
        else:
            return self.cycles_as_permutations()

            

    cpdef cycles_ordered_as_list(self):
        r"""
        This gives the cycles of self in terms of a permutation, i.e.
        (1 2 3)(6 4 5)  is represented by [1,2,3,6,4,5]
        """
        cdef int i,j
        cdef list cycles
        #p = MyPermutation('(1)(2 9 4)(3)(5)(6 8 10)(7)')
        #p0 = MyPermutation('(1)(2 4 9)(3)(5)(6 8 10)(7)')
        #if self==p or self==p0:
        #    print "self=",self
        #    print "self cycles=",self._cycles_list
        #    print "ordered cycles=",self._cycles_ordered_as_list
        if self._cycles_ordered_as_list==[]:
            cycles = self.cycles(order=1,verbose=0)
            #print "cycles=",cycles
            self._cycles_ordered_as_list = []
            for i in range(len(cycles)):
                for j in range(len(cycles[i])):
                    self._cycles_ordered_as_list.append(cycles[i][j])
        return self._cycles_ordered_as_list
            



    cpdef set_cycle_info(self,int reset=0):
        r""" Collect information about all cycles of self.
         Allocate and initialise the vecotr cycles which contain the
         cycles of self in one row, i.e. [[1,2,3],[5],[4,2]] would be represented by the
         vector [1,2,3,5,4,2]


         """
        if self._num_cycles>0 and self._cycles<>NULL and reset==0:
            return ## already set
        if reset==1:
            if self._cycles<>NULL:
                sage_free(self._cycles)
                self._cycles=NULL
            self._num_cycles=0
            self._cycles_list = []
            self._cycles_ordered_as_list = []
            self._cycles_permutations = []
        self.set_cycle_info_c()

    cdef int set_cycle_info_c(self):       
        cdef int i,j,bd_old,ii,cycle_bd,tmp
        self._cycles = NULL
        self._cycles = <int*>sage_malloc(sizeof(int)*(self._N))
        if self._cycles == NULL:
            raise MemoryError
        for i in range(self._N):
            self._cycles[i]=0        
        self._num_cycles = 0
        self._cycle_lens = <int*>sage_malloc(sizeof(int)*(self._N))
        if self._cycle_lens == NULL:
            raise MemoryError
        for i in range(self._N):
            self._cycle_lens[i]=0
        cycle_bd=0
        ii=0
        for i in range(1,self._N+1):
            if _is_in_list(self._cycles,i,self._N):
                continue
            bd_old = cycle_bd
            cycle_bd+=1
            if cycle_bd > self._N:
                break
            self._cycles[cycle_bd-1]=i
            if self._cycles[cycle_bd-1] > self._N:
                break
            tmp = self._entries[self._cycles[cycle_bd-1]-1]
            while tmp<>i and tmp <= self._N and cycle_bd <= self._N:
                cycle_bd+=1
                #print "cycle)(",ii,cycle_bd-1,")=",tmp
                self._cycles[cycle_bd-1]=tmp
                tmp = self._entries[tmp-1]
            cycle_len = cycle_bd - bd_old
            if  ii >= self._N:
                break
            self._cycle_lens[ii]=cycle_len
            #print "len(",ii,")=",cycle_len
            ii+=1
            if cycle_bd > self._N or ii >= self._N:
                break
        self._num_cycles = ii
        #for i in range(0,self._N):
        #    print "cycles[",i,"]=",self._cycles[i]
        return 0

    def cycle_lens(self):
        res = []
        if self._num_cycles == 0:
            self.set_cycle_info()
        for i in range(self._num_cycles):
            res.append(self._cycle_lens[i])
        return res


    # def to_cycles(self,ordered=0):
    #     r"""
    #     Gives the cycles of self as list of lists.
    #     """
    #     if not hasattr(self,"_cycles"):
    #         self._cycles_list=[]
    #     if self._cycles_list == []:
    #         if not self._list or self._entries==NULL:
    #             self.list()
    #         #self._cycles = perm_to_cycle_c(self._N,self._entries)
    #         self._cycles_list = self._perm_to_cycles()
    #     if ordered==1:
    #         self._cycles_list.sort(key = lambda x: len(x))
    #     return self._cycles_list


    cpdef num_cycles(self):
        r"""
        Return the number of cycles of self.
        """
        if self._num_cycles>0:
            return self._num_cycles
        else:
            self._num_cycles = num_cycles_c(self._N,self._entries)
        return self._num_cycles
        
    
    def cycle_type(self):
        r"""
        Return the cycle type of self.
        
        """
        if self._cycle_type == []:
            if self._cycles == NULL or self._num_cycles==0:
                self.set_cycle_info()            
            lens = []
            for i in range(self._num_cycles):                
                lens.append(self._cycle_lens[i])
            lens.sort()
            self._cycle_type = lens
        return self._cycle_type
    
    def cycle_tuples(self):
        r"""
        Gives the cycles of self as list of tuples.
        """
        return map(tuple,self.cycles())

    def cycles_non_triv(self):
        r"""
        Returns non-trivial cycles only, as in the standard SymmetricGroup etc.
        """
        res=list()
        for c in map(tuple,self.cycles()):
            if len(c)>1:
                res.append(c)
        return res


def EndOfList(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


cdef class MyPermutationIterator(SageObject):
    r"""
    Creates an iterator over a set of permutations. 

    INPUT:

     - 'order' -- only consider permutations of this order.
     - 'num_fixed' -- only consider permutations with this number of fixed elements.
     - ''verbose' -- verbosity
     
    """
    
    #    def __cinit__(self,N,verbose=0,order=0,num_fixed=-1,fixed_pts=None,**kwds):
    def __cinit__(self,N,int order=0,int num_fixed=-1,list fixed_pts=[],int verbose=0,is_subiterator=0,map_from=0,map_at_most_to=0):
        r"""

        If a=map_from and b=map_at_least_to with a,b,>0 then we iterate over R s.t. R(a)>b
        
        """
        self._N=N
        self._cur = Integer(0)
#        mpz_set_ui(self._cur.value,0)
        self._verbose=verbose
        self._current_piv=0
        self._order=order
        self._num_fixed = num_fixed 
        self._got_list=0
        self._num=0 
        self._max_num=Integer(0)
        self._max_num=Integer(0) 
        self._current_state_c=NULL
        self._current_state_o=NULL
        self._list_of_perms=NULL
        self._current_perm=NULL
        self._fixed_pts=NULL
        if fixed_pts<>[]:
            self._num_fixed=len(fixed_pts)
        if num_fixed>0:
            self._num_fixed=num_fixed
        self._verbose=0
        if verbose>0:
            # Check for verbose = a*4 + b*2 + c*1
            for i in range(3):
                if verbose % 2 == 1:
                    self._verbose+=1
                    verbose = verbose - 1
                verbose = verbose / 2            
            
        if self._verbose>0:
            print "in cinit! verbose=",self._verbose
            print "fixdpts=",fixed_pts
        self._map_from = <int>map_from
        self._map_at_most_to =<int> map_at_most_to
        if self._verbose>0:
            print "map_from=",map_from,self._map_from
            print "map_at_most_to=",map_at_most_to,self._map_at_most_to   
        if self._map_from>0 and self._map_at_most_to<=0:
            raise ValueError,"Inconsistent condition on R's! {0}: {1}".format(map_from,map_at_most_to)


        self._c_new()

    def _c_new(self):
        ## By default we initialize the list of permutations to just one permutation, the identity
        cdef int i
        self._c_alloc_list_of_perms(1)
        self._current_perm = self._list_of_perms[0]   # the identity perm. used later to store current perm
        self._current_state_c = NULL
        self._current_state_c= <int *> sage_malloc(2*sizeof(int)*self._N) # used for swap
        if self._current_state_c == NULL: raise MemoryError
        self._current_state_o=   self._current_state_c + self._N # used for swap
        self._fixed_pt_free_iterator_labels = NULL
        self._fixed_pts=NULL
        mpz_set_ui(self._cur.value,0) # The index of the current permutation.
        for i in range(self._N):
            self._list_of_perms[0][i]=i+1
            self._current_state_c[i]=-1
            self._current_state_o[i]=1
        if self._verbose>2:
            print "in cinit!"
            print "num_fixedpts=",self._num_fixed
        if self._num_fixed>=0:
            self._fixed_pts=<int*>sage_malloc(sizeof(int)*self._num_fixed)
            if self._fixed_pts==NULL: raise MemoryError
            
        
    def __init__(self,N,order=0,num_fixed=-1,fixed_pts=None,verbose=0,is_subiterator=0,map_from=0,map_at_most_to=0):
        r"""

        INPUT:
        - verbose = binary additive.
           -- 2 => output from iterators
           -- 4 => output from misc.
           -- 8 => output from constructors / destructors
           -- 16 => extremely detailed output 
        """
        cdef int i,n
        ## We adjust the verbosity to desired level
        if self._verbose>2:
            print "in init! N={0}".format(self._N)
            print "fixed_pts=",fixed_pts
        self._has_fixed_pt_free_iterator = 0
        self._fixed_pt_free_iterator = None
        #        self._fixed_pt_free_iterator_labels=list()
        self._is_subiterator=is_subiterator
        if fixed_pts<>None and fixed_pts<>[]:
            if isinstance(fixed_pts,list) and self._num_fixed<>len(fixed_pts):
                raise ValueError,"got inconsistent fixed point data! fixpts={0} and num_fixed_pts={1}".format(fixed_pts,self._num_fixed)                     
            self._max_num=num_permutations_without_fixed_pts(N-self._num_fixed)
            if self._max_num==0:
                self._num=0
            for i in range(self._num_fixed):
                n = fixed_pts[i]
                if n<1 or n>N:
                    raise ValueError,"got inconsistent fixed point data! fixpts={0} and num_fixed={1}".format(fixed_pts,self._num_fixed)
                self._fixed_pts[i]=n
                if self._verbose>0:
                    print "setting fixed_pt[",i,"]=",self._fixed_pts[i]
        elif self._num_fixed>=0:
            self._max_num=copy(num_permutations_without_fixed_pts(N-self._num_fixed))
            for i from 0 <= i < self._num_fixed:
                self._fixed_pts[i]=i+1
        elif fixed_pts==[]:
            self._max_num=copy(num_permutations_without_fixed_pts(N))
            self._num_fixed = 0
        else:
            self._num_fixed = -1
            self._max_num = Integer(factorial(N))
        if self._verbose>2:
            print "fixed_pts in init after=",fixed_pts
            print "self._fixed_pts=",self.fixed_pts()
            print "num=",self.num()
        if self._num_fixed>0 and self._order>0:
            if (N - max(self._num_fixed,0)) % self._order <>0:
                self._max_num = Integer(0)
                self._num = 0
#                    raise ValueError,"got inconsistent order and fixed point data! num. fixpts={0} order={1}".format(self._num_fixed,self._order)
                if verbose>0:
                    print "Got inconsistent order and fixed point data! num. fixpts={0} order={1}".format(self._num_fixed,self._order)
        self.set_fixed_pt_free_iterator()
        #if self._num_fixed<>0:
            #    # Since we start with the identity permutation by default we have to go to the next
        #if self._verbose>0:
        #        print "go to next!"
        if not self._is_subiterator:
            self._goto_first()  
        if self._verbose>2:
            print "self._fixed_pts=",self.fixed_pts()


    def set_fixed_pt_free_iterator(self):
        if self._num_fixed<=0:
            return 
        self._fixed_pt_free_iterator = MyPermutationIterator(self._N-self._num_fixed,num_fixed=0,order=self._order,verbose=self._verbose,is_subiterator=1,map_from=self._map_from,map_at_most_to=self._map_at_most_to)
        #self._fpfree_labels=list()
        self._has_fixed_pt_free_iterator=1
        if self._fixed_pt_free_iterator_labels==NULL:
            self._fixed_pt_free_iterator_labels=<int*>sage_malloc(self._N*sizeof(int))
        if self._fixed_pt_free_iterator_labels==NULL:
            raise MemoryError
        #self._fixed_pt_free_iterator_labels=list()
        cdef int i ,j,addi,ni
        ni = 0
        for i in range(1,self._N+1):
            addi=1
            for j in range(self._num_fixed):
                if i==self._fixed_pts[j]:
                    addi=0
                    break
            if addi==1:
                self._fixed_pt_free_iterator_labels[ni]=i
                ni+=1
        if ni<self._N-self._num_fixed:
            raise ArithmeticError,"Too few elements in fixed point free iterator!"
            #for i in range(ni,self._N-self._num_fixed):
            #    self._fixed_pt_free_iterator_labels[i]=0
        self._fixed_pt_free_iterator.set_verbosity(self._verbose)
         #self._fixed_pt_free_iterator.set_labels(l)
         
    def free_iterator(self):
        return self._fixed_pt_free_iterator

    def free_iterator_labels(self):
        res = []
        for i in range(self._fixed_pt_free_iterator.N()):
            res.append(self._fixed_pt_free_iterator_labels[i])
        return res

    def __iter__(self):
        self._rewind()
        return self

    def __getitem__(self,index):
        cdef MyPermutation p
        cdef int i
        if self._list_of_perms<>NULL and self._num>=index:
            p=MyPermutation(length=self._N)
            p.set_entries(self._list_of_perms[index])
            return p
        for i in range(1,index+1): #from 1<=i<=index:
            self.next()
        return self.current_perm()

    

    def __dealloc__(self):
        if self._verbose>2:
            print "in __dealloc__"
        self._c_dealloc()

    def _c_dealloc(self):
        r"""
        Accessible from the outside...
        """
        cdef int i
        if self._verbose>2:
            print "in c_dealloc!"
        if self._list_of_perms<>NULL:
            if self._verbose>2:
                print "num_in_deal=",self._num
                print "list_o_p="
                printf("%p ", self._list_of_perms)
            for i in range(self._num):
                if self._verbose>2:
                    print "list_o_p[",i,"]="
                    printf("%p ", self._list_of_perms[i])
                if self._list_of_perms[i]<>NULL:
                    sage_free(self._list_of_perms[i])
                    self._list_of_perms[i]=NULL
            sage_free(self._list_of_perms)
            self._list_of_perms=NULL
        if self._fixed_pts<>NULL:
            if self._verbose>2:
                print "fixed_pts="
                printf("%p ", self._fixed_pts)
            sage_free(self._fixed_pts)
            self._fixed_pts=NULL
        if self._current_state_c<>NULL:
            self._current_state_o=NULL
            sage_free(self._current_state_c)
            self._current_state_c=NULL
        if self._fixed_pt_free_iterator_labels<> NULL:
            sage_free(self._fixed_pt_free_iterator_labels)
            self._fixed_pt_free_iterator_labels=NULL
        if self._verbose>2:
            print "end of dealloc!"
            print "self_l_o_p=",printf("%p ",self._list_of_perms)
        
    def _c_alloc_list_of_perms(self,int num=0):
        r"""
        We allocate self._list_of_perms in this routine *only*.
        This makes it easier to control memory errors...
        This routine and the next are also the *only* place where self._num is to be set after initialization!

        """
        if self._verbose>2:
            print "-----------------------------------------"
            print "want to allocate list_perm with num=",num
        cdef int i
        if num==0:
            self._list_of_perms=NULL
            return
        # first deallocate if we have already allocated
        if self._num>0:
            if self._list_of_perms<>NULL: 
                if self._verbose>2:
                    print "we already have list of length:",self._num
                for i from 0 <= i <= self._num-1:
                    if self._list_of_perms[i]:
                        sage_free(self._list_of_perms[i])
                        self._list_of_perms[i]=NULL
                sage_free(self._list_of_perms)
                self._list_of_perms=NULL
        if self._max_num>0 and num > self._max_num:
            self._num = self._max_num ## We never need more than this
            if self._verbose>2:
                print "Set to maxnum:",self._max_num
        else:
            self._num=num
        self._list_of_perms = <int**> sage_malloc(sizeof(int*)*self._num)
        if not self._list_of_perms: raise MemoryError
        if self._verbose>2:
            print "we will allocate list at address:",printf("%p ", self._list_of_perms)
        for i in range(self._num):
            self._list_of_perms[i]=NULL
            self._list_of_perms[i] = <int*> sage_malloc(sizeof(int)*self._N)
            if not self._list_of_perms[i]: raise MemoryError
            if self._verbose>2:
                print " allocated l_o_p[",i,"]=",printf("%p ", self._list_of_perms[i])
        if self._verbose>2:
            print "we now have allocated list of length:",self._num            
            print " at address:",printf("%p ", self._list_of_perms)

    def _c_dealloc_last_perm(self):
        if self._num<=0:
            return 
        if self._list_of_perms:
            if self._list_of_perms[self._num-1]<>NULL:
                sage_free(self._list_of_perms[self._num-1])
                self._list_of_perms[self._num-1]=NULL
            self._num=self._num-1
            if self._verbose>2:
                print "We have decreased the list to length:",self._num

    def verbose(self):
        return self._verbose

    def set_verbosity(self,verbose):
        self._verbose = verbose
    
    def current_piv(self):
        return self._current_piv

    def __repr__(self):
        s="Permutation iterator over %s elements " % self._N
        if self._num_fixed>0:
            s+="\n, which fixes the points "+str(self.fixed_pts())+"."
        elif self._num_fixed==0:
            s+=" without fixed points."
        if self._list_of_perms <>NULL:
            s+="\n List of %s permutations available "% self.num()
        if self._max_num==0:
            s+="This is an empty iterator!"
        else:
            s+="\n Current permutation: %s " % self.current_perm()
        return s
    
    def fixed_pts(self):
        l=list()
        cdef int i
        for i from 0 <= i < self._num_fixed:
            l.append(self._fixed_pts[i])
        return l

    def num_fixed_pts(self):
        return self._num_fixed

    cpdef long num(self):
        r"""
        Returns self._num -- the index of the last permutation in self._list_of_perms.

        """
        if self._verbose>4:
            print "pointer_to_list=",printf("%p ", self._list_of_perms)
            print "num is now=",self._num
        return self._num

    def max_num(self):
        return self._max_num
    
    def current_perm(self):
        r"""
        Returns current permutation of self.
        """
        if self._current_perm == NULL:
            raise ArithmeticError," We do not have a current permutation!"
        res = []
        if mpz_cmp_si(self._max_num.value,0)==0:
            return None
        for i in range(self._N):
            res.append(self._current_perm[i])
        return MyPermutation(res)

    def cur(self):
        return self._cur

    def show_current_state(self):
        s=" self_current_state_c:("
        for i in range(self._N):
            s+=str(self._current_state_c[i])+" "
        s+=")\n self_current_state_o:("
        for i in range(self._N):
            s+=str(self._current_state_o[i])+" "
        s+=")"
        print s

    

    def order(self,set_to=0):
        if set_to<=0:
            return self._order
        else:
            self._order = set_to
    ### Iterator methods

    def __contains__(self,x):
        l=x.list()
        if set(l)<>set(range(1,self._N+1)):
            if self._verbose>0:
                print "x is not a permutation on 1,...,{0} but on {1}".format(self._N,set(l))
            return False
        if fixed_elements(l)<>self.fixed_pts():
            if self._verbose>0:
                print "x has fixed points: {0}. self has fix pts:{1}".format(fixed_elements(l)<>self.fixed_pts())
            return False
        if self.order()>0 and self.order()<>x.order():
            if self._verbose>0:
                print "x has order: {0}. self has order:{1}".format(x.order(),self.order())
            return False
        
        return True

    

    cpdef _goto_first(self):
        r"""
        Go to the first permitted permutation.
        """
        if self._verbose>0:
            print "in _goto_first! N=",self._N," cur=",self._cur
            print "map_from0=",self._map_from,self._map_at_most_to
        cdef int done = 1
        cdef int t=0
        if mpz_cmp_ui(self._cur.value,0)<>0:
            self._rewind()
        if self._current_perm<>NULL:
            if self._order>0:
                done = _is_of_order(self._N,self._current_perm,self._order)
            if done==1 and self._num_fixed>=0:
                done = _num_fixed_eq(self._current_perm,self._N,self._num_fixed)
            if done==1 and self._num_fixed>=0:
                done = _set_fixed_eq(self._current_perm,self._N,self._fixed_pts,self._num_fixed)
            if done==1 and self._map_from>0:
                if self._current_perm[self._map_from-1]>self._map_at_most_to:
                    done = 0
        if done==0: # If the current perm is not in the desired resulting list we go forward
            if self._fixed_pt_free_iterator:
                self._fixed_pt_free_iterator._goto_first()
                self.set_current_from_ffree()                
                mpz_add_ui(self._cur.value,self._cur.value,1)
            else:
                t = self._next()
                if t<>0:
                    return t
        else:
            mpz_set_ui(self._cur.value,1)
        if self._verbose>1:
            print "done=",done
            print "First allowed one is: ",self.current_perm()
        return 0
    ##                 
    def __next__(self):
        if self._verbose>0:
            print "__next__ N=",self._N," cur=",self._cur
        if mpz_cmp(self._max_num.value,self._cur.value)<0: #self._cur > self._max_num:
            raise StopIteration
        if mpz_cmp_ui(self._cur.value,0)==0: #self._cur > self._max_num:
            self._rewind()
            self._goto_first()
        res = self.current_perm()
        if self._verbose>0:
            print "cur =",res
        t = self._next()
        if t<>0:
            mpz_add_ui(self._cur.value,self._max_num.value,1)
            #self._cur = self._max_num+1
        if self._verbose>0:
            print "t =",t    
        return  res #self.current_perm()
        

    cpdef _next(self):
        cdef int* fpc
        cdef int i
        if self._verbose>0:
            print "in _next! N=",self._N," cur=",self._cur
            print "self._num_fixed=",self._num_fixed
            print "self._fixed=",self.fixed_pts()
        i = mpz_cmp(self._max_num.value,self._cur.value)
        if i==0:
            mpz_add_ui(self._cur.value,self._cur.value,1)   
            return 1 # If we are at the last element we simply return 1
        elif i < 0: 
            # if we are after the last element we raise StopIteration #if self._cur>self._max_num:
            if self._verbose>0:
                print "raising stop iteration!"
            raise StopIteration #IndexError,"End of list!"
        if self._has_fixed_pt_free_iterator==0:
            return self._next_free()
        else:
            if self._verbose>0:
                print "Have fixed point free iterator: ",self._fixed_pt_free_iterator," cur=",self._fixed_pt_free_iterator._cur
            mpz_add_ui(self._cur.value,self._cur.value,1) # self._fixed_pt_free_iterator._cur
            i = self._fixed_pt_free_iterator._next() #_free()
            if i<>0:
                return i
            self.set_current_from_ffree()                
        if self._verbose>0:
            for i in range(self._N):
                print "cur_p[",i,"]=",self._current_perm[i]
            print "cur =",self._cur
        return 0


    cdef set_current_from_ffree(self):
        cdef int i,j,k
        cdef int* fpc = self._fixed_pt_free_iterator._current_perm
        if self._verbose>0:
            print "setting current_perm from fixed-point free!"
        for i in range(self._num_fixed):
            self._current_perm[self._fixed_pts[i]-1]=self._fixed_pts[i]
            if self._verbose>0:
                print "fixed_pt[",i,"]=",self._current_perm[self._fixed_pts[i]-1]
        for i in range(self._fixed_pt_free_iterator._N):
            j = self._fixed_pt_free_iterator_labels[i]-1
            k = fpc[i]-1
            self._current_perm[j]=self._fixed_pt_free_iterator_labels[k]
            if self._verbose>0:
                print "cur_perm[",self._fixed_pt_free_iterator_labels[i]-1,"]=",self._current_perm[self._fixed_pt_free_iterator_labels[i]-1]
            
    ###
    ### Functions which iterate over the permutations one by one.
    ###
    cpdef _next_free(self):
        r"""
        Sets the current perm to the next perm.
        Used for either a set of permutations without fixed points or a list with all different fixed points.
        """
        cdef int *lista
        cdef int start,tmp2,num,i,k,t,done
        cdef long tmp,n,f
        cdef Integer mmax
        cdef mpz_t mpzmmax
        cdef unsigned long ii
        #cdef unsigned long lmmax
        if mpz_cmp_si(self._max_num.value,0)==0:
            return 1
        if self._verbose>0:
            print "in _next_free! N=",self._N," cur=",self._cur
        #if self._cur>=self._max_num:
        if mpz_cmp(self._max_num.value,self._cur.value)<=0:
            mpz_add_ui(self._cur.value,self._cur.value,1)
            return 1
            #raise StopIteration
        if self._list_of_perms<>NULL and self._num>1:
            # If we are iterating over a whole list
            if self._verbose>0:
                print "fetching next  from list: num=",self.num()
            mpz_add_ui(self._cur.value,self._cur.value,1) #            self._cur += 1
            self._current_perm=self._list_of_perms[self._cur]
        else:
            #if self._verbose>0:
            #    print "computing next! N=",self._N," cur=",self._cur
            #    print "max_num=",self._max_num
            #    print "current_perm=",self.current_perm()
            if self._cur==0:
                #if self._cur==1 and mpz_cmp_ui(self._max_num.value,1)==0:
                #    return
                if debug>3: 
                    print "before reset!"
                    self.show_current_state()
                reset_swap(self._N,self._current_state_c,self._current_state_o)
                if self._verbose > 3:
                    print "after reset!"
                    self.show_current_state()
                    #reset_swap(self._N,self._current_state_c,self._current_state_o)
            if self._order<=1 and self._num_fixed<0 and self._map_from<=0:
                if debug>3:
                    print "before next swap"
                    self.show_current_state()
                t=next_swap(self._N,self._current_state_c,self._current_state_o)
                if debug>3:
                    print "after next swap"
                    self.show_current_state()
                if debug>0:
                    print "swap t=",t
                mpz_add_ui(self._cur.value,self._cur.value,1) #self._cur = self._cur+1
                if t<0:
                    if self._verbose>2:
                        print "current_perm0=",self.current_perm()
                     #mpz_clear(mpzmmax)
                    return -1 
                #raise ArithmeticError,"Got index <0 from swap!"
                tmp = self._current_perm[t]
                self._current_perm[t] = self._current_perm[t+1]
                self._current_perm[t+1] = tmp
            else:
                
                done=0
                ii=1
                mpz_init(mpzmmax)
                mmax=my_factorial(self._N)
                mpz_set(mpzmmax,mmax.value)
                while done == 0: # and ii<=mmax:
                    if self._verbose>2: 
                        print "before swap:"
                        self.show_current_state()
                    t=next_swap(self._N,self._current_state_c,self._current_state_o)
                    if self._verbose>2: 
                        print "after swap:"
                        self.show_current_state()
                        print "swap t,ii=",t,ii
                    if t<0:
                        done = 0 # raise ArithmeticError,"Got index <0 from swap!"
                        break
                    tmp = self._current_perm[t]
                    self._current_perm[t] = self._current_perm[t+1]
                    self._current_perm[t+1] = tmp                    
                    if debug>0: 
                        if self._verbose>2 and self._current_perm[0]==1 and self._current_perm[1]==2:
                            print "current_perm1=",self.current_perm()
                    done = 1
                    if self._verbose > 0:
                        print "map_from1=",self._map_from
                        print "map at most to =",self._map_at_most_to                     
                    if self._map_from>0:
                        if self._verbose > 1:
                            print "map_from2=",self._map_from
                            print "map at most to =",self._map_at_most_to
                        if self._current_perm[self._map_from-1]>self._map_at_most_to:
                            done = 0
                            if self._verbose > 1:
                                print "Continue since R({0})={1}>{2}".format(self._map_from,self._current_perm[self._map_from-1],self._map_at_most_to)
                    if done==1 and self._order>0:
                        done = _is_of_order(self._N,self._current_perm,self._order)
                        #else:
                        #    done=1
                    if self._verbose > 1:
                        print "current_perm2=",self.current_perm()
                    #    print "of correct order=",done
                    if self._num_fixed==0 and done==1:
                        done = _num_fixed_eq(self._current_perm,self._N,self._num_fixed)
                    if done==0:
                        #if ii>mmax:
                        if mpz_cmp_si(mpzmmax,ii)<0:
                            break
                        ii=ii+1
                mpz_clear(mpzmmax)
                if done==1:
                    #print "cur=",self._cur
                    mpz_add_ui(self._cur.value,self._cur.value,1) #self._cur = self._cur+1
                else: ## We have reached the end. Should start over from beginning.
                    #print "PI=",self
                    return 1
                    #raise StopIteration  #,"You are at the end of the list!"
                    #raise ArithmeticError,"Could not find permutation with correct properties!"
                #if self._verbose > 3:
                #    print "perm is of order ",self._order
        return 0
        #if self._verbose > 0:
        #    print "end of next!",self._order        
        

    ## def _next_old(self):
    ##     r"""
    ##     Sets the current perm to the next perm.
    ##     """
    ##     cdef int *lista
    ##     cdef int start,tmp2,num,i,k,ii,t,done
    ##     cdef long tmp,n
    ##     cdef Integer mmax
    ##     #cdef int *fixed_pts
    ##     if self._max_num==0:
    ##         return
    ##     if self._cur>=self._max_num:
    ##         raise IndexError
    ##     if self._list_of_perms<>NULL and self._num>1:
    ##         # If we are iterating over a whole list
    ##         if self._verbose>0:
    ##             print "fetching next  from list: num=",self.num()
    ##         self._cur += 1
    ##         self._current_perm=self._list_of_perms[self._cur]

                
    ##     else:
    ##         if self._verbose>0:
    ##             print "computing next!"
    ##             print "cur=",self._cur
    ##             print "max_num=",self._max_num
    ##         if self._cur==0:
    ##             ##
    ##             if self._max_num==1 and self._cur==1:
    ##                 # Then we do not need to iterate more...
    ##                 return
    ##             if self._verbose > 3:
    ##                 print "before reset!"
    ##                 self.show_current_state()
    ##             reset_swap(self._N,self._current_state_c,self._current_state_o)
    ##             if self._verbose > 3:
    ##                 print "after reset!"
    ##                 self.show_current_state()
    ##                 #reset_swap(self._N,self._current_state_c,self._current_state_o)
    ##         if self._order<=1 and self._num_fixed<0:
    ##             if self._verbose > 3:
    ##                 print "before next swap"
    ##                 self.show_current_state()
    ##             t=next_swap(self._N,self._current_state_c,self._current_state_o)
    ##             if self._verbose > 3:
    ##                 tmplist=list();tmplist1=list()
    ##                 print "after next swap"
    ##                 self.show_current_state()
    ##             if self._verbose > 3:
    ##                 print "swap t=",t
    ##             if t<0:
    ##                 raise ArithmeticError,"Got index <0 from swap!"
    ##             tmp = self._current_perm[t]
    ##             self._current_perm[t] = self._current_perm[t+1]
    ##             self._current_perm[t+1] = tmp
    ##             self._cur = self._cur+1
    ##         else:
    ##             done=0
    ##             ii=0
    ##             mmax = my_factorial(self._N)
    ##             #print "here: ii=",ii
    ##             while done == 0 and ii<=mmax:
    ##                 if self._verbose > 3:
    ##                     print "before swap:"
    ##                     self.show_current_state()
    ##                 t=next_swap(self._N,self._current_state_c,self._current_state_o)
    ##                 if self._verbose > 3:
    ##                     print "after swap:"
    ##                     self.show_current_state()
    ##                 #print "swap t,ii=",t,ii
    ##                 ii=ii+1
    ##                 if t<0:
    ##                     done = 0 # raise ArithmeticError,"Got index <0 from swap!"
    ##                     break
    ##                 tmp = self._current_perm[t]
    ##                 self._current_perm[t] = self._current_perm[t+1]
    ##                 self._current_perm[t+1] = tmp                    
    ##                 #if self._verbose > 2:
    ##                 if self._verbose>2 and self._current_perm[0]==1 and self._current_perm[1]==2:
    ##                     print "current_perm=",self.current_perm()
    ##                 done = _is_of_order(self._N,self._current_perm,self._order)
    ##                 if self._verbose > 3:
    ##                     print "of correct order=",done
    ##                 if self._num_fixed>=0 and done==1:
    ##                     if self._num_fixed>0:
    ##                         done = _set_fixed_eq(self._current_perm,self._N,self._fixed_pts,self._num_fixed)
    ##                         if (self._verbose > 2 and done) or self._verbose>3:
    ##                             print "fixed are equal=",done
    ##                         #done = _num_fixed_eq(self._current_perm,self._N,self._num_fixed)
    ##                     else:
    ##                         done = _num_fixed_eq(self._current_perm,self._N,self._num_fixed)
    ##                         if (self._verbose > 2 and done) or  self._verbose>3:
    ##                             print "number of fixed points are {0}! done={1}".format(self._num_fixed,done)
    ##             if done==1:
    ##                 self._cur = self._cur+1
                    
    ##             else: ## We have reached the end. Should start over from beginning.
    ##                 print "PI=",self
    ##                 raise IndexError  #,"You are at the end of the list!"
    ##                 #raise ArithmeticError,"Could not find permutation with correct properties!"
    ##             if self._verbose > 3:
    ##                 print "perm is of order ",self._order
    ##     if self._verbose > 0:
    ##         print "end of next!",self._order   
            


    cpdef _rewind(self):
        r"""
        Sets current permutation to the identity and resets the states.
        If we have fixed points or oorder we also iterate forward to the first
        correct permutation.
        """
        if self._verbose>0:
            print "rewind!"
        mpz_set_ui(self._cur.value,0)
        if self._fixed_pt_free_iterator:
            self._fixed_pt_free_iterator._rewind()
        reset_swap(self._N,self._current_state_c,self._current_state_o)
        self._c_alloc_list_of_perms(1) # self._num=1
        cdef int i
        for i in range(self._N):
            self._list_of_perms[0][i]=i+1
        self._current_perm = self._list_of_perms[0]
        if self._order>0 or self._num_fixed>=0:
            self._goto_first()
            if self._verbose>0:
                print "Current=",self.current_perm()
                
    ###
    ### Functions which generates a whole list of permutations
    ###
    def __list__(self):
        return self.list()

    cpdef list_new(self):
        cdef list res=[]
        cdef MyPermutation x
        for x in self:
            res.append(x)
        return res


    def list_new2(self):
        cdef list res=[]
        cdef MyPermutation x
        for x in self:
            res.append(x)
        return res


    cpdef list(self):
        cdef int i,j
        cdef list l
        if self._max_num==0:            
            return []
        if self._got_list==0:
            self._get_list_of_perms()
        if self._list_of_perms<>NULL:
            res=range(self._num)
            for i in range(self._num):
                if self._list_of_perms[i]<>NULL:
                    l = range(self._N) #list()
                    for j in range(self._N):
                        #l.append(self._list_of_perms[i][j])
                        l[j]=self._list_of_perms[i][j]
                    res[i]=l
                else:
                    res[i]=[]
            return res
        else:
            return []
        
    cpdef _get_list_of_perms(self):
        r"""
        Populates self._list_perms with all permutations. Should only be used for small indices.
        """
        if self._N>12:
                raise ValueError, "Do not use for N>12!"
        cdef int *lista
        cdef int i
        cdef long num
        if self._got_list==1: #list_of_perms<>NULL and self._num>2:            
            if self._verbose>0:
                print "we already have a list of length : %s " %  self.num()
            return 
        self._got_list=1
        self._c_alloc_list_of_perms(factorial(self._N))
        if self._verbose > 1:
            print "allocated num=",self._num
        lista = <int*> sage_malloc(sizeof(int)*self._N)
        if not lista: raise MemoryError
        for i in range(self._N): 
            lista[i]=i+1
        mpz_set_ui(self._cur.value,0)
        num = factorial(self._N) #deepcopy(self.num())
        if self._verbose>1:
            print "will get next_perm recursive!"
        self._get_next_permutation_recursive(0,0,num,lista)        
        # make sure we got the correct number of permutations
        if  mpz_cmp(self._max_num.value,self._cur.value) <0 and self._order<=1:
            if self._verbose > 1:
                print "number of perms=",self.num()
            raise ArithmeticError,"Something wrong with generation of permutations!"
        mpz_set_ui(self._cur.value,0)
        if self._verbose > 1:
            print "num in get list=",self.num()
        sage_free(lista)
        lista=NULL

    cdef _get_next_permutation_recursive(self,int tmp, int start,long num,int *lista):        
        r"""
        Recursive function which will produce all permutations.
        TODO: Figure out how to use this to just get the next permutation.
        """
        cdef int i,j,k,test_o,test_f
        cdef int *fixed_pts
        if self._verbose>1:
            print "lista=",printf("%p ", lista)
            tmpl=list()
            for i from 0 <= i <= self._N-1:
                tmpl.append(lista[i])
            print "current perm=",tmpl
            print "tmp=",tmp
            print "start=",start
            print "num=",num
            print "cur=",self._cur
            print "map_from3=",self._map_from,self._map_at_most_to
        if self._order<=1 and self._num_fixed<0:
            for i in range(self._N):
                self._list_of_perms[self._cur][i]=lista[i]            
            mpz_add_ui(self._cur.value,self._cur.value,1)
        else:
            # First 
            if self._map_from>0:
                test_o  = 1
                if self._current_perm[self._map_from-1]>self._map_at_most_to:
                    test_o = 0
                if self._verbose > 1:
                    print "Continue since R({0})={1}>{2}".format(self._map_from,self._current_perm[self._map_from-1],self._map_at_most_to)
            # Then test if of correct order.
            if test_o==1:
                test_o =_is_of_order(self._N,lista,self._order)
                if self._verbose>3:
                    print "test_o=",test_o
                    print "num_fixed=",self._num_fixed
            if self._num_fixed>0 and test_o==1:
                if self._verbose>3:
                    print "compare fixed points of perm with ",self.fixed_pts()
                test_f  = _set_fixed_eq(lista,self._N,self._fixed_pts,self._num_fixed)
                if self._verbose>3:
                    print "fixed points agree :",test_f
            elif self._num_fixed==0 and test_o==1:
                test_f =_num_fixed_eq(lista,self._N,self._num_fixed)
                if self._verbose>3:
                    print "test_f=",test_f
            if test_f==1 and test_o==1:
                for i in range(self._N):
                    if self._verbose>3:
                        print "list_of_perm[",self._cur,"][",i,"]=",lista[i]
                    self._list_of_perms[self._cur][i]=lista[i]                            
                mpz_add_ui(self._cur.value,self._cur.value,1)
                if self._verbose>3:
                    print "self._cur=",self._cur
                    print "self._num=",self._num
                    print "self._max_num=",self._max_num
                # If we have enough we stop
                if mpz_cmp(self._max_num.value,self._cur.value) <= 0:
                    num=0
            elif test_o<>1:
                ## If we have correct number of fixed points but the order does not match then we allocated
                ## a too long list to begin with and we now reduce the list
                test_f  = _set_fixed_eq(lista,self._N,self._fixed_pts,self._num_fixed)
                if test_f==1:
                    if self._verbose>3:
                        print "have num to:",self.num()
                    self._c_dealloc_last_perm()
                    if self._verbose>3:
                        print "decreased num to:",self.num()

        if num>0:
            num=num-1
            if self._verbose>3:
                print "start=",start
            if start < self._N-1:
                i = self._N-2
                while i>=start:
                    for j from i+1 <= j < self._N:
                        # !!! swap lista(i) och lista(j) for i+1<=j<=self._N-1
                        if self._verbose>3:
                            print "lista_i[",i,"]=",lista[i]
                            print "lista_j[",j,"]=",lista[j]
                        tmp=lista[i]
                        lista[i]=lista[j]
                        lista[j]=tmp
                        if self._verbose>3:
                            print "lista_i[",i,"]=",lista[i]
                            print "lista_j[",j,"]=",lista[j]
                            print "tmp=",tmp
                            for k from 0 <= k <= self._N-1:
                                print "lista[",k,"]=",lista[k]
                            print "num_before=",num        
                        self._get_next_permutation_recursive(tmp,i+1,num,lista)
                        if self._verbose>3:
                            print "num_after=",num        
                        if num<=0:
                            #print "setting start to :",i-1
                            self._current_piv = i-1
                            return
                        # !! rotate lista(i:n) to the left
                    tmp=lista[i]
                    for k from i <= k <= self._N-2:
                        lista[k]=lista[k+1]
                    lista[self._N-1]=tmp
                    i = i-1
                # return tmp
        else:
            pass
        #    print "start==",start


    cpdef get_perms_test(self):
        r"""
        More or less directly copied from permutation_cython.pyx
        A less efficient way of generating all permutations, but we can use 'next_swap' to iterate over permutations.
        """
        cdef int *c, *o, N, m,t,i
        cdef int *perm,perm_new
        cdef list T
        if self._N <= 1:
            return []
        if self._N > 12:
            raise ValueError, "Cowardly refusing to enumerate the permutations on more than 12 letters."
        if not self._current_state_c: raise MemoryError
        if not self._current_state_o: raise MemoryError
        self._c_alloc_list_of_perms(factorial(self._N))
        reset_swap(self._N,self._current_state_c,self._current_state_o)
        for i in range(self._N):
            self._list_of_perms[0][i]=i+1
        if self._verbose>1:
            for i in range(self._N):                
                print "start_list[",0,"]=",self._list_of_perms[0][i]            
        if self._verbose>1:
            s=""; ss=""
            for i in range(self._N):
                s+=" ",self._current_state_o[i]
                ss+=" ",self._current_state_c[i]
            print "o[0]="+s
            print "c[0]="+ss
        for m in range(1,self._num):
            t=next_swap(self._N,self._current_state_c,self._current_state_o)
            if self._verbose>1:
                s=""; ss=""
                for i from 0<= i< self._N:
                    s+=" ",self._current_state_o[i]
                    ss+=" ",self._current_state_c[i]
                print "o[0]="+s
                print "c[0]="+ss

            for i in range(self._N):
                if i==t:
                    self._list_of_perms[m][t]=self._list_of_perms[m-1][t+1]
                elif i==t+1:
                    self._list_of_perms[m][t+1]=self._list_of_perms[m-1][t]
                else:
                    self._list_of_perms[m][i]=self._list_of_perms[m-1][i]



cdef class CycleCombinationIterator(Parent):
    r"""
    Class implementing a sequence of permutations given by combinations of cycles in a given permutation.
    This should be a small number since we construct all permutations.

    NOTE: Currently only works with cycles of the same length.
    
    """

    def __cinit__(self,MyPermutation R):
        r"""
        C-Init self.
        """
        self._cycles = []
        self._cycle_types = NULL
        self._N = R.N()
        self._length = 1
        self._num_cycles=0
        self._dbase = 1
        self._cache = []
        cdef int j
        self._got = 0
        cdef list tmp
        cdef int jj
        self._c_new()
        for c in R.cycles_as_permutations():
            j = c.order()
            if j==1:
                continue
            tmp = []
            for jj in range(1,j):
                if jj == 1:
                    tmp.append(c)
                elif jj==2:
                    tmp.append(c.square())
                else:
                    tmp.append(c.pow(jj))
            self._cycles.append(tmp)
            self._cycle_types[self._num_cycles]= j
            self._length = self._length*j
            self._num_cycles = self._num_cycles + 1
            if j > self._dbase:
                self._dbase = j

    def _c_new(self):
        r"""
        Allocate stuff.
        """
        if self._cycle_types<>NULL:
            sage_free(self._cycle_types)
        self._cycle_types = <int*>sage_malloc(sizeof(int)*self._N)
        if self._cycle_types == NULL:
            raise MemoryError
                
    
    def __init__(self,MyPermutation R):
        r"""
        Initialize self.
        """
        pass

    def __dealloc__(self):
        if self._cycle_types <> NULL:
            sage_free(self._cycle_types)
        self._cycle_types=NULL
    
    
    def __repr__(self):
        r""" print self (as a list)"""
        s=""
        for i in range(self._length):
            s+=str(self.permutation_nr(i))+' '
        return s
    
    cpdef length(self):
        r"""
        The number of elements in self.
        """
        return self._length

    cpdef list(self):
        r"""
        A list of elements of self.
        """
        cdef int i
        for i in range(len(self._cache),self._length):
            self._cache.append(self.permutation_nr_c(i))
            self._got = self._length
        return self._cache

    
    def get_permutation_nr(self,M):
        r"""
        The M'th element of self.
        """
        return self.permutation_nr(M)
    
    cpdef MyPermutation permutation_nr(self,int M):
        r"""
        The M'th element of self.
        """
        M = M % (self._length + 1)
        cdef int i
        if M >= self._got:
            for i in range(self._got,M+1):
                self._cache.append(self.permutation_nr_c(i))
            self._got = M+1
        return self._cache[M]

    @cython.cdivision(True)
    cdef MyPermutation permutation_nr_c(self,int M): #,MyPermutation p):
        r"""
        The M'th element of self.
        """
        cdef MyPermutation res
        res = MyPermutation(length=self._N)
        cdef int i, M1,j1
        for i in range(self._num_cycles):
            M1 = M % self._dbase            
            #res = res._mult_perm(<MyPermutation>(self._cycles[i]).pow(M1))
            j1 = M1 % self._cycle_types[i]
            if j1>0:
                #res = res._mult_perm(<MyPermutation>(self._cycles[i][j1-1]))
                _mult_perm_unsafe(self._N,res._entries,(<MyPermutation>(self._cycles[i][j1-1]))._entries,res._entries)                
            M = M / self._dbase # integer division, i.e. M = (M-M1)/self._dbase
        return res

    @cython.cdivision(True)
    cdef int permutation_nr_c_ptr(self,int M,int* res): #,MyPermutation p):
        r"""
        The M'th element of self.
        """
        #cdef MyPermutation res
        #res = MyPermutation(length=self._N)        
        cdef int i, M1,j1
        for i in range(self._N):
            res[i]=i+1
#            print "res[i]=",res[i]
        for i in range(self._num_cycles):
            M1 = M % self._dbase            
            #res = res._mult_perm(<MyPermutation>(self._cycles[i]).pow(M1))
            j1 = M1 % self._cycle_types[i]
            if j1>0:
                #res = res._mult_perm(<MyPermutation>(self._cycles[i][j1-1]))
                _mult_perm_unsafe(self._N,res,(<MyPermutation>(self._cycles[i][j1-1]))._entries,res)                
            M = M / self._dbase # integer division, i.e. M = (M-M1)/self._dbase
        return 0
        #return res

    
### Helper functions

cdef int _is_of_order(int n,int *plist, int order):
    r"""
    Check if list is a permutation of given order.
    """
    cdef int i,j,tmp
    cdef list l
    if order<1:
        return 1
    elif order==3:
        for i in range(n):
            if plist[plist[plist[i]-1]-1]<>i+1:
                return 0
        return 1
    elif order==2:
        for i in range(n):
            if plist[plist[i]-1]<>i+1:
                return 0
        return 1
    for i in range(n):
        tmp = plist[i]-1
        if tmp>=n or tmp<0:
            l=list()
            for j in range(n):
                l.append(plist[j])
            raise ValueError,"Not a valid permutation in {0}".format(l)
        for j in range(2,order+1): #from 2 <= j <= order:        
            tmp = plist[tmp]-1
        if tmp<>i:
            return 0
    return 1

cdef int _num_fixed_eq(int *list,int n, int num_fixed):
    r"""
    Check if list is a permutation with given number of fixed-points.
    """
    cdef int i,tmp
    if num_fixed < 0:
        return 1
    tmp = 0
    for i from 0 <= i < n:
        if list[i]==i+1:
            tmp+=1
        if tmp > num_fixed:
            return 0
    if tmp<>num_fixed:
        return 0
    return 1


cdef int _set_fixed_eq(int *list,int n, int *fixed_pts, int num_fixed_pts):
    r"""
    Check if list is a permutation with given number of fixed-points.
    """
    cdef int i,tmp
    if num_fixed_pts<0:
        return 1
    tmp = 0
    for i from 0 <= i < n:
        if list[i]==i+1:
            if _is_in_list(fixed_pts,i+1,num_fixed_pts)==0:
                return 0
            tmp+=1
        if tmp > num_fixed_pts:
            return 0
    if tmp<>num_fixed_pts:
        return 0
    return 1

cpdef are_conjugate_perm(MyPermutation A,MyPermutation B):
   r"""
   Check if the permutation A is conjugate to B
   """
   cdef int N = A._N
   if N<>B._N:
       return 0
   cdef list ctA,ctB
#   A.set_cycle_info(reset=1)
#   B.set_cycle_info(reset=1)   
   ctA = A.cycle_type()
   ctB = B.cycle_type()
   if ctA<>ctB:
       return 0
   #print "A=",A
   #print "B=",B   
   lA = A.cycles_ordered_as_list()
   lB = B.cycles_ordered_as_list()
   #print "lA=",lA
   #print "lB=",lB
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

cdef MyPermutation  get_conjugating_perm_ptr_unsafe(int mu, int* Al,int* Bl):
    r"""
    Find permutation which conjugates Al to Bl, where Al and Bl are lists given by cycle structures.
    """
    cdef int* cperm=NULL
    cdef MyPermutation res
    cdef int i
    cperm = <int*> sage_malloc(sizeof(int)*mu)
    if cperm==NULL:
        raise MemoryError
    for i in range(mu):
        cperm[Al[i]-1]=Bl[i]
    res = MyPermutation(length=mu,init=0,check=0)
    res.set_entries(cperm)
    if cperm <> NULL:
        sage_free(cperm)
    return res

##
## Algorithms for working with cycles in terms of pointers (and lists)
##


cpdef perm_to_cycle(perm):
    r"""
    Writes a permutation, given as a list of integers, in terms of cycles 
    """
    cdef int N,i,last,l
    cdef int *cycle_lens=NULL
    cdef int *cycle=NULL
    cdef int *permv=NULL
    N = len(perm)
    cycle = <int*> sage_malloc(sizeof(int)*N)
    if cycle==NULL: raise MemoryError
    cycle_lens = <int*> sage_malloc(sizeof(int)*N)
    if cycle_lens==NULL: raise MemoryError
    permv = <int*> sage_malloc(sizeof(int)*N)
    if permv==NULL: raise MemoryError
    for i  in range(N):
        permv[i]=perm[i]
        cycle[i]=0
    _to_cycles(N,permv, cycle, cycle_lens)
    res = list()
    cy = list()
    last = 0
    for i in range(N):
        l = cycle_lens[i]
        if l==0:
            continue
        cy=list()
        for j in range(last,last+l):
            cy.append(cycle[j])
        last = last+l
        res.append(cy)
        if last>=N:
            break
    if permv<>NULL:
        sage_free(permv)
    if cycle_lens<>NULL:
        sage_free(cycle_lens)
    if cycle<>NULL:
        sage_free(cycle)
    return res


   
cdef int num_cycles_c(int N,int *perm):
    cdef int i,n,t0,t1
    cdef int* used = NULL
    used = <int*> sage_malloc(sizeof(int)*N)
    if used==NULL: raise MemoryError
    for i in range(N):
        used[i]=0
    n = 0
    for i in range(N):
        if used[i]==1:
            continue
        t0 = perm[i]
        t1 = t0
        used[i]=1
        for j in range(N):
            used[t1-1]=1
            t1 = perm[t1-1]
            if t1==t0:
                break
        n+=1
    if used<>NULL:
        sage_free(used)
    return n

cdef int perm_to_cycle_c(int N,int *perm,int *cycle,int *cycle_lens):
    cdef int *permv=NULL
    cdef list res,cy
    cdef int i,j,k,last,l
    permv = <int*> sage_malloc(sizeof(int)*N)
    if permv==NULL: raise MemoryError
    for i in range(N):
        permv[i]=perm[i]
        cycle[i]=0
    _to_cycles(N,permv, cycle, cycle_lens)
    if permv<>NULL:
        sage_free(permv)
    return 1

cdef _to_cycles(int N, int *perm, int *cycle, int *cycle_lens):
    r"""
    Returns cycles and cycle lengths of a permutation.

    INPUT:
      - N integer
      - perm -- allocated array of integers of length N 
      - cycle -- allocated int* of length N
      - cycle_lens -- allocated int* of length N
    """
    cdef int i,bd_old,ii
    cdef int cycle_bd,tmp
    cycle_bd=0
    ii=0
    for i in range(N):
        cycle_lens[i]=0 
    for i in range(1,N+1):
        if _is_in_list(cycle,i,N):
            continue
        bd_old = cycle_bd
        cycle_bd+=1
        if cycle_bd > N:
            break
        cycle[cycle_bd-1]=i
        if cycle[cycle_bd-1] > N:
            break
        tmp = perm[cycle[cycle_bd-1]-1]
        while tmp<>i and tmp <= N and cycle_bd <=N:
            cycle_bd+=1
            cycle[cycle_bd-1]=tmp
            tmp = perm[tmp-1]
        cycle_len = cycle_bd - bd_old
        if  ii >= N:
            break
        cycle_lens[ii]=cycle_len
        ii+=1
        if cycle_bd > N or ii >= N:
            break



        
cdef void _mult_perm_unsafe(int N,int *left,int *right,int* res):
        r"""
        This is the standard order of multiplication of permutations.
        res = left*right : i-> right(left(i))
        This is 'unsafe' since no check of validity of data or arrays is performed.
        """
        cdef int i
        for i in range(N):
            res[i]=right[left[i]-1]

cdef void _mult_perm_right_unsafe(int N,int *left,int *right):
        r"""
        This is the standard order of multiplication of permutations.
        left = left*right : i-> right(left(i))
        This is 'unsafe' since no check of validity of data or arrays is performed.
        """
        cdef int i
        for i in range(N):
            left[i]=right[left[i]-1]



cpdef conjugate_perm(perma,permb):
    r"""
    Conjugate perma with permb

    """
    cdef int N = len(perma)
    cdef int *a, *b, *c
    cdef int i
    a = <int*> sage_malloc(sizeof(int)*N)
    if not a: raise MemoryError
    b = <int*> sage_malloc(sizeof(int)*N)
    if not b: raise MemoryError
    c = <int*> sage_malloc(sizeof(int)*N)
    if not c: raise MemoryError    
    for i from 0 <= i <N:
        a[i]=perma[i]
        b[i]=permb[i]
    _conjugate_perm(N,c,a,b)
    res = list()
    for i from 0 <= i <N:
        res.append(c[i])
    sage_free(a)
    sage_free(b)
    sage_free(c)
    return res

cdef void _conjugate_perm(int N,int* res,int *a,int* b):
        r"""
        res = b*a*b^-1
        """
        cdef int i,j
        for i from 0 <= i < N:
            for j from 0 <= j < N:
                if b[j]<>i+1:
                    continue
                else:
                    break
            # now i = b^-1(j+1)
            if a[j]<1 or a[j]>N:
                raise IndexError,"Array index out of bounds! a[{0}]={1}".format(j,a[j])
            res[i]=b[a[j]-1]

cdef void _conjugate_perm_list(int N,int *res,int *a,list b):
        r"""
        res = b*a*b^-1
        """
        cdef int i,j
        
        for i from 0 <= i < N:
            for j from 0 <= j < N:
                if b[j]<>i+1:
                    continue
                else:
                    break
            # now i = b^-1(j+1)
            res[i]=b[a[j]-1]
            

cdef void square_my_perm_to_pt(MyPermutation perm, int *res):
    cdef int i
    for i from 0 <= i < perm._N:
        res[i]=perm._entries[perm._entries[i]-1]

### Factorial ###

cdef dict dict_of_factorials={}
        
cpdef Integer my_factorial(int n):
    cdef mpz_t res
    cdef Integer e
    cdef long nn = n
    e = dict_of_factorials.get(n)
    if not e:
        e = Integer(1)
        if n>20:
            mpz_init(res)
            mpz_fac_ui(res,nn) #factorial_mpz(n,res)
            mpz_set(e.value, res)
            mpz_clear(res)
            #return e #mpz_get_d(res)
        else:
            nn = factorial(n) 
            if isinstance(nn,int):
                mpz_set_si(e.value,nn)
            elif isinstance(nn,long): #PyLong_Check(nn):                
                mpz_set_pylong(e.value,nn)
        dict_of_factorials[n]=e
    return e

cdef void factorial_mpz(int n,mpz_t result):
    # result should be initialized
    #mpz_init_set_ui(result,1)
    mpz_set_ui(result,1)
    cdef int i
    for i from 1<= i <= n:
        mpz_mul_ui(result,result,i)

cdef long factorial(int n):
    if n>20:
        return -1
    cdef long i,result
    result = 1
    for i from 1<= i <= n:
        result*=i
    return result
###
cdef MyPermutation perm_power_c(MyPermutation perm,int k):
    cdef MyPermutation res,fac
    #cdef int* entries_res 
    #res = MyPermutation(length=self._N)
    cdef int p,m,j,o
    if k==0:
        return MyPermutation(length=p._N)
    if k==1:
        return perm
    if k==-1:
        return perm.inverse()
    if k<0:
        k = -k
        fac = perm.inverse()
    else:
        fac = perm
    o = perm.order()
    k = k % o
    res = MyPermutation(length=perm._N)
    for p,m in ZZ(k).factor():
        if m==1:
            if p==2:
                res=res*fac.square()
            else:
                for j from 0<=j<p:
                    res=res*fac
        else:
            for j from 0<=j<m:
                res=res*fac.pow(p)
    return res
# cdef _power(int k, int N, int o=0, int* entries,int* res):
#     cdef int i
#     if k==0:
#         for i in range(N):
#             res[i]=i+1
#     elif k==1:
#         for i in range(N):
#             res[i]=entries[i]
#     else:
#         if o>0:
#             if k >1:
#                 inv = 0 
#             if k > o:
#                 k = k % o
#                 inv = 1
#          if k%2==0:
#              tmp=self.square()
#              res=tmp
#             #for i in range(k/2-1):
#              for i from 0<=i<k/2-1:
#                  res=res*tmp
#              else:
#                  tmp=self.square()
#                  res=self
#                  for i in xrange((k-1)/2):
#                      res=res*tmp
#              return res
#          elif k<0:
#              k = -k
#              o = self.order()
#              if k>o:
#                  k = k % o
#              return self.inverse().__pow__(k)
#          elif k==0:
#              return MyPermutation(length=self.N(),init=1)  


### Fixed points ###
cpdef num_fixed(list x):
    cdef int i,num_fixed
    num_fixed=0
    for i from 1 <= i < len(x)+1:
        if x[i-1]==i:
            num_fixed+=1
    return num_fixed


cpdef fixed_elements(list x):
    cdef int i,num_fixed
    res=list()
    for i from 1 <= i < len(x)+1:
        if x[i-1]==i:
            res.append(i)
    return res

###

### Check transitivity

cpdef are_transitive_permutations(MyPermutation pS,MyPermutation pR,int ret_maps=0,int verbose=0):
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

         sage: E=Permutations(4)([1,2,4,3]); E.cycles()
         [(1,), (2,), (3, 4)]
         sage: R=Permutations(4)([2,1,3,4]); R.cycles()
         [(1, 2), (3,), (4,)]
         sage: are_transitive_permutations(E,R)
         False
         sage: R=Permutations(4)([2,3,1,4]); R.cycles()
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
    cdef int *gotten=NULL
    cdef int N = pS.N()
    cdef int *maps=NULL
    cdef int nmaps
    gotten = <int *>sage_malloc(sizeof(int)*N)
    if gotten==NULL: raise MemoryError
    #Sl = <int *>sage_malloc(sizeof(int)*N)
    #if not Sl: raise MemoryError
    #Rl = <int *>sage_malloc(sizeof(int)*N)
    #if not Rl: raise MemoryError    
    #for i in range(N):
    #    Sl[i]=pS._entries[i]; Rl[i]=pR._entries[i]        
    res = are_transitive_perm_c(pR._entries,pS._entries,gotten,N,verbose)
    #sage_free(Rl)
    #sage_free(Sl)
    if gotten <> NULL:
        sage_free(gotten)
    return res


cpdef are_transitive_permutations_wmaps(MyPermutation pS,MyPermutation pR,int verbose=0):
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

         sage: E=Permutations(4)([1,2,4,3]); E.cycles()
         [(1,), (2,), (3, 4)]
         sage: R=Permutations(4)([2,1,3,4]); R.cycles()
         [(1, 2), (3,), (4,)]
         sage: are_transitive_permutations(E,R)
         False
         sage: R=Permutations(4)([2,3,1,4]); R.cycles()
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
    cdef int *gotten
    cdef int N=pS.N()
    cdef int *maps=NULL
    cdef int nmaps
    cdef MyPermutation pT=pS*pR
    if pS.order()<>2 or pR.order()<>3:
        raise ValueError,"Need pS of order 2 and pR of order 3"
    gotten = <int *>sage_malloc(sizeof(int)*N)
    if not gotten: raise MemoryError
    cdef int num,num_old
    cdef int i,j,k,x
    cdef dict maps_list={}
    cdef int* cycle_lens=NULL
    pT.cycles(order=0)  ## make sure we have the cycle containing 1 first
    cdef int numc = pT.num_cycles()
#    cdef list cycle_lens = pT.cycle_lens() #<int*>sage_malloc(sizeof(int)*numc)
    if cycle_lens==NULL:
        raise MemoryError
    for i in range(numc):
        cycle_lens[i]=pT.cycle_lens()[i]
    for i in range(N):
        gotten[i]=0
    cdef SL2Z_elt A,R,S,T
    A = SL2Z_elt(1,0,0,1)
    R = SL2Z_elt(0,-1,1,1)
    S = SL2Z_elt(0,-1,1,0)
    T = SL2Z_elt(1,1,0,1)
    x = 1
    maps_list[1]=A
    gotten[0]=1
    for i in range(1,cycle_lens[0]):
        x = pT._entries[x-1]
        gotten[x-1]=x
        if verbose>0:
            print "gotten[{0}]={1}".format(x-1,x)
        maps_list[x]=SL2Z_elt(1,i,0,1)
        #print "maps{0}={1}".format(x,maps_list[x])
    num = cycle_lens[0]

    if verbose>0:
        print "gotten=",print_vec(N,gotten)
        for i in range(N):
            print "gotten[{0}]={1}".format(i,gotten[i])
    cdef int cycle_bd=cycle_lens[0]
    for j in range(1,numc):
        if verbose>0:
            print "cycle[{0}]={1}".format(j,pT.cycles()[j])
        ## want to connect with the next cycle of pT
        num_old = num
        for i in range(N):
            if gotten[i]==0:
                continue
            if verbose>0:
                print "Checking {0}".format(i+1)
                print "S({0})={1}".format(i+1,pS._entries[i])
                print "R({0})={1}".format(i+1,pR._entries[i])
                print "R^2({0})={1}".format(i+1,pR._entries[pR._entries[i]-1])                
            if _is_in_list(pT._cycles+cycle_bd,pS._entries[i],cycle_lens[j]):
                x = pS._entries[i]
                A = maps_list[i+1]*S
            elif _is_in_list(pT._cycles+cycle_bd,pR._entries[i],cycle_lens[j]):
                x = pS._entries[i]
                A = maps_list[i+1]*R
            elif _is_in_list(pT._cycles+cycle_bd,pR._entries[pR._entries[i]-1],cycle_lens[j]):
                x =  pR._entries[pR._entries[i]-1]
                A = maps_list[i+1]*R*R
            else:
                continue
            if verbose>0:
                print "Here A=",A
            gotten[x-1]=x
            maps_list[x]=A
            num+=1
            for k in range(1,cycle_lens[j]):
                A = A*T
                x = pT._entries[x-1]
                gotten[x-1]=x
                maps_list[x]=A
                num+=1
            if num==N:
                break
        if num==N:
            break
        if num_old==num:
            # we didn't get to the next cycle so we are not transitive
            if cycle_lens<>NULL:
                sage_free(cycle_lens)
            return 0,{}
        cycle_bd+=cycle_lens[j]
    if cycle_lens<>NULL:
        sage_free(cycle_lens)
    return 1,maps_list
                
    #     num_old = num
        
    #     for k in range(num_old):
    #         x = gotten[k]
    #         if verbose>0:
    #             print "gotten[{0}]={1}".format(k,x)
    #             print "Sl[x-1]=",pS._entries[x-1]
    #         A = S*A #tmp_list.append('S')
    #         if _is_in_list(gotten,pS._entries[x-1],num)==0:
    #             gotten[num]=pS._entries[x-1];  num+=1
    #             maps_list[pS._entries[x-1]]=A #copy(tmp_list)
    #         x = pS._entries[x-1]
    #         if verbose>0:
    #             print "gotten=",print_vec(N,gotten)
    #             print "x=",x
    #             print "Rl[x-1]=",pR._entries[x-1]
    #         A = R*A #tmp_list.append('R')
    #         if _is_in_list(gotten,pR._entries[x-1],num)==0:
    #             gotten[num]=pR._entries[x-1];  num+=1
    #             maps_list[pR._entries[x-1]]=A #copy(tmp_list)
    #             #tmp_list.append('R')
    #         A = R*R*A
    #         x=pR._entries[x-1]
    #         if verbose>0:
    #             print "x1=",x
    #         if _is_in_list(gotten,pR._entries[x-1],num)==0:
    #             gotten[num]=pR._entries[x-1];  num+=1
    #             maps_list[pR._entries[x-1]]=A #copy(tmp_list)
    #         if verbose>0:
    #             print "num=",num
    #             print "gotten[end]=",print_vec(N,gotten)
    #             print "maps_list = ",maps_list
    #     if num == num_old:
    #         if num<>N:
    #             sage_free(gotten)
    #             return 0,{}
    #         else:
    #             sage_free(gotten)
    #             return 1,maps_list
    # sage_free(gotten)
    # return 0,{}


cdef int are_transitive_perm_c(int *Rl,int *Sl,int *gotten, int n,int verbose=0):
    cdef int num,num_old
    cdef int j,x
    if gotten==NULL: raise MemoryError
    gotten[0]=1; num=1
    if Rl[0]<>1:
        gotten[1]=Rl[0]
        num=2        
    if verbose>0:
        for j in range(n):
            print "Rl[{0}]={1}".format(j,Rl[j])
        for j in range(n):
            print "Sl[{0}]={1}".format(j,Sl[j])
    for j in range(n):
        num_old = num
        for k in range(num_old):
            x = gotten[k]
            if verbose>0:
                print "gotten[{0}]={1}".format(k,x)
                print "Sl[x-1]=",Sl[x-1]
            if _is_in_list(gotten,Sl[x-1],num)==0:
                gotten[num]=Sl[x-1];  num+=1
            x = Sl[x-1]
            if verbose>0:
                print "gotten=",print_vec(n,gotten)
                print "x=",x
                print "Rl[x-1]=",Rl[x-1]
            if _is_in_list(gotten,Rl[x-1],num)==0:
                gotten[num]=Rl[x-1];  num+=1
            x=Rl[x-1]
            if verbose>0:
                print "x1=",x
            if _is_in_list(gotten,Rl[x-1],num)==0:
                gotten[num]=Rl[x-1];  num+=1
            if verbose>0:
                print "num=",num
                print "gotten[end]=",print_vec(n,gotten)
        if num == num_old:
            if num<>n:
                return 0
            else:
                return 1
    return 0

###

cdef int _is_in_list(int *lista,int y,int num):
    cdef int i
    for i from 0 <= i < num:
        if lista[i]==y:
            return 1
    return 0

### Printing ###

cdef print_vec(int n,int *list):
    s="["+str(list[0])
    cdef int i
    for i in range(1,n):
        s+=","+str(list[i])
    s+="]"
    return s

cdef void _print_vec(int n,int* l):
    cdef int i
    cdef str s #char* str
    s=""
    for i from 0 <= i < n:
        s+=str(l[i])+" "
    print s
    


# def print_cycles(perm):
#     if isinstance(perm,str):
#         cycles = perm_to_cycle(eval(perm))
#     else:
#         cycles = perm_to_cycle(perm)
#     s =str(cycles)
#     return s

# cdef _print_cycles(int n,int *perm):
#     cdef MyPermutation tmp
#     tmp=MyPermutation(length=n)
#     tmp.set_entries(perm)
#     s =str(tmp.cycles())
#     return s

# def print_ppair(a,b):
#     s=print_cycles(a)
#     s+=" "+print_cycles(b)
#     return s

###


cdef _are_eq_vec(int n,int *a,int *b):
    cdef int i,res
    for i from 0 <= i< n:
        if a[i]<>b[i]:
            return 0
    return 1

def verbosity(v,level):
    return has_binkey(v,level)

def has_binkey(x,n):
    r"""
    Test if x has the binary digit number n or not
    """
    if x==0:
        return 0
    if not isinstance(x,basestring):
        s=bin(x)
    else:
        s=x
    d=s.split('b')[1]
    if len(d)<=n:
        return 0
    if d[-1-n]=='1':
        return 1
    else:
        return 0
    




cdef dict num_perm_dict={}
# Calculate number of permutations 
#@staticmethod
cpdef num_permutations_without_fixed_pts(N):
    r"""
    Returns the number of permutations of N letters without fixed points.
    Uses a recursive formula.
    REFERENCE: www.math.umn.edu/~garrett/crypto/Overheads/06_perms_otp.pdf
    """
    cdef Integer res
    res = num_perm_dict.get(N)
    if not res:
        res=Integer(0)
        _num_permutations_without_fixed_pts(N,res.value)
        num_perm_dict[N]=res
    return res

cdef void _num_permutations_without_fixed_pts(long N,mpz_t res):
    cdef mpz_t tmp,b
    cdef int k
    if N==0:
        mpz_set_ui(res,1)
    else:
        mpz_init(tmp)
        mpz_init(b)
        factorial_mpz(N,res)
        for k from 1 <= k<=N:
            mpz_bin_uiui(b,N,k)
            _num_permutations_without_fixed_pts(N-k,tmp)
            mpz_mul(tmp,tmp,b)
            mpz_sub(res,res,tmp)
        mpz_clear(tmp)
        mpz_clear(b)






def sort_perms(a,b):
    if a<>b:
        sa=str(a)
        sb=str(b)
        if sa<sb:
            return -1
        else:
            return 1
    else:
        return 0



cpdef test_p_1(int N):
    cdef int j
    #cdef list fixpts
    fixpts=vector(ZZ,N)
    for j from 1<=j<=N:
        fixpts[j-1]=j




cpdef test_p_2(N):
    cdef int j
    fixpts=list()
    for j from 1<=j<=N:
        fixpts.append(j)


cpdef test_p_3(N):
    cdef int j
    fixpts=range(1,N+1)

cpdef test_p_4(int N):
    test_p_4_c(N)


cdef void test_p_4_c(int N):
    cdef int j
    cdef list fixpts
    fixpts=list()
    for j from 1<=j<=N:
        fixpts.append(j)

#cpdef print_vec(n,l):
#    _print_vec_(n,l)


cpdef test_pi(int n,list fix=[]):
    cdef int i
    cdef MyPermutationIterator PONEI
    for i from 0<= i <= 100:
        PONEI = MyPermutationIterator(n,fixed_pts=fix,verbose=3)
    return 0


cpdef test_p(int n):
    cdef int i
    cdef MyPermutation P
    for i from 0<= i <= 100:
        print ":::::::::::::::::::::::::::::::::::::::::::::::::::i=",i
        P = MyPermutation(length=n)
        print "P(",i,")=",P
    return 0

cpdef test_is_of_order(perml,o):
    cdef int* perm
    cdef int i,N,r
    N = len(perml)
    perm=<int*>sage_malloc(sizeof(int)*N)
    for i from 0 <= i <N:
        perm[i]=perml[i]
    r=_is_of_order(N,perm, o)
    sage_free(perm)
    return r
    

cpdef test_list(int N):
    cdef int i
    cdef list l
    l=[]
    for i from 0 <= i <=N:
        l.append(i)
    return l


cpdef transposition(int N,int i,int j):
    cdef MyPermutation p
    p=MyPermutation(length=N,init=1)
    p._entries[i-1]=j
    p._entries[j-1]=i
    p._list=[]
    p._str=''
    return p


#def _cmp_list_of_lists(l1,l2):
#    return len(l1) < len(l2)

## cimport cython
## from cython.operator cimport dereference as deref, preincrement as inc #dereference and increment operators
## from libcpp.vector cimport vector

## cdef extern from "<vector>" namespace "std":
##     cdef cppclass vector[T]:
##         #cppclass iterator:
##         #    T operator*()
##         #    iterator operator++()
##         #    bint operator==(iterator)
##         #    bint operator!=(iterator)
##         vector()
##         void push_back(T&)
##         T& operator[](int)
##         T& at(int)
##         #iterator begin()
##         #iterator end()


## cpdef test_list1(int N):
##     cdef int i
##     cdef vector[int] *v = new vector[int]()
##     for i from 0 <= i <=N:
##        v.push_back(i)
##     return v


cdef int gcd(int m, int n):
    cdef int tmp
    while m<>0:
        tmp = m
        m = n % m
        n = tmp       
    return n

cdef int lcm(int m, int n):
    return m / gcd(m, n) * n
