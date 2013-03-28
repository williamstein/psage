include "sage/ext/interrupt.pxi" 
include "sage/ext/stdsage.pxi"  
include "sage/ext/cdefs.pxi"
include "sage/rings/mpc.pxi"
#include "../ext/gmp.pxi"
include "sage/ext/python_int.pxi"
from sage.structure.sage_object cimport SageObject
from sage.rings.integer cimport Integer

cdef class MyPermutation(SageObject):
    cdef int *_entries
    cdef int _N
    cdef long _hash
    cdef void c_new(self,int* entries)
    cdef list _list
    cdef str _str
    cdef int _init
    cdef int _order
    cdef int _verbose
    cdef int _rep  # determines how the permutation is printed
    cpdef pow(self,int k)
    #cpdef mult(self,MyPermutation other)
    #cdef _mult(self,MyPermutation other,MyPermutation res)
    #cdef _mult_list(self,list other,MyPermutation res)
    cdef void set_entries(self,int *entries)
    cdef int get_unsafe(self,int i)
    cpdef int _get_unsafe(self,int i)
    cpdef conjugate(self,other)
    cpdef is_order_eq(self,int o)
    cpdef _get_order(self)
    cpdef _inverse(self)
    cpdef str export_as_string(self,str sep=*)
    cdef void _conj_w_transp(self,int a, int b)
    cpdef conjugate_with_transposition(self,int a, int b)
    cdef void _conjugate_ptr(self,int *entries,int *other)
    cdef MyPermutation _conjugate(self,MyPermutation other)
    cdef MyPermutation _conjugate_list(self,list other)
    cdef MyPermutation _conjugate_vec(self,int * other)
    cdef _pickle_version0(self)
    cdef _unpickle_version0(self, data)
    cpdef fixed_elements(self)
    cdef int num_fixed_c(self)
    cpdef num_fixed(self)
    cpdef _mult_perm(self,MyPermutation other)
            
    cdef int eq_c(self,MyPermutation other)
    cpdef int eq(self,MyPermutation other)
    cdef int N_c(self)
    cpdef int N(self)
     
cdef class MyPermutationIterator(SageObject):
    cdef int _N
    cdef long _num
    cdef Integer _max_num
    cdef long _cur
    cdef int _verbose
    cdef int _current_piv
    cdef int * _current_state_c
    cdef int * _current_state_o
    cdef int * _current_perm 
    cdef int * _fixed_pts   # -- make into a python object instead -- due to malloc problems...
    cdef int ** _list_of_perms
    cpdef int _order
    cdef int _is_subiterator
    cdef int _got_list
    cdef int _num_fixed
    cdef MyPermutationIterator _fixed_pt_free_iterator
    cdef list _fixed_pt_free_iterator_labels
    cpdef list(self)
    cpdef list_new(self)
    cpdef _get_list_of_perms(self)
    cpdef _next(self)
    cpdef _goto_first(self)
    cpdef _next_free(self)
    cpdef long num(self)
    cpdef _rewind(self)
    cdef set_current_from_ffree(self)
    cdef _get_next_permutation_recursive(self,int tmp, int start,long num,int *lista)
    cpdef get_perms_test(self)



cdef void _mult_perm(int N, int* res,int *left,int* right)
    
cdef print_vec(int n,int *list)
cdef _are_eq_vec(int n,int *a,int *b)
cdef void _conjugate_perm(int N,int* res,int *a,int* b)
cpdef transposition(int N,int i,int j)
cdef int are_transitive_perm_c(int *Rl,int *Sl, int *gotten, int n,int verbose=?)
cdef list perm_to_cycle_c(int N,int *perm)
cdef print_vec(int n,int *list)
