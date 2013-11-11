include "sage/ext/interrupt.pxi" 
include "sage/ext/stdsage.pxi"  
include "sage/ext/cdefs.pxi"
include "sage/rings/mpc.pxi"
#include "../ext/gmp.pxi"
#include "sage/ext/python_int.pxi"
from sage.structure.sage_object cimport SageObject
from sage.structure.parent cimport Parent
from sage.rings.integer cimport Integer

cdef class MyPermutation(SageObject):
    cdef int *_entries
    cdef int _N
    cdef long _hash
    cdef void c_new(self,int* entries)
    cdef list _list
    cdef list _cycles_list
    cdef list _cycles_permutations
    cdef int _cycles_list_is_ordered 
    cdef list _cycle_type
    cdef list _cycles_ordered_as_list
    cdef int* _cycles # matrix containing all cycles and their lenghts
    cdef int* _cycle_lens
    cdef int _num_cycles
    cdef str _str
    cdef int _init
    cdef int _order
    cdef int _verbose
    cdef int _rep  # determines how the permutation is printed
    cpdef pow(self,int k)
    #cpdef mult(self,MyPermutation other)
    #cdef _mult(self,MyPermutation other,MyPermutation res)
    #cdef _mult_list(self,list other,MyPermutation res)
    cdef int set_entries(self,int *entries)
    cdef int get_unsafe(self,int i)
    cpdef int _get_unsafe(self,int i)
    cpdef conjugate(self,other)
#    Routines dealing with the cycles
    cpdef cycles_as_lists(self)
    cpdef cycles_as_permutations(self)
    cpdef cycles_ordered_as_list(self)
    cpdef set_cycle_info(self,int reset=*)
    cdef int set_cycle_info_c(self)
    cpdef  list _list_from_cycles(self,list cycles)
    cpdef num_cycles(self)

    
    cpdef is_order_eq(self,int o)
    cdef int is_order_eq_c(self,int o)
    cpdef _get_order(self)
    cpdef _inverse(self)
    cpdef _copy_c(self)
    cpdef str export_as_string(self,str sep=*)
    cpdef list list_from_string(self,s)
    cdef void _conj_w_transp(self,int a, int b,int verbose=*)
    cpdef conjugate_with_transposition(self,int a, int b,int verbose=*)
    cpdef is_conjugate_to(self,MyPermutation other,int ret_perm=*)
    cdef void _conjugate_ptr(self,int *entries,int *other)
    cdef MyPermutation _conjugate(self,MyPermutation other)
    cdef MyPermutation _conjugate_list(self,list other)
    cdef MyPermutation _conjugate_vec(self,int * other)
    cdef _pickle_version0(self)
    cdef _unpickle_version0(self, data)
    cpdef fixed_elements(self)
    cpdef non_fixed_elements(self)
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
    cdef Integer _cur
    cdef int _verbose
    cdef int _current_piv
    cdef int _map_from
    cdef int _map_at_most_to
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
    cdef int _has_fixed_pt_free_iterator
    cdef int * _fixed_pt_free_iterator_labels
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

cdef class CycleCombinationIterator(Parent):
    cdef list _cycles
    cdef int* _cycle_types
    cdef int _N
    cdef int _length
    cdef int _num_cycles
    cdef int _dbase
    cdef list _cache
    cdef int _got
    cpdef list(self)
    #    cdef dict __cached_methods
    cpdef length(self)
    cpdef  MyPermutation permutation_nr(self,int M)
    cdef MyPermutation permutation_nr_c(self,int M)
    cdef int permutation_nr_c_ptr(self,int M,int* res)    
#    cdef int permutation_nr_c(self,int M,MyPermutation p)

cdef void _mult_perm_unsafe(int N, int* res,int *left,int* right)
cpdef are_conjugate_perm(MyPermutation A,MyPermutation B)
cpdef get_conjugating_perm_list(list Al,list Bl)
cdef print_vec(int n,int *list)
cdef _are_eq_vec(int n,int *a,int *b)
cdef void _conjugate_perm(int N,int* res,int *a,int* b)
cpdef transposition(int N,int i,int j)
cdef int are_transitive_perm_c(int *Rl,int *Sl, int *gotten, int n,int verbose=?)
cdef int perm_to_cycle_c(int N,int *perm,int* cycle,int* cycle_lens)
cdef print_vec(int n,int *list)
cdef int num_cycles_c(int N,int *perm)
cdef MyPermutation  get_conjugating_perm_ptr_unsafe(int mu, int* Al,int* Bl)
