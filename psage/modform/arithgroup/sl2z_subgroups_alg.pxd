from psage.groups.permutation_alg cimport MyPermutation 

cpdef are_mod1_equivalent(MyPermutation R1,MyPermutation S1, MyPermutation R2,MyPermutation S2,int verbose=?)

cdef int are_mod1_equivalent_c(int N,MyPermutation S1,MyPermutation R1,MyPermutation S2,MyPermutation R2,int* pres,int verbose=?,int check=?)
