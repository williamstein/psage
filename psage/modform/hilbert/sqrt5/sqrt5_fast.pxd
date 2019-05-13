from sage.rings.ring cimport CommutativeRing
from sage.rings.number_field.number_field_element_quadratic cimport NumberFieldElement_quadratic

# Residue element
ctypedef long residue_element[2]

cdef class ResidueRingElement

cdef class ResidueRing_abstract(CommutativeRing):
    cdef object P, F
    cdef public object element_class, _residue_field
    cdef long n0, n1, p, e, _cardinality
    cdef long im_gen0
    cdef object element_to_residue_field(self, residue_element x)
    cdef void add(self, residue_element rop, residue_element op0, residue_element op1)
    cdef void sub(self, residue_element rop, residue_element op0, residue_element op1)
    cdef void mul(self, residue_element rop, residue_element op0, residue_element op1)
    cdef int inv(self, residue_element rop, residue_element op) except -1
    cdef bint is_unit(self, residue_element op)
    cdef void neg(self, residue_element rop, residue_element op)
    cdef ResidueRingElement new_element(self)
    cdef int coefficients(self, long* v0, long* v1, NumberFieldElement_quadratic x) except -1
    cdef int coerce_from_nf(self, residue_element r, NumberFieldElement_quadratic x) except -1
    cdef bint element_is_1(self, residue_element op)
    cdef bint element_is_0(self, residue_element op)
    cdef void set_element_to_1(self, residue_element op)
    cdef void set_element_to_0(self, residue_element op)
    cdef void set_element(self, residue_element rop, residue_element op)
    cdef int set_element_from_tuple(self, residue_element rop, x) except -1
    cdef int cmp_element(self, residue_element left, residue_element right)
    cdef int pow(self, residue_element rop, residue_element op, long e) except -1
    cdef bint is_square(self, residue_element op) except -2
    cdef int sqrt(self, residue_element rop, residue_element op) except -1
    cdef int ith_element(self, residue_element rop, long i) except -1
    cpdef long cardinality(self)
    cdef void unsafe_ith_element(self, residue_element rop, long i)
    cdef int next_element(self, residue_element rop, residue_element op) except -1
    cdef bint is_last_element(self, residue_element op)
    cdef long index_of_element(self, residue_element op) except -1
    cdef long index_of_element_in_P(self, residue_element op) except -1
    cdef int next_element_in_P(self, residue_element rop, residue_element op) except -1
    cdef bint is_last_element_in_P(self, residue_element op)
    cdef element_to_str(self, residue_element op)


cdef class ResidueRingElement:
    cdef residue_element x
    cdef ResidueRing_abstract _parent
    cpdef long index(self)
    cpdef parent(self)
    cdef new_element(self)
    cpdef bint is_unit(self)
    cpdef bint is_square(self)
    cpdef sqrt(self)




