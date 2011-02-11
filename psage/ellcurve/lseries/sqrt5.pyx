#################################################################################
#
# (c) Copyright 2011 William Stein
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
#################################################################################

"""
Fast Cython code needed to compute L-series of elliptic curves over F = Q(sqrt(5)).

USES:

   - The fast residue class rings defined in
     psage.modform.hilbert.sqrt5.sqrt5_fast for naive enumeration.

   - Drew Sutherlands smalljac for point counting

   - Lcalc for evaluating the L-series

   - Dokchitser as well.

   - Computes the *root number* in addition to the L-series.
   
"""

####################################################################
# Straightforward Elliptic Curve Pointcount
####################################################################

from psage.modform.hilbert.sqrt5.sqrt5_fast cimport (
    ResidueRingElement, ResidueRing_abstract, residue_element
    )

cpdef long ap_via_enumeration(ResidueRingElement a4, ResidueRingElement a6) except -1:
    """
    Compute the trace of Frobenius `a_p` on the elliptic curve defined
    by `y^2 = x^3 + a_4 x + a_6` using a straightforward enumeration
    algorithm.  Here `p` must be a prime of good reduction for the
    equation of the curve.

    EXAMPLES::

        sage: import psage.ellcurve.lseries.sqrt5 as s
        sage: K.<a> = NumberField(x^2-x-1)
        sage: E = EllipticCurve([1,-a,a,5*a-4,-3*a+4])
        sage: _,_,_,a4,a6 = E.short_weierstrass_model().a_invariants()
        sage: P = K.ideal(163)
        sage: import psage.modform.hilbert.sqrt5.sqrt5_fast as t
        sage: R = t.ResidueRing(P, 1)
        sage: s.ap_via_enumeration(R(a4), R(a6))
        120
        sage: k = P.residue_field(); E0 = E.change_ring(k); k.cardinality() + 1 - E0.cardinality()
        120
    """
    cdef long i, cnt = 1  # point at infinity
    cdef ResidueRing_abstract R = a4.parent()

    if R.e != 1:
        raise ValueError, "residue ring must be a field"
    
    cdef residue_element x, z, w
    for i in range(R.cardinality()):
        R.unsafe_ith_element(x, i)
        R.mul(z, x, x)   # z = x*x
        R.mul(z, z, x)   # z = x^3
        R.mul(w, a4.x, x)  # w = a4*x
        R.add(z, z, w)     # z = z + w = x^3 + a4*x
        R.add(z, z, a6.x)  # z = x^3 + a4*x + a6
        if R.element_is_0(z):
            cnt += 1
        elif R.is_square(z):
            cnt += 2
    return R.cardinality() + 1 - cnt


