#################################################################################
#
# (c) Copyright 2010 William Stein
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
Bases of newforms for classical GL2 modular forms over QQ.
"""

def degrees(N, k, eps=None):
    """
    Return the degrees of the newforms of level N, weight k, with character eps.

    INPUT:

        - N -- level; positive integer or Dirichlet character
        - k -- weight; integer at least 2
        - eps -- None or Dirichlet character; if specified N is ignored
        
    EXAMPLES::

        sage: import psage
        sage: psage.modform.rational.degrees(11,2)
        [1]
        sage: psage.modform.rational.degrees(37,2)
        [1, 1]
        sage: psage.modform.rational.degrees(43,2)
        [1, 2]
        sage: psage.modform.rational.degrees(DirichletGroup(13).0^2,2)
        [1]
        sage: psage.modform.rational.degrees(13,2,DirichletGroup(13).0^2)
        [1]
        sage: psage.modform.rational.degrees(13,2)
        []
    """
    group = eps if eps else N
    from sage.all import ModularSymbols, dimension_new_cusp_forms
    d = dimension_new_cusp_forms(group, k)
    if d == 0:
        # A useful optimization!
        return []
    M = ModularSymbols(group=group, weight=k, sign=1).cuspidal_subspace()
    N = M.new_subspace()
    D = N.decomposition()
    # TODO: put in a consistency check.
    degs = [f.dimension() for f in D]
    assert sum(degs) == d, "major consistency check failed in degrees"
    return degs

def eigenvalue_fields(N, k, eps=None):
    """
    Return Hecke eigenvalue fields of the newforms in the space with
    given level, weight, and character.

    Note that computing this field involves taking random linear
    combinations of Hecke eigenvalues, so is not deterministic.  Set
    the random seed first (set_random_seed(0)) if you want this to be
    deterministic.

    INPUT:

        - N -- level; positive integer
        - k -- weight; integer at least 2
        - eps -- None or Dirichlet character; if specified N is ignored
        
    EXAMPLES::
    
        sage: import psage
        sage: psage.modform.rational.eigenvalue_fields(11,2)
        [Rational Field]
        sage: psage.modform.rational.eigenvalue_fields(43,2)
        [Rational Field, Number Field in alpha with defining polynomial x^2 - 2]
        sage: psage.modform.rational.eigenvalue_fields(DirichletGroup(13).0^2,2)
        [Cyclotomic Field of order 6 and degree 2]
        sage: psage.modform.rational.eigenvalue_fields(13,2,DirichletGroup(13).0^2)
        [Cyclotomic Field of order 6 and degree 2]
    """
    group = eps if eps else N
    from sage.all import ModularSymbols
    M = ModularSymbols(group=group, weight=k, sign=1).cuspidal_subspace()
    N = M.new_subspace()
    D = N.decomposition()
    X = [f.compact_system_of_eigenvalues([2])[1].base_ring() for f in D]
    return X






