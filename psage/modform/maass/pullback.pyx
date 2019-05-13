# -*- coding=utf-8 -*-
#*****************************************************************************
#  Copyright (C) 2011 Fredrik Strömberg <stroemberg@mathematik.tu-darmstadt.de>
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
#*****************************************************************************
r"""
Various algorithms for pulling back points to a fundamental domain.

"""
include 'sage/ext/stdsage.pxi'
from sage.libs.gmp.all cimport *
from cysignals.signals cimport *
include "sage/ext/gmp.pxi"
include "sage/rings/mpc.pxi"

from sage.libs.mpfr cimport *

from psage.modforms.mysubgroups_alg cimport *

### Toplevel functions.


def pullback_to_G(G,x,y,mp_ctx=None):
    r""" Mpmath version of general pullback alg

    INPUT:
    
     - ``G`` -- MySubgroup 
     - ``x`` -- real
     - ``y`` -- real > 0
     - ``mp_ctx`` -- mpmath context to use for pullback (default None uses cython double)

    EXAMPLES::

    
        sage: G=MySubgroup(SL2Z)
        sage: [x,y,A]=pullback_to_G_mp(G,mpmath.mpf(0.1),mpmath.mpf(0.1))
        sage: x,y
        (mpf('6.9388939039072274e-16'), mpf('4.9999999999999991'))
        sage: A
        [ 5 -1]
        [ 1  0]
        sage: G=MySubgroup(Gamma0(3))
        sage: [x,y,A]=pullback_to_G_mp(G,mpmath.mpf(0.1),mpmath.mpf(0.1))
        sage: x,y
        (mpf('-0.038461538461538491'), mpf('0.19230769230769232'))
        sage: A
        [-1  0]
        [ 6 -1]

    

    """
    import mpmath
    try:
        reps=G._coset_reps
    except AttributeError:
        reps=list(G.coset_reps())

    if(mp_ctx==mpmath.mp):
        try:
            x.ae; y.ae
        except AttributeError:
            raise Exception,"Need x,y in mpmath format!"
        [x1,y1,[a,b,c,d]]=pullback_to_psl2z_mp(x,y)
    else:
        [x1,y1,[a,b,c,d]]=pullback_to_psl2z_dble1(x,y)
    A=SL2Z([a,b,c,d])
    try:
        for V in reps:
            B=V*A
            if(B in G):
                raise StopIteration
    except StopIteration:            
        pass
    else:
        raise Exception,"Did not find coset rep. for A=%s" % A
    [xpb,ypb]=apply_sl2_map(x,y,B,mp_ctx=mp_ctx)
    return [xpb,ypb,B]
    


