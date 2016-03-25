r"""
The diffential operator involved in the definition of the Satoh bracket.

AUTHORS:

- Martin Raum (2010 - 03 - 16) Initial version.
"""

#===============================================================================
# 
# Copyright (C) 2010 Martin Raum
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================

include "cysignals/signals.pxi"
include "stdsage.pxi"
include "cdefs.pxi"

cpdef satoh_dz(f, R) :
    ## Satoh's definition of the operation on Siegel modular forms yields
    ## a polynomial substitution x -> dx - cy, y -> -b + ay.
    ## We will need x -> ax + by, y -> cx + dy.
    ## We hence calculate the operator matrix(2, [d z_2, -d z_3 / 2, *, d z_1]) 
    
    x = R.gen(0)
    y = R.gen(1)
    xsq = x*x
    xy  = x*y
    ysq = y*y
        
    res = PY_NEW(dict)
    res2 = PY_NEW(dict)
    res3 = PY_NEW(dict)
    for k in f :
        res[k] = xsq * (k[2] * f[k]) + xy * (k[1] * f[k]) + ysq * (k[0] * f[k])

    return res
