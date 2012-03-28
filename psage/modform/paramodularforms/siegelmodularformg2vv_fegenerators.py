r"""
Generator functions for the Fourier expansion of vector valued Siegel modular forms. 

AUTHORS:

- Martin Raum (2010 - 05 - 12) Initial version
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

from psage.modform.paramodularforms.siegelmodularformg2_fourierexpansion import SiegelModularFormG2VVFourierExpansionRing
from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_ambient import EquivariantMonoidPowerSeriesAmbient_abstract
from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_lazyelement import EquivariantMonoidPowerSeries_lazy
from psage.modform.paramodularforms.siegelmodularformg2vv_fegenerators_cython import satoh_dz

#===============================================================================
# SiegelModularFormG2SatohBracket
#===============================================================================

def SiegelModularFormG2SatohBracket(f, g, f_weight = None, g_weight = None) :
    r"""
    INPUT:
    
    - `f` -- Fourier expansion of a classical Siegel modular form.
    
    - `g` -- Fourier expansion of a classical Siegel modular form.
    
    - ``f_weight`` -- the weight of `f` (default: ``None``).
     
    - ``g_weight`` -- the weight of `g` (default: ``None``).
    
    OUTPUT:
    
    - The Fourier expansion of a vector-valued Siegel modular form.
    """        
    if f_weight is None :
        f_weight = f.weight()
    if g_weight is None :
        g_weight = g.weight()
    
    if not isinstance(f.parent(), EquivariantMonoidPowerSeriesAmbient_abstract) :
        f = f.fourier_expansion()
    if not isinstance(g.parent(), EquivariantMonoidPowerSeriesAmbient_abstract) :
        g = g.fourier_expansion()
        
    if f.parent() != g.parent() :
        if f.parent().has_coerce_map_from(g.parent()) :
            g = f.parent()(g)
        elif g.parent().has_coerce_map_from(f.parent()) :
            f = g.parent()(f)
        else :
            from sage.categories.pushout import pushout
            parent = pushout(f.parent(), g.parent())
            f = parent(f)
            g = parent(g)
            
    precision = min(f.precision(), g.precision())
    
    expansion_ring = SiegelModularFormG2VVFourierExpansionRing(
                                            f.parent().coefficient_domain().fraction_field() )
    coefficients_factory = DelayedFactory_SMFG2_satohbracket(
                             f, g, f_weight, g_weight, expansion_ring.coefficient_domain() )

    return EquivariantMonoidPowerSeries_lazy( expansion_ring, precision,
                                              coefficients_factory.getcoeff )

#===============================================================================
# DelayedFactory_SMFG2_satohbracket
#===============================================================================

class DelayedFactory_SMFG2_satohbracket :
    def __init__( self, f, g, f_weight, g_weight, coefficient_domain ) :
        self.__f = f
        self.__g = g
        self.__f_weight = f_weight
        self.__g_weight = g_weight
        self.__coefficient_domain = coefficient_domain
        
        self.__series = None
    
    def getcoeff( self, key, **kwds ) :
        (_, k) = key
        # for speed we ignore the character 
        if self.__series is None :
            self.__series = \
             _satoh_bracket( self.__f, self.__g,
                             self.__f_weight, self.__g_weight )
        
        return self.__coefficient_domain( self.__series[k] )

def _satoh_bracket(f, g, f_weight, g_weight) :
    r"""
    We suppose that `f` and `g` are either contained in a ring of a scalar
    valued Siegel modular forms or that they are equivariant monoid
    power series.
    `f` and `g` must have the same parent after conversion to fourier
    expansions.
    
    OUTPUT:
    
    - The return value is a Equivariant monoid power series.
    """
    if not f.parent() == g.parent() :
        raise ValueError, "The fourier expansions of f and g must" + \
                          " have the same parent"

    base_ring = f.parent().coefficient_domain()
    if len(f.non_zero_components()) != 1 :
        raise ValueError, "f must have only one non-vanishing character"
    if len(g.non_zero_components()) != 1 :
        raise ValueError, "g must have only one non-vanishing character"
    
    fe_ring = SiegelModularFormG2VVFourierExpansionRing(base_ring)
    R = fe_ring.coefficient_domain()
            
    dzf = fe_ring(satoh_dz(f.coefficients(False), R))
    dzg = fe_ring(satoh_dz(g.coefficients(False), R))
    
    return (g * dzf) / f_weight - (f * dzg) / g_weight
