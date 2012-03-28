r"""
Functions to create and use basis of graded modules which are products of
Maass lifts.

AUTHORS:

- Martin Raum (2009 - 08 - 03) Initial version
"""

#===============================================================================
# 
# Copyright (C) 2009 Martin Raum
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

from sage.matrix.constructor import matrix
from sage.misc.misc import prod
from sage.rings.arith import random_prime
from sage.rings.integer_ring import ZZ
from sage.rings.padics.factory import Qp
from sage.rings.rational_field import QQ
from psage.modform.paramodularforms.siegelmodularformg2_fegenerators import SiegelModularFormG2MaassLift
from psage.modform.fourier_expansion_framework.monoidpowerseries.monoidpowerseries_lazyelement import EquivariantMonoidPowerSeries_LazyMultiplication
from psage.modform.paramodularforms.siegelmodularformg2_submodule import SiegelModularFormG2Submodule_maassspace
from sage.interfaces.magma import magma

def spanning_maass_products(ring, weight, subspaces = None, lazy_rank_check = False) :
    r"""This function is intended to quickly find good generators. This may be
    calculated in low precision and then ported to higher.
    
    the maass spaces of all subspaces are expected to be ordered the same way as
    ``ring.type()._maass_generator_preimages``.
    """
    ## Outline :
    ##  - get the basis as polynomials form the ring
    ##  - get the maass products, convert them to polynomials and try to 
    ##  construct a basis
    ##  - get the maass space associated to this weight and get its coordinates
    ##  - return a basis containing of maass forms and the products

    if subspaces is None :
        subspaces = dict()
        
    assert ring.base_ring() is QQ or ring.base_ring() is ZZ, \
           "%s doesn't have rational base field" % ring
    
    dim = ring.graded_submodule(weight).dimension()
    precision = ring.fourier_expansion_precision()
    maass_lift_func = ring.type()._maass_generators
    maass_lift_preimage_func = ring.type()._maass_generator_preimages
                
    lift_polys = []; preims = []
    monomials = set()
    
    for k in xrange(min((weight // 2) - ((weight // 2) % 2), weight - 4), 3, -2) :
        if k not in subspaces :
            subspaces[k] = ring.graded_submodule(k)
        if weight - k not in subspaces :
            subspaces[weight - k] = ring.graded_submodule(weight - k)
        
        try :
            kgs = zip(subspaces[k].maass_space()._provided_basis(), maass_lift_preimage_func(k))
            ogs = zip(subspaces[weight - k].maass_space()._provided_basis(), maass_lift_preimage_func(weight - k))
        except ValueError :
            continue
        
        for f_coords,f_pre in kgs :
            for g_coords,g_pre in ogs :
                assert f_coords.parent().graded_ambient() is ring
                assert g_coords.parent().graded_ambient() is ring
 
                f = ring(f_coords)
                g = ring(g_coords)
                    
                fg_poly = (f*g)._reduce_polynomial()
                monomials = monomials.union(set(fg_poly.monomials()))
                
                M = matrix( QQ, len(lift_polys) + 1,
                            [ poly.monomial_coefficient(m)
                              for poly in lift_polys + [fg_poly]
                              for m in monomials] )
                
                # TODO : Use linbox
                try :
                    if magma(M).Rank() <= len(lift_polys) :
                        continue
                except TypeError :
                    for i in xrange(10) :
                        Mp = matrix(Qp(random_prime(10**10), 10), M)
                        if Mp.rank() > len(lift_polys) : break
                    else :
                        if lazy_rank_check :
                            continue
                        elif M.rank() <= len(lift_polys) :
                            continue
                
                lift_polys.append(fg_poly)
                preims.append([f_pre, g_pre])
                        
                if len(lift_polys) == dim :
                    break
                
    if len(lift_polys) == dim :
        ss = construct_from_maass_products(ring, weight, zip(preims, lift_polys), True, False, True)
        poly_coords = [[0]*i + [1] + (len(lift_polys) - i - 1)*[0] for i in xrange(len(lift_polys))]
    else :
        # The products of two Maass lifts don't span this space
        # we hence have to consider the Maass lifts
        ss = ring.graded_submodule(weight)
        poly_coords = [ss.coordinates(ring(p)) for p in lift_polys]

    maass_lifts = maass_lift_func(weight, precision)
    preims[0:0] = [[p] for p in maass_lift_preimage_func(weight)]
    all_coords = map(ss.coordinates, maass_lifts) + poly_coords
    
    M = matrix(QQ, all_coords).transpose()

    nmb_mls = len(maass_lifts)
    for i in xrange(10) :
        Mp = matrix(Qp(random_prime(10**10), 10), M)
        pvs = Mp.pivots()
        if pvs[:nmb_mls] == range(nmb_mls) and \
           len(pvs) == dim :
            break
    else :
        pvs = M.pivots()
        

    if len(pvs) == dim :
        return [(preims[i], ring(ss(all_coords[i])).polynomial()) for i in pvs]
    else :        
        raise RuntimeError, "The products of at most two Maass lifts don't span " + \
                            "this space"

def construct_from_maass_products(ring, weight, products, is_basis = True,
                                  provides_maass_spezialschar = False,
                                  is_integral = False,
                                  lazy_rank_check = True) :
    r"""
    Pass the return value of spanning_maass_products of a space of Siegel modular
    forms of same type. This will return a space using these forms.
    Whenever ``is_basis`` is False the the products are first filtered to yield a
    basis.
    """
    assert QQ.has_coerce_map_from(ring.base_ring()) or ring.base_ring() is ZZ, \
           "%s doesn't have rational base field" % ring
           
    dim = ring.graded_submodule(weight).dimension()

    ## if the products don't provide a basis, we have to choose one
    if not is_basis :
        ## we prefer to use Maass lifts, since they are very cheap
        maass_forms = []; non_maass_forms = []
        for p in products :
            if len(p[0]) == 1 :
                maass_forms.append(p)
            else :
                non_maass_forms.append(p)
        
        
        monomials = set()
        lift_polys = []
        products = []
        for lifts, lift_poly in maass_forms + non_maass_forms :
            red_poly = ring(ring.relations().ring()(lift_poly))._reduce_polynomial()
            monomials = monomials.union(set(red_poly.monomials()))
                
            M = matrix( QQ, len(lift_polys) + 1,
                        [ poly.monomial_coefficient(m)
                          for poly in lift_polys + [lift_poly]
                          for m in monomials] )                
                                    
            # TODO : Use linbox
            try :
                if magma(M).Rank() > len(lift_polys) :
                    break
            except TypeError :
                for i in xrange(10) :
                    Mp = matrix(Qp(random_prime(10**10), 10), M)
                    if Mp.rank() > len(lift_polys) : break
                else :
                    if lazy_rank_check :
                        continue
                    elif M.rank() <= len(lift_polys) :
                        continue
                    
            lift_polys.append(red_poly)
            products.append((lifts, red_poly))

            if len(products) == dim :
                break
        else :
            raise ValueError, "products don't provide a basis"
    
    basis = []
     
    if provides_maass_spezialschar :
        maass_form_indices = []
    for i, (lifts, lift_poly) in enumerate(products) :
        e = ring(ring.relations().ring()(lift_poly))

        if len(lifts) == 1 :
            l = lifts[0]
            e._set_fourier_expansion(
                SiegelModularFormG2MaassLift(l[0], l[1], ring.fourier_expansion_precision(),
                                             is_integral = is_integral) )
        elif len(lifts) == 2 :
            (l0, l1) = tuple(lifts)
            e._set_fourier_expansion(
                EquivariantMonoidPowerSeries_LazyMultiplication(
                    SiegelModularFormG2MaassLift(l0[0], l0[1], ring.fourier_expansion_precision(),
                                                 is_integral = is_integral),
                    SiegelModularFormG2MaassLift(l1[0], l1[1], ring.fourier_expansion_precision(),
                                                 is_integral = is_integral) ) ) 

        else :
            e._set_fourier_expansion( 
                prod( SiegelModularFormG2MaassLift(l[0], l[1], ring.precision(),
                                                   is_integral = is_integral)
                      for l in lifts) )            
        basis.append(e)
        
        if provides_maass_spezialschar and len(lifts) == 1 :
            maass_form_indices.append(i)
    
    ss = ring._submodule(basis, grading_indices = (weight,), is_heckeinvariant = True)
    if provides_maass_spezialschar : 
        maass_coords = [ ss([0]*i + [1] + [0]*(dim-i-1))
                         for i in maass_form_indices ]
        ss.maass_space.set_cache( 
          SiegelModularFormG2Submodule_maassspace(ss, map(ss, maass_coords)) )
                                 
    return ss
