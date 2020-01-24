r"""
Gradings of a polynomial quotient ring.

AUTHOR :
    -- Martin Raum (2009 - 07 - 27) Initial version
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

from past.builtins import cmp
from builtins import str
from builtins import map
from builtins import range
from operator import mul
from operator import xor
from sage.misc.latex import latex
from sage.modules.free_module import FreeModule
from sage.rings.infinity import infinity
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.structure.all import Sequence
from sage.structure.sage_object import SageObject
from functools import reduce

class Grading_abstract ( SageObject ) :
    r"""
    A common interface for monomial gradings of polynomial rings
    `R[x_1, .., x_n]`.
    """
    
    def ngens(self) :
        r"""
        The number of generators of the polynomial ring.
        
        OUTPUT:
            An integer.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import *
            sage: Grading_abstract().ngens()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError
        
    def gen(self, i) :
        r"""
        The grading associated to the `i`-th generator of the polynomial ring.
        
        INPUT:
            - `i` -- An integer.
        
        OUTPUT:
            A grading index.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import *
            sage: Grading_abstract().gen(1)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError
        
    def gens(self) :
        r"""
        The gradings of the generators of the polynomial ring.
        
        OUTPUT:
            A tuple of grading indices.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import *
            sage: Grading_abstract().gens()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError
        
    def index(self, x) :
        r"""
        The grading value of `x` with respect to this grading. 

        INPUT:
            - `x` -- A tuple of length `n`.

        OUTPUT:
            A grading index.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import *
            sage: Grading_abstract().index((1))
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError
        
    def basis(self, index, vars = None) :
        r"""
        All monomials that are have given grading index involving only the
        given variables.
        
        INPUT:
            - ``index`` -- A grading index.
            - ``vars``  -- A list of integers from `0` to `n - 1` or ``None``
                           (default: ``None``); If ``None`` all variables
                           will be considered.
        
        OUTPUT:
            A list of tuples of integers, each of which has length `n`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import *
            sage: Grading_abstract().basis(1)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError
    
    def subgrading(self, gens) :
        r"""
        The grading of same type for the ring with only the variables given by
        ``gens``.
        
        INPUT:
            - ``gens`` - A list of integers.
        
        OUTPUT:
            An instance of :class:~`.Grading_abstract`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import *
            sage: Grading_abstract().subgrading([1,2])
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError
    
    def __contains__(self, x) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import *
            sage: 1 in Grading_abstract()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError
    
    def _repr_(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import *
            sage: Grading_abstract()
            A grading
        """
        return "A grading"
    
    def _latex_(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import *
            sage: latex( Grading_abstract() )
            \text{A grading}
        """
        return r"\text{A grading}"
                
#===============================================================================
# DegreeGrading
#===============================================================================

class DegreeGrading( Grading_abstract ) :
    r"""
    This class implements a monomial grading for a polynomial ring
    `R[x_1, .., x_n]`.
    """
  
    def __init__( self, degrees ) :
        r"""
        INPUT:
            - ``degrees`` -- A list or tuple of `n` positive integers.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import *
            sage: g = DegreeGrading((1,2))
            sage: g = DegreeGrading([5,2,3])
        """
        
        self.__degrees = tuple(degrees)
        self.__module = FreeModule(ZZ, len(degrees))
        self.__module_basis = self.__module.basis()
        
    def ngens(self) :
        r"""
        The number of generators of the polynomial ring.
        
        OUTPUT:
            An integer.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import *
            sage: DegreeGrading((1,2)).ngens()
            2
            sage: DegreeGrading([5,2,3]).ngens()
            3
        """
        return len(self.__degrees)
     
    def gen(self, i) :
        r"""
        The number of generators of the polynomial ring.
        
        OUTPUT:
            An integer.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import *
            sage: DegreeGrading([5,2,3]).gen(0)
            5
            sage: DegreeGrading([5,2,3]).gen(3)
            Traceback (most recent call last):
            ...
            ValueError: Generator 3 does not exist.
        """
        if i < len(self.__degrees) :
            return self.__degrees[i]
        
        raise ValueError("Generator %s does not exist." % (i,))
    
    def gens(self) :
        r"""
        The gradings of the generators of the polynomial ring.
        
        OUTPUT:
            A tuple of integers.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import *
            sage: DegreeGrading((1,2)).gens()
            (1, 2)
            sage: DegreeGrading([5,2,3]).gens()
            (5, 2, 3)
        """
        return self.__degrees
                
    def index(self, x) :
        r"""
        The grading value of `x` with respect to this grading. 

        INPUT:
            - `x` -- A tuple of length `n`.

        OUTPUT:
            An integer.
                
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import *
            sage: g = DegreeGrading([5,2,3])
            sage: g.index((2,4,5))
            33
            sage: g.index((2,3))
            Traceback (most recent call last):
            ...
            ValueError: Tuple must have length 3.
        """
        ## TODO: We shouldn't need this.
        #if len(x) == 0 : return infinity
        
        if len(x) != len(self.__degrees) :
            raise ValueError( "Tuple must have length %s." % (len(self.__degrees),))
        
        return sum( map(mul, x, self.__degrees) )
        
    def basis(self, index, vars = None) :
        r"""
        All monomials that are have given grading index involving only the
        given variables.
        
        INPUT:
            - ``index`` -- A grading index.
            - ``vars``  -- A list of integers from `0` to `n - 1` or ``None``
                           (default: ``None``); If ``None`` all variables
                           will be considered.
        
        OUTPUT:
            A list of elements in `\Z^n`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import *
            sage: g = DegreeGrading([5,2,3])
            sage: g.basis(2)
            [(0, 1, 0)]
            sage: g.basis(10, vars = [0])
            [(2, 0, 0)]
            sage: g.basis(20, vars = [1,2])
            [(0, 10, 0), (0, 7, 2), (0, 4, 4), (0, 1, 6)]
            sage: g.basis(-1)
            []
            sage: g.basis(0)
            [(0, 0, 0)]
        """
        if index < 0 :
            return []
        elif index == 0 :
            return [self.__module.zero_vector()]
        
        if vars == None :
            vars = [m for (m,d) in enumerate(self.__degrees) if d <= index]
        if len(vars) == 0 :
            return []
            
        res = [ self.__module_basis[vars[0]] + m
                for m in self.basis(index - self.__degrees[vars[0]], vars) ]
        res += self.basis(index, vars[1:])
        
        return res

    def subgrading(self, gens) :
        r"""
        The grading of same type for the ring with only the variables given by
        ``gens``.
        
        INPUT:
            - ``gens`` - A list of integers.
        
        OUTPUT:
            An instance of :class:~`.DegreeGrading`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import *
            sage: g = DegreeGrading([5,2,3])
            sage: g.subgrading([2])
            Degree grading (3,)
            sage: g.subgrading([])
            Degree grading ()
        """
        return DegreeGrading([self.__degrees[i] for i in gens])
        
    def __contains__(self, x) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import *
            sage: g = DegreeGrading([5,2,3])
            sage: 5 in g
            True
            sage: "t" in g
            False
        """
        return isinstance(x, (int, Integer)) 
    
    def __cmp__(self, other) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import *
            sage: g = DegreeGrading([5,2,3])
            sage: g == DegreeGrading([5,2,3])
            True
            sage: g == DegreeGrading((2,4))
            False
        """
        c = cmp(type(self), type(other))
        
        if c == 0 :
            c = cmp(self.__degrees, other.__degrees)

        return c
    
    def __hash__(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import *
            sage: hash( DegreeGrading([5,2,3]) )
            -1302562269         # 32-bit
            7573306103633312291 # 64-bit
        """
        return hash(self.__degrees)
    
    def _repr_(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import *
            sage: DegreeGrading([5,2,3])
            Degree grading (5, 2, 3)
        """
        return "Degree grading %s" % str(self.__degrees)
    
    def _latex_(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import *
            sage: latex( DegreeGrading([5,2,3]) )
            \text{Degree grading } \left(5, 2, 3\right)
        """
        return r"\text{Degree grading }" + latex(self.__degrees)

class TrivialGrading ( Grading_abstract ) :
    r"""
    A grading for a polynomial ring `R[x_1, .., x_n]` assigning to every
    element the same, arbitrary index.
    """
    
    def __init__(self, nmb_generators, index) :
        r"""
        INPUT:
            - ``nmb_generators`` -- A positive integer; The number `n` of
                                    variables of the graded polynomial ring.
            - ``index``          -- An arbitrary object; The index assigned to
                                    all monomials.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import *
            sage: g = TrivialGrading( 3, "t" )
            sage: g = TrivialGrading( 3, None )
        """
        self.__ngens = nmb_generators
        self.__index = index
        
    def ngens(self) :
        r"""
        The number of generators of the polynomial ring.
        
        OUTPUT:
            An integer.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import *
            sage: TrivialGrading( 3, "t" ).ngens()
            3
        """
        return self.__ngens
    
    def gen(self, i) :
        r"""
        The number of generators of the polynomial ring.
        
        OUTPUT:
            An arbitrary object.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import *
            sage: TrivialGrading( 3, "t" ).gen(2)
            't'
            sage: TrivialGrading( 3, "t" ).gen(3)
            Traceback (most recent call last):
            ...
            ValueError: Generator 3 does not exist.
        """
        if i < self.__ngens :
            return self.__index
        
        raise ValueError("Generator %s does not exist." % (i,))
    
    def gens(self) :
        r"""
        The gradings of the generators of the polynomial ring.
        
        OUTPUT:
            A tuple of objects of arbitrary but equal type.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import *
            sage: TrivialGrading( 3, "t" ).gens()
            ('t', 't', 't')
        """
        return tuple(self.__ngens*[self.__index])
    
    def index(self, x) :
        r"""
        The grading value of `x` with respect to this grading. 

        INPUT:
            - `x` -- A tuple of length `n`.

        OUTPUT:
            A grading index.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import *
            sage: TrivialGrading( 3, "t" ).index((1,0,1))
            't'
            sage: TrivialGrading( 3, "t" ).index((1,0))
            Traceback (most recent call last):
            ...
            ValueError: Tuple must have length 3.
        """
        if len(x) != self.__ngens :
            raise ValueError( "Tuple must have length %s." % (self.__ngens,))

        return self.__index
    
    def basis(self, index, vars = None) :
        r"""
        All degree one monomials that are have given grading index involving only the
        given variables.
        
        INPUT:
            - ``index`` -- A grading index.
            - ``vars``  -- A list of integers from `0` to `n - 1` or ``None``
                           (default: ``None``); If ``None`` all variables
                           will be considered.
        
        OUTPUT:
            A list of tuples of integers, each of which has length `n`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import *
            sage: TrivialGrading( 3, "t" ).basis('t')
            [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
            sage: TrivialGrading( 3, "t" ).basis( None )
            []
        """
        if index == self.__index :
            if vars is None :
                vars = list(range(self.__ngens))
            
            return [ tuple(i*[0] + [1] + (self.__ngens - i - 1)*[0])
                     for i in vars ]
        else :
            return [] 
        
    def subgrading(self, gens) :
        r"""
        The grading of same type for the ring with only the variables given by
        ``gens``.
        
        INPUT:
            - ``gens`` - A list of integers.
        
        OUTPUT:
            An instance of :class:~`.TrivialGrading`.
        
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import *
            sage: TrivialGrading( 3, "t" ).subgrading([1,2])
            Trivial grading on 2 generators with index 't'
            sage: TrivialGrading( 3, "t" ).subgrading([])
            Trivial grading on 0 generators with index 't'
        """
        return TrivialGrading(len(gens), self.__index)
    
    def __contains__(self, x) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import *
            sage: g = TrivialGrading( 3, "t" )
            sage: "t" in g
            True
            sage: (1,2,1) in g
            False
            sage: None in g
            False
        """
        return x == self.__index
    
    def __cmp__(self, other) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import *
            sage: g = TrivialGrading( 3, "t" )
            sage: g == TrivialGrading( 3, "t" )
            True
            sage: g == TrivialGrading( 2, "t" )
            False
            sage: g == TrivialGrading( 3, None )
            False
        """
        c = cmp(type(self), type(other))
        
        if c == 0 :
            c = cmp(self.__ngens, other.__ngens)
        if c == 0 :
            c = cmp(self.__index, other.__index)
            
        return c
    
    def __hash__(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import *
            sage: hash( TrivialGrading( 3, "t" ) )
            1963142774     # 32-bit
            14848044662    # 64-bit
        """
        return reduce(xor, list(map(hash, [self.__ngens, self.__index])))
    
    def _repr_(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import *
            sage: TrivialGrading( 3, "t" )
            Trivial grading on 3 generators with index 't'
        """
        return "Trivial grading on %s generators with index %s" \
                % (self.__ngens, repr(self.__index))
    
    def _latex_(self) :
        r"""
        TESTS::
            sage: from psage.modform.fourier_expansion_framework.gradedexpansions.gradedexpansion_grading import *
            sage: latex( TrivialGrading( 3, "t" ) )
            \text{Trivial grading on $3$ generators with index $\verb|t|$}
        """
        return r"\text{Trivial grading on $%s$ generators with index $%s$}" \
                % (latex(self.__ngens), latex(self.__index))
