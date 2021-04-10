"""
Precision for Fourier expansions of Siegel modular forms
"""
from __future__ import absolute_import
from __future__ import division


from builtins import str
from builtins import range
from past.utils import old_div
from copy import deepcopy
import operator

from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.functions.other import floor
from sage.functions.other import ceil
from sage.misc.functional import isqrt
from sage.structure.sage_object import SageObject
from sage.structure.element import RingElement
from sage.rings.all import infinity
import sage.structure.element
from sage.misc.latex import latex

#from copy import deepcopy
#import operator



#from sage.rings.integer import Integer
#from sage.functions.other import (floor, ceil)
#from sage.misc.functional import isqrt
#from sage.structure.sage_object import SageObject
#from sage.rings.all import infinity
#from fastmult import reduce_GL


#load 'fastmult.spyx'
#from fastmult import reduce_GL



class SiegelModularFormPrecision (SageObject):
    r"""
    Stores information on the precision of the Fourier expansion of a Siegel
    modular form and provides checking.
    """
    def __init__(self, prec):
        r"""
        INPUT
            prec -- an integer or infinity or a triple (aprec, bprec, cprec) of integers.
                    or an instance of class SiegelModularFormPrecision
        NOTE
            We distinguish two types of precision bounds:
            'disc' or 'box'. prec indicates that all terms q^(a,b,c)
            with GL(2,Z)-reduced (a,b,c) are available such that
            either 4ac-b^2 < prec ('disc'-precision) or else
            one has componentwise
            (a,b,c) < (aprec, bprec, cprec)
            ('box'-precision).

        EXAMPLES::

            sage: from sage.modular.siegel.siegel_modular_form_prec import SiegelModularFormPrecision
            sage: prec = SiegelModularFormPrecision(101)
            sage: prec2 = SiegelModularFormPrecision(prec)
            sage: prec
            Discriminant precision for Siegel modular form with bound 101
            sage: prec2
            Discriminant precision for Siegel modular form with bound 101
            sage: prec == prec2
            True
            sage: prec3 = SiegelModularFormPrecision((5,5,5))
            sage: prec3
            Box precision for Siegel modular form with bound (5, 5, 5)
            sage: prec4 = SiegelModularFormPrecision(infinity)
            sage: prec4
            Infinity precision for Siegel modular form with bound +Infinity
            sage: prec5 = SiegelModularFormPrecision((5,0,5))
            sage: prec5.prec()
            (0, 0, 0)

        TESTS::

            sage: prec = SiegelModularFormPrecision(101)
            sage: TestSuite(prec).run()
        """
        ## copy constructor
        if isinstance(prec, SiegelModularFormPrecision):
            self.__type = deepcopy(prec.type())
            self.__prec = deepcopy(prec.prec())

        elif prec is infinity:
            self.__type = 'infinity'
            self.__prec = infinity
            
        elif isinstance(prec, (int, Integer)):
            self.__type = 'disc'
            prec = max(0, prec)
            Dmod = prec % 4
            if Dmod >= 2:
                self.__prec = Integer(prec - Dmod + 1)
            else:
                self.__prec = Integer(prec)
                
        elif isinstance(prec, tuple) and len(prec) == 3 and all([isinstance(e, (int, Integer)) for e in list(prec)]):
            self.__type = 'box'
            a,b,c = prec
            a = min(a, c)
            b = min(b, a)
            if b <= 0:
                self.__prec = (0,0,0)
            else:
                self.__prec = (Integer(a), Integer(b), Integer(c))
                
        else:
            raise TypeError("incorrect type {0} for prec".format(type(prec)))
    
    def _repr_(self):
        r"""
        Returns the repr of self

        EXAMPLES::

            sage: from sage.modular.siegel.siegel_modular_form_prec import SiegelModularFormPrecision
            sage: prec = SiegelModularFormPrecision(101)
            sage: prec._repr_()
            'Discriminant precision for Siegel modular form with bound 101'
        """

        if self.type() == 'infinity':
            n = "Infinity"
        elif self.type() == 'disc':
            n = "Discriminant"
        elif self.type() == 'box':
            n = "Box"
        else:
            raise RuntimeError("Unexpected value of self.__type")
            
        return "{0} precision for Siegel modular form with bound {1}".format(n, repr(self.__prec))

    
    def _latex_(self):
        r"""
        Returns the latex of self

        EXAMPLES::

            sage: from sage.modular.siegel.siegel_modular_form_prec import SiegelModularFormPrecision
            sage: prec = SiegelModularFormPrecision(101)
            sage: prec._latex_()
            'Discriminant precision for Siegel modular form with bound $101$'

        """

        if self.__type == 'infinity':
            n = "Infinity"
        elif self.__type == 'disc':
            n = "Discriminant"
        elif self.__type == 'box':
            n = "Box"
        else:
            raise ValueError("Unexpected value of self.__type")

        return "{0} precision for Siegel modular form with bound ${1}$".format(str(n),str(latex(self.__prec)))


    def type(self):
        r"""
        Returns the type of self: box, disc or infinity

        EXAMPLES::

            sage: from sage.modular.siegel.siegel_modular_form_prec import SiegelModularFormPrecision
            sage: prec = SiegelModularFormPrecision((2,2,2))
            sage: prec = SiegelModularFormPrecision(infinity)
            sage: prec.type()
            'infinity'
            sage: prec = SiegelModularFormPrecision(101)
            sage: prec.type()
            'disc'

        """
        return self.__type

    
    def prec(self):
        r"""
        Returns the tuple of the maximal box, the integer for the maximal discriminant
        and infinity otherwise.  If self is of type disc then it will return the 
        largest fundamental discriminant just below.
    
        EXAMPLES::
     
            sage: from sage.modular.siegel.siegel_modular_form_prec import SiegelModularFormPrecision
            sage: prec = SiegelModularFormPrecision(11)
            sage: prec.prec()
            9
            sage: prec = SiegelModularFormPrecision(8)
            sage: prec.prec()
            8


        """
        return self.__prec

    
    def is_in_bound(self, t):
        r"""
        Return True or False accordingly as the GL(2,Z)-reduction of
        the tuple (a,b,c)=t is within the bounds of this precision.

        EXAMPLES::

            sage: from sage.modular.siegel.siegel_modular_form_prec import SiegelModularFormPrecision
            sage: prec = SiegelModularFormPrecision(101)
            sage: prec.is_in_bound((1,0,1))
            True
            sage: prec.is_in_bound((0,0,25))
            True
            sage: prec.is_in_bound((0,0,26))
            False
            sage: prec.is_in_bound((6,0,6))
            False
            sage: prec = SiegelModularFormPrecision((2,2,2))
            sage: prec.is_in_bound((1,0,1))
            True
            sage: prec.is_in_bound((2,2,2))
            False


        TODO 
            fix this

        NOTE
            It is assumed that (a,b,c) is semi positive definite
        """
        (a, b, c) = t

        from .fastmult import reduce_GL
        if self.__type == 'infinity':
            return True        
        elif self.__type == 'disc':
            (a,b,c) = reduce_GL(a, b, c)
            if a == 0 and b == 0:
                return c < self.get_contents_bound_for_semi_definite_forms()
            D = 4*a*c-b**2
            if D < 0:
                return False
            if D > 0:
                return D < self.__prec
        elif self.__type == 'box':
            (a, b, c) = reduce_GL(a, b, c)
            return a < self.__prec[0] and b < self.__prec[1] and c < self.__prec[2]
        else:
            raise RuntimeError("Unexpected value of self.__type")        


    def get_contents_bound_for_semi_definite_forms(self):
        r"""
        Return the bound for the contents of a semi positive definite
        form falling within this precision.

        EXAMPLES::

            sage: from sage.modular.siegel.siegel_modular_form_prec import SiegelModularFormPrecision
            sage: prec = SiegelModularFormPrecision(7)
            sage: prec.get_contents_bound_for_semi_definite_forms()
            2
            sage: prec = SiegelModularFormPrecision(101)
            sage: prec.get_contents_bound_for_semi_definite_forms()
            26


        NOTE
            If (a,b,c) is semi positive definite, then it is GL(2,Z)-equivalent
            to a the form (0,0,g) where g is the contents (gcd) of a,b,c.
   
            I don't know what this does.  It seems to only do it for singular forms-- NR
        """
        if self.__type == 'infinity':
            return infinity
        elif self.__type == 'disc':
            return ceil((self.__prec+1)/4)
        elif self.__type == 'box':
            return self.__prec[2]
        else:
            raise RuntimeError("Unexpected value of self.__type")        


    def __lt__(self, other):
        r"""
        Return True if the set of GL(2,Z)-reduced forms within the precision self
        is contained in but not equal to the corresponding set of forms within the
        precision other.

        EXAMPLES::
            
            sage: from sage.modular.siegel.siegel_modular_form_prec import SiegelModularFormPrecision
            sage: prec1 = SiegelModularFormPrecision(101)
            sage: prec2 = SiegelModularFormPrecision(100)
            sage: prec3 = SiegelModularFormPrecision(201)
            sage: prec1 
            Discriminant precision for Siegel modular form with bound 101
            sage: prec2
            Discriminant precision for Siegel modular form with bound 100
            sage: prec2 < prec1
            True
            sage: prec3 < prec1
            False
            sage: prec1.__lt__(prec3)
            True

        """
        if not isinstance(other, SiegelModularFormPrecision):
            raise NotImplementedError("can only compare with elements of the same class")

        ## TODO: implement comparison of precisions of different types
        if self.__type != other.__type:
            return False

        if self.__type == 'box':
            a,b,c = self.__prec
            ao,bo,co = other.__prec
            return a <= ao and b <= bo and c <= co \
                    and any([a < ao, b < bo, c < co])

        elif self.__type == 'disc':
            return self.__prec < other.__prec

        else:
            raise RuntimeError("Unexpected value of self.__type")

        
    def __le__(self, other):
        r"""
        Returns whether self <= other

        EXAMPLES::

            sage: from sage.modular.siegel.siegel_modular_form_prec import SiegelModularFormPrecision
            sage: prec1 = SiegelModularFormPrecision(101)
            sage: prec2 = SiegelModularFormPrecision(100)
            sage: prec3 = SiegelModularFormPrecision(101)
            sage: prec1 <= prec3
            True
            sage: prec1 <= prec2
            False
            sage: prec1 <= prec1
            True
            sage: prec4 = SiegelModularFormPrecision(infinity)
            sage: prec5 = SiegelModularFormPrecision((5,5,5))
            sage: prec5.__le__(prec4)
            False

        TO DO:

            It seems eevrything should be <= infinity

        """
        return self.__lt__(other) or self.__eq__(other)

    
    def __eq__(self, other):
        r"""
        Returns self == other as defined by having the same prec and type

        EXAMPLES::

            sage: from sage.modular.siegel.siegel_modular_form_prec import SiegelModularFormPrecision
            sage: prec1 = SiegelModularFormPrecision(101)
            sage: prec2 = SiegelModularFormPrecision(100)
            sage: prec3 = SiegelModularFormPrecision(101)
            sage: prec4 = SiegelModularFormPrecision(infinity)
            sage: prec5 = SiegelModularFormPrecision((5,5,5))
            sage: prec1 == prec2
            False
            sage: prec1 == prec3
            True
            sage: prec4 == prec2
            False
            sage: prec4 == prec4
            True
            sage: prec5.__eq__(prec5)
            True

        """

        if not isinstance(other, SiegelModularFormPrecision):
            raise NotImplementedError("can only compare with elements of the same class")
        
        return (self.__type == other.type() and self.__prec == other.prec())

    
    def __ne__(self, other):       
        r"""
        Returns self != other

        EXAMPLES::

            sage: from sage.modular.siegel.siegel_modular_form_prec import SiegelModularFormPrecision
            sage: prec1 = SiegelModularFormPrecision(101)
            sage: prec2 = SiegelModularFormPrecision(100)
            sage: prec3 = SiegelModularFormPrecision(101)
            sage: prec1.__ne__(prec1)
            False
            sage: prec1 != prec2
            True
            sage: prec1 != prec3
            False
        """

        return not self.__eq__(other)

    
    def __gt__(self, other):
        r"""
        Returns self > other

        EXAMPLES::

            sage: from sage.modular.siegel.siegel_modular_form_prec import SiegelModularFormPrecision
            sage: prec = SiegelModularFormPrecision(101)
            sage: prec2 = SiegelModularFormPrecision(100)
            sage: prec > prec2
            True 
            sage: prec2.__gt__(prec) 
            False
        """

        return other.__lt__(self)

    
    def __ge__(self, other):
        r"""
        Returns self >= other

        EXAMPLES::

            sage: from sage.modular.siegel.siegel_modular_form_prec import SiegelModularFormPrecision
            sage: prec = SiegelModularFormPrecision(101)
            sage: prec2 = SiegelModularFormPrecision(100)
            sage: prec >= prec2
            True
            sage: prec.__ge__(prec)
            True
        """

        return self.__eq__(other) or other.__lt__(self)

    
    def __iter__(self):
        r"""
        Iterate over all GL(2,Z)-reduced  semi positive forms
        which are within the bounds of this precision.

        EXAMPLES::

            sage: from sage.modular.siegel.siegel_modular_form_prec import SiegelModularFormPrecision
            sage: prec = SiegelModularFormPrecision(11)
            sage: for k in prec.__iter__(): print k
            (0, 0, 0)
            (0, 0, 1)
            (0, 0, 2)
            (1, 0, 1)
            (1, 0, 2)
            (1, 1, 1)
            (1, 1, 2)


        NOTE
            The forms are enumerated in lexicographic order.
        """
        if self.__type == 'disc':
            bound = self.get_contents_bound_for_semi_definite_forms()
            for c in range(0, bound):
                yield (0,0,c)

            atop = isqrt(self.__prec // 3)
            if 3*atop*atop == self.__prec: atop -= 1
            for a in range(1,atop + 1):
                for b in range(a+1):
                    for c in range(a, ceil((b**2 + self.__prec)/(4*a))):
                        yield (a,b,c)

        elif 'box' == self.__type:
            (am, bm, cm) = self.__prec
            for a in range(am):
                for b in range(min(bm,a+1)):
                    for c in range(a, cm):
                        yield (a,b,c)

        else:
            raise RuntimeError("Unexpected value of self.__type")
            
        raise StopIteration

        
    def positive_forms(self):
        r"""
        Iterate over all GL(2,Z)-reduced  strictly positive forms
        which are within the bounds of this precision.

        EXAMPLES::

            sage: from sage.modular.siegel.siegel_modular_form_prec import SiegelModularFormPrecision
            sage: prec = SiegelModularFormPrecision(11)
            sage: for k in prec.positive_forms(): print k
            (1, 0, 1)
            (1, 0, 2)
            (1, 1, 1)
            (1, 1, 2)

        NOTE
            The forms are enumerated in lexicographic order.
        """
        if self.__type == 'disc':
            atop = isqrt(self.__prec // 3)
            if 3*atop*atop == self.__prec: atop -= 1
            for a in range(1,atop + 1):
                for b in range(a+1):
                    for c in range(a, ceil((b**2 + self.__prec)/(4*a))):
                        yield (a,b,c)

        elif 'box' == self.__type:
            (am, bm, cm) = self.__prec
            for a in range(am):
                for b in range(min(bm,a+1)):
                    for c in range(a, cm):
                        yield (a,b,c)

        else:
            raise RuntimeError("Unexpected value of self.__type")
            
        raise StopIteration
