#*****************************************************************************
#       Copyright (C) 2008 Nils-Peter Skoruppa <nils.skoruppa@uni-siegen.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

r"""
Implementation of the category of finite quadratic modules.

A quadratic module is a pair $(M,Q)$, where $M$ is a finite abelian
group and where $Q:M\rightarrow \Q/Z$ is a quadratic form. The latter
means that $Q(tx)=t^2Q(x)$ for all integers $t$ and all $x$ in $M$,
and that $B(x,y):=Q(x+y)=Q(x)-Q(y)$ defines a bilinear map
$B: M \times M \rightarrow \Q/Z$. A morphism $f:A \rightarrow B$
between quadratic modules $A=(M,Q) and $B=(N,R)$ is a homomorphism
of abelian groups $F:M \rightarrow N$ such that $R \circ f = Q$. 

Quadratic modules show up naturally in the theory of Weil
representations, which in turn are needed e.g. in the theory of Jacobi
forms or elliptic modular forms of integral or half integral
weight. This is due to the fact, that every irreducible representation of $SL(2,\ZZ)$
whose kernel is a congruence subgroup is contained in a Weil
representation.


TECHNICAL NOTE

A finite abelian group $M$ is given as a set of $n$
generators $a,b,c..,$ and an $n \times n$-matrix $R$ encoding the
relations among the generators: $(a,b,c,...)R = 0$. A (finite) quadratic form
on $M$ is then given by its Gram matrix
$$
    G \equiv \frac 12
    \begin{pmatrix}
    B(a,a) & B(a,b) & B(a,c) & \dots
    \\ B(b,a) & B(b,b) & B(b,c) \dots
    \\
    \dots \end{pmatrix}
    \bmod \ZZ^{n \times n}.
$$
The finite quadratic module is thus isomorphic to
$$
    (\ZZ^n/R\ZZ^n, x+R\ZZ^n \mapsto x^tGx)
$$
via the map$(x_a, x_b, \dots) + R\ZZ^n \mapsto x_a a + x_b b + \cdots.$
Accordingly a typical initialization of a finite quadratic module would
be to provide the integer matrix $R$ and the rational matrix $G$.


REMARK

Many of the mathematically more meaningful methods of the class FiniteQuadraticModule_ambient assume that the
represented finite quadratic module is nondegenerate (i.e. $B(x,y)=0$
for all $y$ in $M$ is only possible for $x=0$). Applying such a method
to a degenerate module will raise an exeption (TODO: what exception?).


REFERENCES

    [Sko] Nils-Peter Skoruppa, Finite quadratic modules and Weil representations,
          in preparation 2008

    TODO: find other references


    
AUTHORS:
    -- Hatice Boylan,      <boylan@mathematik.uni-siegen.de>
    -- Martin Frick,       <frick@mathematik.uni-siegen.de>
    -- Lars Fischer,       <lars.fischer@student.uni-siegen.de>
    -- Shuichi Hayashida,  <hayashida@mathematik.uni-siegen.de>
    -- Nils-Peter Skoruppa <nils.skoruppa@uni-siegen.de>

    -- Fredrik Stroemberg <fredrik314@gmail.com>
    -- Stephan Ehlen <ehlen@mathematik.tu-darmstadt.de>
    -- Sebastian Opitz <opitz@mathematik.tu-darmstadt.de>
    
    The CN-Group started the development of this package as a seminar
    project at the university of Siegen. Its initial members have been:
    Hatice Boylan, Martin Frick, Lars Fischer, Shuichi Hayashida,
    Saber Mbarek, Nils-Peter Skoruppa

\section{Tutorial}

TODO: Lots and lots of examples. 
"""

from sage.functions.other                 import floor
from sage.arith.all                       import divisors, is_prime, kronecker, lcm, gcd, prime_divisors, primitive_root, is_square, is_prime_power, inverse_mod, binomial
from sage.rings.all                       import ZZ, QQ, Integer, PolynomialRing,CC
from sage.groups.group                    import AbelianGroup 
from sage.modules.free_module_element     import vector
from sage.matrix.matrix_space             import MatrixSpace
from sage.matrix.constructor              import matrix, diagonal_matrix, identity_matrix
from sage.rings.number_field.number_field import CyclotomicField
from sage.structure.sage_object           import SageObject
from sage.structure.element               import AdditiveGroupElement
from sage.structure.sequence              import Sequence_generic
from sage.structure.all                   import Sequence
from sage.all                             import copy,cached_method,is_even,is_odd,Sequence,prod,uniq,valuation,randrange,is_fundamental_discriminant,xmrange,QuadraticField,xgcd,cartesian_product
from sage.structure.category_object import normalize_names
from sage.graphs.graph import DiGraph
from sage.rings.number_field.number_field_element import NumberFieldElement

###################################
## CLASS QUAD_MODULE
###################################


class FiniteQuadraticModule_ambient (AbelianGroup):
    r"""
    Describes a finite quadratic module $(A,Q)$. The abelian group $A$
    is given by generators \code{(e_0,e_1,\dots)} and a matrix $R$ of relations
    satisfied by the generators (hence $(a,b,\dots)\cdot R = 0$). The
    quadratic form $Q$ is given by the Gram matrix w.r.t. the generators;
    more precisely, the $(i,j)$-th entry $g_{i,j}$ is a rational number
    such that $Q(e_i+e_j)-Q(r_i)-Q(e_j) = g_{i,j} + \ZZ$..

    NOTES ::
        The abelian group may also be thought of as $\ZZ^n/R\ZZ^n$,
        and the generators as $e_i = c_i + R\ZZ^n$, where $c_i$ is
        the canonical basis for $\ZZ^n$.

        In our implementation we think of elements of $\ZZ^n$
        as column vectors.

        Use FiniteQuadraticModule() for more flexibility when creating
        FiniteQuadraticModule_ambient objects.

    EXAMPLES ::
        sage: R = matrix(2,2,[2,1,1,2])
        sage: G = 1/2 * R^(-1)
        sage: A.<a,b> = FiniteQuadraticModule_ambient( R, G); A
        Finite quadratic module in 2 generators:
         gens: a, b
         form: 1/3*x0^2 + 2/3*x0*x1 + 1/3*x1^2
        sage: a == b
        True
        sage: a is b
        False
        sage: a + b in A
        True
    """

    def __init__(self, R, G, check = True, names = None):
        r"""
        Initialize a quadratic module from R and G.
        
        INPUT
            R     -- an integral non-degenerate square matrix of size $n$
                     representing the finite abelien group $\ZZ^n/R\ZZ^n$.
                    
            G     -- a symmetric rational matrix of same size as $R$ and such
                     that $R^tGR$ is half integral and $2*R*M$ is an
                     integer matrix, representing the quadratic
                     form $x + R\ZZ^n \mapsto G[x] + \ZZ$.
                    
            check -- True or False indicating whether R and G should be
                     checked to be a valid set of data for the definition of
                     a quadratic module.
                     
            names -- a string used to name the generators of
                     the underlying abelian group.
        """
        if check:
            #TODO: check if R, G are matrices over ZZ, QQ, admit nonsquare R with rank == size G
            if False == hasattr( R, '_matrix_') or False == hasattr( G, '_matrix_'):
                raise TypeError, "%s: not a matrix" %R
            if False == R.is_square() or 0 == R.det() or R.denominator() != 1:
                raise ValueError, "%s: list not a regular integral square matrix" %R
            if False == G.is_square() or False == G.is_symmetric():
                raise ValueError, "%s: list not a symmetric square matrix" %G
            C0 = G*R; C = R.transpose()*C0; v = vector([C[i,i] for i in range(C.nrows())])
            if C0.denominator() > 2 or C.denominator() > 2 or v.denominator() > 1:
                raise ValueError, "(%s,%s): not a compatible pair" %(R,G)
        AbelianGroup.__init__( self)
        # There is a sort of bug:
##         sage: M = matrix(ZZ,1,[1])       
##         sage: M1 = matrix(ZZ,1,1,{(0,0):1})
##         sage: M.block_sum (M1)
##         ---------------------------------------------------------------------------
##         TypeError                                 Traceback (most recent call last)
##         ..............
##         TypeError: Cannot convert sage.matrix.matrix_integer_sparse.Matrix_integer_sparse to sage.matrix.matrix_integer_dense.Matrix_integer_dense
        # which we cure for the moment as follows:
        MS = MatrixSpace (ZZ,R.nrows(),R.ncols())
        R = MS(R)
        # We keep the initializing pair $(R,G)$
        
        self.__iM = R
        self.__iG = self._reduce_mat(G);

        # We replace $__iM$ by the unique $__R$ in $__iM * GL(n,\ZZ)$ which is
        # in lower Hermite normal form (i.e. is lower triangular and the rows
        # are reduced modulo their rightmost nonzero element).
        
        self.__R = matrix( ZZ, self.__iM).transpose().echelon_form().transpose()

        # For simplicity and effectiveness in various internal computations
        # we use an equivalent form $(__E,__J)$ of our quadratic module,
        # where $__E$ is the diagonal matrix formed from the elementary divisors of $__R$
        # in descending order, and where superfluous $1$'s are thrown out.
        # The system of generators $e_i + __E\ZZ^m$, where $e_i$ is the standard basis of $\ZZ^m$
        # are in the sequel called 'the fundamental system of generators'. 
        # TODO: In addition, $J$ should be put in Jordan form
        
        D,U,V = matrix( ZZ, self.__R).dense_matrix().smith_form();

        # Hence we have $D = U * __R * V$
        
        mask = []
        for n in range(D.nrows()):
            if D[n,n] > 1:
                mask.append(n)
        if 0 == len(mask):
            mask.append(0)
        self.__E = D.matrix_from_rows_and_columns( mask, mask).sparse_matrix()
        T  = U**(-1)
        self.__J = self._reduce_mat(
            (T.transpose()
             *self.__iG
             *T).matrix_from_rows_and_columns( mask, mask))
        n = self.__E.nrows()
        self.__elementary_divisors = tuple( [ self.__E[j,j] for j in range( n) ])

        # Transformation matrices:
        # can_sys = fun_gen * C2F, fun_sys = can_sys * F2C

        self.__C2F = U.matrix_from_rows( mask)
        self.__F2C = T.matrix_from_columns( mask)
        
        # Set the relation, Gram matrix and ngens to be used for the output

        # self.__R = self.__R
        self.__G = self.__iG
        self.__ngens = self.__R.ncols()
        
        # define inerited ngens attribute 
        if None == names:
            names = "e"
            names = normalize_names(self.ngens(),names)
        self._assign_names(names)

        ## Silly class identifier needed since our class does not keep its name in sage....
        self._is_FiniteQuadraticModule_ambient=True

        # zero of self
        self._zero = FiniteQuadraticModuleElement(self, 0, can_coords = True)
        # list of possible x_c's
        self._xcs={}
        
        
    
    ###################################
    ## Introduce myself ...
    ###################################


    def _latex_( self):
        r"""
        EXAMPLES
            sage: A = FiniteQuadraticModule( '3^2.27^-3'); latex(A)
            \left(\left\langle \Z^{5}/\left(\begin{array}{rrrrr}
            3 & 0 & 0 & 0 & 0 \\
            0 & 3 & 0 & 0 & 0 \\
            0 & 0 & 27 & 0 & 0 \\
            0 & 0 & 0 & 27 & 0 \\
            0 & 0 & 0 & 0 & 27
            \end{array}\right)\Z^{5} \right\rangle, \frac{1}{3} x_{0}^{2} + \frac{1}{3} x_{1}^{2} + \frac{2}{27} x_{2}^{2} + \frac{1}{27} x_{3}^{2} + \frac{1}{27} x_{4}^{2}\right)
        """
        n = self.ngens()
        v = vector( PolynomialRing(QQ, 'x', n).gens())
        form = v.dot_product( self.gram() * v)
        return '\\left(\\left\\langle \\Z^{%s}/%s\\Z^{%s} \\right\\rangle, %s\\right)'\
               %(latex(n), latex(self.__R), latex(n), latex(form))


    def _repr_(self):
        r"""
        EXAMPLES
            sage: A = FiniteQuadraticModule('2^2'); A
            Finite quadratic module in 2 generators:
             gens: e0, e1
             form: 1/2*x0*x1
        """
        n = self.ngens()
        v = vector( PolynomialRing(QQ, 'x', n).gens())
        gens = ', '.join([x for x in self._names])
        form = v.dot_product( self.gram() * v)
        return 'Finite quadratic module in %s generators:\n gens: %s\n form: %s' \
               %(n, gens, form)
 

    ###################################
    ## Providing struct. defining items
    ###################################


    def ngens( self):
        r"""
        Return the number of generators of the underlying abelian group.

        EXAMPLES
            sage: F = matrix( QQ, 3, 3, [ 2, 1, 5, 1, 34, 19, 5, 19, 6]); F
            [ 2  1  5]
            [ 1 34 19]
            [ 5 19  6]
            sage: A = FiniteQuadraticModule( F); A
            Finite quadratic module in 3 generators:
             gens: e0, e1, e2
             form: 157/1960*x0^2 + 891/980*x0*x1 + 13/1960*x1^2 + 151/980*x0*x2 + 33/980*x1*x2 + 1893/1960*x2^2
            sage: A.ngens()
            3
        """
        return self.__ngens


    def gen( self, i=0):
        r"""
        Return the $i$-th generator of the underlying abelian group.
        
        EXAMPLES
            sage: A = FiniteQuadraticModule( [2, 4, 8]); A
            Finite quadratic module in 3 generators:
             gens: e0, e1, e2
             form: 1/8*x0^2 + 1/16*x1^2 + 1/32*x2^2
            sage: A.gens()
            (e0, e1, e2)
            sage: A.1
            e1
        """
        x = [0]*self.ngens()
        x[int(i)] = 1
        return FiniteQuadraticModuleElement(self, x, can_coords = True)

    def gens(self):
        r"""
        Return a tuple of generators for self.
        """
        return tuple([ self.gen(i) for i in range(self.ngens())])

    def relations( self):
        r"""
        Return a matrix in Hermite normal form describing the relations
        satisfied by the generators (see class dos string for details).

        EXAMPLES
            sage: A = FiniteQuadraticModule('2^-2'); A
            Finite quadratic module in 2 generators:
             gens: e0, e1
             form: 1/2*x0^2 + 1/2*x0*x1 + 1/2*x1^2
            sage: A.relations()
            [2 0]
            [0 2]
        """
        return self.__R


    def gram( self):
        r"""
        Return the Gram matrix with respect to the generators
        (as rational matrix).
        
        EXAMPLES NONE
        """
        
        return self.__G

    
    def elementary_divisors( self):
        r"""
        Return the orders of the generators of the underlying group.

        EXAMPLES
            sage: A = FiniteQuadraticModule([11,33]); A  
            Finite quadratic module in 2 generators:
             gens: e0, e1
             form: 1/44*x0^2 + 1/132*x1^2
            sage: A.elementary_divisors ()
            (66, 22)
        """
        return self.__elementary_divisors


    def fgens( self):
        r"""
        Return a fundamental system for the underlying abelian group.

        NOTES
            A fundamental system of a finite abelian group $A$ is a
            set of generators $a_i$ such that $A$ equals the direct sum of
            the cyclic subgroups $\langle a_i \rangle$ generated by the
            $a_i$, and if, for each $i$, the order of $a_i$ equals the
            $i$-th elementary divisor of $A$.

            This method returns a fundamental system (which is, in fact,
            the one which was chosen
            when the quadratic module was initialized, and with respect to
            which all internal computations are actually performed).

        EXAMPLES NONE
        """
        return tuple([ self( list(x), can_coords = False) for x in identity_matrix( ZZ, len(self.elementary_divisors())) ])


    ###################################
    ## Coercion
    ###################################


    def __call__( self, x, can_coords = False):
        r"""
        Coerce object into an appopriate child object
        of self if possible.

        We ccoerce
        - an element of this module,
        - a list of coordinates with respect to the
          fundamental generators,
        - the integer $0$.

        EXAMPLES NONE
        """
        if isinstance( x, FiniteQuadraticModuleElement):
            if x.parent() is self:
                return x
        if isinstance( x, list):
            return FiniteQuadraticModuleElement( self, x, can_coords)
        if isinstance( x, (Integer, int, long)) and 0 == x:
            return FiniteQuadraticModuleElement( self, x)
        raise TypeError, "cannot coerce %s to an element of %s" %(x, self)


    ###################################
    ## Invariants
    ###################################


    def order( self):
        r"""
        If self is the quadratic module $(M,Q)$, return the order of $M$.

        EXAMPLES
            sage: A = FiniteQuadraticModule([11,33]);
            sage: A.order()
            1452
        """
        return prod( e for e in self.elementary_divisors())

    
    def exponent( self):
        r"""
        If self is the quadratic module $(M,Q)$, then return the exponent of the
        abelian group $M$.

        EXAMPLES
            sage: A = FiniteQuadraticModule([11,33]);
            sage: A.exponent()
            66
        """
        return max( self.elementary_divisors())
    
    
    def level( self):
        r"""
        If self is the quadratic module $(M,Q)$, then return the smallest positive integer $l$
        such that $l\cdotQ = 0$.

        EXAMPLES
            sage: A = FiniteQuadraticModule([11,33]);
            sage: A.level()
            132
        """
        H = copy(self.__J)
        for i in range(H.ncols()):
            for j in range( i+1, H.ncols()):
                H[i,j] = 2*H[i,j]
                H[j,i] = H[i,j]
        return H.denominator()


    def tau_invariant( self, p = None):
        r"""
        Return +1  or -1 accordingly  as the order of the underlying abelian group
        (resp. the largest power of $p$ dividing this order)
        is a perfect square or not.
        
        EXAMPLES NONE
        """
        q = self.order()
        if p is None:
            return +1 if q.is_square() else -1
        return +1 if is_even(q.valuation(p)) else -1

    @cached_method
    def sigma_invariant( self, p = None):
        r"""
        If this quadratic module equals $A=(M,Q)$, return
        $\sigma(A) = \sqrt{|M|}^{-1/2}\sum_{x\in M} \exp(-2\pi i Q(x))$

        EXAMPLES NONE
        """
        return self.char_invariant( -1, p)[0]

        
    def witt_invariants(self):
        r"""
        Return the family $\{sigma( A(p)\}_p$ as dictionary,
        where $A$ is this module, $A(p)$ its $p$-part,
        and $p$ is running through the divisors of the exponent of $A$
        (see also A.sigma_invariant()).
        
        EXAMPLES NONE
##             sage: A = FiniteQuadraticModule([11,33]);
##             sage: A.witt_class()[3]
##             zeta8^2
        """
        P = prime_divisors( self.exponent())
        d = dict()
        for p in P:
            s = self.sigma_invariant(p)
            t = self.tau_invariant(p)
            d[p] = (s, t);
        return d


    def char_invariant(self, s, p=None, debug=0):
        r"""
        If this quadratic module equals $A = (M,Q)$, return
        the characteristic function of $A$ (or $A(p)$ if $p$ is a prime)
        at $s$, i.e. return
        $$\chi_A (s)= |M|^{-1}\sum_{x\in M} \exp(2\pi i s Q(x))).$$

        NOTE
            We apply the formula in [Sko, Second Proof of Theorem 1.4.1].

        EXAMPLES NONE
        """
        s = s % self.level()
        if s == 0:
            return 1, 1
        if not p is None and not is_prime(p):
            raise TypeError
        if p and 0 != self.order() % p:
            return 1,1
        _p = p
        K = CyclotomicField (8)
        z = K.gen()
        jd = self.jordan_decomposition()
        ci = ci1 = 1
        for c in jd:
            # c: basis, ( prime p,  valuation of p-power n, dimension r, determinant d over p [, oddity o]) 
            p,n,r,d = c[1][:4]
            if debug > 0: print "c=",c
            if debug > 0: print "p={0}, n={1}, r={2}, d={3}".format(p,n,r,d)
            if _p and p != _p:
                continue
            o = None if 4 == len(c[1]) else c[1][4]
            if debug > 0: print "o=",o
            k = valuation( s, p)
            s1 = Integer(s/p**k)
            h = max(n-k,0)
            if debug > 0: print "h={0}".format(h)
            q = p**h
            if p!=2:
                lci  = z**((r*(1-q))%8) * d**(h%2) if h > 0 else 1
                lci1 = q**(-r) if h > 0 else 1
            elif k == n and o != None:
                #print "t!"
                return 0,0
            else:
                f = z**o if o else 1 
                lci = f * d**(h%2) if h > 0 else 1
                lci1 = q**(-r) if h > 0 else 1
                if debug > 0: print f, d, lci 
            if 2 == p: lci = lci**s1
            if debug > 0: print "lci=",lci
            if debug > 0: print "lci1=", lci1
            ci *= lci * kronecker( s1, 1/lci1)
            ci1 *= lci1
        return ci, QQ(ci1).sqrt()

    def signature(self,p=-1):
        r"""
        Compute the p-signature of self.
        p=-1 is the real signature.
        """
        if p == -1:
            p = None
        inv = self.char_invariant(1,p)
        inv = inv[0].list()
        if inv.count(1)>0:
            return inv.index(1)
        else:
            return inv.index(-1) + 4
        
    ###################################
    ## Deriving quadratic modules
    ###################################


    def __add__( self, B):
        r"""
        Return the orthogonal sum of quadratic modules $A + B$,
        where $A$ is this quadratic module.

        INPUT
            B -- quadratic module

        EXAMPLES NONE
##             sage: A = FiniteQuadraticModule([11,33]);
##             sage: B = A + A; B

        """
        return FiniteQuadraticModule_ambient( self.relations().block_sum( B.relations()), \
                            self.gram().block_sum( B.gram()), \
                            False)


    def __mul__( self, _n):
        r"""
        Return the $n$ fold orthogonal sum of this quadratic module..

        EXAMPLES NONE
        """
        n = int(_n)
        if n != _n:
            raise TypeError, "Argument n (= %s) must be an integer."%n
        if n < 0:
            raise ValueError, "Argument n (= %s) must be nonnegative."%n
        if 0 == n:
            return FiniteQuadraticModule()
        if n > 0:
            # TODO: Understand why
            # 1) return sum( self for j in range( n))
            #    does not work, though it works for FiniteQuadraticModuleElements
            #    Explanation for TODO above:
            #      LF: it seems that sum works like tmp=int(0) (with a real int, not Integer),
            #      and then aggregate the sum in tmp
            #      that means the first addition is 0 + quadmodule_object,
            #      which was undefined until now ( implemented __radd__) 
            #     _ For elements there must be some coercion (through inheritance?) or
            #      an inherited addition with ints
            
            # 2) 3*A does not work, but 3*e does work for elements e.
            #      LF: thanks to _r_action 3*A works now, the elements seem to have
            #      an inherited _r_action()

            # old implementation without using sum
            #A = self
            #for i in range(1, n): 
            #    A = A + self
            #return A
            
            return sum( self for j in range( n))
            # sum works now, but depends on an existing __radd__ implementation


    def __radd__(self, thing_on_the_left):
        r"""
        TODO: what is this here?
        Implements 0 + self = self.

        Nothing more. This addition seems to be needed in the sum inside the
        implementation of for example 3*FiniteQuadraticModule().

        EXAMPLES NONE
        """
        
        #print "in __radd__ with a", type(thing_on_the_left), thing_on_the_left
        if (thing_on_the_left == 0):
            return self # 0 + self should be self, right?
        else:
            return ValueError, "Argument n on the left (= %s) must be 0."%n
            # what is for example 1 + FiniteQuadraticModule() ?


    # neither _rmul_ nor __rmul__ work for 3*A
    def _r_action(self, n): 
        r"""
        TODO: whaT is this ?
        Return the $n$ fold orthogonal sum of self.

        This method is used for an expression like 3*A, where A is \emph{on the right}.
        It calls self.__mul__.

        EXAMPLES:
            ::
            sage: A=FiniteQuadraticModule();
            sage: 3*A == A*3
            True
        """
        
        # HINT:
        # !less /usr/local/sage/devel/sage/sage/structure/element.pyx:
        #    __mul__ uses  return coercion_model.bin_op_c(left, right, mul)
        #    with coercion_model = sage.structure.coerce.CoercionModel_cache_maps()
        # then coercion_model??, look for bin_op_c(, then  if op is mul, ....
       
        if type(n) is Integer:
            return self.__mul__(n)
        else:
            raise TypeError, "Argument n (= %s) must be an integer."%n

        
    def __cmp__( self, other):
        r"""
        Return 1 if other is a quadratic module having the same generator names,
        satisfying the same relations and having the same Gram matrix as this module.

        EXAMPLES NONE

        TODO: compare names
        """
        # compare two quadmodules via their relations and Gram-Matrix
        if type(other) is type(self): 
            if (self.__R == other.__R) and (self.gram() == other.gram()):
                return 0
        return -1
        #Why is 0 returned? in the explanation ist says return 1

    def twist( self, s):
        r"""
        If self is the quadratic module $A = (M,Q)$, return the twisted module $A^s = (M,s*G)$.

        INPUT
            s -- an integer

        OUTPUT
            quadratic module
        
        EXAMPLES
            sage: A = FiniteQuadraticModule( '11^-1'); A
            Finite quadratic module in 1 generators:
             gens: e
             form: 10/11*x^2
            sage: B = A^-1; B                           
            Finite quadratic module in 1 generators:
             gens: e
             form: 1/11*x^2
            sage: C = B^11; C                           
            Finite quadratic module in 1 generators:
             gens: e
             form: 0
            sage: D = FiniteQuadraticModule( '4^-2.5^2'); D
            Finite quadratic module in 4 generators:
             gens: e0, e1, e2, e3
             form: 1/4*x0^2 + 1/4*x0*x1 + 1/4*x1^2 + 1/5*x2^2 + 1/5*x3^2
            sage: E = D^2; E
            Finite quadratic module in 4 generators:
             gens: e0, e1, e2, e3
             form: 1/2*x0^2 + 1/2*x0*x1 + 1/2*x1^2 + 2/5*x2^2 + 2/5*x3^2
            sage: E.is_nondegenerate ()
            False
        """
        a = Integer(s)
        if a != s:
            raise TypeError, "Argument a (= %s) must be an integer."%a
        return FiniteQuadraticModule_ambient( self.relations(), a*self.gram())


    def orthogonal_basis( self, p = None):
        r"""
        Return an orthogonal system of generators for the
        underlying group of this quadratic module, if $p$ is None,
        respectively for the $p$-Sylow subgroup if $p$ is a prime.

        NOTES
            See FiniteQuadraticModule_subgroup.orthogonal_basis()
            for detailed explanation.
            
        EXAMPLES
            sage: A.<a,b,c,d,e,f,g,h,j> = FiniteQuadraticModule( '11^-7.2^-2')
            sage: A.orthogonal_basis (11)
            [2*a, 2*b, c, d, e, f, g]
            sage: A.orthogonal_basis (2) 
            [j, h]

            sage: R.<X>= ZZ['X']
            sage: K.<x> = NumberField( X^10 + X^9 - X^7  - X^6 - X^5 - X^4  - X^3 + X + 1)
            sage: L = FiniteQuadraticModule( (1-x)/1001); L

            Finite quadratic module in 10 generators:
             gens: e0, e1, e2, e3, e4, e5, e6, e7, e8, e9
             form: 1/91*x0^2 + 997/1001*x0*x1 + 1000/1001*x1^2 + 999/1001*x0*x2 + 2/1001*x1*x2 + 998/1001*x2^2 + 2/1001*x0*x3 + 995/1001*x1*x3 + 999/1001*x3^2 + 995/1001*x0*x4 + 997/1001*x2*x4 + 10/1001*x3*x4 + 1000/1001*x4^2 + 997/1001*x1*x5 + 10/1001*x2*x5 + 999/1001*x3*x5 + 993/1001*x4*x5 + 997/1001*x5^2 + 997/1001*x0*x6 + 10/1001*x1*x6 + 999/1001*x2*x6 + 993/1001*x3*x6 + 993/1001*x4*x6 + 12/1001*x5*x6 + 993/1001*x6^2 + 10/1001*x0*x7 + 999/1001*x1*x7 + 993/1001*x2*x7 + 993/1001*x3*x7 + 12/1001*x4*x7 + 985/1001*x5*x7 + 8/1001*x6*x7 + 1/1001*x7^2 + 999/1001*x0*x8 + 993/1001*x1*x8 + 993/1001*x2*x8 + 12/1001*x3*x8 + 985/1001*x4*x8 + 8/1001*x5*x8 + 2/1001*x6*x8 + 981/1001*x7*x8 + 1/1001*x8^2 + 993/1001*x0*x9 + 993/1001*x1*x9 + 12/1001*x2*x9 + 985/1001*x3*x9 + 8/1001*x4*x9 + 2/1001*x5*x9 + 981/1001*x6*x9 + 2/1001*x7*x9 + 989/1001*x8*x9 + 4/1001*x9^2
            sage: og_b = L.orthogonal_basis(); og_b long time
            [77*e0,
             77*e0 + 462*e9,
             77*e0 + 77*e8 + 693*e9,
             77*e0 + 693*e7 + 385*e8 + 231*e9,
             77*e0 + 693*e5 + 154*e7 + 308*e8 + 77*e9,
             77*e0 + 693*e5 + 616*e6 + 231*e7 + 462*e8 + 77*e9,
             77*e0 + 308*e4 + 462*e5 + 539*e6 + 616*e7 + 77*e8 + 462*e9,
             77*e0 + 308*e3 + 154*e4 + 154*e5 + 924*e6 + 770*e7 + 462*e8 + 308*e9,
             77*e0 + 770*e2 + 770*e3 + 231*e4 + 77*e5 + 693*e6 + 924*e7 + 77*e8 + 77*e9,
             77*e0 + 77*e1 + 231*e2 + 308*e3 + 539*e4 + 847*e5 + 231*e6 + 539*e7 + 847*e8 + 385*e9,
             91*e1,
             91*e0 + 455*e9,
             91*e0 + 819*e8 + 637*e9,
             91*e0 + 819*e7 + 455*e8 + 546*e9,
             91*e0 + 273*e6 + 637*e7 + 182*e8 + 455*e9,
             91*e0 + 364*e5 + 364*e6 + 273*e7 + 455*e8 + 455*e9,
             91*e0 + 364*e4 + 91*e5 + 364*e6 + 910*e8 + 455*e9,
             91*e0 + 637*e3 + 182*e5 + 364*e6 + 728*e7 + 728*e8 + 182*e9,
             91*e0 + 91*e2 + 455*e3 + 910*e4 + 728*e5 + 91*e6 + 455*e7 + 455*e8 + 819*e9,
             91*e0 + 455*e1 + 455*e3 + 637*e4 + 364*e5 + 273*e6 + 455*e7 + 728*e8 + 819*e9,
             143*e0,
             143*e1 + 429*e9,
             143*e0 + 715*e8 + 715*e9,
             143*e0 + 858*e7 + 143*e8 + 429*e9,
             143*e0 + 858*e4 + 286*e7 + 572*e8 + 715*e9,
             143*e0 + 858*e4 + 858*e6 + 572*e7 + 429*e8 + 429*e9,
             143*e0 + 429*e3 + 858*e4 + 572*e6 + 143*e7 + 572*e8 + 858*e9,
             143*e0 + 429*e3 + 858*e4 + 143*e5 + 858*e6 + 572*e7 + 572*e8,
             143*e0 + 858*e2 + 429*e3 + 858*e4 + 858*e5 + 429*e6 + 286*e7 + 286*e9,
             143*e0 + 143*e1 + 286*e2 + 143*e3 + 572*e5 + 715*e6 + 143*e7 + 858*e8 + 143*e9]
         """
        if not (self.is_nondegenerate()):
            raise ValueError
        if p is None:
            U = self.subgroup( self.gens())
        elif is_prime(p):
            U = self.subgroup( p)
        else:
            raise TypeError
        return U.orthogonal_basis()


    def jordan_decomposition( self):
        r"""
        """
        try:
            return self.__jd
        except AttributeError:
            self.__jd = JordanDecomposition( self)
        return self.__jd

    
    def spawn( self, gens, names = None):
        r"""
        Spawn the subgroup generated by the elements of the list
        gens equipped with the quadratic form induced by this module as finite
        quadratic module.

        EXAMPLES NONE
        """
        return self.subgroup( gens).as_ambient( names)
        

    def quotient( self, U):
        r"""
        Return the quotient module $self/U$ for the isotropic subgroup
        $U$ of self. If $U$ is not isotropic an excepion is thrown.

        INPUT
            U -- a subgroup of self
        
        OUTPUT
            quadratic module
          
        EXAMPLES NONE
##             sage: A = FiniteQuadraticModule([11,33]);
##             sage: A2=A.quotient(list(A.isotropic_subgroups())[0]); A2
##             Finite quadratic module ([33, 11], 1/33*x0^2 + 1/11*x0*x1 + 1/11*x1^2) with generators (e0, e1)

        NOTES
            Let $U^\sharp = K\ZZ^n/M\ZZ^n$ the dual of the subgroup $U =
            H\ZZ^n/M\ZZ^n$ of $M$. The quotient module
            $(U^\sharp/U, x + U \mapsto Q(x) + \ZZ)$
            is then isomorphic to $(\ZZ^n/K^{-1}H\ZZ^n, G[Kx])$.

        """
        if not isinstance(U, FiniteQuadraticModule_subgroup) or U.ambience() is not self or not U.is_isotropic():
            raise ValueError, "%s: not an isotropic subgroup" %U
        V = U.dual()
        K = matrix( V)
        return FiniteQuadraticModule( K**(-1)*matrix(U), K.transpose()*self.__J*K)


    def __div__( self, U):
        r"""
        Return the quotient of $A/U$.

        EAMPLES NONE
        """
        return self.quotient( U)

    
    def __pow__( self, s):
        r"""
        Return the twist $A^s$.

        EXAMPLES NONE
        """
        return self.twist( s)
    

    def anisotropic_kernel( self):
        r"""
        Return the anisotropic quotient of this qudaratic module,
        i.e. if this module is $A$ then return $A/U$, where $U$
        is a maximal isotropic subgroup.

        OUTPUT
            (A/U, f, g), where $U$ is a maximal isotropic subgroup,
            and, where $f:(U^#,Q) \rightarrow A$ and $g:(U^#,Q) \rightarrow A/U$
            ($Q$ denotes the quadratic form of $A$) are the natural morohisms of quadratic modules.
            
        EXAMPLES NONE
        
        TODO:
        Just find a maximal isotropic subgroup
        and return the quotient.
        """
        raise NotImplementedError


    ###################################
    ## Deriving subgroups
    ###################################

    
    def subgroup( self, arg = []):
        r"""
        Return a subgroup of the underlying abelian group $U$ of this quadratic module.
        Return the subgroup of $U$ generated by the elements in arg if arg is a list or tuple.
        Return the $p$-Sylow subgroup if $arg$ is a prime $p$.
        
        INPUT
           arg -- a list of elements of this quadratic module or a prime number

        EXAMPLES
            sage: A.<a,b,c,d,e,f,g> = FiniteQuadraticModule( '11^-3.2_2^4')
            sage: A2 = A.subgroup (2); A2                                     
            < g, f, e, d >
            sage: A3 = A.subgroup (3); A3
            < 0 >
            sage: A11 = A.subgroup (11); A11                                  
            < 2*a, 2*b, 2*c >
            sage: A11.order()
            1331
            sage: A2.order() 
            16
        """
        if isinstance( arg, (list, tuple)):
            if [] == list(arg):
                arg = [self(0)]
            return FiniteQuadraticModule_subgroup( list(arg))
        p = Integer(arg)
        if is_prime(p):
            U = FiniteQuadraticModule_subgroup( list( self.gens()))
            U = U.split(p)[0] if 0 == U.order()%p else self.subgroup([self(0)])
            return U
        raise ValueError
    

    def kernel( self):
        r"""
        Return the dual subgroup of the underlying group of this module,
        i.e. return $\{y \in A : B(y,x) = 0 \text{ for all } x \in A  \}$,
        for this module $(A,Q)$.

        EXAMPLES
            sage: A.<a,b> = FiniteQuadraticModule( [3,3], [1/3,1/3,1/3]); A
            Finite quadratic module in 2 generators:
             gens: a, b
             form: 1/3*x0^2 + 1/3*x0*x1 + 1/3*x1^2
            sage: U = A.kernel(); U
            < a + b >
            sage: B = A.quotient(U); B
            Finite quadratic module in 2 generators:
             gens: e0, e1
             form: 1/3*x0^2 + 1/3*x0*x1 + 1/3*x1^2
            sage: e0,e1 = B.gens()
            sage: e0
            2*e1
            sage: B.jordan_decomposition().genus_symbol()
            '3^-1'
        """
        return self.dual_group( self.subgroup(self.gens()))


    def dual_group( self, U):
        r"""
        Return the dual subgroup
        $U^\sharp = \{y \in A : B(y,x) = 0 \text{ for all } x \in U  \}$
        of the subgroup $U$ of $self = (A,Q)$

        INPUT
            U -- a subgroup of this quadratic module
          
        EXAMPLES NONE

        NOTES
            If  the  dual  group  (w.r.t. the fundamental system) is  given
            as $K\ZZ^r/E\ZZ^n$ then the
            columns of $K$ form a basis for the integral solutions  of
            $2H^tJx \in \ZZ^n$. We solve this  by the trick of augmenting
            $2H^tJx$ by the unit matrix and solving the  corresonding
            system of linear equations over $\ZZ$
        """
        H = matrix( U)
        X = 2*H.transpose()*self.__J
        n = len(self.elementary_divisors())
        Y = X.augment( MatrixSpace(QQ,n).identity_matrix())
        K0 = matrix(ZZ, Y.transpose().integer_kernel().matrix().transpose())
        K = K0.matrix_from_rows( range( n ))        
        l= [ FiniteQuadraticModuleElement(self, x.list() , can_coords=False ) for x in K.columns()  ]
        return self.subgroup(l)


    def kernel_subgroup(self,c):
        r"""
        Return the subgroup D_c={ x in D | cx=0}
        
        """
        if not c in ZZ:
            raise ValueError("c has to be an integer.")
        if gcd(c,self.order())==1:
            return self.subgroup([])
        l=[]
        for x in self:
            y = c*x
            if y==self(0):
                l.append(x)
        return self.subgroup(l)

    def power_subgroup(self,c):
        r"""
        Compute the subgroup D^c={c*x | x in D}
        
        """
        l=[]
        for x in self:
            y = c*x
            if y not in l:
                l.append(y)
        return self.subgroup(l)

    def power_subset_star(self,c):
        r"""    
        Compute the subset D^c*={x in D | cQ(y)+(x,y)=0, for all y in D_c}
        Using D^c* = x_c + D^c
        """
        xc = self.xc(c)
        Dc = self.power_subgroup(c)
        res=[]
        if xc==self._zero:
            return list(Dc)
        for x in Dc:
            res.append(x+xc)
        return res
        
    
    def xc(self,c):
        r"""
        Compute all non-zero values of the element x_c in the group D

        INPUT:

        - 'c'~- integer 

        OUTPUT:

        - 'x_c' element of self such that D^{c*} = x_c + D^{c} 
        
        """
        
        if is_odd(c):
            return self._zero
        k = valuation(c,2)
        return self._xc(k)
        
    def _xc(self,k):
        r"""
        Compute $x_c such that D^{c*} = x_c + D^{c} for c s.t. 2^k || c$
        """
        if self._xcs=={}:
            self._compute_all_xcs()
        if self._xcs.has_key(k):
            return self._xcs[k]
        return self._zero
            
    def  _compute_all_xcs(self):
        r"""
        Computes all non-zero values of the element x_c in the group D
        OUPUT: dictionary k => x_c  where 2^k || c
        """
        J=self.jordan_decomposition()
        res=dict()
        res[0 ]=0 
        for c in J:
            xc=0 
            p,k,r,d = c[1 ][:4 ]
            t = None if 4  == len(c[1 ]) else c[1 ][4 ]
            if(p<>2  or t==None): 
                continue
            q=2**k # = 2^n
            JC=J.constituent(q)
            CL=JC[0 ]
            # print "CL(",q,")=",JC[0]
            HOM=JC[1 ]
            # print "t=",t
            # We have an odd component
            for y in JC[0 ].orthogonal_basis():
                # print "basis = ",y
                z=JC[1 ](y)
                # print "basis = ",z
                xc=xc+(2 **(k-1 ))*z
            res[k]=xc
        self._xcs = res
        
    
    ###################################
    ## Predicates
    ###################################


    def is_multiplicative( self):
        r"""
        Return False since finite quadratic modules are modules.
        
        EXAMPLES
            sage: var('x')
            x
            sage: K.<a> = NumberField( x^3-x+1)
            sage: A = FiniteQuadraticModule(1/(a+11))
            sage: jd = A.jordan_decomposition()
            sage: jd.genus_symbol()
            '1319'
            sage: A.char_invariant(1)
            (-zeta8^2, 1/sqrt(1319))
            sage: A.is_multiplicative()
            False
        """
        return False

    
    def is_nondegenerate( self):
        r"""
        Return True or False accordingly if this module is
        non degenerate or not.
        NOTE: Is a trivial module non-degenerate (according to def)?
        EXAMPLES
            sage: N = FiniteQuadraticModule(); N

            Finite quadratic module in 1 generators:
             gens: 0
             form: 0
            sage: N.is_nondegenerate()
            True
        """
        if self.kernel().order() == 1:
            return True
        else:
            return False
        #return True if 1 == self.kernel().order() else False

    def non_trivial(self):
        r"""
        Check if self is trivial.
        """
        return self.order>1
        
        
    def is_isomorphic( self, A):
        r"""
        Return True or False accordingly as this quadratic module
        is isomorphic to $A$ or not.

        EXAMPLES
            sage: A = FiniteQuadraticModule( '2_2^4')
            sage: B = FiniteQuadraticModule( '2_6^-4')
            sage: A.is_isomorphic(B)
            True

        TODO
            Extend this function so to include also degenerate modules.

            Maybe wo do not need to check all divisors?
        """
        if False == self.is_nondegenerate() or False == A.is_nondegenerate():
            raise TypeError, 'the quadratic modules to compare must be non degenerate'
        if self.elementary_divisors() != A.elementary_divisors():
            return False
        divs = self.level().divisors()
        for t in divs:
            if self.char_invariant(t) != A.char_invariant(t):
                return False
        return True


    def is_witt_equivalent( self, A):
        r"""
        Return True or False accordingly as this quadratic module
        is Witt equivalent to $A$ or not.

        INPUT
            A -- quadratic module

        EXAMPLES NONE

        TODO
            Extend this function so to include also degenerate modules.       
        """
        if False == self.is_nondegenerate() or False == A.is_nondegenerate():
            raise TypeError, 'the quadratic modules to compare must be non degenerate'
        if self.tau_invariant() == A.tau_invariant():
            a = self.witt_invariants()
            b = A.witt_invariants()
            if a == b:
                return True
        return False


    ###################################
    ## Deriving other structures
    ###################################


    def as_discriminant_module( self):
        r"""
        Return a half-integral matrix $F$ such that
        this quadratic module is isomorphic to the
        discrimnant module
        $$D_F = (\ZZ^n/2F\ZZ^n, x + 2F\ZZ^n \mapsto \frac 14 F^{-1}[x] + \ZZ).$$

        EXAMPLES NONE

        TODO: This will be fun !!! (Look up Wall if his poof is effective ...)

        NOTE
            If $D_F$ and $D_G$ are isomorphic then there exist even unimodular matrices $U$ and $V$ such that
            $F+U$ and $G+V$ are equivalent (viewed as quadratic forms) over $\Z$ (see [Sko]). 
        """
        pass
    

    def orthogonal_group( self):
        r"""
        Returns the orthogonal group of this
        quadratic module.

        EXAMPLES NONE

        TODO
            Nice topic for a master thesis.
        """
        pass


    ###################################
    ## Iterators
    ###################################


    def __iter__( self):
        r"""
        Return a generator over the elements of self.

        EXAMPLES NONE return ( x for x in self if self.Q(x) == 0 )
        """
        return( self(x, can_coords = False) for x in xmrange( self.elementary_divisors()) ) 


    def values( self):
        r"""
        If this is $(M,Q)$, return the values of $Q(x)$ ($x \in M$) as a dictionary d.

        OUTPUT
            dictionary -- the mapping Q(x) --> the number of elements x with the same value Q(x)

        EXAMPLES NONE
##             sage: A=FiniteQuadraticModule([1,3])
##             sage: A.values()[7/12]
##             2

        TODO
            Replace the code by theoretical formulas
            using the Jordan decomposition.
            DONE: See the values function in the class JordanDecomposition  
        """
        valueDict = {}
        
        for x in self:
            v= self.Q( x )
            if valueDict.has_key( v ):
                valueDict[v] +=1
            else:
                valueDict[v] = 1
                
        return valueDict
    

    def subgroups( self, d = None  ):
        r"""
        Return a list of all subgroups of $M$ of order $d$, where $M$
        is the underlying group of self, or of all subgroups if d is not set.

        INPUT
            d -- integer

        OUTPUT
            generator for a list of FiniteQuadraticModule_subgroup of order d

        EXAMPLES
            sage: A = FiniteQuadraticModule([1,3]); A      
            Finite quadratic module in 2 generators:
             gens: e0, e1
             form: 1/4*x0^2 + 1/12*x1^2
            sage: list(A.subgroups())             
            [< e1, e0 >,
             < e1 >,
             < e0 + e1 >,
             < 2*e1, e0 >,
             < 2*e1 >,
             < 3*e1, e0 >,
             < 3*e1 >,
             < e0 + 3*e1 >,
             < e0 >,
             < 0 >]
             sage: B.<a> = FiniteQuadraticModule( '25^-1'); B        
             Finite quadratic module in 1 generators:
              gens: a
              form: 2/25*x^2
             sage: I = [U for U in B.subgroups() if U.is_isotropic()]; I
             [< 5*a >, < 0 >]

        NOTES
            Subgroups are  internally represented by lower triangular
            matrices in Hermite normal form dividing $M$.

            Any   subgroup  $U$   of  a finite  abelian group $\ZZ^n/R\ZZ^n$   is  of   the  form
            $N\ZZ^n/R\ZZ^n$,  where  $N$  is  a regular  integral  square
            matrix of  size $n$ such  that $N$ left divides $M$,  (i.e. such
            that  $N^{-1}M$  is  an   integer  matrix).  The  coset  $N
            GL(n,\ZZ)$ is  uniquely determined by $U$.   Recall that any
            such  coset contains  exactly one lower triangular matrix
            $H = (h_[ij})$  in Hermite  normal form  (i.e.  $h_{ii}>0$ and
            $0\le h_{ij}$ < h_{ii}$ for $j < i$).

            Accordingly, this function sets up an iterator over all
            matrices in lower Hermite normal form left dividing the diagonal
            matrix whose diagonal entries are the elementary divisors of
            the underlying abelian group of this module
            (and whose determinants equal $self.order()/d$).

        TODO
            Find a more effective implementation.

            Introduce optional arguments which allow to iterate in addition effectively
            over all subgroups contained in or containig a certain subgroup, being isotropic etc.

            One can use short cuts using e.g. that $N$ defines
            an isotropic subgroups if and only if $N^tJN$ is half integral (where $J$ is the Gram matrix
            w.r.t the fundamental generators of this module.
        """
        elementary_divisors = self.elementary_divisors()
        N = len( elementary_divisors)

        Mat = ZZ_N_times_N = MatrixSpace(ZZ, N )

        def __enumerate_echelon_forms():
            """
            Return an iterator over all matrices in HNF
            left dividing self.__E.
            """
            # If e1, e_2, ...  denote the elemetary divisors of self,
            # we define a string like
            # ( [ [d_0,0], [n_10,d_1], ...]
            #     for d_0 in divisors(e1)
            #     for d_1 in divisors(e2)
            #     for n_10 in range(d_1), ... )
            # and evaluate it to obtain our generator.
            genString = '['
            forStringN = ''
            forStringD = ''

            for r in range(N):
                genString += ' ['
                for c in range(N):
                    if c < r:
                        v= 'n_'+str(r)+str(c)
                        genString+= v+','
                        forStringN+= 'for '+v+' in '+'xrange(d_'+str(r)+') '
                    elif c==r:
                        v = 'd_'+str(r)
                        genString+= v+','
                        forStringD+= 'for '+v+' in '+'divisors('+str(elementary_divisors[r])+') '
                    else:
                        genString+= '0,'
    
                genString = genString[:-1]+ '] ,'
            genString = genString[:-1] + ' ] '
            
            genExpression='( ' + genString + forStringD + forStringN +') '
            #print genExpression
            #print "evaluating:", genExpression
            return eval( genExpression )

        for h in __enumerate_echelon_forms():
            h1 = Mat(h)
            if FiniteQuadraticModule_subgroup._divides( h1, self.__E):
                f= FiniteQuadraticModule_subgroup( [FiniteQuadraticModuleElement( self, list(x), can_coords = False) for x in h1.transpose()])
                if d:
                    if f.order() == d:
                        yield f
                else:
                    # d == None means we return every subgroup
                    yield f
        
    
    ###################################
    ## Auxiliary functions
    ###################################


    @staticmethod
    def _reduce( r):
        r"""
        Return the fractional part of $x$.

        EXAMPLES
        sage: FiniteQuadraticModule_ambient._reduce( 297/100)
        97/100
        """
        return r - r.floor()
    

    @staticmethod
    def _reduce_mat( A):
        r"""
        Return the fractional part of the symmetric matrix $A$.

        EXAMPLES
            sage: F = matrix(3,3,[1,1/2,3/2,1/2,2,1/9,3/2,1/9,1]); F
            [  1 1/2 3/2]
            [1/2   2 1/9]
            [3/2 1/9   1]
            sage: FiniteQuadraticModule_ambient._reduce_mat(F)      
            [  0   0   0]
            [  0   0 1/9]
            [  0 1/9   0]
        """
        B = matrix(QQ, A.nrows(), A.ncols())
        for i in range(A.nrows()):
            for j in range(A.ncols()):
                if i == j:
                    B[i,j] = FiniteQuadraticModule_ambient._reduce( A[i,j])
                else:
                    B[i,j] = FiniteQuadraticModule_ambient._reduce( 2*A[i,j])/2
        return B


    @staticmethod
    def _rdc( R, x):
        r"""
        Returns the $y \equiv x \bmod R\,ZZ^n$ in the fundamental
        mesh $\{R\zeta : 0 \le \zeta_i < 1\}$.
    
        INPUT
            x --- an integral vector of length n
            R --- a regular nxn matrix in lower Hermite normal form

        NOTE
            It is absolutely necessary that $R$ comes in lower Hermite
            normal form.

        EXAMPLES NONE
        """
        y=[0]*R.nrows()
        k=[0]*R.nrows()
        k[0],y[0] = divmod( Integer( x[0]), Integer( R[0,0]))
        for i in range(1,R.nrows()):
            k[i],y[i] = divmod( Integer( x[i] - sum(R[i,j]*k[j] for j in range(i))), Integer( R[i,i]))
        return y


    def _f2c( self, x):
        r"""
        Transform coordinates w.r.t. to the internal fundamental system to the coordinates w.r.t. the generators.

        EXAMPLES NONE
        """
        v = self.__F2C*vector(ZZ,x)
        return list( self._rdc( self.__R, v))


    def _c2f( self, x):
        r"""
        Transform coordinates w.r.t. the generators to coordinates w.r.t. the internal fundamental system.

        EXAMPLES NONE
        """
        v = self.__C2F*vector(ZZ,x)
        return list( self._rdc( self.__E, v))


    def Q( self, x):
        r"""
        Return the value $Q(x)$ (as rational reduced mod $\ZZ$).

        INPUT
            x -- an FiniteQuadraticModuleElement of self (e.g. x = self.0*17 )

        OUTPUT
            rational number -- the value of the quadratic form on x
            
        EXAMPLES NONE
        """
        c = vector(x.list())
        return self._reduce( c.dot_product( self.__J * c))

    def B( self, x, y):
        r"""
        Return the value $B(x,y) = Q(x+y)-Q(x)-Q(y)$ (as rational reduced mod $\ZZ$).

        INPUT
            x, y -- FiniteQuadraticModuleElements of self

        OUTPUT
            rational number -- the value of the bilinear form on x and y
        

        EXAMPLES NONE
        """
        c = vector( x.list()); d = vector( y.list())
        return self._reduce( 2 * c.dot_product( self.__J * d))


    # TODO: Adapt Shuichi's Jordan decomposition to implement this
    @staticmethod
    def _diagonalize( G, n):
        r"""
        Diagonalizes the matrix $G$ over the localisation $Z_(n)$.
        
        INPUT
        G --- a matrix with entries in $Z$
        n --- an integer
        
        OUTPUT
        A matrix $H$ and a matrix $T$ in $GL(r,Z_(n))$
        such that $H = T^tGT$ is in {\em Jordan decomposition form}.

        NOTES
        Here {\em $H$ is in Jordan decomposition form} means the following:

        o $H$ is a block sum of matrices $B_j$ of size 1 or 2
        o For each prime power $q=p^l$ dividing $n$, one has
          - $H \bmod q$ is diagonal,
          - the sequence $gcd(q,H[i,i])$ is increasing,
          - for $i > 1$ one has $H[i,i] \equiv p^k \bmod q$ for some $k$
            unless $gcd(q, H[i-1,i-1]) < gcd(q, H[i,i])$,
        o if $q=2^l$ is the exact 2-power dividing $n$, then
          - $B_j$ is scalar or
            $B_j \equiv p^k [2,1;1,2] \bmod q$
            or $B_j \equiv p^k [0,1;1,0] \bmod q$ for some $k$.
          (- maybe later also a sign walk and oddity fusion normalisation)

        o An implementation could follow:
        [Sko 1] Nils-Peter Skoruppa, Reduction mod $\ell$ of Theta Series of Level $\ell^n$,
                arXiv:0807.4694
        """
        # matrices over $Z/nZ$: matrix( IntegerModRing(n), r, [...])
        # see also: http://www.sagemath.org/doc/html/ref/module-sage.rings.polynomial-quotient-ring-element.html
        pass
    

    ###################################
    ## Misc
    ###################################


    def an_element(self):
        r"""
        Returns an element of this finite quadratic module.
        """
        return sum(self.gens())
    
    def random_element(self, bound=None):
        """
        Return a random element of this group.

        EXAMPLES NONE
        """
        L= [ randrange(0, ed ) for ed in self.__elementary_divisors ]
        return FiniteQuadraticModuleElement( self, L )

    
    def cayley_graph(self):
        """
        Returns the cayley graph for this finite group, as a SAGE
        DiGraph object. To plot the graph with different colors

        EXAMPLES
            sage: A.<a,b> = FiniteQuadraticModule( matrix( QQ, 2, [2,1,1,6])); A
            Finite quadratic module in 2 generators:
             gens: 5*b, b
             form: 3/11*x0^2 + 10/11*x0*x1 + 1/11*x1^2
            sage: D = A.cayley_graph (); D                                     
            Digraph on 11 vertices
            sage: g = D.plot(color_by_label=True, edge_labels=True)              

        TODO
            Make the following work:
            D.show3d(color_by_label=True, edge_size=0.01, edge_size2=0.02, vertex_size=0.03)
            
            Adjust by options so that images is less cluttered
        
        NOTE
            Copied from sage.groups.group.FiniteGroup.cayley_graph().
        """
        #from sage.graphs.graph import DiGraph
        arrows = {}
        for x in self:
            arrows[x] = {}
            for g in self.gens():
                xg = x+g
                if not xg == x:
                    arrows[x][xg] = g
        return DiGraph(arrows, implementation='networkx')


    def _is_valid_homomorphism_(self, codomain, im_gens):
        r"""
        Return True if \code{im_gens} defines a valid homomorphism
        from self to codomain; otherwise return False.

        EXAMPLES NONE

        NOTE
            If $A = \ZZ^n/R$ and $B = \ZZ^m/S$, then \code{im_gens}
            defines a valid homomorphism if $S^-1MS$ is an integral
            matrix, where $M$ is the matrix whose colums are the coordinates
            of the elements of \code{im_gens} w.r.t. the can. system.
        """
        if not isinstance( im_gens, (tuple, list)):
            x = Sequence( [im_gens])
        x = Sequence( im_gens)
        if not x.universe() is codomain:
            raise TypeError, 'im_gens (=%s) must belong to %s' %( im_gens, codomain) 
        m = matrix([ x[i].c_list() for i in range(len(im_gens))]).transpose()
        return FiniteQuadraticModule_subgroup._divides( codomain.relations(), m *self.relations())


    def _Hom_(self, B, cat = None):
        r"""
        Return the set of morphism from this quadratic module to $B$.

        EXAMPLES
            sage: A = FiniteQuadraticModule( matrix( QQ, 2, [2,1,1,2]))
            sage: B = 2*A
            sage: S = A._Hom_(B)
            Set of Homomorphisms from Finite quadratic module in 2 generators:
            gens: e1, e1
            form: 1/3*x0^2 + 2/3*x0*x1 + 1/3*x1^2 to Finite quadratic module in 4 generators:
            gens: e1, e1, e3, e3
            form: 1/3*x0^2 + 2/3*x0*x1 + 1/3*x1^2 + 1/3*x2^2 + 2/3*x2*x3 + 1/3*x3^2
        """
        if not (cat is None or (cat is B.category() and cat is self.category())):
            raise NotImplementedError
        # if not isinstance( B, FiniteQuadraticModule_ambient):
        if not hasattr(B,'_is_FiniteQuadraticModule_ambient'):
            raise TypeError, "B (=%s) must be finte quadratic module."%B
        return FiniteQuadraticModuleHomset(self, B)

    
    def hom(self, x):
        r"""
        Return the homomorphism from this module to the parent module of the
        elements in the list which maps the $i$-th generator of this module
        to the element \code{x[i]}.

        INPUT
            x -- a list of $n$ elements of a quadratic module, where $n$ is the number
                 of generators of this quadratic module.

        EXAMPLES NONE
        """
        v = Sequence(x)
        B = v.universe()
        #if not isinstance( B, FiniteQuadraticModule_ambient):
        if not hasattr(B,'_is_FiniteQuadraticModule_ambient'):
            raise TypeError, "B (=%s) must have universe a FiniteQuadraticModule."%B
        return self.Hom(B)(x)



###################################
## MORPHISMS
###################################
    

from sage.groups.group_homset import GroupHomset_generic
from sage.categories.homset import HomsetWithBase
from sage.categories.commutative_additive_groups import CommutativeAdditiveGroups
from sage.categories.morphism import Morphism

class FiniteQuadraticModuleHomset (GroupHomset_generic):
    r"""
    Implements the set of morphisms of a quadratic module into another.
    
    EXAMPLES NONE
    """
    def __init__( self, A, B):
        r"""
        INPUT
            A, B -- finite quadratic modules

        EXAMPLES NONE
        """
        HomsetWithBase.__init__(self, A, B, CommutativeAdditiveGroups(), base = A.base_ring())

    
    def __call__(self, im_gens, check=True):
        """
        EXAMPLES NONE
        """
        try:
            return FiniteQuadraticModuleHomomorphism_im_gens( self, im_gens, check = check)
        except (NotImplementedError, ValueError), err:
            print err
            raise TypeError, "images (=%s) do not define a valid homomorphism"%im_gens



class FiniteQuadraticModuleHomomorphism_im_gens (Morphism):
    r"""
    Implements elements of the set of morphisms between two quadratic modules.
    
    EXAMPLES
        sage: A = FiniteQuadraticModule( matrix( QQ, 2, [2,1,1,2])); B = 2*A; S=A._Hom_(B); u,v = A.gens(); a,b,c,d = B.gens();
        sage: f = S([a,b]); f                             
        Homomorphism : Finite quadratic module in 2 generators:
         gens: e1, e1
         form: 1/3*x0^2 + 2/3*x0*x1 + 1/3*x1^2 --> Finite quadratic module in 4 generators:
         gens: e1, e1, e3, e3
         form: 1/3*x0^2 + 2/3*x0*x1 + 1/3*x1^2 + 1/3*x2^2 + 2/3*x2*x3 + 1/3*x3^2
        sage: f(u),f(v),f(2*v),f(0)                                                                               
        (e1, e1, 2*e1, 0)
    """
    
    def __init__( self, homset, im_gens, check = True):
        r"""
        INPUT
            homset  -- a set of modphisms between two quadratic modules  
            im_gens -- a list of elements of the codomain.
            
        EXAMPLES NONE
        """
        Morphism.__init__( self, homset) # sets the parent
        if not isinstance(im_gens, Sequence_generic):
            if not isinstance(im_gens, (tuple, list)):
                im_gens = [im_gens]
            #print im_gens[0].parent() is homset.codomain() 
            im_gens = Sequence( im_gens, homset.codomain())
        if check:
            if len(im_gens) != homset.domain().ngens():
                raise ValueError, "number of images must equal number of generators"
            t = homset.domain()._is_valid_homomorphism_( homset.codomain(), im_gens)
            if not t:
                raise ValueError, "relations do not all (canonically) map to 0 under map determined by images of generators."
        self.__im_gens = im_gens
        # For effeciency, we compute the images of the fundamental generators of the domain
        n = len( im_gens)
        self.__c_im_gens = []
        for x in homset.domain().fgens():
            cos = x.c_list()
            y = sum( im_gens[i]*cos[i] for i in  range( n))
            self.__c_im_gens.append( y)
        
        
    def kernel(self):
        r"""
        Return the kernel of this morphism.

        EXAMPLES NONE
        """
        pass


    def image(self, J):
        r"""
        Return the image of this morphism.

        EXAMPLES NONE
        """
        pass

    
    def _repr_(self):
        r"""
        EXAMPLES NONE
        """
        return "Homomorphism : %s --> %s"%(self.domain(),self.codomain())


    def _latex_(self):
        r"""
        EXAMPLES NONE
        """
        return "%s \\rightarrow{} %s"%(latex(self.domain()), latex(self.codomain()))


    def __call__( self, x):
        r"""
        EXAMPLES NONE
        """
        if not x in self.domain():
            raise TypeError, 'x (=%s) must be in %s' %(x,self.domain())
        n = len(self.__c_im_gens)
        cos = x.list()
        return sum( self.__c_im_gens[i]*cos[i] for i in range(n))



###################################
## CONVENIENCE FUNCTION
###################################


def _A( q, s = 1, **args):
    r"""
    Return the quadratic module $A_q^s$ as defined in [Sko].
    
    EXAMPLES
        sage: A = _A(3); A      
        Finite quadratic module in 1 generators:
         gens: e
         form: 1/3*x^2
        sage: B = _A(5,-2); B   
        Finite quadratic module in 1 generators:
         gens: e
         form: 3/5*x^2
        sage: C.<a> = _A(2^5); C
        Finite quadratic module in 1 generators:
         gens: a
         form: 1/64*x^2
        sage: a.order()
        32
    """
    q = Integer(q)
    s = Integer(s)
    if q > 1 and q.is_prime_power() and 1 == q.gcd(s):
        if is_odd(q):
            return FiniteQuadraticModule( [q], [s/q], **args)
        else:
            return FiniteQuadraticModule( [q], [s/(2*q)], **args)
    raise ValueError, 'q (=%s) must be prime power and s (=%s) relatively prime to q' %(q,s)


def _B( q, **args):
    r"""
    Return the quadratic module $C_q$ as defined in [Sko].
    
    EXAMPLES
        sage: B.<a,b> = _B(2^3); B
        Finite quadratic module in 2 generators:
         gens: a, b
         form: 1/8*x0^2 + 1/8*x0*x1 + 1/8*x1^2
        sage: (a+b).order()       
        8
    """
    q = Integer(q)
    if is_even(q) and q.is_prime_power():
        return FiniteQuadraticModule( [q,q], matrix( QQ, 2, [1/q, 1/(2*q), 1/(2*q), 1/q]), **args)
    raise ValueError, 'q = (%s) must be q power of 2'%q


def _C( q, **args):
    r"""
    Return the quadratic module $C_q$ as defined in [Sko].
    
    EXAMPLES
        sage: C.<a,b> = _C(2); C
        Finite quadratic module in 2 generators:
         gens: a, b
         form: 1/2*x0*x1
        sage: C.exponent()
        2
    """
    q = Integer(q)
    if is_even(q) and q.is_prime_power():
        return FiniteQuadraticModule( [q,q], matrix( QQ, 2, [0, 1/(2*q), 1/(2*q), 0]), **args)
    raise ValueError, 'q = (%s) must be q power of 2'%q


def _FiniteQuadraticModule_from_string( S, **args ):
    r"""
    Return the quadratic module described by the string S.
       
    INPUT
        S -- a string representing Jordan constituents of a finite quadratic modules

    NOTES
        The strings which will be accepted have the form
        $$
        'a^{\pm k}.b^{\pm l}.c^{\pm m}. \dots',
        $$
        where $a$, $b$, $c$, \dots are prime powers, and where $k$, $l$, etc.
        are positive integers (if an exponent $\pm k$ equals 1 it can be ommitted).
        If the $a$, $b$, \dots are powers of $2$ we admit also subscripts $t$, i.e. symbols
        of the form $a_t^{\pm k}$, where $0\le t < 8$ is an integer.

        The dot corresponds to the direct sum of quadratic modules, and the symbols $a^{\pm k}$, \dots
        have the following meaning:

        For a $2$-power $a$, the symbol $a^{+k}$, indicates the $k/2$-fold direct sum of the module $B = (\ZZ^2/a\ZZ^2, xy/a)$,
        whereas the symbol $a^{-k}$ denotes the module $(\ZZ^2/a\ZZ^2, (x^2+xy+y^2)/a)$ plus the $k/2-1$-fold sum of $(\ZZ^2/a\ZZ^2, xy/a)$.
        
        A symbol $a^{\pm k}$, for an odd prime power $a$, indicates the quadratic module
        $A_a^e+(k-1)A_a$, where $e$ is an integer such that the Legendre symbol $2^k e$ over $p$ equals $\pm 1$.

        A symbold $a_t^{\pm k}$, for a $2$-power $a$, indicates the quadratic module
        $A_a^{c_1}+\cdots + A_a^{c_k}$, where the $c_i$ are odd integers such that
        $D := c_1 \cdots c_k$ is a quadratic residue or non-residue modulo $8$ according to
        the sign $\pm$, and where $c_1 + \cdots + c_k \equiv t \bmod 8$. Here, for even $a$, we
        use $A_a = (\ZZ/a\ZZ, x^2/2a)$.

        Note that, for a symol $2^{\pm k}$, the $k$ must be even.
        Furthermore, a solution $(c_1,\dots,c_k)$ of the equations $\sum c_i \equiv t \bmod 8$ and 
        legendre symbol of 8 over $\prod c_i$ equal to $\pm 1$ exists if and only if
        $t \equiv k \bmod 2$, legendre symbol of $8$ over $t$ equal to $\pm 1$ for $k=1$,
        if $t\equiv 0 \bmod 8$ then $\pm 1 = +1$ and if $t\equiv 4 \bmod 8$ then \pm 1 = -1$
        for $k=2$. If any of these conditions is not fullfilled an error is raised.

    EXAMPLES
        sage: A.<a,b,c,d> =_FiniteQuadraticModule_from_string ('3^-1.3.5^2'); A
        Finite quadratic module in 4 generators:
         gens: a, b, c, d
         form: 2/3*x0^2 + 1/3*x1^2 + 1/5*x2^2 + 1/5*x3^2
        sage: A =_FiniteQuadraticModule_from_string ('8^+6'); A
        Finite quadratic module in 6 generators:
         gens: e0, e1, e2, e3, e4, e5
         form: 1/8*x0*x1 + 1/8*x2*x3 + 1/8*x4*x5
        sage: A.elementary_divisors ()
        (8, 8, 8, 8, 8, 8)
        sage: A =_FiniteQuadraticModule_from_string ('8_1^+3'); A
        Finite quadratic module in 3 generators:
         gens: e0, e1, e2
         form: 1/16*x0^2 + 1/16*x1^2 + 15/16*x2^2
        sage: A.elementary_divisors ()                           
        (8, 8, 8)
        sage: D.<a,b,c,d,e> = _FiniteQuadraticModule_from_string ('8_1^3.4^-2'); D
        Finite quadratic module in 5 generators:
         gens: a, b, c, d, e
         form: 1/16*x0^2 + 1/16*x1^2 + 15/16*x2^2 + 1/4*x3^2 + 1/4*x3*x4 + 1/4*x4^2
        sage: D.level(), D.exponent(), D.order(), D.elementary_divisors ()
        (16, 8, 8192, (8, 8, 8, 4, 4))
        sage: E.<a,b,c,d,e,f,g> = _FiniteQuadraticModule_from_string ('8_1^3.4^-2.3^-1.11^-1'); E
        Finite quadratic module in 7 generators:
         gens: a, b, c, d, e, f, g
         form: 1/16*x0^2 + 1/16*x1^2 + 15/16*x2^2 + 1/4*x3^2 + 1/4*x3*x4 + 1/4*x4^2 + 1/3*x5^2 + 1/11*x6^2
        sage: E.elementary_divisors ()
        (264, 8, 8, 4, 4)

    TECHNICAL NOTES
        The accepted strings can be described in BNF as follows:
        
        \begin{verbatim}
        <S>            ::= <symbol_list>
        <symbol_list>  ::= <symbol_list> "." <symbol> | <symbol>
        <symbol>       ::= <p-power> | <p-power>"^"<exponent> |
                           <2-power> "_" type | <2-power> "_" <type> "^" <exponent>
        <p-power>      ::= number
        <2-power>      ::= number
        <type>         ::= number
        <exponent>     ::= number | "+"number | "-"number
        \end{verbatim}

        Of course, we impose  the additional requirements that
        number is a positive integer, and number and type satisfy the above requierements.
    """
                           
    S= S.replace(' ','') # filter out spaces 
    List= S.split('.')
    
    ElementList = []
    for item in List:
        L1 = item.split("^")
        if len(L1) > 2: # more than one ^ in item
            raise ValueError
        elif len(L1) == 2:
            k = Integer(L1[1])
        else:
            k = 1
        L1= L1[0].split("_")
        a = Integer(L1[0])
        if len(L1) > 2: # more than one _ in item
            raise ValueError
        elif len(L1) == 2:
            if Integer(L1[1]) in range(8):
                t = Integer(L1[1])
            else:
                raise ValueError, "Type given, which ist not in 0..7: %s"%(L1[1] )
        else:
            t = None
        if not (k != 0 and a != 1 and a.is_prime_power()
                and ( None == t or (is_even(a) and t%2 == k%2))
                and ( not (None == t and is_even(a)) or 0 == k%2)
                ):
            raise ValueError,"{0} is not a valid signature!".format(S)
        c = None
        if is_odd(a):
            c = [1]*abs(k)
            p =  a.factor()[0][0]
            s = kronecker(2,p)**k
            if s*k < 0:
                c[0] = -1 if 3 == p%4 else primitive_root(p)
        if is_even(a) and t != None:
            if 1 == abs(k):
                if k == kronecker(t,2):
                    c = [t]
                else:
                    raise ValueError
            if abs(k) > 1:
                CP= eval( "cartesian_product([" + "[1,3,5,7]," *(abs(k)-1) + "])" )
                # TODO: find better algorithm
                e = 1 if k > 0 else -1
                for x in CP: 
                    s = sum(x)%8
                    if kronecker( prod(x)*(t-s),2) == e:
                        c = list(x)
                        c.append(t-s)
                        break
                if not c:
                    raise ValueError
        entry = {'a':a, 'k':k, 't':t, 'c':c}  
        ElementList.append( entry )

    names = args.pop('names', None)
    # TODO: Once the 0-module is cached replace the next 6 lines by: A = FiniteQuadraticModule()
    sym = ElementList[0]
    q = sym['a']; t = sym['t']; k = sym['k']
    if is_odd(q) or t != None:
        A = sum( _A(q,s,**args) for s in sym['c'])
    if is_even(q) and None == t:
        A = _C(q, **args)*(k//2) if k > 0 else _B(q, **args)
        if (-k)//2 > 1:
            A += _C(q, **args)*((-k)//2 - 1)
    for sym in ElementList[1:]:
        q = sym['a']; t = sym['t']; k = sym['k']
        if is_odd(q) or t != None:
            A += sum( _A(q,s,**args) for s in sym['c'])
        if is_even(q) and None == t:
            A += _C(q, **args)*(k//2) if k > 0 else _B(q, **args)
            if (-k)//2 > 1:
                A += _C(q, **args)*((-k)//2 - 1)
    A = FiniteQuadraticModule_ambient( A.relations(), A.gram(), names = names, **args)
    return A



def FiniteQuadraticModule( arg0=None, arg1=None, **args):
    r"""
    Create an instance of the class FiniteQuadraticModule_ambient.

    INPUT
        Supported formats:

        N.  FiniteQuadraticModule():
                the trivial quadratic module

        S.  FiniteQuadraticModule( string):
                the quadratic module $(L^#/L, B(x,x)/2)$, where $(L,B)$
                is a $\Z_p$-lattice encoded by the string as described
                in Conway-Sloane, p.???. TODO: look up ???

        L.  FiniteQuadraticModule( list):
                discriminant module constructed from the diagonal matrix
                with $2*x$ and $x$ running through list on its diagonal.
        
        M.  FiniteQuadraticModule( matrix):
                discriminant   module   constructed   from  a   regular
                symmetric even integral matrix.

        F.  FiniteQuadraticModule( number_field_element):
                For a nonzero $\omega$ in a numberfield $K$,
                the quadratic module $(\ZZ_K/A, x+A \mapsto tr( \omega x^2) + \ZZ)$,
                where $A$ is determined by
                $\omega D = B/A$ with relatively prime ideals $A$, $B$,
                and with $D$ denoting the different of $K$.
    
        LL. FiniteQuadraticModule( list_of_orders, list_of_coeffs):
                for a list of orders $[e_i]$  of size $n$ and a list of
                coeficients $[a_{ij}]$, the quadratic module
                $(\ZZ/e_1\times\cdots\times\ZZ/e_n,class(x)\mapsto\sum_{i\le j} a_{ij} x_i x_j)$.
        
        LM. FiniteQuadraticModule( list_of_orders, Gram_matrix):
                for  a  list  of  orders  $[e_i]$ of  size  $n$  and  a
                symmetric matric $G$, the quadratic module
                $(\ZZ/e_1 \times \cdots \times \ZZ/e_n, class(x) \mapsto G[x] + \ZZ)$.

        ML. FiniteQuadraticModule( matrix, list_of_coeffs):
                for a matrix $R$ of  size $n$ and a list of coefficients
                $[a_{ij}]$, the quadratic module
                $(\ZZ^n/R\ZZ^n, x+R\ZZ^n \mapsto \sum_{i\le j} a_{ij} x_i x_j)$.

        MM. FiniteQuadraticModule( matrix, Gram_matrix):
                for  a  matrix $R$  and  a  symmetric  matric $G$,  the
                quadratic module $(\ZZ^n/R\ZZ^n, x+R\ZZ^n \mapsto G[x] + \ZZ)$.


    EXAMPLES
        sage: N.<n> = FiniteQuadraticModule(); N                               
        Finite quadratic module in 1 generators:
         gens: 0
         form: 0
        sage: n.order()
        1

        sage: S.<x,y,z> = FiniteQuadraticModule( '7^-1.3.2_3^-1'); S
        Finite quadratic module in 3 generators:
         gens: x, y, z
         form: 6/7*x0^2 + 1/3*x1^2 + 3/4*x2^2

        sage: L.<w> = FiniteQuadraticModule( [13]); L          
        Finite quadratic module in 1 generators:
         gens: w
         form: 1/52*x^2

        sage: E8 = matrix( ZZ, 8, [4,-2,0,0,0,0,0,1,-2,2,-1,0,0,0,0,0,0,-1,2,-1,0,0,0,0,0,0,-1,2,-1,0,0,0,0,0,0,-1,2,-1,0,0,0,0,0,0,-1,2,-1,0,0,0,0,0,0,-1,2,0,1,0,0,0,0,0,0,2]); E8
        [ 4 -2  0  0  0  0  0  1]
        [-2  2 -1  0  0  0  0  0]
        [ 0 -1  2 -1  0  0  0  0]
        [ 0  0 -1  2 -1  0  0  0]
        [ 0  0  0 -1  2 -1  0  0]
        [ 0  0  0  0 -1  2 -1  0]
        [ 0  0  0  0  0 -1  2  0]
        [ 1  0  0  0  0  0  0  2]
        sage: M.<a,b,c,d,e,f,g,h> = FiniteQuadraticModule( 3*F); M
        Finite quadratic module in 8 generators:
         gens: a, b, c, d, e, f, g, h
         form: 1/3*x0^2 + 2/3*x0*x2 + 2/3*x1*x2 + 1/3*x0*x3 + 1/3*x1*x3 + 1/3*x3^2 + 2/3*x0*x5 + 2/3*x1*x5 + 1/3*x3*x5 + 2/3*x4*x5 + 1/3*x0*x6 + 1/3*x1*x6 + 2/3*x3*x6 + 1/3*x4*x6 + 1/3*x6^2 + 2/3*x0*x7 + 2/3*x2*x7 + 1/3*x3*x7 + 2/3*x5*x7 + 1/3*x6*x7 + 2/3*x7^2

        sage: X = QQ['X'].0                          
        sage: K.<x> = NumberField(X^4-8*X^3+1); K    
        Number Field in x with defining polynomial X^4 - 8*X^3 + 1
        sage: F.<a,b,c,d> = FiniteQuadraticModule((x^2-4)/7); F
        Finite quadratic module in 4 generators:
         gens: a, b, c, d
         form: 6/7*x0^2 + 1/7*x0*x1 + 5/7*x1*x2 + 5/7*x0*x3 + 6/7*x2*x3 + 3/7*x3^2

        sage: LL = FiniteQuadraticModule([3,4,30],[1/3,0,1/3,1/8,5/2,7/60]); LL
        Finite quadratic module in 3 generators:
         gens: e0, e1, e2
         form: 1/3*x0^2 + 1/8*x1^2 + 1/3*x0*x2 + 1/2*x1*x2 + 7/60*x2^2
        sage: LL.elementary_divisors ()                                        
        (60, 6)

        sage:  LL2.<u,v> = FiniteQuadraticModule( [5,5], [3/5,1/5,4/5]); LL2
        Finite quadratic module in 2 generators:
         gens: u, v
         form: 3/5*x0^2 + 1/5*x0*x1 + 4/5*x1^2
        sage: LL2.is_nondegenerate ()
        True

        sage: G = matrix(3,3,[1,1/2,3/2,1/2,2,1/9,3/2,1/9,1]); G
        [  1 1/2 3/2]
        [1/2   2 1/9]
        [3/2 1/9   1]
        sage: LM = FiniteQuadraticModule([4,9,18],G); LM        
        Finite quadratic module in 3 generators:
         gens: e0, e1, e2
         form: 2/9*x1*x2
        sage: LM.is_nondegenerate ()                            
        False

        sage: M = matrix( 2, [4,1,1,6]); M                               
        [4 1]
        [1 6]
        sage: ML.<s,t> = FiniteQuadraticModule( M, [3/23,-1/23,2/23]); ML
        Finite quadratic module in 2 generators:
         gens: 17*t, t
         form: 3/23*x0^2 + 22/23*x0*x1 + 2/23*x1^2
         
        sage: E = matrix( 2, [8,3,3,10]); E                       
        [ 8  3]
        [ 3 10]
        sage: MM.<x,y> = FiniteQuadraticModule( E, 1/2 * E^-1); MM
        Finite quadratic module in 2 generators:
         gens: 44*y, y
         form: 5/71*x0^2 + 68/71*x0*x1 + 4/71*x1^2
    """
    if arg0 is None:

        #N. FiniteQuadraticModule():
        if not 'check' in args:
            args['check'] = False
        return FiniteQuadraticModule_ambient( matrix(1,[1]), matrix(1,[0]), **args)
    elif hasattr(arg0,"_is_FiniteQuadraticModule_ambient"):
        return copy(arg0)

    if arg1 is None:

        if isinstance( arg0, str):
            #S FiniteQuadraticModule( string)
            if not 'check' in args:
                args['check'] = False 
            return _FiniteQuadraticModule_from_string( arg0, **args )

        elif isinstance(arg0, list):
            #L. FiniteQuadraticModule( list_of_orders):
            M = matrix( ZZ, len(arg0), len(arg0), \
                        dict([((i,i),2*arg0[i]) for i in range(len(arg0))]))
            G= QQ(1)/QQ(2) * M**(-1) 
        elif isinstance(arg0,(int,Integer)):
            M = matrix(ZZ,1,1,[arg0])
            G= QQ(1)/QQ(2) * M**(-1) 
            
        elif hasattr(arg0, '_matrix_'):
            #M. FiniteQuadraticModule( matrix):
            M = matrix( ZZ, arg0)     
            G = QQ(1)/QQ(2) * M**(-1)

        elif isinstance( arg0,NumberFieldElement):
            #F. FiniteQuadraticModule( number_field_element):
            if arg0.is_zero():
                raise ValueError, "%s: must be nonzero" %arg0    
            K = arg0.parent()
            d = K.different() 
            n = K.degree()
            l = K.integral_basis()        
            # Find G:
            G = matrix(QQ,n,n)
            for i in range(n):
               for j in range(n):
                 G[i,j]=(arg0*l[j]*l[i]).trace()

            # Compute the denominator ideal of omega*different:
            p = arg0*d
            s = p.factor()         
            A=K.ideal([1])
            for i in range(len(s)):
               s_factor=s[i]
               if s_factor[1] < 0:
                   A=A*s_factor[0]**(-s_factor[1])
                 
            # Compute M as the product of two matrices:                 
            L=matrix(QQ,n,n)
            for j in range(n):
                    for i in range(n):        
                      L[j,i]=A.integral_basis()[i].list()[j]
            E=matrix(QQ,n,n)
            for j in range(n):
                   for i in range(n):
                      E[j,i]=l[i].list()[j]
            M=E**(-1)*L      
        else:
            raise ValueError,"Can not construct finite quadratic module from {0}".format(arg0)
    else:

        if isinstance(arg0, list):
            M = matrix( ZZ, len(arg0), len(arg0), dict([((i,i),arg0[i]) for i in range(len(arg0))]))
        elif hasattr(arg0, '_matrix_'):
            M = matrix( ZZ, arg0);
        else:
            raise TypeError, "%s: should be None, a matrix or a list" %arg0           
       
        if isinstance(arg1, list):
            #LL./ML. FiniteQuadraticModule( list_of_orders/matrix, list_of_coeffs):
            n = M.nrows()
            G = matrix(QQ,n,n)
            i = j = 0
            for g in arg1:
                if i == j:
                    G[i,j] = g - floor(g);
                else:
                    G[i,j] = G[j,i] = QQ(g - floor(g))/QQ(2);
                j += 1
                if n == j:
                    i += 1; j = i
        elif hasattr(arg1, '_matrix_'):
            #LM./.MM FiniteQuadraticModule( list_of_orders/matrix, Gram matrix):
            G = matrix(QQ, arg1);
        else:
            raise TypeError, "%s: should be None, a matrix or a list" %arg1
    return FiniteQuadraticModule_ambient( M, G, **args)

# comments reviewed up to here

###################################
## CLASS QUAD_MODULE_ELEMENT 
###################################


class FiniteQuadraticModuleElement(AdditiveGroupElement):
    r"""
    Describes an element of a quadratic module.

    EXAMPLES NONE
##         sage: p = FiniteQuadraticModule([3,3])
##         sage: g0,g1 = p.gens()
##         sage: x = g0 + g1*5; x
##         e0 + 5*e1
##         sage: -x
##         5*e0 + e1
##         sage: x*1001 + g1
##         5*e0 + 2*e1
##         sage: p.Q(x)
##         1/6

    NOTES
        Code partly grabbed from sage.structure.element.AbelianGroupElement.
    """

    def __init__( self, A, x, can_coords = False):
        r"""
        Create the element of the FiniteQuadraticModule_ambient A whose coordinates
        with respect to the generators of $A$ (resp. fundamental
        generators if can_ccords = False)
        are given by the list x; create the zero element if $x=0$.

        INPUT
            A -- quadratic module
            x -- list or 0

        OUTPUT
            an element of the quadratic module A
        
        EXAMPLES NONE
##             sage: q = FiniteQuadraticModule( [2,4]);
##             sage: x = FiniteQuadraticModuleElement( q, [1,1])
##             sage: x
##             e0 + e1
        """
        AdditiveGroupElement.__init__( self, A)
        self.__repr = None
        ed = A.elementary_divisors()
        n = len(ed)
        if isinstance(x, (int, Integer)) and 0 == x:
            self.__intl_rep = [ 0 for i in range(n) ]
        elif isinstance(x, list):
            y = A._c2f( x) if can_coords == True else x
            self.__intl_rep = [ y[i]%ed[i] for i in range(n) ]
        else:
            raise TypeError, "Argument x (= %s) is of wrong type."%x


    ###################################
    ## Introduce myself ...
    ###################################


    def _latex_( self):
        r"""
        EXAMPLES NONE
        """
        s = ""
        A = self.parent()
        x = A.variable_names()
        n = len( A.variable_names())
        v = A._f2c( self.list())
        for i in range(n):
            if v[i] == 0:
                continue
            elif v[i] == 1:
                if len(s) > 0: s += " + "
                s += "%s"%x[i]
            else:
                if len(s) > 0: s += " + "
                s += "%s \cdot %s"%(latex(v[i]),x[i])
        if len(s) == 0: s = "0"
        return s

    
    def _repr_( self):
        r"""
        EXAMPLES NONE
        """
        s = ""
        A = self.parent()
        x = A.variable_names()
        n = len( A.variable_names())
        v = A._f2c( self.list())
        for i in range(n):
            if v[i] == 0:
                continue
            elif v[i] == 1:
                if len(s) > 0: s += " + "
                s += "%s"%x[i]
            else:
                if len(s) > 0: s += " + "
                s += "%s*%s"%(v[i],x[i])
        if len(s) == 0: s = "0"
        return s


    ###################################
    ## Providing struct. defining items
    ###################################


    def list( self):
        r"""
        Return the cordinates of self w.r.t. the fundamental
        generators of self.parent() as a list.

        EXAMPLES NONE
        """
        return self.__intl_rep

    def c_list( self):
        r"""
        Return the cordinates of self w.r.t. the fundamental
        generators of self.parent() as a list.

        EXAMPLES NONE
        """
        return self.parent()._f2c( self.__intl_rep)

    def _vector_( self):
        r"""
        Return the cordinates of self w.r.t. the fundamental
        generators of self.parent() as a vector.

        EXAMPLES NONE
        """
        return vector( self.list())


    ###################################
    ## Operations
    ###################################

    
    def _add_( self, y):
        r"""
        EXAMPLES NONE
        """
        # Same as _mul_ in FreeAbelianMonoidElement except that the
        # exponents get reduced mod the invariant.
        
        invs = self.parent().elementary_divisors()
        n = len( invs)
        z = FiniteQuadraticModuleElement( self.parent(), 0)
        xelt = self.__intl_rep
        yelt = y.__intl_rep
        zelt = [ xelt[i]+yelt[i] for i in range(len(xelt)) ]
        if len(invs) >= n:
            L =  []
            for i in range(len(xelt)):
                if invs[i]!=0:
                    L.append(zelt[i]%invs[i])
                if invs[i]==0:
                    L.append(zelt[i])
            z.__intl_rep = L
        if len(invs) < n:
            L1 =  []
            for i in range(len(invs)):
                if invs[i]!=0:
                    L1.append(zelt[i]%invs[i])
                if invs[i]==0:
                    L1.append(zelt[i])
            L2 =  [ zelt[i] for i in range(len(invs),len(xelt)) ]
            z.__intl_rep = L1+L2
        return z
            
    #TODO use l_action, r_action and _mule_ for the scalar product
    def __mul__( self, _n):
        r"""
        requires that len(invs) = n
        EXAMPLES NONE
        """
        n = int(_n)
        if n != _n:
            raise TypeError, "Argument n (= %s) must be an integer."%n
        invs = self.parent().elementary_divisors()
        N = len( invs)
        z = FiniteQuadraticModuleElement( self.parent(), 0)
        if n < 0:
            L =[(n*self.__intl_rep[i])%invs[i] for i in range(N)]
            z.__intl_rep = L
            return z
        elif n == 0:
            return z
        elif n == 1:
            return self
        elif n == 2:
            return self + self
        k = n//2
        return self*k + self*(n-k)


    def __neg__( self):
        r"""
        EXAMPLES NONE
        """
        return self.__mul__(-1)


    def __sub__( self, x):
        r"""
        EXAMPLES NONE
        """
        return self._add_( -x)
    

    def __cmp__( self, other,):
        r"""
        EXAMPLES NONE
        """
        if(not isinstance(other,type(self))): return False
        return self.__intl_rep == other.__intl_rep

    def __eq__(self,other):
        r"""
        Test if two FQMElements are equal.
        """
        if not isinstance(other,type(self)): return False
        if not self.parent()==other.parent(): return False
        return self.__intl_rep == other.__intl_rep

    def __ne__(self,other):
        r"""
        Test if two FQMElements are equal.
        """
        return not self.__eq__(other)


    
    ###################################
    ## Associated quantities
    ###################################


    def order( self):
        r"""
        Returns the order of this element.
        
        EXAMPLES
            sage: F.<a,b,c> = FiniteQuadraticModule([2,12,34]); F
            Finite quadratic module in 3 generators:
             gens: a, b, c
             form: 1/8*x0^2 + 1/48*x1^2 + 1/136*x2^2
            sage: x = a - b                                      
            sage: x.order()                                      
            24
        """
        A = self.parent()
        if self == FiniteQuadraticModuleElement( A, 0):
            return Integer(1)
        invs = A.elementary_divisors()
        L = list(self.__intl_rep)
        N = lcm([invs[i]/gcd(invs[i],L[i]) for i in range(len(invs)) if L[i]!=0])
        return N


    def norm( self):
        r"""
        If this element is $a$ and belongs to the module $(M,Q)$ then
        return $Q(a)$.
        
        EXAMPLES NONE
        """
        return self.parent().Q( self)


    def dot( self, b):
        r"""
        If this element is $a$ and belongs a module with associated
        bilinear form $B$ then return $B(a,b)$.
        
        EXAMPLES NONE
        """
        return self.parent().B( self, b)
    


###################################
## CLASS QUAD_MODULE_SUBGROUP
###################################


class FiniteQuadraticModule_subgroup(AbelianGroup):
    r"""
    Descibes a subgroup of the underlying group of a finite quadratic module.

    EXAMPLES NONE
##         sage: p=FiniteQuadraticModule([2,3,10])
##         sage: U = FiniteQuadraticModule_subgroup( [p.0*30])
##         sage: U.is_isotropic()
##         True
##         sage: p1 = U.quotient(); p1
##         Finite quadratic module ([30, 2, 2], 7/30*x0^2 + 1/2*x0*x1 + 1/2*x1^2 + 1/4*x2^2) with generators (e0, e1, e2)
##         sage: V = p1( matrix( ZZ, 3, 2, range(6))); V
##         Subgroup generated by [e0 + e1 + e2] of the Finite quadratic module ([30, 2, 2], 7/30*x0^2 + 1/2*x0*x1 + 1/2*x1^2 + 1/4*x2^2) with generators (e0, e1, e2).
##         sage: W = V.dual(); W
##         Subgroup generated by [15*e0 + e2, e1 + e2] of the Finite quadratic module ([30, 2, 2], 7/30*x0^2 + 1/2*x0*x1 + 1/2*x1^2 + 1/4*x2^2) with generators (e0, e1, e2).
##         sage: V.order(), W.order(), p1.order()
##         (30, 4, 120)
        """

    def __init__( self, gens):
        r"""
        Construct the subgroup generated by the list gens of
        elements of the quadratic module class.

        INPUT
            gens -- nonempty list of elements from the class quadratic module

        EXAMPLES NONE
        
        """
        try:
            x = Sequence( gens)
            ambience = self.__ambient_module = x.universe()
            # Does not work: if not isinstance( ambience, __main__.FiniteQuadraticModule_ambient):
            if not hasattr(ambience,'_is_FiniteQuadraticModule_ambient'):
                raise TypeError
        except TypeError:
            raise TypeError, "%s: must be a list of elements of a quadratic module."%gens[0]
        AbelianGroup.__init__( self)
        # save fun. coordinates of initializing list of gens as columns in the matrix __iMatrix
        self.__iMatrix = matrix( ZZ, [ gens[i].list() for i in range(len(gens))]).transpose()
        # Get lattice generated by the columns of __iMatrix and the relations of the fundamental system.
        self.__lattice = self.__iMatrix.column_module() + diagonal_matrix( ZZ, list(ambience.elementary_divisors())).column_module()
        self.__hnf_gens = [ ambience( list( x), can_coords = False) for x in self.__lattice.matrix()]
        Mat = MatrixSpace( ZZ, len( self.__hnf_gens))
        # TODO: why is this coercion necessary? (SAGE returns rat. matrices for ZZ-lattices)
        self.__hnf_matrix = Mat( self.__lattice.matrix().transpose())
        # throw out 0's, if list becomes empty set to [ambience(0)]
        z = ambience(0); self.__hnf_gens = [ g for g in self.__hnf_gens if g != z ]
        if 0 == len( self.__hnf_gens): self.__hnf_gens = [z]


    ###################################
    ## Iterators
    ###################################


    def __iter__( self):
        r"""
        Return a generator over the elements of self.

        TODO: Smarter implementation(?)
        """
        
        orders = [x.order() for x in self.gens()]
        res = [ sum( self.gens()[i]*x[i] for i in range(len(self.gens()))) for x in xmrange(orders)]
        return (x for x in uniq(res))


    ###################################
    ## Introduce myself
    ###################################


    def _latex_( self):
        r"""
        EXAMPLES NONE
        """
        gens = ', '.join([latex(x) for x in self.gens()])
        return '\\langle %s \\rangle' %gens


    def _repr_( self):
        r"""
        EXAMPLES NONE
        """
        gens = ', '.join([x._repr_() for x in self.gens()])
        return '< %s >' %gens


    ###################################
    ## Providing struct. defining items
    ###################################

    
    def _matrix_( self):
        r"""
        EXAMPLES NONE
        """
        return self.__hnf_matrix


    def ngens( self):
        r"""
        EXAMPLES NONE
        """
        return len( self.__hnf_gens)


    def gen( self, i):
        r"""
        EXAMPLES NONE
        """
        return self.__hnf_gens[i]

    def gens(self):
        r"""
        Return a tuple of generators for self.
        """
        return tuple([ self.gen(i) for i in range(self.ngens())])


    
    def ambience( self):
        r"""
        Return the ambient finite quadratic module.

        EXAMPLES NONE
        """
        return self.__ambient_module


    def _relations( self):
        r"""
        Return the $nxn$-matrix $R$ in lower HNF such that
        $\\Z^n \riightarrow \langle \code{self.ngens()} \rangle$,
        $x \mapsto \code{self.ngens()}\cdot x$ defines an isomorphosm
        of abelian groups.

        EXAMPLES NONE        
        """
        H = matrix( [x.list() for x in self.gens()]).transpose()
        Hp = H.augment( diagonal_matrix( list(self.ambience().elementary_divisors())))
        Rp = matrix(ZZ, Hp.transpose().integer_kernel().matrix().transpose())
        R = Rp.matrix_from_rows( range( self.ngens()))
        return R.transpose().echelon_form().transpose()
    
        
    def as_ambient( self, names = None):
        r"""
        Return a pair $B, f$, where $B$ is a finite quadratic module whose underlying
        group is isomorphic to self and whose quadratic form is
        the one induced by the ambient finite quadratic module $A$,
        and where $f: B \rightarrow A$ is the morphism of finite quadratic modules
        corresponding to the inclusion of self in $A$.

        EXAMPLES
            sage: A = FiniteQuadraticModule( matrix( QQ, 2, [2,1,1,2]))
            sage: B = 2*A                                              
            sage: a,b,c,d = B.gens()                                   
            sage: U = B.subgroup( [b+d])
            sage: C,f = U.as_ambient(); f
            Homomorphism : Finite quadratic module in 1 generators:
             gens: e
             form: 2/3*x^2 --> Finite quadratic module in 4 generators:
             gens: e1, e1, e3, e3
             form: 1/3*x0^2 + 2/3*x0*x1 + 1/3*x1^2 + 1/3*x2^2 + 2/3*x2*x3 + 1/3*x3^2
            sage: f(C.0)                                               
            e1 + e3
        """
        R = self._relations()
        def g(x,y):
            return x.norm() if x == y else x.dot(y)/2
        G = matrix( QQ, self.ngens(), [ g(a,b) for a in self.gens() for b in self.gens()])
        B = FiniteQuadraticModule_ambient( R, G, names = names)
        f = B.hom( self.gens())
        return B, f


    ###################################
    ## Associated quantities
    ###################################

    
    def order( self):
        r"""
        Return the order of this subgroup.
        
        EXAMPLES NONE        

        BUG
        determinant is buggy (workaround: set option 'proof = False'
        """
        return Integer(self.ambience().order()/matrix(self).determinant( proof = False))


    def level( self):
        r"""
        Return the level of this subgroup (viewed as quadratic module w.r.t.
        to the quadratic form inherited from its ambient module).

        EXAMPLES NONE        
        """
        n = self.ngens()
        gens = self.gens()
        v = [ x.norm() for x in gens]
        w = [ gens[i].dot( gens[j]) for i in range(n) for j in range(i+1,n)] 
        return vector( v + w).denominator()


    ###################################
    ## Operations
    ###################################


    def dual( self):
        r"""
        Return the dual subgroup in the ambient module.
        
        EXAMPLES NONE                
        """
        return self.ambience().dual_group( self)


    def __add__( self, other):
        r"""
        Return the sum of this module and the other.
        
        EXAMPLES NONE
        """
        return FiniteQuadraticModule_subgroup( self.gens() + other.gens())


    def cap( self, V):
        r"""
        Return the intersection of this subgroup with the subgroup $V$.

        EXAMPLES
            sage: A.<a,b,c,d> = FiniteQuadraticModule( '3^-4')
            sage: U = A.subgroup( [a+b,c])                    
            sage: V = A.subgroup( [a+b,d])
            sage: W = U.cap(V); W
            < a + b >

            sage: A.<a,b,c,d> = FiniteQuadraticModule( '2_2^-4'); A
            Finite quadratic module in 4 generators:
             gens: a, b, c, d
             form: 1/4*x0^2 + 3/4*x1^2 + 3/4*x2^2 + 3/4*x3^2
            sage: U = A.subgroup( [a,b])
            sage: V = U.dual(); V
            < c, d >
            sage: U.cap(V)              
            < 0 >
        """
        ambience = self.ambience()
        if not ambience is V.ambience():
            raise ValueError
        lat0 = diagonal_matrix( ZZ, list(ambience.elementary_divisors())).column_module()
        lat1 = matrix( self).column_module() + lat0
        lat2 = matrix( V).column_module() + lat0
        lat = lat1.intersection( lat2)
        return ambience.subgroup( [ ambience( list( x), can_coords = False) for x in lat.matrix()])


    def quotient( self):
        r"""
        If this is $U$ and the ambient module $A$ then return $A/U$.
        
        EXAMPLES NONE
        """
        return self.ambience().quotient( self)


    def split( self, n):
        r"""
        Return the splitting $U+V$ of this subgroup
        where the exponent of $U$ divides $n^\infty$ and the exponent of
        $V$ is relatively prime to $n$.

        INPUT
            n -- an integer
            
        NOTE
            $U$ and $V$ are orthogonal to each other.
            If $n$ is a prime then $U$ is the Sylow $p$-subgroup of self.

        EXAMPLES
            sage: A.<a,b,c,d,e,f,g,h,j> = FiniteQuadraticModule( '23^4.2_2^4.3'); A
            Finite quadratic module in 9 generators:
             gens: a, b, c, d, e, f, g, h, j
             form: 1/23*x0^2 + 1/23*x1^2 + 1/23*x2^2 + 1/23*x3^2 + 1/4*x4^2 + 1/4*x5^2 + 1/4*x6^2 + 3/4*x7^2 + 1/3*x8^2
            sage: U = A.subgroup( A.gens()); U
            < a + h + 2*j, b + g, c + f, d + e >
            sage: U2 = U.split(2); U2   
            (< h, g, f, e >, < 2*a + j, 2*b, 2*c, 2*d >)
            sage: U3 = U.split(3); U3
            (< 2*j >, < 3*a + h, b + g, c + f, d + e >)
            sage: U23 = U.split(23); U23
            (< 6*a, 2*b, 2*c, 2*d >, < h + j, g, f, e >)
            sage: V = U2[0] + U3[0] + U23[0]
            sage: U == V
            True
            sage: U is V
            False
        """
        # Let e be the exponent of self, write e=e1e2 with gcd(e1,n)=1, e2 | n^infty,
        # choose x, y such that 1 = e1x+e2y, set u=e1x and v=e2y. If a in A,
        # then a = u*a + v*a is the decomposition of a with respect to self = U + V.
        # In fact, if q dnotes the order of a, let q=q1q2 with gcd(q1,n)=1, q2 | n^infty.
        # Then the order of u*a is q/gcd(q,u) = q/gcd(q,e1x) = q/q1 = q2, and
        # the order of v*a is q/gcd(q,v) = = q/gcd(q,e2y) = q/q2 = q1.
        e = self.level()
        e2 = gcd( e, n); e1 = Integer(e/e2)
        while not 1 == e1.gcd(e2):
            n = n**2; e2 = gcd( e,n); e1 = Integer(e/e2)
        g, x, y = xgcd( e1, e2)
        u = x*e1; v = y*e2
        ul =[]; vl =[]
        for a in self.gens():
            ul.append(u*a); vl.append(v*a)
        return  FiniteQuadraticModule_subgroup(ul), FiniteQuadraticModule_subgroup(vl) 


    ###################################
    ## Predicates
    ###################################


    def is_multiplicative( self):
        r"""
        EXAMPLES NONE        
        """
        return False

    
    def is_isotropic( self):
        r"""
        EXAMPLES NONE
        """
        return 1 == self.level()


    ###################################
    ## Relations
    ###################################
    

    def __lt__( self, other):
        return self.__le( other) and not self == other

    
    def __le__( self, other):
        return self._divides( matrix( other), matrix( self))


    def __eq__( self, other):
        return  matrix( self) == matrix( other)


    def __ne__( self, other):
        return  matrix( self) != matrix( other)
    

    def __gt__( self, other):
        return self.__ge( other) and not self == other


    def __ge__( self, other):
        r"""
        EXAMPLES NONE        
        Test here all of the above
        """
        return self._divides( matrix( self), matrix( other))


    ###################################
    ## Misc
    ###################################


    def orthogonal_basis( self):
        r"""
        Return an orthogonal system for this subgroup
        if this subgroup is nondegenerate
        w.r.t. the scalar product induced from the ambient module.
        Otherwise raise an exception.

        NOTE
            An orthonormal system for a subgroup is a system of pairwise
            orthogonal generators $a_i$ of $p$-order, possibly extended by pairs
            of genertors $b_j,c_j$ which are orthogonal to the $a_i$ and such that
            the subgroups $\langle b_j,c_j \rangle$ are pairwise orthogonal and do no possess
            orthogonal generators.
            
        EXAMPLES:
            ::
            sage: A.<a,b,c,d,e,f,g> = FiniteQuadraticModule( '11^-7')
            sage: U = A.subgroup( [a+b,b+c,c+d,d+e,e+f,f+g]); U
            < a + 10*g, b + g, c + 10*g, d + g, e + 10*g, f + g >
            sage: matrix( len(U.gens()), [ x.dot(y) for x in U.gens() for y in U.gens()])
            [   0 9/11 2/11 9/11 2/11 9/11]
            [9/11 4/11 9/11 2/11 9/11 2/11]
            [2/11 9/11 4/11 9/11 2/11 9/11]
            [9/11 2/11 9/11 4/11 9/11 2/11]
            [2/11 9/11 2/11 9/11 4/11 9/11]
            [9/11 2/11 9/11 2/11 9/11 4/11]
            sage: og_b = U.orthogonal_basis(); og_b
            [b + g,
             b + 9*f + 10*g,
             b + 3*e + f + 10*g,
             b + 7*d + 10*e + f + 10*g,
             b + 5*c + d + 10*e + f + 10*g,
             a + 2*b + 9*c + 2*d + 9*e + 2*f + 9*g]
            sage: matrix( len(og_b), [ x.dot(y) for x in og_b for y in og_b])
            [4/11    0    0    0    0    0]
            [   0 1/11    0    0    0    0]
            [   0    0 2/11    0    0    0]
            [   0    0    0 7/11    0    0]
            [   0    0    0    0 5/11    0]
            [   0    0    0    0    0 2/11]
        """
        if not self.as_ambient()[0].is_nondegenerate():
            raise TypeError
        pl = prime_divisors( self.order())
        og_b = []
        V = self
        while pl != []:
            p = pl.pop()
            if pl != []:
                U,V = V.split(p)
            else:
                U = V
            og_b += U._orthogonal_basis()
        return og_b

    
    def _orthogonal_basis( self):
        r"""
        See FiniteQuadraticModule_subgroup.orthogonal_basis().
        Apply only for nondegenerate $p$-groups.
        Do not use this method directly, you may run in an infinite
        loop for degenerate subgroups. Use orthogonal_basis()
        instead.

        NOTE
            The proof that the algorithm implemented here works is as
            follows (see also [Sko]): Let $U$ be a $p$-subgroup of a
            quadratic module $A$ such that the the scalar product of
            $A$ induces a nondegenerate on on $U$. The value subgroup,
            i.e. the subgroup generated by the values $B(x,y)$ ($x,y$
            in $U$), is (as every subgroup of $\QQQ/\ZZ$) cyclic,
            hence generated by a fixed value $B(x,y)$. We may assume
            that $x$ and $y$ are in $U.gens()$ (since $B(x,y)$ beeing
            a generator is equivalent to the statement that $B(x,y)$
            has the largest denominator among all values).  If $x=y$
            then $U = \langle x \rangle + V$, where $V$ is the dual of
            $\langle x \rangle$ in $U$: in fact, $z - tx$, with an
            integer $t$ such that $B(z,x) = tB(x,x)$, is in $V$ for
            every $z$ in $U$. Assume now that $x != y$ and $B(z,z)$
            never generates the value subgroup for some $z$ in
            $U.gens()$.  Then there are two cases. If $p$ is odd then
            $B(x+y,x+y)$ generates the value subgroup and we can do
            the same argument as before with $x$ replaced by $x+y$.
            However, if $p=2$ then $U = \langle x,y \rangle + V$,
            where $V$ is the orthogonal complement of $\langle x,y
            \rangle$ in $U$. Namely, if $z$ in $U$ then we can find
            integers $t$ and $s$ such that $z - tx -sy$ is in $V$,
            i.e. such that $B(z,x) = t B(x,x) - s B(x,y)$ and $B(z,y)
            = t B(x,y) - s B(y,y)$. Since $U$ is nondegenerate we know
            that the sums here are all direct (since $|X||Y| = |U|$
            for any subgroup $X$ in $U$ with $Y$ denoting its
            dual). But then $V$ is nondegenerate too and we can
            proceed by induuction.
        """
        if 1 == self.order():
            return []
        g = self.gens()
        d = 0
        for x in g:
            e = x.dot(x).denominator()
            if e > d:
                d = e
                og_b = [x]
        for x in g:
            i = list(g).index(x)
            for y in g[i+1:]:
                e = x.dot(y).denominator()
                if e > d:
                    d = e
                    og_b = [x+y] if is_odd(d) else [x,y]
        V = self.cap( self.ambience().subgroup( og_b).dual())
        return og_b + V._orthogonal_basis()


        
    ###################################
    ## Auxiliary functions
    ###################################


    @staticmethod
    def _normal_form( gens):
        r"""
        Assuming that $gens = [x,y]$ and $x$ and $y$ span a type II quadratic module
        return a normalized basis $u$, $v$ for the subgroup generated by $x$ and $y$,
        i.e. generators such that either $x.norm() = x.dot(y) = y.norm()$ or
        $x.norm() = y.norm() = 0$ and $x.dot(y) =1/2^n$.
        """
        M = matrix( [x.norm()])
        pass


    @staticmethod
    def _divides( A, B):
        r"""
        Return  True  or  False  accordingly as  the  rational  matrix
        $A^{-1}*B$ is an integer matrix or not.

        EXAMPLES NONE
        
        """
        ZZ_N_times_N = MatrixSpace(ZZ, A.nrows(), B.ncols())
        return A**(-1)*B in ZZ_N_times_N




###################################
## CLASS JordanDecomposition
###################################


class JordanDecomposition( SageObject):
    r"""
    A container class for the Jordan constituents of a
    finite quadratic module.

    EXAMPLES NONE
    """
    def __init__( self, A):
        r"""
        INPUT
            A -- a finite quadratic module
        """
        self.__A = A
        if not A.is_nondegenerate():
            raise TypeError
        U = A.subgroup( A.gens())
        og_b = U.orthogonal_basis()
        jd = dict()
        ol = []
        primary_comps = uniq( map(lambda x: x.order(), og_b))
        for q in primary_comps:
            basis = tuple( [x for x in og_b if x.order() == q])
            p = q.factor()[0][0]
            n = valuation( q, p)
            r = len( basis)
            def f(x,y):
                return x.norm()*2*q if x == y else x.dot(y)*q
            F = matrix( ZZ, r,r, [ f(x,y) for x in basis for y in basis])
            eps = kronecker( F.det(), p)
            genus = [p, n, r, eps]
            if 2 == p and self.is_type_I(F):
                t = sum(filter(lambda x: is_odd(x), F.diagonal())) % 8
                genus.append( t)
            jd[q] = ( basis, tuple(genus))
            ol.append( (p,n))
        self.__jd = jd
        self.__ol = ol


    def _repr_( self):
        r"""
        EXAMPLES NONE
        """
        return 'Jordan decomposition'


    def _latex_( self):
        r"""
        EXAMPLES NONE
        """
        return 'Jordan decomposition'


    def __iter__( self):
        r"""
        Return the Jordan decomposition as iterator.

        NOTE
            The following is guaranteed. Returned is a list of pairs
            basis, ( prime p,  valuation of p-power n, dimension r, determinant e over p[, oddity o]),
            where $n > 0$, ordered lexicographically by $p$, $n$.

        EXAMPLES NONE
        """
        return ( self.__jd[p**n] for p,n in self.__ol )

        
    def genus_symbol( self, p = None):
        r"""
        Return the genus symbol of the Jordan constituents
        whose exponent is a power of the prime $p$.
        Return the concatenation of all local genus symbols
        if no argument is given.
        
        EXAMPLES:

        We check that the calculation of the genus symbol is correct
        for 2-adic symbols.
            ::
            sage: M = FiniteQuadraticModule('2^-2.2_1^1'); M
            Finite quadratic module in 3 generators:
             gens: e0, e1, e2
             form: 1/2*x0^2 + 1/2*x0*x1 + 1/2*x1^2 + 1/4*x2^2
            sage: M.jordan_decomposition().genus_symbol()
            '2_1^-3'
            sage: N = FiniteQuadraticModule('2_1^-3')
            sage: N.is_isomorphic(M)
            True
            sage: N = FiniteQuadraticModule('2_5^-3')
            sage: N.is_isomorphic(M)
            False
        """
        n = self.__A.order()
        if not p:
            _P = n.prime_divisors()
            _P.sort( reverse = True)
        elif is_prime(p):
            _P = [p] if p.divides(n) else []
        else:
            raise TypeError
        s = ''
        while [] != _P:
            p = _P.pop()
            s +=  self._genus_symbol( p)
            if [] != _P:
                s += '.'
        return s
    
        
    def _genus_symbol( self, p):
        r"""
        Return the genus symbol of the Jordan constituent
        whose exponent is a power of the prime $p$.
        Do not use directly, use genus_symbol() instead.
        
        EXAMPLES NONE
        """
        l = [q for q in self.__jd.keys() if 0 == q%p]
        if [] == l:
            return ''
        l.sort( reverse = True)
        s = ''
        while l != []:
            q = l.pop()
            s += str(q)
            gs = self.__jd[q][1]
            e = gs[2]*gs[3]
            if len(gs) > 4:
                s += '_' + str(gs[4])
            if 1 != e:
                s += '^' + str(e)
            if l != []:
                s += '.'
        return s

        
    def orbit_list(self, p = None, short = False):
        r"""
        If this is the Jordan decomposition for $(M,Q)$, return the dictionary of
        dictionaries of orbits corresponding to the p-groups of $M$.
        If a prime p is given only the dictionary of orbits for the p-group is returned. 
        OUTPUT:
            dictionary -- the mapping p --> (dictionary -- the mapping orbit --> the number of elements in this orbit)

        EXAMPLES:
            sage: A = FiniteQuadraticModule('3^-3.25^2')
            sage: J = JordanDecomposition(A)
            sage: J.orbit_list() == \
                  {3: \
                      {(1,): 1, \
                       (3, 1, 0): 8, \
                       (3, 1, 1/3): 6, \
                       (3, 1, 2/3): 12}, \
                   5: {(1,): 1, \
                       (5, 5, 0): 8, \
                       (5, 5, 1/5): 4, \
                       (5, 5, 2/5): 4, \
                       (5, 5, 3/5): 4, \
                       (5, 5, 4/5): 4, \
                       (25, 1, 1, 0, 0): 40, \
                       (25, 1, 1, 1/25, 1/5): 20, \
                       (25, 1, 1, 2/25, 2/5): 20, \
                       (25, 1, 1, 3/25, 3/5): 20, \
                       (25, 1, 1, 4/25, 4/5): 20, \
                       (25, 1, 1, 1/5, 0): 40, \
                       (25, 1, 1, 6/25, 1/5): 20, \
                       (25, 1, 1, 7/25, 2/5): 20, \
                       (25, 1, 1, 8/25, 3/5): 20, \
                       (25, 1, 1, 9/25, 4/5): 20, \
                       (25, 1, 1, 2/5, 0): 40, \
                       (25, 1, 1, 11/25, 1/5): 20, \
                       (25, 1, 1, 12/25, 2/5): 20, \
                       (25, 1, 1, 13/25, 3/5): 20, \
                       (25, 1, 1, 14/25, 4/5): 20, \
                       (25, 1, 1, 3/5, 0): 40, \
                       (25, 1, 1, 16/25, 1/5): 20, \
                       (25, 1, 1, 17/25, 2/5): 20, \
                       (25, 1, 1, 18/25, 3/5): 20, \
                       (25, 1, 1, 19/25, 4/5): 20, \
                       (25, 1, 1, 4/5, 0): 40, \
                       (25, 1, 1, 21/25, 1/5): 20, \
                       (25, 1, 1, 22/25, 2/5): 20, \
                       (25, 1, 1, 23/25, 3/5): 20, \
                       (25, 1, 1, 24/25, 4/5): 20}}
            True
            sage: A = FiniteQuadraticModule('3^-3.27^2')
            sage: J = JordanDecomposition(A)
            sage: J.orbit_list(3) == \
                  {(1,): 1, \
                  (3, 1, 0): 72, \
                  (3, 1, 1/3): 54, \
                  (3, 1, 2/3): 108, \
                  (3, 9, 1/3): 4, \
                  (3, 9, 2/3): 4, \
                  (9, 1, 3, 0, 1/3): 432, \
                  (9, 1, 3, 0, 2/3): 216, \
                  (9, 1, 3, 1/3, 1/3): 288, \
                  (9, 1, 3, 1/3, 2/3): 432, \
                  (9, 1, 3, 2/3, 1/3): 216, \
                  (9, 1, 3, 2/3, 2/3): 288, \
                  (9, 3, 3, 1/9, 1/3): 12, \
                  (9, 3, 3, 2/9, 2/3): 12, \
                  (9, 3, 3, 4/9, 1/3): 12, \
                  (9, 3, 3, 5/9, 2/3): 12, \
                  (9, 3, 3, 7/9, 1/3): 12, \
                  (9, 3, 3, 8/9, 2/3): 12, \
                  (27, 1, 1, 1, 1/27, 1/9, 1/3): 972, \
                  (27, 1, 1, 1, 2/27, 2/9, 2/3): 972, \
                  (27, 1, 1, 1, 4/27, 4/9, 1/3): 972, \
                  (27, 1, 1, 1, 5/27, 5/9, 2/3): 972, \
                  (27, 1, 1, 1, 7/27, 7/9, 1/3): 972, \
                  (27, 1, 1, 1, 8/27, 8/9, 2/3): 972, \
                  (27, 1, 1, 1, 10/27, 1/9, 1/3): 972, \
                  (27, 1, 1, 1, 11/27, 2/9, 2/3): 972, \
                  (27, 1, 1, 1, 13/27, 4/9, 1/3): 972, \
                  (27, 1, 1, 1, 14/27, 5/9, 2/3): 972, \
                  (27, 1, 1, 1, 16/27, 7/9, 1/3): 972, \
                  (27, 1, 1, 1, 17/27, 8/9, 2/3): 972, \
                  (27, 1, 1, 1, 19/27, 1/9, 1/3): 972, \
                  (27, 1, 1, 1, 20/27, 2/9, 2/3): 972, \
                  (27, 1, 1, 1, 22/27, 4/9, 1/3): 972, \
                  (27, 1, 1, 1, 23/27, 5/9, 2/3): 972, \
                  (27, 1, 1, 1, 25/27, 7/9, 1/3): 972, \
                  (27, 1, 1, 1, 26/27, 8/9, 2/3): 972} 
            True
        """
        n = self.__A.order()
        if not p:
            _P = n.prime_divisors()
            if 2 in _P:
                _P.remove(2)
            _P.sort( reverse = True)
        elif is_prime(p):
            return self._orbit_list(p, short = short)
        else:
            raise TypeError
        orbit_dict = dict()
        while [] != _P:
            p = _P.pop()
            orbit_dict[p] = self._orbit_list(p, short = short)
        return orbit_dict

        
    def _orbit_list(self, p, short = False, debug = 0):
        r"""
        If this is the Jordan decomposition for $(M,Q)$, return the dictionary of
        orbits corresponding to the p-group of $M$.

        OUTPUT:
            dictionary -- the mapping orbit --> the number of elements in this orbit

        EXAMPLES:
            sage: A = FiniteQuadraticModule('3^-3.27^2')
            sage: J = JordanDecomposition(A)
            sage: J._orbit_list(3) == \
                  {(1,): 1, \
                  (3, 1, 0): 72, \
                  (3, 1, 1/3): 54, \
                  (3, 1, 2/3): 108, \
                  (3, 9, 1/3): 4, \
                  (3, 9, 2/3): 4, \
                  (9, 1, 3, 0, 1/3): 432, \
                  (9, 1, 3, 0, 2/3): 216, \
                  (9, 1, 3, 1/3, 1/3): 288, \
                  (9, 1, 3, 1/3, 2/3): 432, \
                  (9, 1, 3, 2/3, 1/3): 216, \
                  (9, 1, 3, 2/3, 2/3): 288, \
                  (9, 3, 3, 1/9, 1/3): 12, \
                  (9, 3, 3, 2/9, 2/3): 12, \
                  (9, 3, 3, 4/9, 1/3): 12, \
                  (9, 3, 3, 5/9, 2/3): 12, \
                  (9, 3, 3, 7/9, 1/3): 12, \
                  (9, 3, 3, 8/9, 2/3): 12, \
                  (27, 1, 1, 1, 1/27, 1/9, 1/3): 972, \
                  (27, 1, 1, 1, 2/27, 2/9, 2/3): 972, \
                  (27, 1, 1, 1, 4/27, 4/9, 1/3): 972, \
                  (27, 1, 1, 1, 5/27, 5/9, 2/3): 972, \
                  (27, 1, 1, 1, 7/27, 7/9, 1/3): 972, \
                  (27, 1, 1, 1, 8/27, 8/9, 2/3): 972, \
                  (27, 1, 1, 1, 10/27, 1/9, 1/3): 972, \
                  (27, 1, 1, 1, 11/27, 2/9, 2/3): 972, \
                  (27, 1, 1, 1, 13/27, 4/9, 1/3): 972, \
                  (27, 1, 1, 1, 14/27, 5/9, 2/3): 972, \
                  (27, 1, 1, 1, 16/27, 7/9, 1/3): 972, \
                  (27, 1, 1, 1, 17/27, 8/9, 2/3): 972, \
                  (27, 1, 1, 1, 19/27, 1/9, 1/3): 972, \
                  (27, 1, 1, 1, 20/27, 2/9, 2/3): 972, \
                  (27, 1, 1, 1, 22/27, 4/9, 1/3): 972, \
                  (27, 1, 1, 1, 23/27, 5/9, 2/3): 972, \
                  (27, 1, 1, 1, 25/27, 7/9, 1/3): 972, \
                  (27, 1, 1, 1, 26/27, 8/9, 2/3): 972} 
            True
        """
        if not is_prime(p):
            raise TypeError
        ppowers = [q for q in self.__jd.keys() if 0 == q % p]
        ppowers.sort()
        if debug > 0: print "ppowers=", ppowers
        if [] == ppowers:
            return dict()
        orbitdict = {(1,) : 1}

        levelpower = ppowers[len(ppowers)-1]

        def _orbit_length( r, eps, t):

            if is_even(r):
                n = Integer(r/2)
                if t.is_integral():
                    return p**(r-1) + eps * kronecker(-1,p)**n * (p**n - p**(n-1)) -1
                else:
                    return p**(r-1) - eps * kronecker(-1,p)**n * p**(n-1)
            else:
                if t.is_integral():
                    return p**(r-1)-1
                else:
                    n = Integer((r-1)/2)
                    return p**(r-1) + eps * kronecker(-1,p)**n * kronecker(2,p) * kronecker(Integer(p*t),p) * p**n

        def _multiplicitieslist():
            r"""
            Create a dictionary of possible combinations of
            reduced multiplicities
            """
            
            completemultiplicitieslist = dict()
            
            for k in range(0, valuation( levelpower, p)):

                m = p**(k+1)
                idash = len([x for x in ppowers if x < m])
                completemultiplicitieslist[k] = []
                usedrankslist = []

                for j in range(idash, len(ppowers)):
                    usedranks = [0 for x in ppowers]
                    usedranks[j] = 1
                    completemultiplicitieslist[k].append([Integer(ppowers[j]/m)])
                    usedrankslist.append(usedranks)

                for j in range(k-1,-1,-1):

                    for x in range(0,len(completemultiplicitieslist[k])):

                        multiplicities = completemultiplicitieslist[k][x]
                        usedranks = usedrankslist[x]
                        idash = len([xx for xx in ppowers if xx <= p**j])

                        completemultiplicitieslist[k][x] = [multiplicities[0]] + multiplicities

                        for jj in [xx for xx in range(idash, len(ppowers)) if usedranks[xx] < self.__jd[ppowers[xx]][1][2] and ppowers[xx] < p**(j+1)*multiplicities[0]]:

                            newusedranks = usedranks[:]
                            newusedranks[jj] += 1
                            newmultiplicities = [Integer(ppowers[jj] / p**(j+1))] + multiplicities
                            completemultiplicitieslist[k].append(newmultiplicities)
                            usedrankslist.append(newusedranks)

            multiplicitieslist = []
            for k in completemultiplicitieslist.keys():
                multiplicitieslist += sorted(completemultiplicitieslist[k])
        
            return multiplicitieslist

        multiplicitieslist =  _multiplicitieslist()

        multiplicitieslist.reverse()

        tconstant = Integer(2)
        while kronecker(tconstant, p) != -1 and tconstant < p:
            tconstant += 1

        ranks = [self.__jd[x][1][2] for x in ppowers]
        epsilons = [self.__jd[x][1][3] for x in ppowers]
        
        while multiplicitieslist != []:

            if debug > 0: print "multiplicitieslist=", multiplicitieslist

            multiplicities = multiplicitieslist.pop()
            k = len(multiplicities)-1
            pk = p**k
            m = p*pk

            if debug > 0: print "pk={0}, m={1}, k={2}".format(pk, m, k)

            if multiplicities[0] == multiplicities[k]:

                ordersDv1 = [Integer(x / multiplicities[0]) for x in ppowers if x > multiplicities[0]]
                ranksDv1 = ranks[len(ppowers) - len(ordersDv1):]
                ordersDv1pk = [Integer(x / pk) for x in ordersDv1 if x > pk]
                ranksDv1pk = ranksDv1[len(ordersDv1)-len(ordersDv1pk):]
                if debug > 0: print "ordersDv1 = {0}, ranksDv1={1}".format(ordersDv1, ranksDv1)
                if debug > 0: print "ordersDv1pk = {0}, ranksDv1pk={1}".format(ordersDv1pk, ranksDv1pk)

                if len(ordersDv1pk) != 0 and ordersDv1pk[0] == p:

                    constantfactor = Integer(prod([min(pk, ordersDv1[j])**ranksDv1[j] for j in range(0, len(ordersDv1))]) / pk)
                    if debug > 0: print "constantfactor=", constantfactor
                    constantfactor *= prod(map(lambda x: p**x, ranksDv1pk[1:]))
                    rank = ranksDv1pk[0]
                    eps = epsilons[len(ranks)-len(ranksDv1pk)]
                    values = [constantfactor * _orbit_length(rank, eps, tconstant / p)]
                    values += [constantfactor * _orbit_length(rank, eps, Integer(0))]
                    values += [constantfactor * _orbit_length(rank, eps, Integer(1) / p)]

                    if short == True:

                        orbitdict[tuple([m] + multiplicities)] = tuple(values)

                    else:

                        nonzeros = [t for t in range(0,p) if values[kronecker(t,p) +1] != 0]

                        for tdash in range(0,m,p):
                            
                            for tdashdash in nonzeros:

                                tt = tdash + tdashdash
                                orbitlength = values[kronecker(tdashdash, p)+1]
                                orbit = tuple([m] + multiplicities + map(lambda x: x - floor(x), [tt * p**j / m for j in range(0, k+1)]))
                                orbitdict[orbit] = orbitlength

            else:

                maxdenominators = [p for x in multiplicities]
                for j in range(k-1,-1,-1):
                    maxdenominators[j] = Integer(max(p*multiplicities[j]/multiplicities[j+1]*maxdenominators[j+1],p))

                skips = [0] + [j+1 for j in range(0,k) if multiplicities[j+1] > multiplicities[j]]
                noskips = [j for j in range(1,k+1) if not j in skips]
                ranklist = []
                epslist = []
                constantfactor = p**(len(skips)-1-k)
                ordersD = [Integer(x / multiplicities[0]) for x in ppowers if x > multiplicities[0]]
                ranksD = ranks[len(ppowers)-len(ordersD):]

                for skip in skips[1:]:

                    quotient = Integer(multiplicities[skip] / multiplicities[skip - 1])
                    ordersA = [x for x in ordersD if x <= quotient * p**skip]
                    ranksA = ranksD[:len(ordersA)]
                    ordersB = ordersD[len(ranksA):]
                    ranksB = ranksD[len(ranksA):]
                    pj = p**(skip - 1)
                    constantfactor *= prod([min(pj, ordersA[j])**ranksA[j] for j in range(0,len(ranksA))])
                    ordersApj = [Integer(x/pj) for x in ordersA if x > pj]
                    ranksApj = ranksA[len(ranksA)-len(ordersApj):]

                    if ordersApj != [] and ordersApj[0] == p:

                        ranklist.append(ranksApj[0])
                        epslist.append(epsilons[len(epsilons)-len(ranksD)])
                        constantfactor *= p**sum(ranksApj[1:])
                        ordersD = [Integer(x / quotient) for x in ordersB if x > quotient]
                        ranksD = ranksB[len(ranksB)-len(ordersD):]

                    else:

                        constantfactor = 0
                        break

                if constantfactor:

                    constantfactor *= prod([min(pk, ordersD[j])**ranksD[j] for j in range(0,len(ordersD))])
                    ordersDpk = [Integer(x / pk) for x in ordersD if x > pk]
                    ranksDpk = ranksD[len(ranksD)-len(ordersDpk):]

                    if ordersDpk != [] and ordersDpk[0] == p:

                        ranklist.append(ranksDpk[0])
                        epslist.append(epsilons[len(epsilons)-len(ranksDpk)])
                        constantfactor *= p**sum(ranksDpk[1:])

                    else:

                        constantfactor = 0

                if constantfactor:

                    product = [constantfactor] + [0 for j in skips[2:]]
                    skipindex = 0
                    tpjvalues = [0 for j in skips]
                    tpjs = [[x / maxdenominators[0] for x in range(0,maxdenominators[0])]] + [[] for j in skips[1:]]

                    # if debug > 0: print "tpjs", tpjs
                    
                    while tpjs[0] != []:

                        if tpjs[skipindex] == []:

                            skipindex -= 1
                            tpjs[skipindex].remove(tpjvalues[skipindex])

                        else:

                            if skipindex == len(skips)-1:

                                for tpj in tpjs[skipindex]:

                                    tpjvalues[skipindex] = tpj
                                    tpk = p**(k - skips[skipindex]) * tpj
                                    factor = product[len(product) - 1]
                                    t = p**(skips[skipindex] - skips[skipindex - 1] - 1) * tpjvalues[skipindex - 1]
                                    t -= multiplicities[skips[skipindex]] / multiplicities[skips[skipindex] - 1] / p * tpj
                                    orbitlength1 = _orbit_length(ranklist[skipindex - 1], epslist[skipindex - 1], t)
                                    orbitlength2 = _orbit_length(ranklist[skipindex], epslist[skipindex], tpk)
                                    orbitlengths = orbitlength1 * orbitlength2

                                    if orbitlengths != 0:

                                        reducednorms = [0 for j in range(0,k+1)]
                                        for j in range(0,len(skips)):
                                            reducednorms[skips[j]] = tpjvalues[j]
                                        for j in noskips:
                                            t = p*reducednorms[j-1]
                                            reducednorms[j] = t - floor(t)

                                        orbitdict[tuple([m] + multiplicities + reducednorms)] = factor * orbitlengths

                                skipindex -= 1
                                tpjs[skipindex].remove(tpjvalues[skipindex])

                            else:

                                tpjvalues[skipindex] = tpjs[skipindex][0]

                                if skipindex > 0:
                                    
                                    t = p**(skips[skipindex] - skips[skipindex - 1] - 1) * tpjvalues[skipindex - 1]
                                    t -= multiplicities[skips[skipindex]] / multiplicities[skips[skipindex] - 1] / p * tpjvalues[skipindex]
                                    
                                    orbitlength = _orbit_length(ranklist[skipindex - 1], epslist[skipindex - 1], t)

                                    if orbitlength == 0:

                                        tpjs[skipindex].remove(tpjvalues[skipindex])

                                    else:
                                        maxdenominator = maxdenominators[skips[skipindex + 1]]
                                        tcandidates = [t/maxdenominator for t in range(0,maxdenominator)]
                                        difference = skips[skipindex + 1] - skips[skipindex]
                                        quotient = multiplicities[skips[skipindex + 1]] / multiplicities[skips[skipindex]]
                                        tpjs[skipindex + 1] = [t for t in tcandidates if (quotient * t - p**difference * tpjvalues[skipindex]).is_integral()]
                                        product[skipindex] = orbitlength * product[skipindex - 1]
                                        skipindex += 1

                                else:

                                    maxdenominator = maxdenominators[skips[skipindex + 1]]
                                    tcandidates = [t/maxdenominator for t in range(0,maxdenominator)]
                                    difference = skips[skipindex + 1] - skips[skipindex]
                                    quotient = multiplicities[skips[skipindex + 1]] / multiplicities[skips[skipindex]]
                                    tpjs[skipindex + 1] = [t for t in tcandidates if (quotient * t - p**difference * tpjvalues[skipindex]).is_integral()]

                                    skipindex += 1
                                    
        return orbitdict

    def values( self, debug = 0):
        r"""
        If this is the Jordan decomposition for $(M,Q)$, return the values of $Q(x)$ ($x \in M$) as a dictionary d.

        OUTPUT:
            dictionary -- the mapping Q(x) --> the number of elements x with the same value Q(x)

        EXAMPLES:
            sage: A = FiniteQuadraticModule('3^-3.27^2')
            sage: J = JordanDecomposition(A)
            sage: J.values() == \
                  {0: 729, \
                  1/27: 972, \
                  2/27: 972, \
                  4/27: 972, \
                  5/27: 972, \
                  7/27: 972, \
                  8/27: 972, \
                  1/3: 810, \
                  10/27: 972, \
                  11/27: 972, \
                  13/27: 972, \
                  14/27: 972, \
                  16/27: 972, \
                  17/27: 972, \
                  2/3: 648, \
                  19/27: 972, \
                  20/27: 972, \
                  22/27: 972, \
                  23/27: 972, \
                  25/27: 972, \
                  26/27: 972}
            True
            sage: A = FiniteQuadraticModule('2_3^-3.4^2.8_2^-2')
            sage: J = JordanDecomposition(A)
            sage: J.values() == \
                  {0: 480, \
                  1/8: 512, \
                  3/16: 1024, \
                  1/4: 480, \
                  3/8: 512, \
                  7/16: 1024, \
                  1/2: 544, \
                  5/8: 512, \
                  11/16: 1024, \
                  3/4: 544, \
                  7/8: 512, \
                  15/16: 1024}
            True
        """
        n = self.__A.order()

        values = [1]
        
        def combine_lists( list1, list2):
            N1 = len(list1)
            N2 = len(list2)
            newlength = lcm(N1,N2)
            newlist = [0 for j1 in range(0,newlength)]
            for j1 in range(0,N1):
                n1 = list1[j1]
                if n1 != 0:
                    for j2 in range(0,N2):
                        n2 = list2[j2]
                        if n2 != 0:
                            newlist[(j1 * Integer(newlength / N1) + j2 * Integer(newlength / N2)) % newlength] += (n1 * n2)
            return newlist        

        def values_even2adic(gs):
            p, l, n, eps = gs
            n /= 2
            factor = 2**((l-1)*n)
            if n == 1 and eps == 1:
                return [factor * (l+2)] + [factor * (valuation(j,2)+1) for j in range(1,2**l)]
            else:
                quotient = Integer(2**n-eps)/Integer(2**(n-1)-eps)
                return [factor * (quotient * (2**((n-1)*(l+1))-eps**(l+1)) + eps**(l+1))] + [factor * quotient * (2**((n-1)*(l+1)) - eps**(valuation(j,2)+1) * 2**((n-1)*(l-valuation(j,2)))) for j in range(1,2**l)] 
            
        def values_odd2adic(gs):
            p, l, n, eps, t = gs
            # t = t%8
            if eps == +1:
                if t == 0:
                    tvalues = (1,7)
                elif t == 1:
                    tvalues = (1,)
                elif t == 2:
                    tvalues = (1,1)
                elif t == 3:
                    tvalues = (1,1,1)
                elif t == 4:
                    tvalues = (1,1,1,1)
                elif t == 5:
                    tvalues = (7,7,7)
                elif t == 6:
                    tvalues = (7,7)
                elif t == 7:
                    tvalues = (7,)
                else:
                    raise TypeError
            elif eps == -1:
                if t == 0:
                    tvalues = (5,1,1,1)
                elif t == 1:
                    tvalues = (3,7,7)
                elif t == 2:
                    tvalues = (3,7)
                elif t == 3:
                    tvalues = (3,)
                elif t == 4:
                    tvalues = (3,1)
                elif t == 5:
                    tvalues = (5,)
                elif t == 6:
                    tvalues = (5,1)
                elif t == 7:
                    tvalues = (5,1,1)
            else:
                raise TypeError

            level = 2**(l+1)
            squarerepresentationlist = [0 for j in range(0,level)]
            if is_even(l+1):
                squarerepresentationlist[0] = squarerepresentationlist[2**(l-1)] = 2**((l+1)/2)
            else:
                squarerepresentationlist[0] = squarerepresentationlist[2**l] = 2**(l/2)
            for k in range(0,l-1,2):
                # if debug > 0: print "k:", k
                for a in range(1, 2**(l+1-k),8):
                    # if debug > 0: print "a:", a
                    squarerepresentationlist[2**k * a] = 2**(k/2 + 2)
            # if debug > 0: print "Test the squarelist:", sum(squarerepresentationlist) == level, squarerepresentationlist, level

            # if debug > 0: print "tvalues", tvalues
            
            t1inverse = inverse_mod(tvalues[0], level)
            values = [squarerepresentationlist[(j * t1inverse)%level]/2 for j in range(0,level)]

            if len(tvalues) > 1: #The following works only for tvalues where the last elements coincide

                t2inverse = inverse_mod(tvalues[1], level)
                newvalues = [squarerepresentationlist[(j * t2inverse)%level]/2 for j in range(0,level)]

                for j in range(1,len(tvalues)):
                    values = combine_lists(values, newvalues)

            neven = n - len(tvalues)
            if neven > 0:
                values = combine_lists(values, values_even2adic((p, l, neven, +1)))
            
            return values

                
        _P = n.prime_divisors()
        if 2 in _P:

            _P.remove(2)
            
            l = sorted([q for q in self.__jd.keys() if 0 == q%2])
            if debug > 0: print l
            while l != []:
                q = l.pop()
                gs = self.__jd[q][1]
                if len(gs) > 4:
                    values = combine_lists(values, values_odd2adic(gs))
                else:
                    values = combine_lists(values, values_even2adic(gs))
                if debug > 0: print values

        _P.sort( reverse = True)

        while [] != _P:
            p = _P.pop()
            if debug > 0: print "p = ", p
            shortorbitdict = self.orbit_list(p, short = True)
            level = max(q for q in self.__jd.keys() if 0 == q%p)
            newvalues = [0 for j in range(0,level)]
            newvalues[0] = 1
            
            for orbit in shortorbitdict.keys():

                if orbit != (1,):
                    
                    k = Integer(valuation(orbit[0],p)-1)
                    # if debug > 0: print orbit
                    v1 = orbit[1]
                    if v1 == orbit[k+1]:
                    
                        orbitlengthsbykronecker = shortorbitdict[orbit]
                        for t1 in range(0,p**(k+1)):
                            newvalues[Integer(v1 * t1 * level / p**(k+1)) % level] += orbitlengthsbykronecker[kronecker(t1,p) + 1]

                    else:

                        newvalues[Integer(v1 * orbit[k+2] * level) % level] += shortorbitdict[orbit]

                    # if debug > 0: print "Position1:", sum(newvalues)
            # if debug > 0: print "1:", values
            # if debug > 0: print "2:", newvalues
            values = combine_lists(values, newvalues)
            # if debug > 0: print "3:", values
            # if debug > 0: print "Position2:", values, _P, p

        # if debug > 0: print "Position3:", values
            
        valuesdict = {Integer(j)/len(values) : values[j] for j in range(0,len(values)) if values[j] != 0}

        # if debug > 0: print "Position4:", values, valuesdict, "END"
            
        return valuesdict #, values

    def two_torsion_values( self):
        r"""
        If this is the Jordan decomposition for $(M,Q)$, return the values of $Q(x)$
        for $x \in M_2=\{x \in M | 2*x = 0\}$, the subgroup of two-torsion elements as a dictionary.

        OUTPUT:
            dictionary -- the mapping Q(x) --> the number two-torsion elements x with the same value Q(x)

        EXAMPLES:
            sage: A = FiniteQuadraticModule('2_3^-3.4^2.8_2^-2')
            sage: J = JordanDecomposition(A)
            sage: J.two_torsion_values() == {0: 48, 1/4: 16, 1/2: 16, 3/4: 48}
            True
        """
        n = self.__A.order()

        values = [1]
        
        def combine_lists( list1, list2):
            N1 = len(list1)
            N2 = len(list2)
            newlength = lcm(N1,N2)
            newlist = [0 for j1 in range(0,newlength)]
            for j1 in range(0,N1):
                n1 = list1[j1]
                if n1 != 0:
                    for j2 in range(0,N2):
                        n2 = list2[j2]
                        if n2 != 0:
                            newlist[(j1 * Integer(newlength / N1) + j2 * Integer(newlength / N2)) % newlength] += (n1 * n2)
            return newlist        

        def two_torsion_values_even2adic(gs):
            p, l, n, eps = gs
            n /= 2
            fourn = 4**n
            if l == 1:
                epstwon = eps * 2**n
                return [(fourn + epstwon)/2, (fourn - epstwon)/2]
            else:
                return [fourn]
                
        def two_torsion_values_odd2adic(gs):
            p, l, n, eps, t = gs
            if l == 1:
                # print "n:", n, "eps:", eps, "t:", t, "n-t:", n-t, (n-t)/2
                if eps == -1:
                    t = (t + 4) % 8
                n2 = ((n - t)/2) % 4
                n1 = n - n2
                # print "n1:", n1, "n2:", n2, "t:", t
                list1 = [sum([binomial(n1,k) for k in range(j,n1+1,4)]) for j in range(0,4)]
                # print list1
                if n2 == 0:
                    return list1
                else:
                    list2 = [[1,0,0,1],[1,0,1,2],[1,1,3,3]][n2 - 1]
                    # print list2
                    return combine_lists(list1, list2)
            elif l == 2:
                twonminusone = 2**(n-1)
                return [twonminusone, twonminusone]
            else:
                return [2**n]
            
        l = sorted([q for q in self.__jd.keys() if 0 == q%2])
        while l != []:
            q = l.pop()
            gs = self.__jd[q][1]
            if len(gs) > 4:
                values = combine_lists(values, two_torsion_values_odd2adic(gs))
            else:
                values = combine_lists(values, two_torsion_values_even2adic(gs))
            
        valuesdict = {Integer(j)/len(values) : values[j] for j in range(0,len(values)) if values[j] != 0}

        return valuesdict #, values

        
    def constituent( self, q, names = None):
        r"""
        Return the Jordan constituent whose exponent is the
        prime power "q".
        
        EXAMPLES NONE
        """
        if not is_prime_power(q):
            raise TypeError
        gens = self.__jd.get( q, ((),()))[0]
        return self.__A.spawn( gens, names)


    def finite_quadratic_module(self):
        r"""
        Return the finite quadratic module who initialized
        this Jordan decomposition.
        
        EXAMPLES NONE
        """
        return self.__A


    def basis( p =None):
        r"""
        """
        pass

    
    @staticmethod
    def is_type_I(F):
        r"""
        EXAMPLES NONE
        """
        for i in range( F.nrows()):
            if is_odd( F[i,i]):
                return True
        return False

    
    
r"""
o: For elements return tuples instad of lists. Are their immutable matrices?

o A.quotient(U) should return B, f, g where  B is V/U (V = dual of U), f:V-->A and g:V-->B are
  the natural mophims.

o Implementation details should appear as NOTES not as comments

o Trivial madule + A should return A.

o If U is of type I the U.orthogonal_basis() should really give an OG basis

o Type II should return hyp's or hyp's + a single non hyp

o Less confusing for use: if a basis is already orthogonal the return it.

o The og basis should correspond to a 'normalized' genus_symbol

BUGS:

"""

###  Testing routines ###

def FiniteQuadraticModuleRandom(discbound=100,normbound=100,verbose=0):
    """
    Returns a random finite quaratic module with discriminant within the discriminant bound.
    
    """
    D = 0
    while False == is_fundamental_discriminant(D):
        D = ZZ.random_element( -discbound, discbound)
    K=QuadraticField(D,names='a')
    a = K.gens()[0]
    beta = 0
    while 0 == beta:
        beta = K.random_element()
    alpha = 0
    while 0 == alpha:
        alpha = K.random_element()
    om = beta/alpha
    if verbose>0:
        print "D=",D
        print "K=",K
        print "om=",om
    if om.denominator_ideal().absolute_norm()>normbound:
        return FiniteQuadraticModuleRandom(discbound,normbound,verbose)
    #A.<x,y> =
    A = FiniteQuadraticModule( om/om.parent().gen())
    x,y = A.gens()
    N = A.kernel()
    if verbose>0:
        print "A=",A
        print "|A.list()|=",len(list(A))
        print "N=",N
    if False == N.is_isotropic():
        return FiniteQuadraticModuleRandom(discbound,normbound,verbose)
    if len(list(A)) == 1:
        return FiniteQuadraticModuleRandom(discbound,normbound,verbose)
    if  max( map(max,A.gram().rows()))==0 and min( map(min,A.gram().rows()))==0:
        return FiniteQuadraticModuleRandom(discbound,normbound,verbose)
#if False == A.is_nondegenerate():
    #    return FiniteQuadraticModuleRandom(bound)
    B = A.quotient(N)
    if False == B.is_nondegenerate():
        if verbose>0:
            print "A / N is nondegenerate!"
            print "A=",A
            print "Gram(A)=",A.gram()
            print "B=",B
            print "N=A.kernel()=",N
        raise ValueError,"Quotient is nondegenerate!"
    return B



def test_fqm_random(fqbound=100,nbound=10,cbound=10,size_bd=50,verbose=0):
    r"""
    Check nbound different random finite quadratic modules

    The tests we do are:
    1) Compute Jordan block
    2) Compute sigma invariant
    """
    # first get a random FQM
    ntest=0
    for i in range(nbound):
        l=size_bd+1        
        while l>size_bd:
            FQ=FiniteQuadraticModuleRandom(fqbound,nbound,verbose-1)
            l=len(list(FQ))
            #info=get_factors2(FQ.jordan_decomposition())
        t = test_one_F(FQ)
        if t<>True:
            return t
        ntest+=1
        #assert s0 == si[0]/si[1]
    if verbose>0:
        print "Tested {0} modules!".format(ntest)
    return True

def test_one_F(FQ='4_1'):
    if not hasattr(FQ,"Q"):
        FQ = FiniteQuadraticModule(FQ)
    N = FQ.level()
    z = CyclotomicField(N).gens()[0]
    for a in range(1,N):
        s0 = FQ.char_invariant(a)
        s1 = naive_Gauss_sum(FQ,a)
        #print s0,s1
        if abs(CC(s0[0])*CC(s0[1])-CC(s1[0])/CC(s1[1]**2))>1e-10:
                if verbose>0:
                    print 'a=', a
                    print "s0=",s0,CC(s0[0]*s0[1])
                    print "s1=",s1,CC(s1[0]/s1[1]**2)
                return False,a,FQ
    return True


def test_fqm_from_signature(num_tests=10,num_comp=4,prime_bd=5,pow_bd=3,verbose=0):
    r"""
    Check that the genus symbol determines the finite quadratic module.
    """
    from sage.all import prime_range,random,random_prime
    for n in range(num_tests):
        s=""
        for i in range(num_comp):
            p = random_prime(prime_bd)
            k = ZZ.random_element(1,pow_bd+1)
            q = p**k
            if p == 2:
                if random()<0.5: ## Make a type 1 factor
                    t = 2*ZZ.random_element(4)+1
                    if kronecker(t,2) == 1:
                        s+="{q}_{t}^1".format(t=t,q=q)
                    else:
                        s+="{q}_{t}^-1".format(t=t,q=q)
                else:
                    if random()<0.5:
                        s+="{q}^2".format(q=q) # B
                    else:
                        s+="{q}^-2".format(q=q) # C
            else:
                if random()<0.5:
                    s+="{q}^1".format(q=q)
                else:
                    s+="{q}^-1".format(q=q)                    
            if i<num_comp-1:
                s+="."
        M = FiniteQuadraticModule(s)
        s0 = M.jordan_decomposition().genus_symbol()
        N = FiniteQuadraticModule(s0)
        if verbose>0:
            print "s=",s
            print "s0=",s0
        if verbose>1:
            print "M=",M
            print "N=",N            
            
        if not M.is_isomorphic(N):
            raise ArithmeticError,"{0} and {1} with symbol {2} are not isomorphic!".format(M,N,s)
        # if(FQ.level() % 4 <> 0 or is_odd(info['sign'])):
        #     continue
        # if verbose>0:
        #     print "FQ=",FQ," size=",l," level=",FQ.level()," sigma^4=",FQ.sigma_invariant()^4
        #     print "Genus symbol=",FQ.jordan_decomposition().genus_symbol()
        
        # if test_epd_one_fqm(FQ,cbound):
        #     print "ok for FQ=",i
        # else:
        #     print "*Not* ok for FQ=",i
        #     return False
    ## Functions for testing
def naive_Gauss_sum(FQ,a):
    r"""
    If this quadratic module equals $A = (M,Q)$, return
    the characteristic function of $A$ at $s$, i.e. return
    $$\chi_A (s)= |M|^{-1}\sum_{x\in M} \exp(2\pi i s Q(x)))$$
    computed naively by summing.
    NOTE: This is slow and should only be used for testing purposes.
    """
    N = FQ.level()
    z = CyclotomicField(N).gens()[0]
    res = 0
    for x in list(FQ):
        res += z**(a*(FQ.Q(x)*N))
    return res,ZZ(len(list(FQ))).sqrt()


def orbitlist_test(str = None):
    r"""
    testing if all orbits sum up to the size of the
    finite quadratic module
    """
    if str:
        A = FiniteQuadraticModule(str)
    else:
        A = FiniteQuadraticModuleRandom(discbound=10000,normbound=10000,verbose=0)
    J = JordanDecomposition(A)
    str = J.genus_symbol()
    print str
    olist = J.orbit_list()
    testpassed = True
    for p in olist.keys():
        st = J.genus_symbol(p)
        order = p**valuation(A.order(),p)
        orbitsum = sum(olist[p][key] for key in olist[p].keys())
        # print "   order of A:", order
        # print "sum of orbits:", orbitsum
        testpassed = testpassed and (order == orbitsum)
        print p.str() + ": # of elements in the computed orbits sum up to the order of " + st + ":", order == orbitsum, order, orbitsum
    return testpassed
    
def values_test(str):
    r"""
    testing if the computed dictionary of values sums up
    to the size of the finite quadratic module
    """
    A = FiniteQuadraticModule(str)
    J = JordanDecomposition(A)

    print str

    # valuesdict, values = J.values()
    valuesdict = J.values()
    
    # print "Position1"

    Avalues = A.values()
    b1 = valuesdict == Avalues
    print "Test A.values() == J.values():", b1
    
    # print "Position2"
    
    b2 = sum([valuesdict[key] for key in valuesdict.keys()]) == A.order()
    print "Test sum(values) == A.order():", b2
    
    # print "Position3"

    Atwotorsionvalues = A.kernel_subgroup(2).as_ambient()[0].values()
    Jtwotorsionvalues = J.two_torsion_values()

    b3 = Atwotorsionvalues == Jtwotorsionvalues
    print "Test two_torsion_values():    ", b3

    # print Avalues
    # print valuesdict    

    if not b3:
        print "A:", Atwotorsionvalues
        print "J:", Jtwotorsionvalues

    return b1 and b2 and b3

    
def testing_routine(p):
    r"""
    testing discriminant only with p components
    """
    k = Integer(3)
    p = Integer(p)
    p1 = p.str()
    p2 = (p**2).str()
    p3 = (p**3).str()
    p4 = (p**4).str()
    p5 = (p**5).str()
    str = ''
    for a in range(-k,k+1):
        if a != 0:
            astr = Integer(a).str()
            str1 = str +  p1+'^'+astr+'.'
        else:
            str1 = str
        print "str1:", str1
        for b in range(-k+1,k):
            if b != 0:
                bstr = Integer(b).str()
                str2 = str1 + p2+'^'+bstr+'.'
            else:
                str2 = str1
            print "    str2:", str2
            for c in range(-k+3,k-2):
                if c != 0:
                    cstr = Integer(c).str()
                    str3 = str2 + p3+'^'+cstr+'.'
                else:
                    str3 = str2
                print "        str3:", str3
                for d in range(0,1): # range(-k+4,k-3)
                    if d != 0:
                        dstr = Integer(d).str()
                        str4 = str3 + p4+'^'+dstr+'.'
                    else:
                        str4 = str3
                    print "            str4:", str4
                    for e in range(0,1): # range(-k+4,k-3)
                        if e != 0:
                            estr = Integer(e).str()
                            str5 = str4 + p5+'^'+estr+'.'
                        else:
                            str5 = str4
                        print "                str5:", str5
                        str5 = str5[:-1]
                        if str5 != '':
                            # A = FiniteQuadraticModule(str5)
                            # J = JordanDecomposition(A)
                            if values_test(str5):
                                print str, True
                            else:
                                return str, False
    return True, "All tests successfull"
                                
def testing_routine_odd2adic():

    for p in [3,5,7,9,11,13,17,19,25]:

        q = Integer(2)
        
        while q < 2**4:
            
            for oddstr in ['_0^2', '_0^-4', '_1^1', '_1^-3', '_2^2', '_2^-2', '_3^3',
                           '_3^-1', '_4^4', '_4^-2', '_5^3', '_5^-1', '_6^2', '_6^-2',
                           '_7^1', '_7^-3', '_0^4', '_0^-6', '_1^3', '_1^-5', '_2^4',
                           '_2^-4', '_3^5', '_3^-3', '_4^6', '_4^-4', '_5^5', '_5^-3',
                           '_6^4', '_6^-4', '_7^3', '_7^-5']:

                oddprimestr = '.' + Integer(p).str() + '^' + (-1 + 2 * floor(2 * random())).str() + '.27^-1'
                if not values_test(q.str() + oddstr + oddprimestr):

                    return "Test not passed:", q.str() + oddstr + oddprimestr
            
            q *= 2
