# -*- coding: iso-8859-1 -*-
#*****************************************************************************
#       Copyright (C) 2009 Nils-Peter Skoruppa <nils.skoruppa@uni-siegen.de>
#                          Fredrik Stroemberg <stroemberg@mathematik.tu-darmstadt.de>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

r"""
Weil representations for finite quadratic modules.
Implements a class for working with Weil representations and also the explicit formula for general elements of SL(2,Z).


Note: We did not want to implement the metaplectic group so the Weil representation for odd signature is given by the canonical section rho(A)=rho(A,j_A(z))  where j_A(z)=(cz+d)**-1/2 if A = ( a b // c d )
This means that the representation will not be multiplicative, but
$\rho(A)\rho(B)=\sigma(A,B)\rho(AB)$ where the cocycle $\sigma(A,B)$ is implemented as the function sigma_cocycle().





REFERENCES:
 - [St] Fredrik Strömberg, "Weil representations associated to finite quadratic modules", arXiv:1108.0202
 


 AUTHORS:
   - Fredrik Strömberg
   - Nils-Peter Skoruppa
   - Stephan Ehlen


EXAMPLES::


    sage: F=FiniteQuadraticModule('5^1')


   
"""

#from sage.all_cmdline import *   # import sage library

from sage.all import Integer,RR,CC,QQ,ZZ,sgn,cached_method,copy,CyclotomicField,lcm,is_square,matrix,SL2Z,MatrixSpace,floor,ceil,is_odd,is_even,hilbert_symbol,sqrt,inverse_mod
from sage.structure.formal_sum import FormalSums
from sage.structure.formal_sum import FormalSum
from psage.modform.maass.mysubgroups_alg import factor_matrix_in_sl2z
from sage.modules.vector_integer_dense import Vector_integer_dense
from weil_module_alg import *
from finite_quadratic_module import FiniteQuadraticModuleElement,FiniteQuadraticModule

#_sage_const_3 = Integer(3);
#_sage_const_2 = Integer(2);
#_sage_const_1 = Integer(1);
#_sage_const_0 = Integer(0);
#_sage_const_4 = Integer(4);
#_sage_const_8 = Integer(8);
#_sage_const_10000 = Integer(10000)

###################################
## CLASS WEIL_MODULE
###################################


class WeilModule (FormalSums):
    r"""
    Implements the Weil representation of the metaplectic
    cover $Mp(2,Z)$ or $SL(2,Z)$ of associated to a finite
    quadratic module $A$.
    More pecisely, it implements the $K$-vector space $K[A]$
    as $Mp(2,Z)$-module, where $K$ is the $lcm(l,8)$-th cyclotomic field
    if $l$ denotes the level of $A$.
    """


    def __init__(self, A,**kwds):
        r"""
        Initialize the Weil representation associated to the
        nondegenerate finite quadratic Module A.

        INPUT
            A -- an instance of class  FiniteQuadraticModule_ambient.

        """
        #    raise ValueError, "the finite quadratic module is degenerate"
        if not hasattr(A,"is_nondegenerate"):
            ## If we don't have a finite quadratic module we try to creat one from A
            if isinstance(A,(int,Integer)):
                A = [A]
            A = FiniteQuadraticModule(A)
        try:
            if(A.is_nondegenerate()):
               self._QM = copy(A)
               #self._QM = A
            else:
                raise TypeError," is degenerate!"
        except AttributeError:
            raise TypeError," Argument must be a nondegenerate Finite Quadratic module. Got \n %s" %(A)
        except TypeError:
            raise TypeError," Argument must be a nondegenerate Finite Quadratic module. Got *degenerate* module:\n %s" %(A)
        # Recall that we have elements in CyclotomicField (A.level())
        # times the sigma invariant
        self._verbose = 0 
        self._K = CyclotomicField( lcm(A.level(),8 ))
        FormalSums.__init__( self, base_ring = self._K)        
        #l = copy(self._QM.list())
        #self._list=[] ## Ordered list of elements, in the way that -list[i]=list[D-i]
        #for i in range(len(l)):
        #    self._list.append() 
        # for i in range(len(l)):
        #     x = l[i]
        #     try:
        #         yi = l.index(-x)
        #     except ValueError:
        #         raise ArithmeticError,"Could not find {0} in {1}!".format(-self._L[i],selt._L)
        #     y = l[yi]    
        #     self._neg_ix.append(self._L.index(-self._L[i]))
        
      
        # Silly variable used because
        # issinstance(W,WeilModule)
        # is not reliable... 
        self._is_WeilModule=True

        ### Get the generators Figure out which field we need for the Weil matrices
        #self.__minimal_rho_gen=self.__K.gen()
        # note that the generator of the new R is zeta_ro=zeta_mgm**nm

        self._zl = self._K(CyclotomicField(self._QM.level()).gen())  # e(1/level)
        self._z8 = CyclotomicField(8).gen()                   # e(1/8)
                
        self._n=self._QM.order()
        self._sqn = self._n.sqrt()
        self._gens=self._QM.gens()
        self._gen_orders=[]
        for jj,g in enumerate(self._gens):
            self._gen_orders.append(Integer(g.order()))
        #self._minus_element=[]
        #self._minus_element=self._get_negative_indices()

        self._level= self._QM.level()
        self._xcs  = self._compute_all_xcs()
        # Pre-compute invariants
        #print "self=",type(self)
        self._inv = self._get_invariants(self)
        self._zero = WeilModuleElement([(0,self._QM(0))],parent=self)
        self._even_submodule = None
        self._odd_submodule = None
        self._even_submodule_ix = None
        self._odd_submodule_ix = None
        self._basis = []
        self._dim_param={}
        self._dim_param_numeric={}
        #self._L = list()
        #e=WeilModuleElement(self,self._QM.list()[0])
        #for ii in range(0,self._n):
        #    self._L.append(WeilModuleElement(self,e))

    def basis(self):
        r"""
        Gives a basis of self as a vector space of dimension |D|
        """
        if self._basis==[]:
            for x in list(self._QM):
                self._basis.append(WeilModuleElement([(1,x)],parent=self))
        return self._basis
        
    def _get_negative_indices(self):
        return cython_neg_indices(self._n,self._gen_orders)
        
    def _get_negative_indices_python(self):
        l=[]
        for ii in range(self._n):
            l.append(self._neg_index_python(ii))
        return l

    #@cached_method
    def _el_index(self,c):        
        if not isinstance(c,(list,Vector_integer_dense)):
            raise ValueError,"Need element of list form! Got c={0} of type={1}".format(c,type(c))
        if not len(c)==len(self._gen_orders):
            raise ValueError,"Need element of list of the same length as orders={0}! Got c={1}".format(self._gen_orders,c)
        return cython_el_index(c,self._gen_orders)
    
    @cached_method
    def _neg_index(self,ii):
        return cython_neg_index(ii,self._gen_orders)

    def neg_index(self,x):
        r""" Return the index of minus x in the list of basis elements in self."""
        return self._neg_index(x)
    
    @cached_method
    def _elt(self,ii):
        return cython_elt(ii,self._gen_orders)

    def zero(self):
        return self._zero

    def an_element(self):
        return WeilModuleElement(self._QM.an_element(),parent=self)

    def random_element(self):
        return WeilModuleElement(self._QM.an_element(),parent=self)
        
    ###################################
    ## Introduce myself ...
    ###################################


    def _latex_( self):
        r"""
        EXAMPLES
        """
        return 'Weil module associated to %s' % latex(self._QM)


    def _repr_(self):
        r"""
        EXAMPLES
        """
        return "Weil module associated to %s" % self._QM
 

    ###################################
    ## Coercion
    ###################################


    def __call__( self, x):
        r"""
        Coerce object into an appopriate child object
        of self if possible.

        We coerce
        - an element of $A$ into an elemnt of $K[A]$.

        EXAMPLES
        """
        #print "x=",x
        #print "par=",x.parent()
        if isinstance( x, FiniteQuadraticModuleElement):
            if x.parent() is self._QM:
                return WeilModuleElement([(1 ,x)],parent=self)

        # otherwise we assume that x is a list of pairs (a,k),
        # where a is in self.__QM and k is in self.__K
        #print "x=",x
        return  WeilModuleElement(x,parent=self)

        raise TypeError, "argument must be an element of %s" %self._QM


    ###################################
    ## Basic Methods ...
    ###################################

    def matrix(self, A,filter=None,by_factoring=False,**kwds):
        r"""
        Return the matrix rho(A) giving the action of self on $K[A]$.

        INPUT:
        
        -''A''
        -''filter'' -- integer matrix of the same size as the Weil representation. Describes which elements to compute. Only relevant if using formula
        -''by_factoring'' -- logical (default False). If set to True we use the definition of the Weil representation together with factoring the matrix A.  
        - 'kwds' set for example:

            -'numeric' -- integer. set to 1 to return a matrix_complex_dense

        OUTPUT:
        --''[r,f]'' -- r =  n x n matrix over self._K describing
                       f =  sqrt(len(Dc))/sqrt(len(D))
                       r*f = rho(A)

        EXAMPLES::

        
        sage: F = FiniteQuadraticModule([3,3],[0,1/3,2/3])
        sage: W = WeilModule(F)
        sage: A = SL2Z([1,2,1,3])
        sage: [r,f]=W.matrix(A)
        """
        prec = kwds.get('prec',0)
        if prec > 0:
            return  weil_rep_matrix_mpc(self,A[0,0],A[0,1],A[1,0],A[1,1],filter=filter,prec=prec,verbose=self._verbose)
        # We only need the diagonal elements of rho(A)
        n=len(list(self._QM))
        # Need a WeilModuleElement to compute the matrix
        e=WeilModuleElement(self._QM.gens()[0 ],self,verbose=self._verbose)
        #print "e=",e
        if(by_factoring==False):
            if filter<>None:
                [r,fac]=e._action_of_SL2Z_formula(A,filter,**kwds)
            else:
                [r,fac]=e._action_of_SL2Z_formula(A,**kwds)
        else:
            [r,fac]=e._action_of_SL2Z_factor(A)
        return [r,fac]


    def trace(self, A):
        r"""
        Return the trace of the matrix A in Mp(2,Z) or SL(2,Z)
        as endomorphism of $K[A]$.

        EXAMPLES
        F = FiniteQuadraticModule([3,3],[0,1/3,2/3])
        W = WeilModule(F)
        A = SL2Z([1,2,1,3])
        """
        # We only need the diagonal elements of rho(A)
        n=len(list(self._QM))
        filter=MatrixSpace(ZZ,n).identity_matrix()
        # Need a WeilModuleElement to compute the matrix
        e=self._zero
        [r,fac]=e._action_of_SL2Z_formula(A,filter)
        s=0
        for j in range(0 ,n):
            s=s+r[j,j]
        #print "s=",s,fac
        return [s,fac]


    def trace_all(self):
        r"""
        Return a list of "all" traces
        I.e. W.trace([a,b,c,d]) where 0<=c<=Q.level() and d mod c
        """
        l=list()
        for c in range(self._QM.level()):
            for d in range(c):
                # we only want d=0 if c=1
                if(d==0  and c<>1 ): 
                    continue
                elif(d==0  and c==1 ):
                    a=0 ; b=1 
                else:
                    [g,b,a]=xgcd(c,d)
                    if(g<>1 ):
                        continue
                [t,f]=self.trace([a,-b,c,d])
                #print "a,b,c,d:trace=",a,-b,c,d,':',t
                l.append(t)
        return [l,f]

    def _compute_all_xcs(self):
        r"""
        Computes all non-zero values of the element x_c in the group D
        OUPUT: dictionsry k => x_c  where 2^k || c
        This routine is called at initialization of WeilModuleElement
        """
        J=self._QM.jordan_decomposition()
        res=dict()
        res[0 ]=0 
        for c in J:
            xc=0 
            p,k,r,d = c[1 ][:4 ]
            t = None if 4  == len(c[1 ]) else c[1 ][4 ]
            if(p<>2  or t==None): 
                continue
            q=2 **k # = 2^n
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
        return res

    @staticmethod
    # Compute total invariants of a Jordan decomposition
    def _get_invariants(self):
        r"""
        INPUT: JD is a Jordan decomposition returned from the FQM package
        OUPUT: dictionary res with entries:
        res['total p-excess']= sum of all p-excesses
        res['total oddity']  = sum of all oddities
        res['oddity']        = list of oddities of all components
        res['p-excess']=     = list of p-excesses of all components
        res['signature']     = signature := (odt-pexc_tot) % 8
        res['type']          = dictionary q=>'I' or 'II' telling whether the 2-adic component q is odd or even
        """
        JD=self._QM.jordan_decomposition()
        pexc=dict()
        odts=dict()
        res=dict()
        types=dict()
        sigma_inv=dict()
        sigma_inv[0]=self._QM.sigma_invariant()
        odt=0 
        for c in JD:
            [p,n,r,ep]=c[1 ][:4 ]
            t = 0  if (len(c[1 ])==4 ) else c[1 ][4 ]
            q=p**n
            if( (not is_square(q)) and (ep==-1 )):
                k=1 
            else:
                k=0 
            if p!=2:
                pexc[q]=(r*(q-1 )+4 *k) % 8 
            else:
                odts[q]=((t+4 *k) % 8 )
                if t<>0:
                    types[q]="I" # odd
                else:
                    types[q]="II" # even
            sigma_inv[p]=self._QM.sigma_invariant(p)
        pexc_tot=0 
        odt=0 
        for q in pexc.keys():
            pexc_tot=pexc_tot+pexc[q]
        for q in odts.keys():
            odt=odt+odts[q]        
        res['total p-excess']=pexc_tot % 8 
        res['total oddity']=odt % 8 
        res['oddity']=odts
        res['p-excess']=pexc
        res['signature']=(odt-pexc_tot) % 8 
        res['type']=types
        res['sigma']=sigma_inv
        return res

    def signature(self):
        return self._inv['signature']

    def invariant(self,s):
        if self._inv.has_key(s):
            return  self._inv[s]
        raise ValueError,"Invariant {0} is not defined! Got:{1}".format(s,self._inv.keys())

    def finite_quadratic_module(self):
        return self._QM
    
    def rank(self):
        return self._n
    
    def signature(self):
        return self._inv['signature']

    def oddity(self):
        return self._inv['oddity']

    def level(self):
        return self._level
    
    def pexcess(self,p):
        return self._inv['p-excess'].get(p,0)

    def odd_submodule(self,indices=0):
        if indices==0:
            return self._symmetric_submodule(sign=-1)
        else:
            return self._symmetric_submodule_ix(sign=-1)

    def even_submodule(self,indices=0):
        if indices==0:
            return self._symmetric_submodule(sign=1)
        else:
            return self._symmetric_submodule_ix(sign=1)
        
    def _symmetric_submodule(self,sign=1):
        r"""
        Compute the submodule of self which is C-spanned by
        e_{gamma} + sign*e_{-gamma} for gamma in self._QM
        """
        if sign==1 and self._even_submodule<>None:
            return self._even_submodule
        if sign==-1 and self._odd_submodule<>None:
            return self._odd_submodule
        if sign==1:
            basis = [self(self.finite_quadratic_module()(0))]            
        else:
            basis = []
        elts=[self.finite_quadratic_module()(0)]            
        for x in self.finite_quadratic_module():
            if -x in elts or x in elts: continue
            elts.append(x)
            w1 = self(x); w2=self(-x)
            if sign==1:
                if w1==w2:
                    f = w1
                else:
                    f = w1+w2
            else:
                f = w1-w2
            #print "f0=",f,type(f),f==self.zero()
            if f == self.zero(): continue
            #print "f1=",f,type(f),f==self.zero()
            if f not in basis:
                basis.append(f)
        if sign==1:
            self._even_submodule = basis
        else:
            self._odd_submodule = basis
        return basis

    def _symmetric_submodule_ix(self,sign=1):
        r"""
        Compute the submodule of self which is C-spanned by
        e_{gamma} + sign*e_{-gamma} for gamma in self._QM
        """
        if sign==1 and self._even_submodule_ix<>None:
            return self._even_submodule_ix
        if sign==-1 and self._odd_submodule_ix<>None:
            return self._odd_submodule_ix
        if sign==1:
            basis = [0]
        else:
            basis = []
        el   = list(self.finite_quadratic_module())
        elts = [el[0]]  
        for x in self.finite_quadratic_module():
            if -x in elts or x in elts: continue
            elts.append(x)
            w1 = self(x); w2=self(-x)
            if sign==1:
                if w1==w2:
                    f = w1
                else:
                    f = w1+w2
            else:
                f = w1-w2
            #print "f0=",f,type(f),f==self.zero()
            if f == self.zero(): continue
            #print "f1=",f,type(f),f==self.zero()
            fi = el.index(x)
            if fi not in basis:
                basis.append(fi)
        if sign==1:
            self._even_submodule_ix = basis
        else:
            self._odd_submodule_ix = basis
        return basis


    def dim_parameters(self):
        r"""
        Compute and store parameters which are used in the dimension formula.
        """
        if self._dim_param<>{}:
            return self._dim_param
        res = {}
        F = self.finite_quadratic_module()
        K = CyclotomicField(24)
        s2,d2 = F.char_invariant(2)
        s1,d1 = F.char_invariant(-1)
        s3,d3 = F.char_invariant(3)
        #d1 = K(d1**2).sqrt()
        d2 = K(d2**2).sqrt()
        d3 = K(d3**2).sqrt()
        sq3 = K(3).sqrt()
        sqd = K(self.rank()).sqrt()        
        #if CC(d1).real()<0: d1 = -d1
        if CC(d2).real()<0: d2 = -d2
        if CC(d3).real()<0: d3 = -d3
        if CC(sq3).real()<0: sq3 = -sq3
        if CC(sqd).real()<0: sqd = -sqd
        #s1 = s1
        s2 = s2 #*d2
        s3 = s3 #*d3
        ## The factors which are not contained in the ground field are added later
        res['f1']=1
        res['f2']=1
        res['f3']=1
        if d2 in K:
            s2 = s2*d2
        else:
            res['f2']=res['f2']*d2
        if d3 in K:
            s3 = s3*d3
        else:
            res['f3']=res['f3']*d3
        if sq3 in K:
            s1 = s1/sq3
            s3 = s3/sq3
        else:
            res['f1']=res['f1']/sq3
            res['f3']=res['f3']/sq3
        if sqd in K:
            s2 = s2*sqd
            s3 = s3*sqd
        else:
            res['f2']=res['f2']*sqd
            res['f3']=res['f3']*sqd
        res['sq3']=sq3
        res['s1']=s1 #/sq3
        res['s2']=s2 #*sqd
        res['s3']=s3 # *sqd/sq3

        
        for ep in [-1,1]:
            res[ep]={}
            W = self._symmetric_submodule(ep)
            dim = len(W)
            res[ep]['dim']=dim
            Qv = []; Qv2=[]
            for x in W:
                qx = F.Q(x[0][1]) ## We only take one representative
                if qx>=1 or qx<0:
                    qx = qx - floor(qx)
                qx2 = 1 - qx   ## For the dual rep.
                if qx2>=1 or qx2<0:
                    qx2 = qx2 - floor(qx2)
                Qv.append(qx)
                Qv2.append(qx2)
            K0 = Qv.count(0)
            res[ep]['K0'] =  K0
            res[ep]['Qv']={1:Qv,-1:Qv2}
            res[ep]['parabolic'] =  {1:K0 + sum(Qv),-1:K0 + sum(Qv2)}
        self._dim_param = res
        return res

    def dim_parameters_numeric(self):
        r"""
        Compute and store parameters which are used in the dimension formula.
        """
        if self._dim_param_numeric<>{}:
            return self._dim_param_numeric
        res = {}
        F = self.finite_quadratic_module()
        #K = CyclotomicField(24)
        s2,d2 = F.char_invariant(2)
        s1,d1 = F.char_invariant(-1)
        s3,d3 = F.char_invariant(3)
        s1 = CC(s1)
        s2 = CC(s2)
        s3 = CC(s3)
        #d1 = K(d1**2).sqrt()
        d2 = RR(d2); d3 = RR(d3); sq3 = RR(3).sqrt()
        sqd = RR(self.rank()).sqrt()        
        #s2 = s2 #*d2
        #s3 = s3 #*d3

        res['f1']=1;  res['f2']=1;  res['f3']=1
        res['s1']=s1/sq3
        res['s2']=s2*d2*sqd
        res['s3']=s3*d3/sq3*sqd
        res['sq3']=sq3

        for ep in [-1,1]:
            res[ep]={}
            W = self._symmetric_submodule(ep)
            dim = len(W)
            res[ep]['dim']=dim
            Qv = []; Qv2=[]
            for x in W:
                qx = F.Q(x[0][1]) ## We only take one representative
                if qx>=1 or qx<0:
                    qx = qx - floor(qx)
                qx2 = 1 - qx   ## For the dual rep.
                if qx2>=1 or qx2<0:
                    qx2 = qx2 - floor(qx2)
                Qv.append(qx)
                Qv2.append(qx2)
            K0 = Qv.count(0)
            res[ep]['Qv']={1:Qv,-1:Qv2}
            res[ep]['K0'] =  K0
            res[ep]['parabolic'] =  {1:K0 + sum(Qv),-1:K0 + sum(Qv2)}
        self._dim_param_numeric = res
        return res
            
    def dimension_mod_forms(self,k,sgn=0,verbose=0,numeric=False):
        d,ep = self.dimension_cusp_forms(k=k,sgn=sgn,verbose=verbose,numeric=numeric)
        if numeric:
            K0 = self.dim_parameters_numeric()[ep]['K0']
        else:
            K0 = self.dim_parameters()[ep]['K0']
        return d+K0,ep
        
    def dimension_cusp_forms(self,k,sgn=0,verbose=0,numeric=False):
        r"""
        Compute the dimension of cusp forms weight k with representation self (if sgn=1) or the dual of self (if sgn=-1)
        If numeric = False we work with algebraic numbers.

        INPUT:

        - `k`      -- integer. the weight
        - `sgn`     -- integer. sgn=1 for the Weil representation and -1 for the dual
        - `verbose` -- integer
        - `numeric` -- bool. set to True to use floating-point numbers (much faster) instead of algebraic numbers.

        OUTPUT:

        -`d,ep` -- tuple. dimension of S^{ep}(rho,k) where rho is the Weil representation or its dual and ep is the symmetry of the space.

        """
        if k<=2:
            raise ValueError,"Only weight k>2 implemented! Got k={0}".format(k)
        s = (sgn*self.signature()) % 8
        t = (2*k-s) % 4
        if verbose>0:
            print "k = {0}, (2k-{1}) % 4 = {2}".format(k,s,t)
        ep = 0
        if t == 0: ep = 1
        if t == 2: ep = -1
        if ep==0: return 0,0
        s = self.signature()
        if numeric:
            par = self.dim_parameters_numeric()
        else:
            par = self.dim_parameters()
        
        s1 = par['s1']; s2 = par['s2'];  s3 = par['s3']
        if sgn==-1:
            s1 = s1.conjugate(); s2 = s2.conjugate(); s3 = s3.conjugate(); s=-s
        F = self.finite_quadratic_module()
        if not numeric:
            K = CyclotomicField(24)
            z24 = K.gens()[0]
            z8 = z24**3
            identity_term = QQ(par[ep]['dim']*(k+5))/QQ(12)
            six = QQ(6); eight=QQ(8)
        else:
            twopi_24=RR.pi()/RR(12)
            z24=CC(0,twopi_24).exp()
            twopi_8=RR.pi()/RR(4)
            z8=CC(0,twopi_8).exp()
            identity_term = RR(par[ep]['dim']*(k+5))/RR(12)
            six = RR(6); eight=RR(8)
        arg = (2*k-s) % 8
        elliptic_term_e2 = (s2*ep+s2.conjugate())*z8**arg/eight
        elliptic_term_e2 = elliptic_term_e2 * par['f2'] 
        arg = (4*k+2-3*s) % 24
        if verbose>1:
            print "z8=",z8,CC(z8)
            print "s1=",CC(s1)
            print "s2=",CC(s2)
            print "s3=",CC(s3)
        zf = z24**arg
        e31 = s1*zf
        e32 = s3*zf
        e311 = e31 + e31.conjugate()
        e322 = e32 + e32.conjugate()
        e311*=par['f1']
        e322*=par['f3']
        e3 = (e311+ep*e322)/six
        #e3 = (s1+ep*s3) * zf        
        elliptic_term_e3 = e3 #(e3+e3.conjugate())/6

        parabolic_term = - par[ep]['parabolic'][sgn]
        if verbose>0:
            print "ep=",ep
            print "sgn=",sgn
            print "signature=",s
            print "Id = ",identity_term
            print "E(2)=",elliptic_term_e2,CC(elliptic_term_e2)
            print "E(3)=",elliptic_term_e3,CC(elliptic_term_e3)
            print "P=",parabolic_term
        res = identity_term + elliptic_term_e2 + elliptic_term_e3 + parabolic_term
        ## We want to simplify as much as possible but "res.simplify is too slow"
        if hasattr(res,"coefficients") and not numeric:
            co = res.coefficients()
            if len(co)>1 or co[0][1]<>0:
                raise ArithmeticError,"Got non-real dimension: d={0} = {1}".format(res,res.simplify())
            return co[0][0],ep
        else:
            if not numeric:
                return res,ep
            if abs(res.imag())>1e-10:
                raise ArithmeticError,"Got non-real dimension: {0}".format(res)
            d1 = ceil(res.real()-1e-5); d2 = floor(res.real()+1e-5)
            if d1<>d2:
                raise ArithmeticError,"Got non-integral dimension: d={0}".format(res)
            return d1,ep



        
    
        
###################################
## CLASS WEIL_MODULE_ELEMENT
###################################



class WeilModuleElement (FormalSum):
    r"""
ss    Describes an element of a Weil module $K[A]$.
    """

    def __init__(self,d,parent=None, check = True,verbose=0):
        r"""
        INPUT
            W -- a Weil module $K[A]$
            d -- a dictionary of pairs $a:k$, where $a$ is an element of
                 a finite quadratic module $A$ of level $l$ and $k$ an element
                 of the field $K$ of the $l$-th roots of unity.

        EXAMPLES:

        
           F = FiniteQuadraticModule([3,3],[0,1/3,2/3])
           W = WeilModule(F)
           a,b = F.gens()
           z = W(a+5*b)
           
        """
        ## If d is just a single group element we make a list
        #print "isinst=",isinstance(parent,(WeilModule,WeilModuleElement))
        #if(not isinstance(parent,WeilModule) and ):
        self._coordinates = []
        self._verbose = verbose
        if parent==None:
            parent = WeilModule(d.parent())
        if not hasattr(parent,'_is_WeilModule'):
            raise TypeError, "Call as WeilModuleElement(W,d) where W=WeilModule. Got W=%s" % parent
        if not isinstance(d,list):            
            if hasattr(d.parent(),"_is_FiniteQuadraticModule_ambient"):
                # print "Is instance!"
                x = d
                d=[(1 ,x)]
            elif d==0:
                d=[(1,0)]
            else:
                s="d={0} Is not instance! Is:{1}".format(d,type(d))
                raise TypeError,s
        elif isinstance(d[0],tuple):
            try:
                par = d[0 ][1 ].parent()
                if(not parent._QM):
                    raise TypeError, "Need element of self._QM! got %s" % (d[0 ][1 ].parent())
            except :
                raise  TypeError, "Need element of self._QM! to construct WeilModuleElement! got %s" % (d)
        elif len(d)==parent.rank():
            ## d should be a list of coordinates.
            dd = []
            self._coordinates = d
            for i in range(len(d)):
                dd.append((d[i],list(parent.finite_quadratic_module())[i]))
            d = dd
        else:
            raise ValueError,"Could not interpret {0} as a Weil module element!".format(d)
        ## Init with a formal sum
        
        FormalSum.__init__(self,d, parent, check, True)
        if self._coordinates == []:
            l = list(parent.finite_quadratic_module())
            self._coordinates = [0 for i in range(parent.rank())]
            for i,x in self._data:
                ix = l.index(x)
                self._coordinates[ix]=i

        self._W=parent
        self._parent = parent
        # BaseField including both FQM.level() roots and eight-roots of unity
        self._K  = self._W._K                                          # Base field
        self._QM  = parent._QM                                         # f.q.m.
#ifdef NEW
        self._zl = parent._zl  # e(1/level)
        self._z8 = parent._z8                   # e(1/8)
        self._n  = parent._n                                       
        self._sqn = parent._sqn
                          #    self._minus_element.append(self._QM.list().index(-self._QM.list()[i]))
#endif /* not NEW */
        self._level= self._QM.level()
        self._xcs  = parent._xcs
        # Pre-compute invariants
        #print "self=",type(self)
        self._inv = parent._inv
        #self._B={}  ## stores a 
        # TODO: some checking here and/or study sage.structure.formal_sum.py

        
    ###################################
    ## Operations
    ###################################


    def __repr__(self):
        s=""
        i=0
        for j,x in self:
            if j==1 and i==0:
                s+="{0}".format(x)
            elif j==1 and i>0:
                s+=" + {0}".format(x)
            elif j==-1 and i==0:
                s+="-({0})".format(x)
            elif j==-1 and i>0:
                s+=" - ({0})".format(x)
            elif i==0:
                s+="{0}*({1})".format(j,x)
            else:
                s+=" + {0}*({1})".format(j,x)
            i+=1
        return s




        
    def parent(self):
        return self._parent

    def coordinates(self):
        return self._coordinates
    
    @cached_method
    def Bi(self,i,j):
#ifdef NEW
        return self._QM.B(self._QM(list(self._W._elt(i))),self._QM(list(self._W._elt(j))))
#else /* not NEW */
        #return self._QM.B(self._L[i],self._L[j])
#endif /* not NEW */

    def B(self,other):
        v1 = self._coordinates
        v2 = other._coordinates
        s = 0
        for i in len(v1):
            for j in len(v2):
                s+=self._QM.B(list(self._QM)[i],list(self._QM)[j])
        return s



    @cached_method
    def Q(self,i):
        #return self._QM.Q(self._L[i])
        return self._QM.Q(self._QM(list(self._W._elt(i))))

    @cached_method
    def _minus_element(self,ii):
        return self._W._neg_index(ii)

    def action( self, A):
        r"""
        Return the result of applying the element $A$ of $SL2Z$
        to this element.

        EXAMPLES
        """
        if A not in SL2Z:
            raise TypeError, "%s must be an element of SL2Z" %A

        #if A[1,0] % self._level == 0:
        #    return self._action_of_Gamma0(A)
        #else:
        return self.action_of_SL2Z(A)
        
        #def __rmul__(self,obj):
        #return self.__mul__( obj)

    ## Note that we want matrices from SL2Z to act from the left
    ## i.e. with an _rmul_
    
    def _rmul_(self,obj,act='l'):
        r"""
        Return the result of applying obj to this element from the left. I.e
        INPUT: obj
        OUPUT =
              obj*self = a*self  if obj=a in CyclotomicField(lcm(level,8))
              obj*self = rho_{Q}(A)*self if obj= A in SL2Z
        TODO: Implement the orthogonal group of $A$
        EXAMPLES

        """
        if self._K.has_coerce_map_from(obj.parent()):
            #print "case 1"
            d = list()
            for (j,x) in self:
                #print "(j,x)=",j,x
                d.append( (j*obj,x) )
            return WeilModuleElement(d,self.parent())
        # multiply with the group element instead
        if self._QM.base_ring().has_coerce_map_from(obj.parent()) :
            #print "case 2"
            d = list()
            for (j,x) in self:
                d.append((j,obj*x))
            return WeilModuleElement(self.parent(),d)

        if obj in SL2Z:
            return self.action_of_SL2Z(obj,act)
        # #         # TODO: obj can also be an element of the orthogonal group of $A$
        raise TypeError, "action of %s on %s not defined" %(obj, self)
    

    # def __eq__(self,other):
    #     r"""
    #     Check if self is equal to other
    #     TODO: check reduced words
    #     """
    #     if not hasattr(other,"parent"): return False
    #     if self.parent()<>other.parent(): return False
    #     v1 = self._coordinates
    #     v2 = other._coordinates
    #     if 
    #     return self._data == other._data
    def __add__(self,other):
        r"""
        Add self and other. We use this since FormalSum does not reduce properly
        """
        if other.parent()<>self.parent(): raise ValueError,"Can not add {0} to {1}".format(self,other)
        v1 = self._coordinates
        v2 = other._coordinates
        v_new = [v1[i]+v2[i] for i in range(len(v1))]
        return  WeilModuleElement(v_new,self.parent())

    def __sub__(self,other):
        r"""
        Add self and other. We use this since FormalSum does not reduce properly
        """
        if other.parent()<>self.parent(): raise ValueError,"Can not add {0} to {1}".format(self,other)
        v1 = self._coordinates
        v2 = other._coordinates
        v_new = [v1[i]-v2[i] for i in range(len(v1))]
        return  WeilModuleElement(v_new,self.parent())


    def __mul__( self, obj):
        r"""
        Return the result of applying obj to this element from the right. I.e
        INPUT: obj
        OUPUT =
              obj*self = a*self  if obj=a in CyclotomicField(lcm(level,8))

        NOTE:
        We have only implemented elements in SL2Z to act from the left        
        TODO: Implement the orthogonal group of $A$ and transposed PSL2Z action

        EXAMPLES
        
        """
        return self._rmul_( obj,act='r')
        #if obj.parent() == self.parent().base_ring():
        # Multiply the coefficient of x in the FormalSum


    def action_of_SL2Z(self,M,act='l',mode=0 ):
        r"""
          INPUT:
            M = element of SL2Z
            mode : determines how the Weil rpresentation is calculated
              0 => use formula (default)
              1 => use factorization
            act : determines if we act from the right or left
              'l' => rho(A)*x  
              'r' => x^t * rho(A)
          OUTPUT: [s,fac]
             where s is a Formal Sum over  
             K=CyclotomicField(lcm(8,level) and 
             fact is a scalar and such that
               rho_{Q}(M)*self = fact*s
             where 
          EXAMPLE:
            FQ = FiniteQuadraticModule
            A  = Element of SL2Z
            W = WeilModule(FQ)
            x=FQ.gens()[0]
            WE=W(x)
            [r,s] = x.action_of_SL2Z(A)
            fact*x = rho_{Q}(M) 
        """
        # Since we only want to act on this specific element we first
        # figure out which matrix-entries we need to compute:
        if mode==0:
            filter=matrix(ZZ,self._n)
            for (k,x) in self:
                jj = self._parent._el_index(x.list())
                #jj = self._L.index(x)
                for ii in range(0 ,self._n):
                    if(act=='l'):
                        filter[ii,jj]=1 
                    else:
                        filter[jj,ii]=1 
            [r,fact] = self._action_of_SL2Z_formula(M,filter)
        else:
            # when we use factorization we must compute all elements
            [r,fact] = self._action_of_SL2Z_factor(M)
        # Compute rho_Q(A)*self
        res = FormalSum([(0 ,0 )],self._W)
        for (k,x) in self:
            #print "k,x=",k,x
            jj = self._parent._el_index(x.list())
            #jj = self._L.index(x)
            for ii in range(0,self._n):
                if(act=='l'):
                    res=res+FormalSum([(r[ii,jj],self._QM(self._W._elt(ii)))],self._W)
                    #res=res+FormalSum([(r[ii,jj],self._L[ii])],self._W)
                else:
                    res=res+FormalSum([(r[jj,ii],self._QM(self._W._elt(ii)))],self._W)
                    #res=res+FormalSum([(r[jj,ii],self._L[ii])],self._W)
        return [res,fact]


    ## Action by special (simple) elements of SL(2,Z)   
    def _action_of_T(self,b=1,sign=1,filter=None):
        r""" Action by the generator sign*T^pow=[[a,b],[0,d]]
        where a=d=sign
        """
        r = matrix(self._K,self._n)
        if sign==-1:
            si=self._QM.sigma_invariant()**2 
        else:
            si=1
        for ii in range(0 ,self._n):
            if filter<>None and filter[ii,ii]<>1:
                continue
            if sign==1:
                r[ii,ii] = self._zl**(ZZ(b)*self._level*self.Q(ii))
            else:
                #r[self._n-1-ii,ii] = self._zl**(b*self._level*self._QM.Q(self._L[ii]))
                r[self._minus_element(ii),ii] = si*self._zl**(b*self._level*self.Q(ii))
        return [r,1]

    def _action_of_S(self,filter=None,sign=1,mult_by_fact=False):
        r"""
        Action by the generator S=[[0,-1],[1,0]]
        """
        r = matrix(self._K,self._n)
        if sign==-1:
            si = self._K(self._QM.sigma_invariant()**3)
            if is_odd(self.parent().signature()):
                si = -si # sigma(Z,A)
        else:
            si = self._K(self._QM.sigma_invariant())
        for ii in range(0 ,self._n):
            for jj in range(0 ,self._n):
                if filter<>None and filter[ii,jj]<>1:
                    continue
                arg = -sign*self._level*self.Bi(ii,jj)
                #arg = -self._level*self._QM.B(self._L[ii],self._L[jj])
                r[ii,jj] = si*self._zl**arg
        #r = r*
        fac = self._parent._sqn**-1
        return [r,fac]

    def _action_of_STn(self,pow=1,sign=1,filter=None):
        r""" Action by  ST^pow or -ST^pow
        NOTE: we do not divide by |D|
        """
        ## Have to find a basefield that also contains the sigma invariant
        if pow==0:
            return self._action_of_S(filter,sign)
        r  = matrix(self._K,self._n)
        if sign==-1:
            si = self._K(self._QM.sigma_invariant()**3)
            if is_odd(self.parent().signature()):
                si = -si  # sigma(Z,A)
        else:
            si = self._K(self._QM.sigma_invariant())
        for ii in range(self._n):
            for jj in range(self._n):
                argl=self._level*(pow*self.Q(jj)-sign*self.Bi(ii,jj))
                #ii = self._L.index(x); jj= self._L.index(j)
                if filter<>None and filter[ii,jj]<>1:
                    continue
                r[ii,jj] = si*self._zl**argl
        fac = self._parent._sqn**-1
        return [r,fac]

    def _action_of_Z(self,filter=None):
        r""" Action by  Z=-Id
        NOTE: we do not divide by |D|
        """
        ## Have to find a basefield that also contains the sigma invariant
        r = matrix(self._K,self._n)
        for ii in range(0 ,self._n):
            if filter<>None and filter[ii,ii]<>1:
                continue
            jj=self._W._neg_index(ii)
            r[ii,jj]=1 
        r = r*self._QM.sigma_invariant()**2 
        return [r,1]
    
    def _action_of_Id(self,filter=None):
        r""" Action by  Z=-Id
        NOTE: we do not divide by |D|
        """
        ## Have to find a basefield that also contains the sigma invariant
        r = matrix(self._K,self._n)
        for ii in range(0 ,self._n):
            if filter<>None and filter[ii,ii]<>1:
                continue
            r[ii,ii]=1 
        #r = r*self._QM.sigma_invariant()**2 
        return [r,1]


    # Action by Gamma_0(N) through formula
    def _action_of_Gamma0(self,A,filter=None):
        r"""
        Action by A in Gamma_0(l) 
        where l is the level of the FQM
        INPUT:
           A in SL2Z with A[1,0] == 0 mod l
           act ='r' or 'l' : do we act from left or right'
        filter = |D|x|D| integer matrix with entries 0 or 1
                         where 1 means that we compute this entry
                         of the matrix rho_{Q}(A) 
        """
        #print "A_in_acton_of=",A
        a=A[0 ,0 ]; b=A[0 ,1 ]; c=A[1 ,0 ]; d=A[1 ,1 ] #[a,b,c,d]=elts(A)
        if(c % self._level <>0 ):
            raise ValueError, "Must be called with Gamma0(l) matrix! not A=" %(A)
        r = matrix(self._K,self._n)
        z4 = CyclotomicField(4).gens()[0]
        F = self.parent().finite_quadratic_module()
        L = F.list()
        for ii in range(0 ,self._n):
            for jj in range(0 ,self._n):
                #if(self._L[ii]==d*self._L[jj] and (filter==None or filter[ii,jj]==1 )):
                if (F[ii]==d*F[jj] and (filter==None or filter[ii,jj]==1) ):
                    argl=self._level*b*d*self.Q(jj)
                    r[ii,jj]=self._zl**argl
        # Compute the character 
        signature = self._inv['signature']
        if( self._level % 4  == 0 ):
            test = (signature + kronecker(-1 ,self._n)) % 4
            if(is_even(test)):
                if(test==0 ):
                    power=1 
                elif(test==2 ):
                    power=-1 
                if( d % 4  == 1 ):
                    chi = 1 
                else:
                    chi=z4**power
                chi=chi*kronecker(c,d)
            else:
                if(test==3 ):
                    chi= kronecker(-1 ,d)
                else:
                    chi=1 
            chi = chi*kronecker(d,self._n*2**signature)
        else:
            chi = kronecker(d,self._n*2**signature)
        r=r*chi
        return [r,1 ]
        
    # Now we want the general action

    def _action_of_SL2Z_formula(self,A,filter=None,**kwds):
        r"""
        The Action of A in SL2(Z) given by a matrix rho_{Q}(A)
        as given by the formula
        filter = |D|x|D| integer matrix with entries 0 or 1
                         where 1 means that we compute this entry
                         of the matrix rho_{Q}(A) 
        """
        ## extract eleements from A
        [a,b,c,d]=_entries(A)
        # check that A is in SL(2,Z)
        #if(A not in SL2Z):
        #()    raise  TypeError,"Matrix must be in SL(2,Z)!"
        ## Check if we have a generator
        sign=1
        if c==0:
            if b==0:
                if a<0:
                    return self._action_of_Z(filter)
                else:
                    return self._action_of_Id(filter)
            if a<1:
                sign=-1
            else:
                sign=1
            return self._action_of_T(b,sign,filter)
        if c % self._level == 0 :
            #print "A1=",A
            return self._action_of_Gamma0(A)
        if abs(c)==1 and a==0:
            if self._verbose>0:
                print "call STn with pos={0} and sign={1}".format(abs(d),sgn(c))
            sic = sgn(c)
            return self._action_of_STn(pow=d*sic,sign=sic,filter=filter)        
        # These are all known easy cases
        # recall we assumed the formula
        if c<0 or (c==0 and d<0): # change to -A
            a=-a; b=-b; c=-c; d=-d
            A=SL2Z(matrix(ZZ,2,2,[a,b,c,d]))
            sign=-1 
        else:
            sign=1
        xis=self._get_xis(A)
        xi=1 
        for q in xis.keys():
            xi=xi*xis[q]
        norms_c=self._get_all_norm_alpha_cs(c)
        #norms_c_old=self._get_all_norm_alpha_cs_old(c)
        if self._verbose>0:
            print "c=",c
            print "xis=",xis
            print "xi=",xi,CC(xi)
            print "norms=",norms_c
        #print "11"
        r = matrix(self._K,self._n)
        if sign==-1:
            #r=r*self._QM.sigma_invariant()**2 
            si = self._QM.sigma_invariant()**2
            if is_odd(self.parent().signature()):
                if c>0 or (c==0 and d<0):
                    si = -si ## sigma(Z,A)
        else:
            si=1
        if self._verbose>0:
            print "si=",si
            print "sign=",sign
        for na in range(0 ,self._n):
            for nb in range(0 ,self._n):
                if filter <> None and filter[na,nb]==0:
                    continue
                if sign==-1:
                    #print type(nb)
                    nbm=self._minus_element(nb) #-alpha
                else:
                    nbm=nb
                #beta=self._L[nb]
                #gamma=alpha-d*beta
                # c*alpha' = 
                gi=self.lin_comb(na,-d,nbm)
                try:
                    ngamma_c=norms_c[gi]
                    #ngamma_c_old=norms_c_old[gamma]
                except KeyError:
                    #print alpha," not in D^c*"                    
                    continue
                #ngamma_c_old=self._norm_alpha_c(gamma,c)
                #arg_old=a*ngamma_c_old+b*self._QM.B(gamma,beta)+b*d*self._QM.Q(beta)
                #gi = self._L.index(gamma)
                #CHECK: + or - ? arg=a*ngamma_c+b*self.B(gi,nbm)+b*d*self.Q(nbm)
                arg=a*ngamma_c+b*self.Bi(na,nbm)-b*d*self.Q(nbm)
                larg=arg*self._level
                if self._verbose>0 and na==0 and nb==1:
                    print "na,nb,nbm=",na,nb,nbm
                    print "gi=",gi
                    print "ngamma_c[{0}]={1}".format(gi,ngamma_c)
                    print "b*B(na,nbm)=",b*self.Bi(na,nbm)
                    print "arg=",               a,"*",ngamma_c,"+",b,"*",self.Bi(na,nbm),"-",b,"*",d,"*",self.Q(nbm)
                    print "arg=",arg
                    print "e(arg)=",CC(0,arg*RR.pi()*2).exp()
                    print "e_L(arg)=",CC(self._zl**(larg))
                #if na==nb:
                #    print "arg[",na,"]=",a*ngamma_c,'+',b*self.B(gi,nbm),'-',b*d*self.Q(nbm),'=',arg

                r[na,nb]=si*xi*self._zl**(larg)
                if self._verbose>0 and na==0 and nb==1:
                    print "r[",na,nb,"]=",r[na,nb]
                    print "r[",na,nb,"]=",r[na,nb].complex_embedding(53)
                #print "xi=",xi
                #print "zl=",self._zl
        fac = self._get_lenDc(c)
        #print "12"
        return [r,(QQ(fac)/QQ(self._n)).sqrt()]


    def _get_lenDc(self,c):
        r"""
        compute the number of elements in the group of elements of order c
        """
        # Check to see if everything is precomputed if so we fetch it
        g = gcd(c,self._n)
        try:
            n=self._All_len_Dcs[g]
        except AttributeError:
            # if the dictionary doesn't exist, create it and store the correct value
            n=self._get_one_lenDc(c)
            self._All_len_Dcs=dict()
            self._All_len_Dcs[g]=n
        except KeyError:
            # The dictionary exist but this value does not 
            n=self._get_one_lenDc(c)
            self._All_len_Dcs[g]=n
        return n

    @cached_method
    def _get_one_lenDc(self,c):
        r"""
        compute the number of elements in the group of elements of order c
        """
        n=0 
        for ii in range(0,self._n):
            x=self._QM(self._W._elt(ii))
            if(c*x==self._QM(0)):
                n=n+1 
        return n

    def _get_all_lenDc(self):
        r"""
        Compute the number of elements in the group of elements of order c
        for all c dividing |D|
        """
        res=dict()
        divs = divisors(self._n)
        for c in divs:
            res[c]=self._get_one_lenDc(c)
        # set all lengths
        self._All_len_Dcs=res


    # Setup the functions for computing the Weil representation
        

    def get_xis(self,A,B=None,C=None,D=None,pset=None):
        return  self._get_xis(A,B=B,C=C,D=D,pset=pset)
    
    def _get_xis(self,A,B=None,C=None,D=None,pset=None):
        r"""
        compute the p-adic factors: \xi_0, \xi_p, p | |D|

        if pset = p we only return the p-factor of xi.
        """
        JD=self._QM.jordan_decomposition()
        absD = self._n
        if D==None:
            [a,b,c,d]=_entries(A)
        else:
            a=A; b=B; c=C; d=D
        if self._verbose>0:
            print "pset=",pset
            print "JD=",JD
        #info=get_factors2(JD)
        oddity=self._inv['total oddity']
        oddities=self._inv['oddity']
        pexcesses=self._inv['p-excess']
        sign=self._inv['signature']
        z8 = self._z8
        gammafactor=1 
        xis=dict()
        for comp in JD:
            p=comp[1][0]
            xis[p]=1
            if self._verbose>0:
                print "p=",p
        if a*d<>0:
            if(is_even(sign)):
                argl=- 2 *sign
                xis[0]=z8**argl
            else:
                dc=kronecker(-a,c)
                argl=-2 *sign
                if(is_even(c)):
                    argl=argl+(a+1 )*(odd_part(c)+1 )
                xis[0]=z8**argl*dc
        else:
            argl=-sign
            xis[0]=z8**argl
            if pset==0:
                return {0:xis[0]}
            elif pset<>None:
                return {0:1}
            return xis
        if pset==-1 or pset==0:
            return {0:xis[0]}
        if(xis.keys().count(2 )>0 ):
            if(is_odd(c)):
                argl=(c*oddity) % 8 
            else:
                argl=(-(1+a)*oddity) % 8 
            xis[2]=z8**argl
            if self._verbose>0:
                print "oddity(Q)=",oddity
                if(is_odd(c)):                
                    print "c*oddity=",argl
                else:
                    print "-(a+1)oddity=",argl
            if pset==2:
                return {2:xis[2]}
        for comp in JD:
            [p,n,r,ep]=comp[1 ][:4 ]
            if self._verbose>0:
                print "comp:p,n,r,ep=",p,n,r,ep
            t = None if (len(comp[1 ])==4 ) else comp[1 ][4 ]
            q=p**n
            qc=gcd(q,c)
            qqc=Integer(q/qc)
            ccq=Integer(c/qc)
            nq=valuation(qqc,p)
            gammaf=1 
            dc=1 
            #if(c % q == 0): # This only contributes a trivial factor
            #    continue
            if p==2:
                if is_even(c) :
                    if c % q <> 0:
                        odt=self._get_oddity(p,nq,r,ep,t)
                        argl=-a*ccq*odt
                        gammaf=z8**argl
                        if self._verbose>0:
                            print "odt({q})_{t}^{e},{r}={o}".format(t=t,o=odt,e=ep,r=r,q=p**nq)
                            print "argl=",argl
                        dc=kronecker(-a*ccq,qqc**r)
                        if self._verbose>0:
                            print "kron(-ac_q/q_c^r)=",dc
                    dc=dc*kronecker(-a,q**r)
                    xis[ 2 ]=xis[ 2 ]*gammaf*dc
                    if self._verbose>0:
                        print "xis2=",xis[2]
                else:
                    dc=kronecker(c,q**r)
                    xis[ 2 ]=xis[ 2 ]*dc
            else:
                if(c%q <>0 ):
                    exc=self._get_pexcess(p,nq,r,ep)
                    argl=-exc
                    gammaf=z8**argl
                    dc=kronecker(ccq,qqc**r)
                    if self._verbose>0:
                        print "gamma_{p}={g}".format(p=p,g=gammaf)
                dc=dc*kronecker(-d,qc**r)
                xis[p]=xis[p]*gammaf*dc
            if pset==p:
                return {p:xis[p]}
        return xis


    def _get_xc(self,c):
        r"""
        Returns x_c which is precomputed at initialization
        """
        k=valuation(c,2 )
        if(k==0 ):
            xc=0 
        else:
            try:
                xc=self._xcs[k]
            except KeyError:
                xc=0 
        return xc
    
    def _compute_all_xcs(self):
        r"""
        Computes all non-zero values of the element x_c in the group D
        OUPUT: dictionsry k => x_c  where 2^k || c
        This routine is called at initialization of WeilModuleElement
        """
        J=self._QM.jordan_decomposition()
        res=dict()
        res[0 ]=0 
        for c in J:
            xc=0 
            p,k,r,d = c[1 ][:4 ]
            t = None if 4  == len(c[1 ]) else c[1 ][4 ]
            if(p<>2  or t==None): 
                continue
            q=2 **k # = 2^n
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
        return res

    @cached_method
    def _get_all_norm_alpha_cs(self,c):
        r"""
        Computes a vector of all Q(alpha_c)
        for alpha in D  (==0 unless alpha_c is in D^c*)
        """
        res=dict()
        for ai in range(self._n):
            #for alpha in self._L:
            nc=self._get_norm_alpha_c(ai,c)
            if nc<>None:
                res[ai]=nc
        return res
    
    @cached_method
    def _get_norm_alpha_c(self,ai,c):
        r"""
        FQM = Finite Quadratic Module
        Test before that alpha is in D^c*!!
        """
        alpha=self._QM(self._W._elt(ai))
        xc=self._get_xc(c)
        if xc<>0:
            gammatmp=alpha-xc
        else:
            gammatmp=alpha
        # We need to find its inverse mod c
        # i.e. gammatmp/c
        # print alpha,gammatmp
        if gcd(c,self._level)==1:
            cc = inverse_mod(c,self._level)
            gamma = (cc*gammatmp).list()
        else:
            gamma=[]
            for jj,g in enumerate(self._QM.gens()):
                for x in range(g.order()):
                    if (c*x - gammatmp.list()[jj]) % g.order() == 0:
                        gamma.append(x)
                        break
            if len(gamma)<len(self._QM.gens()):
                if self._verbose>1:
                    print "c=",c
                    print "x_c=",xc
                    print "gammatmp=",gammatmp
                    print "y=gamma/c=",gamma
                return None
        #gamma=y #vector(y)
        #    # raise ValueError, "Found no inverse (alpha-xc)/c: alpha=%s, xc=%s, c=%s !" %(alpha,xc,c) 
        if self._verbose>0:
            print "xc=",xc
            #print "orders=",self._W._gen_orders
            #print "gamma=",gamma
        if len(gamma)<>len(self._W._gen_orders):
            print "W=",self._W
            print "F=",list(self._W._QM)
            print "F.gens=",self._W._QM.gens()
            print "F.gram=",self._W._QM.gram()
            print "is_nondeg=",self._W._QM.is_nondegenerate()
            print "ai=",ai
            print "c=",c
        gi = self._W._el_index(gamma)
        if self._verbose>0:
            print "gi=",gi
        #res=c*self._QM.Q(gamma)
        res=c*self.Q(gi)
        if xc<>0:
            #res=res+self._QM.B(xc,gamma)
            if self._verbose>0:
                print "xc=",xc
                print "xc.list=",xc.list()
                print "orders=",self._W._gen_orders
            xci = self._W._el_index(xc.list())
            res=res+self.Bi(xci,gi)
        return res


    def _get_all_norm_alpha_cs_old(self,c):
        r"""
        Computes a vector of all Q(alpha_c)
        for alpha in D  (==0 unless alpha_c is in D^c*)
        """
        res=dict()
        for alpha in self._L:
            nc=self._get_norm_alpha_c_old(alpha,c)
            if(nc<>None):
                res[alpha]=nc
        return res
    
    def _get_norm_alpha_c_old(self,alpha,c):
        r"""
        FQM = Finite Quadratic Module
        Test before that alpha is in D^c*!!
        """
        ## Make an extra test. This should be removed in an efficient version
        #Dcs=D_upper_c_star(FQM,c)
        #if(Dcs.count(alpha)==0):
        #    raise ValueError, "alpha=%s is not in D^c* for c=%s, D^c*=%s" %(alpha,c,Dcs)
        xc=self._get_xc(c)
        # print "xc=",xc
        if xc<>0:
            gammatmp=alpha-xc
        else:
            gammatmp=alpha
        # first a simple test to see if gammatmp is not in D^c
        test = gammatmp*Integer(self._n/gcd(self._n,c))
        if test<>self._QM(0):
            return None
        # We need to find its inverse mod c
        # i.e. gammatmp/c
        # print alpha,gammatmp
        try:
            for y in self._L:
                yy=c*y
                #print "yy=",yy
                if yy==gammatmp:
                    gamma=y    #(alpha-xc)/c
                    raise StopIteration
        except StopIteration:
            pass
        else:
            return None
            # raise ValueError, "Found no inverse (alpha-xc)/c: alpha=%s, xc=%s, c=%s !" %(alpha,xc,c)
        res=c*self._QM.Q(gamma)
        if xc<>0:
            res=res+self._QM.Bi(xc,gamma)
        return res
    
    # Defines the action of the Weil representation using factorization into product of ST^k_j
    def _action_of_SL2Z_factor(self,A):
        r""" Computes the action of A in SL2Z on self
             Using the factorization of A into S and T's
             and the definition of rho on these elements.
             This method works for any WeilModule but is (obviously) very slow
             INPUT : A in SL2Z
             OUTPUT: [r,fact]
                     r = matrix over CyclotomicField(lcm(level,8))
                     fact = sqrt(integer) 
                     rho(A) = fact*r
        """
        if A not in SL2Z:
            raise TypeError, "%s must be an element of SL2Z" %A
        #[fak,sgn]=self._factor_sl2z(A)
        sgn,n,fak=factor_matrix_in_sl2z(A)
        #[fak,sgn]=self._factor_sl2z(A)
        fak.insert(0,n) 
        [r,fact,M]=self._weil_matrix_from_list(fak,sgn)
        for i in range(2):
            for j in range(2):
                if A[i,j]<>M[i,j] and A[i,j]<>-M[i,j]:
                    raise ValueError, "\n A=\n%s <> M=\n%s, \n factor=%s " % (A,M,fak)
        return [r,fact]

    def _weil_matrix_from_list(self,fak,sgn):
        r"""
        INPUT: fak = output from _factorsl2z
                   = [a0,a1,...,ak]
               sgn = +-1
        OUTPUT [r,fact,M]
               M       = sgn*T^a0*S*T^a1*...*S*T^ak
               r*fact  = rho(M)
               r       = matrix in CyclotomicField(lcm(level,8))
               
        """
        [S,T]=SL2Z.gens()
        M=SL2Z([1 ,0 ,0 ,1 ])
        Z=SL2Z([-1 ,0 ,0 ,-1 ])
        r = matrix(self._K,self._n)
        ss = self._QM.sigma_invariant()
        ## Do the first T^a
        if sgn==-1:
            r,fac=self._action_of_Z()
            rt,fact=self._action_of_T(fak[0])
            r=r*rt
            fac=fac*fact
            M=Z*T**fak[0 ]
        else:
            r,fac=self._action_of_T(fak[0 ])
            M=T**fak[0 ]
        if self._verbose>2:
            print "M=",M
            print "r=",r
        A=SL2Z([1 ,0 ,0 ,1 ]) # Id
        for j in range(1 ,len(fak)):
            A=A*S*T**fak[j]
        if sgn==-1 and sigma_cocycle(Z,A)==-1:
            sfak=ss**4
        else:
            sfak=1 
        # Now A should be the starting matrix except the first T factor
        fact=1 
        #print "A=\n",A
        for j in range(1 ,len(fak)):
            rN,fac=self._action_of_STn(fak[j])
            Mtmp=S*T**fak[j]   
            M=M*Mtmp
            A=(Mtmp**-1 )*A
            if(j<len(fak)-1 ):
                si=sigma_cocycle(Mtmp,A)
                if si==-1:
                    #sig=Integer(si)
                    sfak=sfak*ss**4
            r=r*rN #*sfak
         #   print "fact=",fact
            fact=fact*self._n
        r=r*sfak
        if self._verbose>0:
            print "sfak=",sfak
        #print "|D|=",self._n
        #print "fact=",fact
        if hasattr(fact,"sqrt"):
            fact = fact.sqrt()**-1
        else:
            fact = 1/sqrt(fact)
        return [r, fact,M]



    @staticmethod
    # Compute total invariants of a Jordan decomposition
    def _get_invariants(self):
        r"""
        INPUT: JD is a Jordan decomposition returned from the FQM package
        OUPUT: dictionary res with entries:
        res['total p-excess']= sum of all p-excesses
        res['total oddity']  = sum of all oddities
        res['oddity']        = list of oddities of all components
        res['p-excess']=     = list of p-excesses of all components
        res['signature']     = signature := (odt-pexc_tot) % 8
        res['type']          = dictionary q=>'I' or 'II' telling whether the 2-adic component q is odd or even
        """
        JD=self._QM.jordan_decomposition()
        pexc=dict()
        odts=dict()
        res=dict()
        types=dict()
        odt=0 
        for c in JD:
            [p,n,r,ep]=c[1 ][:4 ]
            t = 0  if (len(c[1 ])==4 ) else c[1 ][4 ]
            q=p**n
            if( (not is_square(q)) and (ep==-1 )):
                k=1 
            else:
                k=0 
            if(p!=2 ):
                pexc[q]=(r*(q-1 )+4 *k) % 8 
            else:
                odts[q]=((t+4 *k) % 8 ) # Can we have more than one 2-adic component at one time?
                if(t<>0 ):
                    types[q]="I" # odd
                else:
                    types[q]="II" # even
        pexc_tot=0 
        odt=0 
        for q in pexc.keys():
            pexc_tot=pexc_tot+pexc[q]
        for q in odts.keys():
            odt=odt+odts[q]        
        res['total p-excess']=pexc_tot % 8 
        res['total oddity']=odt % 8 
        res['oddity']=odts
        res['p-excess']=pexc
        res['signature']=(odt-pexc_tot) % 8 
        res['type']=types
        return res

    # computes oddity and p-excess of scaled Jordan blocks
    def _get_oddity(self,p,n,r,ep,t):
        r"""
        return the oddity of the Jordan block q_t^(r*ep) where q = p^n  
        """
        if n==0  or p==1:
            return 0 
        k=0 
        if n % 2 <>0  and ep==-1:  # q not a square and sign=-1
            k=1 
        else:
            k=0 
        if(t):
            odt=(t+k*4 ) % 8 
        else:
            odt=(k*4 ) % 8 
        return odt


    def _get_pexcess(self,p,n,r,ep):
        r"""
        return the oddity of the corresponding Jordan block
        q = p^n  and the module is q^(r*ep)
        """
        if(n==0  or p==1 ):
            return 0 
        #if( n % 2 <>0  and ep==-1 ):  # q is a square
        if( is_odd(n) and ep==-1 ):  # q is a square
            k=1 
        else:
            k=0 
        exc=(r*(p**n-1 )+ 4 *k) % 8
        return exc


    @cached_method
    def lin_comb(self,a,d,b):
        x = self._QM(self._W._elt(a))+d*self._QM(self._W._elt(b))
        x=vector(x.list())
        x.set_immutable()
        return self._W._el_index(x)


# def WeilRepresentation(FQM):
#     r"""
#     Construct a dummy element of WeilModule. This is useful to extract informaiton about the Weil representation without constructing an actual element.
#     """

    
#### End of WeilModule Element
def _entries(A):
    r"""
    Returns the entries of A
    where A is one of:
    1) Element of SL2Z
    2) 2x2 integer matrix with determinant 1
    3) list [a,b,c,d] with ad-bc=1
    """
    try: 
        if(A in MatrixSpace(ZZ,2,2)):
            a=A[0 ,0 ]; b=A[0 ,1 ]; c=A[1 ,0 ]; d=A[1 ,1 ]
        else:
            [a,b,c,d]=A
    except TypeError:
        [a,b,c,d]=A
    #if(a*d-b*c<>1 ):
    #    return None
    #else:
    return [a,b,c,d]

def sigma_cocycle(A,B):
    r"""
    Computing the cocycle sigma(A,B) using the Theorem and Hilbert symbols
    
    INPUT:
    
    -''A'' -- matrix in SL(2,Z)
    -''B'' -- matrix in SL(2,Z)
        
    OUTPUT:
    
    -''s'' -- sigma(A,B) \in {1,-1}
    
    EXAMPLES::
    
        
    sage: S,T=SL2Z.gens()     
    
    
    """
    [a1,b1,c1,d1]=_entries(A)
    [a2,b2,c2,d2]=_entries(B)
    if c2*c1<>0:
        C=A*B
        [a3,b3,c3,d3]=_entries(C)
        if c3<>0:
            # print "here",c3*c1,c3*c2
            return hilbert_symbol(c3*c1,c3*c2,-1 )
        else:
            return hilbert_symbol(c2,d3,-1 )
    elif c1<>0:
        return hilbert_symbol(-c1,d2,-1 )
    elif c2<>0:
        return hilbert_symbol(-c2,d1,-1 )
    else:
        return hilbert_symbol(d1,d2,-1 )

    

def kubota_cocycle(A,B):
    r"""
    Computing the cocycle sigma(A,B) using the Theorem and Hilbert symbols
    
    INPUT:
    
    -''A'' -- matrix in SL(2,Z)
    -''B'' -- matrix in SL(2,Z)
        
    OUTPUT:
    
    -''s'' -- sigma(A,B) \in {1,-1}
    
    EXAMPLES::
    
        
    sage: S,T=SL2Z.gens()     
    
    
    """
    [a1,b1,c1,d1]=_entries(A)
    [a2,b2,c2,d2]=_entries(B)
    C=A*B
    [a3,b3,c3,d3]=_entries(C)
    sA = kubota_sigma_symbol(c1,d1)
    sB = kubota_sigma_symbol(c2,d2)
    sC = kubota_sigma_symbol(c3,d3)
    res = hilbert_symbol(sA,sB,-1)*hilbert_symbol(-sA*sB,sAB,-1)
    return res

    


def kubota_sigma_symbol(c,d):
    r"""
    Compute sigma_A=sigma(c,d) for A = (a b // c d)
    given by sigma_A = c if c<>0 else = d
    """
    if c<>0:
        return c
    else:
        return d

from finite_quadratic_module import FiniteQuadraticModuleRandom
    
def test_dimensions_1(fqbound=100,nbound=10,cbound=10,size_bd=50,kmin=2,kmax=10,verbose=0,test=1,numeric=False):
    r"""
    Run tests on a set of random quadratic modules
    test >= 0: make sure only exactly one of dim(r) and dim(r*) arenon-zero
    Test >= 1: check that dimensions are integers 

    
    """
    for i in range(nbound):
        l=size_bd+1        
        while l>size_bd:
            FQ=FiniteQuadraticModuleRandom(fqbound,nbound,verbose-1)
            l=len(list(FQ))
            W = WeilModule(FQ)
            g = FQ.jordan_decomposition().genus_symbol()
            s = W.signature()
            if verbose>0:                
                print "Check genus=",g,"sign=",s
            for twok in range(2*kmin+1,2*kmax+1):
                k = QQ(twok)/QQ(2)
                if verbose>1:
                    print "k=",k
                ## There is a built-in integrality test already
                try:
                    d1,ep1 = W.dimension_cusp_forms(k,sgn=1,verbose=verbose-2,numeric=numeric)
                    d2,ep2 = W.dimension_cusp_forms(k,sgn=-1,verbose=verbose-2,numeric=numeric)
                except ArithmeticError:
                    print "Fail 2"
                    print "d=",d
                    print "genus:",g
                    print "FQ=",FQ
                    print "dimS(rho,{0})={1}".format(k,CC(d1))
                    print "dimS(rho*,{0})={1}".format(k,CC(d2))
                    return False,W,d1,d2
    return True
### testing functions
def test_formula(fqbound=100,nbound=10,cbound=10,size_bd=50,kmin=2,kmax=10,verbose=0,test=1,numeric=False,gamma0_test=0):
    r"""
    Run tests on a set of random quadratic modules
    test >= 0: make sure only exactly one of dim(r) and dim(r*) arenon-zero
    Test >= 1: check that dimensions are integers 

    
    """
    for i in range(nbound):
        l=size_bd+1        
        while l>size_bd:
            FQ=FiniteQuadraticModuleRandom(fqbound,nbound,verbose-1)
            if list(FQ)==[FQ(0)]:
                return FQ
            l=len(list(FQ))
            W = WeilModule(FQ)
            g = FQ.jordan_decomposition().genus_symbol()
            s = W.signature()
        if verbose>0:
            print "signature=",s
            print "genus=",g
        w = W.an_element()
        i = 0
        j = 0
        while (i<cbound):
            A = SL2Z.random_element()
            if gamma0_test==1 and A[1,0] % W.level() <>0:
                j+=1
                if j>2000:
                    raise  ArithmeticError,"Error in random elements of SL2Z!"
                continue
            j = 0
            i = i + 1
            t = compare_formula_for_one_matrix(W,A,verbose,gamma0_test)
            if not t:
                if verbose>0:
                    print "A=",A
                    return W,A
            
            
                        #raise ArithmeticError,"Got different matrix! r1={0} and r2={1}".format(r1,r2)
    return True

def compare_formula_for_one_matrix(W,A,verbose=0):
    r"""
    Compare the action of A on W using the formula and factoring into generators.
    """
    r1 = W.matrix(A,by_factoring=False)
    r2 = W.matrix(A,by_factoring=True)
    f1 = CC(r1[1])
    f2 = CC(r2[1])
    for i in range(r1[0].nrows()):
        for j in range(r1[0].ncols()):
            t1 = f1*r1[0][i,j].complex_embedding()
            t2 = f2*r2[0][i,j].complex_embedding()
            if abs(t1-t2)>1e-10:
                if verbose>0:
                    print "i,j=",i,j
                    print "t1=",t1
                    print "t2=",t2
                return False
    return True

def compare_formula_for_one_module(W,nmats=10,verbose=0):
    r"""
    Compare the action of A on W using the formula and factoring into generators.
    """
    for i in range(nmats):
        A = SL2Z.random_element()
        t = compare_formula_for_one_matrix(W,A,verbose)
        if t==False:
            return False
    return True
