r"""
Classes and routines for working with spaces of (vector-valued) period polynomials of modular forms. 


"""
from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

from builtins import map
from builtins import range
from sage.all import SL2Z,SageObject,is_squarefree,RealField,ComplexField,Matrix,conjugate,binomial,Gamma0,is_prime,ZZ,Cusp,gcd,ceil,floor,log_b,CyclotomicField,divisors


from .period_polynomials_algs import *

class PeriodPolynomialSpace(SageObject):
    r"""
    Space of vector-valued period polynomials for the induced representation
    from $\Gamma_0(N)$
    """
    def __init__(self,N,k,prec=53):
        self._N = N
        self._k = k


    def basis(self):
        r"""
        Compute a basis fof self.
        """





class PolynomialModuleElement(SageObject):
    r"""
    An element of the space of polynomials of degree w equipped with
    the weight -w action of  SL(2,Z) 
    An element of the ambient spacePeriodPolynomialSpace
    """
    def __init__(self,w,coefficients=[]):
        self._w = w
        assert isinstance(coefficients,(list,dict))
        if len(coefficients)==w and isinstance(coefficients,list):
            self._coefficients = coefficients
        elif isinstance(coefficients,list):
            self._coefficients = coefficients
            for k in range(len(coefficients),w):
                self._coefficients.append(0)
        else:
            for k in range(w):
                self._coefficients.append(coefficients.get(k,0))
            
class PeriodPolynomial(SageObject):
    r"""
    Class for (numerical approximations to) period polynomials of modular forms
    
    """
    def __init__(self,f,prec=53,verbose=0,data={}):
        r"""
        Initialize the period polynomial of f.

        INPUT:
         - `f` -- Newform
         - ``
        
        EXAMPLES:

        # First elliptic curve

        sage: f=Newforms(11,2)[0]             
        sage: pp=PeriodPolynomial(f,prec=103)
        sage: L=f.cuspform_lseries()
        sage: pp.petersson_norm()
        0.00390834565612459898524738548154 - 5.53300833748482053164300189100e-35*I
        sage: pp.special_value(1)            
        0.253841860855910684337758923267
        sage: L(1)
        0.253841860855911

        # We can also study rationality of the coefficients of the 
        sage: pp.test_rationality()
        P(f 0)^+/omega_+=1.00000000000000000000000000000
        P(f 0)^-/omega_-=0
        P(f 1)^+/omega_+=-1.00000000000000000000000000000
        P(f 1)^-/omega_-=0
        P(f 2)^+/omega_+=0
        P(f 2)^-/omega_-=1.27412157006452456511355616309e-31 + 1.01489446868406102099789763444e-28*I
        P(f 3)^+/omega_+=-5.00000000000000000000000000133
        P(f 3)^-/omega_-=-1.16793587723237638257725820700e-60 + 8.24993716616779655911027615609e-29*I
        P(f 4)^+/omega_+=-2.50000000000000000000000000073
        P(f 4)^-/omega_-=-3.39765752017206550696948310188e-32 + 2.39999999999999999999999999940*I
        P(f 5)^+/omega_+=2.50000000000000000000000000073
        P(f 5)^-/omega_-=-3.39765752017206550696948310188e-32 + 2.39999999999999999999999999940*I
        P(f 6)^+/omega_+=5.00000000000000000000000000133
        P(f 6)^-/omega_-=1.16793587723237638257725820700e-60 - 8.24993716616779655911027615609e-29*I
        P(f 7)^+/omega_+=5.00000000000000000000000000133
        P(f 7)^-/omega_-=-1.16793587723237638257725820700e-60 + 8.24993716616779655911027615609e-29*I
        P(f 8)^+/omega_+=2.50000000000000000000000000073
        P(f 8)^-/omega_-=3.39765752017206550696948310188e-32 - 2.39999999999999999999999999940*I
        P(f 9)^+/omega_+=-2.50000000000000000000000000073
        P(f 9)^-/omega_-=3.39765752017206550696948310188e-32 - 2.39999999999999999999999999940*I
        P(f10)^+/omega_+=-5.00000000000000000000000000133
        P(f10)^-/omega_-=1.16793587723237638257725820700e-60 - 8.24993716616779655911027615609e-29*I
        P(f11)^+/omega_+=0
        P(f11)^-/omega_-=-1.27412157006452456511355616309e-31 - 1.01489446868406102099789763444e-28*I
        o1= 0.0404001869188632792144191985239*I
        o2= -1.36955018267536771772869542629e-33 - 0.0967407815209822856536332369020*I

        
        sage: f=Newforms(5,8,names='a')[0];f
        q - 14*q^2 - 48*q^3 + 68*q^4 + 125*q^5 + O(q^6)
        sage: L=f.cuspform_lseries()  
        sage: pp=PeriodPolynomial(f,prec=103)
        sage: pp.special_value(4)
        -3.34183629241748201816194829811e-29 + 8.45531852674217908762910051272e-60*I
        sage: pp.special_value(3)
        -0.729970592499451187790964533256 + 2.56020879251697431818575640846e-31*I
        sage: L(3)                  
        -0.729970592499451
        sage: L(4)
        0.000000000000000
        
        



        """
        self._f = f
        self._level = None
        if data != {}:
            self.init_from_dict(data)
        if self._level==None:
            self._level = f.level()
            self._k = f.weight()
            self._w = self._k - 2
            self._verbose = verbose
            self._polynomials = {}
            self._prec = prec
            self._rks = {}
            self._coset_reps={}
            self._dim = Gamma0(self._level).index()
            self._coefficients_at_coset = {}
            self._base_coeffs=[]
            self._base_coeffs_embedding=[]
            self._M = {}
            self._M0 = 0
            self._width={}
            self._shift={}
            self._atkin_lehner={}
            self._integrals={}
            self._pplus = {}
            self._pminus={}
            self._peterson_norm=0
            self._canonical_periods=[]
        # Check
        assert is_squarefree(self._level)
        ## Some numerical definitions
        self._RF = RealField(prec)
        self._eps = self._RF(2)**-self._RF(prec)
        self._pi = self._RF.pi()


    def __repr__(self):
        s="Period Polynomial associated to the weight {0} and level {1} MF {2}".format(self._k,self._level,self._f)
        return s

    def __reduce__(self):
        return (PeriodPolynomial,(self._f,self._prec,self._verbose,self.__dict__))
    

    ## If the dict is non-empty we initialize from it
    def init_from_dict(self,data={}):
        self._level = data.get("_level")
        self._k = data.get("_k",None)
        self._w = data.get("_w")
        self._verbose = data.get("_verbose")
        self._polynomials = data.get("_polynomials")
        self._prec = data.get("_prec")
        self._rks = data.get("_rks")
        self._coset_reps=data.get("_coset_reps")
        self._dim = data.get("_dim")
        self._coefficients_at_coset =data.get("_coefficients_at_coset")
        self._base_coeffs=data.get("_base_coeffs")
        self._base_coeffs_embedding=data.get("_base_coeffs_embedding")
        self._M = data.get("_M")
        self._M0 = data.get("_M0")
        self._width=data.get("_width")
        self._shift=data.get("_shift")
        self._atkin_lehner=data.get("_atkin_lehner")
        self._integrals=data.get("_integrals")
        self._pplus =data.get("_pplus")
        self._pminus=data.get("_pminus")
        self._peterson_norm=data.get("_peterson_norm")
        self._canonical_periods=data.get("_canonical_periods")
        
        

    def level(self):
        return self._level

    def weight(self):
        return self._k

    def function(self):
        return self._f

    def polynomial(self,n):
        if not self._polynomials:
            self._get_polynomials()
        if hasattr(n,"matrix"):
            if n in SL2Z:
                n = self._get_coset_n(n)
            else:
                raise ValueError("{0} not in SL(2,Z)".format(n))
        return self._polynomials.get(n,None)

    def polynomials(self):
        r"""
        Return a list of the polynomials  P(f|A_n)(X), 0<=n<dim
        """
        if not self._polynomials:
            self._get_polynomials()
        return self._polynomials
    

    def special_value(self,n):
        r"""
        Compute L(f,n)
        Recall that L(f,n+1)=(2pii)^(n+1)*(-1)^(n+1)/Gamma(n+1) * r_n(f)
        """
        RF = RealField(self._prec)
        CF = ComplexField(self._prec)
        if n < 0 or n>self._w+1:
            print("{0} is not a special point for this f!".format(n))
            return 
        if not self._polynomials:
            self._get_polynomials()
        rk = self.get_rk(0,n-1)
        return CF(0,2*RF.pi())**(n)*(-1)**(n)/RF(n).gamma()*rk
        #print "Warning: r_{0}(f) not computed!".format(n)


    def polynomial_plus(self,n):
        if n not in self._pplus:
            E = Matrix(ZZ,2,2,[-1,0,0,1])
            p1 = self.polynomial(n)
            p2 = self.slash_action(n,E)
            self._pplus[n]=(p1+p2)/2
            self._pminus[n]=(p1-p2)/2
        return self._pplus[n]

    def polynomial_minus(self,n):
        if n not in self._pminus:
            E = Matrix(ZZ,2,2,[-1,0,0,1])
            p1 = self.polynomial(n)
            p2 = self.slash_action(n,E)
            self._pplus[n]=(p1+p2)/2
            self._pminus[n]=(p1-p2)/2
        return self._pminus[n]        


    def petersson_norm(self,type=1):
        r"""
        Compute the Petersson norm of f.

        ALGORITHM:

        Uses
        (f,f) = < rho_{f}^{+} | T - T^-1 , rho_f^{-}>
        for k even and 
        (f,f) = < rho_{f}^{+} | T - T^-1 , rho_f^{+}>
        for k odd.
        See e.g. Thm 3.3. in Pasol - Popa "Modular forms and period polynomials"
        """
        if self._peterson_norm != 0:
            return self._peterson_norm
        T=SL2Z([1,1,0,1])
        Ti=SL2Z([1,-1,0,1])
        norm = 0
        for n in range(self._dim):
            if type==1:
                p1 = self.slash_action(n,T,sym='plus')
                p2 = self.slash_action(n,Ti,sym='plus')
            else:
                p1 = self.slash_action(n,T,sym='minus')
                p2 = self.slash_action(n,Ti,sym='minus')
            pL = p1 - p2
            #print "Pl=",pL
            if self.function().weight() % 2 == 1:
                if type==1:
                    pR = self.polynomial_plus(n)
                else:
                    pR = self.polynomial_minus(n)
            else:
                if type==1:
                    pR = self.polynomial_minus(n)
                else:
                    pR = self.polynomial_plus(n)
            #print "PR=",pR
            t = self.pair_two_pols(pL,pR)
           
            #print "t=",t
            norm+=t
        CF=ComplexField(self._prec)
        c2 = CF(3)*CF(0,2)**(self._k-1)
        self._peterson_norm = -norm/c2
        return self._peterson_norm



    
    def omega_plus_minus(self):
        r"""
        Compute two  periods $\omega_+$ and $\omega_{-}$
        where we take $\omega_{+}$ to be the highest (non-vanishing)
        coefficient of $\rho(f)$ (corresponding to the trivial coset)
        and $\omega_{-}=||f||/\omega_{+}$.

        OUTPUT:
        - '\omega_{+},\omega_{-}'
        """
        if self._canonical_periods != []:
            return self._canonical_periods
        norm = self.petersson_norm()
        o1 = self.get_rk(0,0)
        o2 = norm/o1
        self._o1 = o1
        self._o2 = o2
        self._canonical_periods = o1,o2
        return o1,o2
    
    def pair_two_pols(self,p1,p2):
        res = 0
        cp1 = p1.coefficients(sparse=False)
        cp2 = p2.coefficients(sparse=False)
        for n in range(self._w+1):
            if n < len(cp1):
                c1 = cp1[n]
            else:
                c1 = 0
            k = self._w -n
            if k < len(cp2):
                c2 = cp2[k]
            else:
                c2 = 0
            term = c1*conjugate(c2)*(-1)**(self._w-n)/binomial(self._w,n)
            res = res + term
        return res/self._dim
            
    def _get_polynomials(self,prec=None):
        r"""
        Compute all P(f|A_n)(X), 0<=n<dim
        """
        if prec==None:
            prec = self._prec
        for n in range(self._dim):
            for k in range(self._w+1):
                self.get_rk(n,k)
            self._polynomials[n]=self._get_polynomial(n)
        return self._polynomials
    
    def _get_polynomial(self,n):
        r"""
        Compute P(f|A_n)(X)
        """
        if n in self._polynomials:
            return self._polynomials[n]
        CF = ComplexField(self._prec)
        X=CF['X'].gens()[0]
        k = self._k
        pols = {}
        ## For each component we get a polynomial
        p = X.parent().zero()
        for l in range(self._w+1):
            rk = self.get_rk(n,self._w-l)
            if self._verbose>0:
                print("rk={0}".format(rk))
            fak = binomial(self._w,l)*(-1)**l
            p+=CF(fak)*CF(rk)*X**l
        self._polynomials[n]=p
        return p
#        self._polynomials = #get_period_pol_from_rks(self._f,self._rks,prec=prec,verbose=self._verbose)

        
    def coset_rep(self,n):
        r"""
        Return coset rep nr. n
        """
        if not self._coset_reps:
            self._set_coset_reps()
        return self._coset_reps[n]


    def _set_coset_reps(self):
        if is_prime(self._level):
            self._coset_reps[0] = SL2Z([1,0,0,1])
            for n in range(self._level):
                self._coset_reps[n+1] = SL2Z([0,-1,1,n])
        else:
            G = Gamma0(self._level)
            n = 0
            for A in G.coset_reps():
                self._coset_reps[n] = A
                n+=1
        
    def _get_coset_n(self,A):
        r"""
        Find n s.t. A in Gamma_0(N)A_n 
        """
        G = Gamma0(self._level)
        if A in G:
            return 0  ## Always use 0 for the trivial coset
        else:
            n=0
            for n in range(self._dim):
                An  = self.coset_rep(n)
                if A*An**-1 in G:
                    return n
                n+=1
        raise ArithmeticError("Can not find coset of A={0}".format(A))


    def _get_shifted_coset_m(self,n,A):
        r"""
        Get the index m s.t. A_n*A**-1 is in Gamma0(N)A_m
        """
        An = self.coset_rep(n)
        if A.det()==1:
            m = self._get_coset_n(An*A**-1)
        elif A.det()==-1 and A[0,1]==0 and A[1,0]==0:
            AE = SL2Z([An[0,0],-An[0,1],-An[1,0],An[1,1]])
            m = self._get_coset_n(AE)
        else:
            raise ValueError("Call with SL2Z element or [-1,0,1,0]. Got:{0}".format(A))
        return m
        
    def slash_action(self,n,gamma,sym='none'):
        r"""
        Act of P^{+/-}(f|A_n) with an element of SL2Z
        """        
        # First shift the coset
        m = self._get_shifted_coset_m(n,gamma)
        #print "m=",m
        if sym=='plus':
            P0 = self.polynomial_plus(m)
        elif sym=='minus':
            P0 = self.polynomial_minus(m)
        else:
            P0 = self.polynomial(m)
        if P0==None:
            return
        X = P0.parent().gens()[0]
        coeffs = P0.coefficients(sparse=False)
        p = P0.parent().zero()
        if hasattr(gamma,"matrix"):
            a,b,c,d=gamma
        else:
            a,b,c,d=gamma.list()
            
        denom = (c*X+d)
        numer = (a*X+b)
        w = self._w
        for k in P0.exponents():
            ak = coeffs[k]
            p+=ak*numer**k*denom**(w-k)
        return p


    def slash_action_vector(self,gamma,sym=0):
        res = {}
        for n in range(self._dim):
            res[n] = self.polynomial_slash_action(n,gamma,sym=sym)
        

    def __mul__(self,l):
        r"""
        Acts on all polynomials of self by l
        """

        res = {}
        for n in range(self._dim):
            if l in SL2Z:
                res[n]=self.slash_action(n,l)
            elif isinstance(l,list):
                res[n] = self.action_by_lin_comb(n,l)
            else:
                raise NotImplementedError

        return res
    
    def action_by_lin_comb(self,n,l,sym=0):
        r"""
        Act on polynomial n of self by the element of the group ring c_1*A_1+...+c_k*A_k
        INPUT:
        -`l` -- list of pairs (c,A) with c complex and A in SL2Z
        """
        p = self.polynomial(n).parent().zero()
        try:
            for c,A in l:            
                p+=c*self.slash_action(n,A,sym=sym)
        except:
            raise ValueError("Need a list of tuples! Got:{0}".format(l))
        return p


    def maxM(self):
        r"""
        Compute the maximum truncation point.
        """
        if self._M0 == 0:
            self.get_shifts_and_widths()
        return self._M0

    def set_M(self,Mset):
        r"""
        Sets the truncation point to at most a fixed integer Mset
        """
        M0 = self.maxM()
        for n in range(self._dim):
            m = self._M[n]
            if m > Mset:
                self._M[n] = Mset
        if M0>Mset:
            self._M0 = Mset
        if self._verbose>0:
            print("Warning: the precision of results might be lower than expected!")
        
    
    def truncation_M(self,n):
        r"""
        Compute the truncation point at coset nr. n.
        """        
        if n not in self._M:
            h,w = self.get_shift_and_width_for_coset(n)
            #self._M[n]=self._get_truncation(w,self._k,self._prec)
            self._M[n]=get_truncation(self._k,w,self._prec)
        return self._M[n]


    
    def coefficients_f(self):
        r"""
        Compute the maximum number of coefficients we need to use (with current precision).
        """
        if not self._base_coeffs:
            kf = 0.5*(self._k-1.0)
            M0 = self.maxM()
            prec1 =  self._prec + ceil(log_b(M0*6,2))
            if hasattr(self._f,"O"):
                self._base_coeffs = self._f.coefficients()
                if len(self._base_coeffs)<M0:
                    raise ValueError("Too few coefficients available!")
            else:
                self._base_coeffs = self._f.coefficients(ZZ(M0))
            self._base_coeffs_embedding=[]
            #precn = prec1
            if hasattr(self._base_coeffs[0],"complex_embedding"):
                n=1
                #print "kf=",kf
                #print "kfprec1=",prec1
                for a in self._base_coeffs:
                    ## Use precision high enough that c(n)/n^((k-1)/2) have self._prec correct bits
                    precn = prec1 + kf*ceil(log_b(n,2))
                    self._base_coeffs_embedding.append(a.complex_embedding(precn))
                    n+=1
            else:
                ## Here we have a rational number so we can ue the default number of bits
                n=1
                for a in self._base_coeffs:
                    precn = prec1 + kf*ceil(log_b(n,2))
                    aa = RealField(precn)(a)
                    self._base_coeffs_embedding.append(a) #self._RF(a))
                    n+=1
        return self._base_coeffs
    
    def coefficient_at_coset(self,n):
        r"""
        Compute Fourier coefficients of f|A_n 
        """
        if hasattr(n,"matrix"):
            n = self._get_coset_n(n)
        if not self._base_coeffs:
            self.coefficients_f()
        CF = ComplexField(self._prec)
        if n not in self._coefficients_at_coset:
            M = self.truncation_M(n)
            h,w = self.get_shift_and_width_for_coset(n)
            zN = CyclotomicField(w).gens()[0]
            prec1 = self._prec + ceil(log_b(M,2))
            RF = RealField(prec1)
            coeffs=[]
            fak  = RF(w)**-(RF(self._k)/RF(2))
            #print "f=",f,n
            for i in range(M):
                al = CF(self.atkin_lehner(n))
                c0 = self._base_coeffs_embedding[i]
                #if hasattr(c0,"complex_embeddings"):
                #    c0 = c0.complex_embedding(prec1)
                #else:
                #    c0 = CF(c0)
                qN = zN**((i+1)*h)
                #print "qN=",qN,type(qN)
                an = fak*al*c0*qN.complex_embedding(prec1)
                #an = self.atkin_lehner(n)*self._base_coeffs[i]*zN**((i+1)*h)
                #an = f*an.complex_embedding(prec1)
                coeffs.append(an)
            self._coefficients_at_coset[n] = coeffs
        return self._coefficients_at_coset[n]

    def get_rk(self,n,k):
        r"""
        Compute the coefficient r_k((f|A_n))
        """
        if (n,k) not in self._rks:
            S,T = SL2Z.gens()
            nstar = self._get_shifted_coset_m(n,S)
            kstar = self._k-2-k
            i1 = self.get_integral_from_1_to_oo(n,k)
            i2 = self.get_integral_from_1_to_oo(nstar,kstar)
            if self._verbose>0:
                print("i1={0}".format(i1))
                print("i2={0}".format(i2))
            self._rks[(n,k)] = i1+i2*(-1)**(k+1)
        return self._rks[(n,k)]

    def get_integral_from_1_to_oo(self,n,k):
        r"""
        Compute \int_{1}^{+Infinity} f|A_n(it) t^k dt 
        """
        CF = ComplexField(self._prec)
        if (n,k) not in self._integrals:
            c = self.coefficient_at_coset(n)
            w = self.width(n)
            my_integral = 0
            if self._prec<=53:
                if w==1:
                    l = 1
                    for a in c: 
                        my_integral+= a*exp_z_integral(CF(0,1),l,k)
                        l+=1
                else:
                    l = 1
                    for a in c: 
                        my_integral+= a*exp_z_integral_w(CF(0,1),l,w,k)
                        l+=1
            else:
                l = 1
                for a in c: 
                    my_integral+= a*exp_z_integral_wmp(CF(0,1),l,w,k)
                    l+=1
            self._integrals[(n,k)]=my_integral
        return self._integrals[(n,k)]

            
    def width(self,n):
        r"""
        Compute the width of the cusp at coset nr. n.
        """
        if n not in self._width:
            self.get_shifts_and_widths()
        return self._width[n]
        
    def atkin_lehner(self,n):
        r"""
        Compute the Atkin-Lehner eigenvalue at the cusp of coset nr. n.
        """
        if n not in self._atkin_lehner:
            A = self.coset_rep(n)
            c = Cusp(A[0,0],A[1,0])
            G=Gamma0(self._level)
            for Q in divisors(self._level):
                cQ = Cusp(Q,self._level)
                if G.are_equivalent(c,cQ):
                    self._atkin_lehner[n]=self._f.atkin_lehner_eigenvalue(Q)
                    break
        return  self._atkin_lehner[n]
    
    def get_shifts_and_widths(self):
        r"""
        Compute the shifts and widths of all cosets.
        """
        self._width[0]=1
        self._shift[0]=0
        for n in range(0,self._dim):
            h,w = self.get_shift_and_width_for_coset(n)            
            self._M[n] = get_truncation(self._k,w,self._prec)
        self._M0 = max(self._M.values())
            
    def get_shift_and_width_for_coset(self,n):
        r"""
        Compute the shift and width of coset nr. n.
        """
        if n in self._shift and n in self._width:
            return self._shift[n],self._width[n]
        if not is_squarefree(self._level):
            raise NotImplementedError("Only square-free levels implemented")
        G = Gamma0(self._level)
        A = self.coset_rep(n)
        a = A[0,0]; c=A[1,0]
        width = G.cusp_width(Cusp(a,c))
        if self._verbose>0:
            print("width={0}".format(width))
        Amat = Matrix(ZZ,2,2,A.matrix().list())
        h = -1 
        for j in range(width):
            Th = Matrix(ZZ,2,2,[width,-j,0,1])
            W = Amat*Th
            if self._verbose>1:
                print("W={0}".format(W))
                
            if self.is_involution(W):
                h = j
                break
        if h<0:
            raise ArithmeticError("Could not find shift!")
        self._shift[n]=h
        self._width[n]=width
        return h,width        

    def is_involution(self,W,verbose=0):
        r"""
        Explicit test if W is an involution of Gamma0(self._level)
        """
        G = Gamma0(self._level)
        for g in G.gens():
            gg = Matrix(ZZ,2,2,g.matrix().list())
            g1 = W*gg*W**-1
            if verbose>0:
                print("WgW^-1={0}".format(g1))
            if g1 not in G:
                return False
        W2 = W*W
        c = gcd(W2.list())
        if c>1:
            W2 = W2/c
        if verbose>0:
            print("W^2={0}".format(W2))
        if W2 not in G:
            return False
        return True
    
    def test_relations(self):
        
        I = SL2Z([1,0,0,1])
        S = SL2Z([0,-1,1,0])
        U = SL2Z([1,-1,1,0])
        U2 = SL2Z([0,-1,1,-1])
        lS = [(1,I),(1,S)]
        lU = [(1,I),(1,U),(1,U2)]
        pS = self*lS
        pU = self*lU 
        maxS = 0; maxU=0
        if self._verbose>0:
            print("pS={0}".format(pS))
            print("pU={0}".format(pU))
        for p in pS.values():
            if hasattr(p,"coeffs"):
                if p.coefficients(sparse=False) != []:
                    tmp = max(list(map(abs,p.coefficients(sparse=False))))
                else:
                    tmp = 0
            else:
                tmp = abs(p)
            if tmp>maxS:
                maxS = tmp
        for p in pU.values():
            if hasattr(p,"coefficients"):
                if p.coefficients(sparse=False) != []:
                    tmp = max(list(map(abs,p.coefficients(sparse=False))))
                else:
                    tmp = 0
            else:
                tmp = abs(p)
            if tmp>maxU:
                maxU = tmp                
        print("Max coeff of self|(1+S)={0}".format(maxS))
        print("Max coeff of self|(1+U+U^2)={0}".format(maxU))


    def test_relations_plusminus(self):
        
        I = SL2Z([1,0,0,1])
        S = SL2Z([0,-1,1,0])
        U = SL2Z([1,-1,1,0])
        U2 = SL2Z([0,-1,1,-1])
        lS = [(1,I),(1,S)]
        lU = [(1,I),(1,U),(1,U2)]
        ppp1={};ppp2={};ppp3={};ppp4={}
        for n in range(self._dim):
            p1 = self.action_by_lin_comb(n,lS,sym=1)
            p2 = self.action_by_lin_comb(n,lU,sym=1)
            ppp1[n]=p1
            ppp2[n]=p2
            p1 = self.action_by_lin_comb(n,lS,sym=-1)
            p2 = self.action_by_lin_comb(n,lU,sym=-1)
            ppp3[n]=p1
            ppp4[n]=p2

        for pp in [ppp1,ppp2,ppp3,ppp4]:
            maxv=0
            for p in pp.values():
                if hasattr(p,"coeffs"):
                    if p.coefficients(sparse=False)!=[]:
                        tmp = max(list(map(abs,p.coefficients(sparse=False))))
                    else:
                        tmp = 0
                else:
                    tmp = abs(p)
                if tmp>maxv:
                    maxv = tmp
            print("maxv={0}".format(maxv))

        #print "Max coeff of self|(1+S)=",maxS
        #print "Max coeff of self|(1+U+U^2)=",maxU

    def test_rationality(self):        
        omega1 = self.canonical_periods()[0]
        omega2 = self.canonical_periods()[1]
        for n in range(self._dim):          
            pp = self.polynomial_plus(n)/omega1
            pm = self.polynomial_minus(n)/omega2
            print("P(f{0:2d})^+/omega_+={1}".format(n,pp))
            print("P(f{0:2d})^-/omega_-={1}".format(n,pm))
        print("o1={0}".format(omega1))
        print("o2={0}".format(omega2))


        
         
