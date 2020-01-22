# -*- coding: utf-8 -*-
#*****************************************************************************
#  Copyright (C) 2010 Fredrik Strömberg <stroemberg@mathematik.tu-darmstadt.de>,
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

The Weil representation corresponding to the discriminant form $D$ of a rank-one lattice $L=(ZZ, q:x-> Nx**2)$, i.e. $D=(Z/NZ,q mod 1)$.

"""
from __future__ import print_function

from sage.all import Parent,QQ,ZZ,Integer,SL2Z,CyclotomicField,lcm,odd_part,kronecker,gcd,IntegerModRing,matrix,is_odd,\
    valuation,sqrt,MatrixSpace,CC,powerset,squarefree_part,is_even,floor,QuadraticField,is_fundamental_discriminant,\
    is_square,latex,numerator,denominator,prime_divisors
from sage.arith.all import fundamental_discriminant
from sage.rings.complex_mpc import MPComplexField
from sage.misc.cachefunc import cached_function,cached_method
from psage.matrix.matrix_complex_dense import Matrix_complex_dense
from sage.misc.cachefunc import cached_function,cached_method


class WeilRepDiscriminantForm(Parent):
    r""" An elementary version of the Weil representation of the finite quadratic module
    given by D=Z/2NZ.
    """
    def __init__(self,N,k=None,dual=False,sym_type=0,verbose=0):
        r""" Creates a Weil representation (or its dual) of the discriminant form given by D=Z/2NZ.

        
        EXAMPLES::


            sage: WR=WeilRepDiscriminantForm(1,dual=True)
            sage: WR.D
            [0, 1/2]
            sage: WR.D_as_integers
            [0, 1]
            sage: WR.Qv
            [0, -1/4]
            sage: WR=WeilRepDiscriminantForm(1,dual=False)
            sage: WR.D
            [0, 1/2]
            sage: WR.D_as_integers
            [0, 1]
            sage: WR.Qv
            [0, 1/4]

            
        
        """
        ## If N<0 we use |N| and set dual rep. to true
        self._verbose = verbose
        if N<0:
            self._N=-N
            self.dual = not dual
            self._is_dual_rep= not dual # do we use dual representation or not
        else:
            self._N=N
            self._is_dual_rep=dual
        
        N2=Integer(2*self._N)
        self.group=SL2Z
        self._level=4*self._N
        self._D_as_integers=range(0,N2)
        self._even_submodule=[]
        self._odd_submodule=[]
        self._D=list()
        for x in range(0,N2):
            y=QQ(x)/QQ(N2)
            self._D.append(y)
        self.Qv=list()              # List of +- q(x) for x in D
        self.Qv_times_level=list()      # List of +- 4N*q(x) for x in D
        if self._is_dual_rep: # we add this already here for efficiency
            sig=-1
        else:
            sig=1
        for x in self._D:
            y=sig*self.Q(x)
            self.Qv.append(y)
            self.Qv_times_level.append(self._level*y)

        self._signature = sig
        self._sigma_invariant = CyclotomicField(8).gens()[0]**-self._signature
        self._rank = N2
        self._weight = None
        self._sym_type = sym_type
        if sym_type==0 and k != None: # Then we set it
            self._weight = QQ(k)
            if ((self._weight-QQ(1/2)) % 2) == 0:
                sym_type = sig
            elif ((self._weight-QQ(3/2)) % 2) == 0:
                sym_type = -sig
            else:
                raise ValueError("Got incompatible weight and signature!")
        elif sym_type != 0 and k == None:  ## Set the weight
            if sig==sym_type:
                self._weight = QQ(1)/QQ(2)
            elif sig==-sym_type:
                self._weight = QQ(3)/QQ(2)
            else:
                raise ValueError("Got incompatible symmetry type and signature!")
        elif sym_type==0 and k==None:  ## Set the weight
            ## We make a choice
            self._sym_type = sig
            self._weight = QQ(1)/QQ(2)
        else:
            ## Check consistency
            self._weight = QQ(k)
            if ((self._weight-QQ(1/2)) % 2) == 0 and self._sym_type == sig:
                pass
            elif ((self._weight-QQ(3/2)) % 2) == 0 and  self._sym_type == -sig:
                pass

            else:
                raise ValueError("Need either sym type or weight!")
            


    def list(self):
        return self._D

    def N(self):
        return self._N

    def basis(self):
        return self._D
    
    def rank(self):
        return self._rank

    def even_submodule(self,indices=0):        
        if self._even_submodule==[]:
            self._even_submodule = range(0,self._N+1)
        return self._even_submodule

    def odd_submodule(self,indices=0):
        if self._odd_submodule==[]:
            self._odd_submodule = range(1,self._N)
        return self._odd_submodule
    
    def __reduce__(self):
        r""" Used for pickling.
        """
        return(WeilRepDiscriminantForm,(self._N,self._is_dual_rep))

    def __eq__(self,other):
        if(not isinstance(other,WeilRepDiscriminantForm)):
            return False
        return (self._N==other._N) and (self._is_dual_rep==other._is_dual_rep)

    def __ne__(self,other):
        if(not isinstance(other,WeilRepDiscriminantForm)):
            return True
        return not ((self._N==other._N) and (self._is_dual_rep==other._is_dual_rep))
    
    
    def signature(self):
        return self._signature

    def _repr_(self):
        r"""
        Returns string representation of self.

        EXAMPLES::


            sage: WR=WeilRepDiscriminantForm(1,dual=False);WR
            Weil representation of the discriminant form given by ZZ/2ZZ with quadratic form Q(x)=1*x**2 mod 1.
            sage: WR=WeilRepDiscriminantForm(1,dual=True);WR
            Dual of Weil representation of the discriminant form given by ZZ/2ZZ with quadratic form Q(x)=1*x**2 mod 1.
            
       
        
        """
        if self._is_dual_rep:
            s="Dual of "
        else:
            s=""
        s+="Weil representation of the discriminant form given by ZZ/"+str(2*self._N)+"ZZ with quadratic form Q(x)="+str(self._N)+"*x**2 mod 1."
            
        return s

    def _latex_(self):
        r""" Returns LaTeX string representation of self. 

        EXAMPLES::


            sage: WR=WeilRepDiscriminantForm(2,dual=False)
            sage: latex(WR)
            Weil representation of the discriminant form given by $\mathbb{Z}/4\mathbb{Z}$ with quadratic form $Q(x)=2\,x^{2} \mathrm{mod} 1$.


        """
        
        s="\\begin{verbatim}\\end{verbatim}"
        if self._is_dual_rep:
            s+="Dual of "
        else:
            s+=""
            #        s+="Weil representation of the discriminant form given by  $\\mathbb{Z}/"+str(2*self._N)+"\\mathbb{Z}$ \\text{ with quadratic form } Q(x)="+latex(self._N)+"\\,x^{2}\\, \\mathrm{mod}\\, 1$ .\end{verbatim}}"
        s+="Weil representation of the discriminant form given by  $\\mathbb{Z}/"+str(2*self._N)+"\\mathbb{Z}$"
        s+=" with quadratic form  $Q(x)="+latex(self._N)+"\\,x^{2}\\, \\mathrm{mod}\\, 1$."
            
        return s

    def is_dual(self):
        r"""
        Returns True if we have the dual Weil representation, otherwise False.

        EXAMPLES::


            sage: WR=WeilRepDiscriminantForm(1,dual=True);WR.is_dual()
            True
            sage: WR=WeilRepDiscriminantForm(1,dual=False);WR.is_dual()
            False

        
        """
        return self._is_dual_rep
        
    def Q(self,x):
        r"""
        Quadratic form on x, Q(x) mod 1

        INPUT:
          -''x'' -- rational

        OUTPUT:
        -''Q(x'' -- rational

        
        EXAMPLES::

            sage: DF=DiscriminantForm(1,False)
            sage: DF.Q(1/2)
            1/4

        """
        r=self._N*x*x
        p=r.numerator()
        q=r.denominator()
        res=QQ(p % q)/QQ(q)
        return res

    def D(self):
        return self._D_as_integers

    

    @cached_method
    def neg_index(self,a):
        if a not in self._D_as_integers:
            raise ValueError
        ma = (2*self._N - a) % (2*self._N)
        return ma
    
    def B(self,x,y):
        r"""
        Bilinear form B(x,y) mod 1, givenby the quadratic form Q

        INPUT:
        
          -''x'' -- rational
          -''y'' -- rational

        OUTPUT:

        -''B(x,y)'' -- rational

        
        EXAMPLES::

            sage: WR=WeilRepDiscriminantForm(3,dual=True)
            sage: WR.B(1/6,1/2)
            1/2
            sage: WR.B(1/6,1/6)
            1/6
            sage: WR.B(1/6,-1+1/6)
            1/6


        """
        #print "N=",self._N,x,y
        r=Integer(2)*self._N*x*y
        p=r.numerator()
        q=r.denominator()
        res=QQ(p % q)/QQ(q)
        return res

    def sigma_invariant(self):
        return self._sigma_invariant

    def negative_element(self,r):
        r"""
        Return the negative of r in the abelian group of self.
        """
        if r in self._D:
            minus_r = QQ(1 - r)
        elif r in self._D_as_integers:
            minus_r = self._N*2 - r 
        else:
            raise ValueError("Need element in the abelian group of self! Got {0}".format(r))
        return minus_r
            
    def Qc(self,c,x):
        r""" compute Q_c(x)  for x in D^c*
        """
        Dcstar=self._D_times_c_star(c)
        if (not x in Dcstar):
            raise ValueError(" Call only for x in D^c*! Got x={0} and D^c*={1}".format(x,Dcstar))
        xc=0
        if(valuation(c,2)==valuation(2*self._N,2)):
            xc=QQ(1)/QQ(2)
        cy=x-xc
        Dc=self._D_times_c(c)
        for y in Dc:
            p=numerator(y*c)
            q=denominator(y*c)
            if( QQ(p%q)/QQ(q) == QQ(cy)):
                Qc=c*self.Q(y)+self.B(xc,y)
                return Qc
        return ArithmeticError," Could not find y s.t. x=x_c+cy! x=%s and c=%s " %(x,c)
    
    ###  We now add functions for computing the corresponding Weil representation
    def xi(self,A):
        r""" The eight-root of unity in front of the Weil representation.

        INPUT:
        
        -''N'' -- integer
        -''A'' -- element of PSL(2,Z)

        EXAMPLES::

        
            sage: A=SL2Z([41,77,33,62])
            sage: WR.xi(A)
            -zeta8^3]
            sage: S,T=SL2Z.gens()
            sage: WR.xi(S)
            -zeta8^3
            sage: WR.xi(T)
            1
            sage: A=SL2Z([-1,1,-4,3])
            sage: WR.xi(A)
            -zeta8^2
            sage: A=SL2Z([0,1,-1,0])
            sage: WR.xi(A)
            -zeta8

        """
        a=Integer(A[0,0]); b=Integer(A[0,1])
        c=Integer(A[1,0]); d=Integer(A[1,1])
        if(c==0):
            return 1
        z=CyclotomicField(8).gen()    
        N=self._N
        N2=odd_part(N)
        Neven=ZZ(2*N).divide_knowing_divisible_by(N2)
        c2=odd_part(c)
        Nc=gcd(Integer(2*N),Integer(c))
        cNc=ZZ(c).divide_knowing_divisible_by(Nc)
        f1=kronecker(-a,cNc)
        f2=kronecker(cNc,ZZ(2*N).divide_knowing_divisible_by(Nc))
        if(is_odd(c)):
            s=c*N2
        elif( c % Neven == 0):
            s=(c2+1-N2)*(a+1)
        else:
            s=(c2+1-N2)*(a+1)-N2*a*c2
        r=-1-QQ(N2)/QQ(gcd(c,N2))+s
        xi=f1*f2*z**r
        return xi

    def __call__(self,A):
        return self.rho(A)
    
    def matrix(self,M,silent=0,numeric=0,prec=-1):
        return self.rho(M,silent,numeric,prec)

    def rho(self,M,silent=0,numeric=0,prec=-1):
        r""" The Weil representation acting on SL(2,Z).

        INPUT::

        -``M`` -- element of SL2Z
        - ''numeric'' -- set to 1 to return a Matrix_complex_dense with prec=prec instead of exact
        - ''prec'' -- precision
        EXAMPLES::
        
            sage: WR=WeilRepDiscriminantForm(1,dual=False)
            sage: S,T=SL2Z.gens()
            sage: WR.rho(S)
            [
            [-zeta8^3 -zeta8^3]
            [-zeta8^3  zeta8^3], sqrt(1/2)
            ]
            sage: WR.rho(T)
            [
            [       1        0]
            [       0 -zeta8^2], 1
            ]
            sage: A=SL2Z([-1,1,-4,3]); WR.rho(A)
            [
            [zeta8^2       0]
            [      0       1], 1
            ]
            sage: A=SL2Z([41,77,33,62]); WR.rho(A)
            [
            [-zeta8^3  zeta8^3]
            [   zeta8    zeta8], sqrt(1/2)
            ]

        """
        N=self._N; D=2*N; D2=2*D
        if numeric==0:
            K=CyclotomicField (lcm(4*self._N,8))
            z=K(CyclotomicField(4*self._N).gen())
            rho=matrix(K,D)
        else:
            CF = MPComplexField(prec)
            RF = CF.base()
            MS = MatrixSpace(CF,int(D),int(D))
            rho = Matrix_complex_dense(MS)
            #arg = RF(2)*RF.pi()/RF(4*self._N)
            z = CF(0,RF(2)*RF.pi()/RF(4*self._N)).exp()
        [a,b,c,d]=M
        fak=1; sig=1
        if c<0:
            # need to use the reflection 
            # r(-A)=r(Z)r(A)sigma(Z,A)  where sigma(Z,A)=-1 if c>0
            sig=-1
            if numeric==0:
                fz=CyclotomicField(4).gen() # = i
            else:
                fz=CF(0,1)
            # the factor is rho(Z) sigma(Z,-A)
            #if(c < 0 or (c==0 and d>0)):
            #    fak=-fz
            #else:
            #sig=1
            #fz=1
            fak=fz
            a=-a; b=-b; c=-c; d=-d;
            A=SL2Z([a,b,c,d])
            if numeric==0:
                chi=self.xi(A)            
            else:
                chi=CF(self.xi(A).complex_embedding(prec))
            if(silent>0):
                print("fz={0}".format(fz))
                print("chi={0}".format(chi))
        elif c == 0: # then we use the simple formula
            if d < 0:
                sig=-1
                if numeric == 0:
                    fz=CyclotomicField(4).gen()
                else:
                    fz=CF(0,1)
                fak=fz
                a=-a; b=-b; c=-c; d=-d;
            else:
                fak=1
            for alpha in range(D):
                arg=(b*alpha*alpha ) % D2
                if(sig==-1):
                    malpha = (D - alpha) % D
                    rho[malpha,alpha]=fak*z**arg
                else:
                    #print "D2=",D2
                    #print "b=",b
                    #print "arg=",arg
                    rho[alpha,alpha]=z**arg
            return [rho,1]
        else:
            if numeric==0:
                chi=self.xi(M)            
            else:
                chi=CF(self.xi(M).complex_embedding(prec))
        Nc=gcd(Integer(D),Integer(c))
        #chi=chi*sqrt(CF(Nc)/CF(D))
        if( valuation(Integer(c),2)==valuation(Integer(D),2)):
            xc=Integer(N)
        else:
            xc=0
        if silent>0:
            print("c={0}".format(c))
            print("xc={0}".format(xc))
            print("chi={0}".format(chi))
        for alpha in range(D):
            al=QQ(alpha)/QQ(D)
            for beta in range(D):
                be=QQ(beta)/QQ(D)
                c_div=False
                if(xc==0):
                    alpha_minus_dbeta=(alpha-d*beta) % D
                else:
                    alpha_minus_dbeta=(alpha-d*beta-xc) % D
                if silent > 0: # and alpha==7 and beta == 7):
                    print("alpha,beta={0},{1}".format(alpha,beta))
                    print("c,d={0},{1}".format(c,d))
                    print("alpha-d*beta={0}".format(alpha_minus_dbeta))
                invers=0
                for r in range(D):
                    if (r*c - alpha_minus_dbeta) % D ==0:
                        c_div=True
                        invers=r
                        break
                if c_div and silent > 0:
                    print("invers={0}".format(invers))
                    print(" inverse(alpha-d*beta) mod c={0}".format(invers))
                elif(silent>0):
                    print(" no inverse!")
                if(c_div):
                    y=invers
                    if xc==0:
                        argu=a*c*y**2+b*d*beta**2+2*b*c*y*beta
                    else:
                        argu=a*c*y**2+2*xc*(a*y+b*beta)+b*d*beta**2+2*b*c*y*beta
                    argu = argu % D2
                    tmp1=z**argu  # exp(2*pi*I*argu)
                    if silent>0:# and alpha==7 and beta==7):
                        print("a,b,c,d={0},{1},{2},{3}".format(a,b,c,d))
                        print("xc={0}".format(xc))
                        print("argu={0}".format(argu))
                        print("exp(...)={0}".format(tmp1))
                        print("chi={0}".format(chi))
                        print("sig={0}".format(sig))
                    if sig == -1:
                        minus_alpha = (D - alpha) % D
                        rho[minus_alpha,beta]=tmp1*chi
                    else:
                        rho[alpha,beta]=tmp1*chi
        #print "fak=",fak
        if numeric==0:            
            return [fak*rho,sqrt(QQ(Nc)/QQ(D))]
        else:
            return [CF(fak)*rho,RF(sqrt(QQ(Nc)/QQ(D)))]


    def level(self):
        return self._level
        
    def from_discriminant(self,D):
        r"""
        Return the (r,n) s.t. D=n+-q(r).
        """
        ZI=IntegerModRing(self._level)
        if(self.is_dual()):
            x=ZI(-D)
        else:
            x=ZI(D)
        for j in self._D:
            x=self.Qv[j]
            n=QQ(D)/QQ(self._level)-QQ(x)
            if(n % self._level == 0):
                print("D/4N-q(v)={0}".format(n))
                return (self._D[j],ZZ(QQ(n)/QQ(self._level)))


    def _xc(self,c,as_int=False):
        r"""
        Return the element x_c of order 2 (for this Discriminant form x_c=0 or 1/2)

        INPUT:
        -''c'' -- integer
        -''as_int'' -- logical, if true then we return the set D^c as a list of integers
        
        """
        x_c=0
        if(valuation(2*self._N,2)==valuation(c,2)):
            if(as_int):
                x_c=self._N
            else:
                x_c=QQ(1)/QQ(2)
        return x_c


    def _D_times_c(self,c,as_int=False):
        r"""
        Return the set D^c={cx | x in D}
        INPUT:
        -''c'' -- integer
        -''s_int'' -- logical, if true then we return the set D^c as a list of integers
        """
        Dc=list()
        if(as_int):
            setD=self._D_as_integers
        else:
            setD=self._D
        for x in setD:
            if(as_int):
                z=(c*x) % len(self._D)
            else:
                y=c*x
                p=y.numer(); q=y.denom(); z=QQ(p % q)/QQ(q)
            #print "c*",x,"=",z
            Dc.append(z)
        Dc.sort()
        # make unique
        for x in Dc:
            i=Dc.count(x)
            if(i>1):
                for j in range(i-1):
                    Dc.remove(x)
        return Dc


    def _D_lower_c(self,c,as_int=False):
        r"""
        Return the set D_c={x in D| cx = 0}
        INPUT:
        -''c'' -- integer
        -''s_int'' -- logical, if true then we return the set D^c as a list of integers
        """
        Dc=list()
        if(as_int):
            setD=self._D_as_integers
        else:
            setD=self._D
        for x in setD:
            if(as_int):
                z=(c*x) % len(self._D)
            else:
                y=c*x
                p=y.numer(); q=y.denom(); z=QQ(p % q)/QQ(q)
            #print "c*",x,"=",z
            if(z==0):
                Dc.append(x)
        Dc.sort()
        # make unique
        for x in Dc:
            i=Dc.count(x)
            if(i>1):
                for j in range(i-1):
                    Dc.remove(x)
        return Dc


    def _D_times_c_star(self,c,as_int=False):
        r"""
        Return the set D^c*=x_c+{c*x | x in D}, where x_c=0 or 1/2

        INPUT:
        -''c'' -- integer
        -''as_int'' -- logical, if true then we return the set D^c as a list of integers
        
        """
        Dc=self._D_times_c(c,as_int)
        Dcs=list()
        x_c=self._xc(c,as_int)
       
        for x in Dc:
            if(as_int):
                z=(x + x_c) % len(self._D)
            else:
                y=QQ(c*x)+x_c
                p=y.numer(); q=y.denom(); z=QQ(p % q)/QQ(q)
            #print "c*",x,"=",z
            Dcs.append(z)
        Dcs.sort()
        # make unique
        for x in Dcs:
            i=Dcs.count(x)
            if(i>1):
                for j in range(i-1):
                    Dcs.remove(x)
        return Dcs

    def maximal_isotropic_subgroup(self):
        r"""
        Returns the maximal isotropic subgroup of self. 
        """
        S=list()
        for a in self._D:
            if self.Q(a) == 0 and a != 0:
                S.append(a)
        # S is now a list of all isotropic elements except 0
        PS=list(powerset(S))
        PS.reverse()
        # PS now contains all subsets of isotropic elements (except 0)
        # with decreasing sizes. We now need to find the first, which together with 0
        # is a group.
        #print "PS=",PS
        for A in PS:
            A.append(0)
            #print "Test the isotropic set: S=",A
            ok=True
            for x in A:
                for y in A:
                    z=red_mod1(x+y)
                    if(not z in A):
                        #print "S was not a subgroup!"
                        ok=False
            if(ok):
                A.sort()
                return A
        raise ArithmeticError("Could not find maximal isotropic subgroup!")
    

    def dimension_cusp_forms(self,k,eps=1):
        return self._dimension_formula(k,eps=eps,cuspidal=1)

    def dimension_modular_forms(self,k,eps=1):
        return self._dimension_formula(k,eps=eps,cuspidal=0)        
        
    @cached_method
    def _dimension_formula(self,k,eps=1,cuspidal=1):
        ep = 0
        N = self._N
        if (2*k) % 4 == 1: ep = 1
        if (2*k) % 4 == 3: ep = -1
        if ep==0: return 0,0
        if eps==-1:
            ep = -ep
        twok = ZZ(2*k)
        K0 = 1
        sqf = ZZ(N).divide_knowing_divisible_by(squarefree_part(N))
        if sqf>12:
            b2 = max(sqf.divisors())
        else:
            b2 = 1
        b = sqrt(b2)
        if ep==1:
            K0 = floor(QQ(b+2)/QQ(2))
        else:
            # print "b=",b
            K0 = floor(QQ(b-1)/QQ(2))
        if is_even(N):
            e2 = ep*kronecker(2,twok)/QQ(4)
        else:
            e2 = 0
        N2 = odd_part(N)
        N22 = ZZ(N).divide_knowing_divisible_by(N2)
        k3 = kronecker(3,twok)
        if gcd(3,N)>1:
            if eps==1:
                e3 = -ep*kronecker(-3,4*k+ep-1)/QQ(3)
            else:
                e3 = -1*ep*kronecker(-3,4*k+ep+1)/QQ(3)
            #e3 = -1/3*ep
        else:
            f1 = kronecker(3,2*N22)*kronecker(-12,N2) - ep
            f2 = kronecker(-3,twok+1)
            e3 = f1*f2/QQ(6)
        ID = QQ(N+ep)*(k-1)/QQ(12)
        P = 0
        for d in ZZ(4*N).divisors():
            dm4=d % 4
            if dm4== 2 or dm4 == 1:
                h = 0
            elif d == 3:
                h = QQ(1)/QQ(3)
            elif d == 4:
                h = QQ(1)/QQ(2)
            else:
                h = class_nr_pos_def_qf(-d)
            if self._verbose>1:
                print("h({0})={1}".format(d,h))
            if h!=0:
                P= P + h
        P = QQ(P)/QQ(4)
        if self._verbose>0:
            print("P={0}".format(P))
        P=P + QQ(ep)*kronecker(-4,N)/QQ(8)
        if eps==-1:
            P = -P
        if self._verbose>0:
            print("P={0}".format(P))
        # P = -2*N**2 + N*(twok+10-ep*3) +(twok+10)*ep-1
        if self._verbose>0:
            print("ID={0}".format(ID))
        P =  P - QQ(1)/QQ(2*K0)
        # P = QQ(P)/QQ(24) - K0
        # P = P - K0
        res = ID + P + e2 + e3
        if self._verbose>1:
            print("twok={0}".format(twok))
            print("K0={0}".format(K0))
            print("ep={0}".format(ep))
            print("e2={0}".format(e2))
            print("e3={0}".format(e3))
            print("P={0}".format(P))
        if cuspidal==0:
            res = res + K0
        return res   #,ep




### For the dimension formula

@cached_function
def class_nr_pos_def_qf(D):
    r"""
    Compute the class number of positive definite quadratic forms.
    For fundamental discriminants this is the class number of Q(sqrt(D)),
    otherwise it is computed using: Cohen 'A course in Computational Algebraic Number Theory', p. 233
    """
    if D>0:
        return 0
    D4 = D % 4
    if D4 == 3 or D4==2:
        return 0
    K = QuadraticField(D)
    if is_fundamental_discriminant(D):
        return K.class_number()
    else:
        D0 = K.discriminant()
        Df = ZZ(D).divide_knowing_divisible_by(D0)
        if not is_square(Df):
            raise ArithmeticError("Did not get a discrimimant * square! D={0} disc(D)={1}".format(D,D0))
        D2 = sqrt(Df)
        h0 = QuadraticField(D0).class_number()
        w0 = _get_w(D0)
        w = _get_w(D)
        #print "w,w0=",w,w0
        #print "h0=",h0
        h = 1
        for p in prime_divisors(D2):
            h = QQ(h)*(1-kronecker(D0,p)/QQ(p))
        #print "h=",h
        #print "fak=",
        h=QQ(h*h0*D2*w)/QQ(w0)
        return h
    
def _get_w(D):
    """
    The number of roots of unity in a quadratic field with discriminant D
    :param D:
    :return:
    """
    if D == -3:
        return 6
    if D == -4:
        return 4
    else:
        return 2

def red_mod1(z):
    """
    Compute z mod 1 in [0,1)
    :param z:
    :return:
    """
    if z > 0:
        return z - floor(z)
