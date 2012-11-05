r"""

Algorithms for holomorphic and non-holomorphic Poincare series.

"""

from sage.all import ComplexField,inverse_mod
from multiplier_systems import *
from mysubgroup import *
from poincare_series_alg import *

class PoincareSeries(SageObject):
    def __init__(self,G,k=0,prec=53,multiplier=None,verbose=0,holomorphic=True):
        r"""
        Setup the associated space.
        """
        if isinstance(G,(int,Integer)):
            self._N=G
            self._group=MySubgroup(Gamma0(G))
        else:
            self._group=G
            self._N=G.level()
        self._k = k
        self._verbose = verbose
        self._holomorphic = holomorphic
        self._prec=prec
        self._RF=RealField(prec)
        self._CF=ComplexField(prec)
        
    def K(self,m,n,c):
        r"""
        K(m,n,c) = sum_{d (c)} e((md+n\bar{d})/c)
        """
        summa=0
        z=CyclotomicField(c).gen()
        print "z=",z
        for d in range(c):
            if gcd(d,c)>1:
                continue
            try:
                dbar = inverse_mod(d,c)
            except ZeroDivisionError:
                print "c=",c
                print "d=",d
                raise ZeroDivisionError
            arg=m*dbar+n*d
            #print "arg=",arg
            summa=summa+z**arg
        return summa

    

    def C(self,m,n,Nmax=50):
        k = self._k
        if n>0:
            if self._holomorphic:
                res=Aplus_triv(m,k,n,self._N,self._prec,Nmax,verbose=self._verbose)
            else:
                res=Bplus_triv(m,k,n,self._N,self._prec,Nmax,verbose)
        elif n==0:
            if self._holomorphic:
                res=Azero_triv(m,k,n,self._N,self._prec,Nmax,self._verbose)
            else:
                res=Bzero_triv(m,k,n,self._N,self._prec,Nmax,self._verbose)
        else:
            if self._holomorphic:
                res=Aminus_triv(m,k,n,self._N,self._prec,Nmax,self._verbose)
            else:
                res=Bminus_triv(m,k,n,self._N,self._prec,Nmax,self._verbose)
            
        return res
    
    def B0plus(self,m,n,N=10):        
        k = self._k
        res=Bplus_triv(m,k,n,self._N,self._prec,N)
        #if n==0:
        #    summa=0
        #    for c in range(1,N):
        #        cn=c*self._N
        #        #term=self._CF(self.K(m,n,cn))/self._RF(cn)**2
        #        term=Ktriv(self._N,m,n,cn,self._prec)/self._RF(cn)**2
        #        summa=summa+term
        #f=self._RF(4)*self._RF.pi()**2
        #return f*summa
        return res

    def B0plus_old(self,m,n,k=2,N=10):
        if n==0:
            summa=self._CF(0)
            for c in range(1,N):
                cn=c*self._N
                #term=self._CF(self.K(m,n,cn))/self._RF(cn)**2
                term=self.K(m,n,cn) #    Ktriv(self._N,m,n,cn,self._prec)/self._RF(cn)**2
                #print "K,m,n(",cn,")=",term
                term=self._CF(term)/self._RF(cn)**k
                #print "term(",cn,")=",term
                summa=summa+term
        f=self._RF(2**k)*self._RF.pi()**k
        f=f*self._RF(m)**(k-1)
        f=f*self._CF(0,1)**(self._RF(k+2))
        f=f/self._RF(factorial(k-1))
        return f*summa


    def Aminus_triv(self,m,n):
        k = self._k
        N = self._N
        if n==0:
            summa=self._CF(0)
            for c in range(1,N):
                cn=c*self._N
                #term=self._CF(self.K(m,n,cn))/self._RF(cn)**2
                term=self.K(m,n,cn) #    Ktriv(self._N,m,n,cn,self._prec)/self._RF(cn)**2
                #print "K,m,n(",cn,")=",term
                term=self._CF(term)/self._RF(cn)**k
                #print "term(",cn,")=",term
                summa=summa+term
        f=self._RF(2**k)*self._RF.pi()**k
        f=f*self._RF(m)**(k-1)
        f=f*self._CF(0,1)**(self._RF(k+2))
        f=f/self._RF(factorial(k-1))
        return f*summa
