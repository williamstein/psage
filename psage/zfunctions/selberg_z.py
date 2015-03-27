# -*- coding: utf-8 -*-
r"""
Selberg Zeta function and transfer operator for Hecke triangle groups.

This file contains the classes SelbergZeta and TransferOperator and methods to compute the values of the Selberg Z-function on Hecke triangle groups using the transfer operator method.
For an explanation of the algorithm see F. Str\"omberg, 'Computation of Selberg zeta functions on Hecke triangle groups', Arxiv:...

AUTHORS:

 - Fredrik Stroemberg (2012-05-15)


 
NOTE: If used with the implicit error bounds (i.e. iterate until the desired precision is reached) the methods might take *very* long time.
 


EXAMPLES::


sage: Z=SelbergZeta(3,verbose=0,working_prec=249,delta=1e-7)
sage: st8=Z.make_table_phi(prec=249,M0=150,N=2,outprec=66)




"""
#*****************************************************************************
#  Copyright (C) 2012 Fredrik Strömberg <fredrik314@gmail.com>
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


import mpmath as mpmath
from sage.all import  Parent,RR,ZZ,QQ,is_even,matrix,zeta,is_odd,is_even,ceil,log_b,log,gamma,tan,cos,sin,latex,CyclotomicField,MatrixSpace,sign,binomial
from selberg_z_alg import *
from psage.rings import mpc_pochammer
from sage.rings.complex_mpc import MPComplexField
from sage.rings.complex_field import ComplexField
from psage.matrix.matrix_complex_dense import Matrix_complex_dense
from psage.modform.maass.maass_forms import scattering_determinant_sl2z,scattering_determinant_Hecke_triangle

from sage.misc.sage_timeit import sage_timeit
from sage.parallel.decorate import *

def SelbergZ(q,verbose=0,working_prec=103,digits=5,delta=1e-5,**kwds):
    return SelbergZeta(q,verbose,working_prec,digits,delta,**kwds)


class SelbergZeta(Parent):
    r"""
    The Selberg zeta function of a Fuchsian group.
    So far, only implemented for Hecke triangle groups.

    EXAMPLE::

      sage: F=MPComplexNumbers(103)
      sage: R=F(9.5336952613535575543442352359287703238212563951072519823757904641353489912983477817692555099754353664930447678582858545070606)
      sage: s=R+F(0.5)

    """
    
    def __init__(self,q,verbose=0,working_prec=103,digits=5,delta=1e-5,**kwds):
        r""" Initalize the Selberg zeta function for the Hecke triangle  group G_q. 

        INPUT::

        - `q` -- integer (parameter for the Hecke triangle group)
        - `working_prec` -- integer. The number of bits for initial working precision. This might need to be increased to reach desired precision.)
        - `digits` -- integer. The number of digits we want the values to be accurate to. 
        - `delta` -- positive real number. The maximum distance of eigenvalues treated as the same eigenvalue.
        - `verbose` -- integer. Set the verbosity of self.
        """
        self._q=q
        T=TransferOperator(q,verbose=verbose,prec=working_prec)
        self._transfer_operator=T
        self._MAXIT=10000
        self._verbose=verbose
        self._working_prec=working_prec
        self._delta=delta
        self._eps = 10.0**(-digits)
        #if self._eps < 2.0**(1-prec):
        #    print "Not high enough precision to obtain eps={0}".format(self._eps)
        if self._delta < 2.0**(1-working_prec):
            print "Not high enough precision to obtain delta={0}".format(self._delta)
            self._delta = 2.0**(3-working_prec)
        self._deltad=ceil(abs(log_b(self._eps,10)))
        #self._epsd=digits #ceil(abs(log_b(self._eps,10)))
        
    def __repr__(self):
        r"""  Return string representation of self. """
        s="Selberg Zeta function of the Hecke triange group G_%s."%(self._q)
        return s


    def __call__(self,s,**kwds):
        r"""
        Compute the value of Z(s).
        
        INPUT:
          - `s`  -- complex number

          - `M0` -- integer : the starting value of the rank of the approximation.
          - `Nh` -- integer: the starting value of the differens in ranks of the approximations used.
          - `get_digits` -- integer: set to 0 if you just want to run the algorithm once and se which precision you get.
                                    Otherwise set to D if you want to keep iterating until D digits of precision is estimated.
          - `get_err` -- integer : set to 0 if you just want to return the value and no error estimates in the case of get_digits>0
        """
        return self.value(s,**kwds)
    
    def value(self,s,N=0,Nh=0,sym=1,checks=1,get_digits=0, get_eps=0,get_err=1,ret_param=0,prec=0,approx=0,A1_in=None,A2_in=None,get_evs=0,verbose=0):
        r"""
        Compute the value of Z(s).

        INPUT::
        
        - `s`  -- complex number
        - `N` -- integer (default: computed). The starting value of the rank of the approximation.
        - `Nh` -- integer (default 3). The starting value of the differens in ranks of the approximations used.
        - `get_digits` -- integer (default 1). Possible values
            0 -- if you just want to run the algorithm once and se which precision you get.
            D -- integer >0. If you want to keep iterating until D digits of precision is estimated.
        - get_eps -- real. Alternative to get_digits. Specifiy the precsion you want.
        - `get_err` -- integer : set to 0 if you just want to return the value and no error estimates in the case of get_digits>0
        - `checks` -- integer. Possible values:       
            - 1. Return the tuple:
               z,eps
            - 2. Return the tuple:
               z,k,eps,delta
            - >=3: Return the tuple:
               z,k,eps,delta,err
            Here:
             - z =Z(s)
             - k = number of eigenvalues used
             - eps = abs. value of smallest used eigenvalue
             - delta = max relative error in the used eigenvalues
            - err =  |phi_Z(s)-phi_E(s)|
              where phi_Z is the scattering determinant, computed using Z
              and phi_E is the same, but computed using either the explicit
              formula, if q=3, or the Eisenstein series otherwise.
        - `ret_param` -- integer. Set to 1 if you want to return the N and precision used in a second tuple.
        -`get_evs` -- return a list of verified eigenvalues of the transfer operator of self.
        """
        l=3
        self.set_prec(s,prec)
        verbose=max(self._verbose,verbose)
        current_prec=self._working_prec
        CCF = ComplexField(current_prec)
        CF = MPComplexField(current_prec)
        s0 = CCF(s.real(),s.imag())
        s  = CF(s0.real(),s0.imag())
        s0=s
        if get_digits>0 or get_eps>0:
            if get_eps==0:
                eps=10.0**(-get_digits)
            else:
                eps=get_eps
            if verbose>0:
                print "eps=",eps
        else:
            eps = 1.0
        delta_eps = 1
        err = 1
        if N>0:
            M00 = N
        else:
            M00 = ceil(abs(s.imag())*1.5)
        if Nh>0:
            Nh0=Nh
        else:
            Nh0=10
        M=M00
        err_old = 1
        M_new=0
        if verbose>0:
            print "l=",l

        if A1_in<>None or A2_in<>None:
            A_tmp = self._transfer_operator.matrix_approximation(s,1)
        else:
            A_tmp = 0
        for mj in range(1000):
            #self._working_prec=current_prec
            #s_hp = MPComplexField(self._working_prec)(s.real(),s.imag())
            if A1_in<>None and type(A1_in)==type(A_tmp):
                A1 = A1_in
            elif approx==3:
                A1=self._transfer_operator._approximation3(s,M)
            else:
                if sym==1:
                    A1plus=self._transfer_operator.matrix_approximation(s,M,sym=sym,eps=1)
                    A1minus=self._transfer_operator.matrix_approximation(s,M,sym=sym,eps=-1)
                else:
                    A1=self._transfer_operator.matrix_approximation(s,M,sym=sym)

            if A2_in<>None and type(A2_in)==type(A_tmp):
                A2 = A2_in
            elif approx==3:
                A2=self._transfer_operator._approximation3(s,M+l)
            else:
                if sym==1:
                    A2plus=self._transfer_operator.matrix_approximation(s,M+l,sym=sym,eps=1)
                    A2minus=self._transfer_operator.matrix_approximation(s,M+l,sym=sym,eps=-1)
                else:
                    A2=self._transfer_operator.matrix_approximation(s,M+l,sym=sym)
            if sym==1:
                ev1plus=A1plus.eigenvalues(sorted=1)
                ev1minus=A1minus.eigenvalues(sorted=1)
                ev2plus=A2plus.eigenvalues(sorted=1)
                ev2minus=A2minus.eigenvalues(sorted=1)
                ev1 = ev1plus; ev1.extend(ev1minus)
                ev2 = ev2plus; ev2.extend(ev2minus)
            else:
                ev1=A1.eigenvalues(sorted=1)
                ev2=A2.eigenvalues(sorted=1)
            if verbose>2:
                print "ev1=",ev1
                print "ev2=",ev2
            if checks>=1:
                evs,delta_eps = get_close_values(ev1,ev2,self._delta,ret_err=1)
            else:
                evs = get_close_values(ev1,ev2,self._delta,ret_err=0,verbose=verbose)
            if verbose>1:
                print "evs=",evs
            if  verbose>0:
                print "len(evs)=",len(evs)
            z=CF(1)
            for d in evs:
               z*=(1-d)
            if verbose>0:
                print "M=",M
            # We know that the "true" eigenvalues decrease to 0
            # Hence the error can be estimated by the smallest
            # "true" eigenvalue
            if len(evs)>0:
                err = min(map(abs,evs))
            else:
                err=1
            if verbose>0:
                print "Error estimate={0}".format(err)
                print "Desired eps={0}".format(eps)
                print "det(1-L)=",z
                ## Also print the quotient
                k=self._transfer_operator.det_1_minus_K(s,prec=current_prec)
                print "det(1-L)/det(1-K)=",z/k
            if get_evs==1:
                return evs
            if err<eps or (get_digits==0 and get_eps==0):
                if verbose>0:
                    print "exit loop: Error={0} less than eps={1}".format(err,eps)
                break
            ## We now have to decide if we need to raise precision
            ## To choose a precision and appriximation size is the most difficult part of the algorithm
            inc_prec=0; inc_appr=0
            ## The first time we always increase the approximation
            ## so we can get a better estimate of the error
            if err_old==1:
                inc_appr=1; err_old = err
                M_old = M
            elif M_new>0:
                ## If we have already tried to increase the approximation
                ## to an estimated 'good' value without success
                ## we try to increase the precision
                inc_prec=1
                if verbose>0:
                    print "inc_prec=",inc_prec
            else:
                r_test = exp(log(err_old/err)/(M_old-M))
                C_test = err / r_test**M
                if verbose>0:
                    print "r_test=",r_test
                    print "C_test=",C_test
                if r_test<1:
                    M_new = ceil(abs(log(C_test/eps)/log(r_test)))
                if M_new > 1.5*M:
                    M_new = M

                if verbose>0:
                    print "M_new=",M_new
                if M_new > M:
                    Nh0=M_new-M
                    inc_appr=0
                else:
                    inc_prec=1
                if len(evs)<self._deltad:
                    inc_prec=1
                #if (mj % 2) == 1: ## Every 2 times we also increase 
                #    inc_appr=1
            if inc_appr==0 and inc_prec==0:
                inc_appr=1; Nh0=5
            if inc_appr==1:
                if verbose>0:
                    print "Adding to M! M={0} Mplus={1}".format(M,Nh0)
                M += Nh0
            if inc_prec==1:
                current_prec+=20
                CF = MPComplexField(current_prec)
                if verbose>0:
                    print "raising number of bits to {0}".format(current_prec)
                s = CF(s0)
                if verbose>0:
                    print "j={0}".format(mj)
        if verbose>0:
            if mj>=999:
                print "warning: error={0}".format(err)
            print "det(1-L)=",z
        k=self._transfer_operator.det_1_minus_K(s,prec=current_prec)
        if verbose>0:
            print "det(1-K)=",k
            print "get_err=",get_err
            print "prec=",prec
            print "current_prec=",current_prec
            if hasattr(s,"prec"):
                print "s.prec()=",s.prec()
        z = z/k
        if checks==1:
            if ret_param==0:
                return z,len(evs),float(err),float(delta_eps)
            else:
                return (z,len(evs),float(err),float(delta_eps)),(M,current_prec)
        elif checks>=2: #get_err>0:
            #if get_err==1:  # Use error estimate from |lambda_K|
            #    return z,float(err)
            #else:
            ## Use scattering determinant
            if s.real()<>0.5:
                if verbose>0:
                    print "computing L(1-s) with M=",M
                z2 =self.value(1-s,N=M,Nh=Nh,sym=sym,checks=0,get_digits=get_digits,get_eps=get_eps,get_err=0,prec=prec,approx=approx)[0]
            else:
                z2 = z.conjugate()
            p = self._psi(s,prec=prec)
            ph = z2/z/p
            if verbose>0:
                print "psi=",p
                print "ph=",ph
            if self._q in [3,4,6]:
                p1=scattering_determinant_Hecke_triangle(s,self._q,prec=current_prec)
            else:
                ndig = max(16,ceil(abs(log_b(err,10))))
                if get_digits>0:
                    ndig = min(get_digits,ndig)
                if get_eps>0:
                    ndig = min(ndig,ceil(abs(log_b(get_eps,10))))
                pprec=dps_to_prec(ndig)
                if verbose>0:
                    print "Compute Eisen with prec=",pprec
                p1=scattering_determinant_Hecke_triangle(s,self._q,prec=pprec,use_eisenstein=1)
            if verbose>0:
                print "phi=",p1
            er = abs(p1-ph)
            if verbose>0:
                print "er=",er
            lev = len(evs)
            err = float(RR(err))
            delta_eps = float(RR(delta_eps))
            er = float(RR(er))
            if ret_param==0:
                return (z,lev,err,delta_eps,er)
            else:
                return (z,lev,err,delta_eps,er),(M,current_prec)
        elif get_err==0 and get_digits>0:
            if ret_param==0:
                return z
            else:
                return (z,lev,err,delta_eps,er),(M,current_prec)
        else:
            if ret_param==0:
                return z,float(err)
            else:
                return (z,float(err)),(M,current_prec)

    def scattering_determinant(self,s,N=0,Nh=0,sym=1,checks=1,get_digits=0, get_err=1,ret_param=0,prec=0,approx=0,outprec=0,verbose=0):
#prec=0,M0_start=0,Nh_start=0,do_once=0,checks=0,outprec=0):
        r"""
        Use the functional equation of Z(s) to compute the scattering determinant.

        If checks=1 we return the tuple:
        [phi,k,eps,delta]
        where phi=phi(s)
        - k = number of eigenvalues used
        - eps = abs. value of smallest used eigenvalue
        - delta = max relative error in the used eigenvalues

        """
        old_prec = self._working_prec        
        if self._working_prec<prec:
            self._working_prec=prec
        CF = ComplexField(self._working_prec)
        s = CF(s)
        l = self.value(s,N=N,Nh=Nh,sym=sym,checks=1,get_digits=get_digits, get_err=get_err,ret_param=ret_param,prec=prec,approx=approx,verbose=verbose)
        if ret_param==1:
            M,pprec = l[1]
            l = l[0]
            if verbose>0:
                print "M,pprec=",M,pprec
        if verbose>0:
            print "self.value=",l
        p = self._psi(s,prec=self._working_prec)
        if verbose>0:
            print "self._psi={0} with prec={1}".format(p,p.prec())
        if checks>=1:
            z=l[0]; K=l[1];er=l[2]; delta=l[3]
        elif get_err<>0:
            z = l[0]
        else:
            z=l
        if s.real()<>0.5:
            if verbose>0:
                print "Computing L(1-s) with N=",N
            z2 =  self.value(CF(1)-s,N=N,Nh=Nh,sym=sym,checks=checks,get_digits=get_digits, get_err=ge_err,ret_param=0,prec=prec,approx=approx,verbose=verbose)
            #self.value(1-s,M0=M0_start,Mh=Mh_start,checks=0)
        else:
            z2 = z.conjugate()
        if verbose>0:
            print "working_prec=",self._working_prec
        if z.prec()<>self._working_prec:
            z = CF(z.real(),z.imag())
        if z2.prec()<>self._working_prec:
            z2 = CF(z2.real(),z2.imag())
        ph = z2/z/p
        ph = CF(ph.real(),ph.imag())
        self._working_prec = old_prec        
        if checks==0:
            if ret_param==0:
                return ph
            else:
                return ph,(M,pprec)
        else:
            if self._q in [3,4,6]:
                p1=scattering_determinant_Hecke_triangle(s,self._q,prec=prec)
            else:
                p1=0
            if verbose>0:
                print "p1=",p1
            RF = RealField(3)
            er1=RF(abs(p1-ph))

            ndig = max(16,ceil(abs(log_b(er1,10))))
            if get_digits>0:
                ndig = min(get_digits,ndig)
            eisen_prec=dps_to_prec(ndig)
            if eisen_prec<53:
                eisen_prec=53 ## We want to use at least double precision
            if verbose>0:
                print "Need {0} digits.".format(ndig)
                print "Compute Eisen with prec=",eisen_prec
            p2 = scattering_determinant_Hecke_triangle(s,self._q,prec=eisen_prec,use_eisenstein=1)
            if verbose>0:
                print "p2=",p2
            er2=RF(abs(p2-ph))
            er =RF(er)
            delta = RF(delta)
            if outprec>0:
                CF = ComplexField(outprec)
                ph = CF(ph.real(),ph.imag())
            if ret_param==0:
                return ph,K,er,delta,er1,er2
            else:
                return (ph,K,er,delta,er1,er2),(M,pprec)

    def set_prec(self,s=None,prec=0):
        if prec==0:            
            if hasattr(s,"prec"):
                prec = s.prec()
            else:  ## If e.g. s is symbolic s= 1 + I
                prec = self._working_prec
        self._working_prec=prec

    def Theta(self,t,branch=0):
        r"""
        The argument of Z(1/2+it)
        branch -- 0 : Use principal branch, 1: cut at positive axis
        """
        self.set_prec(t)
        CF = ComplexField(self._working_prec)
        RF = CF.base_ring()
        s = CF(0.5,t)
        z = scattering_determinant_Hecke_triangle(s,self._q,use_eisenstein=0)
        z = z*self._psi(s)
        ar = z.argument()
        if branch <> 0:
            x = z.real()
            y = z.imag()
            if y<0:
                ar = ar + RF(2)*RF.pi()
        return ar

    def RealZ(self,t,M0,Mh=3,do_once=1,checks=0,branch=0):
        self.set_prec(t)
        CF = ComplexField(self._working_prec)
        s = CF(0.5,t)
        z=self.value(s,M0_start=M0,Mh_start=Mh,do_once=do_once,checks=checks)
        th = self.Theta(t,branch)
        z = z*CF(0,0.5*th).exp()
        return z


    def values_in_range(self,t0=5,t1=10,N=100,M0=0,prec=0):
        r"""
        Compute the values of Z(1/2+it) for t in a specific range.

        INPUT:

        - 
        """

        if prec>0: self.set_prec(prec,None)
            
        CF = ComplexField(self._working_prec)
        RF = CF.base_ring()
        t0 = RF(t0); t1=RF(t1)
        h = (t1-t0)/RF(N)
        M0 = 0
        zold = CF(1)
        res=[]
        branch = 0
        for i in range(N):
            t = t0+h*i
            znew,M1 = self.RealZ(t,M0,branch=branch)
            res.append(znew)
            if i==0:
                zold=znew
                continue
            # if abs(znew)<0.5*h:
            #     ## We are close to a zero and keep the current branch
            #     branch = 0
            # if abs(zold)<0.5*h:
            #     ## We might have passed through a zero and if so use same branch
            # slope = (znew-zold).real()/h
            
            # if abs(zold-znew)>h:
            zold = znew
        return res

    def make_table_phi(self,prec=0,N=50,ls=1,lf=10,target_dig=0,get_times=0,outprec=63,verbose=0):
        r"""
        Produce a LaTeX table of values and error estimates of self.

        INPUT:

        - ``prec`` -- integer. Set the working precision.
        - `N` -- integer. The approximation size to use.
        - `ls`  -- integer
        - `lf` -- integer
        - `target_dig` -- integer. If >0 we run the algorithm until target_dig (decimal) digits of precision is found
        - `outprec` -- integer (the precision to use for the output)


        OUTPUT:

        A string.

        """
        st=""
        precold = self._working_prec
        if prec>0:
            self._working_prec=prec
        row_args=[]
        rows={}
        for l in range(ls,lf+1):
            row_args.append((l,prec,N,target_dig,get_times,outprec,verbose))
        
            #stl = self.make_table_row(l,prec=prec,N=N,target_dig=target_dig,get_times=get_times,outprec=outprec,verbose=verbose)
        rows = sorted(list(self.make_table_row(row_args)))
        st=""
        print "rows=",rows
        for row in rows:
            if verbose>0:
                print "row=",row
            if len(row[1])>1:
                l,stl = row[1]
            else:
                l = row[1]
            if verbose>0:
                print "l=",l
                print "stl=",stl
            st+=stl
        self._working_prec=precold
        return st

    #    @fork(verbose=True)
    #@parallel(p_iter='fork',ncpus=4)
    def make_table_row(self,l,prec=0,N=50,target_dig=0,get_times=0,outprec=63,verbose=0):
        r"""
        Make one row of the table.
        """
        CF = ComplexField(self._working_prec)
        if isinstance(l,list) and isinstance(l[0],tuple):
            ll,prec,N,target_dig,get_times,outprec,verbose=l[0]
            l = ll
        #if verbose>0:
        print "l=",l
        s = CF(0.5,l)
        if verbose>0:
            print "Making table row for s={0} with N={1} and prec={2}".format(s,N,self._working_prec)
        if target_dig==0:
            ll = self.scattering_determinant(s,N=N,checks=3,outprec=outprec,verbose=verbose-1)
        else:
            ll = self.scattering_determinant(s,N=N,checks=3,get_digits=target_dig,outprec=outprec,ret_param=1,verbose=verbose-1)
            M,pprec = ll[1]
            ll = ll[0]
            if M>N:                    
                N = M
                if verbose>0:
                    print "Increased: Got M ={0} and pprec= {1}.".format(M,pprec)
            if pprec>self._working_prec:
                self._working_prec=pprec
            pprec = latex(pprec)
            if get_times:
                if verbose>0:
                    print "ROW: {0}".format(ll)
                    print "We are now recomputing in order to get the time."
                globals_dict = globals()
                globals_dict['self']=self
                globals_dict['s1']=ComplexField(self._working_prec)(0.5,l)
                globals_dict['N']=N
                globals_dict['outprec']=outprec
                stime = sage_timeit('self.scattering_determinant(s1,N=N,checks=0,outprec=outprec,verbose=0)',
                                    globals_dict=globals_dict,seconds=True,number=1,repeat=1)
            else:
                stime=""
        z,K,er,delta,er1,er2=map(latex,ll)
        if target_dig==0:
            if self._q in [3,4,6]:
                stl="${0}$ & ${1}$ & ${2}$ & ${3}$ & ${4}$ & ${5}$ & ${6}$\\\\ \n".format(l,z,er1,er2,K,er,delta)
            else:
                stl="${0}$ & ${1}$ & ${2}$ & ${3}$ & ${4}$ & ${5}$ & ${6}$\\\\ \n".format(l,z,er2,K,er,delta)
        else:
            if self._q in [3,4,6]:
                stl="${0}$ & ${1}$ & ${2}$ & ${3}$ & ${4}$ & ${5}$ & ${6}$ & ${7}$ & ${8}$ & ${9}$\\\\ \n".format(l,z,er1,er2,K,er,delta,M,pprec,stime)
            else:
                stl="${0}$ & ${1}$ & ${2}$ & ${3}$ & ${4}$ & ${5}$ & ${6}$  & ${7}$ & ${8}$ \\\\ \n".format(l,z,er2,K,er,delta,M,pprec,stime)
        if verbose>0:
            print "stl=",stl
        return l,stl
                    





    def make_table_phi2(self,s,prec_list=[53,83,166,249],N1=10,N2=70,outprec=63,time=1,verbose=0):
        r"""
        Produce a LaTeX table of values and error estimates of self for precision and N in ranges and fixed s.

        INPUT:

        - ``prec`` -- integer. Set the working precision.
        - `N` -- integer. The approximation size to use.
        - `ls`  -- integer
        - `lf` -- integer
        - `outprec` -- integer (the precision to use for the output)


        OUTPUT:

        A string.

        """        
        st=" s = {0} \n".format(s)
        precold = self._working_prec
        for p in prec_list:
            if p <2:
                continue
            self._working_prec = p 
            s1 = ComplexField(p)(s.real(),s.imag())
            for N in range(N1,N2+1,10):
                if verbose>0:
                    print "Calculate for prec={0} and N={1}".format(p,N)
                if N == N1:
                    p0 = p
                else:
                    p0 =""
                ll = self.scattering_determinant(s1,N=N,checks=3,outprec=outprec)
                z,K,er,delta,er1,er2=map(latex,ll)
                if time==1:
                    ## Do a timeit for the proces without checks
                    globals_dict = globals()
                    globals_dict['self']=self
                    globals_dict['s1']=s1
                    globals_dict['N']=N
                    globals_dict['outprec']=outprec
                    stime = sage_timeit('self.scattering_determinant(s1,N=N,checks=0,outprec=outprec,verbose=0)',
                                        globals_dict=globals_dict,seconds=True)
                    #{'self':self,'s1':s1,'N':N,'outprec':outprec},seconds=True)
                else:
                    stime = ""
                if self._q in [3,4,6]:
                    st+="${0}$ & ${1}$ & ${2}$ & ${3}$ & ${4}$ & ${5}$ & ${6}$ & ${7}$ \\\\ \n".format(p0,N,z,er1,K,er,delta,stime)
                else:
                    st+="${0}$ & ${1}$ & ${2}$ & ${3}$ & ${4}$ & ${5}$ & ${6}$ & ${7}$ \\\\ \n".format(p0,N,z,er2,K,er,delta,stime)
        self._working_prec=precold
        return st

    

        
    def test_functional_equation(self,s,prec=0,N=0,Nh=0,do_once=0):
        r"""
        Test the functional equation of self for s--> 1 - s

        """
        lhs = self.scattering_determinant(s,prec,N,Nh)
        phi = scattering_determinant_Hecke_triangle(self._q,s,prec)

        if self._verbose>0:
            print "Z(s)=",z
            print "phi(Zeta)(s)=",lhs
            print "phi(s)=",phi
        return abs(lhs-phi)

            
            

    def _psi(self,s,prec=0):
        r"""  Compute the factor Psi(s)=Psi(1/2)Z(1-s)/Z(s) in the Functional equation for the Selberg zeta function   (note that the factor Phi(1/2)=1 or -1) see e.g. Hejhal STF vol. 2: p. 499-500.

        INPUT::

        - `s` -- complex number
        - `prec` -- integer (working precision)

        OUTPUT::

        - `Psi(s)=+-1*Z(1-s)/Z(s)`
        
        """
        if hasattr(s,"prec"):
            prec = s.prec()
        elif prec>0:
            prec = prec
        else:
            prec = self._working_prec
        mpmath.mp.prec = prec
        CF = ComplexField(prec)
        ss=CF(s.real()-0.5,s.imag()); sigma=ss.real(); T=ss.imag()
        s1=CF(1.5-s.real(),-s.imag()); s2=CF(0.5+s.real(),s.imag())
        fak1=mpmath.mp.gamma(s1)/mpmath.mp.gamma(s2)
        mpi = mpmath.mp.mpc(0,1)
        mp1 = mpmath.mp.mpf(1)
        mppi = mpmath.mp.pi
        mppii = mpmath.mp.pi*mpi
        twopi = mppi*mpmath.mp.mpf(2)
        twopii = mppi*mpi*mpmath.mp.mpf(2)
        A=mpmath.mp.mpf(self._q-2)/mpmath.mp.mpf(self._q)*mppi
        f1 = lambda y: y*mpmath.mp.tan(mppii*y)
        IH1=-mpmath.mp.quad(f1,[0,T])
        f2 = lambda x: mpmath.mp.mpc(x,T)*mpmath.mp.tan(mppi*mpmath.mp.mpc(x,T))
        IH2=mpmath.mp.quad(f2,[0,sigma])
        H1=-A*(IH1+IH2)
        f3 = lambda y : mpmath.mp.cos(mppi*mpi*y)**-1
        IE1=mpi*mpmath.mp.quad(f3,[0,T])
        f4 = lambda x: mpmath.mp.cos(mppi*mpmath.mp.mpc(x,T))**-1
        IE2=  mpmath.mp.quad(f4,[0,sigma])
        m=self._q
        E1=mppi*(IE1+IE2)/mpmath.mp.mpf(2)
        if self._verbose>2:
            print "prec=",prec
            print "mpmath.mp.dps=",mpmath.mp.dps
            print "mpmath.mp.prec=",mpmath.mp.prec
            print "IH1=",IH1
            print "IH2=",IH2
            print "IE1=",IE1
            print "IE2=",IE2
            print "E1=",E1
        for k in range(1,m): #from 1 to m-1 do:
            km = mpmath.mp.mpf(k)/mpmath.mp.mpf(m)
            g1 = lambda t: mpmath.mp.exp(twopi*km*t)  / (mp1+mpmath.mp.exp(twopi*t ))+mpmath.mp.exp(-twopi*km*t)/(mp1+mpmath.mp.exp(-twopi*t))
            IE11 = mpi*mpmath.mp.quad(g1,[0,T])
            g2 = lambda x: mpmath.mp.exp(-twopii*km*mpmath.mp.mpc(x,T))  / (mp1+mpmath.mp.exp(-twopii*mpmath.mp.mpc(x,T) ))+mpmath.mp.exp(twopii*km*mpmath.mp.mpc(x,T))/(mp1+mpmath.mp.exp(twopii*mpmath.mp.mpc(x,T)))

            IE12 = mpmath.mp.quad(g2,[0,sigma])
            E1=E1+mppi*(IE11+IE12)/mpmath.mp.mpf(m)/mpmath.mp.sin(mppi*km)
            if self._verbose>2:
                print "E1[",k,"]=",E1
        P1=CF(1-2*s.real(),-2*s.imag())*mpmath.mp.ln(mpmath.mp.mpf(2))
        P=fak1*mpmath.mp.exp(H1+E1+P1)
        if self._verbose>2:
            print "P1=",P1
            print "E1=",E1
            print "H1=",H1            
        return ComplexField(prec)(P.real,P.imag)


    def RealZ(self,t,M0=0,Mh=3,do_once=0,branch=0,prec=0,tol=0.01,get_err=0):
        r"""  Compute the rotated Selberg zeta function, i.e. Z(1/2+it)*exp(-iArgZ(1/2+it)).   """
        if self._verbose>0:
            print "RealZ at t=",t
        if hasattr(t,"imag"):            
            if t.imag()<>0:
                raise ValueError,"Call with real argument only! Got:{0}".format(t)
            t = t.real()
        self.set_prec(t,prec)
        CF = ComplexField(self._working_prec)
        s =CF(0.5,t)
        z = self.value(s,M0_start=M0,Mh_start=Mh,do_once=do_once,get_err=get_err)
        th= self.Theta(t,branch)
        if get_err>=1:
            er = z[1]
            z = z[0]
        res = z*CF(0,0.5*th).exp()
        if abs(res.imag())>tol:
            print "Warning: function not real:{0}".format(res)
        if get_err==0:
            return res
        else:
            return res,er

    def theta(t,branch=0):
        r"""
        The argument of Z(1/2+it).

        INPUT:

        - `t` -- real number
        - `branch` -- choice of branch cut
        
        """
        prec  = t.prec()
        CF = ComplexField(prec)
        s = CF(0.5,t)
        z=scattering_determinant_hecke_triangle(self._q,s)*self._psi(s,q)
        ar=arg(z)
        if branch<>0: # Do not use principal branch 
            # Use branch cut at positive axis instead 
            x=z.real()
            y=z.imag()
            if y<0:
                ar=ar+CF(2)*CF.base_ring().pi()
        return ar

    def values_on_interval(t1,t2,step_size=0.05,M0=0,tol=1E-4,filename=None):
        r"""
        Computes the Selberg zeta function along an interval
        and writes the values into a file

        """
        CF = ComplexField(self._working_prec)
        RF = RealField(self._working_prec)
        if filename==None:
            filename='selberg_zeta-q{0}-{1}--{2}.txt'.format(self._q,t1,t2)
        if M0==0:
            M0=20
        if t2<t1:
            raise ValueError,"Need t2 > t1!"
        if step_size==0:
            if ns==0:
            
                raise ValueError,"Need either number of steps or stepsize!"
            step_size = RF(t2-t1)/RF(ns)
        if ns==0:
            ns = ceil(RF(t2-t1)/step_size)
        fp = open(filename)
        z_old = CF(1)
        zl=[]; tl=[]
        for i in range(ns):
            t = t1+i*step_size
            if z_old.real()>=0:
                branch=0
            else:
                branch=1
            tl.append(t)
            z,M = self.Realz(t,M0_start=M0,branch=branch,return_M0=1)
            zl.append(z)
            fp.write((t,z))
            M0=M
            if i==0:
                continue
            if i==1:
                if z.real()>0:
                    branch = 0
                else:
                    branch = 1
                continue
            ### Determine (experimentally) the branch to use
            if abs(z)<0.5*h: # We are probably close to or just passed a zero
                slope = (zl[i]-zl[i-2])/(t[i]-t[i-2])
                const = zl[i-1] - slope*tl[i-1]
                app_val = slope*t+ const
                if app_val.real()*z.real()<0:
                    oldsign = -1
                    z = sign(oldsign)*z









class TransferOperator(Parent):
    r"""
    Transfer Operator for the geodesic flow on the Hecke triangle group G_q.
    """

    def __init__(self,q,prec=53,verbose=0):
        r""" Initializes the transfer operator for the geodesic flow on G_q.  """
        self._q=q
        self._is_set=False
        self._prec=prec
        #self._dps=mpmath.mp.dps
        self._verbose=verbose 
        self._R=0
        self._lambdaq=0
        self._lambdaq_r=0
        self._setup_transfer_operator()
        self._MAXIT=100000
        self._contraction_factor = 0

    def __repr__(self):
        r""" Return string representation of self.  """
        s="Transfer operator of Hecke triangle group G_{0}\n".format(self._q)
        #s+=" lambda="+str(self._lambdaq)
        #s+="\n R="+str(self._R)
        #s+="\n h="+str(self._h)
        #s+="\n dim="+str(self._dim)
        return s

    def _setup_transfer_operator(self,prec=0):
        r"""
        Setup the necessary data for the transfer operator.
        """
        
        #ctx=self._ctx
        if self._is_set and (prec==0 or prec>0 and prec<self._prec):
            return
        if prec>0:
            self._prec=prec
        q=self._q
        RF = RealField(self._prec)
        mp1=RF(1); mp2=RF(2); mp4=RF(4)
        qr=RF(q)
        llambda=mp2*(RF.pi()/qr).cos()
        dim=0
        if is_even(q):
            #R=mp1
            h=ZZ(QQ(q-2)/QQ(2))
            dim=ZZ(2*h)
            #            print "dim=",dim
            NIJ=matrix(ZZ,dim)
            NIJ[0,h-1]=2
            for  i in range(1,h):
                NIJ[i,i-1]=1
                NIJ[i,h-1]=2
            for i in range(h,2*h):
                NIJ[i,h-1]=1
            NIJ[2*h-1,h-1]=1
            for j in range(dim):
                for i in range(h):
                    #print "lhs=",j,dim-1-i,"rhs=",dim-1-j,i
                    NIJ[j,dim-1-i]=-NIJ[dim-1-j,i]
        elif q==3:
            R=(RF(5).sqrt()-mp1)/mp2
            dim=ZZ(2)
            h=ZZ(0)
            NIJ=matrix([[3,-2],[2,-3]])
        elif q>=5:
            #R=llambda/mp2-mp1
            #R+=mp1/mp2*((mp2-llambda**2)+mp4).sqrt()
            h=ZZ(QQ(q-3)/QQ(2))
            dim=ZZ(4*h+2)
            NIJ=matrix(ZZ,dim)	
            NIJ[0,2*h-1]=2
            NIJ[0,2*h]=3
            NIJ[1,2*h]=2
            if self._verbose>0:
                print "h=",h
                print "dim=",dim
            for i in range(2,2*h+2):
                NIJ[i,i-2]=1
                NIJ[i,2*h]=2
            for i in range(2*h+2,4*h+1):
                NIJ[i,2*h-1]=1
                NIJ[i,2*h]=2
            for i in range(4*h+1,4*h+2):
                if self._verbose>1:
                    print i,2*h-1,2*h
                NIJ[i,2*h-1]=1
                NIJ[i,2*h]=2
            for j in range(0,dim):
                for i in range(0,2*h+1):
                    #print "lhs=",j,dim-i,"rhs=",dim-j,i
                    NIJ[j,dim-1-i]=-NIJ[dim-1-j,i]
        self._is_set=True
        self._h=h
        self._lambdaq_r=llambda
        self._R=self.R(prec=prec)
        if self._q<>3:
            K = CyclotomicField(2*self._q)
            z = K.gens()[0]
            llambda = z+z.conjugate()
            if self._verbose>0:
                print "lambda=",llambda
                print "RR(lambda)=",llambda.complex_embedding()
                print "RR(z)=",z.complex_embedding()
            if not llambda.is_real_positive():
                if not (-llambda.is_real_positive()):
                    raise ArithmeticError,"Could not get lambda!"
                else:
                    llambda = -llambda
        else:
            llambda=1
        if self._verbose>0:
            print "lambda=",llambda
        self._lambdaq=llambda
        self._dim=dim
        self._Nij=NIJ
        if is_odd(self._q):
            self._numi = 2*self._h+1
        else:
            self._numi = self._h


    def R(self,prec=0,format='float',verbose=0):
        if prec>0:
            RF = RealField(prec)
        else:
            prec = self._prec
            RF = RealField(self._prec)
        if format=='float':
            if self._R<>0 and self._prec>=prec:
                return RF(self._R)
            if is_even(self._q):
                self._R=RF(1)
            elif self._q==3:
                self._R=RF(RF(5).sqrt()-1)/RF(2)
            else:
                lambdaq=self.lambdaq_r(prec)
                self._R = RF((2-lambdaq)**2+4).sqrt()/RF(2) + lambdaq/RF(2)-RF(1)
                if verbose>0:
                    print "R1=",RF((2-lambdaq)**2+4).sqrt()/RF(2) + lambdaq/RF(2)-RF(1)
                    print "R2=",-RF((2-lambdaq)**2+4).sqrt()/RF(2) + lambdaq/RF(2)-RF(1)
            return self._R
        elif format=='alg':            
            if self._R_alg<>0:
                return self._R_alg
            if is_even(self._q):
                self._R_alg = 1
            #elif self._q==3:
            #    
            #    self._R_alg=(ZZ(5).sqrt()-1)/RF(2)
            else:
                lambdaq=self.lambdaq()
                F = lambdaq.parent()
                z = F['z'].gens()[0]
                f = z^2+(2-llambda)*z+1
                print "f=",f
                #self._R = RF((2-lambdaq)**2+4).sqrt()/RF(2) + lambdaq/RF(2)-RF(1)
            
    def lambdaq(self):
        return self._lambdaq

    def lambdaq_r(self,prec=0):
        if prec>self._prec or self._lambdaq_r==0:
            if hasattr(self._lambdaq,"complex_embedding"):
                self._lambdaq_r=self._lambdaq.complex_embedding(prec).real()
            else:
                self._lambdaq_r=RealField(prec)(self._lambdaq)
        return self._lambdaq_r
        
    
    def numi(self):
        return self._numi

    def dim(self):
        return self._dim

    def h(self):
        return self._h

    def Nij(self):
        return self._Nij
    

    def trace(self,s,eps=1e-10,imax=10000):
        r"""
        Compute the trace of self at s.
        """
        trace = 0
        for m in self._Nij.diagonal():
            if m==0:
                continue
            if m<0:
                sg=-1
            else:
                sg=1
            tmptrace = 0
            for n in range(abs(m),imax):
                tmp = self.trace_ns(n*sg,s)
                tmptrace+=tmp
                if abs(tmp)<eps and n>10:
                    break
            if self._verbose>0:
                print "number of iterations: {0}".format(n)
                print "Trace[{0}]={1}".format(m,tmptrace)
                print "Lasterr=",abs(tmp)
            trace+=tmptrace
        return trace

    def trace_ns(self,n,s):
        prec = s.prec()
        nn = self.norm_n(n,prec)**-1
        one = nn.parent()(1)
        return nn**s/(one-nn)
        

#    @cached_method
    def norm_n(self,n,prec=None):
        r"""
        Returns the trace of pi_s(ST^n) = N(ST^n)^-s/(1-N(ST^n)^-1)
        """
        if prec<>None:
            l = self.lambdaq_r(prec)
        else:
            l = self.lambdaq_r()
        nl = n*l
        two = l.parent()(2)
        a = nl**2-two
        b = abs(nl)*(nl**2-two**2).sqrt()
        x = (a+b)/two
        if x<1:
            if self._verbose>0:
                print "n=",n
                print "x=",x
            x = x - b
        assert x>1
        return x
    
    def get_markov_partition(self,itype='alg',prec=0,iformat='ie'):
        r"""
        Get the Markov partition (or rather half of it)

        - format -- set to 'ie' if you want a list of end points x_i
                    set to 'is' if you want a list of [x_i,x_{i+1}]
                    set to 'cr' if you want a list of [c_i,r_i]
                    
        """
        ## First get the algebraic version of lambda
        llambda_2 = self.lambdaq()/QQ(2)
        iphi=[llambda_2]
        numi = self.numi()
        if self._verbose>1:
            print "num int=",numi
            print "lambda/2=",llambda_2
            if hasattr(llambda_2,"complex_embedding"):
                print "lambda/2=",llambda_2.complex_embedding()
        y = llambda_2
        for i in range(numi):
            if self._verbose>1:
                print "y=",y
            y = self.fq(y)
            iphi.append(y)
        if self._verbose>1:
            print "iphi=",iphi
        if itype=='float':
            res = []
            if prec==0:
                prec = self._prec
            RF = RealField(self._prec)
            for x in iphi:
                if hasattr(x,"complex_embedding"):
                    x = x.complex_embedding(prec).real()
                elif hasattr(x,"real"):
                    x = x.real()
                else:
                    x = RF(x)
                res.append(x)
            res.sort()
        else:
            iphi.sort(cmp=my_alg_sort)
            res = iphi
        if iformat=='ie':
            return res
        elif iformat=='cr':
            res2 = []
            for i in range(len(res)-1):
                x1=res[i]; x2=res[i+1]
                if iformat=='is':
                    res2.append((x1,x2))
                else:
                    c = (x1+x2)/2
                    r = abs(x2-x1)/2
                    res2.append((c,r))
            return res2
        elif iformat=='ar':
            alphas={}; rhos={}
            i = 0 
            res.reverse()
            for i in range(len(res)-1):
                x1=res[i]; x2=res[i+1]
                c = (x1+x2)/2
                r = abs(x2-x1)/2
                alphas[i]=-c
                alphas[i+self._dim/2]=c
                rhos[i]=r
                rhos[i+self._dim/2]=r
                i+=1
            return alphas,rhos


    def get_contracting_discs(self,iformat='ar',verbose=0):
        r"""
        Compute intervals of contraction as in Lemma 4.4. of [MMS2012]





        REFERENCES:

            [MMS2012]   Mayer, M\"uhlenbruch, Str\"omberg, 'The Transfer Operator for the Hecke Triangle Groups', Discrete Contin. Dyn. Syst, Vol. 32, No. 7, 2012. 
        """
        intervals={}
        if self._q==3:
            intervals={1:[-1,QQ(1)/QQ(2)],2:[-QQ(1)/QQ(2),1]}
        elif self._q==4:
            intervals={1:[-1,self.lambdaq_r()/4],2:[-self.lambdaq_r()/4,1]}
        elif is_odd(self._q):
            dim_2=ZZ(self._dim).divide_knowing_divisible_by(ZZ(2))
            llambdaq_4=self.lambdaq_r()/4
            for i in range(self._h+1):
                l2ip1=[-1]
                for j in range(i):
                    l2ip1.append(-1)
                l2ip1.append(-2)
                for j in range(self._h):
                    l2ip1.append(-1)
                if verbose>0:
                    print "cf[2*{0}+1]={1}".format(i,l2ip1)
                x0 = self.cont_frac_to_pt(l2ip1)
                intervals[2*i+1]=[x0,llambdaq_4]
                if verbose>0:
                    print "intervals[{0}]={1}".format(2*i+1,intervals[2*i+1])
            for i in range(1,self._h+1):
                l2i=[-1]
                for j in range(i):
                    l2i.append(-1)
                x0 = self.cont_frac_to_pt(l2i)
                intervals[2*i]=[x0,llambdaq_4]
                if verbose>0:
                    print "intervals[{0}]={1}".format(2*i,intervals[2*i])

            for i in range(1,self._h+1):
                x1 = intervals[2*i][0]; x2 = intervals[2*i][1]
                if verbose>1:
                    print "x1,x2[{0}]={1}".format(2*i,(x1,x2))
                intervals[2*i+dim_2]=[-x2,-x1]
            for i in range(self._h+1):
                x1 = intervals[2*i+1][0]; x2 = intervals[2*i+1][1]
                if verbose>1:
                    print "x1,x2[{0}]={1}".format(2*i+1,(x1,x2))
                intervals[2*i+1+dim_2]=[-x2,-x1]

        else:
            raise NotImplementedError
        if verbose>1:
            print "intervals=",intervals
            print "intervals.keys()=",intervals.keys()
        if iformat=='ar':
            alphas={}
            rhos={}
            for i in range(1,len(intervals)+1):
                x1 = intervals[i][0];x2 = intervals[i][1]
                if verbose>1:
                    print "x1,x2[{0}]={1}".format(i,(x1,x2))
                c = (x1+x2)/2
                r = abs(x1-x2)/2
                alphas[i-1]=c
                rhos[i-1]=r
            return alphas,rhos
        return intervals


    def get_contracting_factor(self):
        r"""
        Returns the max quotient of the radii of the contracting disks
        
        """
        if self._contraction_factor<>0:
            return self._contraction_factor
        ## Check the contracting properties so that we are ok.
        d = self._dim
        #centers={}
        ad={}; bd={}
        centers,r1s = self.get_contracting_discs()
        mina={}; maxb = {}
        for i in range(d):
            mina[i]=0;maxb[i]=0
        ## First check that the original intervals contain the attractive and no repelling fixed-points as well as no singularities.
        for i in range(d):        
            aa = centers[i] - r1s[i]
            bb = centers[i] + r1s[i]
            for j in range(d):
                n = self._Nij[i,j]
                if n >= aa and n <= bb:
                    print "singularity at x={0} in [{1},{2}]".format(n,aa,bb)
                t =  self.STn_fixed_pts(n)
                if isinstance(t,tuple):
                    x1,x2 = t
                else:
                    x1=t; x2 = t
                if  x2 >= aa and x2 <= bb:
                    print "Interval contains repelling fixed point! x={0} in [{1},{2}]".format(x2,aa,bb)
                if not x1 > aa and x1 < bb and x1<>x2:
                    print "Interval does not contain attractive fixed point! x={0} in [{1},{2}]".format(x1,aa,bb)    
        if self._verbose>0:
            print "c=",centers
            print "r=",r1s
        for i in range(d):        
            for j in range(d):
                n = Nij[i,j]
                aa = centers[i] - r1s[i]
                bb = centers[i] + r1s[i]
                phia,phib=image_of_interval(T,aa,bb,n)
                if min(phia,phib) < mina[j]:
                    mina[j] = min(phia,phib)
                if max(phia,phib) > maxb[j]:
                    maxb[j] = max(phia,phib)            
                if self._verbose>0:
                    print "ST^[n>={n}]({a},{b})=[{x},{y}]".                 format(n=n,a=aa,b=bb,x=mina[j],y=maxb[j])
        if self._verbose>0:
            print "mina=",mina
            print "maxb=",maxb
            for j in range(d):            
                print "Phi(I_{j}) \subset [{a},{b}]".format(j=j,a=mina[j],b=maxb[j])
        # We now need to find disks with the same center as before
        radii={};rho={}
        for j in range(d):
            tmpr1 = abs(mina[j]-centers[j])
            tmpr2 = abs(maxb[j]-centers[j])
            radii[j]=max(tmpr1,tmpr2)
            if self._verbose>0:
                print "Disk[{0}]=D({1},{2})".format(j,centers[j],radii[j])
        r3s={}; r2s={};rhos={};r4s={}
        rho1s={};rho2s={};rho3s={}
        for i in range(d):
            r3s[i] = radii[i]  ## Get r3
            r2s[i] =(r3s[i]*r1s[i]).sqrt()  # and r2
            rho1s[i]=r2s[i]/r1s[i]
            rho2s[i]=r3s[i]/r2s[i]
            ## Find the optimal r4:
            r4s[i]=0
            rtmp = 0        
            for j in range(d):
                n = self._Nij()[i,j]
                a = centers[i]-r2s[i]; b = centers[i]+r2s[i]
                a1 = self.STn(a,n) ; b1 = self.STn(b,n) 
                rtmp1 = abs(a1-centers[j])
                rtmp2 = abs(b1-centers[j])
                rtmp = max(rtmp1,rtmp2)
                if rtmp>r4s[i]:
                    r4s[i]=rtmp+0.000001
            rho3s[i]=r4s[i]/r3s[i]
        rho1 = max(rho1s.values())
        rho2 = max(rho2s.values())
        rho3 = max(rho3s.values())
        rho = max([rho1,rho2,rho3])
        print "r4=",r4s
        print "r3=",r3s
        print "r2=",r2s
        # Should check that this gives a valid r4
        t =check_mapping_into(T,centers,r2s,centers,r4s,verbose=verbose)
        print "maps into=",t
        self._contraction_factor=rho
        return rho
        
        
    def check_mapping_into(self,c1,r1,c2,r2):
        r"""
        Check that the system of functions in T maps the disks given by c1,r1
        inside the disks of c2,r2, i.e. Phi^{ij}(D_i) in D_j
        """
        d = self._dim
        assert len(c1.keys())==d; assert len(c2.keys())==d
        assert len(r1.keys())==d; assert len(r2.keys())==d    
        for i in range(d):
            a1 = c1[i]-r1[i]; b1 = c1[i]+r1[i]
            if self._verbose>0:
                print "a1=",a1
                print "b1=",b1
            for j in range(d):
                n = T.Nij()[i,j]
                a2 = c2[j]-r2[j]; b2 = c2[j]+r2[j]
                a,b = image_of_interval(T,a1,b1,n)
                if self._verbose>0:
                    print "F^{i}{j}({a1},{b1}) = ({a2},{b2})".format(i=i,j=j,a1=a1,b1=b1,a2=a,b2=b)
                    print "I_{j}=({a},{b})".format(j=j,a=a2,b=b2)
                if a < a2 or b > b2:
                    return False
        return True

    
    def cont_frac_to_pt(self,cf=[],prec=0,verbose=0):

        r"""
        Compute the point corresponding to a nearest lambda continued fraction.

        INPUT:

         - `cf` list. A nearest continued fraction expansion given as a list with format:
                cf = [a_0,a_1,...] where a_0 is the integer part and a_j are integers.
         - `prec` -- integer. The bits of precision of the returned value.
         - `verbose` -- integer. Set the verbosity.

         OUTPUT:

         - real number of precision prec.
        """
        if prec==0:
            prec = self._prec 
        RF = RealField(prec)
        frac_part = cf[1:]; n = len(frac_part)
        x = RF(0)
        verbose = max(verbose,self._verbose)
        for j in range(n):
            x1 = self.STn(x,frac_part[n-j-1])
            if verbose>1:
                print "ST^{0}({1})={2}".format(frac_part[n-j-1],x,x1)
            x = x1
        x = x+self.lambdaq_r(prec)*RF(cf[0])
        return x

    def nearest_lambda_code(self,x,N=0):
        r""" Compute the nearest lambda continued fraction of x up to N symbols if N>0.  """
        if hasattr(x,"prec"):
            prec = x.prec()
            lambdaq_r = self.lambdaq_r(prec)
        elif isinstance(x,Rational):
            prec = 0
            N = -1
            lambdaq = self.lambdaq()
            lambdaq_r = self.lambdaq_r()
        else:
            prec = 53
            lambdaq_r = self.lambdaq_r(prec)
        RF = RealField(prec)
        a0 = self.nearest_lambda_mult(z)
        res = [a0]
        if prec>0:
            x = x - a0*lambdaq_r
        else:
            x = x - a0*lambdaq
        j = 0
        while x<>0:
            x1,n =  self.fq(y,retn=1)
            x = x1
            res.append(n)
            if j>=N and N>0:
                break
            j+=1
        return res

                    
        
    def STn(self,x,n,prec=0):
        r""" Return ST^n_q(x)= -1/(x+n*lambda_q). """
        if prec==0:
            prec = self._prec
        RF = RealField(prec)
        return RF(-1)/RF(x+n*self.lambdaq_r(prec))

    def STn_fixed_pts(self,n,ret="f"):
        r"""
        Compute fixed points of
        ST^n: x -> -1/(x+n)
        """
        x1 = (-n+QQ(n*n-4).sqrt())/QQ(2)
        x2 = (-n-QQ(n*n-4).sqrt())/QQ(2)
        if ret=="f":
            x1 = RR(x1); x2 = RR(x2)
        t1 = abs(x1+n)**2
        t2 = abs(x2+n)**2
        if t1>1 and t2<1:  ## contracting fixed point
            return x1,x2
        if t1<1 and t2>1:  ## contracting fixed point
            return x2,x1
        if x1==x2:
            return x1
        
    def fq(self,y,retn=0):
        r""" Compute f_q(x)=-1/x - n*lambda_q where n = nearest lambda multiple to -1/x. """
        if y==0:
            return 0
        if abs(y) > self.lambdaq_r()/2:
            n = self.nearest_lambda_mult(y)
            y = y - n
        else:
            z = -y**-1
            n = self.nearest_lambda_mult(z)
            y = z - n*self.lambdaq_r()
        if self._verbose>1:
            print "-1/x={0}".format(z)
            print "-1/x-{0}*lambda={1}".format(n,y)
        if retn==1:
            return y,n
        return y
    
    def nearest_lambda_mult(self,x):
        r""" Compute the nearest lambda multiple.  """
        if hasattr(x,"prec"):
            prec = x.prec()
        else:
            prec = 53
        RF = RealField(prec)
        if hasattr(x,"complex_embedding"):
            if x<>x.conjugate():
                raise ArithmeticError,"Call only with real argument!"
            x = x.complex_embedding().real()
        elif hasattr(x,"real"):
            x = x.real()
        else:
            x = RF(x)
        if self._q<>3:
            x = x.real()/self.lambdaq_r(prec)

        y = x+RF(0.5)
        if self._verbose>1:
            print "x/lambda+1/2={0}".format(y)
        ## Need my own defined floor function
        return my_floor(float(y))
        
    def matrix_approximation(self,s,M,sym=1,it=1,eps=1,alpha_in={},rhos_in={}):
        r"""
        Compute the approximation to self at s by a finite rank (matrix) operator acting on the Banach space B of vector-valued holomorphic functions on a product of discs D x ... x D with D centered at 0.

        INPUT:

        - `s` -- complex number.
        - `M` -- integer. The size of the the finite rank approximation. (The point of truncation of the Power series representating functions in B).
        
        """
        if hasattr(s,"prec"):
            prec = s.prec()
        else:
            prec = 53
        CF = MPComplexField(prec)
        RF = CF.base()
        self._setup_transfer_operator(prec)
        dim=self._dim
        B = self._Nij
        ll = self.lambdaq_r(prec) #CF._base(self._lambdaq)
        Nmax = max(map(abs,self._Nij.list()))
        
        if self._verbose>0:
            print "Nmax = ",Nmax
            print "M=",M
            print "dim=",dim
        s = CF(s.real(),s.imag())
        if sym==0:
            ## We use unsymmetrized operator            
            MS=MatrixSpace(CF,dim*(M+1))
            A=Matrix_complex_dense(MS,0)
            trop_approximation(A,self._Nij,M,self._q, dim, Nmax, s,ll,verbose=verbose)
        else:
            ## We use symmetrized operator
            sym_dim=ZZ(dim).divide_knowing_divisible_by(ZZ(2))
            if self._verbose>0:
                print "sym_dim=",sym_dim                
            MS=MatrixSpace(CF,sym_dim*(M+1))
            A=Matrix_complex_dense(MS,0)
            ## The question is now which intervals to use.
            if it==0:
                alphas=[RF(0) for x in range(self._dim)]
                rhos=[RF(1) for x in range(self._dim)]
            elif it==1: ## Use discs from the Markov partition
                alphas,rhos=self.get_markov_partition(itype='float',iformat='ar')                
            else:
                alphas,rhos=self.get_contracting_discs(iformat='ar')
            if alpha_in or rhos_in:
                if isinstance(alphas_in,dict):
                    for j in alphas_in.keys():
                        alphas[j-1]=RF(alphas_in[j])
                elif isinstance(alphas_in,list):
                    for j in len(alphas_in):
                        alphas[j]=RF(alphas_in[j])
                else:
                    raise ValueError,"Got intervals in wrong format! alphas={0}".format(alphas_in)
                if isinstance(rhos_in,dict):
                    for j in rhos_in.keys():
                        rhos[j-1]=RF(rhos_in[j])
                elif isinstance(rhos_in,list):
                    for j in len(rhos_in):
                        rhos[j]=RF(rhos_in[j])
                else:
                    raise ValueError,"Got intervals in wrong format! rhos={0}".format(rhos_in)
                        
            if self._verbose>1:
                print "alphas=",alphas
                print "rhos=",rhos
            trop_approximation(A,self._Nij,M,self._q, dim, Nmax, s,ll,verbose=self._verbose,approx_type=3,eps=eps,alphas=alphas,rhos=rhos)
            
        return A


    def get_eigenvalues(self,s,N=0,h=3,sym=1,delta=1e-7):
        r"""
        Return the 'verified' eigenvalues of the matrix approximation of self.
        
        """
        evs=[]
        if sym==1:
            A1plus=self.matrix_approximation(s,N,sym=sym,eps=1)
            A1minus=self.matrix_approximation(s,N,sym=sym,eps=-1)
            A2plus=self.matrix_approximation(s,N+h,sym=sym,eps=1)
            A2minus=self.matrix_approximation(s,N+h,sym=sym,eps=-1)
            ev1p = A1plus.eigenvalues(sorted=1)
            ev1m = A1minus.eigenvalues(sorted=1)
            ev2p = A2plus.eigenvalues(sorted=1)
            ev2m = A2minus.eigenvalues(sorted=1)
            ev1 = ev1p; ev1.extend(ev1m)
            ev2 = ev2p; ev2.extend(ev2m)
        else:
            A1=self.matrix_approximation(s,N,sym=sym)
            A2=self.matrix_approximation(s,N+h,sym=sym)
            ev1=A1.eigenvalues(sorted=1)
            ev2=A2.eigenvalues(sorted=1)
        evs = get_close_values(ev1,ev2,delta,ret_err=0,verbose=self._verbose)
        return evs


    def get_approximation_Zagier_type(self,s,M,eps=1):
        prec = s.parent().prec()
        CF = MPComplexField(prec)
        CCF = ComplexField(prec)
        RF = RealField(prec)
        dim=self._dim
        MS = MatrixSpace(CF,dim/2*M)
        A = Matrix_complex_dense(MS,0)
        B = Matrix_complex_dense(MS,0)
        xj = {}
        if self.lambdaq()<>1:
            l = self.lambdaq_r(prec)
        else:
            l = 1
        mpmath.mp.prec = prec
        for i in range(M):
            arg = RF(2*(i+1)-1)*RF.pi()/RF(4*M)
            tmp = arg.sin()**2
            #xj[i] = (tmp-RF(0.5))*RF(l)
            xj[i] = -tmp
            #*RF(l)/RF(2)
        #return xj
        twos = RF(2)*s
        for i in range(dim/2):
            pos=1
            if i > dim/2-1:
                pos = -1
            for j in range(dim/2):
                n = abs(self._Nij[i,j])
                do_zeta=1
                if j < dim/2-1 or j > dim/2+1:
                    do_zeta = 0
                for r in range(M):
                    zarg = CCF(twos.real()+r,twos.imag())
                    ri = i*M+r
                    for k in range(M):
                        kj = j*M+k
                        if l<>1:
                            xj_l = xj[k]/l
                            f = l**-zarg
                        else:
                            f=1
                            xj_l = xj[k]
                        if do_zeta==1:
                            z1 = mpmath.mp.zeta(zarg,xj_l+n)*(-1)**r
                            z2 = mpmath.mp.zeta(zarg,-xj_l+n)
                            z  = z1+eps*z2
                            z  = z*f
                        else:
                            z1=((-1)**r)*(xj[k]+n*l)**-zarg
                            z2=(-xj[k]+n*l)**-zarg
                            z = z1+eps*z2
                        B[ri,kj]=CF(z.real,z.imag)
                        A[ri,kj]=xj[k]**(r) #+eps*(-xj[k])**(k)
        return A,B
                    
   


        
    
    def _approximation3(self,s,M,dprec=0,one_i=0,check=0,alphas_in={},rhos_in={},verbose=0):
        r"""
        Compute the approximation to L_s up to order M Using Taylor expansions centered at different discs 

        INPUT:

        - `s` -- complex number
        - `M` -- integer
        - `dprec` -- integer
        - `verbose` -- integer
        - `one_i` -- if =1 we use one common interval (i.e. the same as approximation2)
        - `check` -- include extra checks
        - `alphas_in` -- a list of centers of intervals
        - `rhos_in` -- a list of radii of intervals
        
        NOTE: Not cythonized. This is a very slow algorithm!

        """
        if hasattr(s,"prec"):
            prec = s.prec()
        else:
            prec = 53
        CF = MPComplexField(prec)
        CFF = ComplexField(prec)
        RF = CF.base()
        llambda=self.lambdaq_r(prec)
        maxprec=max(prec,dprec)
        verbose = max(verbose,self._verbose)
        if one_i==1:
            alphas=[RF(0) for x in range(self._dim)]
            #rhos=[llambda/RF(2) for x in range(self._dim)]
            rhos=[RF(1) for x in range(self._dim)]
        elif one_i==-1:
            alphas=[RF(0) for x in range(self._dim)]
            rhos=[llambda/RF(2) for x in range(self._dim)]
        else:
            alphas,rhos=self.get_markov_partition(itype='float',iformat='ar')
        if alphas_in<>{}:
            alphas={}
            for j in alphas_in.keys():
                alphas[j-1]=RF(alphas_in[j])
        if rhos_in<>{}:
            rhos={}
            for j in rhos_in.keys():
                rhos[j-1]=RF(rhos_in[j])            
        ## Checking that the order is ok, i.e. that the 
        if verbose>0:
            print "alphas=",alphas
            print "rhos=",rhos
        if dprec>maxprec:
            print "Warning: Can not get more than %s bits of precision!" %(CF.prec())
            print "Please redefine the Transfer operator!"
        twos=CF(2)*CF(s.real(),s.imag())
        dim=self._dim
        MS = MatrixSpace(CF,dim*(M+1))
        A=Matrix_complex_dense(MS,0) #CF,dim*(M+1))
        Z=dict();nr=dict()
        mpmath.mp.prec=prec
        for i in range(dim):
            Z[i]=dict()
        for l in range(0,2*M+2):
            nr[l]=RF(l)
            zarg = CFF(twos.real()+nr[l],twos.imag())
            if one_i==1:
                for i in range(dim):
                    z = zeta(zarg)
                    Z[i][l]=CF(z.real(),z.imag())
            else:
                for i in range(dim):
                    if alphas[i]<>0:
                        z = mpmath.mp.zeta(zarg,alphas[i]/llambda+1)
                    else:
                        z = mpmath.mp.zeta(zarg,1)
                    Z[i][l]=CF(z.real,z.imag)
                    ## Note: alphas[i]<0 if i<dim/2 and >0 if i>dim/2
                    ## and alphas[dim/2+i]=-alphas[i] for 0<=i<dim/2
        # we only need ctx.mpf(l) another time for 0<=l<=M
        for l in range(M+1,2*M+1):
            lr=CF(l)
        if verbose>0:
            print "dim=",dim
        #llambda=CF(self._lambdaq_r.real(),self._lambdaq_r.imag())
        for n in range(M+1):
            #pow=nr[n]+twos
            fakn=RF.factorial(n)**-1
            for k in range(M+1):
                sg=1
                if is_odd(n+k):
                    sg=-1
                for i in range(dim):
                    for j in range(dim):
                        B=RF(self._Nij[i,j])
                        ai = RF(alphas[i])*RF(sign(B))
                        aj = RF(alphas[j])*RF(sign(B))
                        ri = rhos[i]; rj=rhos[j]
                        abB=B.abs()
                        if verbose>1 and n==0 and i==0:
                            print "B,|B=",B,abB
                        if B==0:
                            AA=CFF(0)
                        else:
                            summa = CFF(0)
                            for l in range(k+1):
                                poc=mpc_pochammer(twos+l,n)
                                tmp = poc*CF(binomial(k,l))
                                tmp = tmp*ai**(k-l)
                                zarg = CF(twos.real()+l+n,twos.imag())
                                if 2*j < dim -2 or 2*j>dim:
                                    #z = mpmath.mp.zeta(twos+l+n,aj+llambda*abB)                                
                                    z = (aj+llambda*abB)**-zarg
                                else:
                                    if aj<0:
                                        jj = ZZ(j).mod(dim/2)
                                    else:
                                        jj = dim/2+ZZ(j).mod(dim/2)
                                    z = Z[jj][l+n]
                                    if aj==0: 
                                        ztmp=CF(1,0)                                    
                                        for a in range(2,abB):
                                            ztmp = ztmp + CF(a,0)**(-zarg)
                                    else:
                                        ztmp=CF(0)                                    
                                        for a in range(1,abB):
                                            ztmp += CF(aj/llambda+RF(a),0)**(-zarg)
                                            if verbose>1:
                                                print "ztmp+[{0}]={1}".format(a,CF(aj/llambda+RF(a),0)**(-zarg))
                                                print "aj/lambda+a=",CF(aj/llambda+RF(a))
                                    z = z - ztmp
                                    if check==1: ## Extra check that we computed the Hurwitz zeta func. correctly
                                        z1 = mpmath.mp.zeta(zarg,aj/llambda+abB)
                                        z1 = CF(z1.real,z1.imag)
                                        if verbose>1:
                                            print "diff=",abs(z-z1)
                                        if abs(z-z1)>1e-10:
                                            print "z0=Z[",jj,"][",l+n,"]=",Z[jj][l+n]
                                            print "z,z1=",z,z1
                                            raise ArithmeticError,"Problem with H-Z!"
                                    z = CF(z.real(),z.imag())
                                    ## test for lambda=1 z = z*llambda**(-zarg)                                    
                                tmp = tmp*z
                                summa+=tmp
                            AA=summa
                        fak=fakn*(ri**-k*rj**n) #*RF.factorial(k).sqrt()/RF.factorial(n).sqrt()   #*ri**-k)*(rj**n/RF.factorial(n))
                        if verbose>2:
                            print "fak[{0}][{1}]={2}".format(n,k,fak)
                        AA = AA *fak
                        AA = CF(AA.real(),AA.imag())
                        if B>0 and sg<0:
                            A[i*(M+1)+n,j*(M+1)+k]=-AA
                        else:
                            A[i*(M+1)+n,j*(M+1)+k]=AA
        return A


    

    

    def _sym_approximation(self,s,M,eps=1,dprec=0,interval_type=1,check=0,alphas_in={},rhos_in={},verbose=0):
        r"""
        Compute the approximation to L_s up to order M Using Taylor expansions centered at different discs 

        INPUT:

        - `s` -- complex number
        - `M` -- integer
        - `sign` -- integer +1 or -1
        - `dprec` -- integer
        - `verbose` -- integer
        - `one_i` -- if =1 we use one common interval (i.e. the same as approximation2)
        - `check` -- include extra checks
        - `alphas_in` -- a list of centers of intervals
        - `rhos_in` -- a list of radii of intervals
        
        NOTE: Not cythonized. This is a very slow algorithm!

        """
        if hasattr(s,"prec"):
            prec = s.prec()
        else:
            prec = 53
        assert eps in [1,-1]        
        CF = MPComplexField(prec)
        CFF = ComplexField(prec)
        RF = CF.base()
        llambda=self.lambdaq_r(prec)
        maxprec=max(prec,dprec)
        verbose = max(verbose,self._verbose)
        alphas={};rhos={}
        if interval_type==0:
            ### Use one interval centered at zero for all functions
            alphas=[RF(0) for x in range(self._dim)]
            rhos=[llambda/RF(2) for x in range(self._dim)]
            #rhos=[RF(1) for x in range(self._dim)]            
        elif interval_type==1:
            ## Use the markov partition.
            alphas,rhos=self.get_markov_partition(itype='float',iformat='ar')

            
        elif interval_type==2:
            alphas,rhos=self.get_contracting_discs(iformat='ar')
            # USe the contracting discs
        elif alphas_in=={} or rhos_in=={}:
            raise ValueError,"Must have valid format for intervals! Got alphas={0} and rhos={1}".format(alphas,rhos)
        ## If we set intervals we apply them now
        if alphas_in<>{}:
            alphas={}
            for j in alphas_in.keys():
                alphas[j]=RF(alphas_in[j])
        if rhos_in<>{}:
            rhos={}
            for j in rhos_in.keys():
                rhos[j]=RF(rhos_in[j])
        if isinstance(alphas,dict):
            alphas=list(alphas.values())
        if isinstance(rhos,dict):
            rhos=list(rhos.values())
        dim=self._dim
        ### Check that alpha and rhos have correct format
        if len(alphas)<>dim or len(rhos)<>dim:
            raise ValueError,"Must have valid format for intervals! Got alphas={0} and rhos={1}".format(alphas,rhos)
        for j in range(dim):            
            if j<dim/2 and alphas[j]>0 or j>=dim/2 and alphas[j]<0:
                raise ValueError,"Must have valid format for intervals! Got alphas={0} and rhos={1}".format(alphas,rhos)
        if verbose>0:
            print "alphas=",alphas
            print "rhos=",rhos
        if dprec>maxprec:
            print "Warning: Can not get more than %s bits of precision!" %(CF.prec())
            print "Please redefine the Transfer operator!"
        twos=CF(2)*CF(s.real(),s.imag())

        sym_dim=ZZ(dim).divide_knowing_divisible_by(ZZ(2))
        MS = MatrixSpace(CF,sym_dim*(M+1))
        A=Matrix_complex_dense(MS,0) 
        Z=dict();nr=dict()
        B={}; AA={}

        mpmath.mp.prec=prec
        #print "A.nrows()=",A.nrows()
        #print "A.cols()=",A.ncols()
        for i in range(dim):
            Z[i]=dict()
            for j in range(dim):
                Z[i][j]=dict()
        lpow={}
        for l in range(0,2*M+2):
            nr[l]=RF(l)
            zarg = CFF(twos.real()+nr[l],twos.imag())
            lpow[l]=llambda**-zarg
            for i in range(sym_dim):
                for j in range(sym_dim):
                    Z[i][j][l]={}
                    B[0] = RF(self._Nij[i,j])
                    B[1] = RF(self._Nij[i,dim-1-j]) 
                    assert B[0]>=0
                    assert B[1]<=0
                    if 2*j < dim -2 or 2*j>dim:
                       z = (alphas[i]/llambda+B[0])**-zarg
                       Z[i][j][l][0]=CF(z.real(),z.imag())
                       z = (-alphas[i]/llambda-B[1])**-zarg
                       Z[i][j][l][1]=CF(z.real(),z.imag())
                    else:
                        z = mpmath.mp.zeta(zarg,alphas[i]/llambda+B[0])
                        Z[i][j][l][0]=CF(z.real,z.imag)
                        z = mpmath.mp.zeta(zarg,-alphas[i]/llambda-B[1])
                        Z[i][j][l][1]=CF(z.real,z.imag)
                    #print i,j,l,0,CC(Z[i][j][l][0].real(),Z[i][j][l][0].imag())
                    #print i,j,l,1,CC(Z[i][j][l][1].real(),Z[i][j][l][1].imag())
        # we only need ctx.mpf(l) another time for 0<=l<=M
        for l in range(M+1,2*M+1):
            lr=CF(l)
        if verbose>0:
            print "dim=",dim
        for n in range(M+1):
            fakn=RF.factorial(n)**-1
            for k in range(M+1):
                sg=1
                if is_odd(n+k):
                    sg=-1
                for i in range(sym_dim):
                    for j in range(sym_dim):
                        ni = i*(M+1)+n; kj = j*(M+1)+k
                        ri = rhos[i]; rj=rhos[j]
                        ## Do the positive and negative term separately
                        B[0] = RF(self._Nij[i,j])
                        B[1] = RF(self._Nij[i,dim-1-j])
                        for ii in range(2):
                            if ii==0:
                                ai = RF(alphas[i]); aj = RF(alphas[j])
                            else:
                                ai = -RF(alphas[i]); aj = RF(alphas[j])
                            if B[ii]==0:
                                AA[ii]=CFF(0)
                            else:
                                summa=CFF(0)
                                for l in range(k+1):
                                    poc=mpc_pochammer(twos+l,n)
                                    tmp = poc*CF(binomial(k,l))
                                    z = Z[i][j][l+n][ii]
                                    #if ni==0 and kj==32:
                                    #    print "z[",i,j,l+n,ii,"]=",z
                                    tmp = tmp*z*lpow[l+n]*aj**(k-l)
                                    summa+=tmp
                                    #if ni==0 and kj==32:
                                    #    print "tmp[",i,j,l+n,ii,"]=",tmp*fakn
                                if B[ii]>0:
                                    AA[ii]=summa*sg
                                else:
                                    AA[ii]=summa
                            #if ni==0 and kj==32:
                            #    print 
                            if verbose>2:
                                print "fak[{0}][{1}]={2}".format(n,k,fak)
                            #AA[ii] = AA[ii] *fak
                            AA[ii] = CF(AA[ii].real(),AA[ii].imag())
                        fak=fakn*(rj**-k*ri**n) #*RF.factorial(k).sqrt                        
                        #if ni==0 and kj==32:
                        #    print "AA[0]=",AA[0]
                        #    print "AA[1]=",AA[1]             
                        A[ni,kj] = AA[0]+eps*AA[1]*(-1)**k
                        #if ni==0 and kj==32:
                        #    print "A[0,32]=",A[ni,kj]
                        #    print "fak=",fak
                        A[ni,kj] = A[ni,kj]*fak

        return A


    def det_1_minus_K(self,s,prec=0):
        r"""
        Compute the Fredholm determinant det(1-K_s) of the auxiliary operator K_s.

        INPUT:

        - `s` -- complex number
        - `prec` -- integer (default: 0) precision in bits. If prec=0 the precision is determined by the argument s.


        ALGORITHM:
        We know that each eigenvalue of K_s is of the form mu_n = L**(-2s-2n) where L is explicitly given.
        Hence det(1-K_s)=Prod(1-L**(-2s-2n)) where the product is truncated when the relative error is <2**-prec.
        
        """
        if hasattr(s,"prec"):
            if prec==0 or prec<s.prec():
                prec = s.prec()
        else:
            if prec==0:
                prec = 53
        self._prec=prec
        self._setup_transfer_operator(prec)
        RF=RealField(prec)
        CF=ComplexField(prec)
        eps = RF(2)**RF(-prec)
        ss = CF(s.real(),s.imag())
        mp1=RF(1);mp2=RF(2);mp4=RF(4)
        twos=mp2*ss
        llambda=self.lambdaq_r(prec)
        if self._q==3:
            L=mp2+self._R
        elif self._q==4:
            L=mp2.sqrt()+mp1
        elif(is_even(self._q)):
            L=(mp2+llambda)/(mp4-llambda*llambda).sqrt()
        else:
            L=(mp2+self._R*llambda)/(mp2-llambda)
        l=dict()
        K=mp1
        for n in range(self._MAXIT):
            pow=-twos-mp2*RF(n)
            #print "pow[",n,"]=",pow
            l[n]=(L**pow)
            K=K*(mp1-l[n])
            if n>1:
                err=max(abs(l[n]),abs(l[n-1]))*abs(K)
                #print n,l[n],err
                if err<eps:
                    break
        if n > self._MAXIT-2:
            raise ArithmeticError,"Could not obtain approximation in %s iteraations. Raise value of self._MAXIT?"
        #ctx.dps=old_prec
        return K

                    


def Gauss_approximation(s,M):
    r"""
    Try the approximation (for the Gauss map) as suggested by Zagier in 'New points of view on Selberg Zeta function' 2002.

    """
    RF = s.parent()
    CF = MPComplexField(RF.prec())
    xj={}
    mpmath.mp.prec = RF.prec()
    for j in range(1,M+1):
        arg = RF(2*j-1)*RF.pi()/RF(4*M)
        tmp = arg.sin()**2
        xj[j] = tmp
    MS = MatrixSpace(CF,M)
    A = Matrix_complex_dense(MS,0)
    B = Matrix_complex_dense(MS,0)
    twos=s*2
    for k in range(1,M+1):
        for j in range(1,M+1):
            A[k-1,j-1]=xj[j]**(k-1)
            #z = mpmath.mp.zeta(xj[j]+1,twos+k-1)
            z = mpmath.mp.zeta(twos+k-1,xj[j]+1)
            B[k-1,j-1]=CF(z.real,z.imag)
    C = A.inverse()*B
    return C

    

def eigenvalues_Gauss(s,M,h=3,delta=1e-7):
    r"""
    Get approximations to the eigenvalues of the Gauss transfer operator
    """
    A1 = self._approximation_GaussZ_one_mat(s,M)
    A2 = self._approximation_GaussZ_one_mat(s,M+h)
    ev1 = A1.eigenvalues(sorted=1)
    ev2 = A2.eigenvalues(sorted=1)
    if self._verbose>0:
        print "ev1=",ev1
        print "min(ev1-1)=",min(map(lambda x:abs(x-1),ev1))
        print "ev2=",ev2
        print "min(ev2-1)=",min(map(lambda x:abs(x-1),ev2))
    evs = get_close_values(ev1,ev2,delta,ret_err=0,verbose=self._verbose)
    return evs




def get_close_values(l1,l2,delta,test='rel',ret_err=0,verbose=0):
    r"""
    Take two lists of numbers and returns a list of those numbers in list nr. 2 which are closer than delta to a number in list 1.

    INPUT:
     - ``l1`` -- lista nr. 1
     - ``l2`` -- lista nr. 2
     - ``delta`` -- desired error
    - ``prec`` -- work with this precision, if elements of l1 do not have a precision.
    -``test`` -- is either 'rel' or 'abs' if we check relative, or absolute error
    """
    res=[];  used=[]
    assert isinstance(l1,list) and isinstance(l2,list) 
    assert hasattr(l1[0],"prec") and hasattr(l2[0],"prec") 
    prec1 = l1[0].prec();prec2 = l2[0].prec()
    prec = min(prec1,prec2)
    eps = 2.0**(1-prec)
    for i in range(len(l2)):
        used.append(0)
    if ret_err==1:
        delta_eps = 0
    for i in range(len(l1)):
        e1=l1[i]
        if abs(e1)<eps:
            continue
        for j in range(len(l2)):
            if used[j]==1:
                continue
            if test=='rel':
                err=abs(l1[i]-l2[j])/(abs(l1[i])+abs(l2[j]))
            elif test=='abs':
                err=abs(l1[i]-l2[j])
            else:
                raise ValueError,"test must be one of 'rel' or 'abs'!"
            if err<delta:
                if verbose>0:
                    print "(rel)err=",err
                    print "l1[{0}]={1}".format(i,l1[i])
                    print "l2=[{0}]={1}".format(j,l2[j])
                if ret_err==1:
                    if err > delta_eps:
                        delta_eps = err
                used[j]=1
                res.append(l2[j])
    if ret_err==1:
        return res,delta_eps
    else:
        return res
    
def dps_to_prec(dps):
    """ Convert number of digits to bits of precision. """
    return int(RR(dps*log(10)/log(2)))

def prec_to_dps(prec):
    """ Convert bits of precision to number of digits. """
    return int(RR(prec*log(2)/log(10)))


def my_alg_sort(x,y):
    r"""  Sort numbers in a number field."""
    if hasattr(x,"complex_embedding"):
        xx = x.complex_embedding().real()
    else:
        xx = x.real()
    if hasattr(y,"complex_embedding"):
        yy = y.complex_embedding().real()
    else:
        yy = y.real()
    return cmp(xx,yy)
