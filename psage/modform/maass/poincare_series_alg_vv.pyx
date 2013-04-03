#-*- coding: utf-8 -*-
r"""
Cython algorithms for computing holomorphic vector-valued Poincaré series.
In particular to find a basis and compute aq gram matrix consisting of the Fouriercoefficients p_{D_i,r_i}(D_j,r_j)

NOTE: At the moment this is only implemented for the simple Weil representation asociated to D=Z/NZ

NOTE: The current implementation doesn't really use Cython...


AUTHOR:
 - Fredrik Strömberg (May 2010)



 
"""


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

include 'sage/ext/stdsage.pxi'
include "sage/ext/cdefs.pxi"
include 'sage/ext/gmp.pxi'
include 'sage/ext/interrupt.pxi'
#include 'sage/structure/coerce.pxi'
from sage.modules.all import vector
import sage.structure.element
cimport sage.structure.element
#from sage.structure.element cimport Element, ModuleElement, RingElement
import operator

from sage.rings.integer_ring import ZZ
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.rings.complex_field import ComplexField
from sage.rings.number_field.number_field import CyclotomicField
from sage.rings.real_mpfr import RealField
from sage.rings.integer import Integer
from sage.rings.number_field.number_field import is_fundamental_discriminant
from sage.all import RealNumber as makeRealNumber
from sage.all import ComplexNumber as makeComplexNumber

from sage.functions.special import bessel_J
from sage.functions.other import gamma,ceil
#from sage.calculus.calculus import ceil
from sage.rings.arith import gcd,xgcd,kronecker,factor,odd_part,is_square,lcm
from sage.misc.functional import is_even

from sage.rings  cimport integer
from sage.rings.ring cimport *
from sage.rings.real_mpfr cimport *
from sage.rings.complex_number cimport *
from sage.libs.mpfr cimport *


global one, minus_one,two,zero,half,four,i,twopii,pii,CF,RF,digs, epsilon_0, twopi,psilent,eight,twopii_eight
#global inited,set_prec
inited=False
set_prec=0
#psilent=0


#cdef Class PoincareSeriesVV()
#def set_silence_level(s=0):
#    r"""
#    Set the level of debugging output from functions in this file.
#    """
#    global psilent
#    psilent=s

cpdef poincareseries_init(prec,verbose=0):
    r"""
    Initializing stuff for the computation of the Poincare series.
    """
    global one, minus_one,two,zero,i,twopii,pii,CF,RF,digs, epsilon_0, twopi
    global inited,set_prec,eight,four,silent,twopii_eight
#    print "want to init:",inited,set_prec,prec
    if(inited and prec==set_prec):
        return
    print "inited=",inited
    #if verbose<silent:
    #    set_silence_level(silent)
    CF=ComplexField(prec)
    RF=RealField(prec)
    one=RF(1.0)
    minus_one=RF(-1.0)
    two=RF(2.0)
    half=RF(0.5)
    four=RF(4.0)
    eight=RF(8.0)
    zero=RF.zero_element()
    i=CF((0,1))
    twopi=two*CF.pi()
    twopii=twopi*i
    pii=CF.pi()*i
    arg=RF(prec-1)/RF(3.5)
    digs=arg.floor()-1
    epsilon_0=RF(10)**RF(-digs)
    set_prec=prec
    twopii_eight=twopii/eight
    inited=True
    
    #Print "epsilon0=",epsilon_0


cdef weil_rep_matrix(int a_in,int b_in,int c_in,int d_in,int b_p, int b_m, int N,int prec=53):
    return 0




cpdef weil_kloosterman_sum(c,r_in,m,rp_in,n,N,kappa,b_m,b_p,prec,verbose=0):
    r"""
    Twisted Kloosterman sum
    """
    poincareseries_init(prec)
    summa=zero
    r=r_in % (2*N)
    rp=rp_in % (2*N)   
    for d in range(abs(c)):
        (g,x,y)=xgcd(d,c)        
        if(g==1):
            a=x
            b=-y
            R=weil_rep_matrix(a,b,c,d,b_p,b_m,N,prec)
            arg=twopii*RF(m*a+n*d)/RF(c)
            tmp=R[r,rp]*arg.exp()
            summa=summa+tmp
    #    H=summa*CF(exp(-pi*i*sgn(c)*kappa/2)/abs(c))
    if(c>0):
        arg=-pii*kappa/two
    else:
        arg=pii*kappa/two
    H=summa*arg.exp()/abs(c)
    return H



cpdef weil_kloosterman_sum_all_old(c,rs,mms,rps,nns,N,kappa,b_m,b_p,prec,verbose=0):
    r"""
    Twisted Kloosterman sum
    """
    # print c,r_in,m,rp_in,n,N,kappa,b_m,b_p,prec
    poincareseries_init(prec)
    H=dict()
    summa=dict()
    for j in rs.keys():
        summa[j]=zero
#    r=r_in % (2*N)
#    rp=rp_in % (2*N)   
    for d in range(abs(c)):
        #x,y=var('x,y')h
        #        if((gcd(d,abs(c))==1) or (d==0 and abs(c)==1)):
        (g,a,b)=xgcd(d,c)        
        if(g==1):
            b=-b
            R=weil_rep_matrix(a,b,c,d,b_p,b_m,N,prec)
            cr=RF(c)
            for j in rs.keys():
                arg=twopii*RF(mms[j]*a+nns[j]*d)/cr
                tmp=R[rs[j],rps[j]]*arg.exp()
                if verbose>1:
                   print "weil_rep_old[",a,b,c,d,"](",rs[j],rps[j],")=",R[rs[j],rps[j]]
#                   print "exp()=",arg.exp()
                summa[j]=summa[j]+tmp
            # endfor
        # end_if:
    # endfor:
    #    H=summa*CF(exp(-pi*i*sgn(c)*kappa/2)/abs(c))
    if(c>0):
        arg=-pii*RF(kappa)/two
    else:
        arg=pii*RF(kappa)/two
    tmp=arg.exp()/abs(RF(c))
    for j in rs.keys():
        H[j]=tmp*summa[j]
    if verbose>1:
        print "tmp=",tmp
        print "kloo sum=",H
    return H
#end_proc:

# x_c=1/2 iff even_part(D) || c
cpdef get_xc(D,c):
    if(c==1):
        return 0
    facD=ZZ(D).factor()
    facc=ZZ(c).factor()
    xc=0
    if(facD[0][0]==2  and facc[0][0]==2):
        if(facD[0][1] == facc[0][1]):
            xc=1
    return xc

cpdef weil_kloosterman_sum_all(c,rs,mms,rps,nns,N,kappa,b_m,b_p,prec,breaks,verbose=0):
    r"""
    Twisted Kloosterman sum
    """
    # print c,r_in,m,rp_in,n,N,kappa,b_m,b_p,prec
    poincareseries_init(prec)
    H=dict()
    summa=dict()
    for j in rs.keys():
        summa[j]=zero
    # Compute everything depending only on c
    D=ZZ(2*N)
#    xc=get_xc(D,c)
    oddD=odd_part(D)
    oddc=odd_part(c)
    evenc=c/oddc
    evenD=Integer(D/oddD)
    # Get xc
    if( evenc==evenD ):
        xc=1
    else:
        xc=0
    D_c=D_lower_c(D,c)
    gcdDc=gcd(Integer(c),D)
    gcdDc2=2*gcdDc
    cD=Integer(c/gcdDc)
    Dc=Integer(D/gcdDc)
    ## We compute the odddity and p-excess 
    odt=0
    fak=ZZ(D).factor()
    exc=0
    kron_c=1
    for [p,m] in fak:
        q=p**m
        gcdqc=gcd(q,c)
        q_c=ZZ(q).divide_knowing_divisible_by(gcdqc)
        if abs(q_c)>1:
            cc = ZZ(c).divide_knowing_divisible_by(gcdqc)
            kron_c=kron_c*cc.kronecker_symbol(q_c)
            if p>2:
                Dq = D.divide_knowing_divisible_by(q)
                if(not is_square(q_c) and Dq.kronecker_symbol(p)==-1):
                    exc=exc+q_c-1+4
                else:
                    exc=exc+q_c-1
            else:
                if(not is_square(q_c) and oddD.kronecker_symbol(2)==-1):
                    odt=oddD+4
                else:
                    odt=oddD
    cr=RF(c)
    lentmp=RF(len(D_c))/RF(D)
    lentmp=lentmp.sqrt()
    twopii_by_4N=twopii/RF(4*N)
    twopii_by_c=twopii/cr
    N4=RF(4*N)
    if verbose>2:
        print "RF=",RF
        print "CF=",CF
        print "odt=",odt
        print "exc=",exc
        print "kron_c=",kron_c
    odt=odt*oddc
    twoDc=evenc*gcdDc2
    od_m_oc_m_2=oddD-oddc-2
    for d in range(abs(c)):
        #x,y=var('x,y')
        #        if((gcd(d,abs(c))==1) or (d==0 and abs(c)==1)):
        (g,a,b)=xgcd(d,c)
        if(g==1):
            b=-b
            ac=a*c
            bd=b*d
            t=rget_t(D,a,b,c,d)
            if( (t % 4)  == 1):
                iarg=((t-1)*od_m_oc_m_2+odt*t+6-exc)
            else:
                iarg=((t-1)*od_m_oc_m_2+odt*t+4-exc)
            iarg=iarg % 8
            arg=twopii_eight*iarg
            kron=ZZ(t).kronecker_symbol(twoDc)
            kron=kron*ZZ(-d).kronecker_symbol(oddc)*kron_c

            chi=CF(kron*arg.exp())
            #chi=kron*tmp2*ep*tmp3*sqrt(len(Dc)/D)
            for j in rs.keys():            
                if(breaks[j]):  # We are already doe with this one
                    continue
                alpha=Integer(rs[j])
                beta=Integer(rps[j])
                c_div=False
                if(xc==0):
                    alpha_minus_dbeta=Integer((alpha-d*beta) % D)
                else:
                    alpha_minus_dbeta=Integer((alpha-d*beta-D/2) % D)
                for r in range(D):
                    if((r*c - alpha_minus_dbeta) % D ==0):
                        c_div=True
                        invers=Integer(r)
                        break
                if(c_div):
                    if(xc==0):
                        argu = a*c*invers**2-b*d*beta**2+2*b*beta*alpha
                    else:
                        argu=ac*invers**2+a*D*invers-bd*beta**2+b*2*beta*alpha
                    arg=twopii_by_4N*RF(argu)
                    tmp1=arg.exp()
                    weil_rep=tmp1*chi
                else:
                    weil_rep=zero
                    # R=weil_rep_matrix(a,b,c,d,b_p,b_m,N,prec)
                arg=twopii_by_c*(mms[j]*a+nns[j]*d)
                if(abs(weil_rep)>zero):
                    tmp=weil_rep*arg.exp()
                    summa[j]=summa[j]+tmp
            # endfor j
        # end_if:
    # endfor: d
    #    H=summa*CF(exp(-pi*i*sgn(c)*kappa/2)/abs(c))
    if(c>0):
        arg=-pii*RF(kappa)/two
    else:
        arg=pii*RF(kappa)/two
    tmp=arg.exp()/abs(cr)
    for j in rs.keys():
        H[j]=CF(tmp*summa[j]*lentmp)
    if verbose>2:
        print "tmp=",tmp
        print "lentmp=",lentmp
        print "kloo sum=",H
    return H
#end_proc:


# D^c = { gamma | gamma=c*delta, delta \in D}
def D_upper_c(D,c):
    """
    D_upper_c(D,c)
    INPUT:
    D = (Integer) Discriminant
    c = Integer
    OUTPUT:
    v = (List) list of k*c mod D for k mod D
    """
    v=list()
    for k in range(D):
        j=c*k % D
        if(v.count(j)==0):  # We only want to add new members
            v.append(j)
    v.sort()
    return v

# D_c = { gamma | gamma=c*delta, delta \in D}
def D_lower_c(D,c):
    """
    D_lower_c(D,c)
    INPUT:
    D = (Integer) Discriminant
    c = Integer
    OUTPUT:
    v = (List) List of k mod D such that c*k=0 mod D
    """
    v=list()
    for k in range(D):
        if(c*k % D==0):
            if(v.count(k)==0):  # We only want to add new members
                v.append(k)
    v.sort()
    return v

def oddity(D):
    r""" Calculates the oddity of the signature 1 lattice
        
    """
    Dodd=odd_part(D)
    Deven=D.divide_knowing_divisible_by(Dodd)
    if( (not(is_square(Deven))) and Dodd.kronecker(2)==-1):
        oddity=Dodd+4
    else:
        oddity=Dodd
    return oddity  

def eps_c_div_eps_m(D,a,b,c,d,t,verbose=0):
    fac=ZZ(D).factor()
    fak=one
    test_arg([D,a,b,c,d,t],ZZ)
    #psilent_bak=psilent
    #psilent=0
    t = ZZ(t)
    #t=get_t(D,a,b,c,d)
    for [p,j] in fac:
        q=ZZ(p**j)
        if( c % q <> 0):
            q_c = q.divide_knowing_divisible_by(gcd(q,c))
            #fak=fak*RF(kronecker(ZZ(t),ZZ(q/gcd(q,c))))
            fak=fak*RF(t.kronecker(q_c))
    if( c % Integer(D/odd_part(D)) <> 0):
        q=Integer(D/odd_part(D))
        arg=((t-1)*oddity(D)*Integer(c/gcd(q,c))) % 8
        arg=twopii*RF(arg)/RF(8)
        tmp=arg.exp()
        #tmp=exp(pi*I*arg/4)
        if verbose>0:
            print "2-adic arg=",(t-1),"*",c,"*",oddity(D),"/",gcd(q,c),"/",8,"=",arg
            print "2-adic part=",tmp
        fak=fak*tmp
    #psilent=psilent_bak
    return fak

def rget_t(D,a,b,c,d,verbose=0):
    N4=2*D
    f=ZZ(N4).factor()
    #    print "f=",f
    g=gcd(c,N4)
    #
    if verbose>1:
        print "(c,4N)=",g
    if(g==1):
        t=1
    else:
        if(is_even(c)):
            if verbose>2:
                print "even c: -a=",-a," g=",g
            #g2=gcd(4*c, lcm(8,N4))
            g2=gcd(4*c, 4*D)
            if verbose>2:
                print "(8c, lcm(8,4N))=",g2
            for j in range(1,10000):
                t= (a*b*c-d) +j*g2
                if(psilent>2):
                    print "t(",j,")=",t
                if(t > 0 and gcd(t,N4)==1):
                    break
        else:
            if(gcd(d,N4)==1):
                    t= -d  % (N4)
            else:
                for j in range(1,10000):
                    t= -d + j*g
                    g2=gcd(t,N4)
                    if( t > 0  and g2==1 ):
                        #print "t>0! and gcd(t,N4)=",g2
                        break
    if(gcd(t,N4)<>1):
        raise ValueError, "Could not find t for D=%s, a,b,c,d=%s,  tried t=%s" %(D,[a,b,c,d],t)
    return t


def holom_poincare_c_vec(l,kappa, N,NN, maxit,b_m,b_p, prec,tol=1E-40,verbose=0):
    r"""
    Coefficient (Delta') of Poincare series (Delta)
    (corr. to the dual rep.)
    l = [ [(D_1,r_1),(D'_1,r'_1)],[(D_2,r_2),(D'_2,r'_2)],...] = list of coefficients to compute
    """
    poincareseries_init(prec)
    if verbose>0:
        print "tol=",tol
        print "l=",l
        print "range(j)=",len(l)
    N4=4*N
    N2=2*N
    ZN4=IntegerModRing(N4)
    ZN2=IntegerModRing(2*N)
    rs=dict()
    Deltas=dict()
    rps=dict()
    Deltaps=dict()
    ms=dict()
    ns=dict()
    mms=dict()
    nns=dict()
    c1s=dict()
    cstops=dict()
    ars=dict()
    errs=dict()
    summas=dict()
    results=dict()
    results['ok']=None
    results['data']=None
    for j in range(len(l)):
        Deltas[j] =l[j][0]
        rs[j]=l[j][1]
        Deltaps[j]=l[j][2]
        rps[j]=l[j][3]
        a=ZN4(Deltas[j])
        b=ZN4(Deltaps[j])
        if verbose>0:
            print "Delta=",Deltas[j]
            print "Delta'=",Deltaps[j]
            print "r=",rs[j]
            print "r''=",rps[j]
            # test that Delta, Delta' are valid Heegner discriminants
            print "D % 4N=",a,type(a)
            print "D' % 4N=",b,type(b)
            print "sqrta=",a.sqrt()
            print "sqrtb=",b.sqrt()
        if( ZN4(rs[j]**2) <> ZN4(Deltas[j]) ):
            raise ValueError, "D=%s <> r^2, r=%s " % (Deltas[j],rs[j])
        if( ZN4(rps[j]**2) <> ZN4(Deltaps[j]) ):				       raise ValueError, "D'=%s <> r'^2, r=%s " % (Deltas[j],rs[j])
        rsq=(rs[j]**2 % N4)
        rpsq=(rps[j]**2 % N4)
        ms[j]=(rsq-Deltas[j])/N4
        ns[j]=(rpsq-Deltaps[j])/N4
        mms[j]=RF(-Deltas[j]/N4)
        nns[j]=RF(-Deltaps[j]/N4)
        if(mms<0):
            raise ValueError,"This is a non-holomorphic Poincare series! m=%s " % mms[j]
        c1s[j]=zero
        if( (Deltas[j] == Deltaps[j])):
            if( ZN2(rs[j]-rps[j])==0):
                c1s[j]=c1s[j]+one
            if( ZN2(rs[j]+rps[j])==0):
                c1s[j]=c1s[j]+one
        # endif:
        summas[j]=zero
        tmp=mms[j]*nns[j]
        ars[j]=two*twopi*tmp.sqrt()
        #    print("ar=",ar)
    # endfor
    terms=dict()
    averages=dict()
    breaks=dict()
    for j in   rs.keys():
        terms[j]=vector(CF,20)
        averages[j]=one
        breaks[j]=False
        cstops[j]=max(NN+1,22)
    kappa_minus_one=RF(kappa)-one
    WK=dict()
    av_kl=zero
    n_kl=0
    cmax=Integer(max(NN+1,22))
    #print "typeNN+1=",type(NN+1)
    #print "typecmax=",type(cmax)
    for c in range(1,cmax):
        cr=RF(c)
        WK=weil_kloosterman_sum_all(c,rs,mms,rps,nns,N,kappa,b_m,b_p,prec,breaks)
        if verbose>3:
            WKold=weil_kloosterman_sum_all_old(c,rs,mms,rps,nns,N,kappa,b_m,b_p,prec)
            for j in rs.keys():
                tmpd=abs(WK[j]-WKold[j])        
                if(tmpd>0.1):
                    print "diffWK[",j,"]=",tmpd

                    print "args=",c,rs,mms,rps,nns,N,kappa,b_m,b_p,prec
                    print "WK_new(",c,")=",WK
                    print "WK_old(",c,")=",WKold
                    raise Exception, "Stop!"
#       # return
        for j in rs.keys():
            av_kl=av_kl+abs(WK[j])
            n_kl=n_kl+1
        for j in rs.keys():
            if(breaks[j]):
                continue
            errs[j]=100
            tmp1=WK[j]
            if(abs(tmp1)<epsilon_0):
                if(c <=20):
                    terms[j][c-1]=tmp1
                do_jbes=False
            else:
                do_jbes=True
            # endif
            if(do_jbes):
                x=ars[j]/cr
                tmpj=RF(bessel_J(kappa_minus_one,x,prec=prec))
                #print "bessel=",tmpj
                tmp=tmpj*tmp1.real()
            else:
                tmp=tmp1.real()
            summas[j]=summas[j]+tmp
            # I want to break out of the loop when the remainder is sufficiently small.
            # Since the terms may be sporadically small I measure an average over the past 10 coefficients.
            if verbose>2:
                print "summas[",j,"]=",summas[j]
                print "c=",c
                print "tmp=",tmp
            if(c<=20):
                terms[j][c-1]=summas[j]
            else:
                nnz=0
                for k in range(19):
                    if(abs(terms[j][k+1])>zero):
                        nnz=nnz+1
                        terms[j][k]=terms[j][k+1]
                if(abs(tmp)>zero):
                    terms[j][19]=summas[j]
                if(abs(terms[j][19])>0):
                    nnz=nnz+1
                averages[j]=sum(terms[j])/RF(nnz)
#                for k in range(20):
#                    print "term[",j,k,"]=",terms[j][k]
#                print "nonzero=",nnz
                #                print "average[",j,c,"]=",averages[j]
                #                if(averages[j]<tol and c>max(50,NN)):
                tmp=abs(averages[j])+abs(summas[j])
                rel_err=abs(summas[j]-averages[j])/tmp
                errs[j]=min(rel_err,tmp)
                if verbose>1:
                    print "summas("+str(c)+")="+str(summas[j])
                    print "averagess("+str(c)+")="+str(averages[j])
                    print "rel_err("+str(c)+")="+str(rel_err)
                if((c>4*N) and ((rel_err<tol) or  tmp < 10*epsilon_0)):
                    breaks[j]=True
                    cstops[j]=c
                    if verbose>0:
                        print "Break for j=",j
                        print "rel_err("+str(c)+")="+str(rel_err)
                        print "tmp("+str(c)+")="+str(tmp)
                        print "summas("+str(c)+")="+str(summas[j])
                        print "averagess("+str(c)+")="+str(averages[j])
                        if(tmp<10*epsilon_0):
                            print "WK_new(",c,")[j]=",WK[j]
                    continue 
                # endif
            # endif
        # endfor j
        # Test if all components are ok up to desired precision
        break_all=True
        for j in rs.keys():
            if(breaks[j]==False):
                break_all=False
                break
        if(break_all):
            break
    # endfor c
    tmp=(kappa_minus_one)/two
    fourpi=two*twopi
    cp=dict()
    #print "c=",c
    true_err=dict()
    for j in rs.keys():
        fak=fourpi*(nns[j]/mms[j])**tmp
        cp[j]=CF(c1s[j]+fak*summas[j])
        ## Check out the proven bound also
        # print "j=",j
        # print "cstops=",cstops
        l=[Deltas[j],rs[j],Deltaps[j],rps[j]]
        tmperr=proven_bound(kappa,prec,N,l,cstops[j])
        # print "proven error=",tmperr.n(20)
        true_err[j]=tmperr
        if(abs(tmperr)<errs[j]):
            errs[j]=abs(tmperr)
        if(errs[j]<tol):
            breaks[j]=True
        if verbose>0:
            print "c1=",c1s[j]
            print "summas=",summas[j]
            print "mm=",mms[j]
            print "nn=",nns[j]
            print "fak=",fourpi,"*",nns[j],"/",mms[j],"**",tmp,'=',fak
        #print "cp[",j,"]=",cp[j]
    if verbose>0:
        print "Average size of Kloosterman sums=",av_kl/RF(n_kl)
    results['data']=cp
    results['ok']=breaks
    results['errs']=true_err
    return results
#enddef



#def holom_poincare_c_vec_old(l,kappa, N,NN, maxit,b_m,b_p, prec,tol=1E-40):
#     r"""
#     Coefficient (rp,Deltap) of Poincare series (r,Delta)
#     (corr. to the dual rep.)
#     l = list of coefficients to compute
#         the format is
#         [Delta,r,Delta',r']
#     """
#     poincareseries_init(prec)
#     if(psilent>0):
#         print "tol=",tol
#         print "l=",l
#         print "range(j)=",len(l)
#     N4=4*N
#     N2=2*N
#     ZN4=IntegerModRing(N4)
#     ZN2=IntegerModRing(2*N)
#     rs=dict()
#     Deltas=dict()
#     rps=dict()
#     Deltaps=dict()
#     ms=dict()
#     ns=dict()
#     mms=dict()
#     nns=dict()
#     c1s=dict()
#     ars=dict()
#     errs=dict()
#     summas=dict()
#     results=dict()
#     results['ok']=None
#     results['data']=None
#     results['keys']=l
#     for j in range(len(l)):
#         #        print "l=",l[j]
#         Deltas[j] =l[j][0]
#         Deltaps[j]=l[j][1]
#         # test that Delta, Delta' are valid Heegner discriminants
        
#         if(not (is_square(ZN4(Deltas[j]))) or not (is_square(ZN4(Deltaps[j])))):
#             print "is_square(-D)=",is_square(ZN4(Deltas[j]))
#             print "is_square(-D')=",is_square(ZN4(Deltaps[j]))
#             print "N=",N
#             print "ZN4=",ZN4
#             s="Input are not Discriminants satisfying a Heegner condition!  D=%s,  D'=%s" %(Deltas[j],Deltaps[j])
#             raise ValueError, s

#         rs[j] =Integer(sqrt(ZN4(Deltas[j])))
#         rps[j]=Integer(sqrt(ZN4(Deltaps[j])))
#         print "r,D=",rs[j],Deltas[j]
#         print "r',D'=",rps[j],Deltaps[j]
#         rsq=(rs[j]**2 % N4)
#         rpsq=(rps[j]**2 % N4)
#         #    !! Must check consistency of Delta,r in input

#         ms[j]=(rsq-Deltas[j])/N4
#         ns[j]=(rpsq-Deltaps[j])/N4
#         #    print "RF2=",RF
#         mms[j]=RF(-Deltas[j]/N4)
#         nns[j]=RF(-Deltaps[j]/N4)
#         if((Deltas[j]<>Deltaps[j])  or (rs[j]<>rps[j])):
#             c1s[j]=zero
#         elif((ZN2(rs[j]-rps[j])==0) and (ZN2(rs[j]+rps[j])==0)):
#             c1s[j]=two
#         else:
#             c1s[j]=one
#         # endif:
#         summas[j]=zero
#         tmp=mms[j]*nns[j]
#         ars[j]=two*twopi*tmp.sqrt()
#         #    print("ar=",ar)
#     # endfor
#     terms=dict()
#     averages=dict()
#     breaks=dict()
#     for j in   rs.keys():
#         terms[j]=vector(CF,20)
#         averages[j]=one
#         breaks[j]=False
#     kappa_minus_one=kappa-one
#     WK=dict()
#     av_kl=zero
#     n_kl=0
#     for c in range(1,max(NN+1,22)):
#         cr=RF(c)
#         WK=weil_kloosterman_sum_all_old(c,rs,mms,rps,nns,N,kappa,b_m,b_p,prec)
#         if(psilent>0):
#             print "WK_old(",c,")=",WK
#         for j in rs.keys():
#             av_kl=av_kl+abs(WK[j])
#             n_kl=n_kl+1
#         for j in rs.keys():
#             errs[j]=100
#             if(breaks[j]):
#                 continue
#             tmp1=WK[j]
#             if(abs(tmp1)<epsilon_0):
#                 if(c <=20):
#                     terms[j][c-1]=tmp1
#                 do_jbes=False
#             else:
#                 do_jbes=True
#             # endif
#             if(do_jbes):
#                 x=ars[j]/cr
#                 tmpj=RF(bessel_J(kappa_minus_one,x,prec=prec))
#                 tmp=tmpj*tmp1.real()
#             else:
#                 tmp=tmp1.real()
#             summas[j]=summas[j]+tmp
#             #print "tmp=",tmp
#             # I want to break out of the loop when the remainder is sufficiently small.
#             # Since the terms may be sporadically small I measure an average over the past 10 coefficients.
#             if(psilent>1):
#                 print "summas[",j,"]=",summas[j]
#                 print "c=",c
#                 print "tmp=",tmp
#             if(c<=20):
#                 terms[j][c-1]=summas[j]
#             else:
#                 nnz=0
#                 for k in range(19):
#                     if(abs(terms[j][k+1])>zero):
#                         nnz=nnz+1
#                         terms[j][k]=terms[j][k+1]
#                 if(abs(tmp)>zero):
#                     terms[j][19]=summas[j]
#                 if(abs(terms[j][19])>zero):
#                     nnz=nnz+1
#                 averages[j]=sum(terms[j])/RF(nnz)
# #                for k in range(20):
# #                    print "terms[",j,k,"]=",terms[j][k]
# #print "num<>0=",nnz
#                 #                print "average[",j,c,"]=",averages[j]
#                 #                if(averages[j]<tol and c>max(50,NN)):
#                 tmp=abs(averages[j])+abs(summas[j])
#                 rel_err=abs(summas[j]-averages[j])/tmp
#                 errs[j]=min(rel_err,tmp)
#                 if(psilent>0):
#                     print "summas("+str(c)+")="+str(summas[j])
#                     print "averagess("+str(c)+")="+str(averages[j])
#                     print "rel_err("+str(c)+")="+str(rel_err)
#                 if((c>4*N) and (rel_err<tol) or  tmp < 10*epsilon_0):
#                     print "Break for j=",j
#                     breaks[j]=True
#                     continue 
#                 # endif
#             # endif
#         # endfor j
#         # Test if all components are ok up to desired precision
#         break_all=True
#         for j in rs.keys():
#             if(breaks[j]==False):
#                 break_all=False
#                 break
#         if(break_all):
#             break
#     # endfor c
#     tmp=(kappa_minus_one)/two
#     fourpi=two*twopi
#     cp=dict()
# #    print "c=",c
#     for j in rs.keys():
#         fak=fourpi*(nns[j]/mms[j])**tmp
#         cp[j]=CF(c1s[j]+fak*summas[j])
#         if(psilent>0):
#             print "c1=",c1s[j]
#             print "summas=",summas[j]
#             print "mm=",mms[j]
#             print "nn=",nns[j]
#             print "fak=",fourpi,"*",nns[j],"/",mms[j],"**",tmp,'=',fak
#             # print "cp[",j,"]=",cp[j]
#     if(psilent>0):
#         print "Average size of Kloosterman sums=",av_kl/RF(n_kl)
#     results['data']=cp
#     results['ok']=breaks
#     results['errs']=errs
#     return results
# #enddef





cpdef get_trunc_bd(N,k,l,prec,tol=1E-40,format='Disc'):
    r""" Gives the truncation point for a specific error in the computation of the holomorphic Poincare series coefficients
    """
    poincareseries_init(prec)        
    half=RF(0.5)
    three=RF(3.0)
    kappa=RF(k)
    if(format=='Disc'):
        try:
            [Delta,r,Deltap,rp]=l
        except ValueError:
            [(Delta,r),(Deltap,rp)]=l
    else:
        [r,n,rp,np]=l
        Delta=(r**2 % (2*N))+4*N*n
        Deltap=(rp**2 % (2*N))+4*N*np
    if(Delta==0):
        raise ValueError," Delta=0! l=%s" %l

    tmp=two*twopi*(RF(4*N))**(one/kappa)/tol**(one/kappa)
    #print "tmp1=",tmp
    #print "Delta=",Delta
    #print "Delta'=",Deltap
    tmp=tmp*(abs(RF(Deltap))**(half+one/kappa))/(abs(RF(Delta))**(half+three/(two*kappa)))
    #print "tmp2=",tmp
    tmp=tmp*(one/kappa/gamma(kappa-one))**(one/kappa)
    #print "tmp3=",tmp
    #print "|tmp|=",abs(tmp)
    #print "type(tmp3)=",type(tmp)
    #print "C>",tmp
    nret=ceil(abs(tmp))
    #print "ceil(|tmp|)=",nret
    return nret
# enddef


cpdef proven_bound(kappa,prec,N,l,NN,format='Disc'):
    r""" Returns the error bound we can prove if we have  summed up to NN
    """
    if(format=='Disc'):
        [Delta,r,Deltap,rp]=l
    else:
        [r,n,rp,np]=l
        Delta=(r**2 % (2*N))+4*N*n
        Deltap=(rp**2 % (2*N))+4*N*np
    tmp1=RF(4*0.04)/RF(3)*twopi
    tmp2=gamma(kappa)*(kappa+two)
    tmp3=Deltap/Delta
#    tmp3=sqrt(tmp3)**(kappa-one)
    tmp4=CF.pi()*sqrt(Delta*Deltap)/RF(N)
    tmp4=tmp4**kappa
    tmp=RF(tmp1/tmp2*tmp4)*RF(sqrt(tmp3)/NN)**(kappa-one)
    #print "tmp=",tmp
    return tmp

cpdef test_arg(l,M):
    r""" Test if the elements of the list l are in the set M or not
    """
    for x in l:
        if x not in M:
            raise TypeError, "%s not a member of %s" %(x,M)
