#-*- coding: utf-8 -*-
r"""
Programs to compute holomorphic vector-valued Poincaré series.
In particular to find a basis and compute aq gram matrix consisting of the Fouriercoefficients p_{D_i,r_i}(D_j,r_j)


AUTHOR:
 - Fredrik Strömberg (May 2010)


EXAMPLES::


    sage: time A=make_gram_matrices(1,5,19.5,100,1E-20,force_prec=True)
    CPU times: user 13.21 s, sys: 0.07 s, total: 13.28 s
    Wall time: 13.31 s
    sage: A[1]['data']
    {0: 1.9997172464628335446644215348, 1: 0.038518031648348148837786590124, 2: 2.0220235164389306047369407568}
    sage: A[1]['maxerr']
    2.5825453423927096265697084347e-24
    sage: A[1]['indices']
    {0: [-3, 1], 1: [-4, 0]}


"""
from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

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

from builtins import str
from builtins import range
import re,os
p= re.compile('MatrixSpace')
from .poincare_series_alg_vv import *
from sage.all import IntegerModRing,RR,RealField,matrix,norm,load,save,Matrix,pi,is_square,set_verbose,ZZ

import mpmath
# Unless the sage-addd-ons are installed 
try:
    d = dimension_jac_cusp_forms(2,1,-1)
    #load 'sage-add-ons/nils/jacobiforms/dimension_jac_forms.sage'
    #load 'sage-add-ons/nils/jacobiforms/sz_trace_formula.sage'
except:
    ## If this does not exist.
    def dimension_jac_cusp_forms(a,b,c):
        raise NotImplementedError("You need the 'dimension_jac_form' package! Please call this routine with a set dimension instead!")
silent=0


def get_discs(N,dmax,sgn):
    r"""
    Computes a list of pairs (D,r) where D is a discriminant with D=sgn*r^2 mod 4N
    and |D|<=dmax
    """
    ZR=IntegerModRing(4*N)
    lp=list()
    ln=list()
    for j in range(1,dmax):
        D=j
        for r in range(1,N):
            if(ZR(D-r*r)==0):
                #print D,r
                lp.append((D,r))
    for j in range(1,dmax):
        D=-j
        for r in range(1,N):
            if(ZR(D-r*r)==0):
                #print D,r
                ln.append((D,r))
    return lp,ln
        
def rn_from_Dn(N,sgn,l):
    r"""  compute (n,r) from (D,r) with D=4Nn+sgn*(r^2 mod 4N)
    """
    [D,r]=l
    N4=4*N
    ZR=IntegerModRing(4*N)
    x= ZR(D -sgn*r*r)
    if x != 0:
        raise ValueError(" Need D=sgn*r^2 mod 4N got N={0}, r={1} D={2}".format(N,r,D))
    n=(D-sgn*r*r)/ZZ(4*N)
    return [n,r]


def rn_from_D(N,sgn,D):
    r"""  compute (n,r) from (D,r) with D=4Nn+sgn*(r^2 mod 4N)
    """
    N4=4*N
    ZR=IntegerModRing(4*N)
    DD=sgn*ZR(D)
    if not DD.is_square():
        raise ValueError(" Need D square mod 4N got N={0},  D={1}".format(N,D))
    for j in range(N):        
        x= ZR(D -sgn*j*j)
        if x == 0:
            r=j
            n=(D-sgn*r*r)/ZZ(4*N)
            break
    return [n,r]




## Function to compute a whole set of Fourier coefficients for Poincaré series at once
def ps_coefficients_holomorphic_vec(N,weight,l,tol=1E-40,prec=501,maxit=10000,force_prec=False):
    r""" Coefficients of vector-valued holomorphic Poincaré series
         for the dual Weil representation corr. to the lattice ZZ with quadratic form q(x)=Nx^2.
         (Note: dual => r(T)_{a,b}=delta_{a,b}*e(-q(a))
    INPUT: l=list of data for coefficients to compute
    format='Disc' if data supplied is a list of tuples
    [([D,r],[D',r']),...] where D,D' are fundamental discriminants
    and D==r^2 mod 4*N and D'==r'^2 mod 4N
    The coefficients we are computing are the (r',D')-th coefficient of the (r,D)-th Poincaré series.
    
    Conversion between the two formats are:
    -D/4N = n - r^2/4N (for the dual representation) and
    D/4N  = n + r^2/4N for the standard weil representation
    """
    NNmax=0
    NN=dict()
    b_m=0
    b_p=1
    if(silent>1):
        print("l={0}".format(l))
    try:
        l.keys()
        ll=l
    except AttributeError:
        ll=dict()
        for j in range(len(l)):
            ll[j]=(l[j][0][0],l[j][0][1],l[j][1][0],l[j][1][1])
    for j in ll.keys():
        NN[j]=get_trunc_bd(N,weight,ll[j],prec,tol=0.1*tol)

        if(NN[j]>NNmax):
            NNmax=NN[j]
    if(NNmax>1E10):
        raise ValueError("Need too many ({0}) terms! Probably to small weight!".format(NN))
    if(silent>1):
        print("NNmax={0}".format(NNmax))
    #for x in NN.keys():
    #    #print "NN(",x,")=",NN[x]
    #print "tol2=",tol
    #print "force_prec(os_coeff)=",force_prec
    res=holom_poincare_c_vec(ll,weight,N,NNmax,maxit,b_m,b_p,prec,tol,force_prec)    
    #print "Data is ok?:",res['ok']
    #print "Data=",res['data']
    return res

def gram_matrix(N,weight,prec=501,tol=1E-40,sv_min=1E-1,sv_max=1E15,bl=None,set_dim=None,force_prec=False):
    r""" Computes a matrix of p_{r,D}(r',D')
    for a basis of P_{r,D}, i.e. dim linearly independent P's
    INPUT: N      = Integer
           weight = Real
    OPTIONAL: 
           tol    = error bound for the Poincaré series
           sv_min = minimal allowed singular value when determining whether a given set is linarly independent or not.
           sv_max = maximally allowed singular value
           bl     = list of pairs (D_i,r_i) from which  we compute a matrix of coeffficients p_{D_i,r_i}(D_j,r_j)
        """
    # If we have supplied a list of D's and r's we make a gram matrix relative to these
    # otherwise we find a basis, i.e. linearly independent forms with correct dimension
    # find the dimension
    wt='%.4f'% weight
    if(N<10):
        stN="0"+str(N)
    else:
        stN=str(N)
    v=dict()
    filename_work="__N"+stN+"-"+wt+"--finding basis.txt"
    fp=open(filename_work,"write")
    fp.write("starting to find basis")
    fp.close()
    if(silent>0):
        print("Forcing precision:{0}".format(force_prec))
    set_verbose(0)
    if(bl!=None): 
        dim=len(bl)
        l=bl
    else:
        if(set_dim!=None and set_dim >0):
            dim=set_dim
        else:
            dim=dimension_jac_cusp_forms(int(weight+0.5),N,-1)
        l=list_of_basis(N,weight,prec,tol,sv_min,sv_max,set_dim=dim)
    j=0
    for [D,r] in l.values():
        for [Dp,rp] in l.values():
            # Recall that the gram matrix is symmetric. We need only compute the upper diagonal
            if(list(v.values()).count([Dp,rp,D,r])==0):
                v[j]=[D,r,Dp,rp]
                j=j+1
    # now v is a list we can get into computing coefficients
    # first we print the "gram data" (list of indices) to the file 
    s=str(N)+": (AI["+str(N)+"],["
    indices=dict()
    for j in range(len(l)):
        Delta=l[j][0]
        r=l[j][1]
        diff=(r*r-Delta) % (4*N)
        if diff != 0:
            raise ValueError("ERROR r^2={0} not congruent to Delta={1} mod {2}!".format(r*r, Delta, 4*N))
        s=s+"("+str(Delta)+","+str(r)+")"
        indices[j]=[Delta,r]
        if j<len(l)-1:
            s=s+","
        else:
            s=s+"]),"
    s=s+"\n"
    if silent>0:
        print(s+"\n")
    filename2="PS_Gramdata"+stN+"-"+wt+".txt"
    fp=open(filename2,"write")
    fp.write(s)
    fp.close()
    try:
        os.remove(filename_work)
    except os.error:
        print("Could not remove file:{0}".format(filename_work))
        pass
    filename_work="__N"+stN+"-"+wt+"--computing_gram_matrix.txt"
    fp=open(filename_work,"write")
    fp.write("")
    fp.close()
    #print "tol=",tol
    #set_verbose(2)
    #print "force_prec(gram_mat)=",force_prec
    res=ps_coefficients_holomorphic_vec(N,weight,v,tol,prec,force_prec=force_prec)
    set_verbose(0)

    res['indices']=indices
    maxerr=0.0
    for j in res['errs'].keys():
        tmperr=abs(res['errs'][j])
        #print "err(",j,")=",tmperr
        if(tmperr>maxerr):
            maxerr=tmperr
        # switch format for easier vewing
        res['errs'][j]=RR(tmperr)
    if silent>0:
        print("maxerr={0}".format(RR(maxerr)))
    res['maxerr']=maxerr
    wt_phalf='%.4f'% (weight+0.5)
    filename3="PS_Gramerr"+stN+"-"+wt+".txt"
    fp=open(filename3,"write")
    wt
    s="MAXERR["+wt_phalf+"]["+stN+"]="+str(RR(maxerr))
    fp.write(s)
    fp.close()
    if(res['ok']):
        Cps=res['data']
    else:
        print("Failed to compute Fourier coefficients!")
        return 0
    RF=RealField(prec)
    A=matrix(RF,dim)
    kappa=weight
    fourpi=RF(4.0)*pi.n(prec)
    one=RF(1.0)
    N4=RF(4*N)
    C=dict()
    if(silent>1):
        print("v={0}".format(v))
        print("dim={0}".format(dim))
    lastix=0
    # First set the upper right part of A
    for j in range(dim):
        ddim=dim-j
        if(silent>1):
            print("j={0} ddim={1} lastix={2]".format(j,ddim,lastix))
        for k in range(0,ddim):
            # need to scale with |D|^(k+0.5)
            if(silent>1):
                print("k={0}".format(k))
                print("lastix+k={0}".format(lastix+k))
            mm=RF(abs(v[lastix+k][0]))/N4
            tmp=RF(mm**(weight-one))
            if(silent>1):
                print("ddim+k={0}".format(ddim+k))
            A[j,j+k]=Cps[lastix+k]*tmp
            C[v[lastix+k][0],v[lastix+k][1]]=Cps[lastix+k]
        lastix=lastix+k+1
    # And add the lower triangular part to mak the matrix symmetric
    for j in range(dim):
        for k in range(0,j):
            A[j,k]=A[k,j]
    # And print the gram matrix
    res['matrix']=A
    dold=mpmath.mp.dps
    mpmath.mp.dps=int(prec/3.3)
    AInt=mpmath.matrix(int(A.nrows()),int(A.ncols()))
    AMp=mpmath.matrix(int(A.nrows()),int(A.ncols()))
    for ir in range(A.nrows()):
        for ik in range(A.ncols()):
            AInt[ir,ik]=mpmath.mpi(A[ir,ik]-tol,A[ir,ik]+tol)
            AMp[ir,ik]=mpmath.mpf(A[ir,ik])
    d=mpmath.det(AMp)
    if(silent>1):
        print("det(A-as-mpmath)={0}".format(d))
    di=mpmath.det(AInt)
    if(silent>1):
        print("det(A-as-interval)={0}".format(di))
    res['det']=(RF(di.a),RF(di.b))
    
    filename="PS_Gram"+stN+"-"+wt+".txt"
    if(silent>1):
        print("printing to file: {0}".format(filename))
    print_matrix_to_file(A,filename,'A['+str(N)+']')
    if(silent>1):
        print("A-A.transpose()={0}".format(norm(A-A.transpose())))
    B=A^-1
    #[d,B]=mat_inverse(A)
    if(silent>1):
        print("A={0}".format(A.n(100)))
        print("det(A)={0}".format(di))
        print("Done making inverse!")
    #res['det']=d
    res['inv']=B
    mpmath.mp.dps=dold
    filename="PS_Gram-inv"+stN+"-"+wt+".txt"        
    print_matrix_to_file(B,filename,' AI['+str(N)+']')
    # first make the filename
    s='%.1e'%tol
    filename3="PS_Coeffs"+stN+"-"+wt+"-"+s+".sobj"
    # If the file already exist we load it and append the new data
    if(silent>0):
        print("saving data to: {0}".format(filename3))
    try:
        f=open(filename3,"read")
    except IOError:
        if(silent>0):
            print("no file before!")
        # do nothing
    else:
        if silent>0:
            print("file: {0} exists!".format(filename3))
        f.close()
        Cold=load(filename3)
        for key in Cold.keys():
            #                print"key:",key
            if key not in C: # then we add it
                print("key:",key," does not exist in the new version!")
                C[key]=Cold[key]
                save(C,filename3)
    ## Save the whole thing
    filename="PS_all_gram"+stN+"-"+wt+".sobj"
    save(res,filename) 
    ## our work is completed and we can remove the file
    try:
        os.remove(filename_work)
    except os.error:
        print("Could not remove file: {0}".format(filename_work))
        pass
    return res

def list_of_basis(N,weight,prec=501,tol=1e-20,sv_min=1E-1,sv_max=1E15,set_dim=None):
    r""" Returns a list of pairs (r,D) forming a basis
    """
    # First we find the smallest Discriminant for each of the components
    if set_dim != None and set_dim >0:
        dim=set_dim
    else:
        dim=dimension_jac_cusp_forms(int(weight+0.5),N,-1)
    basislist=dict()
    num_gotten=0
    co_tmp=dict()
    num_gotten=0
    C0=1
    RF=RealField(prec)
    if(silent>1):
        print("N={0}".format(N))
        print("dim={0}".format(dim))
        print("sv_min={0}".format(sv_min))
        print("sv_max={0}".format(sv_max))
    Aold=Matrix(RF,1)
    tol0=1E-20  #tol
    # we start with the first discriminant, then the second etc.
    Z2N=IntegerModRing(2*N)
    ZZ4N=IntegerModRing(4*N)
    for Dp in [1..max(1000,100*dim)]:
        D=-Dp # we use the dual of the Weil representation
        D4N=ZZ4N(D)
        if(not(is_square(D4N))):
            continue
        for r in my_modsqrt(D4N,N):
            # I want to make sure that P_{(D,r)} is independent from the previously computed functions
            # The only sure way to do this is to compute all submatrices (to a much smaller precision than what we want at the end)
            # The candidate is [D,r] and we need to compute the vector of [D,r,D',r']
            # for all D',r' already in the list
            ltmp1=dict()
            ltmp2=dict()
            j=0
            for [Dp,rp] in basislist.values():
                ltmp1[j]=[D,r,Dp,rp]
                ltmp2[j]=[Dp,rp,D,r]
                j=j+1
            ltmp1[j]=[D,r,D,r]
            #print "Checking: D,r,D,r=",ltmp1
            ctmp1=ps_coefficients_holomorphic_vec(N,weight,ltmp1,tol0)
            # print "ctmp1=",ctmp1
            if(j >0):
                #print "Checking: D,r,Dp,rp=",ltmp2    # Data is ok?: {0: True} 
                ctmp2=ps_coefficients_holomorphic_vec(N,weight,ltmp2,tol0)
                # print "ctmp2=",ctmp2
            #print "num_gotten=",num_gotten
            A=matrix(RF,num_gotten+1)
            # The old matrixc with the elements that are already added to the basis
            # print "Aold=\n",A,"\n"
            # print "num_gotten=",num_gotten
            # print "Aold=\n",Aold,"\n"
            for k in range(Aold.nrows()):
                for l in range(Aold.ncols()):
                    A[k,l]=Aold[k,l]
                    # endfor
                    # print "A set by old=\n",A,"\n"
                    # Add the (D',r',D,r) for each D',r' in the list
            tmp=RF(1.0)
            for l in range(num_gotten):                
                # we do not use  the scaling factor when
                # determining linear independence
                # mm=RF(abs(ltmp2[l][0]))/N4
                # tmp=RF(mm**(weight-one))
                A[num_gotten,l]=ctmp2['data'][l]*tmp
                # Add the (D,r,D',r') for each D',r' in the list
                # print "ctmp1.keys()=",ctmp1.keys()
            for l in range(num_gotten+1):
                #mm=RF(abs(ltmp1[l][2]))/4N
                #tmp=RF(mm**(weight-one))
                # print "scaled with=",tmp.n(200)
                A[l,num_gotten]=ctmp1['data'][l]*tmp
            #[d,B]=mat_inverse(A) # d=det(A) 
            #if(silent>1):
            #d=det(A)
            #print "det A = ",d
            # Now we have to determine whether we have a linearly independent set or not
            dold=mpmath.mp.dps
            mpmath.mp.dps=int(prec/3.3)
            AInt=mpmath.matrix(int(A.nrows()),int(A.ncols()))
            AMp=mpmath.matrix(int(A.nrows()),int(A.ncols()))
            if(silent>0):
                print("tol0={0}".format(tol0))
            for ir in range(A.nrows()):
                for ik in range(A.ncols()):
                    AInt[ir,ik]=mpmath.mp.mpi(A[ir,ik]-tol0,A[ir,ik]+tol0)
                    AMp[ir,ik]=mpmath.mpf(A[ir,ik])

            d=mpmath.det(AMp)
            di=mpmath.mp.mpi(mpmath.mp.det(AInt))
            #for ir in range(A.nrows()):
            #    for ik in range(A.ncols()):
            #        #print "A.d=",AInt[ir,ik].delta
            if(silent>0):
                print("mpmath.mp.dps={0}".format(mpmath.mp.dps))
                print("det(A)={0}".format(d))
                print("det(A-as-interval)={0}".format(di))
                print("d.delta={0}".format(di.delta))
            #if(not mpmath.mpi(d) in di):
            #    raise ArithmeticError," Interval determinant not ok?"
            #ANP=A.numpy()
            #try: 
            #    u,s,vnp=svd(ANP) # s are the singular values
            #    sl=s.tolist()
            #    mins=min(sl)  # the smallest singular value
            #    maxs=max(sl)
            #    if(silent>1):
            #        print "singular values = ",s
            #except LinAlgError:
            #    if(silent>0):
            #        print "could not compute SVD!"
            #        print "using abs(det) instead"
            #   mins=abs(d)
            #    maxs=abs(d)
            #if((mins>sv_min and maxs< sv_max)): 
            zero=mpmath.mpi(0)
            if(zero not in di):
                if(silent>1):
                    print("Adding D,r={0}, {1}".format(D,r))
                basislist[num_gotten]=[D,r]
                num_gotten=num_gotten+1
                if(num_gotten>=dim):
                    return basislist
                else:
                    #print "setting Aold to A"
                    Aold=A
            else:
                if(silent>1):
                    print(" do not use D,r={0}, {1}".format(D,r))
            # endif
            mpmath.mp.dps=dold
    # endfor
    if(num_gotten < dim):
        raise ValueError(" did not find enough good elements for a basis list!")  



def make_gram_matrices(N1,N2,k,prec=501,tol=1E-40,sv_min=1E-1,sv_max=1E15,force_prec=False):
    A=dict()
    for N in range(N1,N2+1):
        A[N]=gram_matrix(N,k,prec,tol,sv_min,sv_max,force_prec=force_prec)
    return A


def print_matrix_to_file(A,filename,varname='A'):
    dimr=A.nrows()
    dimc=A.ncols()
    prec=A.base_ring().prec()
    fu=open(filename,"write")        
    #    fu.write(varname+"=matrix(RealField("+str(prec)+"),[")
    fu.write(varname+"=matrix(RF,[")
    for j in range(dimr):
        fu.write("[")
        for k in range(dimc):
            fu.write(str(A[j,k]))
            if(k<dimc-1):
                fu.write(",")
        if(j<dimr-1):
            fu.write("],")
        else:
            fu.write("]")
    fu.write("])\n")
    fu.close()


# the square root mod 2N is a set-valued function
def my_modsqrt(a,N):
    r""" Returns a list of $r\in{1,...,2N}$ with $r^2 \equiv a \mod 4 N$

    """
    N4=4*N
    amod4N=a % N4 
    #print "amod4N=",amod4N
    sq=[]
    for i in range(0,N+1):
        #print "i*i mod 4N=",((i*i) % N4)
        if( ((i*i) % N4) == amod4N):
            sq.append(i)
    return sq


def upper_submatrix(A,j0,k0):
    r"""
    Gives the upper submatrix
    B=A[j,k], j<=j0, k<=k0 

    """
    if(j0 > A.nrows()):
        j0=A.nrows()
    if( k0 > A.ncols()):
        k0=A.ncols()
    B=matrix(A.base_ring(),j0+1,k0+1)
    for j in range(j0+1):
        for k in range(k0+1):
            B[j,k]=A[j,k]
    return B

  
