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


import mpmath as mpmath
from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.structure.sage_object import SageObject,cPickle
from sage.functions.all import ln,sqrt,floor
from sage.rings.arith import divisors,gcd,inverse_mod
from sage.modular.dirichlet import DirichletGroup
from sage.rings.all import RR
from sage.modular.arithgroup.all import Gamma0
from sage.all import trivial_character,timeit,log,is_squarefree,prime_range,next_prime,deepcopy
from maass_forms_alg import *
from lpkbessel import *
from automorphic_forms import *
from eisenstein_series import Eisenstein_series_one_cusp
#from mysubgroup import is_Hecke_triangle_group
import matplotlib
import warnings

r"""

Maass waveforms for subgroups of the modular group and Eisenstein series for Hecke triangle groups.


AUTHORS:

 - Fredrik Strömberg (March 2010-)

EXAMPLES::


 ?


 TODO:
   - Nontrivial multiplier systems and weights
   - improve eigenvalue finding alorithms


"""


class MaassWaveForms (AutomorphicFormSpace):
    r"""
    Describes a space of Maass waveforms (cuspforms)
    """
    #def __init__(self,G,prec=500,ST=None,character=None,verbose=0,weight=0,**kwds):
    def __init__(self,G,weight=0,multiplier="",ch=0,sym_type=-1,atkin_lehner={},hecke=False,verbose=0,dprec=None,prec=53,**kwds):
        r"""
        Creates an ambient space of Maass waveforms

        INPUT:

            - `'G``    -- subgroup of the modular group
            - ``prec`` -- default working precision in bits 
            - ``dprec`` -- default working precision in digits 
            - ``Hecke`` -- if set to True we assume that we want Hecke eigenforms. In particular this implies that we have to work with complex numbers also for real characters.
            - ``sym_type` -- Integer. Even (0) or Odd (1)
            - ``atkin_lehner` -- Complex. Atkin-Lehner eigenvalues
        EXAMPLES::


            sage: S=MaassWaveForms(Gamma0(1)); S
            Space of Maass waveforms on the group G:
            Arithmetic Subgroup of PSL2(Z) with index 1. Given by: 
            perm(S)=()
            perm(ST)=()
            Constructed from G=Modular Group SL(2,Z)

        """
        if kwds.get('char_norm')=='Conrey' and ch>0:
            raise NotImplementedError
        #    for x in DirichletGroup_conrey():
        self._ch=ch
        self._hecke=hecke
        if dprec==None and prec==None:
            dprec=15; prec=53
        elif dprec==None:
            dprec=floor(RR(prec)/3.4)
        else:
            prec=ceil(3.4*dprec)+1
        AutomorphicFormSpace.__init__(self,G,weight=weight,multiplier=multiplier,character=ch,holomorphic=False,weak=False,cuspidal=True,unitary_action=1,dprec=dprec,verbose=verbose)
        self._eps = 2.0**(3.0-self._dprec)
        if sym_type == None:
            self._sym_type = -1 #self.set_norm()
        else:
            self._sym_type = sym_type
        self._Weyl_law_const=self._Weyl_law_consts()
        maass_forms=dict() # list of members 
        self._use_real=True
        if not self._multiplier.is_trivial():
            if not self._multiplier.is_real() or hecke:
                self._use_real=False
        if self._weight<>0:
            self._use_real=False
        if self._sym_type not in [0,1]:
            self._use_real=False
        self._symmetry=None
        self._even_odd_symmetries={}
        self._cusp_symmetries={}
        self._atkin_lehner_evs={}
        self._cusp_evs=[]
        #if atkin_lehner<>{}:
        self.set_cusp_evs(atkin_lehner)
        self._check_consistent_symmetrization()
        self._smallest_M0=0
        self._is_maass_waveform_space=True

    def _check_consistent_symmetrization(self):
        r"""
        Check that the symmetrization information supplied is consistent.
        """
        # check consistency between cusp eigenvalues and even/odd symmetry
        for i in range(self.group().ncusps()):            
            if self.cusp_evs()[i]<>0 and self.even_odd_symmetries()[i][0]<>1 and self.sym_type()<>-1:
                raise ValueError,"Got incompatible symmetry information!"
        if self.sym_type() in [0,1]:
            ## We can now use symmetries for most cusps (and leave the other with exponentials)
            if not self.group().is_congruence() and self.group().ncusps()>1:
                raise ValueError,"For non-cycloidal non-congruence subgroups we should not use even/odd symmetry!"
        
        ## Also check that the eigenvalues are compatible
        for c in range(self._group.ncusps()):
            o,d=self.cusp_symmetries().get(c,(-1,0))
            ev=self._atkin_lehner_evs.get(c,0)
            if o==0 and ev<>0:
                s = "The cusp nr. {0} does not appear to have an involution!".format(c)
                warnings.warn(s)
                self._cusp_evs[c]=0
            elif o>0:
                if ev<>0 and abs(ev**o-1)>self._eps:
                    s = "The cusp nr. {0} has involution of order {1} and the assumed eigenvalue {2} does not have this order!".format(c,o,ev)
                    warnings.warn(s)
                    self._cusp_evs[c]=0
        # check consistency between cusp eigenvalues and even/odd symmetry
        for i in range(self._group._ncusps):            
            if self.cusp_evs()[i]<>0 and self.even_odd_symmetries()[i][0]<>1 and self._sym_type<>-1:
                if self._verbose>0:
                    print "i=",i
                    print "self._cusp=",self._cusp_evs[i]
                    print "eo_sym=",self.even_odd_symmetries()[i]
                    print "sym_type=",self._sym_type
                raise ValueError,"Got incompatible symmetry information!"

            
        
    def weight(self):
        return self._weight

    def sym_type(self):
        return self._sym_type

    def atkin_lehner_eigenvalue(self,cusp):
        cc = self.group().cusp_representative(cusp)
        return self._atkin_lehner_evs.get(cc,0)
    
    def atkin_lehner_eigenvalues(self):
        return self._atkin_lehner_evs

    def cusp_evs(self):
        if not self._cusp_evs:
            self._cusp_evs=[1]
            for i in range(self._group._ncusps-1):
                self._cusp_evs.append(0)
        return self._cusp_evs

        
    def set_sym_type(self,s):
        if s not in [0,1,-1]:
            raise ValueError,"{0} is not a valid symmetry type!".format(s)
        # check consistency between cusp eigenvalues and even/odd symmetry
        for i in range(self._group._ncusps):            
            if self.cusp_evs()[i]<>0 and self.even_odd_symmetries()[i][0]<>1 and self._sym_type<>-1:
                raise ValueError,"Got incompatible symmetry information!"        
        self._sym_type=s
    
    def set_cusp_evs(self,evs={}):
        r"""
        Set Atkin-Lehner eigenvalues for self. If only on is supplied we assume it to be the Fricke eigenvalue (i.e. Atkin-Lehner at 0)

        INPUT:
          - ''evs'' -- dict or integer. 
        """
        if isinstance(evs,(int,Integer)):
            evs={0:int(eigenvalues)}
        elif not isinstance(evs,(dict,list)):
            raise TypeError,"Could not get cusp eigenvalues from {0}!".format(evs)
        self._atkin_lehner_evs={}
        for c in self.group().cusps():
            if c<>Cusp(1,0):
                self._atkin_lehner_evs[c]=0
            else:
                self._atkin_lehner_evs={Cusp(1,0):1}  # cusp at infinity                
        ### We want to use cusp representatives 
        if isinstance(evs,dict):
            for c in evs.keys():
                ck = self.group().cusp_representative(c)
                self._atkin_lehner_evs[ck]=evs[c]
            self._cusp_evs = []
            for c in self.group().cusps():
                self._cusp_evs.append(self._atkin_lehner_evs[c])
        else:
            self._cusp_evs = evs
            i=0
            assert len(evs)==self.group().ncusps()
            for c in self.group().cusps():
                self._atkin_lehner_evs[c]=evs[i]
                i+=1
        self._check_consistent_symmetrization()

    #return self._cusp_evs

    def even_odd_symmetries(self):
        r"""
        Check even/odd symmetries and behaviour with respect to the character
        """
        if self._even_odd_symmetries<>{}:
            return self._even_odd_symmetries
        res={}
        for j in range(self._group.ncusps()):
            if self._ch==0:
                res[j]=self._group.is_symmetrizable_even_odd(j),1
            else:
                if self._group.is_symmetrizable_even_odd(j)==0:
                    res[j]=0,0
                    continue
                a,b,c,d=self._group._cusp_data[j]['normalizer']
                res[j]=1,self._character(a*d+b*c)
        return res


    def cusp_symmetries(self):
        r"""
        Check cusp symmetries (involutions) and their behaviour with respect to the character
        """
        if self._cusp_symmetries<>{}:
            return self._cusp_symmetries
        res={}
        for j in range(self._group.ncusps()):
            o,d = self._group.cusp_normalizer_is_normalizer(j)
            if o==0:
                res[j]=0,0
            elif self._ch==0:
                res[j]=o,1
            else:
                x,y,z,w=self._group._cusp_data[j]['normalizer']
                l = self._group._cusp_data[j]['width']
                N = self.level()
                q = self._character.modulus()
                #if self._verbose>0:
                #    #print "lz=",l*z
                if (l*z % N) <> 0:
                    res[j]=0,0
                else:
                    ## Check that the character is uniquely determined on sigma_j gamma sigma_j^-1
                    ## for all gamma = (a b ; c d) in Gamma_0(N)
                    vals=[]
                    for a in range(N):
                        if gcd(a,N)>1:
                            continue
                        d = inverse_mod(a,N)
                        dp=x*w*d-y*z*a
                        xi = self._character(dp)
                        xj = self._character(d)
                        xij=xi/xj
                        if xij not in vals:
                            vals.append(xij)
                    # Note that if vals<>[1] then this map is not really an involution
                    # since it sends a character to another character
                    # but 
                    if len(vals)==1:
                        res[j]=o,vals[0]
                    else:
                        res[j]=0,0
        self._cusp_symmetries=res
        return res
    
    def __repr__(self):
        r"""
        Return the string representation of self.

        EXAMPLES::


            sage: M=MaassWaveForms(MySubgroup(Gamma0(1)));M
            Space of Maass waveforms on the group G:
            Arithmetic Subgroup of PSL2(Z) with index 1. Given by:
                perm(S)=()
                perm(ST)=()
            Constructed from G=Modular Group SL(2,Z)


        """
        s="Space of Maass waveforms "
        s+="of weight = "+str(self._weight)+" "
        if str(self._multiplier).find("theta_multiplier")>0:
            s+=" with theta multiplier "
        elif not self._multiplier.is_trivial():
            s+=" with multiplier:\n"+str(self._multiplier)
        else:
            s+=" with trivial multiplier "
        s+=" on "
        if self._group._is_Gamma0:
            
            s+='Gamma0({0})'.format(self.level())
        else:
            s+="the group G:\n"+str(self._group)
        return s


        return s

    def __reduce__(self):
        r""" Used for pickling.
        """
        return(MaassWaveForms,(self._group,self._weight,self._multiplier,self._character,self._sym_type,self._cusp_evs,self._hecke,self._verbose,self._dprec,self._prec))


    def __ne__(self,other):
        if self._verbose>1:
            print "in MaassWaveForms.__ne__"
        return not self.__eq__(other)

    def __eq__(self,other):
        if not hasattr(other,"_is_maass_waveform_space"):
            return False
        if self._verbose>1:
            print "in MaassWaveForms.__eq__"
        l0=self.__reduce__()
        l1=other.__reduce__()
        ## We allow to differ in precision and verbosity
        if l0[0]<>l0[0]:
            return False
        for j in range(0,5):
            #print "comparing A:",l0[1][j]
            #print "comparing B:",l1[1][j]
            if l0[1][j]<>l1[1][j]:
                #print "A<>B!"
                return False
        return True
    #return self.__reduce__() == other.__reduce__()
    
#G,weight=0,multiplier="",ch=0,dprec=None,prec=None,sym_type=None,verbose=0,hecke=True,**kwds):
    
    def __cmp__(self,other):
        r""" Compare self to other
        """
        if not isinstance(other,type(self)):
            return False
        if(self._group <> other._group or self.prec<>other.prec):
            return False
        else:
            return True

    def group(self):
        return self._group

    def is_congruence(self):
        return self._group.is_congruence()
    
    def level(self):
        if not self._group.is_congruence():
            raise ValueError,"Can only call level for a congruence subgroup!"
        return self._group.generalised_level()

    def get_element(self,R,Mset=None,Yset=None,dim=1,ndigs=12,set_c=[],**kwds):
        #if sym_type==None:
        #    sym_type=self._sym_type
        eps = 10**(1-ndigs)
        gr=kwds.get('gr',0)
        param=self.set_default_parameters(R,Mset,Yset,ndigs)
        Y=param['Y']
        Q=param['Q']
        M=param['M']
        oldf = kwds.get('oldforms')
        norm = kwds.get('norm')
        try: 
            if RR(R).is_infinity() or RR(R).is_NaN() or R<0.0:
                raise Exception
        except:
            raise ValueError,"R must be a (finite) non-negative real! Got R:{0}".format(R)
        if dim>1 and self.weight()==0 and oldf==None and norm==None:
            return  self.get_Hecke_basis(R,None,M,Y,dim,ndigs,set_c)
            ## We assume we have a scalar-valued Maass form for the moment
        else:
            if norm <> None:
                NN = norm
            else:
                NN = self.set_norm(dim,set_c=set_c)
            M0=0; Y0 = float(0.0)
            #if Mset<>None: M0 = int(Mset)
            #if Yset<>None: Y0 = float(Yset)
            if self._verbose>0:
                print "Y=",Y
                print "M=",M
                print "NN=",NN
            #if Y0==0.0 and M0==0:
            #    Y0,M0=get_Y_and_M_dp(self,R,eps)
            #if Y0==0.0 and M0<>0:
            #    Y0=get_Y_for_M_dp(S,R,M0,eps)
            #if Y0<>0 and M0==0:
            #    M0=get_M_for_maass_dp(R,Y0,eps)                
            #if self.weight()==0 and self.group().is_congruence():
            #    C = get_coeff_fast_cplx_dp_sym(self,R,Y,M,0,NN,gr=gr)
            #else:
            #    if self._verbose>0:
            #        print "Using routine without symmetry!"
            #        print "R,Y0,M=",R,Y0,M
            C = get_coeff_fast_cplx_dp(self,R,Y,M,0,NN,gr=gr)
            if gr<>0:
                return C
            if self._verbose>0:
                print "C.keys()==",C.keys()
                print "NN=",NN
            ## Make sure that the coefficient have Sage types
            CF = ComplexField(self._prec)
            for i in C.keys(): ## Functions in the space of dim>=1
                if self.multiplier().rank()>1:
                    for j in C[i].keys(): # Components
                        for r in C[i][j].keys(): # Cusps 
                            for n in C[i][j][r].keys(): # Coefficient c(n)
                                c = C[i][j][r][n]
                                C[i][j][r][n] = CF(c)
                else:
                    for r in C[i].keys(): # Cusps 
                        for n in C[i][r].keys(): # Coefficient c(n)
                            c = C[i][r][n]
                            C[i][r][n] = CF(c)
#            if len(C.keys())>1:
#                res = []
#                for i in C.keys():
#                    res.append(Maasswaveform(self,R,C=C[i],compute=#False,Y=Y,norm=NN))
#                return res
#            else:
            return Maasswaveform(self,R,C=C,compute=False,Y=Y,norm=NN,dim=dim)

            #X=coefficients_for_Maass_waveforms(self,R,Y,M,Q,ndigs,cuspidal=True,sym_type=sym_type,dim=dim,set_c=set_c)
            #F._coeffs[0]=X[0]
        #else:
        #    raise ValueError,"Can not compute Maass forms of dimension {0}".format(dim)
        
            
    def get_Hecke_basis(self,R,p=None,Mset=None,Yset=None,dim=1,ndigs=12,set_c=[]):
        if dim==1:
            return self.get_element(R,Mset,Yset,dim,ndigs,set_c)
        #NN = self.set_norm(dim)
        #param=self.set_default_parameters(R,Mset,Yset,ndigs)
        #Y0=param['Y']; Q=param['Q']; M0=param['M']
        NN = self.set_norm(dim); M0=0; Y0 = float(0.0); Q=0
        eps =1e-12
        if Mset<>None: M0 = int(Mset); Q=M0+10
        if Yset<>None: Y0 = float(Yset)
        if Y0==0.0 and M0==0:
            Y0,M0=get_Y_and_M_dp(self,R,eps)
        if Y0==0.0 and M0<>0:
            Y0=get_Y_for_M_dp(self,R,M0,eps)
        if Y0<>0 and M0==0:
            M0=get_M_for_maass_dp(R,Y0,eps)                
    
        if self._verbose>0:
            print "Get Hecke basis with:{0},{1},{2},{3},{4}".format(R,Y0,M0,Q,dim)
        #if self.weight()==0:
        #    X = get_coeff_fast_cplx_dp_sym(self,R,Y0,M0,0,NN)
        #else:
        X = get_coeff_fast_cplx_dp(self,R,Y0,M0,0,NN)
            #X = get_coeff_fast_cplx_dp_sym(self,R,Y0,M0,Q,NN)
        if p==None:
            p = self.get_primitive_p()
        H = self.Hecke_eigenfunction_from_coeffs(X,p)
        res = []
        for i in H.keys(): #range(dim):
            #print "H[",i,"][0][-1]=",H[i][0][-1]
            #C={0:H[i]}
            F = Maasswaveform(self,R,C={0:H[i]},dim=1,compute=False,hecke_p=p)
            res.append(F)
        return res

    
    def get_element_in_range(self,R1,R2,sym_type=None,Mset=None,Yset=None,dim=1,ndigs=12,set_c=None,neps=10):
        r""" Finds element of the space self with R in the interval R1 and R2

        INPUT:

        - ``R1`` -- lower bound (real)
        - ``R1`` -- upper bound (real)
        
        """
        # Dummy element
        F=Maasswaveform(self,R2,sym_type=sym_type,dim=dim,compute=False)
        param=self.set_default_parameters(R2,Mset,Yset,ndigs)
        Y=param['Y']
        Q=param['Q']
        M=param['M']
        if self._verbose>0:
            print "Y=",Y
            print "M=",M
            print "Q=",Q
        l=self.split_interval(R1,R2)
        if self._verbose>1:
            print "Split into intervals:"
            for [r1,r2,y] in l:
                print "[",r1,",",r2,"]:",y
        Rl=list()
        for [r1,r2,y] in l:
            [R,er]=find_single_ev(self,r1,r2,Yset=y,neps=neps)
            Rl.append([R,er])
        if self._verbose>0:
            print "R=",R
            print "er=",er


    def _Weyl_law_consts(self):
        r"""
        Compute constants for the Weyl law on self.group

        OUTPUT:

        - tuple of real numbers

        EXAMPLES::


            sage: M=MaassWaweForms(MySubgroup(Gamma0(1))
            sage: M._Weyl_law_consts  
            (0, 2/pi, (log(pi) - log(2) + 2)/pi, 0, -2)
        """
        import mpmath
        pi=mpmath.fp.pi
        ix=Integer(self._group.index())
        nc=self._group.ncusps()
        if(self._group.is_congruence()):
            lvl=Integer(self.level())
        else:
            lvl=0
        n2=Integer(self._group.nu2())
        n3=Integer(self._group.nu3())
        if is_Hecke_triangle_group(self._group): 
            if self._group._is_Gamma0:
                q=3
            else:
                q=self._group._q
            c1=(q-2)/(4*q)
        else:
            c1=ix/Integer(12)
        c2=Integer(2)*nc/pi
        c3=nc*(Integer(2)-ln(Integer(2))+ln(pi))/pi
        if lvl<>0:
            A=1
            for q in divisors(lvl):
                num_prim_dc=0
                DG=DirichletGroup(q)
                for chi in DG.list():
                    if(chi.is_primitive()):
                        num_prim_dc=num_prim_dc+1
                for m in divisors(lvl):
                    if(lvl % (m*q) == 0   and m % q ==0 ): 
                        fak=(q*lvl)/gcd(m,lvl/m)
                        A=A*Integer(fak)**num_prim_dc        
            c4=-ln(A)/pi
        else:
            c4=Integer(0)
        # constant term
        c5=-ix/144+n2/8+n3*2/9-nc/4-1
        return (c1,c2,c3,c4,c5)

    def Weyl_law_N(self,T,T1=None):
        r"""
        The counting function for this space. N(T)=#{disc. ev.<=T}
        
        INPUT:
        
        -  ``T`` -- double


        EXAMPLES::

            sage: M=MaassWaveForms(MySubgroup(Gamma0(1))
            sage: M.Weyl_law_N(10)
            0.572841337202191
            
        """
        (c1,c2,c3,c4,c5)=self._Weyl_law_const
        cc1=RR(c1); cc2=RR(c2); cc3=RR(c3); cc4=RR(c4); cc5=RR(c5)
        #print "c1,c2,c3,c4,c5=",cc1,cc2,cc3,cc4,cc5
        t=sqrt(T*T+0.25)
        try: 
            lnt=ln(t)
        except TypeError:
            lnt=mpmath.ln(t)
        #print "t,ln(t)=",t,lnt
        NT=cc1*t*t-cc2*t*lnt+cc3*t+cc4*t+cc5
        if(T1<>None):
            t=sqrt(T1*T1+0.25)
            NT1=cc1*(T1*T1+0.25)-cc2*t*ln(t)+cc3*t+cc4*t+cc5
            return RR(abs(NT1-NT))
        else:
            return RR(NT)

    def next_eigenvalue(self,R):
        r"""
        An estimate of where the next eigenvlue will be, i.e. the smallest R1>R so that N(R1)-N(R)>=1

        INPUT:
        - ``R`` -- real > 0

        OUTPUT:
        - real > R

        EXAMPLES::

            sage: M.next_eigenvalue(10.0)
            12.2500000000000

        
        """
        #cdef nmax
        N=self.Weyl_law_N(R)
        try:
            for j in range(1,10000):
                R1=R+j*RR(j)/100.0
                N1=self.Weyl_law_N(R1)
                if(N1-N >= 1.0):
                    raise StopIteration()
        except StopIteration:
            return R1
        else:
            raise ArithmeticError,"Could not find next eigenvalue! in interval: [%s,%s]" %(R,R1)
        
    def Weyl_law_Np(self,T,T1=None):
        r"""
        The derviative of the counting function for this space. N(T)=#{disc. ev.<=T}
        
        INPUT:
        
        -  ``T`` -- double


        EXAMPLES::

            sage: M=MaassWaweForms(MySubgroup(Gamma0(1))
            sage: M.Weyl_law_Np(10)
        
        """
        (c1,c2,c3,c4,c5)=self._Weyl_law_const
        cc1=RR(c1); cc2=RR(c2); cc3=RR(c3); cc4=RR(c4); cc5=RR(c5)
        #print "c1,c2,c3,c4,c5=",c1,c2,c3,c4,c5
        NpT=2.0*cc1*T-cc2*(ln(T)+1.0)+cc3+cc4
        return RR(NpT)


    def set_default_parameters(self,R,Mset=0,Yset=0,ndigs=12):
        r"""
        Try to set default parameters for computing Maass forms.
        """
        res=dict()
        #R=self._R
        eps=RR(10)**RR(-ndigs)
        if Yset>0:
            Y=float(Yset)
            M0=get_M_for_maass_dp(float(R),float(Y),float(eps))
        else:
            M0 = self.smallest_M0()
            i=0; YY = -1.0
            while i<1000 and YY<0:                
                M0+=1
                YY = get_Y_for_M_dp(self,R,M0,eps)
                i+=1
            #YY = float(self._group.minimal_height()*0.95)
            #M0 = get_M_for_maass_dp(R,YY,eps)
            #[YY,M0]=find_Y_and_M(self._group,R,ndigs)
            if ndigs>=15:
                Y =  RealField(ndigs*3.4)(YY)
            else:
                Y=RR(YY)
        if Mset > 0:
            M=Mset
        else:
            M=M0
        Q=M+10 
        res['Q']=Q
        res['M']=M
        res['Y']=Y
        return res  


    def set_norm(self,k=1,cuspidal=True,sym_type=None,set_c=[],atkin_lehner={}):
        r""" Set normalization for computing maass forms.
        
        INPUT:
        
        - ``k`` -- dimension
        - ``cuspidal`` -- cuspidal maass waveforms (default=True)
        

        OUTPUT:
        
        - ``N`` -- normalization (dictionary)
        -- N['comp_dim'] -- dimension of space we are computing in 
        -- N['SetCs']``     -- which coefficients are set
        -- N['Vals']     -- values of set coefficients     
        
        EXAMPLES::
        
        sage: set_norm_maass(1)
        {'Vals': {0: {0: 0, 1: 1}}, 'comp_dim': 1, 'num_set': 2, 'SetCs': [0, 1]}
        sage: set_norm_maass(1,cuspidal=False)
        {'Vals': {0: {0: 1}}, 'comp_dim': 1, 'num_set': 1, 'SetCs': [0]}
         sage: set_norm_maass(2)
         {'Vals': {0: {0: 0, 1: 1, 2: 0}, 1: {0: 0, 1: 0, 2: 1}}, 'comp_dim': 2, 'num_set': 3, 'SetCs': [0, 1, 2]}
         
         """
        if k<1:
            raise ValueError,"Need to compute at least a one-dimensional space!"
        #if set_c<>[] and set_c<>None:
        #    raise NotImplementedError,"We haven't implemented set c yet!"
        C=dict()
        Vals=dict()
        #  set coeffs c(0),c(1),...,c(k-1) if not cuspidal
        #  set coeffs c(0)=0,c(1),...,c(k) if cuspidal 
        SetCs=dict()
        for j in range(k):
            SetCs[j]=[]
            Vals[j]={}
        ### If we have set some c's explicitly then we only set these (plus constant terms if cuspidal):
            
        if set_c<>[]:
            if len(set_c)<>k:
                raise ValueError,"Need to give a complete set of set coefficients!"
            for j in range(k):
                if cuspidal:
                    for c in range(self.group().ncusps()):
                        if self.alpha(c)[0]==0:
                            SetCs[j].append((c,0))
                            Vals[j][(c,0)]=0
                for r,n in set_c[j].keys():
                    if (r,n) not in SetCs[j]:
                        SetCs[j].append((r,n))
                    Vals[j][(r,n)]=set_c[j][(r,n)]

        else:
            if cuspidal:
                for j in range(k):
                    for l in range(0,k+1):
                        SetCs[j].append((0,l))
                #SetCs[j]=range(0,k+1)
                    for i in range(1,self.group().ncusps()):
                        if SetCs[j].count((i,0))==0 and self.alpha(i)[0]==0: 
                            SetCs[j].append((i,0))
            else:
                for j in range(k):
                    SetCs[j]=[]
                    for l in range(0,k):
                        SetCs[j].append((0,l))

            for j in range(k):
                for r,n in SetCs[j]:
                    Vals[j][(r,n)]=0
        ## Set all valued = 0 first
            for j in range(k):            
                #print "Set Vals cuspidal=",cuspidal
                #print "Vals=",Vals
                if cuspidal:
                    #if not Vals[j].has_key((0,j+1)): ## These are values to set if
                    Vals[j][(0,j+1)]=1
                else:
                    #if not Vals[j].has_key((0,j)):
                    Vals[j][(0,j)]=1
        # Make sure that all constant terms are set to 0    
        #if cuspidal:
        #    for i in range(1,self.group().ncusps()):
        #        Vals[j][(i,0)]=0
        C['cuspidal']=cuspidal
        C['comp_dim']=k
        C['SetCs']=SetCs
        C['Vals']=Vals
        if sym_type <> None:
            C['sym_type'] = sym_type
        if atkin_lehner and isinstance(atkin_lehner,dict):
            C['atkin_lehner']=atkin_lehner
        return C

    def set_norm2(self,k=1,cuspidal=True,sym_type=None,atkin_lehner={},use_c=[]):
        r""" Set normalization for computing maass forms.
        
        INPUT:
        
        - ``k`` -- dimension
        - ``cuspidal`` -- cuspidal maass waveforms (default=True)
        - ``use_c`` -- which coefficients to use 

        OUTPUT:
        
        - ``N`` -- normalization (dictionary)
        -- N['comp_dim'] -- dimension of space we are computing in 
        -- N['SetCs']``     -- which coefficients are set
        -- N['Vals']     -- values of set coefficients     
        
        EXAMPLES::
        
        sage: set_norm_maass(1)
        {'Vals': {0: {0: 0, 1: 1}}, 'comp_dim': 1, 'num_set': 2, 'SetCs': [0, 1]}
        sage: set_norm_maass(1,cuspidal=False)
        {'Vals': {0: {0: 1}}, 'comp_dim': 1, 'num_set': 1, 'SetCs': [0]}
         sage: set_norm_maass(2)
         {'Vals': {0: {0: 0, 1: 1, 2: 0}, 1: {0: 0, 1: 0, 2: 1}}, 'comp_dim': 2, 'num_set': 3, 'SetCs': [0, 1, 2]}
         
         """
        C=dict()
        Vals=dict()
        if len(use_c)==0:
            if cuspidal == 0:
                use_c = range(k)
            else:
                use_c = range(1,k+1)
        if len(use_c)<>k:
                raise ArithmeticError,"Need the same number of coefficients to use as the dimension! Got dim={0}, use_c={1}".format(k,use_c)

        #  set coeffs c(0),c(1),...,c(k-1) if not cuspidal
        #  set coeffs c(0)=0,c(1),...,c(k) if cuspidal 
        SetCs=dict()
        for j in range(k):
            SetCs[j]=[]
            for l in range(k):
                SetCs[j].append((0,use_c[l]))            
            if cuspidal:
                for i in range(0,self.group().ncusps()):
                    if SetCs[j].count((i,0))==0 and self.alpha(i)[0]==0:
                        SetCs[j].append((i,0))


        if cuspidal:
            C['cuspidal']=True
        else:
            C['cuspidal']=False
        for j in range(k):
            Vals[j]=dict()
            for r,n in SetCs[j]:
                Vals[j][(r,n)]=0
        ## Set all valued = 0 first
        for j in range(k):
            Vals[j][(0,use_c[j])]=1

        # Make sure that all constant terms are set to 0    
        #if cuspidal:
        #    for i in range(1,self.group().ncusps()):
        #        Vals[j][(i,0)]=0
        C['comp_dim']=k
        C['SetCs']=SetCs
        C['Vals']=Vals
        if sym_type <> None:
            C['sym_type'] = sym_type
        if atkin_lehner and isinstance(atkin_lehner,dict):
            C['atkin_lehner']=atkin_lehner
        return C
    





    #### Split an interv
    def split_interval(self,R1,R2):
        r"""
        Split an interval into pieces, each containing (on average) at most one
        eigenvalue as well as a 0<Y<Y0 s.t. K_IR(Y) has no zero here

        INPUT:
        
            - ''R1'' -- real
            - ''R2'' -- real

        OUPUT:

            - list of triplets (r1,r2,y) where [r1,r2] does not contain a zero of K_ir(y)

        EXAMPLES::


            sage: M._next_kbes
            sage: l=M.split_interval(9.0,11.0)
            sage: print l[0],'\n',l[1],'\n',l[2],'\n',l[3]
            (9.00000000000000, 9.9203192604549457, 0.86169527676551638)
            (9.9203192704549465, 10.135716354265259, 0.86083358148875089)
            (10.13571636426526, 10.903681677771321, 0.86083358148875089)
            (10.903681687771321, 11.0000000000000, 0.85997274790726208)
            
        """
        import mpmath

        # It is enough to work with double precision
        base=mpmath.fp
        pi=base.pi
        # First we find the next zero 
        # First split into intervals having at most one zero
        ivs=list()
        rnew=R1; rold=R1
        while (rnew < R2):
            rnew=min(R2,self.next_eigenvalue(rold))
            if( abs(rold-rnew)==0.0):
                if self._verbose>0:
                    print "ivs=",ivs
                exit
            iv=(rold,rnew)
            ivs.append(iv)
            rold=rnew

        # We now need to split these intervals into pieces with at most one zero of the K-Bessel function
        Y00=base.mpf(0.995)*base.sqrt(base.mpf(3))/base.mpf(2 *self.level())
        new_ivs=list()
        for (r1,r2) in ivs:
            if self._verbose>0:
                print "r1,r2=",r1,r2
            Y0=Y00; r11=r1
            i=0
            while(r11 < r2 and i<1000):
                t=self._next_kbessel_zero(r11,r2,Y0*pi);i=i+1
                if self._verbose>0:
                    print "t=",t

                    oiv=(r11,t,Y0); new_ivs.append(iv)
                # must find Y0 s.t. |besselk(it,Y0)| is large enough
                Y1=Y0
                #k=base.besselk(base.mpc(0,t),Y1).real*mpmath.exp(t*0.5*base.pi)
                k=besselk_dp(RR(t),Y1)*exp(t*0.5*RR.pi())
                j=0
                while(j<1000 and abs(k)<1e-3):
                    Y1=Y1*0.999;j=j+1
                    #k=base.besselk(base.mpc(0,t),Y1).real*mpmath.exp(t*0.5*base.pi)
                    k=besselk_dp(RR(t),Y1)*exp(t*0.5*RR.pi())
                Y0=Y1
                r11=t+1E-08
        return new_ivs


    def _next_kbessel_zero(self,r1,r2,y):
        r"""
        The first zero after r1 i the interval [r1,r2] of K_ir(y),K_ir(2y)
        
        INPUT:

            - ´´r1´´ -- double
            - ´´r2´´ -- double
            - ´´y´´  -- double
    
        OUTPUT:

            - ''double''

        EXAMPLES::

        
            sage: M=MaassWaveForms(MySubgroup(Gamma0(1)))
            sage: Y0=0.995*sqrt(3.0)/2.0
            sage: M._next_kbessel_zero(9.0,15.0,Y0)
            9.9203192604549439
            sage: M._next_kbessel_zero(9.921,15.0,Y0)
            10.139781183668587


        CAVEAT:
            The rootfinding algorithm is not very sophisticated and might miss roots
            
        """
        base=mpmath.fp
        h=(r2-r1)/500.0
        t1=-1.0; t2=-1.0
        r0=base.mpf(r1)
        kd0=my_kbes_diff_r(r0,y,base)
        #print "r0,y=",r0,y,kd0
        while(t1<r1 and r0<r2):

            # Let us first find a R-value for which the derivative changed sign
            kd=my_kbes_diff_r(r0,y,base)
            i=0
            while(kd*kd0>0 and i<500 and r0<r2):
                i=i+1
                r0=r0+h
                kd=my_kbes_diff_r(r0,y,base)
                #print "r0,kd=",r0,kd
                #print "kd*kd0=",kd*kd0
            #print "-r0,y,kd=",r0,y,kd
            #t1=base.findroot(lambda x :  base.besselk(base.mpc(0,x),base.mpf(y),verbose=True).real,r0)
            try:
                t1=base.findroot(lambda x :  my_kbes(x,y,base),r0)
            except ValueError:
                t1=base.findroot(lambda x :  my_kbes(x,y,mpmath.mp),r0)
            r0=r0+h
        if(r0>=r2 or t1>=r2):
            t1=r2
        r0=r1
        kd0=my_kbes_diff_r(r0,y,base)
        while(t2<r1 and r0<r2):
            kd=my_kbes_diff_r(r0,y,base)
            i=0
            while(kd*kd0>0 and i<500 and r0<r2):
                i=i+1
                r0=r0+h
                kd=my_kbes_diff_r(r0,y,base)
            try:
                t2=base.findroot(lambda x :  my_kbes(x,2*y,base),r0)
            except ValueError:
                t2=base.findroot(lambda x :  my_kbes(x,2*y,mpmath.mp),r0)
            #t2=base.findroot(lambda x :  base.besselk(base.mpc(0,x),base.mpf(2*y),verbose=True).real,r0)
            r0=r0+h
        if(r0>=r2 or t2>=r2):
            t2=r2
            #print "zero(besselk,y1,y2)(",r1,r2,")=",t1,t2
        t=min(min(max(r1,t1),max(r1,t2)),r2)
        return t






    def Hecke_matrix(self,F,p):
        r"""
        Make the matrix of T_p with respect to the basis F.
        Here F is assumed to be a 0-1 normalized basis, s.t.
        F[i]=[a[0],a[1],...,a[M0]]
        and
        F[0]=[0,a[1],0,...,0,a[d+1],...]
        F[1]=[0,0,a[2],0,...,0,a[d+1],...]
        ...
        F[d]=[0,0,0,...,0,a[d],a[d+1],...]
        """
        #if Integer(p).divides(self._level):
        #    raise NotImplementedError,"Have only implemented primitive Hecke operators. Got q={0}|{1}".format(p,self._group._level)
        dim =len(F)
        assert self == F[0]._space
        if p*dim>len(F[0]._coeffs[0][0]):
            raise ValueError,"Need smaller p or more coefficients!"
        x=self.multiplier().character()
        if isinstance(F[0]._coeffs[0][0][1],(complex,float)):
            prec=53
        elif hasattr(F[0]._coeffs[0][0][1],"parent"):
            prec=F[0]._coeffs[0][0][1].parent().prec()
        CF = MPComplexField(prec)
        MS=MatrixSpace(CF,dim,dim)
        Tp = Matrix_complex_dense(MS,0)
        for i in range(dim):
            for j in range(dim):
                c=F[i]._coeffs[0][0][p*(j+1)]
                tmp=CF(c)
                if (j+1) % p ==0:
                    c = F[i]._coeffs[0][0][ZZ((j+1)/p)]
                    tmp+=x(p)*CF(c)
                Tp[i,j]=tmp
        return Tp


    

    def Hecke_eigenfunction(self,F,p,fnr=-1,verbose=0):
        r"""
        Construct a Hecke eigenfunction of T_p from the basis vector F
        If coeffs_only=1 then we only return the coefficients of the first component and/or cusp,\
        otherwise we return a MaassWaveformElement_class
        If fnr<0 we return a vector of all eigenfunctions
        """
        assert isinstance(F,(list,dict))
            ## If F is an Hecke eigenform we return it, otherwise return None
            #return F  ### In a one-
        assert self==F[0]._space
        ## Test if we already have Hecke eigenfunctions.
        is_hecke=1
        for j in range(len(F)):
            err = F[0].test(method='Hecke',format='float')
            if err>1E-6:
                is_hecke=0
                break
        if is_hecke==1:
            if verbose>0:
                print "We already have Hecke eigenfunctions!"
            return F
        C={}
        for j in range(len(F)):
            C[j]=F[j]._coeffs[0]  ## We assume scalar-valued functions
        #M=F[0]._space
        
        Cnew=self.Hecke_eigenfunction_from_coeffs(C,p,cusps='all',
                                             fnr=fnr,verbose=verbose)
        #return Cnew
        res=[];
        for j in Cnew.keys():
            FF=copy(F[0])  ## To get the same basic properties
            FF._coeffs={0:Cnew[j]}
            # If F[0] had a symmetry type set we try to find one for this too...
            FF._sym_type=FF.find_sym_type()
            res.append(FF)
        return res
        # #    l=(Tp.transpose()).eigenvectors()[fnr]
        # #print "l=",l
        # ev=l[0]
        # try:
        #     v=l[1][0]
        # except KeyError:
        #     print "l=",l
        #     raise ArithmeticError,"Problem computing eigenvectors!"
        # #print "v=",v
        # C1=dict()
        # if len(v)==2:
        #     res = F[0]._lin_comb(F[1],v[0],v[1])
        # else:
        #     res=F[0]*v[0]
        #     for j in range(1,len(v)):
        #         res=res+F[j]*v[j]
        # res = res*(1/v[0])

        # return res
        #for j in range(F._coeffs[0][0]):


    def Hecke_matrix_from_coeffs(self,C,p):
        r"""
        Make the matrix of T_p with respect to the basis F.
        Here F is assumed to be a 0-1 normalized basis, s.t.
        F[i]=[a[0],a[1],...,a[M0]]
        and
        F[0]=[0,a[1],0,...,0,a[d+1],...]
        F[1]=[0,0,a[2],0,...,0,a[d+1],...]
        ...
        F[d]=[0,0,0,...,0,a[d],a[d+1],...]
        """
        dim = len(C)
        if p*dim>len(C[0][0]):
            raise ValueError,"Need smaller p or more coefficients!\n Got: p={0} dim={1}, len(C)={2}".format(p,dim,len(C[0][0]))
        #assert p>dim:
        x=self.multiplier().character()
        #if Integer(p).divides(self.level()):
        #    raise NotImplementedError,"Have only implemented primitive Hecke operators. Got q={0}|{1}".format(p,self._group._level)
        #print "C001=",C[0][0][1]
        if isinstance(C[0][0][1],(complex,float)):
            prec=53
            #    Tp = Matrix(CC,dim,dim)
        elif hasattr(C[0][0][1],"parent"):
            prec=C[0][0][1].parent().prec()
            #    Tp = Matrix(F[0]._coeffs[0][0][1].parent(),dim,dim)
        else:
            prec=53
        CF = MPComplexField(prec)
        MS=MatrixSpace(CF,dim,dim)
        Tp = Matrix_complex_dense(MS,0)
        if hasattr(x(p),"complex_embedding"):
            xp = x(p).complex_embedding(prec)
        else:
            xp = ComplexField(prec)(x(p))
        for i in range(dim):
            for j in range(dim):
                c=C[i][0][p*(j+1)]
                tmp=CF(c)
                if (j+1) % p ==0:
                    #print "adding c[",ZZ((j+1)/p)
                    c = C[i][0][ZZ((j+1)/p)]
                    tmp+=xp*CF(c)
                Tp[i,j]=tmp
        return Tp



    def Hecke_eigenfunction_from_coeffs(self,C,p,cusps='all',fnr=-1,verbose=0):
        r"""
        Construct a Hecke eigenfunction of T_p from the basis vector F
        If coeffs_only=1 then we only return the coefficients of the first component and/or cusp,\
        otherwise we return a MaassWaveformElement_class
        If cusps = i we only compute coefficients at cusp i (usualle used if we only want coeffs. at infinity)

        """
        assert isinstance(C,(list,dict))
        assert dict_depth(C)>=3
        Tp=self.Hecke_matrix_from_coeffs(C,p)
        #while Tp.det()<Tp.eps():
        #    if verbose>0:
        #        print "Matrix is near singular, change p!"
        #
        if self._character(p)==-1:
            sorting=-1
        else:
            sorting=1
        if verbose>1:
            print "Tp=",Tp
        try:
            l=(Tp.transpose()).eigenvectors(verbose=0,sorted=sorting)
        except ArithmeticError:
            raise ArithmeticError, "Could not compute eigenvectors of the Hecke matrix T_{0}:\n{1}".format(p,Tp)
        if self._verbose>1 or verbose>0:
            for ll in l:
                print ll
        if fnr<0 or fnr>len(l): # compute all eigenfunctions
            fstart=0; fstop=len(l)
        else:
            fstart=fnr; fstop=fnr+1
        res=dict(); jj=0
        for j in range(fstart,fstop):
            ev=l[j][0]
            if self._verbose>1:
                print "ev=",ev
            if not isinstance(l[j][1],list):
                print "Tp=",Tp
                print "l=",l
                raise ArithmeticError,"Problem computing eigenvectors!"
            #if len(l[j][1])>1:
            #    raise ArithmeticError,"Eigenspace seems to be more than one-dimensional! For p={0} in the space {1}. \n eigenvalue={2} and  vector={3}".format(p,self,ev,l[j][1])
            for v in l[j][1]:
                #v=l[j][1][0] # Eigenvector
                if self._verbose>1:
                    print "v=",v
                CF=ComplexField(v[0].prec())
                # Get a normalizing coefficient
                v_norm=0
                for i in range(len(v)):
                    if abs(v[i])>0:
                        v_norm=CF(v[i].real(),v[i].imag())
                        break
                if v_norm==0:
                    ## Presumably we had a too large-dimensional space 
                    continue
                # raise ArithmeticError,"Could not find non-zero Hecke eigenvector! \n Hecke matrix is:{0} \n".format(Tp.transpose())
                res[jj]=dict()
                for i in C[0].keys():
                    if cusps<>'all' and cusps<>i:
                        continue
                    if self.atkin_lehner_eigenvalue(i)<>0 and i>0:
                        res[jj][i]=self.atkin_lehner_eigenvalue(i)
                        continue
                    res[jj][i]=dict()
                    # print "C[0][",i,"]=",C[0][i]
                    for n in C[0][i].keys():
                        res[jj][i][n]=CF(0)
                        for k in range(len(v)): 
                            vj=CF(v[k].real(),v[k].imag())
                            res[jj][i][n]=res[jj][i][n]+vj*C[k][i][n]
                        res[jj][i][n]=res[jj][i][n]/v_norm
                jj+=1
        return res
    
    def max_assumed_dim(self):
        r"""
        If the character has components of order two the dimension
        of a generic type Maass form doubles.
        Similarly if N has square factors the dimension of generic spaces
        might increase by multiples of two.
        """
        d=1
        if not self._group._is_congruence:
            return d
        x = self._multiplier._character
        if not x.is_trivial():
            for xx in x.decomposition():
                if xx.order()<=2:
                    d=d*2
        for p,m in self.level().factor():
            if m>1:
                d=d*2
        return d
    
    def get_primitive_p(self,p0=0,notone=1):
        r"""
        Gives a prime p to use for Hecke operator on M
        p should be relative prime to the level of M._group
        and to the modulus of M._multiplier._character

        INPUT:

        - 'p0' -- return prime greater than p0
        - 'notone' -- if set to one we return a prime with chi(p)<>1
        """
        if not self._group._is_congruence:
            return next_prime(p0)
        m=self._multiplier
        x=m._character
        if hasattr(x,"modulus"):
            modulus=x.modulus()
        else:
            modulus=1
        prim_to=lcm(self.level(),modulus) 
        p00 = next_prime(p0)
        p01 = p00 + prim_to
        if notone:
            if self.level() % 9 ==0 :
                pq=3
                # celif self._group._level % 4 ==0 :
                #    pq=4
            else:
                pq=1
        
        for p in prime_range(p00,p01+1):
            if notone==1 and p%pq==1:
                continue
            if gcd(p,prim_to)==1:
                return p
        raise ArithmeticError," Could not find appropriate p rel. prime to {0}!".format(prim_to)


    def smallest_M0(self):
        r"""
        Smallest M0 which we can use if we want to test using Hecke relations.
        """
        if is_Hecke_triangle_group(self._group):
            if self._group.is_Gamma0():
                self._smallest_M0=int(12)
            else:
                self._smallest_M0=int(12*self._group._lambdaq)
        if self._smallest_M0>0:
            return self._smallest_M0
        a = self.get_primitive_p()
        b = self.get_primitive_p(a)
        c = a*b
        self._smallest_M0=c+3
        return self._smallest_M0


    def test_Hecke_relation(self,C={},a=0,b=0,signed=False):
        r"""Testing Hecke relations for the Fourier coefficients in C

        INPUT:
        -''C'' -- dictionary of complex (Fourier coefficients)
        -''a'' -- integer
        -''b'' -- integer
        -''signed'' -- Boolean. True if we want to return a positive or negative number

        OUTPUT:
        -''diff'' -- real : |C(a)C(b)-C(ab)| if (a,b)=1

        EXAMPLE::


        sage: S=MaassWaveForms(Gamma0(1))
        sage: R=mpmath.mpf(9.53369526135355755434423523592877032382125639510725198237579046413534)
        sage: Y=mpmath.mpf(0.85)
        sage: C=coefficients_for_Maass_waveforms(S,R,Y,10,20,12)
        sage: d=_test_Hecke_relations(C,2,3); mppr(d)
        '9.29e-8'
        sage: C=coefficients_for_Maass_waveforms(S,R,Y,30,50,20)
        sage: d=_test_Hecke_relations(C,2,3); mppr(d)
        '3.83e-43'


        """
        if a*b==0:
            a = self.get_primitive_p()
            b = self.get_primitive_p(a)
        c=gcd(Integer(a),Integer(b))
        if self._verbose>1:
            print "Test Hecke: a={0},b={1},gcd(a,b)={2}".format(a,b,c)
        if not C.has_key(0):
            return -1
        if not hasattr(C[0],"has_key"):
            C = {0:C}
        if C[0].has_key(a) and C[0].has_key(b) and C[0].has_key(a*b): 
            lhs=C[0][a]*C[0][b]
            rhs=0
            for d in divisors(c):
                if self._ch<>0:
                    x = self._character(d)
                else:
                    x = 1
                m = Integer(a*b/d/d)
                if self._verbose>1:
                    print "rhs+=c*C[0][{0}]={1}".format(m,x*C[0][m])
                rhs=rhs+x*C[0][m]
            

            if self._verbose>1:
                print "|rhs|=",abs(rhs)
                print "|lhs|=",abs(lhs)
                print "self._prec=",self._prec
                print "rhs/lhs-1=",rhs/lhs-1.0
            #if max(abs(rhs),abs(lhs))<max(1e-8,2.0**(-0.5*self._prec)):
            #    return -1
            ## We have to return true also for the zero function
            t = rhs-lhs
            if signed:
                return t #rhs/lhs-1
            else:
                return abs(t) #rhs/lhs-1)
        return -1


    def scattering_determinant(self,s):
        r"""
        Computes the scattering determinant, varphi(s), of self at s
        using the non-holomorphic Eisenstein series
        Only implemented for Hecke triangle groups at the moment.
        """        
        if not is_Hecke_triangle_group(self._group): 
            raise NotImplementedError,"Only implemented for Hecke triangle groups"
        E = EisensteinSeries(self,s,verbose=self._verbose)
        return E._coeffs[0][0]







## def my_kbes_diff_r(r,x,mp_ctx=None):
##     r"""

##     Approximation to the derivative with respect to R of the scaled K-Bessel function.

##     INPUT:

##         - ''r'' -- real
##         - ''x'' -- real
##         - ''ctx'' -- mpmath context (default mpmath.mp)

##     OUTPUT:

##         - real -- K_ir(x)*exp(pi*r/2)

##     EXAMPLES::


##         sage: my_kbes_diff_r(9.45,0.861695276766 ,mpmath.fp)
##         -0.31374673969963851
##         sage: my_kbes_diff_r(9.4,0.861695276766 ,mpmath.fp)
##         0.074219541623676832


##     """
##     import mpmath
##     if(mp_ctx==None):
##         mp_ctx=mpmath.mp
##     if(mp_ctx==mpmath.mp):
##         pi=mpmath.mp.pi()
##     else:
##         pi=mpmath.fp.pi
##     try:
##         k=mp_ctx.besselk(mp_ctx.mpc(0,r),ctx.mpf(x))
##         f=k*mp_ctx.exp(r*mp_ctx.mpf(0.5)*pi)
##     except OverflowError:
##         k=mp_ctx.besselk(mp_ctx.mpc(0,r),mp_ctx.mpf(x))
##         f=k*mp_ctx.exp(r*mp_ctx.mpf(0.5)*pi)
##     f1=f.real
##     try:
##         h=mp_ctx.mpf(1e-8)
##         k=mp_ctx.besselk(mp_ctx.mpc(0,r+h),mp_ctx.mpf(x))
##         f=k*mp_ctx.exp((r+h)*mp_ctx.mpf(0.5)*pi)
##     except OverflowError:
##         h=mp_ctx.mpf(1e-8)
##         k=mp_ctx.besselk(mp_ctx.mpc(0,r+h),mp_ctx.mpf(x))
##         f=k*mp_ctx.exp((r+h)*mp_ctx.mpf(0.5)*pi)
##     f2=f.real
##     diff=(f2-f1)/h
##     return diff

#class MaassWaveformElement (SageObject):


def Maasswaveform(space,eigenvalue,**kwds):
    r"""
    Return a Maass waveform as an element of type MaassWaveformElement_class
    """
    data = {'_space':space,'_R':eigenvalue}
    data['_sym_type'] = kwds.pop('sym_type',-1)
    if data['_sym_type']<>-1:
        data['_space'].set_sym_type(data['_sym_type'])
    else:
        data['_sym_type'] = data['_space'].sym_type()
    
    data['_verbose'] = kwds.pop('verbose',0)
    data['_prec'] = kwds.pop('prec',53)
    data['_coeffs']=kwds.pop('C',kwds.pop('coefficients',{}))
    cusp_evs = kwds.pop('cusp_evs',{})
    if cusp_evs<>{}:
        cusp_evs = data['_space'].set_cusp_evs(cusp_evs)
    data['_cusp_evs']=cusp_evs
    data['_atkin_lehner_evs'] = data['_space']._atkin_lehner_evs
    data['_set_c'] = kwds.pop('set_c',{})
    data['_dim'] = kwds.pop('dim',1)
    if data['_dim'] <> data['_space']._rdim:
        data['_space']._rdim = data['_dim']
    
    #    data['_nd'] = kwds.pop('',12)
    data['_compute'] = kwds.pop('compute',1)
#        ,G,R,C=None,nd=12,sym_type=None,cusp_evs={},verbose=None,prec=53,set_c=None,dim=1,compute=False,test=0,data={},**kwds):
    data['_errest'] = kwds.get('errest',0)
    data['_Y'] = kwds.get('Y',0)
    data['_M0'] = kwds.get('M0',0)
    data['_norm'] =  kwds.get('norm',{})
    data['_nd']=kwds.get('nd',12)
    ## If self is constructed as a Hecke eigenform with respect
    ## to T_p we don't want to use p for testing.    
    data['_from_hecke_p']=kwds.get('hecke_p',0)
    data['_test']=kwds.get('test',0)
    data['compute']=kwds.get('compute',0)
    data['_version']=1  ## Might be necessary in the future
    return MaassWaveformElement_class(data)




class MaassWaveformElement_class(AutomorphicFormElement): #(Parent):
    r"""
    An element of a space of Maass waveforms


    EXAMPLES::


        sage: G=MySubgroup(Gamma0(1))
        sage: R=mpmath.mpf(9.53369526135355755434423523592877032382125639510725198237579046413534)
        sage: F=Maasswaveform(G,R)    
        Maass waveform with parameter R=9.5336952613536
        in Space of Maass waveforms on the group G:
        Arithmetic Subgroup of PSL2(Z) with index 1. Given by:
            perm(S)=()
            perm(ST)=()
        Constructed from G=Modular Group SL(2,Z)
        sage: G=MySubgroup(Gamma0(4))
        sage: R=mpmath.mpf(3.70330780121891)
        sage: F=Maasswaveform(G,R);F
        Maass waveform with parameter R=3.70330780121891
        Member of the Space of Maass waveforms on the group G:
        Arithmetic Subgroup of PSL2(Z) with index 6. Given by: 
        perm(S)=(1,2)(3,4)(5,6)
        perm(ST)=(1,3,2)(4,5,6)
        Constructed from G=Congruence Subgroup Gamma0(4)
        sage: F.C(0,-1)
        mpc(real='-1.0000000000014575', imag='5.4476887980094281e-13')
        sage: F.C(0,15)-F.C(0,5)*F.C(0,3)
        mpc(real='-5.938532886679327e-8', imag='-1.0564743382278074e-8')
        sage: F.C(0,3)
        mpc(real='0.53844676975670527', imag='-2.5525466782958545e-13')
        sage: F.C(1,3)
        mpc(real='-0.53844676975666916', imag='2.4484251009604091e-13')
        sage: F.C(2,3)
        mpc(real='-0.53844676975695485', imag='3.3624257152434837e-13')
        
        

        
    
    """
#    def __init__(self,G,R,C=None,nd=12,sym_type=None,cusp_evs={},verbose=None,prec=53,set_c=None,dim=1,compute=False,test=0,data={},**kwds):
    def __init__(self,data={},**kwds):
        r"""
        Construct a Maass waveform on thegroup G with spectral parameter R and coefficients C
        
        INPUT:

        
        - ``G`` -- Group
        - ``R`` -- Spectral parameter
        - ``C`` -- Fourier coefficients (default None)
        - ``nd``-- Number of desired digits (default 15)


        EXAMPLES::
            sage: G=MySubgroup(Gamma0(1))
            sage: R=mpmath.mpf(9.53369526135355755434423523592877032382125639510725198237579046413534)
            sage: F=Maasswaveform(G,R)    
            Maass waveform with parameter R=9.5336952613536
            in Space of Maass waveforms on the group G:
            Arithmetic Subgroup of PSL2(Z) with index 1. Given by:
                perm(S)=()
                perm(ST)=()
            Constructed from G=Modular Group SL(2,Z)
            sage: G=MySubgroup(Gamma0(4))
            sage: R=mpmath.mpf(3.70330780121891)
            sage: F=Maasswaveform(G,R);F
            Maass waveform with parameter R=3.70330780121891
            Member of the Space of Maass waveforms on the group G:
            Arithmetic Subgroup of PSL2(Z) with index 6. Given by: 
            perm(S)=(1,2)(3,4)(5,6)
            perm(ST)=(1,3,2)(4,5,6)
            Constructed from G=Congruence Subgroup Gamma0(4)
            sage: F.C(0,-1)
            mpc(real='-1.0000000000014575', imag='5.4476887980094281e-13')
            sage: F.C(0,15)-F.C(0,5)*F.C(0,3)
            mpc(real='-5.938532886679327e-8', imag='-1.0564743382278074e-8')
            sage: F.C(0,3)
            mpc(real='0.53844676975670527', imag='-2.5525466782958545e-13')
            sage: F.C(1,3)
            mpc(real='-0.53844676975666916', imag='2.4484251009604091e-13')
            sage: F.C(2,3)
            mpc(real='-0.53844676975695485', imag='3.3624257152434837e-13')

        """
        #import mpmath
        self.__dict__.update(data)
        #print "data=",data
        #print "coeffs",self._coeffs
        AutomorphicFormElement.__init__(self,self._space,self._coeffs,prec=self._prec,principal_part={},verbose=self._verbose)
        #if nd>15:
        #    mpmath.mp.dps=nd
        #    self.mp_ctx=mpmath.mp
        #else:
        #    self.mp_ctx=mpmath.fp
        ## We use the Fourier coefficients to verify whether we really have an eigenvalue
        #self._coeffs={}
        #
        #
        if self._test==1 and self._errest<>0:
            self._errest = self.test()
            
        #self._M0 = max(self._coeffs.values().values().keys())
        dprec=2.**(-self._nd)
        self._M0=get_M_for_maass(self._R,
                                 self._space._group.minimal_height(),
                                 dprec) 
        if data.get('compute',False)==True and self._coeffs=={}:
            self._coeffs=self.get_coeffs()
        elif self._coeffs == {}:
            self._coeffs = {0: {}}
            for j in range(self._space._group.ncusps()):
                self._coeffs[0][j]={}
        
                    
    def _repr_(self):
        r""" Returns string representation of self.

        EXAMPLES::


            sage: R=mpmath.mpf(9.53369526135355755434423523592877032382125639510725198237579046413534)
            sage: F=Maasswaveform(Gamma0(1),R,nd=50);F
            Maass waveform with parameter R=9.5336952613535575543442352359287703238212563951073
            Member of the Space of Maass waveforms on the group G:
            Arithmetic Subgroup of PSL2(Z) with index 1.Given by
                perm(S)=()
                perm(ST)=()
            Constructed from G=Modular Group SL(2,Z)

        
        """
        s="Maass waveform with parameter R="+str(self._R)+". Symmetry: "
        if self._sym_type==1:
            sym = "odd "
        elif self._sym_type==0:
            sym = "even "
        else:
            sym = ""
        s+=sym
        if self._cusp_evs<>[]:
            s+="Atkin-Lehner eigenvalues at cusps:"+str(self._atkin_lehner_evs)
        s+="\nMember of the "+str(self._space)

        return s

   

    def __reduce__(self):
        r""" Used for pickling.
        """
        return(MaassWaveformElement_class,(self.__dict__,))
                                     #self._coeffs,self._nd,self._sym_type,self._cusp_evs,self._verbose,self._prec))

    def group(self):
        r"""
        Return self._group 
        """
        return self._space._group

    def level(self):
        return self._space.level()

    def generalised_level(self):
        return self._space._group.generalised_level()

    def eigenvalue(self):
        return self._R  #eigenvalue

    def find_sym_type(self,tol=1e-7):
        cnr=1
        st_old=-1
        c0 = self._coeffs[0][0][1]; cnr=1
        for k in range(self._M0):
            if abs(abs(self._coeffs[0][0][k])-1.0)<tol:
                c0 = self._coeffs[0][0][k]; cnr=k
                break
        if k>=self._M0:
            print "Could not find c[k] close to 1!"
        eosym = self._space.even_odd_symmetries()
        for j in eosym.keys():
            s,d = eosym[j]
            if s==0:
                continue
            c1=self._coeffs[0][j][cnr]
            if abs(c1-c0*d)<tol:
                st=0
            elif abs(c1+c0*d)<tol:
                st=1
            else:
                return -1
            if st<>st_old and st_old<>-1:
                return -1
            st_old=st
    
    def C(self,i,j=None,r=0):
        r"""
        Return the coefficient C(i,j) i.e. coefficient nr. j at cusp i (or if j=none C(1,i))


        EXAMPLES::

        
            sage: G=MySubgroup(Gamma0(1))
            sage: R=mpmath.mpf(9.53369526135355755434423523592877032382125639510725198237579046413534)
            sage: F=Maasswaveform(G,R)    
            sage: F.C(2)
            mpc(real='-1.068333551223568', imag='2.5371356217909904e-17')
            sage: F.C(3)
            mpc(real='-0.45619735450601293', imag='-7.4209294760716175e-16')
            sage: F.C(2)*F.C(3)-F.C(6)
            mpc(real='8.2470016583210667e-8', imag='1.6951583479643061e-9')


        """
        if(j==None):
            cusp=0
            j=i
        else:
            cusp=i
        if(not self._coeffs[r].has_key(cusp)):
                raise ValueError," Need a valid index of a cusp as first argument! I.e in %s" %self._coeffs.keys()
        if(not self._coeffs[r][cusp].has_key(j)):
            return None
        return self._coeffs[r][cusp][j]
        
    def test(self,method='Hecke',up_to_M0=0,format='digits',verbose=0):
        r""" Return the number of digits we believe are correct (at least) 
        INPUT:
        - method -- string: 'Hecke' or 'pcoeff' or 'TwoY'
        - format = 'digits' or 'float'
        EXAMPLES::

        
            sage: G=MySubgroup(Gamma0(1))
            sage: R=mpmath.mpf(9.53369526135355755434423523592877032382125639510725198237579046413534)
            sage: F=Maasswaveform(G,R)    
            sage: F.test()
            7



        """
        # If we have a Gamma_0(N) we can use Hecke operators
        verbose = max(verbose,self._space._verbose)
        if self.level()==1:
            method='Hecke'
        if method=='Hecke' and self._space._group.is_congruence():
            #a = self._space.get_primitive_p()
            #b = self._space.get_primitive_p(a)
            p = self._from_hecke_p
            a = self._space.get_primitive_p(p)
            b = self._space.get_primitive_p(a)
            if verbose>1:
                print "Check Hecke relations! a={0}, b={1}".format(a,b)

            if len(self._coeffs)>b+3:
                if verbose>=0:
                    print "Too few coefficients for test!"
                return 1
            er=self._space.test_Hecke_relation(self._coeffs[0],a=a,b=b)
            #print "er=",er
            if er == 0:
                return 0
            if er==-1:
                return self.test(method='pcoeff',up_to_M0=up_to_M0,format=format)
            if format=='float':
                return er
            d=floor(-log(er,10))
            if verbose>0:
                print "Hecke is ok up to ",d,"digits!"
            return d
        elif method=='pcoeff':
            #d1 = self.test(method='Hecke',format=format)
            if verbose>0:
                print "Testing prime coefficients!"
            N = self.level()
            d1 = 1
            x = self._space._character
            for p in prime_range(N):
                if x.is_trivial():
                    if valuation(N,p)==1:
                        d2 = abs(abs(self._coeffs[0][0][p])**2-RR(1)/RR(p))
                        if verbose>1:
                            print "Checking c({0})".format(p)
                        if d2<d1: d1=d2
                    elif valuation(N,p)==2:
                        d2 = abs(self._coeffs[0][0][p])
                        if d2<d1: d1=d2
                else:
                    p = x.conductor()
                    d2 = abs(abs(self._coeffs[0][0][p])-1)
                    if d2<d1: d1=d2
            return d1
        else:
            # Test two Y's. Since we don't know which Y we had to start with we use two new ones and compare against the coefficients we already have
            nd=self._nd #+5
            Ctmp = deepcopy(self._coeffs)
            if self._Y<=0:
                [M0,Y0]=find_Y_and_M(self._space._group,self._R,nd)
                self.get_coeffs(Mset=M0,Yset=Y0,ndigs=nd,overwrite=1,norm=self._norm)
                C1 = deepcopy(self._coeffs[0][0])
            else:
                Y0 = self._Y
                M0 = get_M_for_maass_dp(self._R,Y0,10**-nd)
                C1 = deepcopy(self._coeffs[0][0])
            Y1=Y0*0.95
            self.get_coeffs(Mset=M0,Yset=Y1,ndigs=nd,overwrite=1,norm=self._norm)
            C2 = deepcopy(self._coeffs[0][0])
            self._coeffs=Ctmp
            er=RR(0)
            #print "C1.keys()=",C1.keys()
            #print "C2.keys()=",C2.keys()
            for j in range(2,floor(max([M0/2,up_to_M0,5]))):
                if self._coeffs[0][0].has_key(j):
                    t1=abs(C1[j]-self.C(j))
                    t2=abs(C2[j]-self.C(j))
                    t = max(t1,t2)
                    if self._verbose>0:
                        if t==t1:
                            print "|C1-C[{0}]|=|{1}-{2}|={3}".format(j,C1[j],self._coeffs[0][j],t)
                        else:
                            print "|C2-C[{0}]|=|{1}-{2}|={3}".format(j,C2[j],self._coeffs[0][j],t)
                elif C1.has_key(j) and C2.has_key(j):
                    t=abs(C1[j]-C2[j])
                    if self._verbose>0:
                        print "|C2-C[{0}]|=|{1}-{2}|={3}".format(j,C1[j],C2[j],t)
                        
                if t>er:
                    er=t
            if er<>0:
                d=floor(-log(er,10))
            else:
                d = 0
            if self._verbose>0:
                print "Hecke is ok up to ",d,"digits!"
            return d

    def eval(self,x,y,prec=1E-10):
        r"""
        Evaluate self.
        """
        return eval_maass_lp(self,RR(x),RR(y))
            

    def plot(self,xlim,ylim,num_pts,**kwds):
        r"""
        Make a plot of self.
        """
        # we evaluate self over a grid, for efficiency
        P=density_plot(self.eval,(-0.5,0.5),(0.01,1.01),plot_points=100,axes=False,**kwds)
        # P=density_plot(f,(-0.5,0.5),(0.01,1.01),plot_points=120,axes=False,**kwds)
        return P
        (xmin,xmax,Nx)=xlim
        (ymin,ymax,Ny)=ylim
        hy = (ymax-ymin)/RR(Ny)
        hx = (xmax-xmin)/RR(Nx)
        yvec=dict()
        #assert self._G == 
        for i in range(Nx):
            y=ymin+i*hy
            for n in range(1,M0):
                yvec[n]=besselk_dp(self.R,RR(2*pi*n))
            for j in range(Ny):
                x=xmin+i*hx                
                w=0
                argx=CC(2*pi*x*I)
                for n in range(M0):
                    term=self.C[n]*exp(argx*N)
        

 

                    
    def get_coeffs(self,Mset=0 ,Yset=None,ndigs=12,twoy=None,overwrite=False,dim=1,norm={}):
        r"""
        Compute M Fourier coefficients (at each cusp) of a Maass (cusp)
        waveform with eigenvalue R for the group G.
        
        
        INPUT:
        
        - ''S''    -- space of Maass waveforms
        - ''R''    -- Real number with precision prec
        - ''Mset'' -- integer : number of desired coefficients
        - ''ST''   -- set symmetry
        - ''Yset'' -- real
        - ''ndigs''-- integer
        -``overwrite`` --  set to True to overwrite old coefficients
        
        OUTPUT:
        
        -''D'' -- dictionary of Fourier coefficients
        
        EXAMPLES:


          sage: R=mpmath.mpf(9.53369526135355755434423523592877032382125639510725198237579046413534)
          sage:         sage: M=MaassWaveForms(Gamma0(1))
          sage:         sage: C=Maassform_coeffs(M,R)
             

        """
        S=self._space
        R=self._R
        G=S._group
        import mpmath
        if(S._verbose>1):
            print "S=",S
        dold=mpmath.mp.dps
        param=self._space.set_default_parameters(R,Mset,Yset,ndigs)
        if Yset:
            Y=Yset
        else:
            Y=param['Y']*0.5
        Q=param['Q']
        M=param['M']
        sym_type=self._sym_type
        #dim=self._dim
        set_c=self._set_c
        if norm == {}:
            norm = S.set_norm(dim)
        if S._verbose>1:
            print "R,Y,M,Q=",R,Y,M,Q
            print "sym_type=",sym_type
            print "Norm=",norm
            # print "nd=",mpmath.mp.dps
        do_cplx=1
        if S.multiplier().is_real() and sym_type in [0,1] and S._use_real:            
            do_cplx=0
        #if self._dim>1:
        #    raise NotImplementedError,"Vector-valued Maass waveforms are currently not implemented!"
        if ndigs<=15:
            if do_cplx:
                if self.weight()==0 and self.group().is_congruence():
                    X=get_coeff_fast_cplx_dp_sym(S,RR(R),RR(Y),int(M),int(Q),norm)
                else:
                    X=get_coeff_fast_cplx_dp_nosym(S,RR(R),RR(Y),int(M),int(Q),norm)
            else:
                X=get_coeff_fast_real_dp_nosym(S,RR(R),RR(Y),int(M),int(Q),norm)
            ## We still want the variables to have Sage types and not primitive python types
            for i in X.keys():
                for j in X[i].keys():
                    if hasattr(X[i][j],"keys"):
                        for n in X[i][j].keys():
                            c = CC(X[i][j][n])
                            X[i][j][n] = c
                    else:
                        try:
                            c = CC(X[i][j])
                            X[i][j] = c
                        except TypeError as te:
                            raise TypeError,"Could not coerce coeficient {0} to CC: {1}".format(X[i][j],te)
        else:
            raise NotImplementedError,"High precision is currently not (efficiently) inplemented!"
        ## The parameters used to compute the current set of coefficients.xs
        self._M0 = M
        self._Y  = Y
        self._norm = norm
        # If we compute more than one Maass form at one time we simply put the coefficients in the first component
        # And rearrange them later in the "get_element" routine.
        if overwrite==1 or dim>1:
            self._coeffs=X
            return
        if not isinstance(self._coeffs,dict):
            self._coeffs=dict()            
            self._coeffs[0]=dict()
        if self._verbose>0:
            print "X.keys()=",X.keys()
        for j in X.keys():
            if not self._coeffs[0].has_key(j):
                self._coeffs[0][j]=dict()
            for n in X[j].keys():
                if self._coeffs[0][j].has_key(n):
                    continue
                else:
                    self._coeffs[0][j][n]=X[j][n]



    def Hecke_action(self,p):
        r"""
        Return T_p(F)
        Note: Only Fourier coefficients at infinity are computed
        """
        res = copy(self)
        c = res._coeffs
        x=self._space._multiplier._character
        for r in res._coeffs.keys():
            res._coeffs[r]=dict()
            res._coeffs[r][0]=dict()
            Ms = min(self._coeffs[r][0].keys())
            Mf = max(self._coeffs[r][0].keys())
            if Ms<0:
                Ms=ceil(RR(Ms)/RR(p))
                Mf=floor(RR(Mf)/RR(p))
            for n in range(Ms,Mf+1):
                tmp=0
                if self._coeffs[r][0].has_key(n*p):
                    tmp += self._coeffs[r][0][n*p]

                if (n%p)==0:
                    m = Integer(n/p)
                    if  self._coeffs[r][0].has_key(m):
                        tmp+=x(p)*self._coeffs[r][0][m]
                res._coeffs[r][0][n]=tmp
        return res




class EisensteinSeries(AutomorphicFormElement):
    r"""
    Non-holomorphic Eisenstein series
    """
    def __init__(self,G,s,nd=12,compute=True,verbose=0):

        if hasattr(G,"_is_maass_waveform_space"):
            self._space=G
        else:
            self._space= MaassWaveForms(G,cuspidal=False)
        ### The working precision is determined by the input
        if hasattr(s,"prec"):
            prec = s.prec()
        else:
            prec = 53        
        self._prec = prec
        self._verbose = verbose
        CF = MPComplexField(self._prec)
        RF = RealField(self._prec)
        self._sigma= RF(s.real())
        self._R= RF(s.imag())
        self._s = CF(self._sigma,self._R)
        self._ndigs = nd
        self._eps = 2.0**(1-nd)
        AutomorphicFormElement.__init__(self,self._space,C=None,prec=prec,principal_part={},verbose=verbose)
        if compute:
            self.get_coefficients()
    
    def get_coefficients(self,Y0=0,M0=0):
        ## At the moment we have only implemented
        ## Eisenstein series for Hecke triangle groups.
        if not is_Hecke_triangle_group(self._space._group):
            raise NotImplementedError
        Rf = float(abs(self._R))
        if M0>0:
            Y = get_Y_for_M_dp(self._space,Rf,M0,self._eps)
        elif Y0>0:
            Y = Y0
            M = get_M_for_maass_dp(Rf,float(Y0),float(self._eps))
        else:
            Y,M = get_Y_and_M_dp(self._space,abs(self._R),self._eps)
        Ymax = self._space._group.minimal_height()/self._space._group._lambdaq
        if Y>Ymax:
            Y=0.99*Ymax
            M = get_M_for_maass_dp(float(abs(self._R)),float(Y),float(self._eps))
        RF = RealField(self._prec)
        Y = RF(Y)
        if self._verbose>0:
            print "Computing coefficients at s={0} with Y={1}, M={2}".format(self._s,Y,M)
        C = Eisenstein_series_one_cusp(self._space,self._sigma,self._R,Y,M,self._verbose)
        self._coeffs = {0: C}
        


                    
from numpy import array

def coefficients_for_Maass_waveforms(S,R,Y,M,Q,ndigs,cuspidal=True,sym_type=None,dim=1,set_c=None):
    r"""
    Compute coefficients of a Maass waveform given a specific M and Y.
     INPUT:
    
    - ''S''    -- Space of Maass waveforms
    - ''R''    -- real : 1/4+R*R is the eigenvalue
    - ''Y''    -- real number > 9
    - ''M''    -- integer
    - ''Q''    -- integer        
    - ''cuspidal''-- logical (default True)
    - ''sym_type'' -- integer (default None)
    - ''ndigs''-- integer : desired precision

    OUTPUT:

    -''D'' -- dictionary of Fourier coefficients


    EXAMPLES::

        sage: S=MaassWaveForms(Gamma0(1))
        sage: R=mpmath.mpf(9.53369526135355755434423523592877032382125639510725198237579046413534)
        sage: Y=mpmath.mpf(0.85)
        sage: C=coefficients_for_Maass_waveforms(S,R,Y,10,20,12)

    

    """
    G=S.group()
    if S._verbose>1:
        print "R,Y,M,Q,sym_type=",R,Y,M,Q,sym_type
    ## Find out which method to use. I.e. real/complex/multiprec. etc.
    #if ndigs<=12:
    #    
    if ndigs<=12:
        W=setup_matrix_for_Maass_waveforms(S,R,Y,M,Q,cuspidal=True,sym_type=sym_type,low_prec=True)
    else:
        W=setup_matrix_for_Maass_waveforms(S,R,Y,M,Q,cuspidal=True,sym_type=sym_type)
    #set_c=S._set_c
    #dim=S._dim  ## Assumed dimension of ambient space / determines how many F-coeffs. we need to set. 
    N=S.set_norm(dim,cuspidal=True)
    #return [W, N]
    dold=mpmath.mp.dps
    mpmath.mp.dps=max(dold,50)  # We work with low precision initially
    if(S._verbose>1):
        deb=True
    else:
        deb=False
    done=False; j=0
    while(done==False and j<=10):
        if(S._verbose>1):
            print "Trying to solve with prec=",mpmath.mp.dps
        try:
            X=solve_system_for_Maass_waveforms(W,N,deb=deb)
        except ZeroDivisionError:
            pass
        if(is_Integer(X) or isinstance(X,int)):
            mpmath.mp.dps=X+5
        elif(isinstance(X,dict)):
            done=True
        else:
            raise ArithmeticError," Could not solve system!"
        j=j+1
    print "X.keys=",X.keys()
    
    if(S._verbose>1):
        for m in X.keys():
            print "Function nr. ",m+1 
            for j in X[m].keys():
                if(sym_type==None):
                    for n in range(M,1 ,-1 ): 
                        print "C[",n,"]=",X[m][j][n]
                    for n in range(M): 
                        print "C[",n,"]=",X[m][j][n]
                else:
                    for n in range(1,M+1): 
                        print "C[",n,"]=",X[m][j][n]
    #print "C2=",X[0][2]
    #print "c2c3-c6=",X[0][2]*X[0][3]-X[0][6]
    mpmath.mp.dps=dold
    return X

# def verify_eigenvalue(S,R,nd=10,ST=None,method='TwoY'):
#     r""" Verify an eigenvalue and give an estimate of the error.

#     INPUT:
#     -''S'' -- Space of Maass waveforms
#     -''R'' -- real: (tentative) eigenvalue = 1/4+R**2 
#     -''nd''-- integer : number of digits we try to get (at a minimum)
#     """
    
#     C=Maassform_coeffs(S,R,ST=ST ,ndigs=nd)





# def  find_single_ev(S,R1in,R2in,Yset=None,neps=10,method='TwoY',verbose=0):
#     r""" Locate a single eigenvalue on G between R1 and R2

#     INPUT:(tentative)

#     - ''S''    -- space of Maass waveforms
#     - ''R1in'' -- real
#     - ''R1in'' -- real
#     - ''Yset'' -- real (use this value of Y to compute coefficients)
#     - ''neps'' -- number of desired digits

#     OUPUT:
    
#     - ''R'' --

    
#     """
#     G=S.group()
#     jmax=1000  # maximal number of interation
#     if(neps>=15):
#         R1=mpmath.mp.mpf(R1in);R3=mpmath.mp.mpf(R2in)
#         print "mpmath.mp.dps=",mpmath.mp.dps
#         print "R1=",R1,type(R1)
#         print "R3=",R3,type(R3)
#     else:
#         R1=mpmath.fp.mpf(R1in);R3=mpmath.fp.mpf(R2in)
#     if(Yset==None):
#         [Y,M]=find_Y_and_M(G,R1,neps)
#     else:
#         [Y,M]=find_Y_and_M(G,R1,neps,Yset=Yset)
#     Y1=Y; Y2=mpmath.mpf(0.995)*Y1
#     tol=mpmath.mpf(10)**mpmath.mpf(-neps)
#     dold=mpmath.mp.dps
#     mpmath.mp.dps=neps+3  # We work with low precision initially
#     h=dict()
#     signs=dict();diffs=dict()
#     c=dict(); h=dict()
#     c[1]=2 ; c[2 ]=3 ; c[3 ]=4 
#     #met='Hecke'

#     met=method
#     [diffs[1 ],h[1 ]]=functional(S,R1,M,Y1,Y2,signs,c,first_time=True,method=met,ndigs=neps)
#     [diffs[3 ],h[3 ]]=functional(S,R3,M,Y1,Y2,signs,c,first_time=True,method=met,ndigs=neps)
#     if S._verbose>1:
#         print "diffs: met=",met
#         print "R1=",R1
#         print "R3=",R3
#         for n in list(c.keys()): #.sort():
#             for j in list(diffs.keys()): #.sort():
#                 print "diff[",j,c[n],"]=",diffs[j][n]
#     # Sset signs and check zeros
#     if met=='Hecke':
#         if(h[1 ]*h[3]>mpmath.eps()):
#                 # We do not have a signchange
#                 return [0 ,0 ]
#     else:
#         var=0.0
#         for j in range(1 ,3 +1 ):
#             var+=abs(diffs[1 ][j])+abs(diffs[3 ][j])
#         print "var=",var
#         for j in range(1 ,3 +1 ):
#             signs[j]=1 
#             if(diffs[1 ][j]*diffs[3 ][j]>mpmath.eps()):
#                 # If we do not have a signchange
#                 # and the absolute values are relatively large
#                 # there is probably no zero here
#                 if(abs(diffs[1][j])+abs(diffs[3][j]) > 0.01*var):
#                     return [0 ,0 ]
#             elif(diffs[1 ][j]>0 ):
#                 signs[j]=-1         
#             # Recompute functionals using the signs
#         if(S._verbose>1):
#             print "h1=",h
#             print "diffs1=",diffs
#             print "signs=",signs
#         for k in [1,3]:
#             h[k]=0 
#             for j in range(1,3+1):
#                 h[k]=h[k]+signs[j]*diffs[k][j]
#     Rnew=prediction(h[1 ],h[3 ],R1,R3)
#     if S._verbose>1:
#         print "h=",h
#         print "Rnew=",Rnew
#     [diffs[2],h[2]]=functional(S,Rnew,M,Y1,Y2,signs,c,first_time=False,method=met,ndigs=neps)
#     zero_in=is_zero_in(h)
#     if(zero_in == -1 ):
#         R3=Rnew; h[3]=h[2 ]; diffs[3 ]=diffs[2 ]; errest=abs(Rnew-R1)
#     else:
#         R1=Rnew; h[1 ]=h[2 ]; diffs[1 ]=diffs[2 ]; errest=abs(Rnew-R3)
#     step=0
#     for j in range(100):
#         Rnew=prediction(h[1 ],h[3 ],R1,R3)
#         errest=max(abs(Rnew-R1),abs(Rnew-R3))
#         if S._verbose>1:
#             print "R1,R3,Rnew,errest=",R1,R3,Rnew,errest
#         if errest<tol:
#             return [Rnew,errest]
#         [diffs[2 ],h[2 ]]=functional(S,Rnew,M,Y1,Y2,signs,c,first_time=False,method=met,ndigs=neps)
#         zero_in=is_zero_in(h)
#         if zero_in==0:
#             return [Rnew,errest]
#         elif zero_in not in [1,-1]:
#             raise StopIteration()
#         if zero_in==-1:
#             stepz=abs(Rnew-R3)
#             R3=Rnew; h[3 ]=h[2 ]; diffs[3 ]=diffs[2 ]; errest=abs(Rnew-R1)
#         elif zero_in==1:
#             stepz=abs(Rnew-R1)
#             R1=Rnew; h[1 ]=h[2 ]; diffs[1 ]=diffs[2 ]; errest=abs(Rnew-R3)
#         # If we have gone in the same direction too many times we need to modify our approach
#         step=step+zero_in
#         if S._verbose>1:
#             print "step=",step
#         if step>2:    # Have gone too many times to the left
#             Rtest=Rnew + mpmath.mpf(0.5)*stepz  # Need to test if this modified R3 work:
#             if S._verbose>1:
#                 print "Rtest(R)=",Rtest
#             [diffs[2 ],h[2 ]]=functional(S,Rtest,M,Y1,Y2,signs,c,False,met,neps)
#             if is_zero_in(h) ==-1: # all is ok
#                 R3=Rtest; h[3]=h[2]; diffs[3]=diffs[2]; step=step-1
#             else: # Test another one
#                 Rtest=Rnew + mpmath.mpf(0.5)*abs(R1-R3)  # Need to test if this modified R3 work:
#                 if S._verbose>1:
#                     print "Rtest(R)=",Rtest
#                 [diffs[2 ],h[2 ]]=functional(S,Rtest,M,Y1,Y2,signs,c,False,met,neps)
#                 if is_zero_in(h) ==-1: # all is ok
#                     R3=Rtest; h[3]=h[2]; diffs[3]=diffs[2]; step=step-1
#         elif step<-2: # Have gone too many times to the right
#             Rtest=Rnew - mpmath.mpf(0.5)*stepz
#             if S._verbose>1:
#                 print "Rtest(L)=",Rtest
#             [diffs[2 ],h[2 ]]=functional(S,Rtest,M,Y1,Y2,signs,c,False,met,neps)
#             if is_zero_in(h) == 1: # all is ok
#                 R1=Rtest; h[1]=h[2]; diffs[1]=diffs[2]; step=step+1
#             else:
#                 Rtest=Rnew - mpmath.mpf(0.5)*abs(R3-R1)
#                 if S._verbose>1:
#                     print "Rtest(L)=",Rtest
#                 [diffs[2 ],h[2 ]]=functional(S,Rtest,M,Y1,Y2,signs,c,False,met,neps)
#                 if is_zero_in(h) == 1: # all is ok
#                     R1=Rtest; h[1]=h[2]; diffs[1]=diffs[2]; step=step+1
                
# ####
            
def is_zero_in(h):
    r"""
    Tells which interval contains changes of sign.
    

    INPUT:
    
    - ''h'' -- dictionary h[j]=[h1,h2,h3]

    OUTPUT:
    
    - integer  (-1 0 1)

    EXAMPLES::
    
    
    """
    zi=dict(); i=0
    if(h[1]*h[2] < 0):
        zi[-1]=1; i=-1
    if(h[3]*h[2] < 0):
        zi[1]=1; i=1
    if(zi.values().count(1) >1 ): # need to split
        return -2 
    #s="Neeed to split! Not implemented!"
    #raise ValueError,s
    return i



def get_character_sqrt(x):
    if isinstance(x,sage.modular.dirichlet.DirichletCharacter):
        if x.is_even():
            for y in x.parent().list():
                if y*y == x:
                    return y
    raise ValueError,"Need an even character to get a square root of a character! Got:{0}".format(x)

def prediction(f0,f1,x0,x1):
    r"""
    Predict zero using the secant method.

    INPUT:

        - ''f0'' -- real
        - ''f1'' -- real
        - ''x0'' -- real
        - ''x1'' -- real

    OUTPUT:

        - real


    EXAMPLES::


        sage: prediction(-1,1,9,10)
        19/2
        sage: prediction(-1.0,1.0,9.0,10.0)
        9.50000000000000

    
    """
    xnew=x0-f0*(x1-x0)/(f1-f0)
    #if(xnew<x0 or xnew>x1):
    #    st= "Secant method ended up outside interval! \n"
    #    st+="input: f0,f1,x0,x1=%s,%s,%s,%s \n xnew=%s"
    #    raise ValueError,st%(f0,f1,x0,x1,xnew)
    return xnew

def prediction_newton(x,f,df):
    r"""
    Predict zero using the secant method.

    INPUT:

        - ''x''  -- real
        - ''f''  -- real, f(x)
        - ''df'' -- real, f'(x)

    OUTPUT:

        - real


    EXAMPLES::


        sage: prediction_newton(1,9,10)
        19/2
        sage: prediction_newton(1.0,9.0,10.0)
        9.50000000000000

    
    """
    if(df==0.0):
        st= "Newtons method failed! \n"
        st+="f'(x)=0. input: f,df,x=%s,%s,%s"
        raise ValueError,st%(f,df,x)
    xnew=x-f/df
    return xnew



def find_Y_and_M(G,R,ndigs=12,Yset=None,Mset=None):
    r"""
    Compute a good value of M and Y for Maass forms on G

    INPUT:

    - ''G'' -- group
    - ''R'' -- real
    - ''ndigs'' -- integer (number of desired digits of precision)
    - ''Yset'' -- real (default None) if set we return M corr. to this Y
    - ''Mset'' -- integer (default None) if set we return Y corr. to this M

    OUTPUT:

    - [Y,M] -- good values of Y (real) and M (integer)

    EXAMPLES::

    

    TODO:
    Better and more effective bound
    """

    import mpmath
    l=G.generalised_level()
    if(Mset <> None):
        # then we get Y corr. to this M
        Y0=RR(3).sqrt()/RR(2*l)
        
    if(Yset==None):
        Y0=RR(3).sqrt()/RR(2*l)
        Y=mpmath.fp.mpf(0.95*Y0)
    else:
        Y=mpmath.fp.mpf(Yset)
    #print "Y=",Y,"Yset=",Yset
    IR=mpmath.mpc(0,R)
    eps= mpmath.fp.mpf(10 **-ndigs)
    twopiY=mpmath.fp.pi*Y*mpmath.fp.mpf(2)

    M0=get_M_for_maass(R,Y,eps) 
    if(M0<10):
        M0=10
    ## Do this in low precision
    dold=mpmath.mp.dps
    #print "Start M=",M0
    #print "dold=",dold
    #mpmath.mp.dps=100
    try:
        for n in range(M0,10000,3):
            X=mpmath.pi()*Y*mpmath.mpf(2*n)
            #print "X,IR=",X,IR
            test=mpmath.fp.besselk(IR,X)
            if(abs(test)<eps):
                raise StopIteration()
    except StopIteration:
        M=n
    else:
        M=n
        raise Exception,"Error: Did not get small enough error:=M=%s gave err=%s" % (M,test)
    mpmath.mp.dps=dold
    return [Y,M]







def _testing_kbes(Rt=[1,10,10],Xt=[1,10,100]):
    [R0,R1,NR]=Rt
    [X0,X1,NX]=Xt
    NRr=mpmath.mpf(NR)
    NXr=mpmath.mpf(NX)
    for j in range(1,NR):
        rj=mpmath.mpf(j)
        R=R0+R1*rj/NRr
        print "r=",R
        iR=mpmath.mpc(0,R)
        for k in range(1,NX):
            rk=mpmath.mpf(k)
            x=X0+X1*rk/NXr
            print "r,x=",R,x
            if(x>R):
                print "kbes_pow="
                timeit( "besselk_dp(%s,%s)" %(R,x),repeat=1)
            #else:
            #    print "kbes_rec="
            #    timeit( "besselk_dp_rec(R,x)",repeat=1)
            print "mpmath.besselk="
            timeit("mpmath.besselk(%s,%s)" %(iR,x),repeat=1)
            

            #print "t1(",R,x,")=",t1
            #print "t2(",R,x,")=",t2
            if(R<15.0):
                if(x<0.3 *R):
                    print "Case 1"
                elif(x<=max(10.0 +1.2*R,2 *R)):
                    print "Case 2"
            elif(R>20  and x>4 *R):
                print "Case 3"
            else:
                print "Case 4"


def _test_Hecke_relations(a=2,b=3,C={}):
    r"""Testing Hecke relations for the Fourier coefficients in C

    INPUT:
    -''C'' -- dictionary of complex (Fourier coefficients)
    -''a'' -- integer
    -''b'' -- integer

    OUTPUT:
    -''diff'' -- real : |C(a)C(b)-C(ab)| if (a,b)=1

    EXAMPLE::

    
    sage: S=MaassWaveForms(Gamma0(1))
    sage: R=mpmath.mpf(9.53369526135355755434423523592877032382125639510725198237579046413534)
    sage: Y=mpmath.mpf(0.85)
    sage: C=coefficients_for_Maass_waveforms(S,R,Y,10,20,12)
    sage: d=_test_Hecke_relations(C,2,3); mppr(d)
    '9.29e-8'
    sage: C=coefficients_for_Maass_waveforms(S,R,Y,30,50,20)
    sage: d=_test_Hecke_relations(C,2,3); mppr(d)
    '3.83e-43'
    
    
    """
    c=gcd(Integer(a),Integer(b))
    if not C.has_key(0):
        return 0
    if C[0].has_key(a) and C[0].has_key(b) and C[0].has_key(a*b): 
        lhs=C[0][a]*C[0][b]
        rhs=0
        for d in divisors(c):
            rhs=rhs+C[0][Integer(a*b/d/d)]
        return abs(rhs-lhs)
    return 0

def _test_Hecke_relations_all(C={}):
    r"""
    Test all possible Hecke relations.


    EXAMPLE::

    
        sage: S=MaassWaveForms(Gamma0(1))
        sage: mpmath.mp.dps=100
        sage: R=mpmath.mpf(9.53369526135355755434423523592877032382125639510725198237579046413534899129834778176925550997543536649304476785828585450706066844381418681978063450078510030977880577576)
        sage: Y=mpmath.mpf(0.85)
        sage: C=coefficients_for_Maass_waveforms(S,R,Y,30,50,20)
        sage: test=_test_Hecke_relations_all(C); test  
        {4: '9.79e-68', 6: '4.11e-63', 9: '4.210e-56', 10: '9.47e-54', 14: '2.110e-44', 15: '4.79e-42', 21: '4.78e-28', 22: '1.02e-25', 25: '9.72e-19', 26: '2.06e-16'}
        

    We can see how the base precision used affects the coefficients

        sage: mpmath.mp.dps=50
        sage: C=coefficients_for_Maass_waveforms(S,R,Y,30,50,20)
        sage: test=_test_Hecke_relations_all(C); test               
        sage: test=_test_Hecke_relations_all(C); test  
{4: '1.83e-48', 6: '4.75e-43', 9: '3.21e-36', 10: '1.24e-33', 14: '4.41e-25', 15: '1.53e-23', 21: '6.41e-8', 22: '4.14e-6', 25: '91.455', 26: '12591.0'}

    

    """
    N=max(C[0].keys())
    test=dict()
    for a in prime_range(N):
        for b in prime_range(N):
            if(a*b <= N):
                test[a*b]=mppr(_test_Hecke_relations(C,a,b))
    return test


def solve_system_for_Maass_waveforms(W,N=None,deb=False,force_type=None):
    r"""
    Choose the correct solver algorithm.
    """
    x=W['V'][0,0]
    if force_type=="mpc" and hasattr(W['V'],"ctx"):
        A=mat_conv_to_mpc(W['V'])
        W['V']=A
    elif force_type=="mpmath" and hasattr(W['V'],"QR"):
        A=mat_conv_to_mpmath(W['V'])
        W['V']=A
    if hasattr(W['V'],"ctx"): ## We have an mpmath matrix
        ## Recall that we may need to adjust the precision for mpmath solutions   
        return solve_system_for_Maass_waveforms_mpmath(W,N,deb)
    elif hasattr(W['V'],"QR"): # Is a Matirx_complex_dense instance
        return solve_system_for_Maass_waveforms_mpc(W,N,deb)
    else:
        raise ValueError,"Unknown type of matrix!: {0}".format(type(W['V']))
        
def solve_system_for_Maass_waveforms_mpmath(W,N=None,deb=False,gr=False):
    r"""
    Solve the linear system to obtain the Fourier coefficients of Maass forms

    INPUT:

    - ``W`` --   (system) dictionary

        - ``W['Ms']``  -- M start
        - ``W['Mf']``  -- M stop
        - ``W['nc']``  -- number of cusps
        - ``W['V']``   -- matrix of size ((Ms-Mf+1)*nc)**2
        - ``W['RHS']`` -- right hand side (for inhomogeneous system) matrix of size ((Ms-Mf+1)*nc)*(dim)

    - ``N`` -- normalisation (dictionary, output from the set_norm_for_maass function)

        - ``N['SetCs']``   -- Which coefficients are set
        - ``N['Vals'] ``   -- To which values are these coefficients set
        - ``N['comp_dim']``-- How large is the assumed dimension of the solution space
        
    - ``deb`` -- print debugging information (default False)

    OUTPUT:
    
    - ``C`` -- Fourier coefficients

    EXAMPLES::

        sage: S=MaassWaveForms(MySubgroup(Gamma0(1)))
        sage: mpmath.mp.dps=20
        sage: R=mpmath.mpf(9.533695261353557554344235235928770323821256395107251982375790464135348991298347781769255509975435366)
        sage: Y=mpmath.mpf(0.5)
        sage: W=setup_matrix_for_Maass_waveforms(S,R,Y,12,22)
        sage: N=S.set_norm_maass(1)
        sage: C=solve_system_for_Maass_waveforms(W,N)
        sage: C[0][2]*C[0][3]-C[0][6]
        mpc(real='-1.8055426724989656270259e-14', imag='1.6658248366482944572967e-19')

    If M is too large and the precision is not high enough the matrix might be numerically singular

        W=setup_matrix_for_Maass_waveforms(G,R,Y,20,40)  
        sage: C=solve_system_for_Maass_waveforms(W,N)
        Traceback (most recent call last)
        ...
        ZeroDivisionError: Need higher precision! Use > 23 digits!

    Increasing the precision helps
    
        sage: mpmath.mp.dps=25
        sage: R=mpmath.mpf(9.533695261353557554344235235928770323821256395107251982375790464135348991298347781769255509975435366)
        sage: C=solve_system_for_Maass_waveforms(W,N)
        sage: C[0][2]*C[0][3]-C[0][6]
        mpc(real='3.780824715556976438911480324e-25', imag='2.114746048869188750991752872e-99')

        
    """
    import mpmath
    V=W['V']
    Ms=W['Ms']
    Mf=W['Mf']
    nc=W['nc']
    M=W['space']
    verbose=M._verbose
    Ml=Mf-Ms+1
    if N==None:
        N = M.set_norm()
    if hasattr(V,'shape'):
        nrows,ncols=V.shape
    else:
        if hasattr(V,"rows"):
            nrows=V.rows
            ncols=V.cols
    if(ncols<>Ml*nc or nrows<>Ml*nc):
        raise Exception," Wrong dimension of input matrix!"
    if M._verbose>0:
        print "Norm=",N
    SetCs=N['SetCs'][0]
    Vals=N['Vals']
    comp_dim=N['comp_dim']
    if(N['cuspidal']):
        for i in range(1,nc):
            if(SetCs.count((i,0))==0):
                SetCs.append((i,Ml))
            for fn_j in range(comp_dim):
                Vals[fn_j][(i,Ml)]=0
    setc_list=list()
    vals_list=dict()
    for j in range(comp_dim):
        vals_list[j]=dict()
    
    for r,n in SetCs:
        if r*Ml+n-Ms<0:
            continue
        setc_list.append(r*Ml+n-Ms)
        for j in range(comp_dim):
            vals_list[j][r*Ml+n-Ms]=Vals[j][(r,n)]
    if verbose>0:
        print "setc_list=",setc_list
        print "vals_list=",vals_list
    if(Ms<0):
        use_sym=0
    else:
        use_sym=1
    #if(use_sym==1 and SetCs.count(0)>0):
    #    num_set=len(N['SetCs'])-1
    #else:
    num_set=len(setc_list)
    t=V[0,0]
    if(isinstance(t,float)):
        mpmath_ctx=mpmath.fp
    else:  
        mpmath_ctx=mpmath.mp
    if(W.has_key('RHS')):
        RHS=W['RHS']
    else:
        RHS=mpmath_ctx.matrix(int(Ml*nc-num_set),int(comp_dim))
    LHS=mpmath_ctx.matrix(int(Ml*nc-num_set),int(Ml*nc-num_set))
    roffs=0

    if(deb):
        print "num_set,use_sym=",num_set,use_sym
        print "SetCs,Vals=",SetCs,Vals
        print "V.rows,cols=",nrows,ncols
        print "LHS.rows,cols=",LHS.rows,LHS.cols
        print "RHS.rows,cols=",RHS.rows,RHS.cols
        print "mpctx=",mpmath_ctx
    for r in range(nrows):
        #cr=r+Ms
        if setc_list.count(r)>0:
            roffs=roffs+1
            continue
        for fn_j in range(comp_dim):
            RHS[r-roffs,fn_j]=mpmath_ctx.mpf(0)
            for cset in setc_list:
                v=vals_list[fn_j][cset]
                if(mpmath_ctx==mpmath.mp):
                    tmp=mpmath_ctx.mpmathify(v)
                elif(isinstance(v,float)):
                    tmp=mpmath_ctx.mpf(v)
                else:
                    tmp=mpmath_ctx.mpc(v)
                tmp=tmp*V[r,cset]
                RHS[r-roffs,fn_j]=RHS[r-roffs,fn_j]-tmp
                #print "RHS[",r-roffs,fn_j,"]=",RHS[r-roffs,fn_j]
                #print "V[",r,",",cset,"]=",V[r,cset]
        coffs=0
        for k in range(ncols):            
            if setc_list.count(k)>0:
                coffs=coffs+1
                continue
            #print "roffs,coffs=",roffs,coffs
            #print "r-roffs,k-coffs=",r-roffs,k-coffs
            LHS[r-roffs,k-coffs]=V[r,k]
            #print "LHS[",r,k,"]=",LHS[r-roffs,k-coffs]
    if gr:
        return LHS,RHS
    done=0; j=0
    oldprec=mpmath.mp.dps
    while done==0 and j<=10:
        if W['space']._verbose>1:
             print "Trying to solve with prec=",mpmath.mp.dps
        try:
            A, p = mpmath_ctx.LU_decomp(LHS)
            done=1
        except ZeroDivisionError:
            t1=smallest_inf_norm(LHS)
            if verbose>0:
                print "n=",smallest_inf_norm(LHS)
            t2=mpmath_ctx.log10(smallest_inf_norm(LHS))
            t3=mpmath_ctx.ceil(-t2)
            isinf=False
            if(isinstance(t3,float)):
                isinf = (t3 == float(infinity))
            if(isinstance(t3,sage.libs.mpmath.ext_main.mpf)):
                isinf = ((t3.ae(mpmath.inf)) or t3==mpmath.inf)
            if(isinstance(t3,sage.rings.real_mpfr.RealLiteral)):
                isinf = t3.is_infinity()
            if(isinf):
                raise ValueError, " element in LHS is infinity! t3=%s" %t3
            t=int(t3)
            mpmath.mp.dps= t + 5
    X=dict()
    for fn_j in range(comp_dim):
        X[fn_j] = dict() #mpmath.matrix(int(Ml),int(1))
        b = mpmath_ctx.L_solve(A, RHS.column(fn_j), p)
        TMP = mpmath_ctx.U_solve(A, b)
        roffs=0
        res = mpmath_ctx.norm(mpmath_ctx.residual(LHS, TMP, RHS.column(fn_j)))
        #print "res(",fn_j,")=",res
        for i in range(nc):
            X[fn_j][i]=dict()
        for n in range(Ml):
            if setc_list.count(n)>0:
                roffs=roffs+1
                #print "X[",fn_j,",",n,",Vals[fn_j][n]
                X[fn_j][0][n+Ms]=vals_list[fn_j][n]
                continue
            X[fn_j][0][n+Ms]=TMP[n-roffs,0]
            #print "X[",fn_j,",",n+Ms,"=",TMP[n-roffs,0]
        for i in range(1,nc):
            for n in range(Ml):
                if setc_list.count(n+i*Ml)>0:
                    #(SetCs.count(n+Ms+i*Ml)>0):
                    roffs=roffs+1
                    # print "X[",fn_j,",",n,",Vals[fn_j][n]
                    X[fn_j][i][n+Ms]=vals_list[fn_j][n+i*Ml]
                    continue
                X[fn_j][i][n+Ms]=TMP[n+i*Ml-roffs,0]
    # return x
    #print "keys:",X.keys()
    mpmath.mp.dps=oldprec
    return X


def solve_system_for_Maass_waveforms_mpc(W,N=None,gr=False,cn=False):
    r"""
    Solve the linear system to obtain the Fourier coefficients of Maass forms

    INPUT:

    - ``H`` -- Space of Maass waveforms
    - ``W`` --   (system) dictionary
        - ``W['Ms']``  -- M start
        - ``W['Mf']``  -- M stop
        - ``W['nc']``  -- number of cusps
        - ``W['V']``   -- matrix of size ((Ms-Mf+1)*nc)**2
        - ``W['RHS']`` -- right hand side (for inhomogeneous system) matrix of size ((Ms-Mf+1)*nc)*(dim)
    - ``N`` -- normalisation (dictionary, output from the set_norm_for_maass function, default None)
        - if N=None we assume that the solution is uniquely determined by the prinicpal part (in the right hand side)
        - ``N['SetCs']``   -- Which coefficients are set
        - ``N['Vals'] ``   -- To which values are these coefficients set
        - ``N['comp_dim']``-- How large is the assumed dimension of the solution space
        - ``N['num_set']`` -- Number of coefficients which are set
        

        
    - ''gr''  -- only return the reduced matrix and right hand side. do not perform the solving .
    - ''cn''  -- logical (default False) set to True to compute the max norm of V^-1
    OUTPUT:
    
    - ``C`` -- Fourier coefficients

    EXAMPLES::

        sage: G=MySubgroup(Gamma0(1))
        sage: mpmath.mp.dps=20
        sage: R=mpmath.mpf(9.533695261353557554344235235928770323821256395107251982375790464135348991298347781769255509975435366)
        sage: Y=mpmath.mpf(0.5)
        sage: W=setup_matrix_for_Maass_waveforms(G,R,Y,12,22)
        sage: N=set_norm_maass(1)
        sage: C=solve_system_for_Maass_waveforms(W,N)
        sage: C[0][2]*C[0][3]-C[0][6]
        mpc(real='-1.8055426724989656270259e-14', imag='1.6658248366482944572967e-19')

    If M is too large and the precision is not high enough the matrix might be numerically singular

        W=setup_matrix_for_Maass_waveforms(G,R,Y,20,40)  
        sage: C=solve_system_for_Maass_waveforms(W,N)
        Traceback (most recent call last)
        ...
        ZeroDivisionError: Need higher precision! Use > 23 digits!

    Increasing the precision helps
    
        sage: mpmath.mp.dps=25
        sage: R=mpmath.mpf(9.533695261353557554344235235928770323821256395107251982375790464135348991298347781769255509975435366)
        sage: C=solve_system_for_Maass_waveforms(W,N)
        sage: C[0][2]*C[0][3]-C[0][6]
        mpc(real='3.780824715556976438911480324e-25', imag='2.114746048869188750991752872e-99')

    """
    V=W['V']
    Ms=W['Ms']
    Mf=W['Mf']
    H=W['space']
    #nc=W['nc']
    Ml=Mf-Ms+1
    get_reduced_matrix=gr
    verbose = H._verbose
    comp_norm=cn
    nc=H.group().ncusps() 
    if V.ncols()<>Ml*nc or V.nrows()<>Ml*nc:
        raise Exception," Wrong dimension of input matrix!"
    if N==None:
        N = H.set_norm(1)
    SetCs=N['SetCs'][0]
    Vals=N['Vals']
    comp_dim=N['comp_dim']
    num_set=len(SetCs[0])
    t=V[0,0]
    CF=MPComplexField(H._prec)
    MS = MatrixSpace(CF,int(Ml*nc-num_set),int(comp_dim))
    RHS=Matrix_complex_dense(MS,0,True,True)
    MS = MatrixSpace(CF,int(Ml*nc-num_set),int(Ml*nc-num_set))
    LHS=Matrix_complex_dense(MS,0,True,True)
    nrows=V.nrows()
    ncols=V.ncols()
    if(N['cuspidal']):
        for i in range(1,nc):
            if(SetCs.count((i,0))==0):
                SetCs.append((i,Ml))
            for fn_j in range(comp_dim):
                Vals[fn_j][(i,Ml)]=0
    if verbose>0:
        print "SetCs=",SetCs
    setc_list=list()
    vals_list=dict()
    for j in range(comp_dim):
        vals_list[j]=dict()    
    for r,n in SetCs:
        if r*Ml+n-Ms<0:
            continue
        setc_list.append(r*Ml+n-Ms)
        for j in range(comp_dim):
            vals_list[j][r*Ml+n-Ms]=Vals[j][(r,n)]

    if verbose>0:
        print "Ml=",Ml
        print "num_set=",num_set
        print "SetCs=",SetCs
        print "Vals=",Vals
        print "setc_list=",setc_list
        print "vals_list=",vals_list
        print "V.rows=",V.nrows()
        print "V.cols=",V.ncols()
        print "LHS.rows=",LHS.nrows()
        print "LHS.cols=",LHS.ncols()
        print "RHS.rows=",RHS.nrows()
        print "RHS.cols=",RHS.ncols()
        print "N=",N

    num_rhs=0
    if(W.has_key('RHS')):
        num_rhs=W['RHS'].ncols()
    if num_rhs>0 and num_rhs<>comp_dim:
        raise ValueError,"Need same number of right hand sides (or just one) as the number of set coefficients!"
    if V.nrows() <> nc*Ml:
        raise ArithmeticError," Matrix does not have correct size!"
    roffs=0
    for r in range(nrows):
        #cr=r+Ms
        if setc_list.count(r)>0:
            roffs=roffs+1
            continue
        for fn_j in range(comp_dim):
            if W.has_key('RHS'):
                RHS[r-roffs,fn_j]=-W['RHS'][r,rhs_j]
            else:
                RHS[r-roffs,fn_j]=CF(0)
            for cset in setc_list:
                tmp=CF(vals_list[fn_j][cset])
                tmp=tmp*V[r,cset]
                RHS[r-roffs,fn_j]=RHS[r-roffs,fn_j]-tmp
        coffs=0
        for k in range(ncols):            
            if setc_list.count(k)>0:
                coffs=coffs+1
                if verbose>1:
                    print "skipping colum:",k
                continue
                if verbose>1 and r-roffs==1:
                    print "Setting LHS[1,",k-coffs
            LHS[r-roffs,k-coffs]=V[r,k]    # for a in range(nc):
    if get_reduced_matrix:
#        return [LHS,RHS]
        return [LHS,RHS]
    maxit=100;i=0
    done=False
    dps0=CF.prec()
    # while (not done and i<=maxit):
    #     try:
    #         Q,R = LHS.qr_decomposition()
    #         done=True
    #     except ZeroDivisionError:
    #         t=int(ceil(-log_b(smallest_inf_norm(LHS),10)))
    #         dps=t+5*i; i=i+1
    #         print "raising number of digits to:",dps
    #         LHS.set_prec(dps)
    #         # raise ZeroDivisionError,"Need higher precision! Use > %s digits!" % t
    # if i>=maxit:
    #     raise ZeroDivisionError,"Can not raise precision enough to solve system! Should need > %s digits! and %s digits was not enough!" % (t,dps)
    if comp_norm:
        max_norm=LHS.norm()
        for j in range(LHS.rows):
            #y=mpmath_ctx.matrix(LHS.rows,int(1)); y[j,0]=1
            y = Vector_complex_dense(vector(CF,LHS.rows).parent(),0)
            y[j]=1
            TMP = LHS.solve(y) #pmath_ctx.U_solve(A, b)
            tmpnorm=max(map(abs,TMP))
            if(tmpnorm>max_norm):
                max_norm=tmpnorm
        print "max norm of V^-1=",max_norm
    X=dict()
    for fn_j in range(comp_dim):
        X[fn_j] = dict() #mpmath.matrix(int(Ml),int(1))
        v = RHS.column(fn_j)
        if verbose>1:
            print "len(B)=",len(v)
        TMP = LHS.solve(v)
        roffs=0
        res = (LHS*TMP-v).norm()
        if verbose>0:
            print "res(",fn_j,")=",res
        for i in range(nc):
            X[fn_j][i]=dict()
        for i in range(nc):
            for n in range(Ml):
                if setc_list.count(n+i*Ml)>0:
                    roffs=roffs+1
                    X[fn_j][i][n+Ms]=vals_list[fn_j][n+i*Ml]
                    continue
                X[fn_j][i][n+Ms]=TMP[n+i*Ml-roffs]
    return X





    def Hecke_action(self,p):
        r"""
        Return T_p(F)
        Note: Only Fourier coefficients at infinity are computed
        """
        res = copy(self)
        c = res._coeffs
        x=self._space._multiplier._character
        for r in res._coeffs.keys():
            res._coeffs[r]=dict()
            res._coeffs[r][0]=dict()
            Ms = min(self._coeffs[r][0].keys())
            Mf = max(self._coeffs[r][0].keys())
            if Ms<0:
                Ms=ceil(RR(Ms)/RR(p))
            Mf=floor(RR(Mf)/RR(p))
            for n in range(Ms,Mf+1):
                tmp=0
                if self._coeffs[r][0].has_key(n*p):
                    tmp += self._coeffs[r][0][n*p]

                if (n%p)==0:
                    m = Integer(n/p)
                    if  self._coeffs[r][0].has_key(m):
                        tmp+=x(p)*self._coeffs[r][0][m]
                res._coeffs[r][0][n]=tmp
        return res





def solve_system_for_Maass_waveforms_GaussElim(W,N=None,gr=False,cn=False):
    r"""
    Solve the linear system to obtain the Fourier coefficients of Maass forms

  
    """
    V=W['V']
    Ms=W['Ms']
    Mf=W['Mf']
    H=W['space']
    Ml=Mf-Ms+1
    get_reduced_matrix=gr
    verbose = H._verbose
    comp_norm=cn
    nc=H.group().ncusps() 
    if V.ncols()<>Ml*nc or V.nrows()<>Ml*nc:
        raise Exception," Wrong dimension of input matrix!"
    if N==None:
        N = H.set_norm(1)
    SetCs=N['SetCs'][0]
    Vals=N['Vals']
    comp_dim=N['comp_dim']
    num_set=len(SetCs[0])
    t=V[0,0]
    CF=MPComplexField(H._prec)
    MS = MatrixSpace(CF,int(Ml*nc-num_set),int(comp_dim))
    RHS=Matrix_complex_dense(MS,0,True,True)
    MS = MatrixSpace(CF,int(Ml*nc-num_set),int(Ml*nc-num_set))
    LHS=Matrix_complex_dense(MS,0,True,True)
    nrows=V.nrows()
    ncols=V.ncols()
    if(N['cuspidal']):
        for i in range(1,nc):
            if(SetCs.count((i,0))==0):
                SetCs.append((i,Ml))
            for fn_j in range(comp_dim):
                Vals[fn_j][(i,Ml)]=0
    if verbose>0:
        print "SetCs=",SetCs
    setc_list=list()
    vals_list=dict()
    for j in range(comp_dim):
        vals_list[j]=dict()    
    for r,n in SetCs:
        if r*Ml+n-Ms<0:
            continue
        setc_list.append(r*Ml+n-Ms)
        for j in range(comp_dim):
            vals_list[j][r*Ml+n-Ms]=Vals[j][(r,n)]
    if verbose>0:
        print "Ml=",Ml
        print "num_set=",num_set
        print "SetCs=",SetCs
        print "Vals=",Vals
        print "setc_list=",setc_list
        print "vals_list=",vals_list
        print "N=",N

    num_rhs=0
    if(W.has_key('RHS')):
        num_rhs=W['RHS'].ncols()
    if num_rhs>0 and num_rhs<>comp_dim:
        raise ValueError,"Need same number of right hand sides (or just one) as the number of set coefficients!"
    if V.nrows() <> nc*Ml:
        raise ArithmeticError," Matrix does not have correct size!"
    roffs=0
    for r in range(nrows):
        #cr=r+Ms
        if setc_list.count(r)>0:
            roffs=roffs+1
            continue
        for fn_j in range(comp_dim):
            if W.has_key('RHS'):
                RHS[r-roffs,fn_j]=-W['RHS'][r,rhs_j]
            else:
                RHS[r-roffs,fn_j]=CF(0)
            for cset in setc_list:
                tmp=CF(vals_list[fn_j][cset])
                tmp=tmp*V[r,cset]
                RHS[r-roffs,fn_j]=RHS[r-roffs,fn_j]-tmp
        coffs=0
        for k in range(ncols):            
            if setc_list.count(k)>0:
                coffs=coffs+1
                continue
            LHS[r-roffs,k-coffs]=V[r,k]    # for a in range(nc):
    if get_reduced_matrix:
#        return [LHS,RHS]
        return [LHS,RHS]
    maxit=100;i=0
    dps0=CF.prec()
    X=dict()
    for fn_j in range(comp_dim):
        X[fn_j] = dict() #mpmath.matrix(int(Ml),int(1))
        b = RHS.column(fn_j)
        done=0
        while (not done and i<=maxit):
            try:
                C=solve_using_Gauss_elem(LHS,b)
                done=1
            except ZeroDivisionError:
                pass
        if verbose>0:
            res = (LHS*C-b).norm()
            print "res(",fn_j,")=",res
        for i in range(nc):
            X[fn_j][i]=dict()
        for i in range(nc):
            roffs=0
            for n in range(Ml):
                if setc_list.count(n+i*Ml)>0:
                    roffs=roffs+1
                    X[fn_j][i][n+Ms]=vals_list[fn_j][n+i*Ml]
                    continue
                #print "C[",n+i*Ml-roffs,"=",C[n+i*Ml-roffs]
                X[fn_j][i][n+Ms]=C[n+i*Ml-roffs]
    return X
 


def mat_conv_to_mpc(A):
    if hasattr(A,"QR"):
        return A
    if hasattr(A,"ctx"):
        m=A.rows
        n=A.cols
        prec=mpmath.mp.prec
    elif hasattr(A,"nrows"):
        m=A.nrows()
        n=A.ncols()
        prec=A[0,0].parent().prec()
    else:
        raise TypeError,"Cabn not convert matrix of type:{0}".format(type(A))
    CF=MPComplexField(prec)
    MS=MatrixSpace(CF,m,n)
    V=Matrix_complex_dense(MS,0)    
    for i in range(m):
        for j in range(n):
            tmp=A[i,j]
            V[i,j]=CF(tmp.real,tmp.imag)
    return V


def mat_conv_to_mpmath(A):
    if hasattr(A,"ctx"):
        return A
    if hasattr(A,"nrows"):
        m=A.nrows()
        n=A.ncols()
    old_prec=mpmath.mp.prec
    mpmath.mp.prec=A[0,0].parent().prec()
    V=mpmath.mp.matrix(int(m),int(n))
    for i in range(m):
        for j in range(n):
            tmp=A[i,j]
            V[i,j]=mpmath.mp.mpc(tmp.real(),tmp.imag())
    mpmath.mp.prec=old_prec
    return V


def mat_conv_to_complex(A):
    if hasattr(A,"ctx"):
        return A
    if hasattr(A,"nrows"):
        m=A.nrows()
        n=A.ncols()
    old_prec=mpmath.mp.prec
    mpmath.mp.prec=A[0,0].parent().prec()
    V=mpmath.mp.matrix(int(m),int(n))
    for i in range(m):
        for j in range(n):
            tmp=A[i,j]
            V[i,j]=mpmath.mp.mpc(tmp.real(),tmp.imag())
    mpmath.mp.prec=old_prec
    return V



def my_kbes(r,x,mp_ctx=None):
    r"""Scaled K-Bessel function with

    INPUT:

        - ''r'' -- real
        - ''x'' -- real
        - ''mp_ctx'' -- mpmath context (default None)

    OUTPUT:

        - real -- K_ir(x)*exp(pi*r/2)

    EXAsMPLES::


        sage: my_kbes(9.0,1.0)
        mpf('-0.71962866121965863')
        sage: my_kbes(9.0,1.0,mpmath.fp)
        -0.71962866121967572

    

    """
    import mpmath
    if mp_ctx==None or mp_ctx==mpmath.fp or mpamth.dps<=15:
        # use fast routine
        return besselk_dp(RR(r),RR(x))
    else:
        pi=mpmath.mp.pi()
        k=mp_cyx.besselk(mp_ctx.mpc(0,r),mp_ctx.mpf(x))
        f=k*mp_ctx.exp(r*mp_ctx.mpf(0.5)*pi)
        return f.real




def my_kbes_diff_r(r,x,mp_ctx=None):
    r"""

    Approximation to the derivative with respect to R of the scaled K-Bessel function.

    INPUT:

        - ''r'' -- real
        - ''x'' -- real
        - ''ctx'' -- mpmath context (default mpmath.mp)

    OUTPUT:
        - real -- K_ir(x)*exp(pi*r/2)

    EXAMPLES::


        sage: my_kbes_diff_r(9.45,0.861695276766 ,mpmath.fp)
        -0.31374673969963851
        sage: my_kbes_diff_r(9.4,0.861695276766 ,mpmath.fp)
        0.074219541623676832


    """
    h=mp_ctx.mpf(1e-8)
    f1 = my_kbes(r,x+h,mp_ctx)
    f2 = my_kbes(r,x,mp_ctx)
    diff=(f2-f1)/h
    return diff

### If we need to figure out the format of input, i.e. how many levels of dictionary we have
## Ex: C is input coefficients and:
##     C[0][0][0] = 0
##     C[0][0][1] = 1

def dict_depth(d,i=0):
    if isinstance(d,dict):
        if d.has_key(0):
            return dict_depth(d[0],i+1)
    return i
            
def scattering_determinant_Hecke_triangle(s,q,prec=0,use_eisenstein=0,**kwds):
    r"""
    Computes the scattering determinant phi(s)
     for the Hecke triangle group G_q


    INPUT:
     - ``s``  -- complex
     - ``prec`` -- precision used if s does not have precision
    
    """

    if q not in [3,4,6] and use_eisenstein==0:
        ## Otherwise we can only use Eisenstein series
        raise NotImplementedError
    elif use_eisenstein == 0:
        z = scattering_determinant_sl2z(s,prec=prec)
        if q==3:
            return z
        if q==4:
            llambda = s.parent()(2).sqrt()
        elif q==6:
            llambda = s.parent()(3).sqrt()
        f1 = llambda**(1-2*s)
        llambda = llambda.log()
        f = f1*( (1-s)*llambda).cosh()/( s*llambda).cosh()
        return f*z
    else:
        G = HeckeTriangleGroup(q)
        M = MaassWaveForms(G)
        if prec>0:
            s = ComplexField(prec)(s)
        return M.scattering_determinant(s,**kwds)


def scattering_determinant_sl2z(s,prec=0,verbose=0):
    r"""
    Computes the scattering determinant :
      phi(s)=sqrt(pi)gamma(s-1/2)zeta(2s-1)/gamma(s)/zeta(2s)
    for PSL2(Z).

    INPUT:
     - ``s``  -- complex
     - ``prec`` -- precision used if s does not have precision
    
    """
    if prec<=0:
        prec = 53 
    if hasattr(s,"prec"):
        if prec<s.prec():
            prec = s.prec()
    if verbose>0:
        print "prec=",prec
    RF=RealField(prec)
    CF = ComplexField(prec)
    s = CF(s.real(),s.imag())
    sqpi=RF.pi().sqrt()
    mp1=RF(1); mp2=RF(2); mp05=RF(1)/RF(2)
    res = sqpi*(s-mp05).gamma()*(mp2*s-mp1).zeta()
    res = res/s.gamma()/(mp2*s).zeta()
    return res

def eisenstein_series_coefficient_sl2z(s,m,prec=0):
    r"""
    Computes the Fourier coefficients of the Eisenstein series E(s;z) for PSL(2,Z) using the explicit formula.

    INPUT:
     - ``s``  -- complex
     - ``m``  -- integer
     - ``prec`` -- precision used if s does not have precision

    """
    if hasattr(s,"prec"):
        prec = s.prec()
    elif prec>0:
        prec = prec
    else:
        prec = 53
    RF=RealField(prec)
    CF = ComplexField(prec)
    s = CF(s)
    mppi=RF.pi()
    mp1=RF(1); mp2=RF(2); mp05=RF(1)/RF(2)
    res = mp2*mppi**s*abs(m)**(s-mp05)
    res = res/s.gamma()/(mp2*s).zeta()
    summa=CF(0)
    for d in divisors(m):
        summa+=RF(d)**(mp1-2*s)
    res = res * summa
    return res
