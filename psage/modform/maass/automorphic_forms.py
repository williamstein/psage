# -*- coding: utf-8 -*-
#*****************************************************************************
#  Copyright (C) 2010 Fredrik Strömberg <fredrik314@gmail.com>
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
Implements Spaces of automorphic forms, for example Harmonic Weak Maass forms.

AUTHORS:

- Fredrik Strömberg


EXAMPLES::


    # Note: need mpmath and mpc installed for multiprecision computations
    #       install mpmath using:
    #       sage -i mpmath
    #       sage -i mpc

    
    sage: H=HarmonicWeakMaassFormSpace(Gamma0(8),3/2)
    sage: PP=[{'+': {(1,0):0,(3,0):0}, '-': { (0, 0): 1}}]
    sage: setc={(0,-1):0,(0,-2):0}
    sage: F=H.get_element(PP,SetM=15,SetC=setc)
    sage: F.list_coefficients(5,norm=True)  

    # Half integral weight forms
    # Classical Example of weight 3/2 Harmonic Maass form (Eisenstein series) on Gamma0(4).
    # See e.g. Hirzebruch-Zagier...
    # First construct a space of modular forms

    sage: M=HalfIntegralWeightForms(Gamma0(4),1/2);M
    Space of Modular Forms of weight = 1/2  with theta multiplier on Congruence Subgroup Gamma0(4)
    # if we have magma installed we can compute the basis quickly...
    sage: M.basis()            
    [1 + 2*q + 2*q^4 + 2*q^9 + O(q^12)]
    # else we use the generic numerical algorithm...

    # Compute F with \xi_k(F)=M.basis()[0]
    sage: F=M.xi_k_inverse_basis()
    sage: M1=HalfIntegralWeightForms(Gamma0(4),3/2);M1
    Space of Modular Forms of weight = 3/2  with theta multiplier on Congruence Subgroup Gamma0(4)
    sage: M1.basis()                                  
    [1 + 6*q + 12*q^2 + 8*q^3 + 6*q^4 + 24*q^5 + 24*q^6 + 12*q^8 + 30*q^9 + 24*q^10 + 24*q^11 + O(q^12)]
    ## the interesting such F is obtained by subtracting off the theta series above
    ## the following function try to do some educated guesses in this direction... 
    sage: F=M.xi_k_inverse_basis() 


    # normalize to get the class number H(n) as c^+(n)


    sage: F1=F*mpmath.mp.mpf(2)/(mpmath.mp.mpf(16)*mpmath.mp.pi())
    sage: F1.list_coefficients(10,norm=False,cusp=0)
    For normalization we use:
    c0=c[ 0 0 ]= (-0.08333333333333333479249476 + 1.374218596157417238432065e-17j)
    c0**4= (0.0000482253086419753120196638 - 3.18106156517920673754516e-20j)
    c0**-4= (20735.99999999999854765577 + 1.367798246876169665343804e-11j)
    c-[ 0 , -1 ]= (0.1410473958869390714828967 + 1.040527440338284773385003e-18j)
    C[ 0 , -10 : -10.0 ]= (-6.048848586346993727605162e-13 + 4.124029470887434151519853e-13j)
    C[ 0 , -9 : -9.0 ]= (0.4231421876609869221276164 - 1.975630542615197649574086e-13j)
    C[ 0 , -8 : -8.0 ]= (3.831516475716355575569604e-14 - 2.458561245383529927285641e-14j)
    C[ 0 , -7 : -7.0 ]= (-9.500735585919241090009618e-15 + 1.110921884248782522963148e-14j)
    C[ 0 , -6 : -6.0 ]= (-2.384259240288897759344736e-15 + 1.403583988174435702911771e-15j)
    C[ 0 , -5 : -5.0 ]= (4.815538765737730623603249e-16 - 5.761338517547018931126834e-16j)
    C[ 0 , -4 : -4.0 ]= (0.2820947917738782846742579 - 7.217635823264244290311293e-17j)
    C[ 0 , -3 : -3.0 ]= 0.0
    C[ 0 , -2 : -2.0 ]= 0.0
    C[ 0 , -1 : -1.0 ]= (0.1410473958869390714828967+ 1.040527440338284773385003e-18j)
    C[ 0 ,- 0 ]= (0.03978873577297383394222094 + 0.0j)
    C[ 0 ,+ 0 ]= (-0.08333333333333333479249476 + 1.374218596157417238432065e-17j)
    C[ 0 , 1 : 1.0 ]= 0.0
    C[ 0 , 2 : 2.0 ]= (-1.798760787601009629745205e-17 + 1.65022799150830524765262e-16j)
    C[ 0 , 3 : 3.0 ]= (0.3333333333333333196776257 + 1.144304356995054422679724e-16j)
    C[ 0 , 4 : 4.0 ]= (0.4999999999999999916474501 + 8.162174201008810051009322e-17j)
    C[ 0 , 5 : 5.0 ]= (-1.764214819638824774085321e-17 + 2.964560603210838161433913e-16j)
    C[ 0 , 6 : 6.0 ]= (-2.968188127647462036121825e-17 + 3.285957770389621651491347e-16j)
    C[ 0 , 7 : 7.0 ]= (0.999999999999999814307206 + 2.96054723015787697893505e-16j)
    C[ 0 , 8 : 8.0 ]= (0.9999999999999998338067601 + 2.976977203293335589862541e-16j)
    C[ 0 , 9 : 9.0 ]= (1.993722392414312843107319e-15 - 2.548391116259674212745512e-15j)
    C[ 0 , 10 : 10.0 ]= (2.24017436548042636561142e-15 - 1.838191079695203924831949e-15j)




    sage: M=HalfIntegralWeightForms(Gamma0(8),3/2);M
    Space of Modular Forms of weight = 3/2  with theta multiplier on Congruence Subgroup Gamma0(8)
    sage: H=HarmonicWeakMaassFormSpace(M);H             
    Space of Harmonic Weak Maass Forms of weight = 1/2  with theta multiplier on Congruence Subgroup Gamma0(8)
    # magma dependent method
    sage: M.basis() 
    [1 + 8*q^3 + 6*q^4 + 12*q^8 + 24*q^11 + O(q^12), q + 2*q^2 + 4*q^5 + 4*q^6 + 5*q^9 + 4*q^10 + O(q^12)]
    # my numerical method
    sage: B=M.basis_numerical();B;[G1,G2]=B
    [Element of Space of Modular Forms of weight = 3/2  with theta multiplier on Congruence Subgroup Gamma0(8), Element of Space of Modular Forms of weight = 3/2  with theta multiplier on Congruence Subgroup Gamma0(8)]
    sage: G1.list_coefficients(5,norm=True,cusp=0)
    For normalization we use:
    c0=c[ 0 0 ]= 1
    c0**4= 1.0
    c0**-4= 1.0
    C[ 0 ,- 0 ]= 0
    C[ 0 ,+ 0 ]= 1
    C[ 0 , 1 : 1.0 ]= 0
    C[ 0 , 2 : 2.0 ]= (3.33584254329255000726713e-16 + 3.828443554530340392128911e-16j)
    C[ 0 , 3 : 3.0 ]= (8.000000000000000026678265 + 2.7994627907636757613313e-17j)
    C[ 0 , 4 : 4.0 ]= (6.000000000000000038839083 + 4.677622975986922840512461e-17j)
    C[ 0 , 5 : 5.0 ]= (6.514117343061377568982691e-16 + 7.582781737715197768359657e-16j)
    sage: G1.list_coefficients(5,norm=True,cusp=1)
    For normalization we use:
    c0=c[ 1 0 ]= (-0.2973017787506802628334689 - 0.2973017787506803222047077j)
    c0**4= (-0.03125000000000001086426115 - 1.24812654886380782517272e-17j)
    c0**-4= (-31.99999999999998887499659 + 1.27808158603653832431049e-14j)
    C[ 1 ,- 0 ]= 0
    C[ 1 ,+ 0 ]= (-0.2973017787506802628334689 - 0.2973017787506803222047077j)
    C[ 1 , 1 : 1.0 ]= (11.99999999999999865794267 - 1.548857790026970449743772e-15j)
    C[ 1 , 2 : 2.0 ]= (5.999999999999999962775556 - 1.804866240517242240183602e-17j)
    C[ 1 , 3 : 3.0 ]= (23.99999999999999732633269 - 3.12728962174064037529157e-15j)
    C[ 1 , 4 : 4.0 ]= (11.99999999999999973191219 - 2.190742938951590285367188e-16j)
    C[ 1 , 5 : 5.0 ]= (23.99999999999999720631503 - 3.29049103498040159747239e-15j)



    

"""

#*****************************************************************************
#  Copyright (C) 2010 Fredrik Strömberg <stroemberg@mathematik.tu-darmstadt.de>,
#
#  Distributed under the terms of the GNU General Public Licensetheta (GPL)
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

import mpmath
from sage.all import SageObject,Parent,ln,latex,random,divisors,ModularForms,prime_divisors,real,imag,PowerSeriesRing,PolynomialRing,CyclotomicField,dimension_cusp_forms,dimension_modular_forms,CuspForms
from mpmath import mpf
from psage.modform.arithgroup.mysubgroup import *
from automorphic_forms_alg import *
from sage.all import I,dumps,loads,ComplexField,LaurentPolynomialRing


from multiplier_systems import *
from psage.matrix.matrix_complex_dense import *
from psage.modform.arithgroup.all import MySubgroup,MySubgroup_class
from sage.all import magma

from vv_harmonic_weak_maass_forms_alg import vv_harmonic_wmwf_setupV_mpc2,vv_holomorphic_setupV_mpc


class AutomorphicFormSpace(Parent):
    r"""
    General class of automorphic forms.
    Subclasses shoul
    d specialize to various types. 
    """
    def __init__(self,G,weight=0,multiplier="",character=0,holomorphic=False,weak=True,almost_holomorphic=False,cuspidal=False,unitary_action=0,dprec=15,prec=53,verbose=0,**kwds):
        r""" Initialize the space of automorphic forms.
        """
        self._from_group = None # try to keep the group used to construct the MyGroup instance
        if isinstance(G,(MySubgroup_class,HeckeTriangleGroup)):
            self._group=G
        elif is_int(G):
            self._from_group = Gamma0(G)
            self._group=MySubgroup(self._from_group)
        elif str(type(G)).find("gamma")>0 or str(type(G)).find("SL2Z")>0:
            self._from_group = G
            try:
                self._group=MySubgroup(G)
            except TypeError:
                raise TypeError,"Incorrect input!! Need subgroup of PSL2Z! Got :%s" %(G)
        else:
            raise TypeError,"Could not convert G:{0} to a group!".format(G)
        self._unitary_action=unitary_action
        self._sym_type=None
        ## Define the character
        # And then the multiplier system, which should be an instance of MultiplierSystem or Subclass thereof, or None.
        if isinstance(character,sage.modular.dirichlet.DirichletCharacter):
            self._character = character
        elif is_int(character) and self._group.is_Gamma0():
            DG = DirichletGroup(self.level())
            if character >0 and character <len(DG): 
                self._character = DG[character]
            else:
                if verbose>0:
                    print "got character={0} as input!".format(character)
                self._character = trivial_character(1)
        elif character==0:
            self._character = trivial_character(1)
        else:
            raise TypeError,"Could not find character {0} on group {1}".format(character,self._group)

        if not multiplier:
            self._multiplier = TrivialMultiplier(self._group,character=self._character)
        elif isinstance(multiplier,MultiplierSystem):
            self._multiplier=multiplier
            self._character = multiplier._character
        else:
            raise TypeError,"Incorrect multiplier! Got: %s" %multiplier
        self._rdim=self._multiplier._dim
        # We assume weights are given as rational (integer of  half-integers)
        try:
            self._weight=QQ(weight)
        except:
            raise TypeError," Need weights as rational numbers! Got:%s" % weight
        # Check consistency of multiplier
        if not self._multiplier.is_consistent(self._weight):
            #print "mul=",self._multiplier
            #print "wt=",self._weight,type(self._weight)
            #print "even=",self._multiplier._character.is_even()
            #print "test=",self._multiplier.is_consistent(self._weight)
            #return self._multiplier
            raise ValueError," The specified multiplier is not compatible with the given weight! \n multiplier:{0}, weight:{1}".format(self._multiplier,self._weight)

        self._dprec=dprec
        self._prec=prec
        if dprec>15:
            self._mp_ctx=mpmath.mp
            if prec==53:
                self._prec=int(3.4*dprec)+1
        else:
            self._mp_ctx=mpmath.fp
        #self._weight=mpmath.mpf(weight)
        self._holomorphic=holomorphic
        self._almost_holomorphic=almost_holomorphic
        self._weak=weak
        self._verbose=verbose
        self._cuspidal=cuspidal
        self._alphas={}
        self._dimension = -1
        # if we are interested in dimension of the corresponding cusp form space
        self._dimension_cusp_forms = -1
        self._dimension_modular_forms = -1
        self._basis_numerical = None
        # A properly working typing would make this unecessary
        self._is_automorphic_form_space=True
        self._scaled=False
        if(self._multiplier.is_trivial()):
            self._rep=False
        else:
            self._rep=True
        # testing for types are tedious when in the interactive setting... 
        self._is_space_of_automorphic_functions=True
        self._do_mpmath = kwds.get("do_mpmath",0)
        
    def __repr__(self):
        r"""
        Return the string representation of self.

        EXAMPLES::


          sage: S=AutomorphicFormSpace(Gamma0(4));S
          Space of Automorphic Forms of weight = 0  with trivial multiplier  and character: Dirichlet character modulo 4 of conductor 1 mapping 3 |--> 1 on the group G:
          Arithmetic Subgroup of PSL2(Z) with index 6. Given by: 
          perm(S)=(1,2)(3,4)(5,6)
          perm(ST)=(1,3,2)(4,5,6)
          Constructed from G=Congruence Subgroup Gamma0(4)
          sage: S=AutomorphicFormSpace(Gamma0(4),multiplier=theta_multiplier,weight=1/2);S
          Space of Automorphic Forms of weight = 1/2  with theta multiplier  on the group G:
          Arithmetic Subgroup of PSL2(Z) with index 6. Given by: 
          perm(S)=(1,2)(3,4)(5,6)
          perm(ST)=(1,3,2)(4,5,6)
          Constructed from G=Congruence Subgroup Gamma0(4)
        
        """
        s="Space of "
        if self.is_holomorphic():
            if self.is_cuspidal():
                s+="Cusp Forms "
            else:
                s+="Modular Forms "
        elif str(type(self)).find("HarmonicWeakMaassFormSpace")>0:
            s+="Harmonic Weak Maass Forms "
        else:
            s+="Automorphic Forms "
        s+="of weight = "+str(self._weight)+" "
        if str(self._multiplier).find("theta_multiplier")>0:
            s+=" with theta multiplier "
        elif not self._multiplier.is_trivial():
            s+=" with multiplier:\n"+str(self._multiplier)
        else:
            s+=" with trivial multiplier "
        #if(self._character<>trivial and self._character<>trivial_character(self._group.level())):            
        #    s+=" and character: "+str(self._character)+" "
        s+=" on "
        if(self._from_group):
            s+=str(self._from_group)
        else:
            if self.group().index()==1:
                s+=" SL(2,Z)"
            else:
                s+="the group G\n"+str(self._group)
        return s


    def __reduce__(self):
        r"""
        """
        # we can't pickle functions so we store the names instead
        #if(isinstance(self._multiplier,type(trivial))):
        #    multiplier_s=self._multiplier.func_name            
        #else:
        #    multiplier_s=self._multiplier
        #if(isinstance(self._character,type(trivial))):
        #    character_s=self._character.func_name            
        #else:
        #    character_s=self._character
        if(self._from_group):
            G = self._from_group
        else:
            G = self._group
        return(AutomorphicFormSpace,(G,self._weight,self.multiplier(),self._holomorphic,self._weak,self._almost_holomorphic,self._cuspidal,self._dprec,self._verbose))

    def __eq__(self,other):
        r"""
        Compare self to other.
        """
        if self._verbose>0:
            print "in AutomorphicFormSpace.__eq__"
        if(not isinstance(other,type(self))):
            return False
        if(self._weight <> other._weight):
            return False
        if(self._group <> other._group):
            return False
        if(self._multiplier <> other._multiplier):
            return False
        if(self._character <> other._character):
            return False
        if(self._holomorphic <> other._holomorphic):
            return False
        if(self._weak <> other._weak):
            return False
        if(self._cuspidal <> other._cuspidal):
            return False
        #eq = eq and (self._dprec == other._weak)
        #    return False
        #print "eq=",eq
        return True

    
    def __ne__(self,other):
        r"""
        Compare self to other.
        """
        return not self.__eq__(other)
    # return various properties of self

    def group(self):
        r""" Return the group of self.
        """
        return self._group

    def weight(self):
        r""" Return the weight of self.
        """
        return self._weight

    def character(self):
        r""" Return the character of self.
        """
        return self._multiplier._character

    def multiplier(self):
        r""" Return the multiplier of self.
        """
        return self._multiplier

    def sym_type(self):
        if self._sym_type==None:
            if hasattr(self._multiplier,"_sym_type"):
                self._sym_type=self._multiplier._sym_type
            else:
                self._sym_type=0
        return self._sym_type

    def prec(self,prec=None):
        if prec<>None:
            self._prec=prec
        return self._prec

    def dprec(self,dprec=None):
        if dprec<>None:
            self._dprec=dprec
        return self._dprec

    def is_holomorphic(self):
        r"""
        Return True if self is holomorphic, otherwise False. 
        """
        return self._holomorphic

    def is_almost_holomorphic(self):
        r"""
        Return True if self is almost holomorphic, otherwise False. 
        """
        return self._almost_holomorphic            

    def is_cuspidal(self):
        r"""
        Return True if self is cuspidal, otherwise False. 
        """
        return self._cuspidal



    def is_weak(self):
        r"""
        Return True if self is a weak form, i.e. has pole(s) at oo, otherwise False. 
        """
        return self._weak

    def is_harmonic(self):
        r"""
            Return True if self is a harmonic weak maassform, i.e. has pole(s) at oo, otherwise False. 
        """
        return self._harmonic


    def level(self):
        r""" Return the level of self (if self.group is a congruence subgroup).
        """
        if not self._group.is_congruence():
            raise ValueError,"Level is only defined for congruence subgroups!"
        return self._group.generalised_level()

        
    def rep(self,A):
        r"""
        Calculate the representation = multiplier (including character) of self on the matrix A
        """
        v= self._multiplier(A)
        # duality and character should be builtin in the multiplier
        #    v = 1/v
        #v=v*self._character(A[1,1])
        return v

    def alpha(self,i):
        r"""
        Compute the translation at cusp nr. i, i.e. v(T_i)=e(alpha(i))
        -''i'' -- index of the cusp
        
        """
        RF=RealField(self._prec)
        CF=ComplexField(self._prec)
        if(not self._alphas.has_key(i)):
            if(self._multiplier == None or self._multiplier.is_trivial()):
                self._alphas[i]=[RF(0),CF(1)]
            elif self.multiplier().ambient_rank()==1:
                #p=self._group._cusps[i]
                A=self._group._cusp_data[i]['stabilizer']
                tmp = self.rep(A)
                if hasattr(tmp,"complex_embedding"):
                    v=tmp.complex_embedding(self._prec)
                else:
                    v=CF(tmp)
                #ar=mpmath.arg(v)/mpmath.pi()/mpmath.mpf(2)
                ar=v.argument()/RF.pi()/RF(2)
                self._alphas[i]=[ar,v]
            else: ## This is nw a matrix-valued representation
                A=self._group._cusp_data[i]['stabilizer']
                mat = self.rep(A)[0]
                tmp = []
                for n in range(mat.nrows()):
                    a = mat[n,n]
                    if hasattr(a,"complex_embedding"):
                        v=a.complex_embedding(self._prec)
                    else:
                        v=CF(tmp)
                    ar=v.argument()/RF.pi()/RF(2)
                    ## Find the exact alpha
                    l = a.multiplicative_order()
                    if l < Infinity:
                        z = CyclotomicField(l).gens()[0]
                        for k in range(l):
                            if z**k == v:
                                break
                        tmp.append((ar,v,k,l))
                    else:
                        tmp.append((ar,v))
                self._alphas[i]=tmp
        return self._alphas[i]

    def alphas(self):
        return self._alphas

    def set_alphas(self,prec=None):
        r""" Compute the vector containing the shifts at the various cusps.
        """
        precold=self._prec
        if prec<>None and prec<>self._prec:
            self._prec=prec
        for i in range( self._group._ncusps):
            self.alpha(i)
        self._prec=precold
    def dimension(self):
        r"""
        Return the dimension of self if we can compute it.
        """
        if self._dimension>=0:
            return self._dimension
        k=self._weight
        if self._holomorphic and self._rdim==1:
            if is_int(k) and self._multiplier.is_trivial():
                # Use the builtin sage routines
                xi = self.character()
                self._dimension_cusp_forms =  dimension_cusp_forms(xi,k)
                self._dimension_mod_forms =  dimension_modular_forms(xi,k)
            else:
                # In weight 1/2 use Serre-Stark
                if k==QQ(1)/QQ(2):
                    [d_mod,d_cusp]=self._dimension_of_weight_one_half_space()
                    if(self.is_cuspidal()):
                        return d_cusp
                    else:
                        return d_mod                
                elif k<QQ(1)/QQ(2):
                    self._dimension = 0
                    # Else use Cohen-Oesterle: (as described in The Web of Modularity by Ono)
                elif k>QQ(3)/QQ(2):
                    dimension_cusp_forms = self._difference_of_dimensions(k)
                    dimension_mod_forms = -self._difference_of_dimensions(QQ(2-k))
                elif k==QQ(3)/QQ(2):
                    [d_mod,d_cusp]=self._dimension_of_weight_one_half_space()                
                    dimension_cusp_forms = self._difference_of_dimensions(k)+d_mod
                    dimension_mod_forms = d_cusp - self._difference_of_dimensions(QQ(1)/QQ(2))
                # print "dim_cusp_forms=",dimension_cusp_forms
                # print "dim_mod_forms=",dimension_mod_forms            
                self._dimension_cusp_forms = dimension_cusp_forms
                self._dimension_mod_forms = dimension_mod_forms
                if(self.is_cuspidal()):
                    self._dimension = self._dimension_cusp_forms
                else:
                    self._dimension = self._dimension_mod_forms
        elif self._holomorphic:
            #if self._group.level()==1:
            #[d_mod,d_cusp]=self._dimension_of_vector_valued_forms()
            self._dimension = self._dimension_of_vector_valued_forms()
            
        else:
            raise NotImplementedError
                #self._dimension=-1
            
        return self._dimension

    def _difference_of_dimensions(self,k):
        r"""
        Use Cohen and Oesterle to compute the difference: dim(S_k)-dim(M_{2-k})
        """
        N = self._group.generalised_level()
        cond = self._character.conductor()
        #k = QQ(RR(self._weight)) 
        kk = ZZ(k - QQ(1)/QQ(2))
        r2 = valuation(N,2)
        s2 = valuation(cond,2)
        if(r2>=4): #   zeta_k_l_chi = lambda_k_l_chi
            if 2*s2 <= r2:
                if is_even(r2):
                    rp = r2/QQ(2)
                    zeta_k_l_chi = 2**(rp) + 2**(rp-1)
                else:
                    rp = (r2-1)/QQ(2)
                    zeta_k_l_chi = 2*2**(rp)
            elif 2*s2 > r2:
                zeta_k_l_chi = 2*2**(rp-sp)            
        elif(r2==3):
            zeta_k_l_chi = 3
        elif(r2==2):
            zeta_k_l_chi = 0
            ## Condition (C)
            for p in prime_divisors(N):
                if( (p % 4) == 3):
                    rp = valuation(N,p)
                    sp = valuation(cond,p)
                    if(is_odd(rp) or (rp>0 and rp < 2*sp)):
                        zeta_k_l_chi = 2
                        break
            if zeta_k_l_chi== 0: # not (C)
                if(is_even(kk)):
                    if(s2==0):
                        zeta_k_l_chi = QQ(3)/QQ(2)
                    elif(s2==2):
                        zeta_k_l_chi = QQ(5)/QQ(2)
                else:
                    if(s2==0):
                        zeta_k_l_chi = QQ(5)/QQ(2)
                    elif(s2==2):
                        zeta_k_l_chi = QQ(3)/QQ(2)
        if(zeta_k_l_chi<=0):
            raise ArithmeticError,"Could not compute zeta(k,l,chi)!"
        fak = QQ(1)
        for p in prime_divisors(N):
            fak = fak* QQ(1+QQ(1)/QQ(p))
        fak2 = QQ(1)
        for p in prime_divisors(N):
            if(p>2):
                rp = valuation(N,p)
                sp = valuation(cond,p)
                if(rp < 2*sp):
                    lam = QQ(2*p**(rp-sp))
                else:
                    if(is_even(rp)):
                        rprim=QQ(rp)/QQ(2)
                        lam = QQ(p**rprim)+QQ(p**(rprim-1))
                    else:
                        rprim=QQ(rp-1)/QQ(2)
                        lam = QQ(2*p**rprim)
                fak2 = fak2 * lam
            # S_k(Gamma0(N),chi) - M_{2-k}(Gamma0(N),chi)
        diff_dims = fak*QQ(k-1)*QQ(N)/QQ(12)-QQ(zeta_k_l_chi)/QQ(2)*fak2
        return diff_dims

    def _dimension_of_weight_one_half_space(self):
        r"""
        Computes the dimension of M and S_{1/2}(4N,chi) where 4N and chi are the level and character of self.
        
        """
        O = self._Omega()
        dim_mod_forms = len(O)
        # check number of totally even characters for the cusp forms
        nn = 0
        for (xi,t) in O:
            l = xi.decomposition()
            for xip in l:
                if(xip(-1)==1):
                    nn = nn+1
        dim_cusp_forms = len(O) - nn
        return [dim_mod_forms,dim_cusp_forms]

    def _Omega(self):
        r"""
        Computes the set of pairs (psi,t) satisfying:
        (1) r**2 * t | self.level()/4
        (2) chi(n)= psi(n)*kronecker(t,n) for (n,self.level())=1

        """
        N = ZZ( QQ(self.level())/QQ(4))
        D = DirichletGroup(self.level())
        chi = self._character
        Omega=[]
        for t in divisors(N):
            for psi in D:
                r = psi.conductor()
                s = ZZ(r*r*t )
                if(not s.divides(N)):
                    continue
                ok = True
                for n in range(1,self.level()):
                    if(psi(n)*kronecker(t,n)<>chi(n)):
                        ok = False
                        break
                if(ok):
                    Omega.append((psi,t))
        return Omega

    def _dimension_of_vector_valued_forms(self):
        r"""
        Calculates the dimension of self is self is a space of automorphic forms  on SL2(Z).
        """
        if self._dimension>=0:
            return self._dimension
        if self._weight < 2 or self.level()<>1:
            return -1
        try:
            if hasattr(self._multiplier,"dimension_cusp_forms"):
                if self._cuspidal:
                    self._dimension = self._multiplier.dimension_cusp_forms(self._weight)
                else:
                    self._dimension = self._multiplier.dimension_modular_forms(self._weight)
                return self._dimension
        except:
            pass
        term0 = QQ(self._rdim*(self._weight-1))/QQ(12)
        S,T=SL2Z.gens()
        R = S*T; R2=R*R
        if self._rdim>1:
            wS=self._multiplier(S)[0].trace()
            wR=self._multiplier(R)[0].trace()        
            wR2=self._multiplier(R2)[0].trace()
            evs = self._multiplier(T)[0].diagonal()
            wT=sum(evs) #self._multiplier(T).trace()
            alphas = map(lambda x:log(CC(x))/CC(2/pi), evs)
        else:
            wS=self._multiplier(S)
            wR=self._multiplier(R)        
            wR2=self._multiplier(R2)            
            wT=self._multiplier(T)
        z2=CyclotomicField(4).gen()
        term1 = wS/QQ(4)*z2*z2**(self._weight-1)
        z3=CyclotomicField(6).gen()        
        term2 = wR/QQ(3)/sqrt(3)*z2*z3**(self._weight-1)
        term3 = wR/QQ(3)/sqrt(3)*z2*z3**(2*(self._weight-1))
        term4=0
        k0=0
        for a in alphas:
            if a<>0:
                term4+=a-0.5
            else:
                k0+=1
        term5 = k0/2*sign(self._weight-1)
        if self._cuspidal:
            term6 = 0
        else:
            term6 = k0
        if k0<>0:
            if self._weight==1:
                raise ArithmeticError,"Need to compute the scattering determinant!"

        dim = term0 + term1 + term2 + term3 + term4 + term5 + term6
        return dim
        
    ## cuspidal subspace
    def cuspidal_subspace(self):
        r"""
        Construct the cuspidal subspace of self.        
        """
        S = copy(self) # copy self
        S._weak = False # not weak
        S._cuspidal = True # and cuspidal
        #S._holmorphic = True # and holmorphic in H
        #S=HalfIntegralWeightForms(G,self._weight,self._multiplier,character=self._character,holomorphic=self._holomorphic,weak=self._weak,cuspidal=True,dprec=self._dprec,verbose=self._verbose,construct=self._construct)
        return S



    def set_normalization(self,C=None):
        r"""
        -''C'' -- a dictionary of set coefficients in the form C[d][(i,n)]=c
        """
        #print "C0=",C
        N=dict()
        N['comp_dim']=1
        if isinstance(C,dict):
            N['comp_dim']=max(1,len(C.keys()))
        else:
            N['comp_dim']=max(1,len(C))
        N['SetCs']=dict()
        N['cuspidal']=self._cuspidal
        N['weak']=self._weak        
        nc=self._group.ncusps()
        for j in range(N['comp_dim']):
            N['SetCs'][j]=dict()
        #if(P.has_key((0,0)) and H._holo):
        #print "holo"
        #for j in range(N['comp_dim']):
        #    N['SetCs'][j][(0,0)]=0            
        if(N['cuspidal']):
            for icusp in range(nc):
                v=self.alpha(icusp)[1]
                if(v==1):
                    for j in range(N['comp_dim']):
                        N['SetCs'][j][(icusp,0)]=0
        if(not N['weak']):
            for icusp in range(nc):
                al=self.alpha(icusp)[0]
                if(al<-mpmath.eps()):
                    for j in range(N['comp_dim']):
                        N['SetCs'][j][(icusp,0)]=0
        if isinstance(C,dict) and C<>{}:
            for i in C.keys():
                for (r,n) in C[i].keys():
                    N['SetCs'][i][(r,n)]=C[i][(r,n)]
        elif isinstance(C,list) and C<>[]:
            for i in range(len(C)):
                for (r,n) in C[i].keys():
                    N['SetCs'][i][(r,n)]=C[i][(r,n)]
        return N

    def set_normalization_vv(self,P={},C={},c_t="pp"):
        r"""
        Set normalization for vector-valued case
        """
        N=dict()
        if isinstance(P,list):
            Pl=P
        else:
            Pl=[P]
        if isinstance(C,list):
            Cl=C
        else:
            Cl=[C]        
        if len(Pl)>0:
            N['comp_dim']=len(Pl)
        elif len(Cl)>0:
            N['comp_dim']=len(Cl)
        else:
            raise ValueError,"Need either principal parts of set coefficients!"
        if len(Cl)>0:
            if len(Cl)<>len(Pl):
                raise ValueError,"Need same number of principal parts and coefficients to set!"
            keys = Cl[0].keys()
            for j in range(1,N['comp_dim']):
                if Cl[j].keys()<>keys:
                    raise ValueError,"Need to set the same coefficients! (or call the method more than once)"
        else:
            Cl=[]
            for j in range(N['comp_dim']):
                Cl.append(C)
        if self._verbose>0:
            print "Pl=",Pl
            print "Cl=",Cl
        N['Vals']=list()
        N['Vals']=list()
        N['SetCs']=list()
        for i in range(N['comp_dim']):
            N['Vals'].append({})
            N['Vals'].append({})
            N['SetCs'].append([])
            ## First look at all zeroth-coefficients and see if we are forced to set some to zero
            for j in range(self.multiplier().weil_module().rank()):
                a=self.multiplier().weil_module().basis()[j]
                x=self.multiplier().Qv[j]
                #N['Vals'][i][(0,j)]=dict()
                
                if x==0:
                    if c_t=="pp" and Pl[i].has_key((0,j)):
                        N['SetCs'][i].append((j,0))
                        N['Vals'][i][(j,0)]=Pl[i][(j,0)]
                    elif self._cuspidal:
                        N['SetCs'][i].append((j,0))
                        N['Vals'][i][(j,0)]=0 #P[(0,0)]
                    
                elif x<0 and self._holomorphic:
                    N['SetCs'][i].append((j,0))
                    N['Vals'][i][(j,0)]=0 #P[(0,0)]            

            if isinstance(Cl[i],dict):
                for (r,n) in Cl[i].keys():
                    if(N['SetCs'][i].count((r,n))==0):
                        N['SetCs'][i].append((r,n))
                        N['Vals'][i][(r,n)]=Cl[i][(r,n)] 
        return N



    def get_Y_and_M(self,digs=10,principal_part=[{}]):
        r"""
        Get good values for Y and M to truncate with an error of prec digits
        """
        ## todo : more alternatives
        ## we use the largest term of the principal part
        ## to estimate the Y and M
        #print "pp0=",principal_part
        if not isinstance(principal_part,list):
            pp = [principal_part['+']]
        elif len(principal_part)==1:
            if principal_part[0].has_key('+'):
                pp = [principal_part[0]['+']]
            else:
                pp = [principal_part[0]]
        else:
            pp=list()
            for p in principal_part:
                pp.append(p['+'])
        maxn=0; maxa=1; maxr=0
        #print "pp1=",pp
        for P in pp: #rincipal_part:
            for (r,n) in P.keys():
                # Remember that the principal part (i,j):c means different things for scalar and vector-valued forms, i.e. i is the cusp in the first case and the component in the second
                a = P[(r,n)]
                if(a>maxa):
                    maxa=a
                if self.multiplier().ambient_rank()>1:
                    if r >=0 and r< self.multiplier().rank():
                        aln = n + self.alpha(0)[r][0]
                        if aln < maxn:
                            maxr = r
                            maxn = aln
                else:
                    aln = n + self.alpha(r)[0]
                    if aln < maxn:
                        maxr = r
                        maxn = aln

        if self._verbose > 1:
            print "maxa=",maxa
            print "maxn=",maxn
            print "digits=",digs
        if not self._holomorphic or self._weak:
            maxa = maxa * len(pp)
            pp_max = {(maxr,maxn):maxa}
            if self._verbose > 1:
                print "pp_max=",pp_max
            #[Y,M]=self.get_Y_and_M(prec,principal_part=pp)
            [Y,M]=get_Y_and_M_for_hwmf(self._group,pp_max,self._weight,digs)
        else:
            Y=self._group.minimal_height()*mpmath.mpf(95)/mpmath.mpf(100)
            M=get_M_for_holom(Y,self._weight,digs)
        return [Y,M]

    ## def get_element(self,principal_part={},C=None,prec=10,M0_in=None,Y_in=None):
    ##     r"""
    ##     Get an element of self given by either principal part
    ##     or normalization of coefficients.


    ##     """
    ##     print "self.type=",type(self)
    ##     F=AutomorphicFormElement(self,principal_part=principal_part)
    ##     if(Y_in<>None and M0_in<>None):
    ##         Y=Y_in
    ##         M=M0_in
    ##     elif(Y_in<>None):
    ##         Y=Y_in
    ##         if(not self._holomorphic):
    ##             M=get_M_for_hwmf(Y,self._weight,prec,principal_part)
    ##         else:
    ##             M=get_M_for_holom(Y,self._weight,prec)               
    ##     elif(M0_in<>None):
    ##         M=M0_in
    ##         Y=self._group.minimal_height()*mpmath.mpf(95)/mpmath.mpf(100)
    ##     else:
    ##         [Y,M]=self.get_Y_and_M(prec,principal_part)
    ##     #Y=Y*0.7
    ##     Q=M+30
    ##     Ymp=mpmath.mp.mpf(Y)
    ##     if(not (principal_part.has_key('+') or principal_part.has_key('-'))):
    ##         raise ValueError,"Need principal part with keys '+' and '-'!"
    ##     PP=principal_part
    ##     V=setup_matrix_for_harmonic_Maass_waveforms_sv(self,Ymp,M,Q,PP)
    ##     V['PP']=PP
    ##     print "Use M,Y=",M,Y
    ##     # recall that the zeroth coefficient is counted in the principal part
    ##     return V
    ##     if(C==None):
    ##         C=dict()
    ##     for j  in range(len(self._group.cusps())):
    ##         al=self.alpha(j)
    ##         print "al(",j,")=",al
    ##         if(al[1]==1):
    ##             if(PP.has_key((j,0))):
    ##                 for r in C.keys():
    ##                     C[r]=dict()
    ##                     C[r][(j,0)]=0
    ##     print "C=",C
    ##     N=self.set_normalization(C)
    ##     print "N=",N
    ##     D=solve_system_for_harmonic_weak_Maass_waveforms(V,N,deb=True)
    ##     F._coeffs=D
    ##     return F
    def _get_element(self,principal_part,digs=10,dbase_prec=None,SetC=None,SetY=None,SetM=None,SetQ=None,do_mpmath=0,get_mat=False,use_sym=1,get_c=False,gr=0,version=0,threads=1,**kwds):
        r"""
        INPUT:
        
        - `principal_part`   -- list of principal parts of the form:
                 RR = { '+' : {(j,n) : c^+(j,n)}     # j is a cusp and n>=0 an index
                        '-' : {(j,n) : c^-(j,n)}     # j is a cusp and n<=0 an index
                     }
                corresponding to principal parts (in notation of Bruinier-Funke):
                    \( \Sum_{n>0} c^+(j,n)q^{-n} +  \Sum_{n<0} c^-(j,n)H(n\tau)
        
                    PP[c,m]=a if the principal at cusp c contains a*q^m
        - `digs` -- integer (default 10): the number of requested digits
        - `dbase_prec` -- integer (default None): if set, use this number of digits for precision in all mpmath calculations
        - `SetC` -- dictionary containing fourier coefficients to keep fixed (and their values)
                      of the form SetC[n][i]=c_i(n)
        """
        from vv_harmonic_weak_maass_forms import solve_system_for_vv_harmonic_weak_Maass_waveforms_new

        ## the principal part and the SetC should be lists if present
        if self._verbose>0:
            print "PP=",principal_part
            print "gr=",gr
        if(not isinstance(principal_part,list)):
            ppart = [principal_part]
        else:
            ppart = principal_part
        ppart1=list()
        for pp in ppart:
            d=dict()
            d['-']=pp.get('-',{}) # By default we have no princ. part
            d['+']=pp.get('+',{}) 
            #print "pp=",pp
            if isinstance(pp.keys()[0],(list,tuple)):
                d['+']=pp # If only one is given we assume it holomorphic
            # If self._holomorphic is True and we have a negative principal part we assume
            # that the only non-holomorphic part is the principal part
            #if self._holomorphic: # Make sure no non-holom. ppart is given for a holom. form
            #    d['-']={}
            ppart1.append(d)
        if self._verbose>0:
            print "PP1=",ppart1
        ppart = ppart1 #principal_part

        ## Check whether the set coefficients are the same for all elements
        ## if thy are the same we only need to solve the system once (using LU decomposition).
        ## otherwise we need to rerun the system solving several times
        
        if SetC<>None and not isinstance(SetC,list):
            setc=list()
            for i in range(len(ppart)):
                setc.append(SetC)
        elif not SetC:
            setc = []
        else:
            setc=SetC
        if len(setc)>0 and len(setc)<>len(ppart):
            raise ValueError,"Inconsistent lengths of principal part and set coefficients!"        
        # recall that we treat 0-coefficients in the principal part
        # as variables.
        if self._verbose>0:
            print "setc0=",setc
            #print "group=",self.group()
        for i in range(len(ppart)):
            for j in range(self.group().ncusps()):
                if ppart[i]['+'].has_key((j,0)):
                    for ii in range(len(setc),i+1):
                        setc.append({})
                    setc[i][(j,0)]=ppart[i]['+'][(j,0)]
        if self._verbose>0:
            print "setc1=",setc

        solve_once = True
        if setc<>None and len(setc)>0:
            d = setc[0].keys()
            for c in setc:
                if c.keys()<>d:
                    solve_once=False
                    break
#                    raise ValueError," Inconsistent set coefficients. Need the same length in all entries! Got: %s" %(setc)
        if self._verbose>0:
            print "solve_once=",solve_once
            print "ppart=",ppart
            print "setc=",setc
        #pos_part=list()
        #for pp in ppart:
        #   pos_part.append(pp['+'])
        if(not (SetY and SetM)):
            [Y,M]=self.get_Y_and_M(digs,ppart)
        # print "dps=",mpmath.mp.dps
        #Y=Y*0.5 #97.
        if(SetY<>None):
            Y=SetY
        if(SetM<>None):
            M=SetM        
        if SetQ<>None and SetQ>M:
            Q = SetQ
        else:
            Q=M+10

        mpmath.mp.prec = self._prec
        
        Ymp=RealField(self._prec)(Y) #mpmath.mp.mpf(Y)
        dpold=mpmath.mp.dps
        #if dbase_prec<>None:
        #    mpmath.mp.dps=max(self._dprec,dbase_prec)
        #else:
        #    mpmath.mp.dps=self._dprec
        sv = 1
        d = self.multiplier().rank()

        if d>1 or  hasattr(self.multiplier(),"D"):
            sv=0
        if self._verbose>0:
            print "dps=",mpmath.mp.dps
            print "setc=",setc
            print "Y=",Ymp
            print "M=",M
            print "Q=",Q
            print "PP=",ppart
            print "do_mpmath=",do_mpmath
            print "alphas=",self.alphas()
            print "dim=",d
            print "scalar=",sv
        C = None

        if sv==1:
            if do_mpmath==1:
                Ymp = mpmath.mp.mpf(Ymp)
                V=setup_matrix_for_harmonic_Maass_waveforms_sv_bak(self,Ymp,M,Q,ppart)
            elif do_mpmath==2:
                Ymp = mpmath.mp.mpf(Ymp)
                V=setup_matrix_for_harmonic_Maass_waveforms_sv_bak_22(self,Ymp,M,Q,ppart)
            elif (version==0 or gr==1):
                V=setup_matrix_for_harmonic_Maass_waveforms(self,Ymp,M,Q,ppart,use_sym=use_sym,threads=threads)
            else:
                C = setup_and_solve_for_harmonic_Maass_waveforms(self,Ymp,M,Q,ppart,cset=setc)
        else:
            ## Only one principal part is implemented in vv-case
            if self._holomorphic:
                V=vv_holomorphic_setupV_mpc(self,Ymp,M,Q)
            else:
                V=vv_harmonic_wmwf_setupV_mpc2(self,ppart[0],Ymp,M,Q)
            V['space']=self
            V['PP']=ppart
        if gr==1:
            return V
        if solve_once and C==None:
            #N=set_norm_harmonic_weak_maass_forms(self,ppart,setc)
            #if isinstance(setc,list):
            if sv==1:                
                if do_mpmath==0:
                    N = self.set_normalization(setc)
                else:
                    N = self.set_norm(setc)
            else:
                N = self.set_normalization_vv(ppart,setc)
            #else:
            #    N = self.set_normalization(setc)
            V['PP']=ppart
            #return V,N
            if self._verbose>0:
                print "N=",N
            if do_mpmath<>0:
                C=solve_system_for_harmonic_weak_Maass_waveforms_mpmath(V,N)
            else:
                if sv==1:
                    C=solve_system_for_harmonic_weak_Maass_waveforms(V,N)
                else:
                     C=solve_system_for_vv_harmonic_weak_Maass_waveforms_new(self,V,N)

        elif C==None:
            C=list()
            RHS=V['RHS']
            if RHS.cols<>len(ppart):
                raise ValueError,"Inconsistent lengths of principal part and right hand sides!"        
            for i in range(len(ppart)):
                pp=[ppart[i]]; cc=[setc[i]]
                V['PP']=pp
                V['RHS']=RHS.column(i)
                if self._verbose>1:
                    print "cc=",cc
                    print "pp=",pp
                N = self.set_norm(ppart,setc)
                #N=set_norm_harmonic_weak_maass_forms(self,pp,cc)
                if self._verbose>1:
                    print "N=",N
                #return V,N
                try:
                    C.append(solve_system_for_harmonic_weak_Maass_waveforms(V,N)[0])
                except:
                    return C.append((V,N))
        #mpmath.mp.dps=dpold
        if self._verbose>0:
            print "C[0][-1]=",C.get(0,{}).get(0,{}).get(-1,None)

        res=list()
        if get_c:
            return C
        prec = mpmath.mp.dps
        if len(C)>0:
            for i in range(len(C)):
                #print "PP=",ppart[i]
                ppf=dict()
                ppf['+']=ppart[i]['+']
                ppf['-']=ppart[i]['-']
                if self._verbose>1:
                    print "type=",type(self)
                if str(type(self)).find("HalfIntegralWeightForms")>0:
                    F=HalfIntegralWeightFormElement(self,C[i],principal_part=ppf)
                elif str(type(self)).find("HarmonicWeakMaassFormSpace")>0:
                    if self._verbose>1:
                        print "Constructing a Harmonic Weak Maassform"
                        print "pp=",ppf
                    F=HarmonicWeakMaassFormElement(self,C[i],prec=prec,principal_part=ppf)
                else:
                    F=AutomorphicFormElement(self,C[i],prec=prec,principal_part=ppf)

                #print "M0=",M
                F._M0 = M
                res.append(F)
                #print "appended"
                #print "res=",res
        if len(res)==1:
            return res[0]
        else:
            return res
### Now define subclasses which specialize the above general space.

class HalfIntegralWeightForms(AutomorphicFormSpace):
    r"""
    Space of half-integral weight modular forms forms (with theta multiplier as default).
    """
    def __init__(self,G,weight=1/2,multiplier="theta",character=0,holomorphic=True,weak=False,cuspidal=False,dprec=25,verbose=0,construct=True,**kwds):
        r""" Initialize the space of automorphic forms.
        construct == False => do not construct the magma space (might take time)
        """
        # we can initialize a half-integral weight space
        # from an integral weight one, i.e. with Shimura correspondence
        self._shimura_image=None
        self._character = character
        self._basis_len=0
        if(hasattr(G,"group") and hasattr(G,"weight")):
            k=G.weight()
            if(is_even(k)):
                weight=QQ(k/2)+QQ(1)/QQ(2)
            else:
                raise ValueError,"Shimura Correspondence only for even weight!"
            if(not G.character().is_trivial):
                raise NotImplementedError, "We only deal with trivial character for now!"
            self._shimura_image=G
            if(hasattr(G,"is_cuspidal") and G.is_cuspidal()):
                cuspidal=True
            else:
                cuspidal=False
            # by default we go to Gamma0(4N)
            multiplier="theta"
            G=MySubgroup(Gamma0(4*G.generalised_level()))
            print "Initializing through Shimura corr!"
        if(isinstance(G,MySubgroup_class)):
            self._group=G
            self._from_group=G._G
        elif is_int(G):
            self._group=MySubgroup(Gamma0(G))
            self._from_group=Gamma0(G)
        elif( hasattr(G,'is_subgroup') and G.is_subgroup(SL2Z)):
            self._group=MySubgroup(G)
            self._from_group=G
        else:
            raise ValueError,"Did not get a group G={0}".format(G)
        if multiplier=="theta":
            if isinstance(self._character,int):
                modulus = self.level(); ch = self._character
            elif isinstance(self._character,tuple):
                modulus,ch=self._character
            multiplier=ThetaMultiplier(self._group,dchar=(modulus,ch),weight=weight)
        character=multiplier._character
        self._character = character #_set_character(character)

        #else:
        #    raise NotImplemented,"Only implemented theta multiplier for half-integral weight!"
        # fix consistency, i.e. change to conjugate if necessary:
        t1 = multiplier.is_consistent(weight)
        if not t1:
            multiplier.set_dual()
            t2 = multiplier.is_consistent(weight)
            if not t2:
                raise ValueError," Could not find consistent multiplier for multiplier: %s and weight %s!" % (multiplier,weight)
        # construct the space with the given multiplier
        AutomorphicFormSpace.__init__(self,G,weight=weight,multiplier=multiplier,holomorphic=holomorphic,weak=weak,cuspidal=cuspidal,dprec=dprec,verbose=verbose)
        
        self._magma_space=None
        # If we have Magma installed we associate Magma's space
        # unfortunately characters are only implemented for Gamma_1(N) in magma
        if construct:
            if character==None:
                try:
                    s="HalfIntegralWeightForms("+str(self.level())+","+str(QQ(weight))+")"
                    self._magma_space=magma.new(s)
                except TypeError,RuntimeError:
                    pass
            else:
                try:
                    ## Want to find the corresponding magma character
                    D1 = DirichletGroup(self.level())
                    s = "DirichletGroup("+str(self.level())+")"
                    Dm = magma.new(s)
                    magma_index_char=-1
                    # print "my vals=",character.values()
                    for j in range(1,len(Dm.Elements())+1):
                        x = Dm.Elements()[j]
                        # print "x.vals=",x.ValueList()
                        try:
                            for i in range(self.level()):
                                if(x(i)<>character(i)):
                                    raise StopIteration()
                                # if we are here we have found the correct characte
                            # print "found match!"
                            magma_index_char=j
                            break
                        except StopIteration:
                            pass
                    if(magma_index_char>=0):
                        i = magma_index_char
                        s="HalfIntegralWeightForms(Elements(DirichletGroup("
                        s=s+str(self.level())+"))["+str(i)+"],"+str(QQ(weight))+")"
                        print "S=",s 
                        self._magma_space=magma.new(s)
                    else:
                        print "Could not construct a corresponding space in Magma!" 
                except TypeError:
                    self._magma_space=None
                pass
        if isinstance(self._multiplier,ThetaMultiplier):
            w_minus_half=QQ(RR(weight))-QQ(1)/QQ(2)
            k=2*ZZ(w_minus_half)
            if(self._shimura_image==None and weight >1):
                N=ZZ(self.level()/4)
                #
                # print "wt=",weight,type(weight)
                if character==None or (character**2).is_trivial():
                    if cuspidal:                
                        self._shimura_image=CuspForms(N,k)
                    else:
                        self._shimura_image=ModularForms(N,k)
                else:
                    DG=DirichletGroup(N)
                    x=DG(character**2)
                    if cuspidal:                
                        self._shimura_image=CuspForms(x,k)
                    else:
                        self._shimura_image=ModularForms(x,k)
            
                

        self._construct=construct
        self._dimension = self.dimension()


    ## 
        

    ## Using Magma we can compute certain things explicitly
    def _magma_dimension(self):
        r"""
        If we have a holomorphic space we can use magma to compute the dimension.       """
        try:
            return self._magma_space.Dimension()
        except:
            raise NotImplementedError,"Currently we need magma for this functionality!"

    def basis(self,prec=None,method=None):
        if method == None:
            # try magma first
            try:
                res = self.q_expansion_basis(prec,method='magma')
            except ValueError:
                res = self.q_expansion_basis(prec,method='numerical')
        else:
            res = self.q_expansion_basis(prec,method=method)

        return res

        
    def q_expansion_basis(self,prec=None,method=None):
        
        r"""
        If we have a holomorphic space we can use magma to compute a basis.
        """
        if method=='magma':
            if self._magma_space == None:
                raise ValueError,"Magma not supported here. Choose different method!"
            else:
                if(prec<>None):
                    return list(self._magma_space.Basis(prec))
                else:
                    return list(self._magma_space.Basis())                
        elif method=='numerical':
            if prec == None:
                Y = mpmath.mpf(0.95)*self._group.minimal_height()
                prec = get_M_for_holom(Y,self._weight,self._dprec)
            if self._verbose > 0:
                print "Using numerical method with %s coefficients " % M
            B = self.basis_numerical(prec)
            res = []
            for f in B:
                l=[]
                for n in range(prec):
                    l.append(f.C(n))
                res.append(l)
            return res
        else:
            raise NotImplementedError,"Currently supported methods are 'magma' and 'numerical'! Got: %s" % method


    def assert_triangular_basis(self):
        r"""
        Check that the basis of self is upper-triangular.
        """ 
        B = self.basis()
        d = len(B)
        try:
            for i in range(d):
                for j in range(d):
                    c = B[i].Coefficient(j)
                    if(i==j and c==0  or i<>j and c<>0):
                        raise StopIteration()
        except StopIteration:
            return False
        return True
                        
    def plus_space_qseries(self,prec=None):
        r"""
        Get a basis for the Kohnen plus space, i.e. with coefficients satisfying n=0 or (-1)^(k-1/2) mod 4.
        """
        kappa=ZZ((2*self._weight-1)/2)
        if(is_even(kappa)):
            sgn=1
        else:
            sgn=-1
        d=self.dimension()
        dd=self._shimura_image.dimension()
        B=self.q_expansion_basis()
        C=dict()
        nmax=self._shimura_image.sturm_bound()**2  # should be enough
        # but in order to represent a Hecke operator T(p^2) with p not dividing the level we might need more
        for p in range(3,1+next_prime(self.level())):
            if( self.level() % p <> 0):
                break
        print "p=",p
        print "nmax=",nmax
        nmax2=nmax*p*p # How many coefficients we need to make a Hecke basis for Tp
        #print "nmax2=",nmax2
        #for j in range(1,d+1):
        #    C[j-1]=dict()
        #    for n in range(nmax2):
        #        C[j-1][n]=B[j].Coefficient(n)
        # Now have all coefficients of the basis elements.
        #print C
        B0=list();    B1=list()
        for j in range(d):
            F=B[j]
            v0=True
            for n in range(nmax):
                tn=(n*sgn) % 4
                if(F.Coefficient(n) <> 0 and tn <> 0 and tn<> 1):
                    v0=False
                    break
            if(v0):
                B0.append(B[j])
            else:
                B1.append(B[j])
        return {'plus':B0,'compl':B1}

    def Hecke_eigen_basis(self,p):
        r"""
        Return a basis of Hecke eigenfunctions of self.
        """
        raise NotImplementedError

    def HeckeMatrix(self,p=None,F=None):
        r"""
        Return the matrix of the Hecke Operator T(p^2) on
        F, which should span a subspace of self.
        """
        if(F==None):
            B=self.basis()
        else:
            B=[F]
        d=len(B)
        print "B=",B
        first_non_zero=-1
        # if we didn't supply a prime we find the first p which doesn't divide the level
        if(p==None):
            for p in prime_range(next_prime(self.level()+1)):
                if(not p.divides(self.level())):
                    break
        print "Using p=",p
        nmax=self._shimura_image.sturm_bound()**2
        
        # Check that we have linearly the basis is in upper triangular forms
        # and if not we try to permute the basis to achieve that
        ok=True
        fnz=dict()
        # if we have to permute the basis we can go back again
        P=identity_matrix(ZZ,d) 
        try:
            for i in range(nmax):
                #print "i=",i
                first_non_zero=-1
                ok=True
                for j in range(d):
                    #print "B[",j,"]=",B[j]
                    #print "first_non_zero=",first_non_zero
                    for n in range(0,first_non_zero+1):
                        cn=B[j].Coefficient(n)
                        if(self._verbose>2):
                            print "B1[",j,"][",n,"]=",cn
                        if(cn<>0):
                            ok=False
                            break
                    if(not ok):
                        # need to swap this with next element
                        if(self._verbose>2):
                            print "swap ",j,"<-->",j-1
                        for k in range(d):
                            if(k<>j-1):
                                P[j,k]=0                            
                                P[k,j]=1
                        P[j,j-1]=1; P[j-1,j]=1
                        F=B[j-1]
                        B[j-1]=B[j]; B[j]=F
                        #print "New B=",B
                        break
                    for n in range(first_non_zero+1,nmax):
                        cn=B[j].Coefficient(n)
                        #print "B2[",j,"][",n,"]=",cn
                        if(cn<>0):
                            first_non_zero=n
                            fnz[j]=n
                            break
                if(not ok):
                    break
                else:
                    raise StopIteration()
            raise ArithmeticError,"Could not bring basis in upper triangular form!"
        except StopIteration:
            pass

        if(self._verbose>2):
            print "B=",B
            print "fnz=",fnz
        V=matrix(QQ,d,d)
        for j in range(d):
            for n in range(d):
                # figure out how many multiples of F[n] is in Tp[F[j]]
                k=fnz[n]                   
                a=self.HeckeOperatorBn(B[j],p,k)
                V[j,n]=QQ(a)/QQ(B[n].Coefficient(k))
        return P*V*P**-1

    def HeckeOperatorBn(self,f,p,n):
        r"""
        
        """
        eps=self._character(-1)
        k=ZZ (QQ(2*self._weight-1)/QQ(2))
        t1=self._character(p)*kronecker(eps*n*(-1)**k,p)*p**(k-1)
        b=f.Coefficient(n*p*p)+t1*f.Coefficient(n)
        if( ZZ(p*p).divides(n)):
            nn=ZZ (QQ(n)/QQ(p*p))
            b=b+p**(2*k-1)*f.Coefficient(nn)
        return b
    

    def Ur2Operator(self,f,r,prec=12):
        r"""
        Compute the action of U(r^2) on the Fourier expansion of f at infinity. 
        """
        c=dict()
        B=self.basis()
        if(prec<len(B)):
            prec=len(B)+1
        for n in range(prec):
            c[n]=f.Coefficient(n*r**2)
        # we can also express U(r^2) in terms of the basis of self
        # assuming we have an upper-triangular basis
        v=dict()
        for d in range(len(B)):
            v[d]=c[d]
        
        S=PowerSeriesRing(QQ,'q')        
        return S(c,prec)

    def basis_numerical(self,digs=15,numc=0,SetM=None,SetY=None,**kwds):
        r"""
        Computes a numerical basis, i.e. using the automorphy method.
        Note that this includes expansions at all cusps.

        INPUT: 
        -`digs` -- number of correct digits wanted
        -`numc` -- number of coefficients requested for each basis element
        -`SetM` -- compute only this number of coefficients
        -`SetY` -- compute using this hieght of the horocycle

        """
        if SetM==None:
            SetM = kwds.get("setM",kwds.get("setm",kwds.get("Setm",None)))
        if SetY==None:
            SetY = kwds.get("setY",kwds.get("sety",kwds.get("Sety",None)))
        #raise NotImplemented,"TODO!"
        # one problem is that unless the character is trivial
        # and magma is installed we don't know the dimension
        # so we have to compute that also
        if(self._basis_numerical and self._basis_len>=numc):
            return self._basis_numerical
        d = self.dimension()
        setc=list(); pp=list()
        offset = 0
        if(self.is_cuspidal()):
            offset = 1
        for i in range(d):
            setc.append({})
            pp.append({'+':{},'-':{}})
            for j in range(d):
                setc[i][(0,j+offset)]=0
            setc[i][(0,i+offset)]=1            
            #if(self.alpha(0)[1]==1):
            pp[i]['+'][(0,0)]=0
            if self.is_holomorphic():
                pp[i]['-'][(0,0)]=0
        if(not self.is_cuspidal()):
            if(len(pp)==0):
                pp.append({'+':{(0,0):1},'-':{}})
                setc.append({})
            else:
                pp[0]['+'][(0,0)]=1
                
        #ppos=[{'+':{(0,0):1}}]
        [Y,M]=self.get_Y_and_M(digs=digs,principal_part=pp)
        if SetM:
            M = SetM
        if M < numc and numc>0:
            M = numc
        if SetY:
            Y = SetY
        # first get the matrx
        Q = M + 10
        if self._verbose>0:
            print "M,Q,Y=",M,Q,Y
        #V=setup_matrix_for_harmonic_Maass_waveforms(self,Y,M,Q,pp)
        # Fix normalizations
        if self._verbose>1:
            print "pp=",pp
            print "setc=",setc
        B = self._get_element(pp,SetC=setc,SetM=M,SetY=Y)
        if(not isinstance(B,list)):
            self._basis_numerical = [B]
        else:
            self._basis_numerical = B
        self._basis_len = M
        return self._basis_numerical 


    def xi_k_inverse_pp_mpmath(self,G,digs=6,true_value=True):
        r"""
        Computes principal parts for a set of Harmonic weak Maass forms {f_1,...,f_d} with the property
        that {f_i,g_j}=delta_{ij} where {g_1,...,g_d} is a basis of self.
        INPUT:

        - true_value -- True if we use e.g. -1 or (32)**(-1/4) etc.
                        instead of the floating point approximations.

        """
        # set the normalization from G
        pp = list()
        eps=mpmath.power(10,-digs)
        if(not isinstance(G,list)):
            GG = [G]
        else:
            GG = G
        for F in GG:
            d = {'-':{},'+':{}}
            for i in range(G._space._group.ncusps()):
                ## If we try to figure out the true value of the frst coefficients at each cusp
                if(true_value):
                    if(self.alpha(i)[1]<>1):
                        continue
                    c = mpmath.mp.mpc(F.C(i,0).conjugate())
                    x = rational_approximation(abs(c)**4,eps)
                    ar = mpmath.arg(c)/mpmath.mp.pi()
                    ar_rat=rational_approximation(ar,eps)
                    ab = mpmath.power(abs(x),mpmath.mpf(0.25))
                    cc = mpmath.mp.exp(mpmath.mpc(0,ar_rat*mpmath.mp.pi()))*ab
                    cx = cc.real
                    cy = cc.imag
                    if(abs(cx)<eps and abs(cy)<eps):
                        cc = mpmath.mpc(0,0)
                    elif(abs(cx)<eps):
                        cc = mpmath.mpc(0,cy)
                    elif(abs(cy)<eps):
                        cc = mpmath.mpc(cx,0)                    
                    d['-'][(i,0)]=cc
                else:
                    d['-'][(i,0)]=F.C(i,0).conjugate()
            pp.append(d)
        return pp

    def xi_k_inverse_pp(self,G,digs=6,true_value=True):
        r"""
        Computes principal parts for a set of Harmonic weak Maass forms {f_1,...,f_d} with the property
        that {f_i,g_j}=delta_{ij} where {g_1,...,g_d} is a basis of self.
        INPUT:

        - true_value -- True if we use e.g. -1 or (32)**(-1/4) etc.
                        instead of the floating point approximations.

        """
        # set the normalization from G
        pp = list()
        RF = RealField(self._prec)
        CF = ComplexField(self._prec)
        eps=RF(10)**-digs
        if(not isinstance(G,list)):
            GG = [G]
        else:
            GG = G
        for F in GG:
            d = {'-':{},'+':{}}
            for i in range(G._space._group.ncusps()):
                ## If we try to figure out the true value of the frst coefficients at each cusp
                if(true_value):
                    if(self.alpha(i)[1]<>1):
                        continue
                    c = F.C(i,0).conjugate()
                    x = rational_approximation(abs(c)**4,eps)
                    if isinstance(c.imag,(int,mpf)):
                        ic = c.imag
                    else:
                        ic = c.imag()
                    if ic==0:
                        ar = 0
                    else:
                        ar = c.argument()/RF.pi()
                    #ar = mpmath.arg(c)/mpmath.mp.pi()
                    ar_rat=rational_approximation(ar,eps)
                    ab = RF(abs(x))**RF(0.25)
                    #ab = mpmath.power(abs(x),mpmath.mpf(0.25))
                    cc = CF(0,ar_rat*RF.pi()).exp()*ab
                    #cc = mpmath.mp.exp(mpmath.mpc(0,ar_rat*mpmath.mp.pi()))*ab
                    cx = cc.real()
                    cy = cc.imag()
                    if(abs(cx)<eps and abs(cy)<eps):
                        cc = CF(0)
                    elif(abs(cx)<eps):
                        cc = CF(0,cy)
                    elif(abs(cy)<eps):
                        cc = CF(cx,0)                    
                    d['-'][(i,0)]=cc
                else:
                    d['-'][(i,0)]=F.C(i,0).conjugate()
            pp.append(d)
        return pp


   
    def xi_k_inverse_basis(self,digs=10,pp_in=None,**kwds):
        r"""
        Computes a set of Harmonic weak Maass forms {f_1,...,f_d} with the property
        that {f_i,g_j}=delta_{ij} where {g_1,...,g_d} is a basis of self.

        INPUT: 
        -`digs`  -- number of digits precision in result
        -`pp_in` -- principal part (if we want a specific p.part set)
        `**kwds` -- extra argument which propagates to algorithms computing coefficients

        """
        B=self.basis_numerical(digs=digs,**kwds)
        #if(not self.assert_triangular_basis()):
        #    raise ArithmeticError,"Basis is not upper triangular!"
        H = HarmonicWeakMaassFormSpace(self)
        M = H.modular_forms_subspace()
        M._weight = H.weight()
        BB = M.basis_numerical() ## we want to "project away" these 
        ## we want to know which coefficients are non-zero
        cb = dict()
        for n in range(1,len(BB)+1):
            cb[n]=0
            for f in BB:
                #cn =  f.Coefficient(n)
                cn = f.C(n)
                if(cn<>0):
                    cb[n]=cb[n]+1
            if(cb[n]<len(BB)):
                cb[n]=0
        if(cb.values().count(0)==len(cb.values())):
            raise ValueError,"Could not find good coefficients!"
        FF=list()
        setc=list()
        pp = list()
        eps=1E-6
        for G in B:
            # get the correct principal part
            p = self.xi_k_inverse_pp(G,digs=6,true_value=True)[0]
            # if we want to set some specific part of the principal parts
            if(pp_in):
                for (r,n) in pp_in.keys():
                    p[(r,n)]=pp_in[(r,n)]
            pp.append(p)                
            ## If G is in the + space we set some c's too...
            sc=dict()
            for j in range(1,4):
                if self._verbose>1:
                    print "C(0",j,")=",G.C(0,j)
                if(abs(G.C(0,j))<eps):
                    sc[(0,-j)]=0
            # We have to guess how to get rid of any holomorphic forms...
            for n in range(1,len(BB)+1):
                if(cb[n]==len(BB)):
                    sc[(0,n)]=0
                    break
            setc.append(sc)
            # need the correct constant-zero part too... 

        H._verbose=2
        if self._verbose>1:
            print "princ_part=",pp
            print "setc=",setc
        FF = H._get_element(pp,digs=digs,SetC=setc,**kwds) #,dbase_prec=prec)
        return FF
    
class AutomorphicFormElement(SageObject):
    r"""
    A member of a generic space of automorphic forms.
    

    """
    def __init__(self,M,C=None,prec=53,principal_part=None,verbose=0,**kwds):
        r"""
        Initialize an automorphic form.
        
        INPUT:

        
        -''M'' -- Space of automorphic forms
        -''k'' -- Weight.
        -''C''-- Fourier coefficients
        -''prec'' -- integer, precision (if given by construction, default None)
        -''principal_part -- principal part in dictionary 
           '+' : principal part of c^{+}
           '-' : principal part of c^{-}

        EXAMPLES:


        """
        # TODO: fix this. It doesn't work for inherited classes...
        #if(not isinstance(M,AutomorphicFormSpace)):
        #    raise TypeError,"Need an element of AutomorphicFormSpace. got %s" %M
        if M._verbose>1:
            print "MM=",M
            print "dim=",M._rdim
        if(not hasattr(M,"_is_automorphic_form_space")):
             raise TypeError,"Need an element of AutomorphicFormSpace. got %s" %M
        d1=M._rdim
        d2=len(M._group.cusps())
        if C <> None:
            # We need the correct length of the coefficient vector
            if len(C.keys()) > d1 or (len(C.keys())<d1 and self._sym_type==None):
                # If we have smaller amount we believe there is a symmetry at work...
                #or (d1==1 and len(M._group._cusps)<>len(C.keys()))):
                raise ValueError,"Coefficient vector of wrong format! Got length=%s" % len(C)
            self._coeffs=C
        else:
            self._coeffs = {i : {j:{} for j in range(d2)} for i in range(d1)}

        self._space=M
        self._prec=prec
        self._base_ring=MPComplexField(prec)
        self._verbose=verbose
        self._maxdigs=prec # the number of digits needed to be displayed to print all digits of the coefficients
        if(principal_part):
            self._principal_part=principal_part
        else:
            self._principal_part={'+':{},'-':{}}
        self._is_automorphic_form=True
        if(not hasattr(self,"_class_name")):
            self._class_name="AutomorphicFormElement"
        

    def __reduce__(self):
        r"""
        Used for pickling self.
        """
        #return(HarmonicWeakMaassFormElement,(self.space,self.principal_part.items(),self.coeffs,self.prec))
        return(AutomorphicFormElement,(self._space,self._coeffs,self._prec,self._principal_part))
    
    def __cmp__(self,other):
        r"""
        Compare self to other
        """
        if(not isinstance(other,type(self))):
            return False
        eq = (self._space == other._space)
        eq = eq and (self._principal_part == other._principal_part)
        if(not eq):
            return False
        # need to check coefficients
        if(self._coeffs.keys() <> other._coeffs.keys()):
            return False
        for r in self._coeffs.keys():
            if(self._coeffs[r].keys() <> other._coeffs[r].keys()):
                return False
            for n in self._coeffs[r].keys():
                if(self._coeffs[r][n].keys() <> other._coeffs[r][n].keys()):
                    return False
        
    def  _repr_(self):
        r""" Return string representation of self.

        EXAMPLES:

            sage: WR=WeilRepDiscriminantForm(11,dual=True)
            sage: M=VVHarmonicWeakMaassFormSpace(WR,0.5,100)
            sage: PP={(7/22,0):1}
            sage: F=M.get_element(PP,12);F
            Element of Space of Vector-Valued harmonic weak Maass forms on Modular Group SL(2,Z) of weight 1/2  and dimension 10.
            Representation is Dual of Weil representation of the discriminant form given by ZZ/22ZZ with quadratic form Q(x)=11*x**2 mod 1. with principal part: q^-5/44
                    
        """
        
        s="Element of "+str(self._space)
        return s

    def _latex_(self):
        r""" Return LaTeX string representation of self.

        EXAMPLES:
        """
        return "Element of "+latex(self.space)


    def __add__(self,G):
        r"""
        Add self to G.
        """
        ## check that -G and self are the same type of modular form
        return self._lin_comb(G,1,1)
        ok=True
        if(not hasattr(G,'_is_automorphic_form')):
            if self._verbose>0:
                print "No autom form!"
            ok = False
        if G._space <> self._space:
            if self._verbose>0:
                print "Not same space as self! L:{0}, R:{1}".format(self._space,G._space)
            ok = False
        if(not ok):
            raise NotImplementedError,"Addition of elements of type: %s and %s are not implemented!" %(type(self),type(G))
        pp1 = self._principal_part
        pp2 = G._principal_part
        c1  = self._coeffs
        c2  = G._coeffs
        p = dict()
        p['+']=dict()
        p['-']=dict()
        c = dict()
        ## we truncate to the smaller number of coefficients
        for r in c1.keys():
            c[r]=dict()
            for j in c1[r].keys():
                c[r][j]=dict()
                for n in c1[r][j].keys():
                    if c2[r][j].has_key(n):
                        if self._verbose>1:
                            print "adding ",r,j,n
                        c[r][j][n]=c1[r][j][n]+c2[r][j][n]
                        if self._verbose > 1 and n==1:
                            print c1[r][j][n],"+",c2[r][j][n],"=",c[r][j][n]
        ## merge the principal parts
        k1 = pp1['+'].keys(); k2 = pp1['+'].keys();  k1.extend(k2)
        for r in k1:
            if(k1.count(r)>0):
                k1.remove(r)
        for (r,n) in k1:
            t=0
            if(pp1['+'].has_key((r,n))):
                t=t+pp1['+'][(r,n)]
            if(pp2['+'].has_key((r,n))):
                t=t+pp2['+'][(r,n)]            
            p['+'][(r,n)]=t
        k1 = pp1['-'].keys(); k2 = pp1['-'].keys();  k1.extend(k2)
        for r in k1:
            if(k1.count(r)>0):
                k1.remove(r)
        for (r,n) in k1:
            t=0
            if(pp1['-'].has_key((r,n))):
                t=t+pp1['-'][(r,n)]
            if(pp2['-'].has_key((r,n))):
                t=t+pp2['-'][(r,n)]            
            p['-'][(r,n)]=t

        ### To actually construct an element of the correct subclass is not so easy, especially if it has to work in attached files as well as in a standard installation...
        F = eval(self._class_name+'(self._space,c,self._prec,p)')
        #return AutomorphicFormElement(self._space,self._coeffs,self._prec,self._principal_part)
        return F


    def __mul__(self,a):
        r"""
        return a*self
        """
        return self._lin_comb(self,a)

    def __copy__(self):
        r"""
        Return copy of self. Probably silly way of copying but it works...
        """
        #print "COpying ",self
        s=dumps(self)
        return loads(s)
    
    def __div__(self,a):
        r"""
        return self/a
        """
        if a==0:
            raise ZeroDivisionError
        return self._lin_comb(self,a**-1)




    

    def _lin_comb(self,G,a,b=0):
        r"""
        Return a*F+b*G
        """
        ok=True
        res=copy(self)
        if not hasattr(G,'_is_automorphic_form'):
            ok = False
        if G._space <> self._space:
            ok = False
        if not ok:
            raise NotImplementedError,"Addition of elements of type: %s and %s are not implemented!" %(type(self),type(G))
        pp1 = self._principal_part
        pp2 = G._principal_part
        c1  = self._coeffs
        c2  = G._coeffs
        p = dict()
        p['+']=dict()
        p['-']=dict()
        c = dict()
        ## we truncate to the smaller number of coefficients

        #print "a,typea=",a,type(a)
        #print "b,typeb=",b,type(b)
        try:
            aa=self._base_ring(a)
        except:
            aa = self._base_ring(a.real(),a.imag())
        try:
            bb=self._base_ring(b)
        except:
            bb = self._base_ring(b.real(),b.imag())


        for r in c1.keys():
            c[r]=dict()
            for j in c1[r].keys():
                c[r][j]=dict()
                for n in c1[r][j].keys():
                    if(c2[r][j].has_key(n)):
                        #print "adding ",r,j,n
                        #print "aa.parent()=",aa.parent(),type(aa)
                        #print "bb.parent()=",bb.parent(),type(bb)
                        #print "c1.parent=",c1[r][j][n].parent(),type(c1[r][j][n])
                        #print "c2.parent=",c2[r][j][n].parent(),type(c2[r][j][n])
                        c[r][j][n]=aa*c1[r][j][n]+bb*c2[r][j][n]
                        #if(n==1):
                        #    print c1[r][j][n],"+",c2[r][j][n],"=",c[r][j][n]
        ## merge the principal part
        k1 = pp1['+'].keys(); k2 = pp2['+'].keys();  k1.extend(k2)
        for r in k1:
            if(k1.count(r)>1):
                k1.remove(r)
        for (r,n) in k1:
            t=0
            if(pp1['+'].has_key((r,n))):
                t=t+a*pp1['+'][(r,n)]
            if(pp2['+'].has_key((r,n))):
                t=t+b*pp2['+'][(r,n)]            
            p['+'][(r,n)]=t
        k1 = pp1['-'].keys(); k2 = pp1['-'].keys();  k1.extend(k2)
        for r in k1:
            if(k1.count(r)>1):
                k1.remove(r)
        for (r,n) in k1:
            t=0
            if(pp1['-'].has_key((r,n))):
                t=t+a*pp1['-'][(r,n)]
            if(pp2['-'].has_key((r,n))):
                t=t+b*pp2['-'][(r,n)]            
            p['-'][(r,n)]=t

        res._coeffs=c
        res._principal_part=p
        ### To actually construct an element of the correct subclass is not so easy, especially if it has to work in attached files as well as in a standard installation...
        #F = eval(self._class_name+'(self._space,c,self._prec,p)')

        #return AutomorphicFormElement(self._space,self._coeffs,self._prec,self._principal_part)
        return res

    def set_verbositiy(self,verbose):
        r"""

        """
        self._verbose=verbose
    
    def prec(self):
        r"""
        return precision (number of coefficients) of self
        """
        return self._prec

    def base_ring(self):
        return self._base_ring
    
    def space(self):
        r""" Return the ambient space of self.
        """
        return self._space
    
    def coeffs(self):
        r""" Return the coefficients of self.
        """
        return self._coeffs
    
    def prec(self):
        r""" Return the precision of self.
        """
        return self._prec

    def weight(self):
        r""" Return the weight of ambient space.
        """
        return self._space._weight

    def character(self):
        r""" Return the character of ambient space.
        """
        return self._space._character

    def multiplier(self):
        r""" Return the multiplier of ambient space.
        """
        return self._space._multiplier

    def principal_part(self):
        r""" Return the principal part  of self.
        """
        return self._principal_part        

    def is_holomorphic(self):
        r"""
        Return True if self is a mamber of a space of holomorphicforms, otherwise False. 
        """
        return self._space._holomorphic

    def is_weak(self):
        r"""
        Return True if self is a weak form, i.e. has pole(s) at oo, otherwise False. 
        """
        return self._space._weak

    def is_harmonic(self):
        r"""
            Return True if self is a harmonic weak maassform, i.e. has pole(s) at oo, otherwise False. 
        """
        return self._space._harmonic

    def C(self,a=None,b=None,c=None):
        r""" Return coefficient nr. n at cusp nr. i of component r of self.
        C(n) = Coeff n at cusp 0
        C(i,n) = Coeff n at cusp i
        C(r,i,n) = Coeff n at cusp i and component r
        """
        if(c<>None):
            r=a; i=b; n=c
        elif(b<>None):
            r=0; i=a; n=b
        elif(a<>None):
            r=0; i=0; n=a
        else:
            raise ValueError,"Need to supply a valid index for the coefficient!"
        C = self._coeffs
        if(not C.has_key(r)):
            return None
        if(not C[r].has_key(i)):
            return None
        if(not C[r][i].has_key(n)):
            return None         
        return C[r][i][n]

    def my_zzcmp(self,a,b):
        r"""
        Compare a and b.
        """
        if(a=='-0' and ZZ(b)==0):
            return -1
        if(b=='-0' and ZZ(a)==0):
            return 1        
        return cmp(ZZ(a),ZZ(b))
    
    def list_coefficients_old(self,nmax,norm=False,plus=True,minus=True,cusp=None):
        r"""
        List coefficients of self.
        """
        have_minus=self._principal_part.get('-')<>{}
        
        for r in self._coeffs.keys():
            print ""
            for j in self._coeffs[r].keys():
                if cusp<>None and j<>cusp:
                    continue
                print ""
                l=self._coeffs[r][j].keys()
                l.sort(cmp=self.my_zzcmp)
                #c0=self._coeffs[r][j][0]
                #c1=self._coeffs[r][j][1]
                ## Scaling factor for positive coefficients
                # If we have set the constant term then C(r)(0) is in fact related to the non-holomorphic part and not the holomorphic one
                a0_neg=False
                c_neg_norm=1 
                c_norm=1 
                if self._principal_part.get('+',{}).has_key((j,0)):
                    if not self._principal_part.get('-',{}).has_key((j,0)):
                        a0_neg=True
                print "a0_neg=",a0_neg
                if norm:
                    for norm_index in range(0,len(l)):
                        if not (norm_index in l):
                            continue
                        if norm_index==0 and a0_neg:
                            continue
                        nal = norm_index + self._space.alpha(j)[0]
                        if nal < 0:
                            continue
                        c_norm=self._coeffs[r][j][norm_index]#/mpmath.sqrt(abs(nal))
                        if abs(c_norm)>mpmath.eps() and abs(c_norm)>0.0:
                            print "For normalization we use:"
                            print "c0=c[",j,norm_index,"]=",c_norm
                            print "c0**4=",c_norm**4 #mpmath.power(c_norm,4)
                            print "c0**-4=",c_norm**-4 #mpmath.power(c_norm,-4)
                            break
                    if self._verbose>1:
                        print "c_norm=",c_norm
                        # If the first non-zero coefficients is too small we don't want to normalize
                    for norm_index in range(0,len(l)):
                        nal = -norm_index + self._space.alpha(j)[0]
                        print "nal(",norm_index,')=',nal
                        if(not (-norm_index in l)):
                            continue
                        if nal>0:
                            continue
                        if nal==0 and not a0_neg:
                            continue  #cm1=self._principal_part['-'][(j,0)]                       
                        cm1=self._coeffs[r][j][-norm_index]
                        if abs(cm1)>mpmath.eps() and abs(cm1)>0.0:
                            if abs(nal)>0:
                                c_neg_norm=cm1*self._base_ring(nal).abs().sqrt()
                            else:
                                c_neg_norm=cm1
                            print "c-[",j,",",-norm_index,"]=",cm1
                            break
                    if self._verbose>1:
                        print "cm_neg_norm=",c_neg_norm
                ## Scaling factor for negative coefficients
                #cm1=self._coeffs[r][j][-1]
                #if(abs(cm1)>mpmath.eps() and cm1<>0):
                for n in l:
                    #if(is_int(n)):
                    al = self._space.alpha(j)[0]
                    nal = ZZ(n) + al
                    st_nal=sci_pretty_print(nal,3)
                    #elif(n=='-0'):
                    #    nal=0
                    if(abs(nal)>nmax):
                        continue
                    if(plus == False and nal>=0):
                        continue
                    if(minus== False and nal<0):
                        continue                    
                    if al<>0:
                        if(norm and abs(c_norm)>mpmath.eps() and nal>0): # and n<>norm_index):
                            c=self._coeffs[r][j][n]/c_norm
                        elif(norm and nal<0 and abs(c_neg_norm)>mpmath.eps()):
                            #c=self._coeffs[r][j][n]*mpmath.sqrt(abs(nal))/c_neg_norm
                            c=self._coeffs[r][j][n]*sqrt(abs(nal))/c_neg_norm
                        else:
                            c=self._coeffs[r][j][n]
                        if(len(self._coeffs.keys())>1):
                            print "C[",r,",",j,",",n,"]=",c
                        else:
                            print "C[",j,",",n,":"+st_nal+"]=",c
                    else:
                        if(self._verbose>1):
                            print "j,n=",j,n,nal,a0_neg
                        if(a0_neg):
                            c_minus=self._coeffs[r][j][n]
                            c_plus=self._principal_part['+'].get((j,0),0)
                        else:
                            c_minus=self._principal_part['-'].get((j,0),0)
                            c_plus=self._coeffs[r][j][n]
                        if len(self._coeffs.keys())>1:
                            if have_minus:
                                print "C[",r,",",j,",-",n,"]=",c_minus
                            print "C[",r,",",j,",+",n,"]=",c_plus
                        else:
                            if have_minus:
                                print "C[",j,",-",n,"]=",c_minus
                            print "C[",j,",+",n,"]=",c_plus
                            


    def list_coefficients(self,N,cusp="all",norm=False,print_part=['+','-'],component=-1,out_prec=53):
        r"""
        List coefficients of self.

        INPUT::

        
        - `N` -- integer, number of coefficients to list
        - `cusp` -- integer, list only coefficients for this cusp
        - `norm` -- Bool, print normalized expansions
        - `print_part` -- '+','-' or '+,-' (default) 
        - `component` -- integer, print only this component (applies to vector-valued forms)
        """
        have_minus=self._principal_part.get('-')<>{}
        ###
        s=""
        if component<0:
            for r in self._coeffs.keys():
                s+=self.list_coefficients(N=N,cusp=cusp,norm=norm,print_part=print_part,component=r,out_prec=out_prec)
            return s
        if cusp=="all":
            for j in self._coeffs[component].keys():
                s+=self.list_coefficients(N=N,cusp=j,norm=norm,print_part=print_part,component=component,out_prec=out_prec)
                #assert isinstance(cusp,(int,Integer))
            return s
        C = {} #copy(self._coeffs.get(component,{}).get(cusp,{}))
        #Cplus = []
        #Cminus = []
        l = self._coeffs.get(component,{}).get(cusp,{}).keys()
        l_plus = []; l_minus=[]
        for n in l:
            if abs(n)>N:
                continue
            if self.my_zzcmp(n,0)>0:
                l_plus.append(n)
            else:
                l_minus.append(n)
        #for (j,n) in  self.principal_part().get('+',{}).keys(): 
        #    if j==cusp:
        #        l_plus.append(n); C[n]=self.principal_part()['+'][(j,n)]
        #for (j,n) in  self.principal_part().get('-',{}).keys(): 
        #    if j==cusp:
        #        l_minus.append(n); C[n]=self.principal_part()['-'][(j,n)]
        l_plus.sort(cmp=self.my_zzcmp)
        l_minus.sort(cmp=self.my_zzcmp)
        scaling_factors=[]
        RF = RealField(self._prec)
        new_prec = self._prec
        if '+' in print_part:  ## Print the positive part:
            if norm:
                norm_plus = self.principal_part()['+'].get((cusp,0),F.C(cusp,0))
                if hasattr(norm_plus,"prec"):
                    new_prec = norm_plus.prec()
                if abs(norm_plus) < RR(2.0)**(-self.prec()/2.0):
                    try:
                        for n in l_plus:
                            if n>0:
                                norm_plus = F.C(cusp,n)
                                if abs(norm_plus) > RR(2.0)**(-self.prec()/2.0):
                                    raise StopIteration()
                        raise ArithmeticError,"Could not find non-zero positive coefficient to normalize with!"
                    except StopIteration:
                        pass                
            else:
                norm_plus = 1 
            if hasattr(norm_plus,"prec"):
                new_prec = norm_plus.prec()
            if out_prec>new_prec:
                new_prec = out_prec
            RF = RealField(new_prec)
            CF = ComplexField(new_prec)
            norm_plus = CF(1)
            al = self._space.alpha(cusp)[0]
            al0 = QQ(al)
            if self._verbose>0:
                print "al=",al
                print "al0=",al0
                print "al0.denom()=",al0.denominator()
            if al0.denominator()<1000: 
                al = al0
            for n in l_plus:
                nal = n + al
                if n<0:
                    c = C.get(component,n,0)
                else:
                    c = self.C(component,cusp,n)
                if isinstance(c,(int,mpf,complex)): 
                    c = CF(c.real,c.imag)
                elif hasattr(c,"real"):
                    c = CF(c.real(),c.imag())
                else:
                    c = CF(c)
                c = c/norm_plus

                s+="C^{{+}}_{{ {0} }} ({1}) = {2} \n ".format(cusp,nal,c)
        if '-' in print_part:
            ## For the negative coefficients we scale with c_(-1) n^(1-k)
            al = self._space.alpha(cusp)[0]
            al0 = QQ(al)
            if al0.denominator()<1000: 
                al = al0
            if norm:
                norm_minus = self.C(component,cusp,-1)  #self.principal_part()['-'].get((cusp,0),F.C(cusp,0))
                if abs(norm_minus) < RF(2.0)**(-self.prec()/2.0):
                    try:
                        for n in l_minus:
                            if n>0:
                                norm_minus = self.C(component,cusp,n)
                                if abs(norm_minus) > RF(2.0)**(-self.prec()/2.0):
                                    raise StopIteration()
                        raise ArithmeticError,"Could not find non-zero negative coefficient to normalize with!"
                    except StopIteration:
                        pass                
            else:
                norm_minus = CF(1)             
            for n in l_minus:    
                nal = n + al
                c = self.C(component,cusp,n)
                if isinstance(c,(int,mpf,complex)): 
                    c = CF(c.real,c.imag)
                else:
                    c = CF(c.real(),c.imag())
                if nal>0:
                    nn = abs(nal)**(1-self.weight())
                    c = c/norm_minus/nn

                s+="C^{{-}}_{{ {0} }} ({1}) = {2} \n ".format(cusp,n,c)
        return s

    def print_table_of_coeffs(self,nmax,norm=False,plus=True,minus=True,cusp=None,table_format="table"):
        r"""
        Print a Latex table of coefficients.
        """
        zl=1E-8  ### smaller coefficients than this are treated as zero
        th="\\begin{"+table_format+"}[h]"
        res=""
        ph="\\hphantom{$-$}"
        prec = self.space().prec()
        eps = 2.0*10.0**(-prec*ln(2.0)/ln(10.0))
        for r in self._coeffs.keys():
            ## Make one table per vector component
            tbl=" $n$ "
            cols="|l"
            for j in self._coeffs[r].keys():
                if(cusp<>None and j<>cusp):
                    continue
                tbl+= "& $p="+latex(self._space._group.cusps()[j])+"$ "
                cols+="|l"
            tbl+="\\\\ \n"
            tbl+="\hline \\noalign{\smallskip} \n"
            rh="\\begin{tabular}{"+cols+"|} \n"
            res=res+th+rh

            # Use the union of indices for all cusps as index set.
            Nlist=[]
            for j in self._coeffs[r].keys():
                if(cusp<>None and j<>cusp):
                    continue
                l=self._coeffs[r][j].keys()
                for n in l:
                    al = self._space.alpha(j)[0]
                    nal = ZZ(n) + al
                    if( ((not plus) and nal>0) or ((not minus) and nal<0)):
                        continue
                    if(abs(n)>nmax):
                        continue
                    if(Nlist.count(n)==0):                        
                        Nlist.append(n)
            Nlist.sort(cmp=self.my_zzcmp)
            ### First get the scaling factors
            pos_factor=dict()
            neg_factor=dict()
            a0_neg=dict()
            for j in self._coeffs[r].keys():
                if(cusp<>None and j<>cusp):
                    continue
                # # Scaling factor for positive coefficients
                # If we have set the constant term then C(r)(0) is in fact related to the non-holomorphic part and not the holomorphic one
                a0_neg[j]=False
                if(self._principal_part['+'].has_key((j,0))):
                    if(not self._principal_part['-'].has_key((j,0))):
                        a0_neg[j]=True
                for norm_index in range(0,max(Nlist)):
                    nal = norm_index + self._space.alpha(j)[0]
                    if( (not (norm_index in l))  or (norm_index==0 and a0_neg[j]) or (nal <0)):
                        continue
                    c_norm=self._coeffs[r][j][norm_index]#/mpmath.sqrt(abs(nal))
                    if(abs(c_norm)>eps and abs(c_norm)>0.0):
                        print "For normalization we use:"
                        print "c0=c[",j,norm_index,"]=",c_norm
                        print "c0**4=",c_norm**4
                        print "c0**-4=",c_norm**-4
                        break
                pos_norm_i=norm_index
                pos_factor[j]=(pos_norm_i,c_norm)
                c_neg_norm=1 #mpmath.mpf(1)
                ## Scaling factor for negative coefficients
                # If the first non-zero coefficients is too small we don't want to normalize
                for norm_index in range(0,max(Nlist)):                    
                    nal = -norm_index + self._space.alpha(j)[0]
                    # print "nal=",nal
                    if(not (-norm_index in l)):
                        continue
                    if(nal>=0):
                        continue
                    if(norm_index==0 and not a0_neg[j]):
                        continue  #cm1=self._principal_part['-'][(j,0)]                       
                    cm1=self._coeffs[r][j][-norm_index]
                    if(abs(cm1)>eps and abs(cm1)>0.0):
                        if(abs(nal)>0):
                            if hasattr(nal,"ae"):
                                c_neg_norm=cm1*mpmath.sqrt(abs(nal))
                            else:
                                c_neg_norm=cm1*sqrt(abs(nal))
                        else:
                            c_neg_norm=cm1
                        print "c-[",j,",",-norm_index,"]=",cm1
                        # print "c_neg_norm=",c_neg_norm
                        break
                neg_norm_i=-norm_index
                print "cm_neg_norm=",c_neg_norm
                neg_factor[j]=(neg_norm_i,c_neg_norm)
                # Now we can finally make the table
            for n in Nlist:
                if(n<>0):
                    row=" $"+str(n)+"$"
                    for j in self._coeffs[r].keys():
                        if(cusp<>None and j<>cusp):
                            continue
                        al = self._space.alpha(j)[0]
                        nal = ZZ(n) + al                
                        if(plus == False and nal>=0):
                            continue
                        if(minus== False and nal<0):
                            continue                    
                        if(norm and abs(c_norm)>eps and nal>0): # and n<>norm_index):
                            c=self._coeffs[r][j][n]/pos_factor[j][1]
                        elif(norm and nal<0 and abs(c_neg_norm)>eps):
                            if hasattr(nal,"ae"):
                                c=self._coeffs[r][j][n]*mpmath.sqrt(abs(nal))/neg_factor[j][1]
                            else:
                                c=self._coeffs[r][j][n]*sqrt(abs(nal))/neg_factor[j][1]                             
                        else:
                            c=self._coeffs[r][j][n]
                        s = norm_sci_pretty_print(c,16,latex_pow=True,zero_lim=zl)
                        if(s.lstrip()[0]<>"-"):                            
                            row+="& "+ph+"$"+s+"$ "
                        else:
                            row+="& $"+s+"$ "
                    row+="\\\\ \n "
                else:
                    ## Have to distinguish between alpha(j)=0 and <>0 
                    ## Add two rows of c(0)
                    print "a0neg=",a0_neg
                    row=""
                    if(minus):
                        row+=" $a^{-}(0)$ "
                        for j in self._coeffs[r].keys():
                            if(cusp<>None and j<>cusp):
                                continue
                            al = self._space.alpha(j)[0]
                            nal = ZZ(n) + al                
                            if(nal<>0):
                                row+=" & "
                            else:
                                if(a0_neg[j]):
                                    c_minus=self._coeffs[r][j][n]
                                else:
                                    c_minus=self._principal_part['-'].get((j,0),0)
                                s = norm_sci_pretty_print(c_minus,16,latex_pow=True,zero_lim=zl) 
                                if(s.lstrip()[0]<>"-"):                            
                                    row+="& "+ph+"$"+s+"$ "
                                else:
                                    row+="& $"+s+"$ "
                            #print "c_minus(",j,")=",c_minus
                            #print "s=",s
                            #row+=s
                        row+="\\\\ \n "
                        row+="$a^{+}(0)$ "
                    else:
                        row+="$0$ "
                    for j in self._coeffs[r].keys():
                        if(cusp<>None and j<>cusp):
                            continue
                        al = self._space.alpha(j)[0]
                        nal = ZZ(n) + al                
                        if(nal<>0):
                            c_plus = self._coeffs[r][j][n]
                        else:
                            if(a0_neg[j]):
                                c_plus=self._principal_part['+'][(j,0)]
                            else:
                                c_plus=self._coeffs[r][j][n]
                        if(norm):
                            c_plus=c_plus/pos_factor[j][1]
                        print "c_plus=",c_plus
                        ## We now add two rows
                        s=norm_sci_pretty_print(c_plus,16,latex_pow=True,zero_lim=zl)
                        if(s.lstrip()[0]<>"-"):                            
                            row+="& "+ph+"$"+s+"$ "
                        else:
                            row+="& $"+s+"$ "
                    row+="\\\\ \n"
                tbl+=row
            tbl+="\end{tabular} \end{"+table_format+"}"
            res+=tbl
        return res

    
    def print_as_q_exp(self,nmax,base_ring=QQ,eps=1E-12):
        r"""
        Print self as q-expansions with coefficients truncated to nearest integer (if error is less than prec)

        INPUT:
        
        - `nmax` -- number of terms in the q-expansion
        - `base_ring` -- Base ring for q-series
        - `eps` -- precision for rational approximation 

        """
        qexps=list()
        ## We try with rationals...
        x=ZZ['x'].gen()
        QQi=QQ.extension(x**2+1,'i')
        if base_ring=='QQi':
            base_ring=QQi
        S=PowerSeriesRing(base_ring,'q')        
        coeff_f_minus=dict()
        coeff_f_plus=dict()
        eps = 2.0*10.0**(-self.space().prec()*ln(2.0)/ln(10.0))
        print "Principal part:"+self.print_principal_part()
        const=dict()
        for r in self._coeffs.keys():
            coeff_f_minus[r]=dict()
            coeff_f_plus[r]=dict()
            const[r]=dict()
            for j in self._coeffs[r].keys():
                coeff_f_minus[r][j]=list()
                coeff_f_plus[r][j]=list()
                #print "c0=",c0
                if(j<>0):
                    c0=self._coeffs[r][j][0]
                    n0=0
                    if(abs(c0)==0):
                        c0=self._coeffs[r][j][1]
                        n0=1
                else:
                    c0=1
                    n0=0
                if base_ring==QQ:
                    cr = rational_approximation(c0**4,eps)
                else:
                    c0=base_ring(c0)
                    cr=base_ring(c0**4)
                if hasattr(c0,"ae"):
                    ar = rational_approximation(mpmath.arg(c0)/mpmath.mp.pi(),eps)
                elif hasattr(c0,"argument"):
                    ar = rational_approximation(c0.argument()/c0.base_ring().pi(),eps)
                else:
                    if c0>0:
                        ar = 0
                    else:
                        ar = RR.pi()
                
                const[r][j]=[cr,ar]
                #print "c0=",c0
                # First f_j^-
                for n in range(0,nmax):
                    al=-n+self._space.alpha(j)[0]
                    if al<eps and self.is_holomorphic():
                        continue
                    if base_ring==QQ:
                        c=rational_approximation(self._coeffs[r][j][-n]/c0,eps)
                    else:
                        c=base_ring(self._coeffs[r][j][-n]/c0)
                    coeff_f_minus[r][j].append(c)
                    #print "f[",j,",",n,"]^-=",c
                # Then f_j^+
                for n in range(0,nmax):
                    al=n+self._space.alpha(j)[0]
                    if al<-eps:
                        continue
                    tmp=self._coeffs[r][j][n]/c0
                    print "tmp=",tmp,type(tmp)
                    if isinstance(tmp,int):
                        rtmp=tmp
                        itmp=0
                    else:
                        rtmp=tmp.real()
                        if hasattr(tmp,"complex_embedding"):
                            itmp=tmp.complex_embedding(self.space().prec()).imag()
                        else:
                            itmp=tmp.imag()
                    if base_ring==QQ:
                        c=rational_approximation(rtmp,eps)
                    if base_ring==QQi:
                        cx=rational_approximation(rtmp,eps)
                        cy=rational_approximation(itmp,eps)
                        c=cx+i*cy
                    else:
                        c=base_ring(tmp)
                    if base_ring==QQ:
                        coeff_f_plus[r][j].append(c.real())                
                    else:
                        coeff_f_plus[r][j].append(c)                
                    #print "f^+[",j,",",n,"]^+=",c
                if len(self._coeffs.keys())>1:
                    if(coeff_f_minus[r][j].count(0)==len(coeff_f_minus[r][j])):
                        fm="0"
                    else:
                        fm=str(coeff_f_minus[r][j])
                    fp0 = str(S(coeff_f_plus[r][j]).add_bigoh(nmax))
                    al=self._space.alpha(j)[0]
                    if(al<>0):
                        fp="q^{"+str(rational_approximation(al))+"}("+fp0+")"
                    else:
                        fp=fp0
                    print "f[",r,",",j,"]^+=",fp
                else:
                    if(coeff_f_minus[r][j].count(0)==len(coeff_f_minus[r][j])):
                        fm="0"
                    else:
                        fm=str(coeff_f_minus[r][j])
                    if not self.is_holomorphic():
                        print "f[",j,"]^-=",fm
                    al=self._space.alpha(j)[0]                    
                    fp0 = str(S(coeff_f_plus[r][j]).add_bigoh(nmax))                        
                    if(al<>0):
                        fp="q^{"+str(rational_approximation(al))+"}("+fp0+")"
                    else:
                        fp=fp0
                    print "f[",j,"]^+=",fp
        for r in self._coeffs.keys():
            for j in self._coeffs[r].keys():
                if const[r][j][0]<>0 and const[r][j][1]<>0:
                    print "c[",r,j,"]^4= %s arg = %s pi" % (const[r][j][0],const[r][j][1])

    def print_principal_part(self):
        r""" Print the principal part of self as q-series.
        """
        s=""
        sp=""
        for (r,n) in self._principal_part['+']:
            a=self._principal_part['+'][(r,n)]
            if a<>0:
                x=QQ(n+self._space.alpha(r)[0])            
            if a<>1:
                if a>0 and len(sp)>0:
                    ast="+"+str(a)
                else:
                    ast=str(a)
                if x<>0:
                    sp=sp+ast+"q^{"+str(x)+"}"
                else:
                    sp=sp+ast
            else:
                if x<>0:
                    sp=sp+"q^{"+str(x)+"}"
                else:
                    sp=sp+"1"
        s=s+sp+""
        for (r,n) in self._principal_part['-']:
            s+=' + non-holomorphic principal part...(todo: print better)'
        return s


    def xi_k_inverse_G(self,prec=53,pp_in=None):
        r"""
        Compute xi_k^-1(g) for self.
        """
        M = self._space
        pp = [M.xi_k_inverse_pp(self)]
        H = HarmonicWeakMaassFormSpace(M) 
        F = H.get_element(pp,prec)
        eps = 2.0*10.0**(-self.space().prec()*ln(2.0)/ln(10.0))
        if(pp_in):
            p = dict()
            for (r,n) in pp_in.keys():
                p[(r,n)]=pp_in[(r,n)]
            pp.append(p)                
        ## If G is in the + space we set some c's too...
        sc=dict()
        for j in range(1,4):
            if self._verbose > 1:
                print "C(0",j,")=",G.C(0,j)
            if(abs(G.C(0,j))<eps):
                sc[(0,-j)]=0
            # We have to guess how to get rid of any holomorphic forms...
        for n in range(1,len(BB)+1):
            if(cb[n]==len(BB)):
                sc[(0,n)]=0
                break
        setc.append(sc)
        F = H.get_element(pp,SetC=setc,prec=prec)
        return F


            
class HalfIntegralWeightFormElement(AutomorphicFormElement):
    r"""
    A member of the space of half-integral weight modular forms forms (with theta multiplier as default).
    """ 
    def __init__(self,M,C=None,prec=53,principal_part=None,**kwds):
        r"""
        Initialize an automorphic form.
        
        INPUT:

        
        -''M'' -- Space of automorphic forms
        -''k'' -- Weight.
        -''C''-- Fourier coefficients
        -''prec'' -- integer, precision (if given by construction, default None)
        -''principal_part -- principal part in dictionary 
           '+' : principal part of c^{+}
           '-' : principal part of c^{-}
        EXAMPLES:
        """
        if(hasattr(C,'get_magma_attribute')):
            ## We can initialize a half integral weight form from a magma form (q-expansion)
            self._magma_form=C
            coeff = dict()
            coeff[0] = C.Coefficients()
        else:
            self._magma_form=None
            coeff = C
        self._class_name = "HalfIntegralWeightFormElement"
        AutomorphicFormElement.__init__(self,M,C=coeff,prec=prec,principal_part=principal_part)

class HarmonicWeakMaassFormElement(AutomorphicFormElement):    
    r"""
    Create an Harmonic Weak Maass form.
    """ 

    def __init__(self,M,C=None,prec=None,principal_part=None,**kwds):
        r"""
        See ``HarmonicWeakMaassFormElement`` for full documentation.

        Initialize an harmonic weak Maass form.
        
        INPUT:

        
            -''M'' -- Space of automorphic forms
            
            -''k'' -- Weight.
            -''C''-- Fourier coefficients
            -''prec'' -- integer, precision (if given by construction, default None)

        EXAMPLES::
        """
        #print "typeM=",type(M)
        self._class_name ="HarmonicWeakMaassFormElement"
        AutomorphicFormElement.__init__(self,M,C,prec=prec,principal_part=principal_part)
        

        


# class AlmostHolomorphicModularForms(AutomorphicFormSpace):
#     r"""
#     Space of Harmonic weak Maass forms

#     """
#     def __init__(self,G,weight,multiplier=None,holomorphic=False,cuspidal=False):
#         holomorphic=False # otherwise we have a holomorphic modular form
class AlmostHolomorphicModularFormSpace(AutomorphicFormSpace):
    pass
        
class HarmonicWeakMaassFormSpace(AlmostHolomorphicModularFormSpace):
    r"""
    Space of Harmonic weak Maass forms.
    """
    def __init__(self,G,weight=0,multiplier=None,holomorphic=False,weak=True,cuspidal=False,verbose=0,**kwds):
        r"""
        Initialize the space of Harmonic weak Maass forms, that is the space of weakly modular
        functions with q-expansions of the form:
        f(x+iy) = c/y^{k-1} + P(q^-1) + Sum_{n} c(n) W_n(y)q^n
        where W_n(y)=1 for n>0 and Gamma(k-1,4piy) for n<0 
        (n - alpha(i)) is an integer where v(T_i)=e(alpha(i)) as usual.
        """
        if(isinstance(G,MySubgroup_class)):
            self._group=G
            self._from_group=G._G
        else:
            if is_int(G):
                self._group=MySubgroup(Gamma0(G))
                self._from_group=Gamma0(G)
            elif( hasattr(G,'is_subgroup') and G.is_subgroup(SL2Z)):
                self._group=MySubgroup(G)
                self._from_group=G
            elif( hasattr(G,'_is_space_of_automorphic_functions')):
                # We can use the inverse of the \xi_k operator to map M_{k,rho} to H_{2-k,\bar{rho}} and 
                self._group=G._group
                self._from_group=G._from_group
                self._verbose = G._verbose
                weight=QQ(2)-QQ(G.weight())
                x=G.character()
                if x.is_trivial():
                    character=trivial_character(self.level())
                elif(isinstance(x,sage.modular.dirichlet.DirichletCharacter)):
                    if(x.order()<=2):
                        character=x
                    else:
                        # the conjugate character
                        character=x.parent()[0]/x   
                elif(isinstance(x,type(trivial))):  # if x is a function
                    character = x
                else:
                    raise ValueError, "Unknown character! x:%s" % x
            else:
                raise ValueError, "Could not initialize space from G:%s" % G
            
            
        # We need a level divisible by 4
        if multiplier==None:
            if is_int(weight):
                multiplier=TrivialMultiplier(self._group)
            else:
                if(self.level() % 4 <>0):
                    raise ValueError," Need level divisible by 4. Got:%s " % self.level()
                #if ( int(2*weight) % 4 == 1):
                multiplier=ThetaMultiplier(self._group,weight=weight)
                #else:
                #    multiplier=ThetaMultiplier(self._group,dual=True)
        self._class_name ="HarmonicWeakMaassFormSpace"
        #AutomorphicFormSpace.__init__(self,GG,weight=weight,multiplier=multiplier,character=character,holomorphic=holomorphic,weak=weak,cuspidal=cuspidal,dprec=dprec,verbose=verbose)
        AutomorphicFormSpace.__init__(self,self._group,weight=weight,multiplier=multiplier,holomorphic=holomorphic,weak=weak,cuspidal=cuspidal,verbose=verbose,**kwds)


    def modular_forms_subspace(self,cuspidal_subspace=False):
        r"""
        Construct the modular forms (i.e. holomorphic and non-weak) subspace of self.
        
        """
        if(not is_int(self.weight()) and is_int(2*self.weight())):
            if(self._from_group):
                G = self._from_group
            else:
                G = self._group
            M=HalfIntegralWeightForms(G,weight=self.weight(),character=self._character,multiplier=self._multiplier,holomorphic=True,weak=False,cuspidal=cuspidal_subspace,dprec=self._dprec,verbose=self._verbose)
        else:
            M=AutomorphicFormspace(G,weight=self.weight(),character=self._character,multiplier=self._multiplier,holomorphic=True,weak=False,cuspidal=cuspidal_subspace,dprec=self._dprec,verbose=self._verbose)
        return M

    def get_element(self,principal_part=None,prec=53,dbase_prec=None,ndig=10,SetC=None,SetY=None,SetM=None,**kwds):
        r"""
        Get an element of the space of HarmonicWeakMaassFormSpace.
        
        INPUT:
        
         - `principal_part`   -- list of principal parts of the form:
                 RR = { '+' : {(j,n) : c^+(j,n)}     # j is a cusp and n>=0 an index
                        '-' : {(j,n) : c^-(j,n)}     # j is a cusp and n<=0 an index
                     }
                corresponding to principal parts (in notation of Bruinier-Funke):
                    \( \Sum_{n>0} c^+(j,n)q^{-n} +  \Sum_{n<0} c^-(j,n)H(n\tau)
        - ''ndig'' -- integer (default 10): the number of requested digits
        - ''dbase_prec'' -- integer (default None): if set, use this number of digits for precision in all mpmath calculations
        - ''SetC'' -- dictionary containing fourier coefficients to keep fixed (and their values)
                      of the form SetC[n][i]=c_i(n)

        """
        pp = principal_part
        if not isinstance(pp,list):
            pp = [pp]
        C = AutomorphicFormSpace._get_element(self,principal_part=principal_part,SetC=SetC,SetY=SetY,SetM=SetM,get_c=True,**kwds)
        if kwds.get('gr',0)<>0:
            return C
        res=list()
        if self._verbose>0:
            print "principal_part=",pp
        if not isinstance(SetC,list):
            numfuns=1
        else:
            numfuns = max(1,len(SetC))
        if len(C)>0:
            for i in range(numfuns):
                F=HarmonicWeakMaassFormElement(self,C[i],prec=prec,principal_part=pp[i])
                res.append(F)
        if len(res)==1:
            res=res[0]
        return res
                
    def set_norm(self,P=None,C=None):
        r"""
        Set normalization for computing Fourier coefficients.
        -''P'' -- principal part = list of dictionaries
        -''C'' -- a list of dictionaries of coefficients in the form SetC[k][(i,n)]=c
              (if only one dictionary is supplied we use the same for all)
        """
        N=dict()
        if isinstance(P,list):
            N['comp_dim']=len(P)
        else:
            N['comp_dim']=1
            P=[P]
        if N['comp_dim']==0:
            raise ValueError,"Need to specify at least one set of conditions!"
        if self._verbose > 0:
            print "comp_dim:=",N['comp_dim']
            print "PP=",P
        N['SetCs']=dict()
        N['SetCs_neg']=dict()
        N['cuspidal']=self._cuspidal
        N['weak']=self._weak        
        nc=self._group.ncusps()
        eps = 2.0*10.0**(-self.prec()*ln(2.0)/ln(10.0))
        if C<>None and len(C)>0:
            C1=dict()
            if isinstance(C,dict):
                C1[0]=[C]
                for i in range(1,N['comp_dim']):
                    C1.append(C)
            elif isinstance(C,list):
                C1 = copy(C)
                if len(C1)<>N['comp_dim']:
                    raise ValueError,"Need the same length of coefficients to set as number of principal parts!" 
            else:
                raise ValueError,"Need the same length of coefficients to set as number of principal parts!"    
        else:
            C1=list()
            for j in range(N['comp_dim']):
                C1.append({})
    
        if self._verbose > 1:
            print "SetC=",C1
        # Set coefficients (e.g. in +-space or similar)
        if C1<>None:
            for j in range(N['comp_dim']):
                N['SetCs'][j]=copy(C1[j]) # dict()
        else:
            for j in range(N['comp_dim']):
                N['SetCs'][j]=dict()
        
        # # Impose further restrictions by cuspidality or non-weakness
        for icusp in range(nc):
            al=self.alpha(icusp)[0] # 
            if self._verbose > 1:
                print "alpha(",icusp,")=",al
            if (al<-eps and not N['weak'])  or (al<=eps and N['cuspidal']):
                for j in range(N['comp_dim']):
                    N['SetCs'][j][(icusp,0)]=0
        ## Set coefficients given by the principal parts
        ## (if applicable, i.e. only the 0th-coefficient)
        for j in range(N['comp_dim']):
            N['SetCs_neg'][j]=dict()
            if(P[j].has_key("-")):
                for (r,n) in P[j]['-']:
                    N['SetCs_neg'][j][(r,n)]=P[j]['-'][(r,n)]
            if P[j].has_key("+"):        
                for (r,n) in P[j]['+']:
                    if n==0: #if H.alpha(r)[0]==0:
                        if not N['SetCs'][j].has_key((r,n)):
                            N['SetCs'][j][(r,n)]=P[j]['+'][(r,n)]
        return N

    

    
class HolomorphicModularForms(AutomorphicFormSpace):
    r"""
    Space of Holomorphic modular forms.

    EXAMPLES:
    
    sage: S=HolomorphicModularForms(Gamma0(1),12,prec=203)
    sage: F=S.get_element(SetM=50)
    sage: F.print_as_q_exp(10,prec=1E-12)
    Principal part:
    f[ 0 ]^+= q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 - 6048*q^6 - 16744*q^7 + 84480*q^8 - 113643*q^9 + O(q^10)


    """
    def __init__(self,G,weight=None,multiplier=None,weak=False,cuspidal=True,dprec=25,verbose=0,**kwds):
        r""" Initialize the space of automorphic forms.
        """
        self._verbose=verbose
        if( hasattr(G,'_is_space_of_automorphic_functions')):
            # Look at the holomorphic space on the group of G with same weight
            # and multiplier as default
            self._group=G._group
            self._from_group=G._from_group
            self._verbose = G._verbose
            if weight==None:
                weight=G.weight()
            if multiplier==None:
                multiplier=G.multiplier()
        else:
            if(isinstance(G,MySubgroup_class)):
                self._group=G
                self._from_group=G._G
            else:
                if(is_int(G)):
                    self._from_group=Gamma0(G)
                    self._group=MySubgroup(self._from_group)
                elif( hasattr(G,'is_subgroup') and G.is_subgroup(SL2Z)):
                    self._group=MySubgroup(G)
                    self._from_group=G
        #AutomorphicFormSpace.__init__(self,GG,weight=weight,multiplier=multiplier,character=character,holomorphic=holomorphic,weak=weak,cuspidal=cuspidal,dprec=dprec,verbose=verbose)
        if hasattr(multiplier,"dimension_cusp_forms") and not weak:
            if weight >2:
                if cuspidal:
                    self._dimension = multiplier.dimension_cusp_forms(weight)
                else:
                    self._dimension = multiplier.dimension_modular_forms(weight)
            else:
                self._dimension=-1
        AutomorphicFormSpace.__init__(self,self._group,weight=weight,multiplier=multiplier,holomorphic=True,weak=weak,cuspidal=cuspidal,dprec=dprec,verbose=verbose,**kwds)
        

        
    def get_element(self,principal_part=None,ndig=10,prec=53,dbase_prec=None,SetC=None,SetY=None,SetM=None,do_mpmath=False,get_mat=False,**kwds):
        if SetC==None:
            if self.is_cuspidal():
                SetC=[{(0,0):0,(0,1):1}]
            else:
                SetC=[{(0,0):1,(0,1):0}]
        if principal_part<>None:
            pp = principal_part            
        else:
            pp=list()
            for i in range(len(SetC)):
                pp.append({'+':{},'-':{}})
        if not isinstance(pp,list):
            pp0=[pp]
        else:
            pp0=pp
        if self._verbose>0:
            print "Using pp=",pp
            print "SetC=",SetC
            print "SetY=",SetY
            print "SetM=",SetM
        C = AutomorphicFormSpace._get_element(self,principal_part=pp,SetC=SetC,SetY=SetY,SetM=SetM,do_mpmath=do_mpmath,get_mat=get_mat,get_c=True,**kwds)
        if kwds.get('gr',0)==1:
            return C
        res=list()
        if len(C)>0:
            for i in range(len(C)):
                F= HolomorphicModularFormElement(self,C[i],prec=prec,principal_part=pp0[i])
                res.append(F)
        if len(res)>1:
            return res
        else:
            return res[0]

        
class HolomorphicModularFormElement(AutomorphicFormElement):    
    r"""
    Create an Harmonic Weak Maass form.
    """ 

    def __init__(self,M,C=None,prec=None,principal_part=None):
        r"""
        See ``HarmonicWeakMaassFormElement`` for full documentation.

        Initialize an harmonic weak Maass form.
        
        INPUT:

        
            -''M'' -- Space of automorphic forms
            
            -''k'' -- Weight.
            -''C''-- Fourier coefficients
            -''prec'' -- integer, precision (if given by construction, default None)

        EXAMPLES::

        """
        #print "typeM=",type(M)
        self._class_name ="HolomorphicModularFormElement"
        AutomorphicFormElement.__init__(self,M,C,prec=prec,principal_part=principal_part)

    def list_coefficients(self,nmax,norm=False,plus=True,minus=True,cusp=None):
        r"""
        List coefficients of self.
        """
        for r in self._coeffs.keys():
            print ""
            for j in self._coeffs[r].keys():
                if cusp<>None and j<>cusp:
                    continue
                print ""
                l=self._coeffs[r][j].keys()
                l.sort(cmp=self.my_zzcmp)
                for n in l:
                    al = self._space.alpha(j)[0]
                    if al<>0:
                        nal = ZZ(n) + al
                        st_nal=sci_pretty_print(nal,3)
                    else:
                        nal=str(n)
                    c=self._coeffs[r][j][n]
                    print "C[",j,",",nal,"]=",c


### Construction routines for specific types of forms
def WeakModularForm(G,weight=0,principal_part="q^-1",**kwds):
    r"""
    Construct a Weaklyp holomorphic modular form on G.

    INPUT:
    - `principal_part` -- Principal part given as a (vector) of q-series with coefficients in QQ.
    - `weight` -- weight
    - `kwds` -- keywords to get passes on to the HolomorphicModulaForms space

    NOTE: If a vector of principal series is given the entries correspond to the cusps in G (with the same order)

    """
    M = HolomorphicModularForms(G,weight=weight,weak=True,**kwds)
    pp = extract_princial_part(M,principal_part)
    F = M.get_element(principal_part={'+':pp})
    return F


def HarmonicWeakMaassForm(G,weight=0,principal_part="q^-1",verbose=0,**kwds):
    r"""
    Construct a Harmonic weak Maass form form on G.

    INPUT:
    - `principal_part` -- Principal part given as a (vector) of q-series with coefficients in QQ.
    - `weight` -- weight
    - `kwds` -- keywords to get passes on to the HolomorphicModulaForms space

    NOTE: If a vector of principal series is given the entries correspond to the cusps in G (with the same order)

    NOTE: can not specify 

    EXAMPLE:
    # To compute E2 -- the non-holomorphic weight 2 Eisenstein series for SL(2,Z)
    sage: pp={'-':{(0,0):RR(-3)/RR.pi()},'+':{(0,0):RR(1)}}
    sage: F=psage.modform.maass.automorphic_forms.HarmonicWeakMaassForm(Gamma0(1),weight=2,principal_part=pp,almost_holomorphic=True,SetM=10,SetY=0.4)
    sage: F.coeffs()[0][0]
{0: 1.00000000000000,
 1: -24.0000000000000,
 2: -72.0000000000000,
 3: -95.9999999999992,
 4: -168.000000000001,
 5: -143.999999999936,
 6: -288.000000000510,
 7: -191.999999986893,
 8: -360.000000118637,
 9: -311.999997918382,
 10: -431.999999938662}
    """
    #print "kwds=",kwds
    M = extract_hwmf_space(G,weight,**kwds)
    #print "kwds=",kwds
    pp = extract_princial_part(M,principal_part)
    M._verbose = verbose
    if pp.get('-',{})<>{}:
        M._almost_holomorphic = True
    #M = HarmonicWeakMaassFormSpace(G,weight=weight,weak=True,verbose=verbose)
    
    F = M.get_element(principal_part=pp,**kwds)
    return F

def extract_hwmf_space(X,weight=0,**kwds):
    multiplier = extract_multiplier(X,weight,**kwds)
    if kwds.get('vv',0)==1:
        return VVHarmonicWeakMaassFormSpace(G=multiplier.group(),weight=weight,multiplier=multiplier,**kwds)
    else:
        character=kwds.pop('character',0)
        holomorphic=kwds.pop('holomorphic',False)
        weak=kwds.pop('weak',True)
        almost_holomorphic=kwds.pop('almost_holomorphic',False)
        cuspidal=kwds.pop('cuspidal',False)
        unitary_action=kwds.pop('unitary_action',0)
        dprec=kwds.pop('dprec',15)
        prec=kwds.pop('prec',53)
        return HarmonicWeakMaassFormSpace(G=multiplier.group(),weight=weight,multiplier=multiplier,character=character,holomorphic=holomorphic,weak=weak,almost_holomorphic=almost_holomorphic,cuspidal=cuspidal,unitary_action=unitary_action,dprec=dprec,prec=prec,**kwds)

def extract_multiplier(X,weight=0,**kwds):
    vv = 0    
    if kwds.get('vector-valued',0)==1 or kwds.get('vv',0)==1:
        vv = 1
    if isinstance(X,MultiplierSystem):
        rank = X.rank()
        multiplier = X
    else: 
        ## Try to construct an appropriate multiplier
        try:
            if isinstance(X,WeilModule):
                multiplier = WeilRepMultiplier(X,weight=weight)
            elif vv == 0:
                if isinstance(X,(int,Integer)):
                    G = MySubgroup(Gamma0(X))                
                elif not isinstance(X,MySubgroup_class):
                    G = MySubgroup(X)
                if is_int(weight):
                    multiplier = TrivialMultiplier(X,weight=weight)
                elif is_int(2*weight):
                    multiplier = ThetaMultiplier(G,weight=weight)
            else:
                multiplier = WeilRepMultiplier(X,weight=weight) ## Weil rep. of Z, x-> X*x^2
        except ArithmeticError:
            raise ValueError,"Could not construct space from {0}".format(X)
    return multiplier

            
def extract_princial_part(M,principal_part):
    r"""
    Compute the principal part in the format we want for the constructors.
    """### Interpret the principal part .
    LP = LaurentPolynomialRing(QQ,name='q')
    q = LP.gens()[0]
    YP = LaurentPolynomialRing(QQ,name='y')
    y = YP.gens()[0]
    ppdict={'-':{},'+':{}}
    ## We start by setting the constant terms to zero
    ## (unless they are explicitly overriden in the given principal part)
    for j in range(M.group().ncusps()):
        if M.alpha(j)[0]<=0:
            ppdict['+'][(j,0)]=0
            ppdict['-'][(j,0)]=0
    if isinstance(principal_part,dict):
        for r,k in principal_part.get('+',{}):
            ppdict['+'][(r,k)]=principal_part['+'][(r,k)]
        ppdict['-']=principal_part.get('-',{})
        #print "dict:",ppdict
        return ppdict
    if not isinstance(principal_part,list):
        principal_part=[principal_part]
    for j in range(len(principal_part)):
        pp = principal_part[j]        
        if isinstance(pp,str):
            print "pp0=",pp
            # # We need to convert to a Laurent polynomial.            
            pp = eval(pp)
        elif not hasattr(pp,"add_bigoh") and pp.base_ring()==QQ:
            raise ValueError,"Did not get principal part of correct type! Got:{0}".format(principal_part)
        print "pp=",pp
        ck=0
        for k in pp.exponents():
            ppdict['+'][(j,k[0])]=pp.coefficients()[ck]
            ck+=1
    #print "ppdict=",ppdict
    return ppdict
    

    
def shifts(FQM):
    v=dict()
    for x in FQM.list():
        v[x]=FQM.Q(x)
    return v


def set_norm_harmonic_weak_maass_forms(H,P=None,C=None):
    r"""
    -''H'' -- space of automorphic forms
    -''P'' -- principal parts = list of dictionaries
    -''C'' -- a list of dictionaries of coefficients in the form SetC[k][(i,n)]=c
              (if only one dictionary is supplied we use the same for all)
    """
    N=dict()
    N['comp_dim']=len(P)
    if(N['comp_dim']==0):
        raise ValueError,"Need to specify at least one set of conditions!"
    if H._verbose > 0:
        print "comp_dim:=",N['comp_dim']
        print "PP=",P
    N['SetCs']=dict()
    N['SetCs_neg']=dict()
    N['cuspidal']=H._cuspidal
    N['weak']=H._weak        
    nc=H._group.ncusps()
    eps = 2.0*10.0**(-H.prec()*ln(2.0)/ln(10.0))
    if C<>None and len(C)>0:
        C1=dict()
        if(isinstance(C,dict)):
            C1[0]=[C]
            for i in range(1,N['comp_dim']):
                C1.append(C)
        elif(isinstance(C,list)):
            C1 = copy(C)
            if(len(C1)<>N['comp_dim']):
                raise ValueError,"Need the same length of coefficients to set as number of principal parts!" 
        else:
            raise ValueError,"Need the same length of coefficients to set as number of principal parts!"    
    else:
        C1=list()
        for j in range(N['comp_dim']):
            C1.append({})
    
    if H._verbose > 1:
        print "SetC=",C1
    # Set coefficients (e.g. in +-space or similar)
    if(C1<>None):
        for j in range(N['comp_dim']):
            N['SetCs'][j]=copy(C1[j]) # dict()
    else:
        for j in range(N['comp_dim']):
            N['SetCs'][j]=dict()
    ## Set coefficients given by the principal parts
    for j in range(N['comp_dim']):
        N['SetCs_neg'][j]=dict()
        if(P[j].has_key("-")):
            for (r,n) in P[j]['-']:
                N['SetCs_neg'][j][(r,n)]=P[j]['-'][(r,n)]
        if P[j].has_key("+"):        
            for (r,n) in P[j]['+']:
                if H.alpha(r)[0]==0:
                    N['SetCs'][j][(r,n)]=P[j]['+'][(r,n)]
    ## Impose further restrictions by cuspidality or non-weakness
    for icusp in range(nc):
        al=H.alpha(icusp)[0]
        if H._verbose > 1:
            print "alpha(",icusp,")=",al
        if (al<-eps and not N['weak'])  or (al<=eps and N['cuspidal']):
            for j in range(N['comp_dim']):
                N['SetCs'][j][(icusp,0)]=0
        if al<=eps:
            for j in range(N['comp_dim']):
                if(P[j]['+'].has_key((icusp,0))):
                    ## We only have to set the zeroth coefficients if we don't use it as a variable for a^{-}
                    if(P[j].has_key('-') and P[j]['-'].has_key((icusp,0))):
                        N['SetCs'][j][(icusp,0)]=P[j]['+'][(icusp,0)]
                #else:
                #    N['SetCs'][j][(icusp,0)]=0
    # remove coefficients that are already set using the principal part
    # and set them using the supplied coefficients otherwise
    #for i in range(len(PP)):
    #    for (r,n) in PP[i].keys():
    #        N['SetCs'][i][(r,n)]=PP[i]
    #print "N end=",N
    return N



def set_norm_harmonic_weak_maass_forms2(H,P,C=None):
    r"""
    -''H'' -- space of automorphic forms
    -''P'' -- principal parts
    -''C'' -- a dictionary of set coefficients in the form SetC[(i,n)]=c
    """
    N=dict()
    N['comp_dim']=1            
    N['SetCs']=dict()
    N['cuspidal']=H._cuspidal
    N['weak']=H._weak        
    nc=H._group.ncusps()
    eps = 2.0*10.0**(-H.prec()*ln(2.0)/ln(10.0))
    if(C<>None and len(C.keys())>0):
        if(is_int(C.keys()[0])):
            SetC=C
        else:
            SetC=dict()
            SetC[0]=C
    else:
        SetC=C
    if H._verbose > 1:
        print "SetC=",SetC
    if(C<>None and len(SetC.keys())>0):
        N['comp_dim']=len(SetC.keys())
    #print "N=",N
    #print "SetC.keys=",SetC.keys()
    for j in range(N['comp_dim']):
        N['SetCs'][j]=dict()
    if(P.has_key((0,0)) and H._holo):
        #print "holo"
        for j in range(N['comp_dim']):
            N['SetCs'][j][(0,0)]=0            
    if(N['cuspidal']):
        for icusp in range(nc):
            v=H.alpha(icusp)[1]
            if(v==1):
                for j in range(N['comp_dim']):
                    N['SetCs'][j][(icusp,0)]=0
    if(not N['weak']):
        for icusp in range(nc):
            al=H.alpha(icusp)[0]
            if al<-eps:
                for j in range(comp_dim):
                    SetCs[j][(icusp,0)]=0
    #print "N=",N
    if(SetC<>None):
        for i in SetC.keys():
            for (r,n) in SetC[i].keys():
                if(P.has_key((r,n))):
                    N['SetCs'][i][(r,n)]=P[(r,n)]
                else:
                    N['SetCs'][i][(r,n)]=SetC[i][(r,n)]
        
    #print "N=",N
    return N



def solve_system_for_harmonic_weak_Maass_waveforms(W,N,gr=0):
    r"""
    Solve the linear system to obtain the Fourier coefficients of Maass forms

    INPUT:

    - ``W`` --   (system) dictionary
        - ``W['Ms']``  -- M start
        - ``W['Mf']``  -- M stop
        - ``W['nc']``  -- number of cusps
        - ``W['space']``  -- space of automorphic forms
        - ``W['V']``   -- matrix of size ((Ms-Mf+1)*nc)**2
        - ``W['RHS']`` -- right hand side (for inhomogeneous system) matrix of size ((Ms-Mf+1)*nc)*(dim)
    - ``N`` -- normalisation (dictionary, output from the set_norm_for_maass function)
        - ``N['SetCs']``   -- Which coefficients are set and their values
        - ``N['comp_dim']``-- How large is the assumed dimension of the solution space
        - ``N['num_set']`` -- Number of coefficients which are set
        

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
    nc=W.get('nc',1)
    PP=W.get('PP',[])
    H = W.get('space',None)
    if not H:
        raise TypeError,"Need a space together with our W!"
    verbose = H._verbose
    #alphas=W['alphas']
    alphas = H.alphas()
    Ml=W['Ml'] #Mf-Ms+1
    variable_a_plus=W['var_a+']
    variable_a_minus=W['var_a-']
    if(V.ncols()<>Ml*nc or V.nrows()<>Ml*nc):
        raise Exception," Wrong dimension of input matrix!"
    # we have to assume that all normalizations use the same coefficients
    maxit=1000
    SetCs=N['SetCs']
    SetCs_neg=N.get('SetCs_neg',{})
    CF = MPComplexField(H.prec())
    zero = CF(0)
    comp_dim=N['comp_dim']
    use_sym=0
    SetClist=dict()
    for j in range(0,comp_dim):
        SetClist[j]=dict()
    if len(PP)>0 and ((comp_dim<>len(SetCs.keys()) and comp_dim<>len(PP))):
        print "comp_dim=",comp_dim
        print "SetC=",SetCs
        print "PP=",PP
        raise ValueError," Inconsistent normalization SetCs:%s" % SetCs
    num_set=0
    for j in range(0,comp_dim):
        # # First we treat set values of coefficients not corresponsing to the principal part
        for (r,n) in SetCs[j].keys():
            nr = r*Ml+n
            if nr>=0 or not H.is_holomorphic():
                SetClist[j][nr]=SetCs[j][(r,n)]
            elif PP[j]['+'].has_key((r,n)) and PP[j]['-'].has_key((r,n)):
                SetClist[j][nr]=0
        if verbose>0:
            print "SetClist_pos=",SetClist
            print "var_a+=",variable_a_plus[j]
            print "var_a-=",variable_a_minus[j]
        ## Then we check the zeroth coefficients
        for r in range(nc):
            if(alphas[r][1]==1):
                if( (not variable_a_plus[j][r]) and (not variable_a_minus[j][r])):
                    nr = r*Ml
                    if(SetCs_neg.get(j,{}).has_key((r,0))):
                        SetClist[j][nr]=CF(SetCs_neg[j][(r,0)]) 
        num_set=len(SetClist[0].keys())
    if verbose>0:
        print "SetClist_tot=",SetClist
    t=V[0,0]
    if(isinstance(t,float)):
        mpmath_ctx=mpmath.fp
    else:  
        mpmath_ctx=mpmath.mp
    if verbose>0:
        print "mpmath_ctx=",mpmath_ctx
    #use_symmetry=False
    MS = MatrixSpace(CF,int(Ml*nc-num_set),int(comp_dim))
    RHS = Matrix_complex_dense(MS,0,True,True)
    # We allow for either a variation of principal parts or of set coefficients
    if(W.has_key('RHS')):
        l=W['RHS'].ncols()
        if(l>1 and l<>comp_dim):
            raise ValueError,"Incorrect number of right hand sides!"
        
    MS2 = MatrixSpace(CF,int(Ml*nc-num_set),int(Ml*nc-num_set))
    LHS = Matrix_complex_dense(MS2,0,True,True)
    #LHS=mpmath_ctx.matrix(int(Ml*nc-num_set),int(Ml*nc-num_set))
    roffs=0
    if verbose>0:
        print "Ml=",Ml
        print "num_set=",num_set
        print "SetCs=",SetCs
        print "SetClist=",SetClist
        #print "Valslist=",Valslist
        print "V.rows=",V.nrows()
        print "V.cols=",V.ncols()
        print "LHS.rows=",LHS.nrows()
        print "LHS.cols=",LHS.ncols()
        print "RHS.rows=",RHS.nrows()
        print "RHS.cols=",RHS.ncols()
        print "use_sym=",use_sym
    for r in range(V.nrows()):
        cr=r+Ms
        if(SetClist[0].keys().count(r+Ms)>0):
            roffs=roffs+1
            continue
        for fn_j in range(comp_dim):
            if(W.has_key('RHS') and W['RHS'].ncols()>fn_j):
                RHS[r-roffs,fn_j]=CF(-W['RHS'][r,fn_j])
            elif(W.has_key('RHS')):
                RHS[r-roffs,fn_j]=CF(-W['RHS'][r,0])
            else:
                RHS[r-roffs,fn_j]=zero
            for c in SetClist[fn_j].keys():
                v=CF(SetClist[fn_j][c])
                tmp=v*V[r,c-Ms]
                RHS[r-roffs,fn_j]=RHS[r-roffs,fn_j]-tmp
        coffs=0
        for k in range(V.ncols()):            
            if(SetClist[0].keys().count(k+Ms)>0):
                coffs=coffs+1
                continue
            try:                
                LHS[r-roffs,k-coffs]=V[r,k]
            except IndexError:
                print "r,k=",r,k
                print "V.rows=",V.nrows()
                print "V.cols=",V.ncols()
                print "roffs,coffs=",roffs,coffs
                print "r-roffs,k-coffs=",r-roffs,k-coffs
                print "LHS.rows=",LHS.nrows()
                print "LHS.cols=",LHS.ncols()                
                raise IndexError,"Matrix / coefficients is set up wrong!"
            #print "LHS[",r,k,"]=",LHS[r-roffs,k-coffs]
    if gr==1:
        return LHS,RHS
    smin=smallest_inf_norm(LHS)
    if verbose>0:
        print "sminfn=",smin
    dps0=CF.prec()
    done=False
    i=1
    while (not done and i<=maxit):
        try:
            Q,R = LHS.qr_decomposition()
            #A, p = mpmath_ctx.LU_decomp(LHS)
            done=True
        except ZeroDivisionError:
            #t=int(mpmath_ctx.ceil(-mpmath_ctx.log10(smallest_inf_norm(LHS))))
            t=int(ceil(-log_b(smallest_inf_norm(LHS),10)))
            dps=t+5*i; i=i+1
            if verbose>-1:
                print "raising number of digits to:",dps
            LHS.set_prec(dps)
            # raise ZeroDivisionError,"Need higher precision! Use > %s digits!" % t
    if(i>=maxit):
        raise ZeroDivisionError,"Can not raise precision enough to solve system! Should need > %s digits! and %s digits was not enough!" % (t,dps)
    X=dict()
    for fn_j in range(comp_dim):
        X[fn_j] = dict() 
        X[fn_j][0] = dict() 
        v = RHS.column(fn_j)
        if verbose>0:
            print "len(B)=",len(v)
            #print "RHS=",v
        #b = mpmath_ctx.L_solve(A, RHS.column(fn_j), p)
        TMP = LHS.solve(v) #mpmath_ctx.U_solve(A, b)
        roffs=0
        res = (LHS*TMP-v).norm()
        if verbose>0:
            print "res(",fn_j,")=",res
        #res = mpmath_ctx.norm(mpmath_ctx.residual(LHS, TMP, RHS.column(fn_j)))
        #print "res(",fn_j,")=",res
        for i in range(0,nc):
            X[fn_j][0][i]=dict()
        for i in range(nc):
            roffs2=0
            for n in range(Ml):
                nn=i*Ml+n+Ms
                key=n+Ms
                #if(i==1):
                #    print n,key
                if(SetClist[fn_j].keys().count(nn)>0):
                    if verbose>1:
                        print "We have set ",nn
                    roffs=roffs+1
                    X[fn_j][0][i][key]=SetClist[fn_j][nn] 
                    if verbose>0:
                        print "X[",fn_j,",",i,",",key,"]=",SetClist[fn_j][nn]
                        print "nn=",nn
                    continue
                try:
                    #X[fn_j][0][i][n-roffs2+Ms]=TMP[nn-Ms-roffs,0]
                    X[fn_j][0][i][key]=TMP[nn-Ms-roffs]
                except IndexError:
                    print "n*Mli-roffs=",n,'+',Ml,'*',i,'-',roffs,"=",n+Ml*i-roffs
                ## We also insert the principal part if it is applicable
    #mpmath.mp.dps=dpold
    # return x
    return X

def smallest_inf_norm(V):
    r"""
    Computes the smallest of the supremum norms of the columns of a matrix.

    INPUT:

        - ``V`` -- matrix (real/complex)

    OUTPUT:

        - ``t`` -- minimum of supremum norms of the columns of V

    EXAMPLE::


        sage: A=mpmath.matrix([['0.1','0.3','1.0'],['7.1','5.5','4.8'],['3.2','4.4','5.6']])
        sage: smallest_inf_norm(A)
        mpf('5.5')

    
    """
    minc=100
    mi=0
    try:
        nc = V.ncols(); nr=V.nrows()
    except AttributeError:
        nc = V.cols; nr=V.rows
        
    for j in range(nc):
        maxr=0
        for k in range(nr):
            t=abs(V[k,j])
            if(t>maxr):
                maxr=t
        if(maxr<minc):
            minc=maxr
            mi=j
    return minc


def smallest_inf_norm_mpmath(V):
    r"""
    Computes the smallest of the supremum norms of the columns of a matrix.

    INPUT:

        - ``V`` -- matrix (real/complex)

    OUTPUT:

        - ``t`` -- minimum of supremum norms of the columns of V

    EXAMPLE::


        sage: A=mpmath.matrix([['0.1','0.3','1.0'],['7.1','5.5','4.8'],['3.2','4.4','5.6']])
        sage: smallest_inf_norm(A)
        mpf('5.5')

    
    """
    minc=100
    mi=0
    for j in range(V.cols):
        maxr=0
        for k in range(V.rows):
            t=abs(V[k,j])
            if(t>maxr):
                maxr=t
        if(maxr<minc):
            minc=maxr
            mi=j
    return minc



def solve_system_for_harmonic_weak_Maass_waveforms_mpmath(W,N):
    r"""
    Solve the linear system to obtain the Fourier coefficients of Maass forms

    INPUT:

    - ``W`` --   (system) dictionary
        - ``W['Ms']``  -- M start
        - ``W['Mf']``  -- M stop
        - ``W['nc']``  -- number of cusps
        - ``W['V']``   -- matrix of size ((Ms-Mf+1)*nc)**2
        - ``W['space']``   -- space of automorphic forms
        - ``W['RHS']`` -- right hand side (for inhomogeneous system) matrix of size ((Ms-Mf+1)*nc)*(dim)
    - ``N`` -- normalisation (dictionary, output from the set_norm_for_maass function)
        - ``N['SetCs']``   -- Which coefficients are set and their values
        - ``N['comp_dim']``-- How large is the assumed dimension of the solution space
        - ``N['num_set']`` -- Number of coefficients which are set
        
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
    nc=W['nc']
    PP=W.get('PP',[])
    alphas=W['alphas']
    H=W.get('space',None)
    if not H:
        raise ValueError," Need a space in W!"
    verbose =H._verbose
    Ml=W['Ml'] #Mf-Ms+1
    variable_a_plus=W['var_a+']
    variable_a_minus=W['var_a-']
    if(V.cols<>Ml*nc or V.rows<>Ml*nc):
        raise Exception," Wrong dimension of input matrix!"
    # we have to assume that all normalizations use the same coefficients
    SetCs=N['SetCs']
    SetCs_neg=N['SetCs_neg']
    zero=mpmath.mp.mpf(0)
    #Vals=N['Vals']
    ### We have to determine whether we have two "constant" terms
    #if(Ms<0 and Ms+Mf<>0):
    #two_terms=1
    #else:
    #    two_terms=0
    #print "Using two constant terms!:",two_terms
    comp_dim=N['comp_dim']
    use_sym=0
    #if(len(SetCs.keys())>0):
    #    num_set=len(SetCs[0])
    #else:
    #    num_set=0
    #if(len(SetCs_neg.keys())>0):
    #    num_set=num_set+len(SetCs_neg[0])
    ## converse the dictionary to a list 
    SetClist=dict()
    for j in range(0,comp_dim):
        SetClist[j]=dict()
    if(comp_dim<>len(SetCs.keys()) and comp_dim<>len(PP)):
        print "comp_dim=",comp_dim
        print "SetC=",SetCs
        print "PP=",PP
        raise ValueError," Inconsistent normalization SetCs:%s" % SetCs
    num_set=0
    for j in range(0,comp_dim):
        # # First we treat set values of coefficients not corresponsing to the principal part
        for (r,n) in SetCs[j].keys():
            for j in range(comp_dim):
                #if(two_terms and n>=0):
                #    nr = r*Ml+n+1
                #else:
                nr = r*Ml+n
                SetClist[j][nr]=SetCs[j][(r,n)]
        if verbose>0:
            print "SetClist_pos=",SetClist
        #if(not two_terms):
        #continue
        #for (r,n) in SetCs_neg[j].keys():
        #    for j in range(comp_dim):
        #        nr = r*Ml+n
        #        SetClist[j][nr]=SetCs_neg[j][(r,n)]
        ## Then we check the zeroth coefficients
        for r in range(nc):
            if(alphas[r][1]==1):
                if( (not variable_a_plus[r]) and (not variable_a_minus[r])):
                    nr = r*Ml
                    #if(SetCs_neg[j].has_key((r,0))):
                    #    
                    #SetClist[j][nr]=zero
        num_set=len(SetClist[0].keys())
    if verbose>0:
        print "SetClist_tot=",SetClist
    t=V[0,0]
    if(isinstance(t,float)):
        mpmath_ctx=mpmath.fp
    else:  
        mpmath_ctx=mpmath.mp
    if verbose>0:
        print "mpmath_ctx=",mpmath_ctx
    #use_symmetry=False
    RHS=mpmath_ctx.matrix(int(Ml*nc-num_set),int(comp_dim))
    # We allow for either a variation of principal parts or of set coefficients
    #
    if(W.has_key('RHS')):
        l=W['RHS'].cols
        if(l>1 and l<>comp_dim):
            raise ValueError,"Incorrect number of right hand sides!"

    LHS=mpmath_ctx.matrix(int(Ml*nc-num_set),int(Ml*nc-num_set))
    roffs=0
    if verbose>0:
        print "Ml=",Ml
        print "num_set=",num_set
        print "SetCs=",SetCs
        print "SetClist=",SetClist
        #print "Valslist=",Valslist
        print "V.rows=",V.rows
        print "V.cols=",V.cols
        print "LHS.rows=",LHS.rows
        print "LHS.cols=",LHS.cols
        print "RHS.rows=",RHS.rows
        print "RHS.cols=",RHS.cols
        print "use_sym=",use_sym
    for r in range(V.rows):
        cr=r+Ms
        if(SetClist[0].keys().count(r+Ms)>0):
            roffs=roffs+1
            continue
        for fn_j in range(comp_dim):
            if(W.has_key('RHS') and W['RHS'].cols>fn_j):
                RHS[r-roffs,fn_j]=-W['RHS'][r,fn_j]
            elif(W.has_key('RHS')):
                RHS[r-roffs,fn_j]=-W['RHS'][r,0]
            for c in SetClist[fn_j].keys():
                v=SetClist[fn_j][c]
                if(mpmath_ctx==mpmath.mp):
                    tmp=mpmath_ctx.mpmathify(v)
                elif(isinstance(v,float)):
                    tmp=mpmath_ctx.mpf(v)
                else:
                    tmp=mpmath_ctx.mpc(v)
                tmp=tmp*V[r,c-Ms]
                RHS[r-roffs,fn_j]=RHS[r-roffs,fn_j]-tmp
        coffs=0
        for k in range(V.cols):            
            if(SetClist[0].keys().count(k+Ms)>0):
                coffs=coffs+1
                continue
            try:                
                LHS[r-roffs,k-coffs]=V[r,k]
            except IndexError:
                print "r,k=",r,k
                print "V.rows=",V.rows
                print "V.cols=",V.cols
                print "roffs,coffs=",roffs,coffs
                print "r-roffs,k-coffs=",r-roffs,k-coffs
                print "LHS.rows=",LHS.rows
                print "LHS.cols=",LHS.cols                
                return
            #print "LHS[",r,k,"]=",LHS[r-roffs,k-coffs]
    #return LHS
    smin=smallest_inf_norm_mpmath(LHS)
    if verbose>0:
        print "sminfn=",smin
    if(smin<>0):
        t=int(mpmath_ctx.ceil(-mpmath_ctx.log10(smin)))
    else:
        raise ValueError,"Something wrong with normalization. Got min norm=0!"
    dpold=mpmath.mp.dps
    mpmath.mp.dps=max(dpold,t+5)
    maxit=100;i=0
    done=False
    if verbose>0:
        print "using number of digits:",mpmath.mp.dps
    while (not done and i<=maxit):
        try:
            A, p = mpmath_ctx.LU_decomp(LHS)
            done=True
        except ZeroDivisionError:
            mpmath.mp.dps=mpmath.mp.dps+5*i; i=i+1
            print "raising number of digits to:",mpmath.mp.dps
            # raise ZeroDivisionError,"Need higher precision! Use > %s digits!" % t
    if(i>=maxit):
        raise ZeroDivisionError,"Can not raise precision enough to solve system! Should need > %s digits! and %s digits was not enough!" % (t,mpmath.mp.dps)

    #try:
    #    A, p = mpmath_ctx.LU_decomp(LHS)
    #except ZeroDivisionError:
    #    t=int(mpmath_ctx.ceil(-mpmath_ctx.log10(smallest_inf_norm(LHS))))
    #    raise ZeroDivisionError,"Need higher precision! Use > %s digits!" % t
    #return A
    #for k in range(RHS.rows):
    #    print "RHS(",k,")=",RHS[k,0]
    X=dict()
    for fn_j in range(comp_dim):
        X[fn_j] = dict() #mpmath.matrix(int(Ml),int(1))
        X[fn_j][0] = dict() #mpmath.matrix(int(Ml),int(1))
        if verbose>0:
            print "len(B)=",len(RHS.column(fn_j))
        b = mpmath_ctx.L_solve(A, RHS.column(fn_j), p)
        TMP = mpmath_ctx.U_solve(A, b)
        roffs=0
        res = mpmath_ctx.norm(mpmath_ctx.residual(LHS, TMP, RHS.column(fn_j)))
        #print "res(",fn_j,")=",res
        for i in range(0,nc):
            X[fn_j][0][i]=dict()
        #for n in range(Ml):
        #    if(SetClist[fn_j].keys().count(i*Ml+n+Ms)>0):
        #        roffs=roffs+1
        #        #print "X[",fn_j,",",n,",Vals[fn_j][n]
        #        X[fn_j][0][n+Ms]=SetClist[fn_j][i*Ml+n+Ms]
        #        continue
        #    X[fn_j][0][n+Ms]=TMP[n-roffs,0]

        for i in range(nc):
            roffs2=0
            for n in range(Ml):
                nn=i*Ml+n+Ms
                ## get the key of the coefficient we set
                #if(two_terms):
                #    if(n+Ms==0):                        
                #        key='-0'  ## the constant 'negative' term
                #    elif(n+Ms>=1):
                #        key=n-1+Ms  ## the constant 'positive' term
                #    else:
                #        key=n+Ms
                #else:
                key=n+Ms
                #if(i==1):
                #    print n,key
                if(SetClist[fn_j].keys().count(nn)>0):
                    if verbose>1:
                        print "We have set ",nn
                    roffs=roffs+1
                    X[fn_j][0][i][key]=SetClist[fn_j][nn] 
                    if verbose>1:
                        print "X[",fn_j,",",i,",",key,"]=",SetClist[fn_j][nn]
                        print "nn=",nn
                    continue
                #if(two_terms and n+Ms==1):
                #    X[fn_j][0][i]['-0']=TMP[nn-Ms-roffs,0]
                #if(two_terms and n+Ms==1):
                #    roffs2=roffs2+1
                try:
                    #X[fn_j][0][i][n-roffs2+Ms]=TMP[nn-Ms-roffs,0]
                    X[fn_j][0][i][key]=TMP[nn-Ms-roffs,0]
                except IndexError:
                    print "n*Mli-roffs=",n,'+',Ml,'*',i,'-',roffs,"=",n+Ml*i-roffs
                ## We also insert the principal part if it is applicable
    mpmath.mp.dps=dpold
    # return x
    return X

def is_int(q):
    r"""
    Find out if the rational number q is an integer.

    INPUT:
    -''q'' -- integer/rational/real

    OUTPUT:
    - logical -- True if q is an integer otherwise False

    EXAMPLES::

        sage: is_int(1)
        True
        sage: is_int(float(1.0))
        True
        sage: is_int(RR(1.0))   
        True
        sage: is_int(6/3) 
        True
        sage: is_int(6/4)
        False
        sage: is_int(1.5)    
        False
        sage: is_int(Gamma0(1))    
        False
        
    """
    if(isinstance(q,sage.rings.integer.Integer) or isinstance(q,int)):
        return True
    if(isinstance(q,sage.rings.rational.Rational)):
        n=q.denominator()
        if(n==1):
            return True
    if(isinstance(q,tuple)):
        return False
    try:
        if(floor(q)==ceil(q)):
            return True
    except:
        pass
    return False


def rational_approximation(x,eps=1E-12):
    r""" Computes an approximation to x in QQ(i).

    INPUT:
    -''x'' -- real or complex number
    -''eps'' -- desired precisionin the approximation
    OUTPUT:
    -''y'' -- rational approximation with error eps
    If |x-[x]|<eps we return [x], the nearest integer to x otherwise |x-ncf(x)| < eps where the ncf(x) is a continued fraction approximation.

    """
    if hasattr(x,"imag"):
        if hasattr(x,"ae"):
            xr = x.real; xi=x.imag 
        else:
            #    print "x=",x,type(x)
            xr = real(x); xi=imag(x)
        if xi<>0:
            xr=rational_approximation(xr,eps)
            xi=rational_approximation(xi,eps)
            return xr+I*xi
        x = xr
    if isinstance(x,mpf):
        x = RealField(mpmath.mp.prec)(x)
        #return rational_approximation(x,eps)
    if isinstance(x,(int,Integer,ZZ)):
        return QQ(x)
    if hasattr(x,"parent"):
        if x.parent() == QQ:
            return x
    if not isinstance(x,(float,type(RR(1)))):
        raise ValueError,"Can not find rational approximation to x:{0} (of type {1}) with eps={2}".format(x,type(x),eps)
    n = nearest_integer(x)
    prec=x.parent().prec()
    dprec=ceil(eps*ln(2.0)/ln(10.0))
    RF = RealField(prec)
    if abs(x-n)<eps:
        return QQ(n)
    else:
        # xr is in [0,1]

        xr = RealField(prec)(x)-RealField(prec)(n)
        # Note: ncf[0]=0
        ncf = nearest_integer_continued_fraction(xr,dprec)
        ## Get approximation from sequence: 
        #print "ncf=",ncf
        for j in range(len(ncf)):
            y = real_from_nearest_integer_continued_fraction(ncf[0:j+1])
            err=RR(abs(y-xr))
            #print "y=",y
            #print "err=",err,type(err)
            #print "eps=",eps,type(eps)
           #print "cmp:",(err<eps)
            if err<RR(eps):
                #print "Do quit!"
                break
        return QQ(y+n)


def real_from_nearest_integer_continued_fraction(ncf):
    r"""
    Take a list of coeffficients in a nearest integer cotinued fraction and return the corresponding real point.
    """
    y = 0
    for j in range(len(ncf),1,-1):
        y = QQ(-1/(y+ncf[j-1]))
    return QQ(y+ncf[0])


## needed for pickling
import __main__
__main__.AutomorphicFormSpace=AutomorphicFormSpace
__main__.HarmonicWeakMaassFormSpace=HarmonicWeakMaassFormSpace
__main__.AutomorphicFormElement=AutomorphicFormElement
__main__.HalfIntegralWeightForms=HalfIntegralWeightForms
#__main__.
#__main__.


###  temporary routines for testing purposes
def _test_Vs(W1,W2):
    Ms=W1['Ms']
    Mf=W1['Mf']
    nc=W1['nc']
    Ml1=W1['Mf']-W1['Ms']+1
    Ml2=W2['Mf']-W2['Ms']+1
    for i in range(nc):
        for j in range(nc):
            for r in range(Mf):
                for k in range(Mf):
                    a=W1['V'][r+i*Ml1-W1['Ms'],k+j*Ml1-W1['Ms']]
                    b=W2['V'][r+i*Ml2-W2['Ms'],k+j*Ml2-W2['Ms']]
                    t=a-b
                    if(abs(t)>0.1):
                        
                        print i,j,':',r,k,':',a,b,'diff:',t
                        print "ia=",r+i*Ml1-W1['Ms'],':',k,'+',j,'*',Ml1,'-',W1['Ms'],'=',k+j*Ml1-W1['Ms']
                        print "ib=",r+i*Ml2-W2['Ms'],k+j*Ml2-W2['Ms']
                        return 





def _test_lin_comb(F,G,x0,x1,N=100):
    r"""
    Test if the coefficients of F+xG are integral for x0<=x <=x1
    """
    prec=F.space().prec()
    RF=RealField(prec)
    h=RF(x1-x0)/RF(N)
    ermin=1.0
    for j in range(N):
        x = x0+j*h
        P = F._lin_comb(G,1,x)
        c0 = P._coeffs[0][0][0]
        er_loc=0
        c = dict()
        for k in range(1,9):
            c[k]=RR((P._coeffs[0][0][k]/c0).real)
            ## remember we might have rational numbers.
            # hopefully they only have 2 or 3 in the denominator...
            er = abs(nearest_integer(c[k]*RF(6))-c[k]*RF(6))
            if abs(c[k]) > 0.01:
                if(er>er_loc):
                    er_loc=er
        if er_loc<0.1:
            print "err loc=",er_loc
            for k in range(1,9):
                er = abs(nearest_integer(c[k])-c[k])
                if abs(c[k]) > 0.01:
                    print x,k,c[k],er
        if er_loc<ermin:
            ermin=er_loc
            xmin=x
            Pmin=P
    print "xmin=",xmin
    print "ermin=",ermin
    return Pmin

def norm_sci_pretty_print(c,nd=0,es='e',latex_pow=False,zero_lim=0):
    if(is_int(c)):
        return str(c)
    if(abs(c)<zero_lim):
        return "0"
    if hasattr(c,"ae"):
        x=c.real; y=c.imag
    else:
        x=c.real(); y=c.imag()
    if(abs(x)>1E-5):
        sx = sci_pretty_print(x,nd,'',latex_pow)
    elif(abs(x)>zero_lim):
        sx = sci_pretty_print(x,2,'',latex_pow)
    else:
        sx=""
    # print x,y
    if(y>0 and sx<>""):
        p="+"
    else:
        p=""
    if(abs(y)>1E-5):
        sy=p+sci_pretty_print(y,nd,'',latex_pow)+"i"
    elif(abs(y)>zero_lim):
        sy=p+sci_pretty_print(y,2,'',latex_pow)+"i"
    else:
        sy=""
    ## un-scientify numbers between 0.1 and 1, i.e. with exponent -01
    if(sx.find("10^{-01}")>0):
        ss=sx.replace("\cdot 10^{-01}","")
        sx=ss.replace(".","")
        if(x<0):
            sx="-0."+sx.replace("-","")
        else:            
            sx="0."+sx
    if(sy.find("10^{-01}")>0):
        ss=sy.replace("\cdot 10^{-01}","")
        sy=ss.replace(".","")
        if(y<0):
            sy="-0."+sy.replace("-","")
        else:            
            sy="0."+sy

    s=sx+sy
    ss=s.replace("\cdot 10^{00}","")
    
    #print c,x,y,"::",ss
    return ss

def sci_pretty_print(s,nd=0,es='e',latex_pow=False):
    r""" Take a string representation of a number and returns it in scientific notation to desired number of digits.
    """
    # remove leading sign #
    #x=mpmath.mp.mpf(s)
    #if(abs(x)<1):
    #    raise NotImplementedError," Only implemented for |x|>1!!  got x=%s" %x
    if(not isinstance(s,str)):
        s=str(s)
        
    s=s.replace("(","")
    s=s.replace(")","")
    s=s.strip()
    if(s.count("I")+s.count("i")+s.count("j")>0):
        # Get a default complex notation
        s=s.replace("*","")
        s=s.replace("I","i")
        s=s.replace("j","i")
        ## We have to find imaginary and real parts 
        l=s.split("+")                
        if len(l)>1:
            (s1,s2)=l
        else:
            (s1,s2)=s.split("-")                
            s2="-"+s2
        if(s1.count("i")>0): # put imaginary part in standard form
            ims=s1.strip("i").lstrip(); res=s2             
        else:
            ims=s2.strip("i").lstrip(); res=s1 
        if(ims==""): ims="1"
        sres=sci_pretty_print(res,nd,es,latex_pow)    
        sims=sci_pretty_print(ims,nd,es,latex_pow)    
        if sres=="0": sres=""
        if sims=="0":
            sims=""
        else:
            sims=sims+"i"

        if sims.count("-")>0:
            return sres+" "+sims.replace(" -"," - ")
        elif sims<>"" and sres<>"":
            return sres+" + "+sims
        elif sres<>"":
            return sres
        elif sims<>"":
            return sims
        else:
            raise ValueError,"Could not find pretty print for s=%s " %s 
    s=s.strip()    
    if len(s.replace(".","").strip("0"))==0:
        return "0"
    if s.count(".")==0:
        s=s+".0"
    if s[0]=='-':
        ss=s.strip("-")
        ss=sci_pretty_print(ss,nd,es,latex_pow)    
        return "-"+ss

    l=s.split(".")
    if len(l)>1:
        (sint,sdigs)=l
    elif len(s)<nd:
        return s
    elif len(l)>0:
        sint=l[0]
        sdigs=""
    else:
        raise ValueError," Can not do pretty print for s=%s" %s

    if sdigs.count("e")>0:
        l=sdigs.split("e")
        sdigs=l[0]
        ex=int(l[1])
    else:
        ex=0
    if len(sint)==1 and sint=="0":
        # find the number of leading zeros
        sss=sdigs.lstrip("0")        
        nz=len(sdigs)-len(sss)+1
        if nz<10:
            ex="-0"+str(nz)
        else:
            ex=str(-nz)
        # Fix correct rounding
        rest=sss[nd:len(sss)]
        ix=nd-1
        if len(rest)>0:
            if int(rest) < 5*10**(len(rest)-1):
                # print " round < : since "+rest+"<"+str(5*10**(len(rest)-1))
                d=int(sss[ix])
            elif int(rest) > 5*10**(len(rest)-1):
                # print " round > : since "+rest+">"+str(5*10**(len(rest)-1))
                d=int(sss[ix])+1
            else:
                # if we have an exact half we round randomly 
                d=int(sss[ix])+int(random())            
            if d<10:
                ssdigs=sss[0:ix]+str(d) # We account for the leading digit too
            else:
                ssdigs=sss[0:ix-1]+str(int(sss[ix-1])+1)+str(d-10) # We account for the leading digit too
            if latex_pow:
                return ssdigs[0]+"."+ssdigs[1:nd]+"\cdot 10^{"+ex+"}"
            else:
                return ssdigs[0]+"."+ssdigs[1:nd]+es+ex            
    ex=int(ex)+len(sint)-1
    if abs(ex)<10:
        ex="0"+str(ex)
    else:
        ex=str(ex)
    #ssdigs=sint[1:len(sint)]+sdigs
    # cut away to nd digits
    if nd>0:
        #ssdigs=sdigs[0:nd-1] # We acount the leading digit too
        # Try to do correct rounding        
        rest=sdigs[nd-len(sint):len(sdigs)]
        #print "sdigs=",sdigs," nd=",nd
        #print "rest=",rest
        ix=nd-len(sint)-1
        if len(rest)>0:
            if int(rest) < 5*10**(len(rest)-1):
                # print " round < : since "+rest+"<"+str(5*10**(len(rest)-1))
                d=int(sdigs[ix])
            elif int(rest) > 5*10**(len(rest)-1):
                # print " round > : since "+rest+">"+str(5*10**(len(rest)-1))
                d=int(sdigs[ix])+1
            else:
                # if we have an exact half we round randomly 
                random.seed()
                d=int(sdigs[ix])+int(random.getrandbits(1))            
            if d<10:
                ssdigs=sdigs[0:ix]+str(d) # We account for the leading digit too
            else:
                if ix>0:
                    ssdigs=sdigs[0:ix-1]+str(int(sdigs[ix-1])+1)+str(d-10) # We account for the leading digit too
                else:
                    ssdigs=str(d-10) # We account for the leading digit too
                    if len(sint)==1:
                        sint=str(int(sint)+1)
                    else:
                        ll=len(sint); stmp=sint;
                        sint=stmp[1:ll-1]+str(int(stmp[ll-1])+1)
        else:
            ssdigs=sdigs[0:ix] 
        #print "rest=",rest,len(rest)
        ssdigs=sint[1:len(sint)]+ssdigs    
    else:
        ssdigs=sint[1:len(sint)]+sdigs    
        
    # print sint[0]+"."+ssdigs+" \cdot 10^{"+ex+"}"
    if latex_pow:
        res=sint[0]+"."+ssdigs+" \cdot 10^{"+ex+"}"
        return res 
    else:
        return sint[0]+"."+ssdigs+es+ex


def _set_character(character):
    if isinstance(character,str):
        if ["","trivial"].count(character)==0:
            raise NotImplemented,"Incorrect character! Got: %s" %multiplier
        character = None
    elif isinstance(character,sage.modular.dirichlet.DirichletCharacter) or isinstance(character,function):
        pass
    elif is_int(character):
        self._character=DirichletGroup(self.level()).list()[character]
    else:
        raise TypeError," Got an unknown character : %s " % character
    return character



### Functions for error estimates
def c_Ax(A,x):
    """ This function shows up in the estimate of the incomplete Gamma function"""
    if A<1:
        return x**(A-1)
    elif A>1:
        return A*x**(A-1)
    else:
        return 1

def incgamma_upper_bound(A,x):
    """  A bound for Gamma(A,x) """
    return exp(-x)*c_Ax(A,x)


def sum_of_gamma_bound(M,alpha,A,X):
    """ A bound of Sum_{M+1} n^{alpha} Gamma(A,nX) """
    if X<0:
        raise ValueError,"Can not bound this sum for X={0}".format(X)
    if M<alpha/X:
        M0 = ceil(alpha/X)
        print "Need to increase M to {0}".format(ceil(alpha/X))
        return -1,ceil(alpha/X)
    f = exp(-X*M)*X**-alpha
    if A<1:
        return f*c_Ax(alpha+A-1,X*M)
    elif A==1:
        return f/X*c_Ax(alpha+1,X*M)
    else:
        return f*A*c_Ax(alpha+A-1,X*M)
    
def error_bound_minus(k,M,Y):
    """
    Bound truncated non-holomorphic part of a Harmonic weak Maass form
    """
    if Y>=1: ## It is the level and not the height
        Y0 = sqrt(3.0)/2/Y
    else:
        Y0=Y
    fourpi = RR.pi()*RR(4)
    c1 = max(1.0-k,1.0)
    c2 = max(1.5*k+0.75,1.0)
    f1 = 6.0/fourpi**2*Y0**(-k-1)
    f2 = M**(0.5*k+0.75)
    f3 = exp(-RR.pi()*2*Y0*M)
    return c1*c2*f1*f2*f3
    
        
