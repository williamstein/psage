# -*- coding: utf-8 -*-
#*****************************************************************************
#  Copyright (C) 2011 Fredrik Strömberg <fredrik314@gmail.com>
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
Implements a class of multiplier systems which can be used to define automorphic forms.

AUTHORS:

- Fredrik Strömberg


EXAMPLES::

sage: K=CyclotomicField(168)
sage: z=K.gen()
sage: rT=matrix(K,3,3,[z**-1,0,0,0,z**-25,0,0,0,z**-121])
sage: rS=matrix(ComplexField(53),3,3) 
sage: for i in range(3):^J    for j in range(3):^J        rS[i,j]=2*sin(pi*(j+1)*(i+1)/7)    
sage: fak=CC(1/sqrt(7*I))


"""
from __future__ import print_function
from __future__ import division

from builtins import str
from builtins import range
from sage.all import SageObject,CyclotomicField,Integer,is_even,ZZ,QQ,Rational,kronecker,is_odd,SL2Z,Gamma0,matrix,floor,ceil,lcm,copy,trivial_character,var,DirichletGroup
from sage.all import kronecker_character,kronecker_character_upside_down,I
from sage.misc.cachefunc import cached_function,cached_method
from sage.modular.arithgroup.arithgroup_element import ArithmeticSubgroupElement
from sage.modular.etaproducts import qexp_eta
from psage.modform.arithgroup.mysubgroups_alg import factor_matrix_in_sl2z,SL2Z_elt
from psage.modform.arithgroup.mysubgroup import MySubgroup
from psage.modules.weil_module import WeilModule
from psage.groups.dirichlet_conrey import DirichletGroup_conrey,DirichletCharacter_conrey
class MultiplierSystem(SageObject):
    r"""
    Base class for multiplier systems.
    A multiplier system is a function:
    v : Gamma - > C^dim
    s.t. there exists a merom. function of weight k f:H->C^dim
    with f|A=v(A)f
    """
    def __init__(self,group,dchar=(1,1),dual=False,is_trivial=False,dimension=1,**kwargs):
        r"""
        if dual is set to true we use the complex conjugate of the representation (we assume the representation is unitary)


        
        The pair dchar = (conductor,number) gives the character in conrey numbering
        If char_nr = -1  = > kronecker_character
        If char_nr = -2  = > kronecker_character_upside_down
        """
        self._group = group
        if not hasattr(self,'_weight'):
            self._weight = None
        self._dim = dimension
        self._ambient_rank=kwargs.get('ambient_rank',None)
        self._kwargs = kwargs
        self._verbose = kwargs.get("verbose",0)
        if self._verbose>0:
             print("Init multiplier system!")
        (conductor,char_nr)=dchar
        self._conductor=conductor
        self._char_nr=char_nr
        self._character = None
        self._level = group.generalised_level()
        if 'character' in kwargs:
            if str(type(kwargs['character'])).find('Dirichlet')>=0:
                self._character = kwargs['character']
                self._conductor=self._character.conductor()
                try:
                    self._char_nr=self._character.number()
                except:
                    for x in DirichletGroup_conrey(self._conductor):
                        if x.sage_character() ==  self._character:
                            self._char_nr=x.number()
        elif group.is_congruence():
            if conductor<=0:
                    self._conductor=group.level(); self._char_nr=1
            if char_nr>=0:
                self._char_nr=char_nr
            
            if self._char_nr==0 or self._char_nr==1:
                self._character = trivial_character(self._conductor)
            elif self._char_nr==-1:                
                if self._conductor % 4 == 3:
                    self._character = kronecker_character(-self._conductor)
                elif self._conductor % 4 == 1:
                    self._character = kronecker_character(self._conductor)                    
                assert self._character.conductor()==self._conductor
            elif self._char_nr<=-2:
                self._character = kronecker_character_upside_down(self._conductor)    
            else:
                self._character = DirichletCharacter_conrey(DirichletGroup_conrey(self._conductor),self._char_nr).sage_character()
        else:
            self._conductor = 1
            self._character = trivial_character(1)
        #if not hasattr(self._character,'is_trivial'):
        #    if isinstance(self._character,(int,Integer)) and group.is_congruence():
        #        j = self._character
        #        self._character = DirichletGroup(group.level())[j]
        ## Extract the class name for the reduce algorithm
        self._class_name=str(type(self))[1:-2].split(".")[-1]
        if not isinstance(dimension,(int,Integer)):
            raise ValueError("Dimension must be integer!")
        self._is_dual = dual
        self._is_trivial=is_trivial and self._character.is_trivial()
        if is_trivial and self._character.order()<=2:
            self._is_real=True
        else:
            self._is_real=False
        self._character_values = [] ## Store for easy access

    def __getinitargs__(self):
        #print "get initargs"
        return (self._group,(self._conductor,self._char_nr),self._is_dual,self._is_trivial,self._dim)
        
    def __reduce__(self):
        #print "reduce!"
        t = self.__getinitargs__()
        return self.__class__,t 
    #return(TrivialMultiplier,(self._group,self._dim,(self._conductor,self._char_nr),self._is_dual,self._is_trivial))
    #    return(type(self),(self._group,self._dim,(self._conductor,self._char_nr),self._is_dual,self._is_trivial))

    def dual_multiplier(self):
        r"""
        Returns the dual multiplier of self.
        """
        m = copy(self)
        m._is_dual = int(not self._is_dual)
        o = self._character.order()
        m._character = self._character**(o-1)
        m._weight = QQ(2) - QQ(self.weight())
        return m
    
    def __repr__(self):
        r"""
        Needs to be defined in subclasses.
        """
        raise NotImplementedError

    def group(self):
        return self._group

    def is_dual(self):
        return int(self._is_dual)
    def level(self):
        return self._level
    def weight(self):
        r"""
        Return (modulo 2) whish weight self is a multiplier system consistent with
        """
        if self._is_trivial:
            self._weight = 0
            return self._weight
        if self._dim == 1:
            if self._weight is None:
                Z=SL2Z([-1,0,0,-1])
                v = self._action(Z)
                if v == 1:
                    self._weight = QQ(0)
                if v == -1:
                    self._weight = QQ(1)
                elif v == -I:
                    self._weight  = QQ(1)/QQ(2)
                elif v == I:
                    self._weight = QQ(1)/QQ(2)
                else:
                    raise NotImplementedError
        return self._weight
    
    def __call__(self,A):
        r"""
        For eficientcy we should also allow to act on lists
        """
        if isinstance(A,(ArithmeticSubgroupElement,SL2Z_elt,list)):
            if A not in self._group:
                raise ValueError("Element {0} is not in {1}! ".format(A,self._group))
            return self._action(A)
        else:
            raise NotImplementedError("Do not know how the multiplier should act on {0}".format(A))

        
        
    def _action(self):
        raise NotImplemented(" Needs to be overridden by subclasses!")
    def is_trivial(self):
        return self._is_trivial
    def is_real(self):
        return self._is_real

    def set_dual(self):
        self._is_dual = True #not self._is_dual

    def character(self):
        return self._character

    def weil_module(self):
        if hasattr(self,"_weil_module"):
            if self._verbose>1:
                print("self has weil_module!")
            return self._weil_module
        return None
            
    
    def rank(self):
        return self._dim

    def ambient_rank(self):
        r"""
        Return the dimension of the unsymmetrized space containing self.
        """
        if self._ambient_rank!=None:
            return self._ambient_rank
        else:
            return self._dim
        
    def __eq__(self,other):
        r"""
        A character is determined by the group it is acting on and its type, which is given by the representation.
        """
        if str(type(other)).find('Multiplier')<0:
            return False
        if self._group != other._group:
            return False
        if self._dim != other._dim:
            return False
        return self.__repr__()==other.__repr__()

    def __ne__(self,other):
        return not self.__eq__(other)
    
    def is_consistent(self,k):
        r"""
        Checks that v(-I)=(-1)^k,
        """
        Z=SL2Z([-1,0,0,-1])
        zi=CyclotomicField(4).gen()
        v = self._action(Z)
        if self._verbose>0: 
            print("test consistency for k={0}".format(k))
            print("v(Z)={0}".format(v))
        if self._dim==1:
            if isinstance(k,Integer) or k.is_integral():
                if is_even(k):
                    v1 = ZZ(1)
                else:
                    v1 = ZZ(-1)
            elif isinstance(k,Rational) and (k.denominator()==2 or k==0):
                v1 = zi**(-QQ(2*k))
                if self._verbose>0: 
                     print("I**(-2k)={0}".format(v1))
            else:
                raise ValueError("Only integral and half-integral weight is currently supported! Got weight:{0} of type:{1}".format(k,type(k)))
        else:
            raise NotImplemented("Override this function for vector-valued multipliers!")
        return v1==v
    
        
class TrivialMultiplier(MultiplierSystem):
    r"""
    The trivial multiplier.
    """
    def __init__(self,group,dchar=(0,0),dual=False,is_trivial=True,dimension=1,**kwargs):
        #print "kwargs0=",kwargs
        MultiplierSystem.__init__(self,group,dchar=dchar,dual=dual,is_trivial=True,dimension=dimension,**kwargs)

    def __repr__(self):
        if self._character != None and not self._character.is_trivial():
            s="Character "+str(self._character)
        else:
            s="Trivial multiplier!"
        return s
    def __getinitargs__(self):
        #print "get initargs"
        return (self._group,(self._conductor,self._char_nr),self._is_dual,self._is_trivial,self._dim)
        
    
    def _action(self,A):
        if self._character != None:
            if isinstance(A,(list,SL2Z_elt)):
                d=A[3]
            else:
                d=A[1,1]
            if len(self._character_values)>0:
                d = d % self._character.modulus()
                return self._character_values[d]
            else:
                return self._character(d)
        else:
            return 1

    def set_prec(self,prec=None):
        if prec in ["float","double"]:
            self._prec=prec
        else:
            prec = None

    def order(self):
        return 1
        
    def set_character_values(self):
        l = list(self._character.values())
        self._character_values=[]
        if self._prec=="float":
            for x in l:
                self._character_values.append(complex(x))
        elif self._prec=="double":
            for x in l:
                self._character_values.append(complex(x))
            
class ThetaMultiplier(MultiplierSystem):
    r"""
    Theta multiplier
    """
    def __init__(self,group,dchar=(0,0),dual=False,is_trivial=False,dimension=1,**kwargs):
        if not ZZ(4).divides(group.level()):
            raise ValueError(" Need level divisible by 4. Got: {0} ".format(self._group.level()))
        MultiplierSystem.__init__(self,group,dchar=dchar,dual=dual,is_trivial=is_trivial,dimension=dimension,**kwargs)
        self._i = CyclotomicField(4).gen()
        self._one = self._i**4
        self._weight= QQ(kwargs.get("weight",QQ(1)/QQ(2)))
        ## We have to make sure that we have the correct multiplier & character 
        ## for the desired weight
        if self._weight!=None:
            if floor(2*self._weight)!=ceil(2*self._weight):
                raise ValueError(" Use ThetaMultiplier for half integral or integral weight only!")
            t = self.is_consistent(self._weight)
            if not t:
                self.set_dual()
                t1 = self.is_consistent(self._weight)    
                if not t1:
                    raise ArithmeticError("Could not find consistent theta multiplier! Try to add a character.")
                
    def __repr__(self):
        s="Theta multiplier"
        if self._is_dual:
            s+=" v^-1 "
        else:
            s+=" v "
        if self._character != None and not self._character.is_trivial():
            s+=" and character "+str(self._character)
        return s
    def __latex__(self):
        if self._is_dual:
            s=" \bar{v_{\theta}} "
        else:
            s=" v_{\theta} "
        if self._character != None and not self._character.is_trivial():
            s+=" \cdot "
            if self._character == kronecker_character(self._conductor):
                s+=" \left( \frac\{\cdot\}\{ {0} \}\right)".format(self._conductor)
            elif self._character == kronecker_character_upside_down(self._conductor):
                s+=" \left( \frac\{ {0} \}\{ \cdot \}\right)".format(self._conductor)
            else:
                s+=" \chi_\{ {0}, {1} \}".format(self._conductor,self._char_nr)
        return s

    def __getinitargs__(self):
        #print "get initargs"
        return (self._group,(self._conductor,self._char_nr),self._is_dual,self._is_trivial,self._dim,self._weight)
    def order(self):
        return 4

    def _action(self,A):
        [a,b,c,d]=A
        v=kronecker(c,d)*self._one  ### I want the type of the result always be a number field element
        if(d % 4 == 3):
            v=-v*self._i
        elif (c % 4 != 0):
            raise ValueError("Only use theta multiplier for 4|c!")
        if self._character!=None:
            v = self._character(d)*v
        if self._is_dual:
            v = v**-1
        return v

class EtaMultiplier(MultiplierSystem):
    r"""
    Eta multiplier. Valid for any (real) weight.
    """
    def __init__(self,G,k=QQ(1)/QQ(2),number=0,ch=None,dual=False,version=1,dimension=1,**kwargs):
        r"""
        Initialize the Eta multiplier system: $\nu_{\eta}^{2(k+r)}$.
        INPUT:

        - G -- Group
        - ch -- character
        - dual -- if we have the dual (in this case conjugate)
        - weight -- Weight (recall that eta has weight 1/2 and eta**2k has weight k. If weight<>k we adjust the power accordingly.
        - number -- we consider eta^power (here power should be an integer so as not to change the weight...)
                
        """
        self._weight=QQ(k)
        if floor(self._weight-QQ(1)/QQ(2))==ceil(self._weight-QQ(1)/QQ(2)):
            self._half_integral_weight=1
        else:
            self._half_integral_weight=0
        MultiplierSystem.__init__(self,G,character=ch,dual=dual,dimension=dimension)
        number = number % 12
        if not is_even(number):
            raise ValueError("Need to have v_eta^(2(k+r)) with r even!")
        self._pow=QQ((self._weight+number)) ## k+r
        self._k_den=self._pow.denominator()
        self._k_num=self._pow.numerator()
        self._K = CyclotomicField(12*self._k_den)
        self._z = self._K.gen()**self._k_num
        self._i = CyclotomicField(4).gen()
        self._fak = CyclotomicField(2*self._k_den).gen()**-self._k_num
        self._version = version
        
        self.is_consistent(k) # test consistency

    def __repr__(self):
        s="Eta multiplier "
        if self._pow != 1:
            s+="to power 2*"+str(self._pow)+" "
        if self._character != None and not self._character.is_trivial():
            s+=" and character "+str(self._character)
        s+="with weight="+str(self._weight)
        return s
        
    def order(self):
        return 12*self._k_den

    def z(self):
        return self._z
    
     
    def _action(self,A):
        if self._version==1:
            return self._action1(A)
        elif self._version==2:
            return self._action2(A)
        else:
            raise ValueError

    def _action1(self,A):
        [a,b,c,d]=A
        return self._action0(a,b,c,d)
    def _action0(self,a,b,c,d):
        r"""
        Recall that the formula is valid only for c>0. Otherwise we have to use:
        v(A)=v((-I)(-A))=sigma(-I,-A)v(-I)v(-A).
        Then note that by the formula for sigma we have:
        sigma(-I,SL2Z[a, b, c, d])=-1 if (c=0 and d<0) or c>0 and other wise it is =1.
        """

        fak=1
        if c<0:
            a=-a; b=-b; c=-c;  d=-d; fak=-self._fak
        if c==0:
            if a>0:
                res = self._z**b
            else:
                res = self._fak*self._z**b
        else:
            if is_even(c):
                arg = (a+d)*c-b*d*(c*c-1)+3*d-3-3*c*d
                v=kronecker(c,d)
            else:
                arg = (a+d)*c-b*d*(c*c-1)-3*c
                v=kronecker(d,c)
            if not self._half_integral_weight:
                # recall that we can use eta for any real weight
                v=v**(2*self._weight)
            arg=arg*(self._k_num)
            res = v*fak*self._z**arg
            if self._character:
                res = res * self._character(d)
        if self._is_dual:
            res=res**-1
        return res


    def _action2(self,A):
        [a,b,c,d]=A
        fak=1
        if c<0:
            a=-a; b=-b; c=-c;  d=-d; fak=-self._fak
        if c==0:
            if a>0:
                res = self._z**b
            else:
                res = self._fak*self._z**b
        else:
            arg = dedekind_sum(-d,c)
            arg = arg+QQ(a+d)/QQ(12*c)-QQ(1)/QQ(4)
            # print "arg=",arg
            arg=arg*QQ(2)
            den = arg.denominator()*self._k_den
            num = arg.numerator()*self._k_num
            K = CyclotomicField(2*den)
            z=K.gen()
            if z.multiplicative_order()>4:
                fak=K(fak)
                # z = CyclotomicField(2*arg.denominator()).gen()
            res = z**num #rg.numerator()
            if self._character:
                ch = self._character(d)
                res=res*ch
            res = res*fak
        if self._is_dual:
            return res**-1
        return res

class TestMultiplier(MultiplierSystem):
    r"""
    Test of multiplier for f(q). As in e.g. the paper of Bringmann and Ono.
    """
    def __init__(self,group,dchar=(0,0),dual=False,weight=QQ(1)/QQ(2),dimension=1,version=1,**kwargs):
        self._weight=QQ(weight)
        MultiplierSystem.__init__(self,group,dchar=dchar,dual=dual,dimension=dimension,**kwargs)
        self._k_den=self._weight.denominator()
        self._k_num=self._weight.numerator()
        self._K = CyclotomicField(12*self._k_den)
        self._z = self._K.gen()**self._k_num
        self._sqrti = CyclotomicField(8).gen()
        self._i = CyclotomicField(4).gen()
        self._fak = CyclotomicField(2*self._k_den).gen()**-self._k_num
        self._fak_arg=QQ(self._weight)/QQ(2)
        self._version = version
        self.is_consistent(weight) # test consistency


    def order(self):
        return 12*self._k_den

    def z(self):
        return self._z

    def __repr__(self):
        s="Test multiplier"
        if self._character != None and not self._character.is_trivial():
            s+="and character "+str(self._character)
        return s


        
    def _action(self,A):
        [a,b,c,d]=A
        fak=0
        if c<0:
            a=-a; b=-b; c=-c;  d=-d; fak=self._fak_arg
        if c==0:
            if a>0:
                res = self._z**-b
            else:
                res = self._fak*self._z**-b
        else:
            arg=-QQ(1)/QQ(8)+QQ(c+a*d+1)/QQ(4)-QQ(a+d)/QQ(24*c)-QQ(a)/QQ(4)+QQ(3*d*c)/QQ(8)
            # print "arg=",arg
            arg = arg-dedekind_sum(-d,c)/QQ(2)+fak #self._fak_arg
            den=arg.denominator()
            num=arg.numerator()
            # print "den=",den
            # print  "num=",num
            res = self._K(CyclotomicField(den).gen())**num
            #res = res*fak
        if self._is_dual:
            return res**-1
        return res


class MultiplierByGenerator(MultiplierSystem):
    r"""
    Give a generator by its values on the generators.
    For now, only implemented for SL2Z.


    EXAMPLE:

    sage: K = CyclotomicField(24)
    sage: z=K.gen()
    sage: rS = matrix(K,3,3,[0, z**-3, 0, z**-3, 0, 0, 0, 0, -z**-3])
    sage: rT = matrix(K,3,3,[z**-1, 0, 0, 0, 0, z**8, 0, z**8, 0])
    sage: rZ = rS*rS
    sage: rho=MultiplierByGenerator(SL2Z,rS,rT,rZ)

    sage: te=TestMultiplier(Gamma0(2),weight=1/2)
    sage: r=InducedRepresentation(Gamma0(2),v=te)



    
    """
    def __init__(self,group,gens=[],vs=[],**kwargs):
        if len(gens)!=len(vs):
            raise ValueError("Need generators and values of the same form!")
        if hasattr(vs[0],'nrows'):
            dim=vs[0].nrows()
            assert vs[0].ncols()==dim and vs[0].nrows()==dim
        else:
            dim=1
        MultiplierSystem.__init__(self,group,dimension=dim)
        for i in range(len(gens)):
            self.vals[gens[i]]=vs[i]
        #if Z<>None:
        #self.T=rT
        #self.S=rS
        #self.Z=rZ
        
    def __repr__(self):
        s="Multiplier system defined by action on the generators:"
        if self._group==SL2Z:
            S,T=SL2Z.gens()
            Z=S*S
            s+="r(S)=",self.vals[S]
            s+="r(T)=",self.vals[T]
            s+="r(Z)=",self.vals[Z]
        else:
            for g in self.vals.keys():
                s+="r(",g,")=",self.vals[g]
        return s
    
    def _action(self,A):
        #if A not in self._group:
        #    raise ValueError,"Action is only defined for {0}!".format(self._group) 
        a,b,c,d=A
        if self._group==SL2Z:
            [z,n,l]=factor_matrix_in_sl2z(int(a),int(b),int(c),int(d))
            # l,ep = factor_matrix_in_sl2z_in_S_and_T(A_in)
            res = self.T.parent().one()
            if z==-1 and self.Z != None:
                res = self.Z #v(SL2Z[-1,0,0,-1])
                # for j in range(self.dim):
                #    res[j,j]=tmp
            if n != 0:
                res = self.T**n
            for i in range(len(l)):
                res = res*self.S*self.T**l[i]
        elif A in self.vals.keys():
            res = self.vals[A]
        else:
            raise NotImplementedError("Can not write as word in generators of {0}".format(self._group))
        return res

class InducedRepresentationMultiplier(MultiplierSystem):
    def __init__(self,G,v=None,**kwargs):
        r"""
        G should be a subgroup of PSL(2,Z).

        EXAMPLE::

        sage: te=TestMultiplier(Gamma0(2),weight=1/2)
        sage: r=InducedRepresentation(Gamma0(2),v=te)

        """
        dim = len(list(G.coset_reps()))
        MultiplierSystem.__init__(self,Gamma0(1),dimension=dim)
        self._induced_from=G
        # setup the action on S and T (should be faster...)
        self.v = v        
        if v != None:

            k = v.order()
            if k>2:
                K = CyclotomicField(k)
            else:
                K=ZZ
            self.S=matrix(K,dim,dim)
            self.T=matrix(K,dim,dim)
        else:
            self.S=matrix(dim,dim)
            self.T=matrix(dim,dim)
        S,T=SL2Z.gens()
        if hasattr(G,"coset_reps"):
            if isinstance(G.coset_reps(),list):
                Vl=G.coset_reps()
            else:
                Vl=list(G.coset_reps())
        elif hasattr(G,"_G"):
            Vl=list(G._G.coset_reps())
        else:
            raise ValueError("Could not get coset representatives from {0}!".format(G))
        self.repsT=dict()
        self.repsS=dict()
        for i in range(dim):
            Vi=Vl[i]
            for j in range(dim):
                Vj=Vl[j]
                BS = Vi*S*Vj**-1
                BT = Vi*T*Vj**-1
                #print "i,j
                #print "ViSVj^-1=",BS
                #print "ViTVj^-1=",BT
                if BS in G:
                    if v != None:
                        vS=v(BS)
                    else:
                        vS=1
                    self.S[i,j]=vS
                    self.repsS[(i,j)]=BS
                if BT in G:
                    if v != None:
                        vT=v(BT)
                    else:
                        vT=1
                    self.T[i,j]=vT
                    self.repsT[(i,j)]=BT


    def __repr__(self):
        s="Induced representation from "
        s+=str(self.v)+" on the group "+str(self._induced_from)
        return s


    def induced_from(self):
        self._induced_from

    def _action(self,A):
        #if A not in self._group:
        #    raise ValueError,"Action is only defined for {0}!".format(self._group) 
        a,b,c,d=A
        [z,n,l]=factor_matrix_in_sl2z(int(a),int(b),int(c),int(d))
        #l,ep = factor_matrix_in_sl2z_in_S_and_T(A_in)
        res = copy(self.T.parent().one())
        if z==-1 and self.v != None:
            tmp = self.v(SL2Z([-1,0,0,-1]))
            for j in range(self._dim):
                res[j,j]=tmp
        if n != 0:
            res = self.T**n
        for i in range(len(l)):
            res = res*self.S*self.T**l[i]
        return res

    

class WeilRepMultiplier(MultiplierSystem):
    def __init__(self,WR,weight=QQ(1)/QQ(2),use_symmetry=True,group=SL2Z,dual=False,**kwargs):
        r"""
        WR should be a Weil representation.
        INPUT:
         - weight -- weight (should be consistent with the signature of self)
         - use_symmetry -- if False we do not symmetrize the functions with respect to Z
                        -- if True we need a compatible weight
        - 'group' -- Group to act on.
        """
        if isinstance(WR,(Integer,int)):
            #self.WR=WeilRepDiscriminantForm(WR)
            #self.WR=WeilModule(WR)
            self._weil_module = WeilModule(WR)
        elif hasattr(WR,"signature"):
            self._weil_module=WR
        else:
            raise ValueError("{0} is not a Weil module!".format(WR))
        self._sym_type = 0

        if group.level() != 1:
            raise NotImplementedError("Only Weil representations of SL2Z implemented!")
        self._group = MySubgroup(Gamma0(1))
        self._dual = int(dual)
        self._sgn = (-1)**self._dual
        self.Qv=[]
        self._weight = weight
        self._use_symmetry = use_symmetry
        self._kwargs = kwargs
        self._sgn = (-1)**int(dual)        
        half = QQ(1)/QQ(2)
        threehalf = QQ(3)/QQ(2)
        if use_symmetry:
            ## First find weight mod 2:
            twok= QQ(2)*QQ(weight)
            if not twok.is_integral():
                raise ValueError("Only integral or half-integral weights implemented!")
            kmod2 = QQ(twok % 4)/QQ(2)
            if dual:
                sig_mod_4 = -self._weil_module.signature() % 4
            else:
                sig_mod_4 = self._weil_module.signature() % 4
            sym_type = 0
            if (kmod2,sig_mod_4) in [(half,1),(threehalf,3),(0,0),(1,2)]:
                sym_type = 1
            elif (kmod2,sig_mod_4) in [(half,3),(threehalf,1),(0,2),(1,0)]:
                sym_type = -1
            else:
                raise ValueError("The Weil module with signature {0} is incompatible with the weight {1}!".format( self._weil_module.signature(),weight))
            ## Deven and Dodd contains the indices for the even and odd basis
            Deven=[]; Dodd=[]
            if sym_type==1:
                self._D = self._weil_module.even_submodule(indices=1)
            else: #sym_type==-1:
                self._D = self._weil_module.odd_submodule(indices=1)
            self._sym_type=sym_type
            dim = len(self._D)  #Dfinish-Dstart+1
        else:
            dim = len(self._weil_module.finite_quadratic_module().list())
            self._D = list(range(dim))
            self._sym_type=0
        if hasattr(self._weil_module,"finite_quadratic_module"):
#            for x in list(self._weil_module.finite_quadratic_module()):
            for x in self._weil_module.finite_quadratic_module():
                self.Qv.append(self._weil_module.finite_quadratic_module().Q(x))
        else:
            self.Qv=self._weil_module.Qv
        ambient_rank = self._weil_module.rank()
        MultiplierSystem.__init__(self,self._group,dual=dual,dimension=dim,ambient_rank=ambient_rank)


    def __repr__(self):
        s="Weil representation corresponding to "+str(self._weil_module)
        return s

    def _action(self,A):
        #return self._weil_module.rho(A)
        return self._weil_module.matrix(A)

    def sym_type(self):
        return self._sym_type
    
    def dual_multiplier(self):
        r"""
        Returns the dual multiplier of self.
        """
        weight = QQ(2) - QQ(self._weight)
        dual = int(not self._is_dual)
        m = WeilRepMultiplier(self._weil_module,weight,self._use_symmetry,self._group,dual=dual,**self._kwargs)
        return m
   
    def D(self):
        r"""
        Returns the indices of a basis of the (symmetrized) discriminant form of self.
        """
        return self._D

    @cached_method
    def neg_index(self,a):
        r"""
        The index of -a in self._D
        """
        return self.weil_module()._neg_index(a)
    
    def is_consistent(self,k):
        r"""
        Return True if the Weil representation is a multiplier of weight k.
        """
        if self._verbose>0:
            print("is_consistent at wr! k={0}".format(k))
        twok =QQ(2)*QQ(k)
        if not twok.is_integral():
            return False
        if self._sym_type != 0:
            if self.is_dual():
                sig_mod_4 = -self._weil_module.signature() % 4
            else:
                sig_mod_4 = self._weil_module.signature() % 4
            if is_odd(self._weil_module.signature()):                
                return (twok % 4 == (self._sym_type*sig_mod_4 %4))
            else:
                if sig_mod_4 == (1 - self._sym_type) % 4:
                    return twok % 4 == 0
                else:
                    return twok % 4 == 2
        if is_even(twok) and  is_even(self._weil_module.signature()):
                return True
        if is_odd(twok) and  is_odd(self._weil_module.signature()):
                return True
        return False

    def dimension_cusp_forms(self,k,eps=0):
        r"""
        Compute the dimension of the space of cusp forms of weight k with respect to
        the weil module of self if eps=1 and its dual if eps=-1
        
        """
        if eps==0:
            eps = self._sgn
        return self._weil_module.dimension_cusp_forms(k,eps)

    def dimension_modular_forms(self,k,eps=0):
        r"""
        Compute the dimension of the space of holomorphic modular forms of weight k with respect to
        the weil module of self if eps=1 and its dual if eps=-1
        
        """
        if eps==0:
            eps = self._sgn
        return self._weil_module.dimension_modular_forms(k,eps)    


class EtaQuotientMultiplier(MultiplierSystem):
    r"""
    Eta multiplier given by eta(Az)^{r}/eta(Bz)^s
    The weight should be r/2-s/2 mod 2.
    The group is Gamma0(lcm(A,B))
    """
    def __init__(self,args=[1],exponents=[1],ch=None,dual=False,version=1,**kwargs):
        r"""
        Initialize the Eta multiplier system: $\nu_{\eta}^{2(k+r)}$.
        INPUT:

        - G -- Group
        - ch -- character
        - dual -- if we have the dual (in this case conjugate)
        - weight -- Weight (recall that eta has weight 1/2 and eta**2k has weight k. If weight<>k we adjust the power accordingly.
        - number -- we consider eta^power (here power should be an integer so as not to change the weight...)

        EXAMPLE:
        
                
        """
        assert len(args) == len(exponents)
        self._level=lcm(args)
        G = Gamma0(self._level)
        k = sum([QQ(x)*QQ(1)/QQ(2) for x in exponents])
        self._weight=QQ(k)
        if floor(self._weight-QQ(1)/QQ(2))==ceil(self._weight-QQ(1)/QQ(2)):
            self._half_integral_weight=1
        else:
            self._half_integral_weight=0
        MultiplierSystem.__init__(self,G,dimension=1,character=ch,dual=dual)
        self._arguments = args
        self._exponents =exponents
        self._pow=QQ((self._weight)) ## k+r
        self._k_den = self._weight.denominator()
        self._k_num = self._weight.numerator()
        self._K = CyclotomicField(12*self._k_den)
        self._z = self._K.gen()**self._k_num
        self._i = CyclotomicField(4).gen()
        self._fak = CyclotomicField(2*self._k_den).gen()**-self._k_num
        self._version = version
        self.is_consistent(k) # test consistency

    def __repr__(self):
        s="Quotient of Eta multipliers :  "
        for i in range(len(self._arguments)):
            n = self._arguments[i]
            e = self._exponents[i]
            s+="eta({0}z)^{1}".format(n,e)
            if i < len(self._arguments)-1:
                s+="*"
        if self._character != None and not self._character.is_trivial():
            s+=" and character "+str(self._character)
        s+=" with weight="+str(self._weight)
        return s

    def level(self):
        return self._level
        
    def order(self):
        return 12*self._k_den

    def z(self):
        return self._z

    def q_shift(self):
        r"""
        Gives the 'shift' at the cusp at infinity of the q-series.
        The 'true' q-expansion of the eta quotient is then q^shift*q_expansion
        """
        num =  sum([self._argument[i]*self._exponent[i] for i in range(len(self._arguments))])

        return QQ(num)/QQ(24)
    
    def q_expansion(self,n=20):
        r"""
        Give the q-expansion of the quotient.
        """
        eta = qexp_eta(ZZ[['q']],n)
        R = eta.parent()
        q = R.gens()[0]
        res = R(1)
        prefak = 0
        for i in range(len(self._arguments)):        
            res = res*eta.subs({q:q**self._arguments[i]})**self._exponents[i]
            prefak = prefak+self._arguments[i]*self._exponents[i]
        if prefak % 24 == 0:
            return res*q**(prefak/QQ(24))
        else:
            return res,prefak/QQ(24)
        #etA= et.subs(q=q**self._arg_num).power_series(ZZ[['q']])
        #etB= et.subs(q=q**self._arg_den).power_series(ZZ[['q']])
        #res = etA**(self._exp_num)/etB**(self._exp_den)
        #return res
    #def _action(self,A):
    #    return self._action(A)
        
    def _action(self,A):
        [a,b,c,d]=A
        if not c % self._level == 0 :
            raise ValueError("Need A in {0}! Got: {1}".format(self.group,A))
        fak=1
        if c<0:
            a=-a; b=-b; c=-c;  d=-d; fak=-self._fak
            #fak = fak*(-1)**(self._exp_num-self._exp_den)
        res = 1
        exp = 0
        for i in range(len(self._exponents)):
            z = CyclotomicField(lcm(12,self._exponents[i].denominator())).gen()
            arg,v = eta_conjugated(a,b,c,d,self._arguments[i])
            #arg2,v2 = eta_conjugated(a,b,c,d,self._arg_den)
            #res=self._z**(arg1*self._exp_num-arg2*self._exp_den)
#            exp += arg*self._exponents[i]
            if v != 1:
                res=res*v**self._exponents[i]
            #if v2<>1:
            #res=res/v2**self._exp_den
            res = res*z**(arg*self._exponents[i].numerator())
#        res = res*self._z**exp
        if fak != 1:
            res=res*fak**exp
        return res

        
class EtaQuotientMultiplier_2(MultiplierSystem):
    r"""
    Eta multiplier given by eta(Az)^{r}/eta(Bz)^s
    The weight should be r/2-s/2 mod 2.
    The group is Gamma0(lcm(A,B))
    """
    def __init__(self,A,B,r,s,k=None,number=0,ch=None,dual=False,version=1,**kwargs):
        r"""
        Initialize the Eta multiplier system: $\nu_{\eta}^{2(k+r)}$.
        INPUT:

        - G -- Group
        - ch -- character
        - dual -- if we have the dual (in this case conjugate)
        - weight -- Weight (recall that eta has weight 1/2 and eta**2k has weight k. If weight<>k we adjust the power accordingly.
        - number -- we consider eta^power (here power should be an integer so as not to change the weight...)

        EXAMPLE:
        
                
        """
        self._level=lcm(A,B)
        G = Gamma0(self._level)
        if k==None:
            k = (QQ(r)-QQ(s))/QQ(2)
        self._weight=QQ(k)
        if floor(self._weight-QQ(1)/QQ(2))==ceil(self._weight-QQ(1)/QQ(2)):
            self._half_integral_weight=1
        else:
            self._half_integral_weight=0
        MultiplierSystem.__init__(self,G,dimension=1,character=ch,dual=dual)
        number = number % 12
        if not is_even(number):
            raise ValueError("Need to have v_eta^(2(k+r)) with r even!")
        self._arg_num = A
        self._arg_den = B
        self._exp_num = r
        self._exp_den = s
        self._pow=QQ((self._weight+number)) ## k+r
        self._k_den=self._pow.denominator()
        self._k_num=self._pow.numerator()
        self._K = CyclotomicField(12*self._k_den)
        self._z = self._K.gen()**self._k_num
        self._i = CyclotomicField(4).gen()
        self._fak = CyclotomicField(2*self._k_den).gen()**-self._k_num
        self._version = version
        self.is_consistent(k) # test consistency

    def __repr__(self):
        s="Quotient of Eta multipliers :  "
        s+="eta({0})^{1}/eta({2})^{3}".format(self._arg_num,self._exp_num,self._arg_den,self._exp_den)
        if self._character != None and not self._character.is_trivial():
            s+=" and character "+str(self._character)
        s+=" with weight="+str(self._weight)
        return s

    def level(self):
        return self._level
        
    def order(self):
        return 12*self._k_den

    def z(self):
        return self._z

    def q_shift(self):
        r"""
        Gives the 'shift' at the cusp at infinity of the q-series.
        The 'true' q-expansion of the eta quotient is then q^shift*q_expansion
        """
        num =  self._arg_num*self._exp_num-self._arg_den*self._exp_den
        return QQ(num)/QQ(24)
    
    def q_expansion(self,n=20):
        r"""
        Give the q-expansion of the quotient.
        """
        q = ZZ[['q']].gen()
        et = qexp_eta(ZZ[['q']],n)
        etA= et.subs(q=q**self._arg_num).power_series(ZZ[['q']])
        etB= et.subs(q=q**self._arg_den).power_series(ZZ[['q']])
        res = etA**(self._exp_num)*etB**(-self._exp_den)
        return res
    #def _action(self,A):
    #    return self._action(A)
        
    def _action(self,A):
        [a,b,c,d]=A
        if not c % self._level == 0 :
            raise ValueError("Need A in {0}! Got: {1}".format(self.group,A))
        fak=1
        if c<0:
            a=-a; b=-b; c=-c;  d=-d; fak=-self._fak
            #fak = fak*(-1)**(self._exp_num-self._exp_den)
        arg1,v1 = eta_conjugated(a,b,c,d,self._arg_num)
        arg2,v2 = eta_conjugated(a,b,c,d,self._arg_den)
        res=self._z**(arg1*self._exp_num-arg2*self._exp_den)
        if v1 != 1:
            res=res*v1**self._exp_num
        if v2 != 1:
            res=res*v2**(-self._exp_den)
        if fak != 1:
            res=res*fak**(self._exp_num-self._exp_den)
        return res
            
def eta_conjugated(a,b,c,d,l):
    r"""
    Gives eta(V_l A V_l^-1) with A=(a,b,c,d) for c>0
    """
    assert c>=0 and (c%l)==0
    if l != 1:
        cp = QQ(c)/QQ(l)
        bp = QQ(b)*QQ(l)
    else:
        cp=c; bp=b
    if c==0:
        l*bp,1
    res = (a+d)*cp-bp*d*(cp*cp-1)
    if is_odd(c):
        return res-3*cp,kronecker(d,cp)
    else:
        return res+3*d-3-3*cp*d,kronecker(cp,d)


def eta_argument(a,b,c,d):
    res = (a+d)*c-b*d*(c*c-1)
    den = 24
    if is_odd(c):
        return res-3*c,den,kronecker(d,c)
    else:
        return res+3*d-3-3*c*d,den,kronecker(d,c)



def saw_tooth_fn(x):
    if floor(x) == ceil(x):
        return 0
    elif x in QQ:
        return QQ(x)-QQ(floor(x))-QQ(1)/QQ(2)
    else:
        return x-floor(x)-0.5

def dedekind_sum(d,c):
    res=0
    for i in range(c):
        tmp = saw_tooth_fn(QQ(i)/QQ(c))
        tmp*= saw_tooth_fn(QQ(d*i)/QQ(c))
        res=res+tmp
    return res
