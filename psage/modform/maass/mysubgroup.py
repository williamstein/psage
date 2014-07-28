# -*- coding: utf-8 -*-
# cython: profile=True
r"""
A general class for subgroups of the (projective) modular group, PSL(2,Z).
Extends the standard classes with methods needed for Maass waveforms.

NOTE: The fact that we deal with subgroups of PSL(2,Z) and not SL(2,Z) is EXTREMELY important.
  More general algorithms may come later. 


AUTHORS:

 - Fredrik Strömberg


EXAMPLES::


   sage: P=Permutations(6)
   sage: pS=P([2,1,4,3,6,5])
   sage: pR=P([3,1,2,5,6,4])
   sage: G
   Arithmetic Subgroup of PSL2(Z) with index 6.
   Given by
   perm(S)=[2, 1, 4, 3, 6, 5]
   perm(ST)=[3, 1, 2, 5, 6, 4]
   Constructed from G=Arithmetic subgroup corresponding to permutations L=(2,3,5,4), R=(1,3,6,4)
   sage: TestSuite.run()


Commutator subgroup
   sage: pR=MyPermutation('(1 3 5)(2 4 6)')
   sage: pS=MyPermutation('(1 2)(3 4)(5 6)')
   sage: G=MySubgroup(o2=pS,o3=pR)
Gamma^3
   sage: pR=MyPermutation('(1 2 3)')
   sage: pS=MyPermutation('(1)(2)(3)')
   sage: G=MySubgroup(o2=pS,o3=pR)
   
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

#from sage.all_cmdline import *   # import sage library

from sage.rings.arith    import xgcd
from sage.rings.all import Integer,CC,ZZ,QQ,RR,RealNumber,I,infinity,Rational,gcd
from sage.rings.real_mpfr import RealNumber as RealNumber_class
#from sage.combinat.permutation import (Permutations,PermutationOptions)
from sage.modular.cusps import Cusp
#from Cusp import is_gamma0_equiv
from sage.modular.arithgroup.all import *
from sage.symbolic.expression import Expression
from sage.modular.modsym.p1list import lift_to_sl2z 
from sage.functions.other import ceil,floor,sqrt
from sage.all import Matrix,SageObject,numerator,denominator,copy,log_b,is_odd
from sage.modular.arithgroup import congroup_gamma0
from sage.modular.arithgroup.arithgroup_generic import ArithmeticSubgroup
from sage.modular.arithgroup.arithgroup_perm import EvenArithmeticSubgroup_Permutation
from sage.modular.arithgroup.congroup_gamma0 import Gamma0_class
from sage.rings.integer import is_Integer
from sage.groups.all import SymmetricGroup
from sage.rings.arith import lcm
from copy import deepcopy
from psage.modform.maass.mysubgroups_alg import * 
from psage.modform.maass.permutation_alg import MyPermutation 
#from psage.modform.maass.permutation_alg import MyPermutation,MyPermutationIterator
from plot_dom import draw_funddom_d,draw_funddom,my_hyperbolic_triangle
from psage.modform.maass.permutation_alg import are_transitive_permutations,num_fixed
from psage.modform.maass.sl2z_subgroups_alg import are_mod1_equivalent

from sage.plot.all import Graphics
from sage.plot.circle import circle
from sage.plot.colors import to_mpl_color
from sage.plot.misc import options, rename_keyword
#from sage.plot.all import hyperbolic_arc, hyperbolic_triangle, text

import types
import warnings
import sys,os
import matplotlib.patches as patches
import matplotlib.path as path

from sage.modular.arithgroup.arithgroup_perm import *
#from subgroups_alg import *
#load "/home/stromberg/sage/subgroups/subgroups_alg.spyx"

def MySubgroup(A=None,B=None,verbose=0,version=0,display_format='short',data={},**kwds):
    r"""
    Create an instance of MySubgroup_class.

    INPUT:

    - A -- can be 
    """
    s2 = None; s3=None; is_Gamma0=None; level = None
    is_symmetric = kwds.get('is_symmetric')
    symmetry_map = kwds.get('symmetry_map')
    if isinstance(A,MySubgroup_class):
        return MySubgroup_class(data=A.__dict__,**kwds)
    if isinstance(A,ArithmeticSubgroup):
        ## If A is not a subgroup of PSL(2,Z) so we have projectivize with .to_even_subgroup
        s2 = MyPermutation(A.as_permutation_group().to_even_subgroup().S2().domain())
        s3 = MyPermutation(A.as_permutation_group().to_even_subgroup().S3().domain())
        if A.is_congruence():
            level = A.level()
        if isinstance(A,Gamma0_class):
            is_Gamma0=True
        else:
            is_Gamma0=False            
    elif A<>None and B<>None:
        if isinstance(A,MyPermutation) and isinstance(B,MyPermutation):
            s2 = A; s3 = B
        else:
            try:
                s2 = MyPermutation(A.domain())
                s3 = MyPermutation(B.domain())
            except ValueError as ve:
                raise ValueError,"Can not construct permutations! {0}".format(ve)
    elif A<>None:
        if hasattr(A,"p2") and hasattr(A,"p3"):
            s2 = MyPermutation(A.p2); s3 = MyPermutation(A.p3)
        elif hasattr(A,"s2") and hasattr(A,"s3"):
            s2 = MyPermutation(A.s2); s3 = MyPermutation(A.s3)
    else:
        s2 = kwds.get("o2",None)
        s3 = kwds.get("o3",None)
    is_symmetric = kwds.get('is_symmetric',None); symmetry_map = kwds.get('symmetry_map',None)
    if hasattr(A,'__dict__'):
        is_symmetric = A.__dict__.get('is_symmetric',is_symmetric)
        symmetry_map = A.__dict__.get('symmetry_map',symmetry_map)
    if s2==None or s3==None:
        s2 = kwds.get("s2",None)
        s3 = kwds.get("s3",None)
    if s2==None or s3==None:
        raise ValueError,"Could not construct subgroup from input!"
    reps_from_farey = kwds.get('farey',None)
    if is_Gamma0:
        return MySubgroup_congruence_class(o2=s2,o3=s3,verbose=verbose,is_Gamma0=is_Gamma0,is_symmetric=1,symmetry_map=SL2Z_elt(1,0,0,1),reps_from_farey=reps_from_farey)
    return MySubgroup_class(o2=s2,o3=s3,verbose=verbose,is_Gamma0=is_Gamma0,level=level,is_symmetric=is_symmetric,symmetry_map=symmetry_map,reps_from_farey=reps_from_farey)


class MySubgroup_class (EvenArithmeticSubgroup_Permutation):
    r"""
    A class for subgroups of the modular group SL(2,Z).
    Extends the standard classes with methods needed for Maass waveforms.
    
    EXAMPLES::
    
    
    sage: G=MySubgroup(Gamma0(5));G
    Arithmetic Subgroup of PSL2(Z) with index 6.
    Given by
    perm(S)=[2, 1, 4, 3, 5, 6]
    perm(ST)=[3, 1, 2, 5, 6, 4]
    Constructed from G=Congruence Subgroup Gamma0(5)
    sage: G.is_subgroup(G5)
    True
    sage: G.cusps()
    [Infinity, 0]
    G.cusp_normalizer(G.cusps()[1])
    [ 0 -1]
    [ 1  0]
    
    
    """
    def __init__(self,o2=None,o3=None,verbose=0,display_format='short',data={},**kwds):
        r""" Init a subgroup in the following forms: 
          1. G = Arithmetic subgroup
          2. (o2,o3) = pair of transitive permutations of order 2 and
          order 3 respectively.
          
          INPUT:
          - 'G'  --  Subgroup of the modular group.
          - 'o2' --  Permutation of order 2. The permutation implementation needs to have method: .cycle_tuples()
          - 'o3' --  Permutation of order 3. 
          - 'str' -- String: either prepresenting permutations or a subgroup.
          - 'verbose' -- integer, set verbosity with positive values.
          - 'display_format' -- 'short' or 'long'
          INPUT TYPES:
          

          
          ATTRIBUTES:

            permS = permutation representating S
            permT = permutation representating T  
            permR = permutation representating R=ST  
            permP = permutation representating P=STS 
            (permR,permT) gives the group as permutation group
            Note:
              The usual permutation group has two parabolic permutations L,R
              L=permT, R=permP

          EXAMPLES::

          
              sage: G=SL2Z
              sage: MG=MySubgroup(G);MG
              Arithmetic Subgroup of PSL2(Z) with index 1.
              Given by
              perm(S)=[1]
              perm(ST)=[1]
              Constructed from G=Modular Group SL(2,Z)
              sage: P=Permutatons(6)
              sage: pS=[2,1,4,3,6,5]
              sage: pR=[3,1,2,5,6,4]
              sage: G=MySubgroup(o2=pS,o3=pR);G
              Arithmetic Subgroup of PSL2(Z) with index 6.
              Given by
              perm(S)=[2, 1, 4, 3, 6, 5]
              perm(ST)=[3, 1, 2, 5, 6, 4]
              Constructed from G=Arithmetic subgroup corresponding to permutations L=(2,3,5,4), R=(1,3,6,4)

          
        """
        self._verbose = verbose
        self._display_format = display_format
        self._level=kwds.get('level',None)
        if self._level<>None:
            self._generalised_level=self._level;  self._is_congruence=True
        else:
            self._generalised_level=None;  self._is_congruence=None
        self._perm_group=None;
        self._coset_rep_perms={}
        self._coset_rep_strings={}
        self._cusps_as_cusps=[]
        self._generators_as_slz_elts=[]
        self._symmetrizable_cusp={}
        self._cusp_normalizer_is_normalizer=None
        self._vertices=None; self._vertex_data=None
        self._cusps=None; self._nvertices=None
        self._cusp_data=[]; self._vertex_widths=[]
        self._vertex_maps=[]; self._cusp_maps=[]
        self.permT=None; self.permP=None; self.permS=None; self.permR=None
        self._coset_reps_v0 = None; self._coset_reps_v1 = None
        self._coset_reps_v2 = None
        self._nu2=None; self._nu3=None
        self._ncusps=None; self._genus=None 
        self._signature=None; self._verbose=verbose
        self._reps_from_farey = kwds.get('reps_from_farey',None)
        self._is_Gamma0=kwds.get('is_Gamma0',None)
        self._is_symmetric=kwds.get('is_symmetric')
        self._symmetry_map = kwds.get('symmetry_map')        
        self._translational_symmetry = None
        self._reflectional_symmetry = None        
        self._vertices_as_cusps = []
        self._cusp_data_sage_format = {}
        self._reflected_group = None
        self._symmetry_type_Ia = []; self._symmetry_type_Ib = []
        self._symmetry_type_IIa = []; self._symmetry_type_IIb = []
        self._checked_symmetry_type_Ia = 0; self._checked_symmetry_type_Ib = 0
        self._checked_symmetry_type_IIa = 0; self._checked_symmetry_type_IIb = 0        
        self._modular_correspondences = {}
        if self._verbose>1:
            print "o2=",o2
            print "o3=",o3
            print "str=",str
            print "is_Gamma0=",self._is_Gamma0
            print "kwds=",kwds
        if data<>{}:
            self.init_group_from_dict(data,**kwds)
        elif o2<>None and o3<>None:
            self.init_group_from_permutations(o2,o3)
        else:
            raise ValueError,"Incorrect input to subgroup! Got G={0}, o2={1} nad o3={2}".format(o2,o3)
        self._display_format = display_format                
        self._uid = self._get_uid()
        self.class_name='MySubgroup_class'            
        self._name = ''
        self._latex_name = ''        

    def _repr_(self):
        r"""
        Return the string representation of self.

        EXAMPLES::


            sage: P=Permutatons(6)
            sage: pS=[2,1,4,3,6,5]
            sage: pR=[3,1,2,5,6,4]
            sage: G=MySubgroup(o2=pS,o3=pR);G._repr_()
            Arithmetic Subgroup of PSL2(Z) with index 6.
            Given by
            perm(S)=[2, 1, 4, 3, 6, 5]
            perm(ST)=[3, 1, 2, 5, 6, 4]
            Constructed from G=Arithmetic subgroup corresponding to permutations L=(2,3,5,4), R=(1,3,6,4)
            sage: G=SL2Z
            sage: MG=MySubgroup(G);MG._repr_()
            Arithmetic Subgroup of PSL2(Z) with index 1.
            Given by
            perm(S)=[1]
            perm(ST)=[1]
            Constructed from G=Modular Group SL(2,Z)            

        """
        if self.index()==1:
            return "SL(2,Z)"
        if self._is_Gamma0==True:
            return "Gamma0({0})".format(self._level)
        s ="Arithmetic Subgroup of SL(2,Z) with index "+str(self._index)+". "
        s+="Given by: \n \t perm(S)="+str(self.permS)+"\n \t perm(ST)="+str(self.permR)
        if hasattr(self,"_display_format") and self._display_format=='long':
            s+="\nConstructed from G="+super(MySubgroup_class,self)._repr_()

        return s


    def _get_uid(self):
        r""" Constructs a unique identifier for the group

        OUTPUT::

            A unique identifier as string.
            The only truly unique identifier is the set of permutations
            so we return those as strings.
            Because of the extra stucture given by e.g. Gamma0(N) we also return
            this information if available.


        EXAMPLES::


            sage: G=MySubgroup(Gamma0(5))
            sage: G._get_uid()
            [2_1_4_3_5_6]-[3_1_2_5_6_4]-Gamma0(5)'
            sage: P=Permutations(6)
            sage: pS=P([2,1,4,3,6,5])
            sage: pR=P([3,1,2,5,6,4])
            sage: G=MySubgroup(o2=pS,o3=pR)       
            sage: G._get_uid()
            '[2_1_4_3_6_5]-[3_1_2_5_6_4]'        
        """
        # If we have a Congruence subgroup it is (more or less)  easy
        s=""
        N=self.generalised_level()
        s=str(self.permS.list())+"-"
        s=s+str(self.permR.list())
        s=s.replace(", ","_")
        return s 


        
    def __reduce__(self):
        r"""
        Used for pickling self.
        
        EXAMPLES::


            sage: G=MySubgroup(Gamma0(5))

            Not implmented!
        """
        data = self.__dict__
        return (self.__class__, (self.permS,self.permR,self._verbose,data))            



    def __cmp__(self,other):
        r""" Compare self to other.

        EXAMPLES::

            sage: G=MySubgroup(Gamma0(5))
            sage: G <> Gamma0(5)
            False
            sage: GG=MySubgroup(None,G.permS,G.permR)
            sage: GG == G

        
        """
        if not isinstance(other,MySubgroup_class):
            return -1
        return super(MySubgroup_class,self).__cmp__(super(MySubgroup_class,self))

    def __ne__(self,G):
        return not self.__eq__(G)
        
        
    def __eq__(self,G):
        r"""
        Test if G is equal to self.

        EXAMPLES::

        
            sage: G=MySubgroup(Gamma0(5))
            sage: G==Gamma0(5)
            False
            sage: GG=MySubgroup(None,G.permS,G.permR)
            sage: GG == G
            True
        """
        ## Perform some fast tests first
        if not isinstance(G,ArithmeticSubgroup):
            return False
        if G.index() <> self.index():
            return False
        if G.generalised_level() <> self.generalised_level():
            return False
        return super(MySubgroup_class,self).__cmp__(G) == 0
        # if not isinstance(G,MySubgroup_class):
        #     S=G.as_permutation_group().S2()
        #     R=G.as_permutation_group().S3()
        #     t,p= are_mod1_equivalent(self.permR,self.permS,R,S)  
        # else:
        #     t,p= are_mod1_equivalent(self.permR,self.permS,G.permR,G.permS)
        # if t==1:
        #     return True
        # return False
        #return self.is_subgroup(G) and G.is_subgroup(self)

    def relabel(self,inplace=True,label_on='R'):
        r"""

        INPUT:
        
        - inplace -- bool
        - label_on    -- 'R', 'S' or 'T (default 'R').
                      Relabel so that either permR,permS or permT are ordered, e.g.
                      R = (1 2 3)(4 5 6)(7 8 9)
                      S = (1)(2 3)(4 5)
                      T = (1 2 3 4 5 6 7 8 9)
        """
        if inplace==False:
            G = self.__class__(o2=self.permS,o3=self.permR,verbose=0)
            G.relabel(inplace=True,label_on=label_on)
            return G
        if label_on == 'R':
            super(MySubgroup_class,self).relabel(inplace=True)
            self.permS=MyPermutation([x+1 for x in self._S2])
            self.permR=MyPermutation([x+1 for x in self._S3])
            ## Relabel the rest as well
            self.permT = self.permS*self.permR
            self.permP = self.permT*self.permS*self.permT
        else:
            if label_on=='S':
                Pold = self.permS
                cycle_lens = self.permS.cycle_type()
            elif label_on=='T':
                cycle_lens = self.permT.cycle_type()
                Pold = self.permT            
            Pnew = []; j=1
            for i in cycle_lens:
                Pnew.append([x for x in range(j,i+j)])
                j = i+j
            if self._verbose>0:
                print "Pnew=",Pnew
            Pnew = MyPermutation(Pnew)
            if self._verbose>0:
                print "Pnew=",Pnew
            t,p = Pold.is_conjugate_to(Pnew,ret_perm=1)
            if t<>1:
                raise ArithmeticError,"Could not conjugate {0} to {1}!".format(Pold,Pnew)
            if self._verbose>0:
                print "p=",p
            self.permS = self.permS.conjugate(p)
            self.permR = self.permR.conjugate(p)
            self.permT = self.permS*self.permR
            self.permP = self.permT*self.permS*self.permT            
        if self._verbose>0:
            print "Snew=",self.permS
            print "Rnew=",self.permR
            print "Tnew=",self.permT            
        ## Finally we have to make sure that the coset representatives have the same order
        ## The simplest way to do this is simply to reset the lists.
        self._coset_reps_v0 = None
        self._coset_reps_v1 = None
        self._coset_reps_v2 = None
        # if self._coset_reps_v0<>None or self._coset_reps_v1<>None:
        #     reps_values = [self.permutation_action(x)(1) for x in self.coset_reps()]
        #     reps_indices= [reps_values.index(x+1) for x in range(self.index())]
        #     if self._coset_reps_v0<>None:
        #         reps = []
        #         for i in range(self.index()):
        #             reps.append(self._coset_reps_v0[reps_indices[i]])
        #         self._coset_reps_v0 = reps
        #     if self._coset_reps_v1<>None:
        #         reps = []
        #         for i in range(self.index()):
        #             reps.append(self._coset_reps_v1[reps_indices[i]])
        #         self._coset_reps_v1 = reps
                
        
    def init_group_from_permutations(self,o2,o3):
        r"""
        Initialize the group using the two permutations of order 2 and 3.
        """
        if not isinstance(o2,MyPermutation):
            o2 = MyPermutation(o2)
        if not isinstance(o3,MyPermutation):
            o3 = MyPermutation(o3)
        o2.set_rep(3)
        o3.set_rep(3)
        if self._verbose>0:
            print "in init_from_perm"
            print "o2=",o2
            print "o3=",o3
            if self._verbose>1:
                print "dict=",self.__dict__
        if isinstance(o2,MyPermutation):            
            self.permS=o2
        else:
            self.permS=MyPermutation(o2)
        if isinstance(o3,MyPermutation):            
            self.permR=o3
        else:
            self.permR=MyPermutation(o3)
        self._index = self.permR.N()
        ## The generators of EvenArithmeticSubgroup_Permutation is corresponding to
        ## S2 = S, S3 = ZST^-1, L=T, R=Z*ST^-1*S where Z = S^2 = [-1,0,0,-1]
        ## Recall that I assume my input is o3 = S*T
        self.permT = self.permS*self.permR
        self.permP = self.permT*self.permS*self.permT
        s2 = [i-1 for i in self.permS.list()]
        s3 = [i-1 for i in self.permR.inverse().conjugate(self.permS).list()]
        l  = [i-1 for i in self.permT.list()]
        r  = [i-1 for i in self.permT.conjugate(self.permS).inverse().list()]
        super(MySubgroup_class,self).__init__(s2,s3,l,r)
        if self._is_congruence == None:
            self._is_congruence = super(MySubgroup_class,self).is_congruence()
#         if self._is_congruence==True:
#             #print "Adding level!"
# #            setattr(MySubgroup_class,'level', types.MethodType(level,self,MySubgroup_class))
#             self.level = types.MethodType(level,self,MySubgroup_class)
#             #self._level = self.level()
#             #print "level=",self._level
        self.get_data_from_group()       

    def init_group_from_dict(self,data,**kwds):
        r"""
        Initalize self from a dictionary.
        """
        if self._verbose>0:
            print "in init_from_dict"
        for key in data:
            if data[key]<>None:
                self.__dict__[key] = data[key]
        for key in kwds:
            if kwds[key]<>None:
                self.__dict__[key] = kwds[key]            

    def get_data_from_group(self):
        if self._verbose>0:
            print "in Get_data_from_Group"
        if self._level <> None:
            self._generalised_level = self._level
        else:
            self._generalised_level = super(MySubgroup_class,self).generalised_level()
            if self._is_congruence:
                self._level = self._generalised_level
            else:
                self._level = None            
        #self._coset_reps_list=copy(self._coset_reps)
        self._test_consistency_perm(self.permS,self.permR)
        self._nu2=num_fixed(self.permS.list())
        self._nu3=num_fixed(self.permR.list())
        
        ## Get information about cusps and vertices
        l=self._get_all_cusp_data(self.coset_reps())
        if self._verbose>0:
            print "coset_reps=",self.coset_reps()
        self._vertices,self._vertex_data,self._cusps,self._cusp_data=l        
        self._nvertices=len(self._vertices)
        #self._vertex_widths=list()
        #self._vertex_maps=list()
        #self._cusp_maps=list()
        for i in range(len(self._vertices)):
            wi = self._cusp_data[self._vertex_data[i]['cusp']]['width']
            self._vertex_widths.append(wi)
            N=self._cusp_data[self._vertex_data[i]['cusp']]['normalizer']
            N = SL2Z_elt(N[0],N[1],N[2],N[3])
            U = self._vertex_data[i]['cusp_map']
            self._cusp_maps.append(U) #[U[0,0],U[0,1],U[1,0],U[1,1]])
            N = N.inverse()*U
            self._vertex_maps.append(N) #[N[0,0],N[0,1],N[1,0],N[1,1]])

        # We might also want to see which cusps are simultaneously
        # symmetrizable with respect to reflection in the imaginary axis 
        #self._symmetrizable_cusp=dict()

        ## Then we chek if the cusps are symmetrizable in the sense that the normalizing maps
        ## are normalizers of the group.
        ## The entries of this dict are pairs: (o,d) where
        ## N^o is in self and has [1,1] element d.
        self._cusp_normalizer_is_normalizer={0: (1,1) }# The first map is just the identity
        for j in range(1,self.ncusps()):
            d=self.cusp_normalizer_is_normalizer(j,1)
        if self._verbose>1:
            print "inited from group, dict=",self.__dict__

    def index(self):
        if self._index == None:
            self._index = self.permR.N() 
        return self._index
       
    def genus(self):
        r""" Genus of self

        EXAMPLES::


            sage: G=MySubgroup(Gamma0(5));
            sage: G.genus()
            5

        """
        if self._genus == None:
            self._genus=1 +QQ(self.index() - 6*self.ncusps()-3*self.nu2()-4 *self.nu3())/QQ(12)
        return int(self._genus)

    def nu2(self):
        r""" Number of elliptic fixed points of order 2 of self.

        EXAMPLES::


            sage: G=MySubgroup(Gamma0(5));
            sage: G.genus()
            5

        """
        if self._nu2==None:
            self._nu2 = self.permS.num_fixed()
        return self._nu2

    def nu3(self):
        r""" Number of elliptic fixed points of order 3 of self.
        EXAMPLES::


            sage: G=MySubgroup(Gamma0(5));
            sage: G.genus()
            5

        """
        if self._nu3==None:
            self._nu3 = self.permR.num_fixed()
        return self._nu3


    def ncusps(self):
        r"""
        Return the number of cusps of self.
        """
        if self._ncusps == None:
            self._ncusps = len(self.permT.cycle_type())
        return self._ncusps
    
    def signature(self):
        r"""
        Returns the signature of self: (index,h,nu2,nu3,g).
        
        """
        if self._signature == None:
            self._signature =  (self.index(),self.ncusps(),self.nu2(),self.nu3(),self.genus())
        return self._signature

    def reflected_group(self):
        r"""
        Return G^* = JGJ where J:z -> -\bar(z)
        """
        if self._reflected_group == None:
            pS = self.permS
            pR = self.permS*self.permR**2*self.permS
            self._reflected_group = MySubgroup(pS,pR)
        return self._reflected_group

    ###
    ## More advanced functions
    ###

    def is_Gamma0(self,test=0):
        if self._is_Gamma0==None:
            self._is_Gamma0=False
            if self.is_congruence():
                N = self.generalised_level() # Is the level in this case
                ## TODO: See which is the quickest test
                G = Gamma0(N)
                if G.index() <> self.index(): # If self is a subgroup of Gamma0(N) and the index is equal...
                    self._is_Gamma0=False
                else:                    
                    self._is_Gamma0=True
                    if test==1:
                        for g in self.gens():
                            if g.c() % N <> 0:
                                self._is_Gamma0 = False
                                break
                    else:
                        if not self.is_subgroup(G):
                            self._is_Gamma0=False

        return self._is_Gamma0
    
    def is_symmetric(self,ret_map=0,recompute=False,force_check=False,verbose=0):
        r"""
        Check if self has a reflectional symmetry, i.e. check that G^* is conjugate to G

        If it is then we set the internal variables A=self._symmetry_map s.t. AG^*A^-1=G
        and a permutation p=self._sym_perm s.t. p^-1*(S,R)p ~ (S*,R*) mod 1
        i.e. so that the pairs represent the same group.
        
        """
        if self._is_symmetric <>None and recompute==False:
            if ret_map==1:
                return self._is_symmetric,self._symmetry_map
            else:
                return self._is_symmetric
        if self.is_Gamma0():
            self._is_symmetric = True
            self._sym_perm = MyPermutation(length=self.index())            
        else:
            self._is_symmetric = False
            self._sym_perm = None
            ## We first check for symmetries of type IIa
            self._is_symmetric = self._has_symmetry_type_IIa(verbose)
            if self._is_symmetric:
                self._sym_perm = MyPermutation(length=self.index())
                A = self._symmetry_type_IIa[0] # = AA*J
                self._symmetry_map =  SL2Z_elt(-A.a(),A.b(),-A.c(),A.d())
            else:
                ## Check more complicated symmetry of type IIb:
                self._is_symmetric = self._has_symmetry_type_IIa(verbose)
                if self._is_symmetric:
                    self._sym_perm = MyPermutation(length=self.index())
                    A = self._symmetry_type_IIa[0][2] 
                    ## A^-1 T^n G* T^-n A = G
                    self._symmetry_map =  SL2Z_elt(-A.a(),A.b(),-A.c(),A.d())                
            if not self._sym_perm and force_check==True:
                if verbose>0:
                    print "Checking symmetry with conjugation of general maps!"
                # Check if we are symmetric with some other map
                # Since this might take long time we only do this if explicitly told to... 
                pS = self.permS
                pR = self.permS*self.permR**2*self.permS
                ## Now have to see if (pS,pR) is conjugate to (self.permS,self.permR)
                ## by a permutation fixing 1
                t,p = are_mod1_equivalent(self.permS,self.permR,pS,pR,verbose=verbose)
                if t==1:
                    self._is_symmetric = True
                    self._sym_perm = p
                else:
                    self._is_symmetric = False
                    self._sym_perm = MyPermutation(length=self.index())
            if self._sym_perm == None and self._symmetry_map<>None:
                self._sym_perm = self.permutation_action(self._symmetry_map)
            
        if ret_map==1:
            return self._is_symmetric,self._symmetry_map
        else:
            return self._is_symmetric


    def has_translational_symmetry(self,verbose=0):
        r"""
        Check if self has a symmetry of the form  or T^k G T^-k = G with T^2k in G
        and T^k not in G.
        
        OUTPUT:

        - 'k' -- integer, the smallest k>0 such that T^kGT^-k=G and T^2k in G.
        
        """
        if self._translational_symmetry == None:
            if self.is_Gamma0():
                self._translational_symmetry = 0  # symmetry z -> -bar(z) given by T^0
            else:
                self._translational_symmetry = -1 
                if self._verbose>0:
                    print "Checking symmetry with conjugation of T^n!"
                for n in range(1,self._cusp_data[0]['width']):
                    t = [ SL2Z_elt(x.a()+n*x.c(),x.b()+n*(x.d()-x.a())-n*n*x.c(),x.c(),x.d()-n*x.c()) in self for x in self.gens()].count(False)
                    if t == 0 and SL2Z_elt(1,2*n,0,1) in self:
                        self._translational_symmetry = n  # symmetry z -> n-bar(z) given by J*T^n
                        break
        return self._translational_symmetry


    ## Locate symmetries of different types.
    def _has_symmetry_type_Ia(self,verbose=0):
        r"""
        Check if self has a symmetry of the form T**-j G T^j = G 

        OUTPUT:

        - k -- integer >=0 if we have T^k G T^-k=G and k=-1 if no such exist.
        
        """
        if self._checked_symmetry_type_Ia == 1:
            return self._symmetry_type_Ia <> []
        for n in range(1,self._cusp_data[0]['width']+1):
            Tn = SL2Z_elt(1,n,0,1)
            Tni = SL2Z_elt(1,-n,0,1)
            t = [ Tn*x*Tni in self for x in self.generators_as_slz_elts()].count(False)
            if t==0 and Tn not in self:
                self._symmetry_type_Ia.append(SL2Z_elt(1,n,0,1))
        self._checked_symmetry_type_Ia = 1
        return self._symmetry_type_Ia <> []

    def _has_symmetry_type_Ib(self,verbose=0):
        r"""
        Check if self has a symmetry of the form T**j*A G A**-1*T**-jj = G 
        where A is a cusp normalizer of G.
        OUTPUT:

        - 't' -- Boolean
        
        """
        if self._checked_symmetry_type_Ib == 1:
            return self._symmetry_type_Ib <> []
        Gs = self.generators_as_slz_elts()
        for n in range(0,self._cusp_data[0]['width']):
            Tn = SL2Z_elt(1,n,0,1)
            for j in range(1,self.ncusps()):
                A = self.cusp_normalizer(j); Ai = A.inverse()
                for AA in [A,Ai]:
                    if AA==Ai: j=-j
                    B = Tn*AA
                    Bi = B.inverse()
                    if verbose>0:
                        print "Check B=T^{0}A_{1} \t\t=\t {2}".format(n,j,B)
                    t = [ B*x*Bi in self for x in Gs].count(False)
                    if verbose>0:
                        print [ B*x*Bi in self for x in Gs]
                    if t==0 and B not in self:
                        self._symmetry_type_IIb.append((j,n,SL2Z_elt(B.a(),B.b(),B.c(),B.d())))
                    # Checking after conjugation by J:
                    B = SL2Z_elt(-B.a(),B.b(),B.c(),-B.d())
                    Bi = B.inverse()
                    if verbose>0:
                        print "Check B=T^{0}A^-1_{1}\t\t =\t {2}".format(n,j,B)                    
                    t = [ B*x*Bi in self for x in Gs].count(False)
                    if verbose>0:
                        print "tests:",t
                    if t==0 and B not in self:
                        s = (j,-n,SL2Z_elt(B.a(),B.b(),B.c(),B.d()))
                        if s not in self._symmetry_type_IIb:
                            self._symmetry_type_IIb.append(s)
        self._checked_symmetry_type_Ib = 1
        return self._symmetry_type_Ib <> []
    
    def _has_symmetry_type_IIa(self,verbose=0):
        r"""
        Check if self has a symmetry of the form T**-j G T^j = JGJ 

        OUTPUT:

        - k -- integer >=0 if we have T^k JGJ T^-k=G and k=-1 if no such exist.
        
        """
        if self._checked_symmetry_type_IIa == 1:
            return self._symmetry_type_IIa <> []
        Gs = self.reflected_group().generators_as_slz_elts()
        for n in range(0,self._cusp_data[0]['width']):
            Tn = SL2Z_elt(1,n,0,1)
            Tni = SL2Z_elt(1,-n,0,1)
            t = [ Tn*x*Tni in self for x in Gs].count(False)
            if t==0:
                A = GL2Z_elt(-1,n,0,1)
                self._symmetry_type_IIa.append(A)
        self._checked_symmetry_type_IIa = 1
        return self._symmetry_type_IIa <> []
            
    def _has_symmetry_type_IIb(self,verbose=0):
        r"""
        Check if self has a symmetry of the form T**-j * A * G * A^-1 * T^j = JGJ 

        OUTPUT:

        - k -- integer >=0 if we have T^k G T^-k=G and k=-1 if no such exist.
        
        """
        if self._checked_symmetry_type_IIb == 1:
            return self._symmetry_type_IIb <> []
        Gs = self.reflected_group().generators_as_slz_elts()
        for j in range(1,self.ncusps()):
            A = self.cusp_normalizer(j)
            for n in range(0,self._cusp_data[0]['width']):
                Tn = SL2Z_elt(1,n,0,1)
                B = Tn*A
                Bi = B.inverse()
                if verbose>0:
                    print "Tn=",Tn
                    print "Check n={0} and j={1} map ={2}".format(n,j,B)
                t = [ B*x*Bi in self for x in Gs].count(False)
                #  print "t0=",t
                if t==0:
                    s = (j,n,GL2Z_elt(-B.a(),B.b(),-B.c(),B.d()))
                    if s not in self._symmetry_type_IIb:
                        self._symmetry_type_IIb.append(s)
                t = [ Bi*x*B in self for x in Gs].count(False)
                # print "t1=",t
                if t==0:
                    s = (-j,n,GL2Z_elt(-Bi.a(),Bi.b(),-Bi.c(),Bi.d()))
                    if s not in self._symmetry_type_IIb:
                        self._symmetry_type_IIb.append(s)
#t = [ SL2Z_elt(x.a()-n*x.c(),-x.b()+n*(x.d()-x.a())+n*n*x.c(),-x.c(),x.d()+n*x.c()) in self for x in self.gens()].count(False)
        self._checked_symmetry_type_IIb = 1
        return self._symmetry_type_IIb <> []


    def has_modular_correspondence(self,ret_map=False,as_2by2=False):
        r"""
        Check if we can find a modular correspondence for self. 


        Note: 

        """
        if self._modular_correspondences == {}:
            self._has_symmetry_type_Ia(); self._has_symmetry_type_Ib()
            self._has_symmetry_type_IIa(); self._has_symmetry_type_IIb()        
            l1a = self._symmetry_type_Ia; l1a = self._symmetry_type_Ib
            l2a = self._symmetry_type_IIa; l2b = self._symmetry_type_IIb        
            ## Check which (if any) of these symmetries preserve cusps
            l1a = filter(lambda x: self.is_modular_correspondence(x) and self.is_involution(x),  self._symmetry_type_Ia)
            l1b = filter(lambda x: self.is_modular_correspondence(x[2]) and self.is_involution(x[2]) ,  self._symmetry_type_Ib)
            l2a = filter(lambda x: self.is_modular_correspondence(x) and self.is_involution(x),  self._symmetry_type_IIa)
            l2b = filter(lambda x: self.is_modular_correspondence(x[2]) and self.is_involution(x[2]),  self._symmetry_type_IIb)
            self._modular_correspondences = {'tIa':l1a,'tIb':l1b,'tIIa':l2a,'tIIb':l2b}
        return sum(map(len,self._modular_correspondences.values()))

    def modular_correspondence(self,t):
        r"""
        
        """
        if t not in ['tIa','tIIa','tIb','tIIb']:
            raise ValueError,"Need t in 'tIa','tIIa','tIb','tIIb']"
        self.has_modular_correspondence()
        return self._modular_correspondences.get(t)

    def modular_correspondence_matrix(self,t=None):
        self.has_modular_correspondence()
        if t in ['tIa','tIIa']:
            return self._modular_correspondences.get(t,None)
        if t in ['tIb','tIIb']:
            return self._modular_correspondences.get(t,None)[0][2]
        ## Otherwise we take as simple as possible:
        if t == 'all':
            res = []
            for t in ['tIa','tIIa','tIb','tIIb']:            
                m = self.modular_correspondence_matrix(t)            
        for t in ['tIa','tIIa','tIb','tIIb']:            
            m = self.modular_correspondence_matrix(t)
            if m<>[]:
                return m[0]
        
    def modular_correspondence_string(self,t=None):
        self.has_modular_correspondence()
        A = self._modular_correspondences.get('tIa',[])
        B = self._modular_correspondences.get('tIIa',[])
        s = ""
        if len(A)>0:
            A = A[0]
            s = "Conjugaction by T^{0}".format(QQ(A.b()))
        elif len(B)>0:
            B = B[0]
            s = "Reflection in Re(z)={0}".format(QQ(B.b())/QQ(2))
        else:
            l = self._modular_correspondences.get('tIb',[])
            if len(l)>0:
                ia,na,A = l[0]
                if A<>None:    
                    s = "Conjugation by T^{0}\sigma_{1}".format(na,ia)
            l = self._modular_correspondences.get('tIIb',[])                       
            if len(l)>0:
                ib,nb,B = l[0]       
                if B<>None:    
                    s = "Reflection with respect to cusp {0} / shifted by T^{1}".format(ib,nb)
        return s
    

    def _modular_correspondence_desc(self,A,ia=None,na=None):
        r"""
        Get description
        """
        if A.determinant() == 1:
            if A.c()==0 and ia==None:
                return "Conjugaction by T^{0}".format(QQ(A.b()))
            elif ia<>0:
                return "Conjugation by T^{0}\sigma_{1}".format(na,ia)
        else:
            if A.c()==0 and ia==None:
                return "Reflection in Re(z)={0}".format(QQ(A.b())/QQ(2))
            elif ia<>0:
                return "Reflection with respect to cusp {0} / shifted by T^{1}".format(ia,na)
        return ""
    
    def is_modular_correspondence(self,A,verbose=0):
        r"""
        Check if A is a modular correspondence for self. That is, if A is in PGL(2,Z) s.t. A^2 in self
        and A preserves cusp classes.


        """
        try:
            if not hasattr(A,"a"):                
                a=A[0,0]; b=A[0,1]; c=A[1,0]; d=A[1,1]
            else:
                a=A.a(); b=A.b(); c=A.c(); d=A.d()
        except TypeError:
            return False
        #if A**2 not in self: ### This is for an involution not a correspondence
        #    return False
        if self._verbose>0 or verbose>0:
            print A
        for cusp in self.cusps():
            x = cusp.numerator(); y = cusp.denominator()
            Ac = Cusp(a*x+b*y,c*x+d*y)
            if self._verbose>1 or verbose>0:
                print "A({0})={1}".format(cusp,Ac)
            #print Ac,c
            if not self.are_equivalent(Ac,cusp):
                return False
        return True

    def is_involution(self,A):
        return A**2 in self
    
    def cusp_normalizer_is_normalizer(self,j,brute_force=0):
        r"""
        The dictionary self._cusp_normalizer_is_normalizer
        consists of key-value pairs:
         j => (o,d)
         where o is the order of the cusp normalizer of cusp nr. j
         and =0 if the normalizer is not a normalizer of the group (or identity)
         and d is factor which the lower right entry in N*A*N^-1 differs from that of A by
         I.e. chi(NAN^-1)=chi(d)chi(A)
         (should mostly be 1 if everything works as hoped)
        """
        if self._verbose>0:
            print "Checking j:",j
            print "have:",self._cusp_normalizer_is_normalizer
        if self._cusp_normalizer_is_normalizer.has_key(j):
            return self._cusp_normalizer_is_normalizer[j]
        if not self.is_Gamma0():
            self._cusp_normalizer_is_normalizer[j]=(0,0)
            return self._cusp_normalizer_is_normalizer[j]
        l = self._level
        self._cusp_normalizer_is_normalizer[j]=(0,0)
        N0=self._cusp_data[j]['normalizer']
        if self._verbose>0:
            print "N={0}".format(N0)
        a,b,c,d=N0
        w = self._cusp_data[j]['width']
        if self._verbose>0:
            print "w=",w
        if self._cusps[j]==(0,1):
            if a==0 and b*c==-1 and d==0:
                self._cusp_normalizer_is_normalizer[j]=(2,1)
            else:
                warnings.warn("\nIt appears that the normalizer of 0 is not w_N! Have:{0}".format(a,b,c,d))
        elif j>0 and j<self._ncusps:
            ## First see if we actually have an Atkin-Lehner involution directly:
            ## Multiply with the scaling matrix.
            aa=a*w
            cc=c*w            
            if cc==l and (l % aa)==0:
                if self._verbose>0:
                    print "possible A-L invol: Q/N=",aa,"/",cc
                if aa*d-b*cc==aa:
                    self._cusp_normalizer_is_normalizer[j]=(2,d)
                    return self._cusp_normalizer_is_normalizer[j]
            # The other Atkin-Lehner involutions are more tricky...
            if self._verbose>1:
                print "cnc1=",self._cusp_normalizer_is_normalizer[j]                        
            p,q=self._cusps[j]
            if q.divides(l*p):
                Q  = (l*p).divide_knowing_divisible_by(q)
                if self._verbose>0:
                    print "Q=",Q
                ## It is necessary that gcd(Q,N/Q)=1 for an A-L involution to exist
                if Q.divides(l):
                    lQ=l.divide_knowing_divisible_by(Q)
                    if self._verbose>0:
                        print "N/Q=",lQ
                        #print "cnc1=",self._cusp_normalizer_is_normalizer[j]                        
                    if gcd(Q,lQ)==1:
                        if self._verbose>0:
                            print "gcd(Q,N/Q)=1 => possible Atkin-Lehner involution here."
                        fak=lcm(a,Q)  ## We try to see if we can writ e the map as an A-L inv.
                        aa=fak*a; bb=fak*b; cc=fak*c; dd=fak*d
                        if Q.divides(aa) and Q.divides(dd) and l.divides(cc) and (aa*dd-bb*cc)==Q:
                            # We now have a normalizer so we have to find its order.
                            for k in range(2,l):
                                N=N*N0
                                if N in self:
                                    self._cusp_normalizer_is_normalizer[j]=(k,N[1,1])
                                    break #return self._cusp_normalizer_is_normalizer[j]
                            if k>=l-1 and self._cusp_normalizer_is_normalizer[j]==(0,0):
                                warnings.warn("It appears that the normalizer does not have finite order! N={0}".format(N0))
            # We also have a brute force method:
            # print "doing brute force!"
            if self._verbose>1:
                print "cnc2[",j,"]=",self._cusp_normalizer_is_normalizer[j]                        
            if brute_force==1:
                if self.is_normalizer(N0):
                    N = SL2Z_elt(N0[0],N0[1],N0[2],N0[3])
                    for k in range(2,l+1):
                        N=N._mul(N0) #mul_list_maps(N,N)
                        # print "N=",N,type(N)
                        if N in self:
                            self._cusp_normalizer_is_normalizer[j]=(k,N[3])
                            break
                    if k>=l-1 and self._cusp_normalizer_is_normalizer[j]==(0,0):
                        warnings.warn("It appears that the normalizer does not have finite order! j={0}, N={1}. Group given by S={2} and R={3}".format(j,N0,self.permS,self.permR))
        return self._cusp_normalizer_is_normalizer[j]                        


    
    def is_symmetrizable_even_odd(self,j):
        r"""
        Returns 1 if this cusp is symmetrizable in the sense that the normalizing map N=A*rho
        satisfies: JNJ^-1=AN where A is a memeber of  self.
        Note: If self is a G amma_0(l) then A[1,1]==1 mod l

        """
        if j not in self._cusp_data:
            raise ValueError,"{0} is not the index of a cusp! ".format(j)
        if not self._symmetrizable_cusp.has_key(j):
            a,b,c,d=self._cusp_data[j]['normalizer']
            self._symmetrizable_cusp[j]=0
            if self._is_Gamma0:
                if self._level.divides(2*d*c):
                    self._symmetrizable_cusp[j]=1
            elif self.is_symmetric():
                if  [a*d+b*c,-2*a*b,-2*d*c,a*d+b*c] in self:
                    self._symmetrizable_cusp[j]=1
            
        return self._symmetrizable_cusp[j]

    def normalizer_order(self,j):
        r"""
        If cusp number j has normalizer N
        satisfies: JNJ^-1=AN where A is a memeber of  self.
        Note: If self is a G amma_0(l) then A[1,1]==1 mod l

        """
        return self._symmetrizable_cusp[j]
            

    def _get_vertices(self,reps):
        r""" Compute vertices of a fundamental domain corresponding
             to coset reps. given in the list reps.
        INPUT:
        - ''reps'' -- a set of matrices in SL2Z (coset representatives)
        OUTPUT:
        - [vertices,vertex_reps]
          ''vertices''    = list of vertices represented as cusps of self
          ''vertex_reps'' = list of matrices corresponding to the vertices
                          (essentially a reordering of the elements of reps)


        EXAMPLES::

        
            sage: G=MySubgroup(Gamma0(5))
            sage: l=G._coset_reps_v0
            sage: G._get_vertices(l)
            sage: G._get_vertices(l)
            [[Infiniy, 0], {0: [ 0 -1]
            [ 1 -2], Infinity: [1 0]
            [0 1]}]
            sage: P=Permutations(6)
            sage: pS=P([2,1,4,3,6,5])
            sage: pR=P([3,1,2,5,6,4])
            sage: G=MySubgroup(o2=pS,o3=pR)   
            sage: G._get_vertices(G._coset_reps_v0)
            [[Infinity, 0, -1/2], {0: [ 0 -1]
            [ 1  2], Infinity: [1 0]
            [0 1], -1/2: [-1  0]
            [ 2 -1]}]
        """
        v=list()
        vr=dict()
        for A in reps:
            if(A[1 ,0 ]<>0  and v.count(Cusp(A[0 ,0 ],A[1 ,0 ]))==0 ):
                c=Cusp(A[0 ,0 ],A[1 ,0 ])
                v.append(c)
            elif(v.count(Cusp(infinity))==0 ):
                c=Cusp(infinity)
                v.append(c)
            vr[c]=A
        # Reorder so that Infinity and 0 are first
        if(v.count(Cusp(0 ))>0 ):
            v.remove(Cusp(0 )); v.remove(Cusp(infinity))
            v.append(Cusp(0 )); v.append(Cusp(infinity))
        else:
            v.remove(Cusp(infinity))
            v.append(Cusp(infinity))
        v.reverse()
        return [v,vr]



    def permutation_action(self,A,type=1):
        r""" The permutation corresponding to the matrix A
        INPUT:
         - ''A'' Matrix in SL2Z
        OUPUT:
         element in self._S given by the presentation of the group self
         - type: 0 gives PermutationGroup element and 1 gives MyPermutation

        EXAMPLES::


            P=Permutations(6)
            sage: pS=P([2,1,4,3,6,5])
            sage: pR=P([3,1,2,5,6,4])
            sage: G=MySubgroup(o2=pS,o3=pR)
            sage: G.permutation_action(S*T).cycle_tuples()
            [(1, 3, 2), (4, 5, 6)]
            sage: G.permutation_action(S).cycle_tuples()
            [(1, 2), (3, 4), (5, 6)]
            sage: pR.cycle_tuples()
            [(1, 3, 2), (4, 5, 6)]
            sage: pS.cycle_tuples()
            [(1, 2), (3, 4), (5, 6)]
            sage: G=MySubgroup(Gamma0(5))
            sage: G.permutation_action(S).cycle_tuples()
            [(1, 2), (3, 4), (5,), (6,)]
            sage: G.permutation_action(S*T).cycle_tuples()
            [(1, 3, 2), (4, 5, 6)]

        """
        [sg,t0,cf]=factor_matrix_in_sl2z(A) #_in_S_and_T(A)
        if self._verbose>1:
            print "sg=",sg
            print "t0=",t0
            print "cf=",cf
            print "S=",self.permS
            print "T=",self.permT
            print "index=",self._index
            # The sign doesn't matter since we are working in practice only with the projective group PSL(2,Z)
        # We now have A=T^a0*S*T^a1*...*S*T^an
        # Hence
        # h(A)=h(T^an)*h(S)*...*h(T^a1)*h(S)*h(T^a0)
        #print "cf=",cf
        n=len(cf)
        
        if t0==0:
            p=MyPermutation(length=self._index)
        else:
            p=self.permT**t0 #p=ppow(self.permT,cf[0 ])
        # print "p(",0,")=",p.cycles()
        #if self._verbose>1:
        #    print "p=",p
        for j in range(n):
            a=cf[j]
            if a<>0:
                Ta = self.permT**a
                #if self._verbose>1:
                #    print "Ta=",Ta
                p=p*self.permS*Ta
            else:
                p=p*self.permS
        if type==0:
            return SymmetricGroup(self._index)(p.list())
        else:
            return p


    def _get_coset_reps_from_G(self,G,string=False,perm=False):
        r"""
        Compute a better/nicer list of right coset representatives [V_j]
        i.e. SL2Z = \cup G V_j
        Use this routine for known congruence subgroups.

        EXAMPLES::


            sage: G=MySubgroup(Gamma0(5))
            sage: G._get_coset_reps_from_G(Gamma0(5))
            [[1 0]
            [0 1], [ 0 -1]
            [ 1  0], [ 0 -1]
            [ 1  1], [ 0 -1]
            [ 1 -1], [ 0 -1]
            [ 1  2], [ 0 -1]
            [ 1 -2]]
    
        """
        if self._verbose>0:
            print "In get_coset_reps from G"
        if self.is_Gamma0():
            return self._get_coset_reps_from_Gamma0N()
        cl=list()
        S=[0,-1,1,0]
        T=[1,1,0,1]
        #S,T=SL2Z.gens()
        lvl=12 #G.generalised_level()
        # Start with identity rep.
        cl.append([1 ,0 ,0 ,1 ])
        if self.is_Gamma0() and self._level>1:
            cl.append(S)
        elif not S in G:
            cl.append(S)
        # If the original group is given as a Gamma0 then
        # the reps are not the one we want
        # I.e. we like to have a fundamental domain in
        # -1/2 <=x <= 1/2 for Gamma0, Gamma1, Gamma
        if self.is_Gamma0():
            N=lvl #self._level #G.level()
            if (N % 2) == 0:
                #reprange=range(-N/2+1,N/2+1)
                reprange=range(0,N)
            else:
                reprange=range(0,N)
                #reprange=range(-(N-1)/2,(N-1)/2+1)
            for j in reprange:
                if j==0:
                    continue
                cl.append([0 ,-1 ,1 ,j])
        else:
            #print "HERE!"
            for j in range(1 , floor(0.5*lvl) + 2):
                for ep in [1 ,-1]:
                    if(len(cl)>=self._index):
                        break
                    # The ones about 0 are all of this form
                    A=[0 ,-1 ,1 ,ep*j]
                    # just make sure they are inequivalent
                    try:
                        for V in cl:
                            t1=A<>V
                            #A*V^-1=[a,b]*[ v11 -v01]
                            #       [c d] [-v10 v00]
                            if self.is_Gamma0():
                                t2= (c*v11-d*v10) % lvl == 0 
                            else:
                                t2 = SL2Z_elt(a*v11-b*v10,-v01*a+v00*b,c*v11-d*v10,-c*v01+d*v00) in G
                            if t1 and t2:
                                #if((A<>V and A*V**-1  in G) or cl.count(A)>0 ):
                                raise StopIteration()
                        cl.append(A)
                    except StopIteration:
                        pass
        # We now addd the rest of the "flips" of these reps.
        # So that we end up with a connected domain
        i=1 
        Ti = [1,-1,0,1]
        while(True):
            lold=len(cl)
            for V in cl:
                for A in [S,T,Ti]:
                    #B=V*A
                    B=mul_list_maps(V,A)
                    try:
                        for W in cl:
                            tmp = mul_list_maps(B,W,inv=2)
                            if self.is_Gamma0():
                                t1= (tmp[2] % lvl)==0 
                            else:
                                #print "HERE!22"
                                t1= mul_list_maps(B,W,inv=2) in G
                            #if( (B*W**-1  in G) or cl.count(B)>0 ):
                            if t1:
                                raise StopIteration()
                        cl.append(B)
                    except StopIteration:
                        pass
            if(len(cl)>=self._index or lold>=len(cl)):
                # If we either did not addd anything or if we addded enough
                # we exit
                break
        # If we missed something (which is unlikely)        
        if len(cl)<>self._index:
            print "cl=",cl
            raise ValueError,"Problem getting coset reps! Need %s and got %s" %(self._index,len(cl))
        return cl

    def _get_coset_reps_from_Gamma0N(self):
        if self._verbose>0:
            print "In get_coset_reps from Gamma0(N)"
        cl=list()
        S=SL2Z_elt(0,-1,1,0)
        T=SL2Z_elt(1,1,0,1)
        lvl=self.level()
        cl.append(SL2Z_elt(1 ,0 ,0 ,1 ))
        if self._index==1:
            return cl
        cl.append(S)
        # If the original group is given as a Gamma0 then
        # the reps are not the one we want
        # I.e. we like to have a fundamental domain in
        # -1/2 <=x <= 1/2 for Gamma0, Gamma1, Gamma
        N=lvl
        if (N % 2) == 0:
            reprange=range(-N/2+1,N/2+1)
            #reprange=range(N)
        else:
            reprange=range(-(N-1)/2,(N-1)/2+1)
        decomp_seq=[[0,0],[0]]

        for j in reprange:
            if j==0:
                continue
            cl.append(SL2Z_elt(0 ,-1 ,1 ,j))
            tmp=[0]
            if j>0:
                for i in range(j):
                    tmp.append(1)
            else:
                for i in range(-j):
                    tmp.append(-1)
            decomp_seq.append(tmp)
        # We now addd the rest of the "flips" of these reps.
        # So that we end up with a connected domain
        i=1 
        Ti = SL2Z_elt(1,-1,0,1)
        gens=[S,T,Ti]
        last = -1
        while(True):
            clold=cl
            lold=len(clold)
            for i in range(1,lold):
                V=clold[i]
                #for V in clold:
                last=decomp_seq[i][-1]
                if self._verbose>0:
                    print "V=",V
                for j in range(2):
                    if last==0 and j==0:
                        continue
                    if last==1 and j==2:
                        continue
                    if last==2 and j==1:
                        continue
                    
                    A=gens[j]
                    B=V*A
                    ## We are getting new elements by multiplying the
                    ## old with S,T or T^-1
                    if self._verbose>0:
                        print "V=",V
                        print "A=",A
                        print "B=",B
                    #B=mul_list_maps(V,A)
                    try:
                        for W in cl:
                            #tmp = mul_list_maps(B,W,inv=2)
                            tmp = B._mul(W,inv=2)
                            #print "tmp[2]=",tmp[2]
                            if tmp[2] % lvl==0:
                                raise StopIteration()
                        cl.append(B)
                        tmpl = decomp_seq[i]
                        tmpl.append(j)
                        decomp_seq.append(tmpl)
                        #last=j
                    except StopIteration:
                        pass
            if len(cl)>=self._index or lold>=len(cl):
                # If we either did not addd anything or if we addded enough
                # we exit
                break
        # If we missed something (which is unlikely)        
        if self._verbose>0:
            print "cosets=",cl
        if len(cl)<>self._index:
            print "cl=",cl
            raise ValueError,"Problem getting coset reps! Need %s and got %s" %(self._index,len(cl))
        return cl
        
    
    def _get_coset_reps_from_G_2(self,G,string=False,perm=False):
        r"""
        Compute a better/nicer list of right coset representatives [V_j]
        i.e. SL2Z = \cup G V_j
        Use this routine for known congruence subgroups.

        EXAMPLES::


            sage: G=MySubgroup(Gamma0(5))
            sage: G._get_coset_reps_from_G(Gamma0(5))
            [[1 0]
            [0 1], [ 0 -1]
            [ 1  0], [ 0 -1]
            [ 1  1], [ 0 -1]
            [ 1 -1], [ 0 -1]
            [ 1  2], [ 0 -1]
            [ 1 -2]]
    
        """
        if self._verbose>0:
            print "In get_coset_reps from G_2"
        cl=list()
        lvl=G.generalised_level()
        # Start with identity rep.
        S = SL2Z_elt(0, -1 ,1 ,0 )
        T = SL2Z_elt(1, 1 ,0 ,1 )
        cl.append(SL2Z_elt(1 ,0 ,0 ,1 ))
        if(not S in G):
            cl.append(S)
        # If the original group is given as a Gamma0 then
        # the reps are not the one we want
        # I.e. we like to have a fundamental domain in
        # -1/2 <=x <= 1/2 for Gamma0, Gamma1, Gamma
        for j in range(1 , floor(0.5*lvl) + 2):
            for ep in [1 ,-1 ]:
                if(len(cl)>=self._index):
                    break
                # The ones about 0 are all of this form
                A=SL2Z_elt(0 ,-1 ,1 ,ep*j)
                # just make sure they are inequivalent
                try:
                    for V in cl:
                        if((A<>V and A*V**-1  in G) or cl.count(A)>0 ):
                            raise StopIteration()
                    cl.append(A)
                except StopIteration:
                    pass
        # We now addd the rest of the "flips" of these reps.
        # So that we end up with a connected domain
        i=1 
        while(True):
            lold=len(cl)
            for V in cl:
                for A in [S,T,T**-1 ]:
                    B=V*A
                    try:
                        for W in cl:
                            if( (B*W**-1  in G) or cl.count(B)>0 ):
                                raise StopIteration()
                        cl.append(B)
                    except StopIteration:
                        pass
            if(len(cl)>=self._index or lold>=len(cl)):
                # If we either did not addd anything or if we addded enough
                # we exit
                break
        # If we missed something (which is unlikely)        
        if(len(cl)<>self._index):
            print "cl=",cl
            raise ValueError,"Problem getting coset reps! Need %s and got %s" %(self._index,len(cl))
        return cl


    
    def permutation_coset_rep(self,j):
        r"""
        Return the permutation corresponding to coset-representative nr. j
        """
        if not self._coset_rep_perms.has_key(j):
            perm = self.permutation_action(self.coset_reps()[j])
            self._coset_rep_perms[j]=perm
        return self._coset_rep_perms[j]


    def string_coset_rep(self,j):
        r"""
        Return the permutation corresponding to coset-representative nr. j
        """
        if not self._coset_rep_strings.has_key(j):
            s = self.element_as_string(self.coset_reps()[j])
            self._coset_rep_strings[j]=s
        return self._coset_rep_strings[j]


    def element_as_string(self,A):
        r"""
        Writes A as a string in generators S and T
        """
        i=0
        z,t0,cf=factor_matrix_in_sl2z(A)
        if z==-1:
            s='-'
        else:
            s=''
        if t0<>0:
            s=s+"T^%" %t0
        for n in range(len(cf)):
            s=s+"ST^{%}" % cf[n]
        return s


    def _get_coset_reps_from_perms(self):
        r"""
        Compute a better/nicer list of right coset representatives
        i.e. SL2Z = \cup G V_j

        INPUT:
        - self 

        OUTPUT:
        - list of (right) coset-representatives of the group given by pS and pR
          the rep V[j] has the property that self.permutation_action(V[j])=j


        EXAMPLES::

            sage: P=SymmetricGroup(7)
            sage: pS=P([1,3,2,5,4,7,6]); pS
            (2,3)(4,5)(6,7)
            sage: pR=P([3,2,4,1,6,7,5]); pR
            (1,3,4)(5,6,7)
            sage: G=MySubgroup(o2=pS,o3=pR)
            sage: G._get_coset_reps_from_perms(G.permS,G.permR,G.permT)
            [[1 0]
            [0 1], [1 2]
            [0 1], [1 1]
            [0 1], [1 3]
            [0 1], [1 5]
            [0 1], [1 4]
            [0 1], [ 4 -1]
            [ 1  0]]

        # Testing that the reps really are correct coset reps

            sage: l=G._get_coset_reps_from_perms(G.permS,G.permR,G.permT)
            sage: for V in l:
            ....:     print G.permutation_action(V)
            ....:
            ()
            (1,2,6)(3,4,5)
            (1,3,2,4,6,5)
            (1,4)(2,5)(3,6)
            (1,5,6,4,2,3)
            (1,6,2)(3,5,4)
            (1,7,6,3,4,2)

        
        """
        pS = self.permS; pR=self.permR; pT=self.permT
        ix=len(pR.list())
        T=SL2Z_elt(1 ,1 ,0 ,1 )
        S=SL2Z_elt(0 ,-1 ,1 ,0 )
        Id=SL2Z_elt(1 ,0 ,0 ,1 )
        R=S*T
        coset_reps=dict()
        coset_reps[1]=Id
        if self._verbose>0:
            print "coset_reps=",coset_reps,len(coset_reps)
            print "T=",T,pT.cycle_tuples()
            print "S=",S,pS.cycle_tuples()
            print "R=",R,pR.cycle_tuples()
        cycT=pT.cycle_tuples()
        next_cycle = cycT[0]
        new_index = 0
        got_cycles=[]
        old_map = Id
        for cyi in range(len(cycT)):
            cy = next_cycle
            r = len(cy)
            i=pT(cy[new_index])
            if self._verbose>0:
                print "cy=",cy
                print "new_index=",new_index
                print "i=pT(cy[0])=",i
                print "cy[new_index]=",cy[new_index]
            # adding the rest of the cusp
            if i<>cy[new_index]:
                for j in range(r):
                    if self._verbose>0:
                        print "i=",i
                        print "j=",j
                    if j==new_index:
                        continue
                    k = (j - new_index)
                    if k<=r/2:
                        #_add_unique(coset_reps,cy[j],old_map*T**k)
                        coset_reps[cy[j]]=old_map*T**k
                    else:
                        #_add_unique(coset_reps,cy[j],old_map*T**(k-r))
                        coset_reps[cy[j]]=old_map*T**(k-r)                    
                    if self._verbose>0:
                        print "k=",k
                        print "coset_reps[",cy[j],"]=",coset_reps[cy[j]]
            got_cycles.append(cycT.index(cy))
            # we have now added all translate inside the same cusp 
            # and we should see if we can connect to anther cusp
            # if there is any left
            if self._verbose>0:
                print "cyi=",cyi
                print "len(cycT)-1=",len(cycT)-1
            if cyi>= len(cycT)-1:
                if self._verbose>0:
                    print "break!"
                break
            # otherwise we use the order two element to connect the next cycle to one of the previous ones
            # since (S,T) are transitive this must be the case.
            old_map = Id
            if self._verbose>0:
                print "got_cycles=",got_cycles
            try:
                for cyii in range(len(cycT)):
                    if cyii in got_cycles:
                        ## If we already treated this cycle
                        continue
                    next_cycle = cycT[cyii]
                    if self._verbose>0:
                        print "next_cycle=",next_cycle
                    for cyj in range(len(cycT)):
                        if self._verbose>0:
                            print "cyj=",cyj
                        if cyj not in got_cycles:
                            # We can only use cycles which are in the list
                            continue
                        cy=cycT[cyj]
                        if self._verbose>0:
                            print "check with cy=",cy
                        for i in cy:
                            j = pS(i)
                            if j in next_cycle:
                                # we have connected the cycles
                                old_map = coset_reps[i]*S # this is the connecting map and we may as well add it
                                #_add_unique(coset_reps,j,old_map)
                                coset_reps[j]=old_map
                                new_index = next_cycle.index(j)
                                if self._verbose>0:
                                    print "connecting: S(",i,")=",j," in new cycle!"
                                    print "ix=",new_index
                                    print "old_map = ",old_map
                                    print "next_cycle=",next_cycle
                                raise StopIteration()
            except StopIteration:
                pass
            if old_map == Id:
                raise ValueError,"Problem getting coset reps! Could not connect %s using  %s" %(cycT,pS)
            if self._verbose>0:
                print "in this step:"
                for j in coset_reps.keys():
                    if(coset_reps[j]<>Id):
                        print "V(",j,")=",coset_reps[j] 
                    else:
                        print "V(",j,")=Id"
        # By construction none of the coset-reps are in self and h(V_j)=j so they are all independent
        # But to make sure we got all we count the keys
        if coset_reps.keys() <> range(1,self._index+1):
            print "ix=",ix
            print "cl=",coset_reps
            raise ValueError,"Problem getting coset reps! Need %s and got %s" %(self._index,len(coset_reps))
        res  = list()
        for i in range(ix):
            res.append(coset_reps[i+1])
        return res





    def are_equivalent(self,x,y,trans=False):
        r"""
        Check whether two cusps are equivalent with respect to self

        The algorithm is just the same as in sage exxcept that we use SL2Z_elt to increase speed.
        EXAMPLES::


            sage: G=MySubgroup(Gamma0(5))
            sage: G.are_equivalent(Cusp(1),Cusp(infinity))
            False
            sage: G.are_equivalent(Cusp(1),Cusp(0))
            True
            
        """
        x = Cusp(x)
        y = Cusp(y)
        if not trans:
            try:
                xr = self.reduce_cusp(x)
                yr = self.reduce_cusp(y)
                if xr != yr:
                    return False
                if xr == yr:
                    return True
            except NotImplementedError:
                pass
        vx = lift_to_sl2z(x.numerator(),x.denominator(), 0)
        dx = SL2Z_elt(vx[2], -vx[0], vx[3], -vx[1])
        vy = lift_to_sl2z(y.numerator(),y.denominator(), 0)
        dy = SL2Z_elt(vy[2], -vy[0], vy[3], -vy[1])
        for i in range(self.index()):
            # Note that the width of any cusp is bounded above by the index of self.
            # If self is congruence, then the level of self is a much better bound, but
            # this method is written to work with non-congruence subgroups as well,
            t = dy * SL2Z_elt(1,i,0,1)*dx.inverse()
            if t in self:
                if trans:
                    return t
                else:
                    return True
            else:
                t = dy * SL2Z_elt(-1,-i,0,-1) * dx.inverse() 
                if (self.is_odd() and t in self):
                    if trans:
                        return t
                    else:
                        return True
        return False

            
    def _get_cusps(self,l):
        r""" Compute a list of inequivalent cusps from the list of vertices

        EXAMPLES::


            sage: S=SymmetricGroup(6)
            sage: pS=S([2,1,4,3,6,5])
            sage: pR=S([3,1,2,5,6,4])
            sage: G=MySubgroup(o2=pS,o3=pR)   
            sage: l=G._get_vertices(G._coset_reps_v0)[0];l
            [Infinity, 0, -1/2]
            sage: G._get_cusps(l)
        
        """
        lc=list()
        for p in l:
            are_eq=False
            for c in lc:
                if self.are_equivalent(p,c):
                    are_eq=True
                    continue
            if(not are_eq):
                lc.append(Cusp(p))
        if(lc.count(Cusp(infinity))==1  and lc.count(0 )==1 ):
            lc.remove(Cusp(infinity))
            lc.remove(0 )
            lc.append(Cusp(0 ))
            lc.append(Cusp(infinity))
        else:
            lc.remove(Cusp(infinity))
            lc.append(Cusp(infinity))
        lc.reverse()        
        return lc


    def is_normalizer(self,A,verbose=0):
        if hasattr(A,"acton"):
            A = matrix(ZZ,2,2,list(A))
        for g in self.gens():
            gg = matrix(ZZ,2,2,list(g))
            AA=A*gg*A**-1
            #print AA[1,0]
            if AA not in self:
                if verbose>0:
                    print "A*{0}A^-1={1}".format(gg,AA)
                if verbose>1:
                    return False,AA
                else:
                    return False
        return True

    
    def coset_reps(self,version=0):
        r""" Returns coset reps of self


        EXAMPLES::

        
            sage: S=SymmetricGroup(6)
            sage: pS=S([2,1,4,3,6,5])
            sage: pR=S([3,1,2,5,6,4])
            sage: pS
            (1,2)(3,4)(5,6)
            sage: pR
            (1,3,2)(4,5,6)
            sage: G=MySubgroup(o2=pS,o3=pR)
            sage: G.coset_reps()
            [[1 0]
            [0 1], [1 2]
            [0 1], [1 1]
            [0 1], [1 3]
            [0 1], [1 5]
            [0 1], [1 4]
            [0 1], [ 4 -1]
            [ 1  0]]

        """
        if self._verbose>0:
            print "in coset_reps: farey = ",self._reps_from_farey
        if self._reps_from_farey:
            version = 2
        if version==0:
            if self._coset_reps_v0==None:
                self._coset_reps_v0 = self._get_coset_reps_from_perms() #self.permS,self.permR)
            return self._coset_reps_v0
        elif version==1:
            if self._coset_reps_v1==None:
                self._coset_reps_v1 = []
                for A in super(MySubgroup_class,self).coset_reps():
                    self._coset_reps_v1.append(SL2Z_elt(A.a(),A.b(),A.c(),A.d()))
   
            return self._coset_reps_v1
        elif version==2:
            if self._coset_reps_v2==None:
                self._coset_reps_v2 = []
                for A in self.farey_symbol().coset_reps():
                    self._coset_reps_v2.append(SL2Z_elt(A.d(),-A.b(),-A.c(),A.a()))
            return self._coset_reps_v2
        else:
            raise NotImplementedError

    
    def _get_perms_from_coset_reps(self):
        r""" Get permutations of order 2 and 3 from the coset representatives

        EXAMPLES::
        
            sage: S=SymmetricGroup(6)
            sage: pS=S([2,1,4,3,6,5])
            sage: pR=S([3,1,2,5,6,4])
            sage: G=MySubgroup(o2=pS,o3=pR)
            sage: G._get_perms_from_coset_reps()
            [(1,2)(3,4)(5,6), (1,3,2)(4,5,6)]
            sage: G=MySubgroup(Gamma0(6))
            sage: p=G._get_perms_from_coset_reps(); p[0]; p[1]
            (1,2)(3,4)(5,8)(6,9)(7,10)(11,12)
            (1,3,2)(4,5,9)(6,11,10)(7,12,8)


        """
        l=self.coset_reps()
        li=list()
        n=len(l)
        ps=range(1 ,self._index+1 )
        pr=range(1 ,self._index+1 )
        S=SL2Z_elt(0,-1,1,0); T=SL2Z_elt(1,1,0,1)
        R=SL2Z_elt(0,-1,1,1)
        #S,T=SL2Z.gens()
        #R=S*T
        level=self.level()
        if isinstance(l[0],SL2Z_elt):
            for i in range(n):
                li.append(l[i].inverse())
                #li.append([l[i][3],-l[i][1],-l[i][2],l[i][0]])
                #li.append( SL2Z(l[i])**-1 )
            ixr=range(n)
        else:
            for i in range(n):
                li.append( l[i].inverse() )
            ixr=range(n)
        ixs=range(n)
        for i in range(n):
            [a,b,c,d]=l[i]
            VS=SL2Z_elt(b,-a,d,-c) # Vi*S
            VR=SL2Z_elt(b,b-a,d,d-c) # Vi*R=Vi*S*T
            for j in ixr:
                Vji=li[j] # Vj^-1
                tmp = VR*Vji
                #tmp = mul_list_maps(VR,Vji)
                if self.is_Gamma0():
                    t = tmp[2] % level == 0 
                else:
                    t = list(tmp) in G
                if t:
                    pr[i]=j+1 
                    ixr.remove(j)
                    break
            for j in ixs:
                Vji=li[j]
                tmp = VS*Vji
                #tmp = mul_list_maps(VS,Vji)
                if self.is_Gamma0():
                    t = tmp[2] % level == 0 
                else:
                    t = list(tmp) in G
                if t: #(VS*Vji in G):
                    ps[i]=j+1 
                    ixs.remove(j)
                    break
        #print "ps=",ps
        #print "pr=",pr
        #print "PS=",MyPermutation(ps)
        #print "PR=",MyPermutation(pr)        
        #return [self._S(ps),self._S(pr)]
        return MyPermutation(ps),MyPermutation(pr)
    


    # Now to public routines

    def coset_rep(self,A):
        r"""
        Indata: A in PSL(2,Z) 
        Returns the coset representative of A in
        PSL(2,Z)/self.G

        EXAMPLES::


            sage: G=MySubgroup(Gamma0(4))        
            sage: A=SL2Z([9,4,-16,-7])
            sage: G.coset_rep(A)
            [1 0]
            [0 1]
            sage: A=SL2Z([3,11,-26,-95])
            sage: G.coset_rep(A)
            [-1  0]
            [ 2 -1]
            sage: A=SL2Z([71,73,35,36])
            sage: G.coset_rep(A)
            [ 0 -1]
            [ 1  0]

        
        """
        if isinstance(A,SL2Z_elt):
            for V in (self.coset_reps()):
                if  A*V.inverse() in self:
                    return V
        else:
            for V in (self.coset_reps()):
                if  V*A**-1  in self:
                    return V
        raise ArithmeticError,"Did not find coset rep. for A=%s" %(A)

    # def pullback(self,x_in,y_in,ret_mat=1,prec=201,**kwds):
    #     r""" Find the pullback of a point in H to the fundamental domain of self
    #     INPUT:

    #      - ''x_in,y_in'' -- x_in+I*y_in is in the upper half-plane
    #      - ''prec''      -- (default 201) precision in bits
    #      - ret_mat  -- set to 0 if you want to return a list instead of a matrix.
    #     OUTPUT:
        
    #      - [xpb,ypb,B]  --  xpb+I*ypb=B(x_in+I*y_in) with B in self
    #                        xpb and ypb are complex numbers with precision prec 
    #     EXAMPLES::


    #         sage: S=SymmetricGroup(6)
    #         sage: pS=S([2,1,4,3,6,5])
    #         sage: pR=S([3,1,2,5,6,4])
    #         sage: G=MySubgroup(o2=pS,o3=pR)
    #         sage: [x,y,B]=G.pullback(0.2,0.5,53); x,y;B
    #         (-0.237623762376238, 0.123762376237624)
    #         [-1  0]
    #         [ 4 -1]
    #         sage: (B**-1).acton(x+I*y)
    #         0.200000000000000 + 0.500000000000000*I


    #     """
    #     #x=deepcopy(x_in); y=deepcopy(y_in)
    #     version = kwds.get('version',0)
    #     if self.is_Gamma0() and version == 0:
    #         if isinstance(x_in,float):
    #             xpb,ypb,a,b,c,d=pullback_to_Gamma0N_dp(self,x_in,y_in,self._verbose)
    #         elif isinstance(x_in,Expression):
    #             prec=round(RR(len(str(x_in).split(".")[1])/log_b(2,10)))
    #             RF=RealField(prec)
    #             x=RF(x_in); y=RF(y_in)
    #             if prec<=53:
    #                 xpb,ypb,a,b,c,d=pullback_to_Gamma0N_dp(self,x,y,self._verbose)
    #             else:
    #                 xpb,ypb,a,b,c,d=pullback_to_Gamma0N_mpfr(self,x,y)
    #         else:
    #             xpb,ypb,a,b,c,d=pullback_to_Gamma0N_mpfr(self,x_in,y_in)
    #         if ret_mat==1:
    #             return xpb,ypb,SL2Z_elt(a,b,c,d)
    #         else:
    #             return xpb,ypb,int(a),int(b),int(c),int(d)
    #     else:
    #         A=pullback_to_psl2z_mat(RR(x_in),RR(y_in))
    #         A=SL2Z_elt(A) #.matrix()
    #         reps = self.coset_reps(version)
    #         try:
    #             for V in reps:
    #                 B=V*A
    #                 if B in self:
    #                     raise StopIteration
    #         except StopIteration:            
    #             pass
    #         else:
    #             raise ArithmeticError,"Did not find coset rep. for A=%s" % A
    #         #if ret_int==1:
    #         #    a,b,c,d=B[0,0],B[0,1],B[1,0],B[1,1]
    #         if isinstance(x_in,float):
    #             xpb,ypb=apply_sl2z_map_dp(x_in,y_in,B[0,0],B[0,1],B[1,0],B[1,1])
    #         else:
    #             xpb,ypb=apply_sl2z_map_mpfr(x_in,y_in,B[0,0],B[0,1],B[1,0],B[1,1])


    #         if ret_mat==1:
    #             return xpb,ypb,B.matrix()
    #         else:
    #             return xpb,ypb,B[0,0],B[0,1],B[1,0],B[1,1]
    def is_congruence(self):
        r""" Is self a congruence subgroup or not?

        EXAMPLES::


            sage: S=SymmetricGroup(6)
            sage: pS=S([2,1,4,3,6,5])
            sage: pR=S([3,1,2,5,6,4])
            sage: G=MySubgroup(o2=pS,o3=pR)
            sage: G.is_congruence()
            True
            sage: S=SymmetricGroup(7)
            sage: pS=S([1,3,2,5,4,7,6]); pS
            (2,3)(4,5)(6,7)
            sage: pR=S([3,2,4,1,6,7,5]); pR
            (1,3,4)(5,6,7)
            sage: G=MySubgroup(o2=pS,o3=pR)
            sage: G.is_congruence()
        False
        
        """
        if self._is_congruence==None:
            self._is_congruence=super(MySubgroup_class,self).is_congruence()
        return self._is_congruence

        
    def generalised_level(self):
        r""" Generalized level of self

        EXAMPLES::y

            sage: S=SymmetricGroup(6)
            sage: pS=S([2,1,4,3,6,5])
            sage: pR=S([3,1,2,5,6,4])
            sage: G=MySubgroup(o2=pS,o3=pR)
            sage: G.generalised_level()
            4
            sage: S=SymmetricGroup(7)
            sage: pS=S([1,3,2,5,4,7,6]); pS
            (2,3)(4,5)(6,7)
            sage: pR=S([3,2,4,1,6,7,5]); pR
            (1,3,4)(5,6,7)
            sage: G=MySubgroup(o2=pS,o3=pR)
            sage: G.generalised_level()
            6
        """
        if self._generalised_level==None:
            # compute the generalized level
            self._generalised_level= lcm(map(len,self.permT.cycle_tuples()))
        return self._generalised_level
    #raise ArithmeticError, "Could not compute generalised level of %s" %(self)



    # def level(self):
    #     r""" Level of self

    #     EXAMPLES::


    #         sage: G=MySubgroup(Gamma0(5));
    #         sage: G.level()
    #         5

    #     """
    #     if(self._is_congruence):
    #         return self.generalised_level()
    #     else:
    #         raise TypeError,"Group is not a congruence group! Use G.generalised_level() instead!"

    
    def generators_as_slz_elts(self):
        if self._generators_as_slz_elts==None or self._generators_as_slz_elts==[]:
            self._generators_as_slz_elts=[]
            for A in self.gens():
                a,b,c,d=A
                self._generators_as_slz_elts.append(SL2Z_elt(a,b,c,d))
        return self._generators_as_slz_elts








    def closest_vertex(self,x,y,as_integers=1):
        r"""
        The closest vertex to the point z=x+iy in the following sense:
        Let sigma_j be the normalized cusp normalizer of the vertex p_j, 
        i.e. sigma_j^-1(p_j)=Infinity and sigma_j*S_j*sigma_j^-1=T, where
        S_j is the generator of the stabiliser of p_j

        The closest vertex is then the one for which Im(sigma_j^-1(p_j))
        is maximal.

        INPUT:

         - ''x,y'' -- x+iy  in the upper half-plane

        OUTPUT:
        
         - ''v'' -- the closest vertex to x+iy
         
        
        EXAMPLES::


        sage: G=MySubgroup(Gamma0(5))
        sage: G.closest_vertex(-0.4,0.2)
        Infinity
        sage: G.closest_vertex(-0.1,0.1)
        0

        """
        ci=closest_vertex(self._vertex_maps,self._vertex_widths,self._nvertices,x,y,self._verbose)
        if as_integers:
            return ci
        else:
            return self._vertices[ci]

    def closest_cusp(self,x,y,vertex=0,as_integers=1):
        r"""
        The closest cusp to the point z=x+iy in the following sense:
        Let sigma_j be the normalized cusp normalizer of the vertex p_j, 
        i.e. sigma_j^-1(p_j)=Infinity and sigma_j*S_j*sigma_j^-1=T, where
        S_j is the generator of the stabiliser of p_j

        The closest vertex is then the one for which Im(sigma_j^-1(p_j))
        is maximal and the closest cusp is the cusp associated to this vertex

        INPUT:

         - ''x,y'' -- x+iy  in the upper half-plane

        OUTPUT:
        
         - ''v'' -- the closest vertex to x+iy
         
        
        EXAMPLES::


        sage: G=MySubgroup(Gamma0(5))
        sage: G.closest_vertex(-0.4,0.2)
        Infinity
        sage: G.closest_vertex(-0.1,0.1)
        0

        """
        vi = closest_vertex(self._vertex_maps,self._vertex_widths,self._nvertices,x,y,self._verbose)
        ci = self._vertex_data[vi]['cusp']
        if vertex==1:
            if as_integers:
                return ci,vi
            else:
                return self._cusps[ci],self._vertices[vi]
        else:
            if as_integers:
                return ci
            else:
                return self._cusps[ci]



    
    
##   def block_systems(self):
##         r"""
##         Return a list of possible Block systems and overgroups for the group G.
##         INDATA: G = subgroup of the modular group
##         OUTDATA res['blocksystems']= Block systems
##         res['overgroups']  = Overgroups corresponding to the blocks
##         in the block systems


##         EXAMPLES::
        

##         """
##         return get_block_systems(self.permS,self.permR,False)

    def dimension_cuspforms(self,k):
        r"""
        The implementation of this in the standard arithmetic subgroup package
        returns a false answer. Therefore it is better to overwrite it...
        """
        return self.dimension_cuspforms(k)

    def dimension_cuspforms(self,k):
        r"""
        Returns the dimension of the space of cuspforms on G of weight k
        where k is an even integer

        EXAMPLES::


            sage: G=MySubgroup(Gamma0(4))
            sage: G.dimension_cuspforms(12)
            4
            sage: S=SymmetricGroup(7)
            sage: pR=S([3,2,4,1,6,7,5]); pR
            (1,3,4)(5,6,7)
            sage: pS=S([1,3,2,5,4,7,6]); pS
            (2,3)(4,5)(6,7)
            sage: G.dimension_cuspforms(4)
            1

        
        """
        #kk=ZZ(k)
        ki = ZZ(k)
        if is_odd(ki):
            raise ValueError, "Use only for even weight k! not k=" %(kk)
        if ki<2:
            dim=0 
        elif ki==2:
            dim=self._genus
        elif ki>=4:
            kk = RR(k)
            dim=ZZ(kk-1)*(self._genus-1)+self._nu2*floor(kk/4.0)+self._nu3*floor(kk/3.0)
            dim+= ZZ(kk/2.0 - 1)*self._ncusps
        return dim

    def dimension_modularforms(self,k):
        r"""
        Returns the dimension of the space of modular forms on G of weight k
        where k is an even integer

        EXAMPLES::

            sage: S=SymmetricGroup(7)
            sage: pR=S([3,2,4,1,6,7,5]); pR
            (1,3,4)(5,6,7)
            sage: pS=S([1,3,2,5,4,7,6]); pS
            (2,3)(4,5)(6,7)
            sage: G.dimension_modularforms(4)
            3
        """
        kk=Integer(k)
        if is_odd(kk):
            raise ValueError, "Use only for even weight k! not k=" %(kk)
        if k==0 :
            dim=1  # the constant functionz 
        elif k<2:
            dim=0 
        else:
            dim=self.dimension_cuspforms(k)+self._ncusps
            #(kk-1.0)*(self._genus-_sage_const_1 )+self._nu2()*int(floor(kk/_sage_const_4 ))+self._nu3*int(floor(kk/_sage_const_3 ))+kk/_sage_const_2 *self._ncusps()
        return dim
    
    ### Overloaded operators
    def __contains__(self,A):
        r"""
        Is A an element of self (this is an ineffective implementation if self is a permutation group)

        EXAMPLES::
        

            sage: G=MySubgroup(Gamma0(5))
            sage: A=SL2Z([-69,-25,-80,-29])
            sage: G.__contains__(A)
            True
            sage: A in G
            True            
            sage: S=SymmetricGroup(7)
            sage: pR=S([3,2,4,1,6,7,5]); pR
            (1,3,4)(5,6,7)
            sage: pS=S([1,3,2,5,4,7,6]); pS
            (2,3)(4,5)(6,7)
            sage: S,T=SL2Z.gens(); R=S*T
            sage: A=S*T^4*S
            sage: A in G
            False
            sage: A=S*T^6*S
            sage: A in G
            True
            
        """
        if hasattr(A,"det"):
            if A.det()<>1:
                if self._verbose>3:
                    print "det(A)=",A.det()
                return False
            if not is_integer(A[0,0]):
                if self._verbose>3:
                    print "A[0,0]=",A[0,0],type(A[0,0])
                return False
            if not is_integer(A[0,1]):
                if self._verbose>3:
                    print "A[0,1]=",A[0,1],type(A[0,1])

                return False
            if not is_integer(A[1,0]):
                if self._verbose>3:
                    print "A[1,0]=",A[1,0],type(A[1,0])
                return False
            if not is_integer(A[1,1]):
                if self._verbose>3:
                    print "A[1,1]=",A[1,1],type(A[1,1])
                return False
        elif isinstance(A,list):
            if len(A)>4:
                return False
            if A[0]*A[3]-A[1]*A[2]<>1:
                return False
            # # See if A has integer entries
            if not is_integer(A[0]):
                return False
            if not is_integer(A[1]):
                return False
            if not is_integer(A[2]):
                return False
            if not is_integer(A[3]):
                return False
        elif not isinstance(A,SL2Z_elt):
            return False
        p=self.permutation_action(A)
        if self._verbose>1:
            print "perm(A)=",p
        if p(1)==1:
            return True
        else:
            return False

    def cusps(self):
        r"""
        Returns the cusps of self as cusps

        EXAMPLES::


            sage: G=MySubgroup(Gamma0(5))
            sage: G.cusps()
            [Infinity, 0]
        
        """
        if len(self._cusps_as_cusps)==0:
            for c in self._cusps:
                if c[1]<>0:
                    self._cusps_as_cusps.append(Cusp(QQ(c[0])/QQ(c[1])))
                else:
                    self._cusps_as_cusps.append(Cusp(1,0))
        return self._cusps_as_cusps

    def vertices(self):
        r"""
        Return the certices of the (current) fundamental domain of self.
        """
        if len(self._vertices_as_cusps)==0:
            for c in self._vertices:
                if c[1]<>0:
                    self._vertices_as_cusps.append(Cusp(QQ(c[0])/QQ(c[1])))
                else:
                    self._vertices_as_cusps.append(Cusp(1,0))
        return self._vertices_as_cusps

    def nvertices(self):
        return self._vertices

    def cusp_width(self,c):
        r"""
        Returns the cusp width of cusp

        EXAMPLES::


            sage: G=MySubgroup(Gamma0(4))
            sage: G.cusp_data(Cusp(1/2))
            ([-1  1]
            [-4  3], 1, 1)
            sage: G.cusp_width(Cusp(-1/2))
            1
            sage: G.cusp_width(Cusp(1/2))
            1
            sage: G.cusp_width(Cusp(0))
            4
        """
        if isinstance(c,tuple):
            p=c[0];q=c[1]; cusp=None
        else:
            cusp = Cusp(c)
            p=cusp.numerator(); q=cusp.denominator()
        if (p,q) in self._cusps:
            j = self._cusps.index((p,q))
            return self._cusp_data[j]['width']
        else:
            if cusp == None:
                cusp = Cusp(p,q)
            if self._verbose>1:
                print "cusp=",cusp
            # if we are here we did not find the cusp in the list so we have to find an equivalent in the list
            c=self.cusp_representative(cusp)
            p = c.numerator(); q = c.denominator()
            j = self._cusps.index((p,q))
            return self._cusp_data[j]['width']
        
        raise ArithmeticError,"Could not find the width of %s" %cusp

    def cusp_data(self,c):
        r""":
        Returns cuspdata in the same format as for the generic Arithmetic subgroup, i.e. a tuple (A,h,s) where A is a generator of the stabiliser of c, h is the width of c and s is the orientation. 

        INPUT:

        - 'c' -- Integer or cusp
        
        
        EXAMPLES::


            sage: G=MySubgroup(Gamma0(4))
            sage: G.cusp_data(Cusp(1/2))
            ([-1  1]
            [-4  3], 1, 1)
            sage: G.cusp_data(Cusp(-1/2))
            ([ 3  1]
            [-4 -1], 1, 1)
            sage: G.cusp_data(Cusp(0))
            ([ 0 -1]
            [ 1  0], 4, 1)

        """
        if self._cusp_data_sage_format.get(c)==None:        
            cusp = Cusp(c)
            p = cusp.numerator(); q = cusp.denominator()
            if (p,q) not in self._cusps:
                 ## Then we compute everything using the same method as in sage but with SL2Z_elt
                ## to make it faster
                w = lift_to_sl2z(c.denominator(), c.numerator(), 0)
                g = SL2Z_elt([w[3], w[1], w[2],w[0]])
                for d in range(1,1+self.index()):
                    t = g * SL2Z_elt(1,d,0,1) * g.inverse()
                    print t
                    if t in self:
                        self._cusp_data_sage_format[c] = t, d, 1
                        break
                    else:
                        t = SL2Z_elt(-t.a(),-t.b(),-t.c(),-t.d())
                        if t in self: # Note that in the current implementation this will never hold since we work exclusively with PSL
                            self._cusp_data_sage_format[c] = t, d, -1
                            break
            else:
                i = self._cusps.index((p,q))
                t = self._cusp_data[i]['stabilizer']
                d = self._cusp_data[i]['width']
                self._cusp_data_sage_format[c] = t, d, 1
            if self._cusp_data_sage_format.get(c)==None:
                raise ArithmeticError,"Could not find stabiliser of {0}".format(c)
        return self._cusp_data_sage_format.get(c)


        ## try:
        ##     # if we give a cusp already in the list we return stored values
        ##     w=self.cusp_normalizer(c)
        ##     d=self._cusp_width[c][0 ]
        ##     e=self._cusp_width[c][1 ]
        ##     g=w * SL2Z([e,d*e,0 ,e]) * (~w)
        ##     return (self.cusp_normalizer(c),d*e,1)
        ## except:
        ##     w = lift_to_sl2z(c.denominator(), c.numerator(), 0 )
        ##     g = SL2Z([w[3 ], w[1 ], w[2 ],w[0 ]])
        ##     for d in xrange(1 ,1 +self.index()):
        ##         if g * SL2Z([1 ,d,0 ,1 ]) * (~g) in self:
        ##             return (g * SL2Z([1 ,d,0 ,1 ]) * (~g), d, 1 )
        ##         elif g * SL2Z([-1 ,-d,0 ,-1 ]) * (~g) in self:
        ##             return (g * SL2Z([-1 ,-d,0 ,-1 ]) * (~g), d, -1 )
        ##     raise ArithmeticError, "Can' t get here!"
        ## #

    def reduce_cusp(self,c,trans=False):
        r"""
        Find a reduced cusp which is equivalent to the given cusp.
        For a congruence subgroup we return the same reduced cusp as usual.
        This is not defined for a non-congruence group.
        
        """
        return super(MySubgroup_class,self).reduce_cusp(c)


    def cusp_representative(self,cusp,trans=False):
        r"""
        Return a cusp in self.cusps() which is equivalent to the given cusp.

        OUTPUT:

         # If trans ==  True
         - 'c,A' -- tuple with c a cusp and A an element of SL2Z_elt
         # If trans =  False
         - 'c' -- a cusp 
        
        """
        cc = Cusp(cusp)
        for x in self.cusps():
            t = self.are_equivalent(x,cc)
            if not trans and t:
                return x
            if trans and t:
                return x,t
        raise ArithmeticError,"Cusp {0} is not equivalent to any representative!".format(cusp)
    
    # def cusp_equivalent_to(self,cusp):
    #     r"""
    #     Find a cusp in self._cusps which is equivalent to cusp
    #     and returns the cusp, a map which maps the given point to the cusp 

        
    #     """
    #     p=None;q=None
    #     if isinstance(cusp,tuple):
    #         p=cusp[0];q=cusp[1]
    #     else:
    #         try:
    #             cusp=Cusp(cusp)
    #             p=cusp.numerator(); q=cusp.denominator()
    #         except:
    #             raise TypeError,"Could not coerce {0} to a cusp!".format(cusp)
                
    #     if (p,q) in self._cusps:
    #         one = SL2Z_elt(int(1),int(0),int(0),int(1))
    #         return (p,q),one,one
    #     #print "p,q=",p,q,type(p),type(q)
    #     w = lift_to_sl2z(q, p, 0 )
    #     V = SL2Z_elt(w[3 ], w[1 ], w[2 ],w[0 ])
    #     permv=self.permutation_action(V)
    #     for i in range(self._ncusps):
    #         W = self._cusp_data[i]['normalizer']
    #         permw=self.permutation_action(W)
    #         testi=(permw**-1)(1)
    #         for k in range(0,self._index):
    #             test=(permv*self.permT**k)(testi)
    #             if test==1:  ## v = C w with C = V*T**k*W**-1
    #                 Tk = SL2Z_elt(1,k,0,1)
    #                 mapping=V*Tk*W.inverse()
    #                 return (p,q),mapping,W
    #     raise ArithmeticError, "Could not find equivalent cusp!"

    def cusp_normalizer(self,cusp):
        r"""
        Return the cusnormalizer of cusp

        EXAMPLES::


            sage: S=SymmetricGroup(6)
            sage: pS=S([2,1,4,3,6,5])
            sage: pR=S([3,1,2,5,6,4])
            sage: G=MySubgroup(o2=pS,o3=pR)
            sage: G.cusps()
            [Infinity, 0, -1/2]
            sage: G.cusp_normalizer(Cusp(0))
            [ 0 -1]
            [ 1  0]
            sage: G.cusp_normalizer(Cusp(-1/2))
            [-1  0]
            [ 2 -1]
            sage: G.cusp_normalizer(Cusp(1,0))
            [1 0]
            [0 1]

        
        """
        if isinstance(cusp,(int,Integer)):
            cusp = self._cusps[cusp]
        if isinstance(cusp,Cusp):        
            p=cusp.numerator(); q=cusp.denominator()
        elif isinstance(cusp,tuple):
            p=cusp[0];q=cusp[1]
        else:
            raise ValueError,"Could not construct a cusp from {0}".format(cusp)
        if (p,q) in self._cusps:
            j = self._cusps.index((p,q))
            return self._cusp_data[j]['normalizer']
        else:
            c,A=self.cusp_representative(cusp,trans='matrix')
            p = c.numerator(); q = c.denominator()
            if A==1:
                j = self._cusps.index((p,q))
                return self._cusp_data[j]['normalizer']
            else:
                w = lift_to_sl2z(q, p, 0 )
                g = SL2Z_elt(w[3 ], w[1 ], w[2 ],w[0 ])
            return g

    def normalize_to_cusp(self,x,y,c,inv=0):
        r"""
        Apply the cusp normalizing map of cusp c to x+iy
        """
        xx=x.parent()(x); yy=y.parent()(y)
        [d,b,c,a]=self.cusp_normalizer(c)
        wi = G.cusp_widt(c)
        if inv==1:
            b=-b; c=-c; atmp=a;
            a=d; d=a
            if wi<>1:
                wmul=self.cusp_widt(c)**-1
        elif wi<>1:
            wmul=x.parent()(self.cusp_widt(c))
            xx=wmul*xx; yy=wmul*yy
        xx,yy=apply_sl2z_map(xx,yy)
        if wi<>1 and inv==1:
            xx=xx/wmul; yy=yy/wmul
        return xx,yy


    def minimal_height(self):
        r""" Computes the minimal (invariant) height of the fundamental region of self.

        Note: Only guaranteed for Gamma_0(N)
        EXAMPLES::

        
            sage: G=MySubgroup(Gamma0(6))
            sage: G.minimal_height()
            0.144337567297406
            sage: S=SymmetricGroup(6)
            sage: pS=S([2,1,4,3,6,5])
            sage: pR=S([3,1,2,5,6,4])
            sage: G=MySubgroup(o2=pS,o3=pR)
            sage: G.minimal_height()   
            0.216506350946110


        """
        if self.is_Gamma0():
            l=self.generalised_level()
            if self.is_Gamma0():
                return RR(sqrt(3.0))/RR(2*l)
        # For all other groups we have have to locate the largest width
        maxw=0
        for i in range(self._ncusps):
            l=self._cusp_data[i]['width']
            if l>maxw:
                maxw=l
        return RR(sqrt(3.0))/RR(2*maxw)

    @rename_keyword(color='rgbcolor')
    @options(alpha=1, fill=True, thickness=1, rgbcolor="black", \
             zorder=2, linestyle='solid', show_pairing=False, \
             show_tesselation=True,
             method='Farey',
             model='H')
    
    def draw_fundamental_domain(self,**options):
        r""" Draw fundamental domain

        OPTIONS:
        
         - ''model'' -- (default ''H'')
             = ''H'' -- Upper halfplane
             = ''D'' -- Disk model
         - ''filename''-- filename to print to
         - ''**kwds''-- additional arguments to matplotlib 
         - ''axes''  -- set geometry of output
             =[x0,x1,y0,y1] -- restrict figure to [x0,x1]x[y0,y1]

        EXAMPLES::

            sage: G=MySubgroup(Gamma0(3))
            sage: G.draw_fundamental_domain()

        """
        from plot_dom import HyperbolicTriangle
        draw_axes=options.pop('draw_axes',1)
        npts = options.pop('npts',10)
        if options['method']=='Farey' and options['model']=='H':
            options.pop('model'); options.pop('method');
            return self.farey_symbol().fundamental_domain(**options)
        if options['show_pairing'] and not options['method']=='Farey':
            raise NotImplementedError,"Pairings are only implemented for Farey symbols."
        draw_circle = options.pop('draw_circle',True)
        from sage.plot.colors import rainbow
        L = 10000
        #if options['method']=='Farey':
        #    version = 2
        #else:
        model = options['model']
        verbose = options.get('verbose',0)
        ret_domain = options.pop('domain',False)
        contour_only= options.pop('contour',False)
        version = options.pop('version',0)
        circle_color = options.pop('circle_color','black')
        conjugate_A = options.pop('conjugate_by',SL2Z_elt(1,0,0,1)) ## We draw a conjugated fundamental domain
        T = SL2Z_elt(1,1,0,1)
        S = SL2Z_elt(0,-1,1,0)
        Z = SL2Z_elt(-1,0,0,-1)        
        if conjugate_A == None:
            coset_reps = self.coset_reps(version=version)
        else:
            if verbose>0:
                print "Conjugating by A=",conjugate_A
            coset_reps = []
            for x in self.coset_reps(version=version):
                coset_reps.append(conjugate_A*x)
        if verbose>0:
            print "options=",options
        if model=="D2":
            g=draw_funddom_d(coset_reps,format,I)
        else:
            g = Graphics()
            cntr = Graphics()
            A0 = CC(-0.5,sqrt(3.)/2)
            B0 = CC(0.5,sqrt(3.)/2)
            if model == 'D':
                C0 = CC(infinity)
            else:
                C0 = CC(0,L)
            for x in coset_reps:
                #a, b, c, d = x[3], -x[1], -x[2], x[0]
                A,B = [x.acton(z) for z in [A0,B0]]
                if x.c()<>0:
                    C = CC(RR(x.a())/RR(x.c()),0)
                else:
                    if model == 'D':
                        C = CC(infinity)
                    else:
                        C = CC( (A.real()+B.real())*0.5,L)
                if verbose>0:                    
                    print "coset rep=",x
                    print "Triangle: ({0},{1},{2})".format(A,B,C)
                sides = [1,2,3]                    
                if contour_only: # and model=='H':
                    ## See which sides we keep...
                    # if x.c()==0: ## Have vertical side
                    #     sides = []                    
                    #     if x*T**-1 not in coset_reps:  
                    #         sides.append(3)
                    #     if x*T not in coset_reps:
                    #         sides.append(2)
                    #     if x*S not in coset_reps:
                    #         sides.append(1)
                    sides = []                    
                    if x.c()==0:
                        Nj = SL2Z_elt(1,0,0,1)
                    else:
                        p = Cusp(x.a()/x.c())
                        Nj = self.cusp_normalizer(p)
                    V = Nj.inverse()*x
                    if verbose>0:
                        print "x=",x
                        print "V=",V
                    if Nj*V*T**-1 not in coset_reps:  
                        sides.append(3)
                    if Nj*V*T not in coset_reps:
                        sides.append(2)
                    if Nj*V*S not in coset_reps and Z*Nj*V*S not in coset_reps:
                        if verbose>1:
                            print "NjVS=",Nj*V*S
                        sides.append(1)
                    if verbose>0:
                            print "sides = ",sides

                if sides<>[]:
                    t = my_hyperbolic_triangle(A, B, C, \
                                               color=options['rgbcolor'], \
                                               fill=False, \
                                               alpha=options['alpha'], \
                                               thickness=options['thickness'], \
                                               model=model, verbose=verbose,npts=npts,sides=sides)
                    cntr += t
                if model=='H' and options['show_tesselation'] and options.get('fill',True)==True:
                    g += my_hyperbolic_triangle(A, B, C, color="lightgray",fill=True,
                                                model=model)
            g+=cntr
        if False and contour_only==True:
            # Remove interior arcs from path...
            #print "Removing interior arcs:"
            if model == 'D':
                n = len(g)
                for j in range(n-1,-1,-1):
                    if not hasattr(g[j],"path"):
                        continue
                    path0 = g[j].path
                    path1 = []
                    path_has_interior_pts = False
                    if verbose>0:
                        print "Checking path=",path0
                    for p in path0:
                        p0 = copy(p)
                        for i in range(len(p)):
                            x = CC(p[i])
                            #print "x(orig)=",x,abs(x)==1.0
                            if model == 'D':
                                if abs(abs(x)-1.0)<1e-12:
                                    continue
                                x = CC(0,1)*(CC(x)+CC(1,0))/(CC(1,0)-CC(x))
                            #print "x=",x
                            if self.is_interior_point(x,version=version,verbose=verbose-1):
                                path_has_interior_pts = True
                                if verbose>0:
                                    print "point {0} is interior!".format(x)
                                break
                                #del(p0[i])
                        #path1.append(p0)
                    if path_has_interior_pts:
                        #new_prim = deepcopy(g[j])
                        #new_prim.path = path1
                        #print "removing path:",g[j]
                        del(g[j])
                        #print "adding path:",new_prim
                        #g.add_primitive(new_prim)
                
        d = g.get_minmax_data()
        if model=='H':
            g.set_axes_range(d['xmin'], d['xmax'], 0, min(d['ymax'],2))
            g.SHOW_OPTIONS['ticks']=[range(int(d['xmin']),int(d['xmax'])+1),[1,2]]
        else:
            if not ret_domain and draw_circle:
                g+=circle((0,0),1,edgecolor=circle_color)
            g.set_axes_range(-1, 1, -1, 1)    
            g.SHOW_OPTIONS['ticks']=[range(int(d['xmin']),int(d['xmax'])+1),[1,2]]        
        if draw_axes == 0:
            g.axes(False)
        return g
    def is_interior_point(self,x,version=1,verbose=0):
        r"""
        Test if x is an interior point of the desired fundamental domain of self.
        """
        ep = 1e-2
        x0 = CC(x) + CC(0,ep)
        x1 = CC(x) + CC(0,-ep)
        x2 = CC(x) + CC(ep,0)
        x3 = CC(x) + CC(-ep,0)
        z0 = self.pullback(x0.real(),x0.imag(),version=version)
        if isinstance(z0,tuple):
            z0 = CC(z0[0],z0[1])
        if verbose>0:
            print "Pb of {0} is {1}; diff={2}".format(x0,z0,abs(x0-z0))
        if z0 <> x0:
            return False           
        z1 = self.pullback(x1.real(),x1.imag(),version=version)
        if isinstance(z1,tuple):
            z1 = CC(z1[0],z1[1])
        if verbose>0:
            print "Pb of {0} is {1}; diff={2}".format(x1,z1,abs(x1-z1))
        if z1 <> x1:
            return False
        try:
            z2 = self.pullback(x2.real(),x2.imag(),version=version)
        except OverflowError as e:
            print "Could not pull back {0}".format(x2)
            return False
        if isinstance(z2,tuple):
            z2 = CC(z2[0],z2[1])
        if verbose>0:
            print "Pb of {0} is {1}; diff={2}".format(x2,z2,abs(x2-z2))
        if z2 <> x2:
            return False
        z3 = self.pullback(x3.real(),x3.imag(),version=version)
        if isinstance(z3,tuple):
            z3 = CC(z3[0],z3[1])
        if verbose>0:
            print "Pb of {0} is {1}; diff={2}".format(x3,z3,abs(x3-z3))
        if z3 <> x3:
            return False
        return True
        
    def show_symmetry_props(self):
        r"""
        Display the symmetry properties of self.
        """
        for i in range(self._ncusps):
            print "Cusp Nr. {0} = {1}:{2}".format(i,self._cusps[i][0],self._cusps[i][1])
            print "Normalizer is normalizer : {0}".format(self.cusp_normalizer_is_normalizer(i))
            print "Can be symmetrized even/odd: {0}".format(self.is_symmetrizable_even_odd(i))




    def test_normalize_pt_to_cusp(self,z,ci):
        N=self._cusp_data[ci]['normalizer']
        l=self._cusp_data[ci]['width']
        print "N,l=",N,l
        z1=N.acton(z*RR(l))
        print "z1=",z1
        x,y,A= self.pullback(RR(z1.real()),RR(z1.imag()))
        z2=x+I*y
        print "z2=",z2
        cj,vj=self.closest_cusp(x,y,1)
        print "cj,vj=",cj,vj
        U = self._cusp_maps[vj]
        z3=U.acton(z2)
        print "z3=",z3
        N2=self._cusp_data[cj]['normalizer']
        l2=self._cusp_data[cj]['width']
        z4 = (N2**-1).acton(z3)/RR(l2)
        return z4


    def _get_all_cusp_data(self,coset_reps,test=0):
        r""" Return lists of vertices, inequivalent cusps, cusp normalizers, stabilisers and widths
             from a list of right coset representatives.
             Note: It is much more efficient to have it all here than using generic subgroutines.

        OUTPUT:
          -- [ns,ss,ws]  - list of dictionaries with cusps p as keys.
             ns : ns[p] = A in SL2Z with A(p)=infinity (cusp normalizer)
             ss : ss[p] = A in SL2Z s.t. A(p)=p   (cusp stabilizer)
             ws : ws[p]= (w,s).
                 w = width of p and s = 1 if p is regular and -1 if not

        EXAMPLES::
        

            sage: S=SymmetricGroup(6)
            sage: pS=S([2,1,4,3,6,5])
            sage: pR=S([3,1,2,5,6,4])
            sage: G=MySubgroup(o2=pS,o3=pR)     
            sage: l=G._get_cusp_data()
            sage: l[0]
            {0: [ 0 -1]
            [ 1  0], Infinity: [1 1]
            [0 1], -1/2: [-1  0]
            [ 2 -1]}
            sage: l[1]
            {0: [ 1  0]
            [-4  1], Infinity: [1 1]
            [0 1], -1/2: [ 3  1]
            [-4 -1]}
            sage: l[2]
            {0: (4, 1), Infinity: (1, 1), -1/2: (1, 1)}

        """
        vertices=list()
        cusps=list()
        widths=dict()
        normalizers=dict()
        stabilizer=dict()
        vertex_data=dict()
        permT=self.permT
        if self._verbose>0:
            print "permT=",self.permT
        lws=self.permT.cycles() #.cycle_tuples()
        if self._verbose>0:
            print "permT.cycles=",lws
        ## This is a set of all cusp widths which can occur
        all_cusp_widths=map(len,lws)
        if self._verbose>0:
            print "all cusp widths=",all_cusp_widths
        cusp_widths={} #
        ## We get the width from the cycles:
        for j in range(len(coset_reps)):
            for c in lws:
                if j+1 in c:
                    cusp_widths[j]=len(c)
                    break
        gen_level=lcm(all_cusp_widths) #cusp_widths.values())
        coset_perms=[MyPermutation(length=self._index)]
        vertex_data=dict()
        cusp_data=dict()
        for j in range(1,len(coset_reps)):
            V = coset_reps[j]
            if self._verbose>0:
                print "V[",j,"]=",V
            perm=self.permutation_action(V,1)
            if self._verbose>0:
                print "perm(",j,")(1)=",perm(1)
            coset_perms.append(perm)
        # Recall that we have permutations for all groups here
        # and that the widths can be read from self.permT
        if self._verbose>0:
             print "widths=",cusp_widths
             print "reps=",coset_reps
        Id=SL2Z_elt(1,0,0,1)
        vi=0
        ## First populate the vertices of the fundamental domain we have choosen
        for j in range(len(coset_reps)):
            if coset_reps[j][2]==0:
                v=1,0
                if v not in vertices:
                    vertices.append(v)
                    vertex_data[vi]={'cusp':0,'cusp_map':Id,'coset':[j],'width':int(cusp_widths[j])}
                    vi+=1
                else:
                    vj = vertices.index(v)
                    vertex_data[vj]['coset'].append(j)
            else:
                v=(QQ(coset_reps[j][0])/QQ(coset_reps[j][2]))
                v0 = v.numerator(); v1=v.denominator()
                v=(v0,v1)
                if v not in vertices:
                    vertices.append(v)
                    vertex_data[vi]={'cusp':-1,'cusp_map':Id,'coset':[j],'width':int(cusp_widths[j])}
                    vi+=1
                else:
                    vj = vertices.index(v)
                    vertex_data[vj]['coset'].append(j)
        if self._verbose>0:
            print "vi=",vi
            print "vertex_data=",vertex_data
        # Then the cusps
        if test==1:
            return
        ci=0
        norm_perm={}
        if (1,0) in vertices:
            cusps=[(1,0)]
            width = vertex_data[0]['width']
            vertex_data[0]['cusp']=0
            cusp_data[0]={'normalizer':Id,'width':int(width),'stabilizer':SL2Z_elt(1,width,0,1),'coset':[0]}
            ci+=1
        for j in range(len(vertices)):
            v = vertices[j]
            if self._verbose>0:
                print "Test:",v
            if v in cusps:
                cii=cusps.index(v)
                vertex_data[j]['cusp']=cii
                #cusps.index(v)
                if cusp_data[cii].has_key('vertices'):
                    cusp_data[cii]['vertices'].append(j)
                else:
                    cusp_data[cii]['vertices']=[j]
                continue
            # Check which "canonical cusp" v is equivalent to.
            W,U,p,q,l=self.get_equivalent_cusp(v[0],v[1])
            if self._verbose>0:
                if W<>0:
                    print "W,U,p,q,l=",list(W),list(U),p,q,l
                else:
                    print "W,U,p,q,l=",W,U,p,q,l
                    print "cusps=",cusps
            if (p,q) in cusps:
                cii=cusps.index((p,q))
                vertex_data[j]['cusp']=cii
                vertex_data[j]['cusp_map']=U
                cusp_data[cii]['vertices'].append(j)
                continue
                    
            elif W==0:  ## Check if this cusp is equivalent to another
                do_cont=0
                for pp,qq in cusps:
                    t,A=Cusp(pp,qq).is_gamma0_equiv(Cusp(v),1,transformation='matrix')
                    if self._verbose>0:
                        print "are eq:",Cusp(pp,qq),Cusp(v)
                        print "A=",A,type(A)
                        
                    if t==1 and A in self:
                        try: 
                            a,b,c,d=A
                        except ValueError as e:
                            a=A[0,0]; b=A[0,1]; c=A[1,0]; d=A[1,1]
                            if self._verbose:
                                print "We called with a sage matrix A!"
                        cii=cusps.index((pp,qq))
                        vertex_data[j]['cusp']=cii
                        vertex_data[j]['cusp_map']=SL2Z_elt(d,-b,-c,a)
                        cusp_data[cii]['vertices'].append(j)
                        do_cont=1
                        break
                if do_cont==1:
                    continue


            # We add this vertex or the canonical rep.
            cusp_data[ci]={}
            cusp_data[ci]['vertices']=[]
            vertex_data[j]['cusp']=ci
            if p==0 and q==0:
                if self._verbose>0:
                    print "1 setting cusp ",v
                cusps.append(v)
                vertex_data[j]['cusp_map']=Id
                cusp_data[ci]['vertices'].append(j)
            else:
                if self._verbose>0:
                    print "2 setting cusp ",p,q
                cusps.append((p,q))
                vertex_data[j]['cusp_map']=U
                cusp_data[ci]['vertices'].append(j)
                if (p,q) not in vertices:
                    for vj in vertex_data[j]['coset']:
                        Vtmp = coset_reps[vj]
                        coset_reps[vj] = U*Vtmp
                        if self._verbose>0:
                            print "Changing rep from {0} to {1}".format(Vtmp,coset_reps[j])
                    vertices[j]=(p,q)
                    vertex_data[j]['cusp_map']=Id
            # Setting the normalizer    
            if W<>0:
                if l==0:
                    l = vertex_data[j]['width']
                elif l<>vertex_data[j]['width']:
                    raise ArithmeticError,"Could not calculate width"
                if self._verbose>0:
                    print "3 setting cusp ",p,q
                    print "width=",l
                #S=W*SL2Z([1,l,0,1])*W**-1
                if Cusp(p,q)==Cusp(-1,2):
                    #cusp_data[ci]['normalizer']=
                    W = SL2Z_elt(1,0,-2,1)
                else:
                    W = U*W
                cusp_data[ci]['normalizer']=W
                S = W*SL2Z_elt(1,l,0,1) #S = mul_list_maps(W,[1,l,0,1])
                S = S._mul(W,2) #S = mul_list_maps(S,W,inv=2)
                cusp_data[ci]['width']=int(l)
                cusp_data[ci]['stabilizer']=S
            else:
                coseti=vertex_data[j]['coset'][0]
                permv=coset_perms[coseti]
                if p==0 and q==0:
                    VV=coset_reps[coseti]
                else:
                    VV=U*coset_reps[coseti]
                cusp_data[ci]['normalizer']=VV                    
                l = vertex_data[j]['width']
                cusp_data[ci]['width']=int(l)
                S = VV*SL2Z_elt(1,l,0,1)
                S = S._mul(VV,2)
                cusp_data[ci]['stabilizer']=S
            ci+=1
                # mul_list_maps(S,VV,2)
                #'vertices':[j]}                
                #vertex_data[j]['cusp']=ci
                #vertex_data[j]['cusp_map']=U #[U[0],U[1],U[2],U[3]]
                #cusps.append((p,q))
                #else:
                # Use the one we have.
                # First check if it is equiv to previous one
                # coseti=vertex_data[j]['coset']
                # permv=coset_perms[coseti]
                # VV=coset_reps[coseti]
                # if self._verbose>0:
                #     print "VV=",VV
                # is_eq=0
                # for cii in range(len(cusps)):
                #     W=cusp_data[cii]['normalizer']
                #     #permw=self.permutation_action(W,1)
                #     if not norm_perm.has_key(cii):
                #         norm_perm[cii]=self.permutation_action(W,1)
                #     permw=norm_perm[cii]
                #     #l=cusp_data[cii]['width']
                #     for k in all_cusp_widths:
                #         for ep in [1,-1]:
                #             if (permw.inverse())(((permT**(ep*k))(permv(1))))==1:
                #                 U = VV*SL2Z_elt(1,ep*k,0,1) 
                #                 U = U._mul(W,2) 
                #                 #U = VV*SL2Z([1,ep*k,0,1])*W**-1
                #                 if self._verbose>0:
                #                     print vertices[j]," is equiv to ",cii
                #                 vertex_data[j]['cusp']=cii
                #                 vertex_data[j]['cusp_map']=U.inverse() 
                #                 is_eq=1
                #                 break
                #         if is_eq==1:
                #             break
                #     if is_eq==1:
                #         break
                # if is_eq==1:
                #     continue
                #if self._verbose>0:
                #    print "2 setting cusp ",p,q
                #    print "width=",k
                #if k in all_cusp_widths:
                #    all_cusp_widths.remove(k)
                # cusps.append(v)
                #print "Adding to cusps",v
                #vertex_data[j]['cusp']=ci
                # Need to get all the data for this cusp
                # testi=permv(1)
                # try:
                #     for k0 in all_cusp_widths: #range(1,gen_level+1):
                #         for ep in [1,-1]:
                #             k = k0*ep
                #             if self._verbose>2:
                #                 print "trying width:",k
                #                 print "(~permw)((self.permT**l)(permv(1)))",(permv*self._my_permT**k*(permv.inverse())).to_cycles()
                #                 print "(~permw)((self.permT**l)(permv(1)))",(permv.inverse()*self._my_permT**k*(permv)).to_cycles()       
                #             test=(permv.inverse())((self._my_permT**k)(testi))
                #             if test==1:
                #                 if self._verbose>0:
                #                     print "permv=",permv
                #                 width=k
                #                 S = VV*SL2Z_elt(1,width,0,1)
                #                 S = S._mul(VV,2)
                #mul_list_maps(S,VV,2)
                #S=VV*SL2Z([1,width,0,1])*VV**-1
                #cusp_data[ci]={'normalizer':VV,'width':int(abs(width)),'stabilizer':S,'vertices':[j]}
                #if self._verbose>0:
                #                    print "Setting cusp:",cusps[ci]
                #                    print "width=",abs(width)
                #                raise StopIteration()
                #    raise ArithmeticError,"Could not compute data for cusp v={0}".format(v)
                #except StopIteration:
                #    pass
                #print "cusp_data[",ci,"]=",cusp_data[ci]
                #if abs(width) in all_cusp_widths:
                #    all_cusp_widths.remove(abs(width))
                #ci+=1
        if self._verbose>0:
            print "vertices",vertices
            print "vertex_data=",vertex_data
            print "cusps=",cusps
            print "cusp_data=",cusp_data
            print "coset+reps=",coset_reps
            #continue
            # Small test:
        #if len(cusps)<>len(cusp_widths):
        #    print "cusps=",cusps
        #    print "widths=",cusp_widths
        #    raise ArithmeticError,"Could not compute cusp data!"
        return vertices,vertex_data,cusps,cusp_data

    def PermGroup(self):
        if self._perm_group==None:
            self._perm_group = SymmetricGroup(self._index)
        return self._perm_group
        

    # # Now to public routines

    def coset_representative(self,A):
        r"""
        Indata: A in PSL(2,Z) 
        Returns the coset representative of A in
        PSL(2,Z)/self.G

        EXAMPLES::


            sage: G=MySubgroup(Gamma0(4))        
            sage: A=SL2Z([9,4,-16,-7])
            sage: G.coset_rep(A)
            [1 0]
            [0 1]
            sage: A=SL2Z([3,11,-26,-95])
            sage: G.coset_rep(A)
            [-1  0]
            [ 2 -1]
            sage: A=SL2Z([71,73,35,36])
            sage: G.coset_rep(A)
            [ 0 -1]
            [ 1  0]

        
        """
        for V in (self.coset_reps()):
            if  V*A**-1 in self:
                return V
        raise ArithmeticError,"Did not find coset rep. for A=%s" %(A)

    def pullback(self,x_in,y_in=None,ret_mat=0,prec=201,version = 0,**kwds):
        r""" Find the pullback of a point in H to the fundamental domain of self
        INPUT:

         - ''x_in,y_in'' -- x_in+I*y_in is in the upper half-plane
         - ''prec''      -- (default 201) precision in bits
         - ret_mat  -- set to 0 if you want to return a list instead of a matrix.
         - 'version' -- Integer, 0, 1, or 2. Determine which verson of the fundamental domain to use.
        OUTPUT:
        
         - [xpb,ypb,B]  --  xpb+I*ypb=B(x_in+I*y_in) with B in self
                           xpb and ypb are complex numbers with precision prec 
        EXAMPLES::


            sage: S=SymmetricGroup(6)
            sage: pS=S([2,1,4,3,6,5])
            sage: pR=S([3,1,2,5,6,4])
            sage: G=MySubgroup(o2=pS,o3=pR)
            sage: [x,y,B]=G.pullback(0.2,0.5,53); x,y;B
            (-0.237623762376238, 0.123762376237624)
            [-1  0]
            [ 4 -1]
            sage: (B**-1).acton(x+I*y)
            0.200000000000000 + 0.500000000000000*I


        """
        use_dp = False; use_mpfr = False; use_mat = False
        # First get the input point
        if isinstance(x_in,float):
            x = RR(x_in); y = RR(y_in); use_dp=True
        elif isinstance(x_in,RealNumber_class):
            RF = x_in.parent(); prec = RF.prec()
            x = RF(x_in); y = RF(y_in)
            if prec<=53:
                use_dp = True
            else:
                use_mpfr = True
        else:
            try:
                if y_in==None:
                    z = x_in
                    y = z.imag(); x=z.real()
                else:
                    z = x_in+I*y_in
                    x = x_in; y = y_in

                x = RR(x); y = RR(y)
                use_mat = True
                if y_in < 0.1:
                    prec = max(53,3*log_b(1.0/y,2))
                    RF = RealField(prec)
                    use_mpfr = True
                else:
                    RF = RR
                    use_dp = True
            except TypeError as er:
                raise TypeError,"Could not get point in upper half-plane from {0}! {1}!".format(x_in,er)

        if self._is_Gamma0 and version == 0:
            if use_dp:
                xpb,ypb,a,b,c,d=pullback_to_Gamma0N_dp(self,x,y,self._verbose)
            elif use_mpfr:
                xpb,ypb,a,b,c,d=pullback_to_Gamma0N_mpfr(self,x,y)
            else:
                raise ValueError," Need to use either dp or mpfr!"
        else:
            a,b,c,d=pullback_to_psl2z_mat(RR(x_in),RR(y_in))
            A=SL2Z_elt(a,b,c,d) #.matrix()
            reps = self.coset_reps(version=version)
            try:
                for V in reps:
                    B=V*A
                    if B in self:
                        raise StopIteration
            except StopIteration:            
                pass
            else:
                raise ArithmeticError,"Did not find coset rep. for A=%s" % A
            if use_dp:
                xpb,ypb=apply_sl2z_map_dp(x_in,y_in,B[0],B[1],B[2],B[3])
            else:
                xpb,ypb=apply_sl2z_map_mpfr(x_in,y_in,B[0],B[1],B[2],B[3])
            a,b,c,d = B[0],B[1],B[2],B[3]

        if use_mat:
            return (a*z+b)/(c*z+d)
        if ret_mat==1:
            return xpb,ypb,SL2Z_elt(a,b,c,d)
        else:
            return xpb,ypb,int(a),int(b),int(c),int(d)


    def is_congruence(self):
        r""" Is self a congruence subgroup or not?

        EXAMPLES::


            sage: S=SymmetricGroup(6)
            sage: pS=S([2,1,4,3,6,5])
            sage: pR=S([3,1,2,5,6,4])
            sage: G=MySubgroup(o2=pS,o3=pR)
            sage: G.is_congruence()
            True
            sage: S=SymmetricGroup(7)
            sage: pS=S([1,3,2,5,4,7,6]); pS
            (2,3)(4,5)(6,7)
            sage: pR=S([3,2,4,1,6,7,5]); pR
            (1,3,4)(5,6,7)
            sage: G=MySubgroup(o2=pS,o3=pR)
            sage: G.is_congruence()
        False
        
        """
        if self._is_congruence == None:
            self._is_congruence=super(MySubgroup_class,self).is_congruence()
        return self._is_congruence


    def is_Hecke_triangle_group(self):
        if self.is_Gamma0() and self.generalised_level()==1:
            return True
        else:
            return False
        
    def generalised_level(self):
        r""" Generalized level of self

        EXAMPLES::y

            sage: S=SymmetricGroup(6)
            sage: pS=S([2,1,4,3,6,5])
            sage: pR=S([3,1,2,5,6,4])
            sage: G=MySubgroup(o2=pS,o3=pR)
            sage: G.generalised_level()
            4
            sage: S=SymmetricGroup(7)
            sage: pS=S([1,3,2,5,4,7,6]); pS
            (2,3)(4,5)(6,7)
            sage: pR=S([3,2,4,1,6,7,5]); pR
            (1,3,4)(5,6,7)
            sage: G=MySubgroup(o2=pS,o3=pR)
            sage: G.generalised_level()
            6
        """
        if self._generalised_level==None:
            # compute the generalized level
            self._generalised_level= lcm(map(len,self.permT.cycle_tuples()))
        return self._generalised_level
    #raise ArithmeticError, "Could not compute generalised level of %s" %(self)


    def as_named_group(self):
        r"""
        Identify self as a group of type Gamma0, Gamma1 or Gamma
        if possible.-
        TODO: Check if we are conjugate to one of these groups
        """
        if not self.is_congruence():
            return None
        conj=1
        N = self.generalised_level()
        # do a preliminary test to see if ST^NS or T^N is in the group
        name,g = self.find_name(get_named_group=True)
        if name<>'':
            self._name = name
        if g<>None:
            return g
        
    def  find_name(self,get_named_group=False):
        r"""
        Try to find the name of self. So far only some easy types are implemented.
    
        """
        name = ''; G=None
        if self.is_congruence():
            N = self.generalised_level()
            mu = self.index()
            if mu == self._gamma0N_index(N):
                G = Gamma0(N)
                if self.is_subgroup(G):
                    name = 'Gamma_0({0})'.format(N)
                    latex_name = '\Gamma_0({0})'.format(N)                    
                else: ## We could have Gamma^0(N) = S Gamma_0(N) S 
                    S = SL2Z_elt(0,-1,1,0)
                    t = [S*x*S in G for x in self.generators_as_slz_elts()].count(False)
                    if t==0:
                        name = 'Gamma^0({0})'.format(N); G=None ## We can't get Gamma^0 in sage
                        latex_name = '\Gamma^0({0})'.format(N)
            if name=='' and mu == self._gammaN_index(N):
                G = Gamma(N).as_permutation_group().to_even_subgroup()
                name = 'Gamma({0})'.format(N)
                latex_name = '\Gamma({0})'.format(N)                
            if name== '' and mu == self._gamma1N_index(N):
                G = Gamma1(N).as_permutation_group().to_even_subgroup()
                t = [x.c() % N == 0 and (x.d()% N == 1 or -x.d() % N ==1) for x in self.generators_as_slz_elts()].count(False)
                if t == 0:
                    name = 'Gamma_1({0})'.format(N)
                    latex_name = '\Gamma_1({0})'.format(N)
            if name == '': # Try some other alternatives
                # Try to find if we have a Gamma_0(N,M)
                bs = map( lambda x:x.b(), self.generators_as_slz_elts())
                if self._verbose>0:
                    print "bs=",bs
                M = gcd(bs)
                # We know that Gamma_0^0(N,M) ~= B_n Gamma_0(NM) B_n^-1
                for m in M.divisors():
                    if mu == self._gamma0N_index(N*m): ## possible
                        #check if x in self satisfy 
                        t = [x.c() % N == 0 and x.b()% m == 0 for x in self.generators_as_slz_elts()].count(False)
                        if t== 0:
                            name = "Gamma_0^0({0},{1})".format(N,m)
                            latex_name = "\Gamma_0^0({0},{1})".format(N,m)
                            break
                
                if name=='' and mu==N and self.ncusps()==1 and mu<=3: # for larger mu there is more than one conjugate group of this type.
                    name='Gamma^{0}'.format(N)
                    latex_name = '\Gamma^{0}'.format(N)
        if name=='':
            G = None
        else:
            self._name = name
            self._latex_name = latex_name            
        if get_named_group:
            return name,G
        return name

    def latex_name(self):
        s = self.find_name()
        return self._latex_name

    
    
    def _gammaN_index(self,N):
        if N==2:
            return 6
        i = QQ(N)**3/QQ(2)
        for p in ZZ(N).prime_divisors():
            i*=(1-p**-2)
        return i
        
    def _gamma1N_index(self,N):
        if N==2:
            return 3
        i = QQ(N)**2/QQ(2)
        for p in ZZ(N).prime_divisors():
            i*=(1-p**-2)
        return i
    
    def _gamma0N_index(self,N):
        if N==2:
            return 3
        i = QQ(N)
        for p in ZZ(N).prime_divisors():
            i*=(1+p**-1)
        return i
 
    def get_equivalent_cusp(self,p,q):
        r"""
        Get an equivalent cusp and a cusp-normalizing map
        which is, if possible, an Atkin-Lehner involution.
        """

        if not self.is_Gamma0():
            return 0,0,0,0,0
        Id = SL2Z_elt(1,0,0,1)
        if q==0:
            WQ=Id
            A =Id
            l = 1
            return (WQ,A,p,q,l)
        if p==0:
            WQ=SL2Z_elt(0,-1,1,0)
            A =Id 
            l = self.generalised_level()
            return (WQ,A,p,q,l)
        N = self.generalised_level()
        if self._verbose>0:
            print "p=",p
            print "q=",q
            print "N=",N
        if N % q==0:
            Q = N.divide_knowing_divisible_by(q)*p
            if Q<0:
                Q=-Q
                #q=abs(q)
                #p=abs(p)
            if Q.divides(N):
                if self._verbose>0:
                    print "Q0=",Q
                N1 = N.divide_knowing_divisible_by(Q)
                g,s,t=xgcd(Q,N1)
                if g==1:
                    A=[1,0,0,1]
                    if self._verbose>0:
                        print "s,t=",s,t
                    # a = Q; b =-t; c=N; d=s*Q                
                    if Q>0:
                        WQ=SL2Z_elt(1,-t,N1,Q*s)
                    else:
                        WQ=SL2Z_elt(-1,t,-N1,-Q*s)
                    if Q<0:
                        return (WQ,A,abs(p),abs(q),abs(Q))
        divs=[]

        #for Q in self._level.divisors():
        #    divs.append(Q)
        for Q in N.divisors():
            divs.append(-Q)
            #divs.append(-Q)
        res=[]
        for Q in divs:
            if self._verbose>0:
                print "Q=",Q
            c1 = Cusp(QQ(Q)/QQ(N)) #Cusp(QQ(Q)/QQ(N))
            c2 = Cusp(QQ(p)/QQ(q))
            t,A=Cusp.is_gamma0_equiv(c1, c2, N, transformation='matrix')
            if self._verbose>0:
                #print "c1=",c1
                #print "c2=",c2
                print "t=",t,A
            if t:
                U=SL2Z_elt(A[1,1],-A[0,1],-A[1,0],A[0,0])
                #U=A**-1
                if self._verbose>0:
                    print "c1=",c1
                    print "c2=",c2
                pp=c1.numerator()
                qq=c1.denominator()
                g,s,t=xgcd( Q,ZZ(QQ(N)/QQ(Q)))
                if g<>1:
                    # We can not get Atkin-Lehner but at least we get
                    # a canonical representative of the form Q/N
                    # and if we have a cusp of the form p/q, q|N
                    # (p,q)=1 then we get the width from Iwaniec "Topics in Classical..., Sect. 2.4"
                    # for Gamma_0(N) at least
                    if self._is_Gamma0:
                        if q % N==0:
                            l = N.divide_knowing_divisible_by(gcd(N,q**2))
                        elif qq % N==0:
                            l = N.divide_knowing_divisible_by(gcd(N,qq**2))
                        else:
                            if self._verbose>0:
                                print "q%N=",(q%N)
                                print "qq%N=",(qq%N)
                            l=0 # self._G.cusp_width(c2)
                    else:
                        l=self.cusp_width(c2)
                    ## Return a non-normalizer (but maybe twisted) map W
                    if pp==1:
                        W=SL2Z_elt(1,0,int(qq),1)
                    else:
                        w = lift_to_sl2z(q, p, 0 )
                        W = SL2Z_elt(w[3 ], w[1 ], w[2 ],w[0 ])
                    return W,U,pp,qq,l
                a = Q; b =-t; c=N; d=s*Q
                if self._verbose>0:
                    print "a,b,c,d=",a,b,c,d
                ## We should see if (a b ; c d) = A rho with A in SL2Z
                for l in gcd(a,c).divisors():
                    aa = a.divide_knowing_divisible_by(l)
                    cc = c.divide_knowing_divisible_by(l)
                    if aa*d-cc*b==1:
                        WQ=SL2Z_elt(aa,b,cc,d)
                        pp=numerator(Q/N)
                        qq=denominator(Q/N)
                        #if c1==c2:
                        return (WQ,U,pp,qq,l)
                        #res.append((WQ,U,pp,qq,l))
                        break
        if len(res)==1:
            return res[0]
        if len(res)>1:
            # Pick the best one
            if self._verbose>0:
                print "res>1=",res
            lista= map(lambda x:sum(map(abs,x[1].matrix().list())),res)
            mina=min(lista)
            i=lista.index(mina)
            return res[i]
        return 0,0,0,0,0


    def show_symmetry_props(self):
        r"""
        Display the symmetry properties of self.
        """
        for i in range(self._ncusps):
            print "Cusp Nr. {0} = {1}:{2}".format(i,self._cusps[i][0],self._cusps[i][1])
            print "Normalizer is normalizer : {0}".format(self.cusp_normalizer_is_normalizer(i))
            print "Can be symmetrized even/odd: {0}".format(self.is_symmetrizable_even_odd(i))



    def test_cusp_data(self):
        for i in range(self._ncusps):
            p,q = self._cusps[i]
            N = self._cusp_data[i]['normalizer']
            for j in range(self._ncusps):
                if i==j:
                    continue
                x=Cusp(self._cusps[j])
                if self.are_equivalent(Cusp(p,q),x):
                    raise ArithmeticError, "Cusps {0} and {1} are equivalent!".format(Cusp(p,q),x)
                
            if q==0:
                c = Infinity
                t1 = N[0]==p and N[2]==0
            else:
                c=QQ(p)/QQ(q)
                t1 = c==QQ(N[0])/QQ(N[2])
            if not t1:
                raise ArithmeticError, "{0} is not normalizer of cusp nr.{1}: {2}".format(N,i,c)
            A = self._cusp_data[i]['stabilizer']
            if q==0:
                t2 = A[0]<>0 and A[2]==0
            else:
                t2 = c== acton_list_maps(A,c) 
            if not t2:
                raise ArithmeticError, "{0} is does not stabilize cusp nr. {1}: {2}".format(A,i,c)
            TT = mul_list_maps(N,A,1)
            TT = mul_list_maps(TT,N)
            #TT = N**-1*A*N
            S = [0,-1,1,0]; T=[1,1,0,1]
            #S,T=SL2Z.gens()
            l = self._cusp_data[i]['width']
            if TT<>[1,l,0,1]:
                raise ArithmeticError, "Width {0} is not correct width for cusp nr. {1}: {2}".format(l,i,c)
        for i in range(self._nvertices):
            vp,vq = self._vertices[i] 
            if vq==0:
                v = Infinity
            else:
                v = QQ(vp)/QQ(vq)
            j = self._vertex_data[i]['cusp']
            p,q = self._cusps[j]
            if q==0:
                c=Infinity
            else:
                c = QQ(p)/QQ(q)
            # Want to see that it maps to the correct cusp
            A = self._vertex_data[i]['cusp_map']
            if vq==0:
                if q==0:
                    t3 = A[0]<>0 and A[2]==0
                else:
                    t3 = QQ(A[0])/QQ(A[2])==c
            else:
                if q==0:
                    t3 = (A[0]*v+A[2])<>0 and  (A[2]*v+A[3])==0
                else:
                    #print "Av=",A.acton(v)
                    #print "c=",c
                    t3 = acton_list_maps(A,v)==c
            if not t3:
                raise ArithmeticError,"The map {0} does not map the vertex nr. {1}: {2} to the cusp nr. {3}: {4}".format(A,i,v,j,c)
            # Also check that the correct number of coset reps maps to this vertex
            cnt=0
            for vj in self._vertex_data[i]['coset']:
                V = self.coset_reps()[vj]
                a = V.a(); c=V.c()
                if Cusp(a,c)<>Cusp(vp,vq):
                    raise ArithmeticError,"Coset rep {0} do not map to the vertex {1}!".format(V,v)
                else:
                    cnt+=1
            if cnt>self._vertex_data[i]['width']:
                raise ArithmeticError,"Too many coset reps. maps to the vertex {0}! ".format(v)

    def test_coset_reps(self,version=0):
        r"""
        Test if the coset representatives are independent.
        """
        for j in range(len(self.coset_reps(version))):
            A=self.coset_reps(version)[j]
            #print "A=",A,type(A)
            for k in range(len(self.coset_reps(version))):
                B=self.coset_reps(version)[k]
                #print "B=",B,type(B)
                #if A==B:
                if A[0]==B[0] and A[1]==B[1] and A[2]==B[2] and A[3]==B[3]:
                    continue
                BB = mul_list_maps(A,B,2)
                if self._is_Gamma0:
                    t = BB[2] % self._level == 0
                else:
                    t = BB in self
                if t:
                    print "A[",j,"]=",A
                    print "B[",k,"]=",B
                    print "AB^-1=",BB
                    raise ArithmeticError,"Coset representatives are not correct!"
            ## We also want sigma(A)(1)=j if we have our own coset reps. (i.e. version 0 or 1)
            k = self.permutation_action(A)(1)
            if k<>j+1:
                raise ArithmeticError,"sigma(V_{0})(1)={1}".format(j,k)
        return True

    def test_normalize_pt_to_cusp(self,z,ci):
        assert isinstance(ci,(int,Integer))
        assert z.imag()>0
        N=self._cusp_data[ci]['normalizer']
        l=self._cusp_data[ci]['width']
        print "N,l=",N,l
        z1=N.acton(z*RR(l))
        print "z1=",z1
        x,y,A= self.pullback(RR(z1.real()),RR(z1.imag()))
        z2=x+I*y
        print "z2=",z2
        cj,vj=self.closest_cusp(x,y,1)
        print "cj,vj=",cj,vj
        U = self._cusp_maps[vj]
        z3=U.acton(z2)
        print "z3=",z3
        N2=self._cusp_data[cj]['normalizer']
        l2=self._cusp_data[cj]['width']
        z4 = (N2**-1).acton(z3)/RR(l2)
        return z4

    def _test_consistency_perm(self,o2,o3):
        r""" Check the consistency of input permutations.

        INPUT:
        - ''o2'' Permutation of N1 letters
        - ''o3'' Permutation of N2 letters
        
        OUTPUT:
        - True : If N1=N2=N, o2**2=o3**3=Identity, <o2,o3>=Permutations(N)
                 (i.e. the pair o2,o3 is transitive) where N=index of self
        - Otherwise an Exception is raised


        EXAMPLES::


            sage: P=Permutations(6)
            sage: pS=P([2,1,4,3,6,5])
            sage: pR=P2([3,1,2,5,6,4])
            sage: G=MySubgroup(o2=pS,o3=pR)     
            sage: G._test_consistency_perm(pS,pR)
            True

        """
        #SG1=MyPermutation(length=self._index) #self._S.identity()
        #o22=o2*o2; o33=o3*o3*o3
        if self._index==1:
            return True

        if o2.order() not in [1,2] or o3.order() not in [1,3]:
            print "o2=",o2,o2.order()
            print "o3=",o3,o3.order()
            s="Input permutations are not of order 2 and 3: \n perm(S)^2={0} \n perm(R)^3={1} " 
            raise ValueError, s.format(o2*o2,o3*o3*o3)
        if not are_transitive_permutations(o2,o3):
            S = SymmetricGroup(self._index)
            G = S.subgroup([S(self.permS.list()),S(self.permR.list())])
            #GS=self._S.subgroup([self.permS,self.permR])
            s="Error in construction of permutation group! Permutations not transitive: <S,R>={0},{1}. Generated group={2}".format(o2,o3,G)
            raise ValueError, s
        return True

class MySubgroup_congruence_class (MySubgroup_class):
    r"""
    Subclass of congruence subgroups.
    """
    def __init__(self,o2=None,o3=None,verbose=0,display_format='short',data={},**kwds):
        super(MySubgroup_congruence_class,self).__init__(o2,o3,verbose,display_format,data,**kwds)
        if not self.is_congruence():
            raise ValueError,"Self is not a congruence subgroup!"
        self.class_name='MySubgroup_congruence_class'            

    def level(self):
        return self.generalised_level()
 
    
### Hecke triangle groups
class HeckeTriangleGroup(SageObject):
    def __init__(self,q=3,prec=53,verbose=0):
        r"""
         Hecke triangle groups.
         NOTE: For now we don't have subgroups.
        """
        if not isinstance(q,(int,Integer)):
            raise ValueError,"Need an integer q!"
        self._q=q
        self._prec=prec
        self._verbose=verbose
        RF=RealField(prec)
        self._lambdaq=RF(2)*(RF.pi()/RF(q)).cos()
        self._coset_reps = [1]
        self._vertices   = [(1,0)]  # List of vertices corresponding to the coset reps.
        self._cusps      = [(1,0)] # Subset of representatives
        self._cusps_as_cusps=[Cusp(1,0)]
        self._S=None
        self._R=None
        self._T=None
        self._index=1
        self._ncusps=1
        self._nvertices=1
        #A = Matrix(ZZ,2,2,[1,0,0,1])
        self._vertex_data=dict()
        self._level=1  # Should be lambdaq?
        self._vertex_data[0]={'cusp':0,'cusp_map':[1,0,0,1],'coset':[0]}
        self._cusp_data=dict()
        Tl=matrix(RR,2,2,[1,self._lambdaq,0,1])
        self._cusp_data[0]={'width':self._lambdaq,'normalizer':[1,0,0,1],'stabilizer':Tl}
        self._vertex_widths=dict()
        self._vertex_widths=[self._lambdaq]
        self._vertex_maps=[[1,0,0,1]]
        self._cusp_maps=[[1,0,0,1]]
        self._is_Hecke_triangle_group=True
        self._is_Gamma0=False
        self._is_congruence=True
        self._symmetrizable_cusp={0:1}
        self._cusp_normalizer_is_normalizer={0:(1,1)}
        self.class_name='HeckeTriangleGroup'       

    def __repr__(self):
        return 'Hecke Triangle group G_{0}'.format(self._q)
        
    def closest_cusp(self,x,y,vertex=0):
        if vertex==1:
            return (1,0),0
        else:
            return (1,0)
        
    def ncusps(self):
        return self._ncusps

    def nu2(self):
        return 1

    def nu3(self):
        if self._q<>3:
            return 0
        else:
            return 1

    def is_congruence(self):
        return True  # It is a congruence group of itself...

    def is_Gamma0(self):
        return False

    def is_Hecke_triangle_group(self):
        return True
    
    def __rep__(self):
        s="Hecke triangle group G_{0}".format(self._q)
        return s

    def level(self):
        return self._level

    def index(self):
        return self._index

    def cusps(self):
        return self._cusps

    def minimal_height(self):
        return (RR.pi()/self._q).sin()
    
    def pullback(self,x_in,y_in,prec=201 ):
        r""" Find the pullback of a point in H to the fundamental domain of self
        INPUT:

         - ''x_in,y_in'' -- x_in+I*y_in is in the upper half-plane
         - ''prec''      -- (default 201) precision in bits

        OUTPUT:
        
         - [xpb,ypb,B]  --  xpb+I*ypb=B(x_in+I*y_in) with B in self
                           xpb and ypb are complex numbers with precision prec 
        EXAMPLES::



        """
        x=deepcopy(x_in); y=deepcopy(y_in)
        RF=RealField(prec)
        a,b,c,d=pullback_to_hecke_triangle_mat_dp(x,y,RR(self._lambdaq))
        xpb=RR(x_in); ypb=RR(y_in)
        #print a,b,c,d
        xpb,ypb=apply_sl2r_map_dp(x_in,y_in,a,b,c,d)
        return [xpb,ypb,[a,b,c,d]] #Matrix(RR,2,2,[a,b,c,d])



# def factor_matrix_in_sl2z_in_S_and_T(A_in):
#         r"""
#         Factor A in SL2Z in generators S=[[0,-1],[1,0]] and T=[[1,1],[0,1]]
#         INPUT:
#          - ''A_in'' -- Matrix in SL2Z
#         OUTPUT:
#          - ''[[a0,a1,...,an],ep] with a0 in ZZ, ai in ZZ \ {0,\pm1}, ep=1,-1
#            Determined by:
#            A=ep*(T^a0*S*T^a1*S*...*S*T^an)
#            The algorithm uses the Nearest integer Continued Fraction expansion. Note that there due to relations in SL2Z the sequence a_j is not unique
#            (see example below)

#         EXAMPLES::
#             sage: nearest_integer_continued_fraction(0.5);cf
#             [1, 2]
#             sage: A=ncf_to_matrix_in_SL2Z(cf); A
#             [1 1]
#             [1 2]
#             sage: factor_matrix_in_sl2z_in_S_and_T(A)
#             [[1, 2], 1]

#         An example where we do not get back the same sequence::

#             sage: cf=nearest_integer_continued_fraction(pi.n(100),nmax=10);cf
#             [3, -7, 16, 294, 3, 4, 5, 15, -3, -2, 2]
#             sage: A=ncf_to_matrix_in_SL2Z(cf); A
#             [ -411557987 -1068966896]
#             [ -131002976  -340262731]
#             sage: factor_matrix_in_sl2z_in_S_and_T(A)
#             [[3, -7, 16, 294, 3, 4, 5, 15, -2, 2, 3], -1]
        
#         """
#         if A_in not in SL2Z:
#             raise TypeError, "%s must be an element of SL2Z" %A_in
#         S,T=SL2Z.gens()
#         # If A is not member of SL(2,Z) but a plain matrix
#         # we cast it to a member of SL2Z
#         A = SL2Z(A_in)
#         if(A.matrix()[1 ,0 ] == 0 ):
#             return [[A.matrix()[0 ,1 ]],1 ]
#         x = Rational(A.matrix()[0 ,0 ] / A.matrix()[1 ,0 ])
#         cf = nearest_integer_continued_fraction(x)
#         B  = ncf_to_matrix_in_SL2Z(cf)
#         # we know that A(oo) = x = B*S (oo) and that A and BS differ at most by a translation T^j
#         Tj = S**-1  * B**-1  * A 
#         sgn=1 
#         if(Tj.matrix()[0 ,0 ]<0 ):
#             j=-Tj.matrix()[0 ,1 ]
#             sgn=-1 
#             #Tj=SL2Z([1,j,0,1])
#         else:
#             j=Tj.matrix()[0 ,1 ]
#         # if Tj = Id i.e. j=0  then A=BS otherwise A=BST^j
#         cf.append(j)
#         # To make sure we test
#         C = B*S*Tj
#         try:
#             for ir in range(2):
#                 for ik in range(2):
#                     if(C[ir,ik]<>A[ir,ik] and C[ir,ik]<>-A[ir,ik]):
#                         raise StopIteration
#         except StopIteration:
#             print "ERROR: factorization failed %s <> %s " %(C,A)
#             print "Cf(A)=",cf
#             raise ArithmeticError," Could not factor matrix A=%s" % A_in
#         if(C.matrix()[0 ,0 ]==-A.matrix()[0 ,0 ] and C.matrix()[0 ,1 ]==-A.matrix()[0 ,1 ]):
#             sgn=-1 
#         return [cf,sgn]





# def list_factor_matrix_in_sl2z_in_S_and_T(A):
#     r"""
#     Factor A in SL2Z in generators S=[[0,-1],[1,0]] and T=[[1,1],[0,1]]
#     INPUT:
#     - ''A_in'' -- Matrix in SL2Z
#     OUTPUT:
#     - ''[[a0,a1,...,an],ep] with a0 in ZZ, ai in ZZ \ {0,\pm1}, ep=1,-1
#     Determined by:
#     A=ep*(T^a0*S*T^a1*S*...*S*T^an)
#     The algorithm uses the Nearest integer Continued Fraction expansion. Note that there due to relations in SL2Z the sequence a_j is not unique
#     (see example below)
    
#     EXAMPLES::
#     sage: nearest_integer_continued_fraction(0.5);cf
#     [1, 2]
#         sage: A=ncf_to_matrix_in_SL2Z(cf); A
#         [1 1]
#         [1 2]
#         sage: factor_matrix_in_sl2z_in_S_and_T(A)
#         [[1, 2], 1]

#     An example where we do not get back the same sequence::

#         sage: cf=nearest_integer_continued_fraction(pi.n(100),nmax=10);cf
#         [3, -7, 16, 294, 3, 4, 5, 15, -3, -2, 2]
#         sage: A=ncf_to_matrix_in_SL2Z(cf); A
#         [ -411557987 -1068966896]
#         [ -131002976  -340262731]
#         sage: factor_matrix_in_sl2z_in_S_and_T(A)
#         [[3, -7, 16, 294, 3, 4, 5, 15, -2, 2, 3], -1]

#     """
#     #if A_in not in SL2Z:
#     if not isinstance(A,list) or not (A[0]*A[3]-A[1]*A[2]==1):
#         raise TypeError, "%s must be an element of SL2Z" %A_in
#     S = [0,-1,1,0]; T=[1,1,0,1]
#     #S,T=SL2Z.gens()
#     # If A is not member of SL(2,Z) but a plain matrix
#     # we cast it to a member of SL2Z
#     #A = SL2Z(A_in)
#     if(A[2] == 0 ):
#         return [[ b ],1 ]
#     x = Rational(A[0] / A[2])
#     cf = nearest_integer_continued_fraction(x)
#     B  = list_ncf_to_matrix_in_SL2Z(cf)
#     # we know that A(oo) = x = B*S (oo) and that A and BS differ at most by a translation T^j
#     Tj = mul_list_maps(B,A,1)
#     Tj = mul_list_maps(S,Tj,1)
#     #Tj = S**-1  * B**-1  * A 
#     sgn=1 
#     if Tj[0]<0:
#         j=-Tj[1]
#         sgn=-1 
#         #Tj=SL2Z([1,j,0,1])
#     else:
#         j=Tj[1]
#     # if Tj = Id i.e. j=0  then A=BS otherwise A=BST^j
#     cf.append(j)
#     # To make sure we test
#     C = mul_list_maps(B,S)
#     C = mul_list_maps(C,Tj)
#     #C = B*S*Tj
#     try:
#         for ir in range(4):
#             if C[ir]<>A[ir] and C[ir]<>-A[ir]:
#                 raise StopIteration
#             #for ik in range(2):
#             #    if(C[ir,ik]<>A[ir,ik] and C[ir,ik]<>-A[ir,ik]):
#             #        raise StopIteration
#     except StopIteration:
#         print "ERROR: factorization failed %s <> %s " %(C,A)
#         print "Cf(A)=",cf
#         raise ArithmeticError," Could not factor matrix A=%s" % A_in
#     if (C[0]==-A[0] and C[1]==-A[1]):
#         sgn=-1 
#         return [cf,sgn]


# def nearest_integer_continued_fraction(x,nmax=None):
#     r""" Nearest integer continued fraction of x
#     where x is a rational number, and at n digits

#     EXAMPLES::

#         sage: nearest_integer_continued_fraction(0.5)
#         [1, 2]
#         nearest_integer_continued_fraction(pi.n(100),nmax=10)
#         [3, -7, 16, 294, 3, 4, 5, 15, -3, -2, 2]

#     """
#     if(nmax == None):
#         if(x in QQ):
#             nmax=10000 
#         else:
#             nmax=100  # For non-rational numbrs  we don't want so many convergents
#     jj=0 
#     cf=list()
#     n=nearest_integer(x)
#     cf.append(n)
#     if(x in QQ):
#         x1=Rational(x-n)
#         while jj<nmax and x1<>0  :
#             n=nearest_integer(-1 /x1)
#             x1=Rational(-1 /x1-n)
#             cf.append(n)
#             jj=jj+1 
#         return cf
#     else:
#         try:
#             RF=x.parent()
#             x1=RF(x-n)
#             while jj<nmax and x1<>0  :
#                 n=nearest_integer(RF(-1 )/x1)
#                 x1=RF(-1 )/x1-RF(n)
#                 cf.append(n)
#                 jj=jj+1 
#             return cf
#         except AttributeError:
#             raise ValueError,"Could not determine type of input x=%s" %x

# def nearest_integer(x):
#     r""" Returns the nearest integer to x: [x]
#     using the convention that 
#     [1/2]=0 and [-1/2]=0


#     EXAMPLES::

#         sage: nearest_integer(0)
#         0
#         sage: nearest_integer(0.5)
#         1
#         sage: nearest_integer(-0.5)
#         0
#         sage: nearest_integer(-0.500001)
#         -1
#     """
#     return floor(x+ 0.5)


def _get_perms_from_str(st):
    r""" Construct permutations frm the string st of the form
    st='[a1_a2_..._an]-[b1_...bn]'
    

    EXAMPLES::


    sage: G=MySubgroup(Gamma0(5))
    sage: s=G._get_uid();s
    '[2_1_4_3_5_6]-[3_1_2_5_6_4]-Gamma0(5)'
    sage: _get_perms_from_str(s)
    [(1,2)(3,4), (1,3,2)(4,5,6)]
    
    sage: S=SymmetricGroup(7)
    sage: pS=S([1,3,2,5,4,7,6]); pS
    (2,3)(4,5)(6,7)
    sage: pR=S([3,2,4,1,6,7,5]); pR
    (1,3,4)(5,6,7)
    sage: G=MySubgroup(Gamma0(4))
    sage: G=MySubgroup(o2=pS,o3=pR)
    sage: s=G._get_uid();s
    '[1_3_2_5_4_7_6]-[3_2_4_1_6_7_5]'
    sage: _get_perms_from_str(s)
    [(2,3)(4,5)(6,7), (1,3,4)(5,6,7)]
        
    """
    (s1,sep,s2)=s.partition("-")
    if(len(sep)==0 ):
        raise ValueError,"Need to supply string of correct format!"
    if(s2.count("-")>0 ):
        # split again
        (s21,sep,s22)=s2.partition("-")
        s2=s21
    ix=s1.count("_")+ZZ(1)
    pS=list(range(ix))
    pR=list(range(ix))
    s1=s1.strip("["); s1=s1.strip("]");
    s2=s2.strip("["); s2=s2.strip("]");
    vs1=s1.split("_")
    vs2=s2.split("_")
    P=SymmetricGroup(ix)
    for n in range(ix):
        pS[n]=int(vs1[n])
        pR[n]=int(vs2[n])
    return [P(pS),P(pR)]

def get_perms_from_cycle_str(s,N,sep=' ',deb=False):
    r""" Construct permutations frm the string st of the form
    st='[(a1;a2;...;an)(b1;...;bn)'
    where ';' can be replaced by any separator.
    
    EXAMPLES::

    """
    # First get the cycles (and remove the beginning and end of the string)
    cycles = s[1:-1].split(")(")
    perm = range(1,N+1)
    if deb:
        print "N=",N
    for c in cycles:
        if deb:
            print "c=",c
        l=map(int,c.split(sep)) # elements of the cycle are separated by "sep"
        nn = len(l)-1
        if nn>0:
            for i in range(nn):
                perm[l[i]-1]=l[i+1]
                perm[l[nn]-1]=l[0]        
    S=SymmetricGroup(N)
    return S(perm)


###
def level(cls):
    return cls.generalised_level()




### Algorithms which deal with lists of groups


def list_valid_signatures(index,nc=None,ne2=None,ne3=None,ng=None):
    r"""
    Returns a list of all signatures , i.e.tuples:
    (index,h,e2,e3,g)
    with
      h   = number of cusps
      e2 = number of order 2 elliptic points
      e2 = number of order 3 elliptic points
      g   = genus
    of subgroups of PSL(2,Z) of the given index.

    INPUT:
    
    - ''index'' -- a positive integer
    - ''nc''  -- number of cusps desired (default None)
    - ''ne2'' -- number of elliptic points of order 2 desired (default None)
    - ''ne3'' -- number of elliptic points of order 2 desired (default None)
    - ''ng''  -- genus desired (default None)
    
    EXAMPLES::


        sage: get_list_of_valid_signatures(2)
        [[2, 1, 0, 2, 0]]
        sage: get_list_of_valid_signatures(5)
        [[5, 1, 1, 2, 0]]
        sage: get_list_of_valid_signatures(6)
        [[6, 1, 0, 0, 1], [6, 1, 0, 3, 0], [6, 1, 4, 0, 0], [6, 2, 2, 0, 0], [6, 3, 0, 0, 0]]


    """
    if(index not in ZZ or index <=0):
        raise TypeError, "index must be a positive integer! Got: %s"%index
    res=list()
    for h in range(1,index+1):
        if(nc<>None and h<>nc):
            continue
        for e2 in range(0,index+1):
            if(ne2<>None and ne2<>e2):
                continue
            for e3 in range(0,index+1):
                if(ne3<>None and ne3<>e3):
                    continue
                g=12+index-6*h-3*e2-4*e3
                if(ng<>None and g<>ng):
                    continue
                if( (g>=0) and ((g % 12)==0)): # we have a valid type
                    g=g/12
                    res.append([index,h,e2,e3,g])
    return res




def perm_cycles(perm):
    c=perm.cycle_tuples()
    #if sum(map(lambda x:len(x),c))<>len(e.list()):
    # If there are fixed points we add them 
    l = perm.list()
    for i in range(len(l)):
        if l[i]==i+1:
            c.append((i+1,))
    return c

def mul_list_maps(l1,l2,inv=0):
    r"""
    A*B=   (a b) (x y)
           (c d) (z w)
    If inv == 1 invert A
    If inv == 2 invert B
    """
    a,b,c,d=l1
    x,y,z,w=l2
    if inv ==1:
        tmp=a
        a=d
        d=tmp
        b=-b; c=-c
    if inv==2:
        tmp=x
        x=w
        w=tmp
        y=-y; z=-z
    aa=a*x+b*z
    bb=a*y+b*w
    cc=c*x+d*z
    dd=c*y+d*w
    return [aa,bb,cc,dd]

def acton_list_maps(A,c):
    return (A[0]*c+A[1])/(A[2]*c+A[3])



def test_gamma0N(N1,N2):
    for N in range(N1,N2):
        print "N=",N
        G=MySubgroup(Gamma0(N))
        G.test_cusp_data()
    return 1

def is_integer(x):
    if isinstance(x,(int,Integer)):
        return True
    elif isinstance(x,Rational):
        if x.denominator()==1:
            return True
    elif isinstance(x,(float,RealNumber)):
        if ceil(x)==floor(x):
            return True
    return False

def is_Hecke_triangle_group(G):
    if not hasattr(G,"is_Hecke_triangle_group"):
        return False
    return G.is_Hecke_triangle_group()
