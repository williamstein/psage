# -*- coding: utf-8 -*-
r"""
A general class for subgroups of the modular group.
Extends the standard classes with methods needed for Maass waveforms.


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

from sage.all_cmdline import *   # import sage library
from sage.rings.all    import (Integer, CC,ZZ,QQ,RR,RealNumber,I,infinity,Rational)
from sage.combinat.permutation import (Permutations,PermutationOptions)
from sage.modular.cusps import Cusp
from sage.modular.arithgroup.all import *
from sage.modular.modsym.p1list import lift_to_sl2z 
from sage.functions.other import (floor,sqrt)
from sage.all import matrix #import constructor
from sage.modular.arithgroup import congroup_gamma0
from sage.groups.all import SymmetricGroup
from sage.rings.arith import lcm
from copy import deepcopy
from mysubgroups_alg import * #are_transitive_permutations



_sage_const_3 = Integer(3); _sage_const_2 = Integer(2);
_sage_const_1 = Integer(1); _sage_const_0 = Integer(0);
_sage_const_6 = Integer(6); _sage_const_4 = Integer(4);
_sage_const_100 = Integer(100);
_sage_const_1En12 = RealNumber('1E-12');
_sage_const_201 = Integer(201); _sage_const_1p0 = RealNumber('1.0');
_sage_const_1p5 = RealNumber('1.5'); _sage_const_0p0 = RealNumber('0.0');
_sage_const_2p0 = RealNumber('2.0'); _sage_const_0p2 = RealNumber('0.2');
_sage_const_0p5 = RealNumber('0.5'); _sage_const_1000 = Integer(1000);
_sage_const_10000 = Integer(10000); _sage_const_1En10 = RealNumber('1E-10');
_sage_const_20 = Integer(20); _sage_const_3p0 = RealNumber('3.0')




import sys,os
import matplotlib.patches as patches
import matplotlib.path as path

from sage.modular.arithgroup.arithgroup_perm import *
#from subgroups_alg import *
#load "/home/stromberg/sage/subgroups/subgroups_alg.spyx"

class MySubgroup (ArithmeticSubgroup):
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
    def __init__(self,G=None,o2=None,o3=None,str=None):
        r""" Init a subgroup in the following forms: 
          1. G = Arithmetic subgroup
          2. (o2,o3) = pair of transitive permutations of order 2 and
          order 3 respectively.
          
          Attributes:

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
        # Get index also check consistency of indata


        self._coset_reps = None
        self._vertices   = None  # List of vertices corresponding to the coset reps.
        self._cusps      = None # Subset of representatives
        self._cusp_representative = None # Associate cusp reps to vertices
        self._vertex_reps=None # The coset rep associated to each vertex
        if(str<>None):
            try:
                G=sage_eval(str)
            except:
                try:
                    [o2,o3]=self._get_perms_from_str(str)
                    o2.to_cycles()
                except:
                    raise ValueError,"Incorrect string as input to constructor! Got s=%s" %(str)

        self._index=self._get_index(G,o2,o3)
        self._level=None
        self._generalised_level=None
        self._is_congruence=None
        ## The underlying permutation group
        self._P=Permutations(self._index)
        #self._S=SymmetricGroup(self._index) 
        ## Check if the multiplication of permutations is set as left
        ## to right, which is what we expect. I.e. we need
        ## (f*g)(x)=g(f(x))
        ## With this conventions the map h:SL2Z->Permutations(N)
        ## is a homomorphis. I.e. h(AB)=h(A)h(B)
        if(PermutationOptions()['mult']<>'l2r'):
            raise ValueError, "Need permutations to multply from left to right! Got PermutationOptions()="  %PermutationOptions()
        PermutationOptions(display='cycle')
        if G <>None:
            self._G=G                 
            if hasattr(G,'permutation_action'):
                # is a permutation subgroup
                self.permS=self._P( (G.L*(~G.R)*G.L).list())
                self.permT=self._P( G.L.list()  )
                self.permR=self.permS*self.permT
                self.permP=self._P( G.R.list()  )
                self._coset_reps=self._get_coset_reps_from_perms(self.permS,self.permR,self.permT)
            elif(hasattr(G,'is_congruence')):
                self._coset_reps=self._get_coset_reps_from_G(G)
                [self.permS,self.permR]=self._get_perms_from_coset_reps()
                self.permT=self.permS*self.permR
                self.permP=self.permT*self.permS*self.permT
            else:
                raise TypeError,"Error input G= %s" % G
            ArithmeticSubgroup.__init__(G)                        
        elif( (o2<>None and o3<>None) or (str<>None)):
            if(o2 in self._P):
                self.permS=o2                
                self.permR=o3
            else:
                self.permS=self._P(o2.list())
                self.permR=self._P(o3.list())
            # Check the consistency of the input permutations
            self._test_consistency_perm(self.permS,self.permR)
            self.permT=self.permS*self.permR
            self.permP=self.permT*self.permS*self.permT
            #print "self.permS=",self.permS
            #print "self.permT=",self.permT
            #print "self.permR=",self.permR                
            ## Need to get rid of the necessity for this
            S=SymmetricGroup(self._index)
            L=S(list(self.permT))
            R=S(list(self.permP))
            G=ArithmeticSubgroup_Permutation(L,R)
            ArithmeticSubgroup.__init__(G)                        
            self._G=G
            self._coset_reps=self._get_coset_reps_from_perms(self.permS,self.permR,self.permT)
            ## Recall that generalized level is just the greatest common divisor of cusp widths
            ## or alternatively gcd of all cycle lenghts of sigma_T
            cycles=self.permT.to_cycles()
            lens=list()
            for c in cycles:
                lens.append(len(c))
            self._generalised_level=lcm(lens)
        else:
            raise TypeError,"Error incorrect input to constuctor!"
            # Figure out if we have permutations in some form
        ## Just to make sure that the permutations are implemented correctly
        self._nu2=len(self.permS.fixed_points())
        self._nu3=len(self.permR.fixed_points())
        #print "To get_vertices"
        [self._vertices,self._vertex_reps]=self._get_vertices(self._coset_reps)
        #print "GOT_vertices"
        self._cusps=self._get_cusps(self._vertices)
        self._ncusps=len(self._cusps)
        self._genus=1 +1 /_sage_const_2 *(self._index/_sage_const_6 -self._ncusps-1 /_sage_const_2 *self._nu2-_sage_const_2 /_sage_const_3 *self._nu3)
        # Recall that the true permutations are given by the coset reps.
        #print "GOT_CUSPS"
        # These are dictionaries with vertices/cusps as keys
        self._vertex_map= None    # Maps each vertex to the corr. cusp 
        self._cusp_normalizer=None # Cusp normalizers (unscaled)
        self._cusp_stabilizer=None # Generators of the stabilizers
        self._cusp_width=None      # Cusp widths
        if(self._generalised_level==None):
            self._generalised_level=self._G.generalised_level()
        # for some reason this doesn't work for SL2Z
        if(len(self.permS)==1 and len(self.permR)==1):
            self._is_congruence=True
        else:
            self._is_congruence=self._G.is_congruence()
        if(self.is_congruence()):
            self._level=self._generalised_level

        #print "1"
        [self._cusp_representative,self._vertex_map]=self._get_vertex_maps()
        #print "2 GOT_vertex_DATA"
        [self._cusp_normalizer,self._cusp_stabilizer,self._cusp_width]=self._get_cusp_data()
        #print "3 GOT_CUSP_DATA"
        # We need some manageble unique (string) identifier
        # Which also contains info in clear text...
        self._uid = self._get_uid()
 
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
        s ="Arithmetic Subgroup of PSL2(Z) with index "+str(self._index)+". "
        s+="Given by: \n \t perm(S)="+str(self.permS)+"\n \t perm(ST)="+str(self.permR)
        s+="\nConstructed from G="+self._G._repr_()

        #self._name=s
        #print "s=",s
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
        domore=False
        s=str(self.permS._list)+"-"
        s=s+str(self.permR._list)
        s=s.replace(", ","_")
        if(self._is_congruence):
            # test:
            if(self._G == Gamma0(N) or self._G == Gamma0(N).as_permutation_group()):
                s+="-Gamma0("+str(N)+")"
            else:
                if(self._G == Gamma(N)):
                    s+="-Gamma("+str(N)+")"
                elif(Gamma(N).is_even()):
                    if(self._G == Gamma(N).as_permutation_group()):
                        s+="-Gamma("+str(N)+")"
                if(self._G == Gamma1(N)):
                    s+="-Gamma1("+str(N)+")"
                elif(Gamma1(N).is_even()):
                    if(self._G == Gamma1(N) or self._G == Gamma1(N).as_permutation_group()):
                        s+="-Gamma1("+str(N)+")"
        return s 

    def __reduce__(self):
        r"""
        Used for pickling self.
        
        EXAMPLES::


            sage: G=MySubgroup(Gamma0(5))

            Not implmented!
        """
        #raise NotImplementedError, "Not (correctly) implemented!"
        return (MySubgroup, (None,self.permS,self.permR))

    def __cmp__(self,other):
        r""" Compare self to other.

        EXAMPLES::

            sage: G=MySubgroup(Gamma0(5))
            sage: G <> Gamma0(5)
            False
            sage: GG=MySubgroup(None,G.permS,G.permR)
            sage: GG == G

        
        """
        return self._G.__cmp__(other._G)

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
        try:
            pR=G.permR
            pS=G.permS
            if(pR==self.permR and pS==self.permS):
                return True
        except:
            pass
        return False

    def _get_index(self,G=None,o2=None,o3=None):
        r""" Index of self.
        
        INPUT:

            - ``G`` --  subgroup
            - ``o2``--  permutation
            - ``o3``--  permutation

        OUTPUT: 

            - integer -- the index of G in SL2Z or the size of the permutation group.

        EXAMPLES::


            sage: G=MySubgroup(Gamma0(5))
            sage: G._get_index(G._G,None,None)
            6
            sage: G._get_index(Gamma0(8))
            12
            sage: pS=P([2,1,4,3,6,5])
            sage: pR=P([3,1,2,5,6,4])
            sage: G._get_index(o2=pS,o3=pR) 
            6
        """
        if(G<>None):
            try:
                if(G.is_subgroup(SL2Z)):
                    ix=G.index()
            except:
                raise TypeError,"Wrong format for input G: Need subgroup of PSLZ! Got: \n %s" %(G) 
        else:
            if(o2==None):
                raise TypeError,"Need to supply one of G,o2,o3 got: \n G=%s,\n o2=%s,\o3=%s"%(G,o2,o3)
            if(hasattr(o2,"to_cycles")):
                ix=len(o2)
            elif(hasattr(o2,"list")): # We had an element of the SymmetricGroup
                ix=len(o2.list())
            else:
                raise TypeError,"Wrong format of input: o2=%s, \n o2.parent=%s"%(o2,o2.parent())
        return ix
    
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
            sage: l=G._coset_reps
            sage: G._get_vertices(l)
            sage: G._get_vertices(l)
            [[Infiniy, 0], {0: [ 0 -1]
            [ 1 -2], Infinity: [1 0]
            [0 1]}]
            sage: P=Permutations(6)
            sage: pS=P([2,1,4,3,6,5])
            sage: pR=P([3,1,2,5,6,4])
            sage: G=MySubgroup(o2=pS,o3=pR)   
            sage: G._get_vertices(G._coset_reps)
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
        SG1=self._P.identity()
        o22=o2*o2; o33=o3*o3*o3
        if(o22 <>SG1 or o33<>SG1):
            s="Input permutations are not of order 2 and 3: \n perm(S)^2=%s \n perm(R)^3=%s " 
            raise ValueError, s % o2*o2,o3*o3*o3
        if(not are_transitive_permutations(o2,o3)):
            GS=S.subgroup([self.permS,self.permR])
            s="Error in construction of permutation group! Permutations not transitive: <S,R>=%s"
            raise ValueError, s % GS.list()
        return True


    def permutation_action(self,A):
        r""" The permutation corresponding to the matrix A
        INPUT:
         - ''A'' Matrix in SL2Z
        OUPUT:
         element in self._P given by the presentation of the group self


        EXAMPLES::


            P=Permutations(6)
            sage: pS=P([2,1,4,3,6,5])
            sage: pR=P([3,1,2,5,6,4])
            sage: G=MySubgroup(o2=pS,o3=pR)
            sage: G.permutation_action(S*T).to_cycles()
            [(1, 3, 2), (4, 5, 6)]
            sage: G.permutation_action(S).to_cycles()
            [(1, 2), (3, 4), (5, 6)]
            sage: pR.to_cycles()
            [(1, 3, 2), (4, 5, 6)]
            sage: pS.to_cycles()
            [(1, 2), (3, 4), (5, 6)]
            sage: G=MySubgroup(Gamma0(5))
            sage: G.permutation_action(S).to_cycles()
            [(1, 2), (3, 4), (5,), (6,)]
            sage: G.permutation_action(S*T).to_cycles()
            [(1, 3, 2), (4, 5, 6)]

        """
        [cf,sg]=factor_matrix_in_sl2z_in_S_and_T(A)
        # The sign doesn't matter since we are working in paractice only with PSL(2,Z)
        # We now have A=T^a0*S*T^a1*...*S*T^an
        # Hence
        # h(A)=h(T^an)*h(S)*...*h(T^a1)*h(S)*h(T^a0)
        #print "cf=",cf
        n=len(cf)
        #print "n=",n    
        def ppow(perm,k):
            if(k>=0 ):
                pp=perm
                for j in range(k-1 ):
                    pp=pp*perm
                return pp
            else:
                permi=perm.inverse()
                pp=permi
                for j in range(abs(k)-1 ):
                    pp=pp*permi
                return pp
        if(cf[0 ]==0 ):
            p=self._P.identity()
        else:
            p=ppow(self.permT,cf[0 ])
        #print "p(",0,")=",p.to_cycles()
        for j in range(1 ,n):
            a=cf[j]
            if(a<>0 ):
                p=p*self.permS*ppow(self.permT,a)
            else:
                p=p*self.permS
            #print "p(",j,")=",p.to_cycles()
        return p


    
    def _get_coset_reps_from_G(self,G):
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
        cl=list()
        S,T=SL2Z.gens()
        lvl=G.generalised_level()
        # Start with identity rep.
        cl.append(SL2Z([1 ,0 ,0 ,1 ]))
        if(not S in G):
            cl.append(S)
        # If the original group is given as a Gamma0 then
        # the reps are not the one we want
        # I.e. we like to have a fundamental domain in
        # -1/2 <=x <= 1/2 for Gamma0, Gamma1, Gamma
        for j in range(1 ,Integer(lvl)/_sage_const_2 +_sage_const_2 ):
            for ep in [1 ,-1 ]:
                if(len(cl)>=self._index):
                    break
                # The ones about 0 are all of this form
                A=SL2Z([0 ,-1 ,1 ,ep*j])
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

    
        

    def _get_coset_reps_from_perms(self,pS,pR,pT):
        r"""
        Compute a better/nicer list of right coset representatives
        i.e. SL2Z = \cup G V_j
        Todo: Check consistency of input

        EXAMPLES::

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
        ix=len(pR)
        S,T=SL2Z.gens()
        R=S*T
        Id=SL2Z([1 ,0 ,0 ,1 ])
        cl=list(range(ix))
        pR2=pR*pR
        for i in range(ix):
            cl[i]=Id
        deb=False# True #False
        if(deb):
            print "T=",T,pT.to_cycles()
            print "S=",S,pS.to_cycles()
            print "R=",R,pR.to_cycles()
        cycT=pT.to_cycles()
        for cy in cycT:
            i=pT(cy[0 ])
            if(i<>cy[0 ]):
                if(cl[i-1 ]==Id and i-1 <>0 ):
                    cl[i-1 ]=cl[cy[0 ]-1 ]*T
                    if(deb):
                        print "VT(",i-1 ,")=cl[",cy[0 ]-1 ,"]*T"
                        print "VT(",i-1 ,")=",cl[cy[0 ]-1 ],"*",T
                j=pT(i)
                while(j<>cy[0 ]):
                    if(cl[j-1 ]==Id and j-1 <>0 ):
                        cl[j-1 ]=cl[i-1 ]*T
                        if(deb):
                            print "VT(",j-1 ,")=cl[",i-1 ,"]*T"
                            print "VT(",j-1 ,")=",cl[i-1 ],"*",T
                    i=j
                    j=pT(j)
            #for A in cl:
            #    print A
            if(pS(i)<>i):
                if(cl[pS(i)-1 ]==Id and pS(i)-1 <>0 ):
                    cl[pS(i)-1 ]=cl[i-1 ]*S
                    if(deb):
                        print "VS(",pS(i)-1 ,")=",cl[i-1 ]*S
                        print "cl[",pS(i)-1 ,"]=cl[",i-1 ,"]*S"
                        print "=",cl[pS(i)-1 ]
            if(pR(i)<>i):
                i1=pR(i)-1 
                i2=pR(pR(i))-1 
                #print "i1=",i1
                #print "i2=",i2
                if(cl[i1]==Id and i1<>0 ):
                    cl[pR(i)-1 ]=cl[i-1 ]*R
                    if(deb):
                        print "VR(",i1,")=",cl[i-1 ]*R
                if(cl[i2]==Id and i2<>0 ):
                    cl[i2]=cl[i-1 ]*R*R
                    if(deb):
                        print "VRR(",i2,")=",cl[i-1 ],"R^2=",cl[i-1 ]*R*R
                        print "p(V(i2))=",self.permutation_action(cl[i2])
                        print "p(V(i-1))*p(R^2)=",self.permutation_action(cl[i-1 ])*pR2
                # If we missed something (which is unlikely)        
        if(cl.count(Id)>1  or len(cl)<>ix):
            #print "cl=",cl
            raise ValueError,"Problem getting coset reps! Need %s and got %s" %(self._index,len(cl))
        ## Test the reps
        #print "cl=",cl
        try:
            for A in cl:
                for B in cl:
                    if(A<>B):
                        C=A*(B**-1 )
                        p=self.permutation_action(C)
                        if(p(1 )==1 ):
                            #print "A=",A
                            #print "B=",B
                            #print "AB^-1=",C
                            #print "p(AB^-1)=",p
                            pA=self.permutation_action(A)
                            pB=self.permutation_action(B)
                            #print " pA=",pA
                            #print "pB=", pB
                            #print "~pB=", pB.inverse()
                            pBi=self.permutation_action(B**-1 )
                            #print " B^-1=",pBi
                            raise StopIteration()
        except StopIteration:
            raise ArithmeticError," Could not get coset reps! \n %s * %s ^-1 in G : p(AB^-1)=%s" % (A,B,p)
        return cl


    def are_equivalent_cusps(self,p,c):
        r"""
        Check whether two cusps are equivalent with respect to self

        EXAMPLES::


            sage: G=MySubgroup(Gamma0(5))
            sage: G.are_equivalent_cusps(Cusp(1),Cusp(infinity))
            False
            sage: G.are_equivalent_cusps(Cusp(1),Cusp(0))
            True
            
        """
        if(p==c):
            return True
        if(p.denominator()<>0  and c.denominator()<>0 ):
            # The builtin equivalent cusp function is
            # slow for the cusp at infinty infinity
            return self._G.are_equivalent_cusps(p,c)
        elif(p.denominator()<>0  and c.denominator()==0 ):
            w = lift_to_sl2z(p.denominator(), p.numerator(), 0 )
            n = SL2Z([w[3 ], w[1 ], w[2 ],w[0 ]])
            #print "n=",n
            for i in range(len(self.permT.to_cycles()[0 ])):
                test=n*SL2Z([1 ,i,0 ,1 ])
                #print "test=",test
                if test in self:
                    return True
            return False
        elif(p.denominator()==0  and c.denominator()<>0 ):
            w = lift_to_sl2z(c.denominator(), c.numerator(), 0 )
            n = SL2Z([w[3 ], w[1 ], w[2 ],w[0 ]])
            for i in range(len(self.permT.to_cycles()[0 ])):
                if n*SL2Z([1 ,i,0 ,1 ]) in self:
                    return True
            return False
        else:
            return True
        
            
    def _get_cusps(self,l):
        r""" Compute a list of inequivalent cusps from the list of vertices

        EXAMPLES::


            sage: P=Permutations(6)
            sage: pS=P([2,1,4,3,6,5])
            sage: pR=P([3,1,2,5,6,4])
            sage: G=MySubgroup(o2=pS,o3=pR)   
            sage: l=G._get_vertices(G._coset_reps)[0];l
            [Infinity, 0, -1/2]
            sage: G._get_cusps(l)
        
        """
        lc=list()
        for p in l:
            are_eq=False
            for c in lc:
                if(self.are_equivalent_cusps(p,c)):
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

    def _get_vertex_maps(self):
        r"""
        OUTPUT:
         -- [vc,lu]
            vc = dictionary of vertices and corresponding cusps
            lu = dictionary of the maps mapping the vertices to the
                 corresponding  cusp representatives

         EXAMPLES::


             sage: G=MySubgroup(Gamma0(5))
             G._get_vertex_maps()
             [{0: 0, Infinity: Infinity}, {0: [1 0]
             [0 1], Infinity: [1 0]
             [0 1]}]
             sage: P=Permutatons(6)
             sage: pS=[2,1,4,3,6,5]
             sage: pR=[3,1,2,5,6,4]
             sage: G=MySubgroup(o2=pS,o3=pR)            
             sage: G._get_vertex_maps()
             [{0: 0, Infinity: Infinity, -1/2: -1/2}, {0: [1 0]
             [0 1], Infinity: [1 0]
             [0 1], -1/2: [1 0]
             [0 1]}]

        """
        lu=dict()
        vc=dict()
        Id=SL2Z([1 ,0 ,0 ,1 ])
        for p in self._vertices:
            if(self._cusps.count(p)>0 ):
                lu[p]=Id; vc[p]=p
                continue
            Vp=self._vertex_reps[p]
            try:
                for c in self._cusps:
                    Vc=self._vertex_reps[c]**-1 
                    for j in range(1 ,self._index):
                        Tj=SL2Z([1 ,j,0 ,1 ])
                        g=Vp*Tj*Vc
                        if(g in self):
                            lu[p]=g**-1 
                            vc[p]=c
                            raise StopIteration()
            except StopIteration:
                pass
        # We reorder the cusps so that the first in the list is Infinity and the second is 0 (if they exist)
        # Recall that there is always a cusp at infinity here
        return [vc,lu]
    
    def _get_cusp_data(self):
        r""" Return dictionaries of cusp normalizers, stabilisers and widths

        OUTPUT:
          -- [ns,ss,ws]  - list of dictionaries with cusps p as keys.
             ns : ns[p] = A in SL2Z with A(p)=infinity (cusp normalizer)
             ss : ss[p] = A in SL2Z s.t. A(p)=p   (cusp stabilizer)
             ws : ws[p]= (w,s).
                 w = width of p and s = 1 if p is regular and -1 if not

        EXAMPLES::
        

            sage: P=Permutations(6)
            sage: pS=P([2,1,4,3,6,5])
            sage: pR=P([3,1,2,5,6,4])
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
        try:
            p=self._cusps[0 ]
        except:
            self._get_cusps(self._vertices)
        ns=dict()
        ss=dict()
        ws=dict()
        lws=self.permT.to_cycles()
        cusp_widths=map(len,lws)
        cusp_widths.sort()
        # Recall that we have permutations for all groups here
        # and that the widths can be read from self.permT
        for p in self._cusps:
            # could use self._G.cusp_data but is more efficient to redo all
            if(p.denominator()==0 ):
                wi=len(lws[0 ])
                ns[p]= SL2Z([1 , 0 , 0 ,1 ])
                ss[p]= SL2Z([1 , wi,0 ,1 ])
                ws[p]=(wi,1 )  
                cusp_widths.remove(wi)
            else:
                w = lift_to_sl2z(p.denominator(), p.numerator(), 0 )
                n = SL2Z([w[3 ], w[1 ], w[2 ],w[0 ]])
                stab=None
                wi=0 
                try:
                    for d in cusp_widths:
                        if n * SL2Z([1 ,d,0 ,1 ]) * (~n) in self:
                            stab = n * SL2Z([1 ,d,0 ,1 ]) * (~n)
                            wi = (d,1 )
                            cusp_widths.remove(d)
                            raise StopIteration()
                        elif n * SL2Z([-1 ,-d,0 ,-1 ]) * (~n) in self:
                            stab = n * SL2Z([-1 ,-d,0 ,-1 ]) * (~n)
                            wi=(d,-1 )
                            cusp_widths.remove(d)
                            raise StopIteration()
                except StopIteration:
                    pass
                if(wi==0 ):
                    raise ArithmeticError, "Can not find cusp stabilizer for cusp:" %p
                ss[p]=stab
                ws[p]=wi
                ns[p]=n
        return [ns,ss,ws]

    def coset_reps(self):
        r""" Returns coset reps of self


        EXAMPLES::

        
            sage: P=Permutations(6)
            sage: pS=P([2,1,4,3,6,5]); pS
            (1,2)(3,4)(5,6)
            sage: pR=P([3,1,2,5,6,4]); pR
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

        return self._coset_reps

    def _get_perms_from_coset_reps(self):
        r""" Get permutations of order 2 and 3 from the coset representatives

        EXAMPLES::
        

            sage: P=Permutations(6)
            sage: pS=P([2,1,4,3,6,5]); pS
            (1,2)(3,4)(5,6)
            sage: pR=P([3,1,2,5,6,4]); pR
            (1,3,2)(4,5,6)
            sage: G=MySubgroup(o2=pS,o3=pR)
            sage: G._get_perms_from_coset_reps()
            [(1,2)(3,4)(5,6), (1,3,2)(4,5,6)]
            sage: G=MySubgroup(Gamma0(6))
            sage: p=G._get_perms_from_coset_reps(); p[0]; p[1]
            (1,2)(3,4)(5,8)(6,9)(7,10)(11,12)
            (1,3,2)(4,5,9)(6,11,10)(7,12,8)


        """
        l=self._coset_reps
        li=list()
        n=len(l)
        G=self._G
        ps=range(1 ,self._index+1 )
        pr=range(1 ,self._index+1 )
        S,T=SL2Z.gens()
        R=S*T
        for i in range(n):
            li.append(l[i]**-1 )
        #used_inS=dict()
        #used_inR=dict()
        #for i in range(len(l)):
        ##    used_inS[i]=False
        #    used_inR[i]=False
        ixr=range(n)
        ixs=range(n)
        for i in range(n):
            [a,b,c,d]=l[i]
            VS=SL2Z([b,-a,d,-c]) # Vi*S
            VR=SL2Z([b,b-a,d,d-c]) # Vi*R=Vi*S*T
            for j in ixr:
                Vji=li[j] # Vj^-1
                if(VR*Vji in G):
                    pr[i]=j+1 
                    ixr.remove(j)
                    break
            for j in ixs:
                Vji=li[j]
                if(VS*Vji in G):
                    ps[i]=j+1 
                    ixs.remove(j)
                    break
        return [self._P(ps),self._P(pr)]

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
        for V in (self._coset_reps):
            if( (A*V**-1 ) in self):
                return V
        raise ArithmeticError,"Did not find coset rep. for A=%s" %(A)

    def pullback(self,x_in,y_in,prec=201 ):
        r""" Find the pullback of a point in H to the fundamental domain of self
        INPUT:

         - ''x_in,y_in'' -- x_in+I*y_in is in the upper half-plane
         - ''prec''      -- (default 201) precision in bits

        OUTPUT:
        
         - [xpb,ypb,B]  --  xpb+I*ypb=B(x_in+I*y_in) with B in self
                           xpb and ypb are complex numbers with precision prec 
        EXAMPLES::


            sage: P=Permutatons(6)
            sage: pS=[2,1,4,3,6,5]
            sage: pR=[3,1,2,5,6,4]
            sage: G=MySubgroup(o2=pS,o3=pR)
            sage: [x,y,B]=G.pullback(0.2,0.5,53); x,y;B
            (-0.237623762376238, 0.123762376237624)
            [-1  0]
            [ 4 -1]
            sage: (B**-1).acton(x+I*y)
            0.200000000000000 + 0.500000000000000*I


        """
        x=deepcopy(x_in); y=deepcopy(y_in)
        A=pullback_to_psl2z_mat(RR(x),RR(y))
        A=SL2Z(A)
        try:
            for V in self._coset_reps:
                B=V*A
                if(B in self):
                    raise StopIteration
        except StopIteration:            
            pass
        else:
            raise ArithmeticError,"Did not find coset rep. for A=%s" % A
        #B=V*A
        z=CC.to_prec(prec)(x_in)+I*CC.to_prec(prec)(y_in)
        [a,b,c,d]=B
        zpb=CC.to_prec(prec)(a*z + b)/CC.to_prec(prec)(c*z + d)
        xpb=zpb.real();ypb=zpb.imag()
        #print "A=",A        
        return [xpb,ypb,B]

    def is_congruence(self):
        r""" Is self a congruence subgroup or not?

        EXAMPLES::


        sage: P=Permutations(6)
        sage: pR=P([3,1,2,5,6,4]); pR
        (1,3,2)(4,5,6)
        sage: pS=P([2,1,4,3,6,5]); pS
        (1,2)(3,4)(5,6)
        sage: G=MySubgroup(o2=pS,o3=pR)
        sage: G.is_congruence()
        True
        sage: P=Permutations(7)
        sage: pS=P([1,3,2,5,4,7,6]); pS
        (2,3)(4,5)(6,7)
        sage: pR=P([3,2,4,1,6,7,5]); pR
        (1,3,4)(5,6,7)
        sage: G=MySubgroup(o2=pS,o3=pR)
        sage: G.is_congruence()
        False
        
        """
        try:
            return self._is_congruence
        except:
            self._is_congruence=self._G.is_congruence()
            return self._is_congruence

    def generalised_level(self):
        r""" Generalized level of self

        EXAMPLES::


        sage: P=Permutations(6)
        sage: pR=P([3,1,2,5,6,4]); pR
        (1,3,2)(4,5,6)
        sage: pS=P([2,1,4,3,6,5]); pS
        (1,2)(3,4)(5,6)
        sage: G=MySubgroup(o2=pS,o3=pR)
        sage: G.generalised_level()
        4
        sage: P=Permutations(7)
        sage: pS=P([1,3,2,5,4,7,6]); pS
        (2,3)(4,5)(6,7)
        sage: pR=P([3,2,4,1,6,7,5]); pR
        (1,3,4)(5,6,7)
        sage: G=MySubgroup(o2=pS,o3=pR)
        sage: G.generalised_level()
        6
        """
        try:
            if(self._generalised_level<>None):
                return self._generalised_level
            else:
                raise ArithmeticError, "Could not compute generalised level of %s" %(self)
        except:
            self._generalised_level=self._G.generalised_level()
            return self._generalised_level


    def level(self):
        r""" Level of self

        EXAMPLES::


            sage: G=MySubgroup(Gamma0(5));
            sage: G.level()
            5

        """
        if(self._is_congruence):
            return self._level
        else:
            raise TypeError,"Group is not a congruence group! Use G.generalised_level() instead!"
            

    def gens(self):
        r""" Generators of self


        EXAMPLES::


            sage: G=MySubgroup(Gamma0(5));
            sage: G.gens()
            ([1 1]
            [0 1], [-1  0]
            [ 0 -1], [ 1 -1]
            [ 0  1], [1 0]
            [5 1], [1 1]
            [0 1], [-2 -1]
            [ 5  2], [-3 -1]
            [10  3], [-1  0]
            [ 5 -1], [ 1  0]
            [-5  1])
            sage: G3=MySubgroup(o2=Permutations(3)([1,2,3]),o3=Permutations(3)([2,3,1]))
            sage: G3.gens()
            ([1 3]
            [0 1], [ 0 -1]
            [ 1  0], [ 1 -2]
            [ 1 -1], [ 2 -5]
            [ 1 -2])

        """
        return self._G.gens()


    def closest_vertex(self,x,y):
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
        ymax=0 
        vmax=None
        for v in self._vertices:
            c=self._cusp_representative[v]
            U=self._vertex_map[v]
            N=self._cusp_normalizer[c]
            w=self._cusp_width[c][0 ]
            A=(N**-1 *U)
            [a,b,c,d]=A
            den=(c*x+d)**2 +(c*y)**2
            y2=y/den/w # We only want the y coordinate
            #print "v=",v
            #print "c=",c
            #print "A=",A
            #print "y=",y
            if(y2>ymax):
                ymax=y2
                vmax=v
        if(vmax==None):
            raise ArithmeticError," Did not find closest vertex to z=%s+i%s," %(x,y)
        return vmax



    
    
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
        Returns the dimension of the space of cuspforms on G of weight k
        where k is an even integer

        EXAMPLES::


            sage: G=MySubgroup(Gamma0(4))
            sage: G.dimension_cuspforms(12)
            4
            sage: P=Permutations(7)
            sage: pR=P([3,2,4,1,6,7,5]); pR
            (1,3,4)(5,6,7)
            sage: pS=P([1,3,2,5,4,7,6]); pS
            (2,3)(4,5)(6,7)
            sage: G.dimension_cuspforms(4)
            1

        
        """
        kk=Integer(k)
        if(is_odd(kk)):
            raise ValueError, "Use only for even weight k! not k=" %(kk)
        if(kk<2 ):
            dim=0 
        elif(kk==2 ):
            dim=self._genus
        elif(kk>=_4 ):
            dim=Integer(kk-1 )*(self._genus-1 )+self._nu2*int(floor(kk/_sage_const_4 ))+self._nu3*int(floor(kk/_sage_const_3 ))+(kk/_sage_const_2 - _sage_const_1)*self._ncusps
        return dim

    def dimension_modularforms(self,k):
        r"""
        Returns the dimension of the space of modular forms on G of weight k
        where k is an even integer

        EXAMPLES::


            sage: P=Permutations(7)
            sage: pR=P([3,2,4,1,6,7,5]); pR
            (1,3,4)(5,6,7)
            sage: pS=P([1,3,2,5,4,7,6]); pS
            (2,3)(4,5)(6,7)
            sage: G.dimension_modularforms(4)
            3
        """
        kk=Integer(k)
        if(is_odd(kk)):
            raise ValueError, "Use only for even weight k! not k=" %(kk)
        if(k==0 ):
            dim=1 
        elif(k<2 ):
            dim=0 
        else:
            dim=(kk-_sage_const_1 )*(self._genus-_sage_const_1 )+self._nu2()*int(floor(kk/_sage_const_4 ))+self._nu3*int(floor(kk/_sage_const_3 ))+kk/_sage_const_2 *self._ncusps()
        return dim
    
    ### Overloaded operators
    def __contains__(self,A):
        r"""
        Is A an element of self (this is an ineffective implementation if self._G is a permutation group)

        EXAMPLES::
        

            sage: G=MySubgroup(Gamma0(5))
            sage: A=SL2Z([-69,-25,-80,-29])
            sage: G.__contains__(A)
            True
            sage: A in G
            True            
            sage: P=Permutations(7)
            sage: pR=P([3,2,4,1,6,7,5]); pR
            (1,3,4)(5,6,7)
            sage: pS=P([1,3,2,5,4,7,6]); pS
            (2,3)(4,5)(6,7)
            sage: S,T=SL2Z.gens(); R=S*T
            sage: A=S*T^4*S
            sage: A in G
            False
            sage: A=S*T^6*S
            sage: A in G
            True
            
        """
        p=self.permutation_action(A)
        if(p(1 )==1 ):
            return True
        else:
            return False

    def cusps(self):
        r"""
        Returns the cusps of self

        EXAMPLES::


            sage: G=MySubgroup(Gamma0(5))
            sage: G.cusps()
            [Infinity, 0]
        
        """
        return self._cusps

    def cusp_width(self,cusp):
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
        try:
            return self._cusp_width[cusp][0 ]
        except KeyError:
            (A,d,e)=self.cusp_data(cusp)
            return d

    def cusp_data(self,c):
        r""":
        Returns cuspdata in the same format as for the generic Arithmetic subgroup

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
        try:
            # if we give a cusp already in the list we return stored values
            w=self.cusp_normalizer(c)
            d=self._cusp_width[c][0 ]
            e=self._cusp_width[c][1 ]
            g=w * SL2Z([e,d*e,0 ,e]) * (~w)
            return (self.cusp_normalizer(c),d*e,1)
        except:
            w = lift_to_sl2z(c.denominator(), c.numerator(), 0 )
            g = SL2Z([w[3 ], w[1 ], w[2 ],w[0 ]])
            for d in xrange(1 ,1 +self.index()):
                if g * SL2Z([1 ,d,0 ,1 ]) * (~g) in self:
                    return (g * SL2Z([1 ,d,0 ,1 ]) * (~g), d, 1 )
                elif g * SL2Z([-1 ,-d,0 ,-1 ]) * (~g) in self:
                    return (g * SL2Z([-1 ,-d,0 ,-1 ]) * (~g), d, -1 )
            raise ArithmeticError, "Can' t get here!"
        #

    def cusp_normalizer(self,cusp):
        r"""
        Return the cuspnormalizer of cusp

        EXAMPLES::


            sage: P=Permutations(6)
            sage: pR=P([3,1,2,5,6,4]); pR
            (1,3,2)(4,5,6)
            sage: pS=P([2,1,4,3,6,5]); pS
            (1,2)(3,4)(5,6)
            sage: G=MySubgroup(o2=pS,o3=pR)
            sage: G.cusps()
            [Infinity, 0, -1/2]
            sage: G.cusp_normalizer(Cusp(0))
            [ 0 -1]
            [ 1  0]
            sage: G.cusp_normalizer(Cusp(-1/2))
            [-1  0]
            [ 2 -1]
            sage: G.cusp_normalizer(Cusp(Infinity))
            [1 0]
            [0 1]

        
        """

        try:
            return self._cusp_normalizer[cusp]
        except KeyError:
            try:
                w = lift_to_sl2z(c.denominator(), c.numerator(), 0 )
                g = SL2Z([w[3 ], w[1 ], w[2 ],w[0 ]])
                self._cusp_normalizer[cusp]=g
                return g
            except:
                raise ValueError,"Supply a cusp as input! Got:%s" %cusp
        
    def draw_fundamental_domain(self,model="H",axes=None,filename=None,**kwds):
        r""" Draw fundamental domain
        INPUT:
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
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        if(model=="D"):
            g=self._draw_funddom_d(format,I)
        else:
            g=self._draw_funddom(format)
        if(axes<>None):
            [x0,x1,y0,y1]=axes
        elif(model=="D"):
            x0=-1 ; x1=1 ; y0=-1 ; y1=1 
        else:
            # find the width of the fundamental domain
            w=0  #self.cusp_width(Cusp(Infinity))
            wmin=0 ; wmax=1 
            for V in self._coset_reps:
                if(V[1 ,0 ]==0  and V[0 ,0 ]==1 ):
                    if(V[0 ,1 ]>wmax):
                        wmax=V[0 ,1 ]
                    if(V[0 ,1 ]<wmin):
                        wmin=V[0 ,1 ]
            #print "wmin,wmax=",wmin,wmax
            #x0=-1; x1=1; y0=-0.2; y1=1.5
            x0=wmin-1 ; x1=wmax+1 ; y0=-_sage_const_0p2 ; y1=_sage_const_1p5 
        g.set_aspect_ratio(1 )
        g.set_axes_range(x0,x1,y0,y1)
        if(filename<>None):
            fig = g.matplotlib()
            fig.set_canvas(FigureCanvasAgg(fig))
            axes = fig.get_axes()[0 ]
            axes.minorticks_off()
            axes.set_yticks([])
            fig.savefig(filename,**kwds)
        else:
            return g
#        g.show(figsize=[5,5])

    def _draw_funddom(self,format="S"):
        r""" Draw a fundamental domain for G.

        INPUT:

         - ``format``  -- (default 'Disp') How to present the f.d.
         -   ``S`` -- Display directly on the screen
               
        EXAMPLES::        


            sage: G=MySubgroup(Gamma0(3))
            sage: G._draw_funddom()
        
        """
        pi=RR.pi()
        from sage.plot.plot import (Graphics,line)
        from sage.functions.trig import (cos,sin)
        g=Graphics()
        x1=-_sage_const_0p5 ; y1=sqrt(_sage_const_3 )/_sage_const_2 
        x2=_sage_const_0p5 ; y2=sqrt(_sage_const_3 )/_sage_const_2 
        xmax=_sage_const_20 
        l1 = line([[x1,y1],[x1,xmax]])
        l2 = line([[x2,y2],[x2,xmax]])
        l3 = line([[x2,xmax],[x1,xmax]]) # This is added to make a closed contour
        c0=_circ_arc(pi/_sage_const_3p0 ,_sage_const_2p0 *pi/_sage_const_3p0 ,0 ,1 ,_sage_const_100 )
        tri=c0+l1+l3+l2
        g+=tri
        for A in self._coset_reps:
            [a,b,c,d]=A
            if(a==1  and b==0  and c==0  and d==1 ):
                continue
            if(a<0 ):
                a=-a; b=-b; c=-c; d=-1 
            if(c==0 ): # then this is easier
                L0 = [[cos(pi/_sage_const_3 *i/_sage_const_100 )+b,sin(pi/_sage_const_3 *i/_sage_const_100 )] for i in range(_sage_const_100 ,_sage_const_201 )]
                L1 = [[x1+b,y1],[x1+b,xmax]]
                L2 = [[x2+b,y2],[x2+b,xmax]]
                L3 = [[x2+b,xmax],[x1+b,xmax]]
                c0=line(L0); l1=line(L1); l2=line(L2); l3=line(L3)
                tri=c0+l1+l3+l2
                g+=tri
            else:
                den=(c*x1+d)**_sage_const_2 +c**_sage_const_2 *y1**_sage_const_2 
                x1_t=(a*c*(x1**_sage_const_2 +y1**_sage_const_2 )+(a*d+b*c)*x1+b*d)/den
                y1_t=y1/den
                den=(c*x2+d)**_sage_const_2 +c**_sage_const_2 *y2**_sage_const_2 
                x2_t=(a*c*(x2**_sage_const_2 +y2**_sage_const_2 )+(a*d+b*c)*x2+b*d)/den
                y2_t=y2/den
                inf_t=a/c
                #print "A=",A
                #print "arg1=",x1_t,y1_t,x2_t,y2_t
                c0=_geodesic_between_two_points(x1_t,y1_t,x2_t,y2_t)
                #print "arg1=",x1_t,y1_t,inf_t
                c1=_geodesic_between_two_points(x1_t,y1_t,inf_t,_sage_const_0p0 )
                #print "arg1=",x2_t,y2_t,inf_t
                c2=_geodesic_between_two_points(x2_t,y2_t,inf_t,_sage_const_0p0 )
                tri=c0+c1+c2
                g+=tri
        return g


    def _draw_funddom_d(self,format="MP",z0=I):
        r""" Draw a fundamental domain for self in the circle model
        INPUT:
         - ''format''  -- (default 'Disp') How to present the f.d.
               =  'S'  -- Display directly on the screen
         - z0          -- (default I) the upper-half plane is mapped to the disk by z-->(z-z0)/(z-z0.conjugate())
        EXAMPLES::
        

            sage: G=MySubgroup(Gamma0(3))
            sage: G._draw_funddom_d()
        
        """
        # The fundamental domain consists of copies of the standard fundamental domain
        pi=RR.pi()
        from sage.plot.plot import (Graphics,line)
        g=Graphics()
        bdcirc=_circ_arc(0 ,_sage_const_2 *pi,0 ,1 ,1000 )
        g+=bdcirc
        # Corners
        x1=-_sage_const_0p5 ; y1=sqrt(_sage_const_3 )/_sage_const_2 
        x2=_sage_const_0p5 ; y2=sqrt(_sage_const_3 )/_sage_const_2 
        z_inf=1 
        l1 = _geodesic_between_two_points_d(x1,y1,x1,infinity)
        l2 = _geodesic_between_two_points_d(x2,y2,x2,infinity)
        c0 = _geodesic_between_two_points_d(x1,y1,x2,y2)
        tri=c0+l1+l2
        g+=tri
        for A in self._coset_reps:
            [a,b,c,d]=A
            if(a==1  and b==0  and c==0  and d==1 ):
                continue
            if(a<0 ):
                a=-a; b=-b; c=-c; d=-1 
            if(c==0 ): # then this is easier
                l1 = _geodesic_between_two_points_d(x1+b,y1,x1+b,infinity)
                l2 = _geodesic_between_two_points_d(x2+b,y2,x2+b,infinity)
                c0 = _geodesic_between_two_points_d(x1+b,y1,x2+b,y2)
                # c0=line(L0); l1=line(L1); l2=line(L2); l3=line(L3)
                tri=c0+l1+l2
                g+=tri
            else:
                den=(c*x1+d)**_sage_const_2 +c**_sage_const_2 *y1**_sage_const_2 
                x1_t=(a*c*(x1**_sage_const_2 +y1**_sage_const_2 )+(a*d+b*c)*x1+b*d)/den
                y1_t=y1/den
                den=(c*x2+d)**_sage_const_2 +c**_sage_const_2 *y2**_sage_const_2 
                x2_t=(a*c*(x2**_sage_const_2 +y2**_sage_const_2 )+(a*d+b*c)*x2+b*d)/den
                y2_t=y2/den
                inf_t=a/c
                c0=_geodesic_between_two_points_d(x1_t,y1_t,x2_t,y2_t)
                c1=_geodesic_between_two_points_d(x1_t,y1_t,inf_t,_sage_const_0p0 )
                c2=_geodesic_between_two_points_d(x2_t,y2_t,inf_t,_sage_const_0p0 )
                tri=c0+c1+c2
                g+=tri
        g.xmax(1 )
        g.ymax(1 )
        g.xmin(-1 )
        g.ymin(-1 )
        g.set_aspect_ratio(1 )
        return g

    def minimal_height(self):
        r""" Computes the minimal (invariant) height of the fundamental region of self.

        EXAMPLES::

        
            sage: G=MySubgroup(Gamma0(6))
            sage: G.minimal_height()
            0.144337567297406
            sage: P=Permutations(6)
            sage: pR=P([3,1,2,5,6,4])
            sage: pS=P([2,1,4,3,6,5])
            sage: G=MySubgroup(o2=pS,o3=pR)
            sage: G.minimal_height()   
            0.216506350946110


        """
        if self._is_congruence:
            l=self.level()
            if(self._G == Gamma0(l)):
                return RR(sqrt(3.0))/RR(2*l)
        # For all other groups we have have to locate the largest width
        maxw=0
        for c in self.cusps():
            l=self.cusp_width(Cusp(c))
            if(l>maxw):
                maxw=l
        return RR(sqrt(3.0))/RR(2*maxw)


#### Methods not dependent explicitly on the group 
def _geodesic_between_two_points(x1,y1,x2,y2):
    r""" Geodesic path between two points hyperbolic upper half-plane

    INPUTS:
    
    - ''(x1,y1)'' -- starting point (0<y1<=infinity)
    - ''(x2,y2)'' -- ending point   (0<y2<=infinity)
    - ''z0''  -- (default I) the point in the upper corresponding
                 to the point 0 in the disc. I.e. the transform is
                 w -> (z-I)/(z+I)
    OUTPUT:

    - ''ca'' -- a polygonal approximation of a circular arc centered
    at c and radius r, starting at t0 and ending at t1

    
    EXAMPLES::


        sage: l=_geodesic_between_two_points(0.1,0.2,0.0,0.5)
    
    """
    pi=RR.pi()
    from sage.plot.plot import line
    from sage.functions.trig import arcsin
    if(x1==x2):
        # The line segment [x=x1, y0<= y <= y1]
        return line([[x1,y1],[x2,y2]])  #[0,0,x0,infinity]
    c=(y1**_sage_const_2 -y2**_sage_const_2 +x1**_sage_const_2 -x2**_sage_const_2 )/(_sage_const_2 *(x1-x2))
    r=sqrt(y1**_sage_const_2 +(x1-c)**_sage_const_2 )
    r1=y1/r; r2=y2/r
    if(abs(r1-1 )<_sage_const_1En12 ):
        r1=_sage_const_1p0 
    elif(abs(r2+_sage_const_1 )<_sage_const_1En12 ):
        r2=-_sage_const_1p0 
    if(abs(r2-_sage_const_1 )<_sage_const_1En12 ):
        r2=_sage_const_1p0 
    elif(abs(r2+1 )<_sage_const_1En12 ):
        r2=-_sage_const_1p0 
    if(x1>=c):
        t1 = arcsin(r1)
    else:
        t1 = pi-arcsin(r1)
    if(x2>=c):
        t2 = arcsin(r2)
    else:
        t2 = pi-arcsin(r2)
    tmid = (t1+t2)*_sage_const_0p5 
    a0=min(t1,t2)
    a1=max(t1,t2)
    ##print "c,r=",c,r
    #print "t1,t2=",t1,t2
    return _circ_arc(t1,t2,c,r)

def _geodesic_between_two_points_d(x1,y1,x2,y2,z0=I):
    r""" Geodesic path between two points represented in the unit disc
         by the map w = (z-I)/(z+I)
    INPUTS:
    - ''(x1,y1)'' -- starting point (0<y1<=infinity)
    - ''(x2,y2)'' -- ending point   (0<y2<=infinity)
    - ''z0''  -- (default I) the point in the upper corresponding
                 to the point 0 in the disc. I.e. the transform is
                 w -> (z-I)/(z+I)
    OUTPUT:
    - ''ca'' -- a polygonal approximation of a circular arc centered
    at c and radius r, starting at t0 and ending at t1

    
    EXAMPLES::

        sage: l=_geodesic_between_two_points_d(0.1,0.2,0.0,0.5)
    
    """
    pi=RR.pi()
    from sage.plot.plot import line
    from sage.functions.trig import (cos,sin)
    # First compute the points
    if(y1<0  or y2<0 ):
        raise ValueError,"Need points in the upper half-plane! Got y1=%s, y2=%s" %(y1,y2)
    if(y1==infinity):
        P1=CC(1 )
    else:
        P1=CC((x1+I*y1-z0)/(x1+I*y1-z0.conjugate()))
    if(y2==infinity):
        P2=CC(1 )
    else:
        P2=CC((x2+I*y2-z0)/(x2+I*y2-z0.conjugate()))
        # First find the endpoints of the completed geodesic in D
    if(x1==x2):
        a=CC((x1-z0)/(x1-z0.conjugate()))
        b=CC(1 )
    else:
        c=(y1**2 -y2**2 +x1**2 -x2**2 )/(2 *(x1-x2))
        r=sqrt(y1**2 +(x1-c)**2 )
        a=c-r
        b=c+r
        a=CC((a-z0)/(a-z0.conjugate()))
        b=CC((b-z0)/(b-z0.conjugate()))
    if( abs(a+b) < _sage_const_1En10 ): # On a diagonal
        return line([[P1.real(),P1.imag()],[P2.real(),P2.imag()]])
    th_a=a.argument()
    th_b=b.argument()
    # Compute the center of the circle in the disc model
    if( min(abs(b-1 ),abs(b+1 ))<_sage_const_1En10  and  min(abs(a-1 ),abs(a+1 ))>_sage_const_1En10 ):
        c=b+I*(1 -b*cos(th_a))/sin(th_a)
    elif( min(abs(b-1 ),abs(b+1 ))>_sage_const_1En10  and  min(abs(a-1 ),abs(a+1 ))<_sage_const_1En10 ):
        c=a+I*(1 -a*cos(th_b))/sin(th_b)
    else:
        cx=(sin(th_b)-sin(th_a))/sin(th_b-th_a)
        c=cx+I*(1 -cx*cos(th_b))/sin(th_b)
    # First find the endpoints of the completed geodesic
    r=abs(c-a)
    t1=CC(P1-c).argument()
    t2=CC(P2-c).argument()
    #print "t1,t2=",t1,t2
    return _circ_arc(t1,t2,c,r)


def _circ_arc(t0,t1,c,r,num_pts=_sage_const_100 ):
    r""" Circular arc
    INPUTS:
    - ''t0'' -- starting parameter
    - ''t1'' -- ending parameter
    - ''c''  -- center point of the circle
    - ''r''  -- radius of circle
    - ''num_pts''  -- (default 100) number of points on polygon
    OUTPUT:
    - ''ca'' -- a polygonal approximation of a circular arc centered
    at c and radius r, starting at t0 and ending at t1

    
    EXAMPLES::

        sage: ca=_circ_arc(0.1,0.2,0.0,1.0,100)
    
    """
    from sage.plot.plot import line
    from sage.functions.trig import (cos,sin)
    t00=t0; t11=t1
    ## To make sure the line is correct we reduce all arguments to the same branch,
    ## e.g. [0,2pi]
    pi=RR.pi()
    while(t00<0.0):
        t00=t00+2.0*pi
    while(t11<0):
        t11=t11+2.0*pi
    while(t00>2*pi):
        t00=t00-2.0*pi
    while(t11>2*pi):
        t11=t11-2.0*pi

    xc=CC(c).real()
    yc=CC(c).imag()
    L0 = [[r*cos(t00+i*(t11-t00)/num_pts)+xc,r*sin(t00+i*(t11-t00)/num_pts)+yc] for i in range(0 ,num_pts)]
    ca=line(L0)
    return ca


def factor_matrix_in_sl2z_in_S_and_T(A_in):
        r"""
        Factor A in SL2Z in generators S=[[0,-1],[1,0]] and T=[[1,1],[0,1]]
        INPUT:
         - ''A_in'' -- Matrix in SL2Z
        OUTPUT:
         - ''[[a0,a1,...,an],ep] with a0 in ZZ, ai in ZZ \ {0,\pm1}, ep=1,-1
           Determined by:
           A=ep*(T^a0*S*T^a1*S*...*S*T^an)
           The algorithm uses the Nearest integer Continued Fraction expansion. Note that there due to relations in SL2Z the sequence a_j is not unique
           (see example below)

        EXAMPLES::
            sage: nearest_integer_continued_fraction(0.5);cf
            [1, 2]
            sage: A=ncf_to_matrix_in_SL2Z(cf); A
            [1 1]
            [1 2]
            sage: factor_matrix_in_sl2z_in_S_and_T(A)
            [[1, 2], 1]

        An example where we do not get back the same sequence::

            sage: cf=nearest_integer_continued_fraction(pi.n(100),nmax=10);cf
            [3, -7, 16, 294, 3, 4, 5, 15, -3, -2, 2]
            sage: A=ncf_to_matrix_in_SL2Z(cf); A
            [ -411557987 -1068966896]
            [ -131002976  -340262731]
            sage: factor_matrix_in_sl2z_in_S_and_T(A)
            [[3, -7, 16, 294, 3, 4, 5, 15, -2, 2, 3], -1]
        
        """
        if A_in not in SL2Z:
            raise TypeError, "%s must be an element of SL2Z" %A_in
        S,T=SL2Z.gens()
        # If A is not member of SL(2,Z) but a plain matrix
        # we cast it to a member of SL2Z
        A = SL2Z(A_in)
        if(A.matrix()[1 ,0 ] == 0 ):
            return [[A.matrix()[0 ,1 ]],1 ]
        x = Rational(A.matrix()[0 ,0 ] / A.matrix()[1 ,0 ])
        cf = nearest_integer_continued_fraction(x)
        B  = ncf_to_matrix_in_SL2Z(cf)
        # we know that A(oo) = x = B*S (oo) and that A and BS differ at most by a translation T^j
        Tj = S**-1  * B**-1  * A 
        sgn=1 
        if(Tj.matrix()[0 ,0 ]<0 ):
            j=-Tj.matrix()[0 ,1 ]
            sgn=-1 
            #Tj=SL2Z([1,j,0,1])
        else:
            j=Tj.matrix()[0 ,1 ]
        # if Tj = Id i.e. j=0  then A=BS otherwise A=BST^j
        cf.append(j)
        # To make sure we test
        C = B*S*Tj
        try:
            for ir in range(_sage_const_2 ):
                for ik in range(_sage_const_2 ):
                    if(C[ir,ik]<>A[ir,ik] and C[ir,ik]<>-A[ir,ik]):
                        raise StopIteration
        except StopIteration:
            print "ERROR: factorization failed %s <> %s " %(C,A)
            print "Cf(A)=",cf
            raise ArithmeticError," Could not factor matrix A=%s" % A_in
        if(C.matrix()[0 ,0 ]==-A.matrix()[0 ,0 ] and C.matrix()[0 ,1 ]==-A.matrix()[0 ,1 ]):
            sgn=-1 
        return [cf,sgn]

def ncf_to_matrix_in_SL2Z(l):
    r""" Convert a nearest integer continued fraction of x to the
    matrix A where x=A(0)
    
    EXAMPLES::
    
        sage: nearest_integer_continued_fraction(0.5);cf
        [1, 2]
        sage: A=ncf_to_matrix_in_SL2Z(cf); A
        [1 1]
        [1 2]
        sage: factor_matrix_in_sl2z_in_S_and_T(A)
        [[1, 2], 1]
        sage: cf=nearest_integer_continued_fraction(pi.n(100),nmax=10);cf
        [3, -7, 16, 294, 3, 4, 5, 15, -3, -2, 2]
        sage: A=ncf_to_matrix_in_SL2Z(cf); A
        [ -411557987 -1068966896]
        [ -131002976  -340262731]
        sage: factor_matrix_in_sl2z_in_S_and_T(A)
        [[3, -7, 16, 294, 3, 4, 5, 15, -2, 2, 3], -1]
        
    """
    S,T=SL2Z.gens()
    A=SL2Z([1,l[0],0,1])
    for j in range(1,len(l)):
        A=A*S*T**l[j]
    return A

def nearest_integer_continued_fraction(x,nmax=None):
    r""" Nearest integer continued fraction of x
    where x is a rational number, and at n digits

    EXAMPLES::

        sage: nearest_integer_continued_fraction(0.5)
        [1, 2]
        nearest_integer_continued_fraction(pi.n(100),nmax=10)
        [3, -7, 16, 294, 3, 4, 5, 15, -3, -2, 2]

    """
    if(nmax == None):
        if(x in QQ):
            nmax=10000 
        else:
            nmax=100  # For non-rational numbrs  we don't want so many convergents
    jj=0 
    cf=list()
    n=nearest_integer(x)
    cf.append(n)
    if(x in QQ):
        x1=Rational(x-n)
        while jj<nmax and x1<>0  :
            n=nearest_integer(-1 /x1)
            x1=Rational(-1 /x1-n)
            cf.append(n)
            jj=jj+1 
        return cf
    else:
        try:
            RF=x.parent()
            x1=RF(x-n)
            while jj<nmax and x1<>0  :
                n=nearest_integer(RF(-1 )/x1)
                x1=RF(-1 )/x1-RF(n)
                cf.append(n)
                jj=jj+1 
            return cf
        except AttributeError:
            raise ValueError,"Could not determine type of input x=%s" %s

def nearest_integer(x):
    r""" Returns the nearest integer to x: [x]
    using the convention that 
    [1/2]=0 and [-1/2]=0


    EXAMPLES::

        sage: nearest_integer(0)
        0
        sage: nearest_integer(0.5)
        1
        sage: nearest_integer(-0.5)
        0
        sage: nearest_integer(-0.500001)
        -1
    """
    return floor(x+_sage_const_1 /_sage_const_2 )


    def _get_perms_from_str(st):
        r""" Construct permutations frm the string st of the form
             st='[a1_a2_..._an]-[b1_...bn]'


        EXAMPLES::


            sage: G=MySubgroup(Gamma0(5))
            sage: s=G._get_uid();s
            '[2_1_4_3_5_6]-[3_1_2_5_6_4]-Gamma0(5)'
            sage: _get_perms_from_str(s)
            [(1,2)(3,4), (1,3,2)(4,5,6)]

            sage: P=Permutations(7)
            sage: pS=P([1,3,2,5,4,7,6]); pS
            (2,3)(4,5)(6,7)
            sage: pR=P([3,2,4,1,6,7,5]); pR
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
        ix=s1.count("_")+_sage_const_1 
        pS=list(range(ix))
        pR=list(range(ix))
        s1=s1.strip("["); s1=s1.strip("]");
        s2=s2.strip("["); s2=s2.strip("]");
        vs1=s1.split("_")
        vs2=s2.split("_")
        P=Permutations(ix)
        for n in range(ix):
            pS[n]=int(vs1[n])
            pR[n]=int(vs2[n])
        return [P(pS),P(pR)]
