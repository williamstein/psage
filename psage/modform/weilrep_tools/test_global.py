from __future__ import print_function
from __future__ import division
from builtins import str
from builtins import range
from sage.quadratic_forms.genera.genus import GenusSymbol_global_ring, Genus_Symbol_p_adic_ring, is_GlobalGenus
from psage.modules.finite_quadratic_module import FiniteQuadraticModule
from sage.matrix.matrix_space import MatrixSpace
from sage.all import ZZ, Zmod, sys, magma, is_fundamental_discriminant, divisors, RR, log, kronecker, kronecker_character,\
    Newforms, is_odd, prime_divisors, sqrt, factor, euler_phi, QuadraticField, Integer, prod, partitions, Set,\
    DirichletGroup, CuspForms, is_even, dimension_new_cusp_forms, pi, sigma, dimension_cusp_forms
from copy import deepcopy
from .dimension import VectorValuedModularForms
#from psage.modform.weilrep import VectorValuedModularForms

def is_global(M,r,s,return_symbol=False):
    r"""
    Test if the FiniteQuadraticModule M can be represented by a Z-lattice
    of signature ``(r,s)``.

    INPUT:

        -``M`` -- FiniteQuadraticModule
        - ``r`` -- positive integer
        - ``s`` -- positive integer
        
    OUTPUT:
        - boolean
    """

    J=M.jordan_decomposition()
    symbols={}
    n=r+s
    sig=r-s
    for A in J:
        p=A[1][0]
        if p not in symbols:
            symbols[p]=list()
        sym=list(A[1][1:len(A[1])])
        if p==2:
            if len(A[1])==4:
                sym.append(0)
                sym.append(0)
            else:
                if sym[3].kronecker(2)==sym[2]:
                    det=sym[3] % 8
                else:
                    if sym[2]==-1:
                        det=3
                    else:
                        det=1
                sym = [sym[0], sym[1], det, 1, sym[3] % 8]
                #print sym
                #if sym[1]==1:
                #    if  sym[2].kronecker(2)==sym[4].kronecker(2):
                #        sym[2]=sym[4]
                #    else:
                #        return False
        #print p, sym
        symbols[p].append(sym)
    D=M.order()*(-1)**s
    for p in symbols.keys():
        prank = sum([sym[1] for sym in symbols[p]])
        v = sum([ sym[0]*sym[1] for sym in symbols[p] ])
        Dp=D//(p**v)
        if prank != n:
            eps= (Dp*Integer(prod([ sym[2] for sym in symbols[p] ]))).kronecker(p)
            if p==2:
                if eps==-1:
                    eps=3
                symbols[p].append([0,n - prank, eps,0,0])
            else:
                if eps==-1:
                    for x in Zmod(p):
                        if not x.is_square():
                            eps=x
                            break
                symbols[p].append([0,n - prank, eps])
    symbol=GenusSymbol_global_ring(MatrixSpace(ZZ,r+s,r+s).one())
    symbol._local_symbols=[Genus_Symbol_p_adic_ring(p,syms) for p,syms in symbols.items()]
    symbol._signature=(r,s)
    #print r,s, symbol._local_symbols
    isglob=is_GlobalGenus(symbol)
    if return_symbol:
        return symbol, isglob
    else:
        return isglob

list_freitag = {
    3: ['2_7^1', '2_3^-1', '2_7^+3', '2^+2.4_7^+1', '2^2.8_3^-1',
        '2_7^+1.4^2', '2^+4.4_7^+1', '2^4.8_3^-1',
        '2_7^+1.4^+4', '2_7^+1.4^-4', '2_7^+1.3^-2', '2_7^+1.3^4',
        '4_7^1', '2_1^+1.3^+1', '2_1^1.3^-1', '8_7^1'
        ],
    4: ['2^+4.3^+1', '2^+6.3^+1', '2^+2.3^+1', '3^+1',
        '3^-3', '3^+5'
        ],
    5: ['2^+4.4_5^-1', '2^2.4_5^-1', '2^6.4_5^-1',
        '2^-6.4_1^+1', '4_5^-1'
        ],
    6: ['2^-2', '2^-4', '2^-6', '2^-8',
        '5^+1'
        ],
    7: ['4_3^-1', '2_1^+1.3^-1'
        ],
    8: ['2_2^+2', '2^-4',
        '3^-1', '7^1'
        ],
    9: ['2_1^+1', '4_1^+1', '8_1^+1'
        ],
    10: ['2^+2'
        ]
}

def test_list_freitag(m=-1):
    all=True if m==-1 else False
    for n,symbols in list_freitag.items():
        if all or n==m:
            for symbol in symbols:
                M=FiniteQuadraticModule(symbol)
                V=VectorValuedModularForms(M)
                k = Integer(2+n)/Integer(2)
                glob = is_global(M,2,n)
                if glob:
                    d=V.dimension_cusp_forms(k)
                else:
                    d="?"
                print("n = {0} {1}: Dimension({2}) = {3}, {4}, |M|={5}".format( n, symbol, str(k), d, glob, M.order()))

def search_global_symbols(n,D):
    rank=2+n
    sign=2-n
    csymbols=list() # a list of canonical symbols to avoid duplicates
    symbols=list()
    #print D
    D=(-1)**n*D
    fac = Integer(D).factor()
    symbols=list()
    for p, v in fac:
        psymbols=list()
        parts=partitions(v)
        Dp=D//(p**v)
        for vs in parts:
            #print "partition:", vs
            l=list() # list of p-symbols corresponding to the partition vs
            if len(vs) <= rank:
                exponents=Set(list(vs))
                # now we set up a list ll for each vv in the partition vs
                # that contains an entry for each possibility
                # and then update l with ll (see below)
                if p==2:
                    for vv in exponents:
                        mult=vs.count(vv)
                        ll=list()
                        for t in [0,1]: # even(0) or odd(1) type
                            for det in [1,3,5,7]: # the possible determinants
                                if mult % 2 == 0 and t==0:
                                    ll.append([vv,mult,det,0,0])
                                if mult==1:
                                    odds=[det]
                                elif mult==2:
                                    if det in [1,7]:
                                        odds=[0,2,6]
                                    else:
                                        odds=[2,4,6]
                                else:
                                    odds=[o for o in range(8) if o%2==mult%2]
                                for oddity in odds:
                                    if t==1:
                                        ll.append([vv,mult,det,1,oddity])
                                    #else:
                                        #ll.append([vv,1,det,0,0])
                                        #if mult % 2 == 0 and mult>2:
                                        #    for x in range(1,Integer(mult)/Integer(2)):
                                        #        if mult-2*x==2 and det in [1,7] and oddity not in [0,2,6]:
                                        #            continue
                                        #        elif mult-2*x==2 and det in [3,5] and oddity not in [2,4,6]:
                                        #            continue
                                        #        ll.append([[vv,2*x,det,0,0],[vv,mult-2*x,det,1,oddity]])                                                
                        #print "ll:\n",ll
                        if len(l)==0:
                            for t in ll:
                                if type(t[0])==list:
                                    l.append({p: t})
                                else:
                                    l.append({p: [t]})
                        else:
                            newl=list()
                            for t in ll:
                                for sym in l:
                                    newsym = deepcopy(sym)
                                    #print newsym
                                    if type(t[0])==list:
                                        newsym[p]=newsym[p]+t
                                    else:
                                        newsym[p].append(t)
                                    #print newsym
                                    newl.append(newsym)
                                    #print l
                            l=newl
                        #print "l:\n",l
                else:
                    for vv in exponents:
                        ll=[[vv,vs.count(vv),1],[vv,vs.count(vv),-1]]
                        if len(l)==0:
                            for t in ll:
                                l.append({p: [t]})
                        else:
                            newl=list()
                            for t in ll:
                                for sym in l:
                                    sym[p].append(t)
                                    newl.append(sym)
                            l=newl
                #print "l=\n",l
                #print "psymbols=\n",psymbols
                #print psymbols+l
                psymbols=psymbols+l
        if len(symbols)==0:
            symbols=psymbols
        else:
            symbols_new=list()
            for sym in symbols:
                for psym in psymbols:
                    newsym=deepcopy(sym)
                    newsym.update(psym)
                    symbols_new.append(newsym)
            symbols=symbols_new
    global_symbols = []
    for sym in symbols:
        #print sym
        for p in sym.keys():
            prank = sum([s[1] for s in sym[p]])
            v = sum([ s[0]*s[1] for s in sym[p] ])
            Dp=D//(p**v)
            if prank != rank:
                eps= (Dp*Integer(prod([ s[2] for s in sym[p] ]))).kronecker(p)
                if p==2:
                    if eps==-1:
                        eps=3
                    sym[p].insert(0,[0,rank - prank, eps,0,0])
                else:
                    if eps==-1:
                        for x in Zmod(p):
                            if not x.is_square():
                                eps=x
                                break
                    sym[p].insert(0,[0,rank - prank, eps])
        symbol=GenusSymbol_global_ring(MatrixSpace(ZZ,rank,rank).one())
        symbol._local_symbols=[Genus_Symbol_p_adic_ring(p,syms) for p,syms in sym.items()]
        symbol._signature=(2,n)
        #print symbol._local_symbols
        isglob=is_GlobalGenus(symbol)
        if isglob:
            #print "GLOBAL SYMBOL:"
            #print symbol._local_symbols
            #return symbol
            #for s in symbol._local_symbols:
            #    s = s.canonical_symbol()
            append=True
            for j,s in enumerate(symbol._local_symbols):
                if s._prime==2:
                    sc=deepcopy(symbol)
                    sc._local_symbols[j]=sc._local_symbols[j].canonical_symbol()
                    if csymbols.count(sc)>0:
                        append=False
                    else:
                        csymbols.append(sc)    
                    break
            if append:
                global_symbols.append(symbol)
    return global_symbols
       

def search_for_simple_lattices(n=3,min_D=2,max_D=100):
    simple_symbols=dict()
    rank=2+n
    sign=2-n
    k = Integer(2+n)/Integer(2)
    for D in range(min_D,max_D+1):
        global_symbols=search_global_symbols(n,D)                
        #global_symbols=Set(global_symbols)
        if len(global_symbols)>0:
            print("Symbols for D ={0} : {1}".format(D,len(global_symbols)))
        for sym in global_symbols:
            symstr=get_symbol_string(sym)
            print(symstr)
            sys.stdout.flush()
            if is_simple(symstr):
                print(symstr+"simple!")
                if D not in simple_symbols:
                    simple_symbols[D]=list()
                simple_symbols[D].append(symstr)
                #print sym._local_symbols[0].canonical_symbol()
    return simple_symbols

def is_simple(symstr):
    M=FiniteQuadraticModule(symstr)
    V=VectorValuedModularForms(M)
    d=V.dimension_cusp_forms(V.weight())
    if d == 0:
        return True
    else:
        return False

def get_symbol_string(sym):
    symstr = ''
    #print sym._local_symbols
    for lsym in sym._local_symbols:
        p=lsym._prime
        for s in lsym.symbol_tuple_list():
            if s[0]==0:
                continue
            if len(symstr)!=0:
                symstr=symstr + '.'
            symstr = symstr + str(p**s[0])
            if p == 2:
                sgn = '+' if (Integer(s[2]).kronecker(2) == 1) else '-' 
                if s[3]==1:
                    symstr = symstr + '_' + str(s[4])
                symstr = symstr + '^' + sgn + str(s[1])
                    #else:
                    #    symstr = symstr + '^' + str(s[1])
            else:
                sgn = '+' if (s[2] == 1)  else '-'
                symstr = symstr + '^' + sgn + str(s[1])
    return symstr

def compare_vv_scalar(V,k):
    dv=V.dimension_cusp_forms(k)
    N=V._level
    m=V._M.order()
    if k in ZZ:
        D=DirichletGroup(N)
        if is_even(k):
            chi = D(kronecker_character(m))
        else:
            chi = D(kronecker_character(-m))
        S=CuspForms(chi,k)
        N=Newforms(chi,k,names='a')
        ds=S.dimension()
        B=N       
    else:
        D=magma.DirichletGroup(V._level)
        chi=magma.D(magma.KroneckerCharacter(2*m))
        M=magma.HalfIntegralWeightForms(chi,k)
        S=M.CuspidalSubspace()
        ds=S.Dimension()
        B=S.Basis()
    return dv,ds,S,B

def compare_spaces(n,min_D,max_D):
    spaces=dict()
    rank=2+n
    sign=2-n
    k = Integer(2+n)/Integer(2)
    bad_symbols=dict()
    for D in range(min_D,max_D+1):
        global_symbols=search_global_symbols(n,D)                
        #global_symbols=Set(global_symbols)
        if len(global_symbols)>0:
            print("Symbols for D ={0} : {1}".format(D,len(global_symbols)))
        for sym in global_symbols:
            symstr=get_symbol_string(sym)
            print(symstr)
            sys.stdout.flush()
            M=FiniteQuadraticModule(symstr)
            V=VectorValuedModularForms(M)
            dv,ds,S,B=compare_vv_scalar(V,k)
            if dv==0 and ds!=0 or dv!=0:
                print("dim S({0},{1},rho*) = {2}, dim S(Gamma0({3}),chi_D,{1})={4}".format(symstr,k,dv,V._level,ds))
                if ds != 0:
                    W=V._W
                    f,index=find_lifting_candidate(W,B)
                    #spaces[D]={"W": W, "dimW": dv, "S": S, "dimS":ds, "new_lifting_candidate": f and f.is_new()}
                    if f: #and f.is_new:
                        print("Newform lifting candidate: {0} \n index={1}".format(f.qexp(index+3), index))
                    #elif f:
                    #    print "Oldform lifting candidate found."
                    else:
                        print("**** No newform lifting candidate")
                        if V._level not in bad_symbols:
                            bad_symbols[V._level]=list()
                        bad_symbols[V._level].append(symstr)
                    #else:
                    #    print "No lifting candidate, basis:"
                    #    S
                elif dv==0:
                    if V._level not in bad_symbols:
                        bad_symbols[V._level]=list()
                    bad_symbols[V._level].append(symstr)
            else:
                if dv==0 and ds==0:
                    print("simple and scalar valued space = 0")
    return bad_symbols

def find_lifting_candidate(W,B,max=2):
    vals=list(W.finite_quadratic_module().values())
    vals=list(vals)
    vals=[1-x for x in vals] # we consider the dual representation!
    vals.sort()
    for n in range(max):
        for m in vals:
            index=n*W.level()+m*W.level()
            for f in B:
                if f[index]!=0:
                    #if f.is_new():
                    return f,index
                    #else:
                    #    f_cand=f
                    #    index_cand=index
    return None, index

def get_values(M,dual=False):
    vals=list(M.values())
    vals=list(vals)
    if dual:
        vals=[1-x for x in vals] # we consider the dual representation!
    vals.sort()
    return vals

def compare_formulas_1(D,k):
    DG=DirichletGroup(abs(D))
    chi=DG(kronecker_character(D))
    d1=dimension_new_cusp_forms(chi,k)
    #if D>0:
    #    lvals=sage.lfunctions.all.lcalc.twist_values(1,2,D)
    #else:
    #    lvals=sage.lfunctions.all.lcalc.twist_values(1,D,0)
    #s1=RR(sum([sqrt(abs(lv[0]))*lv[1]*2**len(prime_factors(D/lv[0])) for lv in lvals if lv[0].divides(D) and Zmod(lv[0])(abs(D/lv[0])).is_square()]))
    #d2=RR(1/pi*s1)
    d2=0
    for d in divisors(D):
        if is_fundamental_discriminant(-d):
            K=QuadraticField(-d)
            DD=old_div(ZZ(D),ZZ(d))
            ep=euler_phi((chi*DG(kronecker_character(-d))).conductor())
            #ep=euler_phi(squarefree_part(abs(D*d)))
            print("ep=", ep, D, d)
            ids=[a for a in K.ideals_of_bdd_norm(-DD)[-DD]]
            eulers1=[]
            for a in ids:
                e=a.euler_phi()
                if e!=1 and ep==1:
                    if K(-1).mod(a)!=K(1).mod(a):
                        e=old_div(e,(2*ep))
                else:
                    e=old_div(e,ep)
                eulers1.append(e)
            print(eulers1, ep)
            s=sum(eulers1)
            if ep==1 and not (d.divides(DD) or abs(DD)==1):
                continue
            print(d, s)
            if len(eulers1)>0:
                d2+=s*K.class_number()
    return d1-d2

def compare_formulas_2(D,k):
    d1=old_div(RR(abs(D)),RR(6))
    if D<0:
        D=-D
    s1=RR(sqrt(abs(D))*sum([log(d) for d in divisors(D) if is_fundamental_discriminant(-d) and kronecker(-d,old_div(D,d))==1]))
    d2=RR((old_div(2,(sqrt(3)*pi)))*s1)
    return d1-d2,d2,RR(2*sqrt(D)*log(D)/pi)

A3=lambda N: sum([kronecker(-N,x) for x in Zmod(N) if x**2+x+Zmod(N)(1) == 0])

def compare_formulas_2a(D,k):
    d1=dimension_new_cusp_forms(kronecker_character(D),k)
    if D<0:
        D=-D
    d2=RR(1/pi*sqrt(D)*sum([log(d)*sigma(old_div(D,d),0) for d in divisors(D) if Zmod(d)(old_div(D,d)).is_square() and is_fundamental_discriminant(-d)]))
    return d1-d2

def compare_formulas_3(D,k):
    d1=dimension_cusp_forms(kronecker_character(D),k)
    d2=RR(sqrt(abs(D))*log(abs(D))/pi + sqrt(abs(D))*log(abs(D))*9/(sqrt(3)*pi))
    return d1-d2

def formtest(minD,maxD,k,eps=-1):
    ds={}
    for D in range(minD,maxD+1):
        if (eps*D)%4==1:
            dif,d,dd=compare_formulas_2(eps*D,k)
            if dif<=1e-6:
                print(D, factor(D))
                print(dif,d,dd)
            #ds[D]=(d[1]-d[2])/RR(D)
    #return ds

def formtest_2(minD,maxD):
    s=0
    for D in range(minD,maxD):
        if is_fundamental_discriminant(-D):
            for d in divisors(D):
                if is_fundamental_discriminant(-d):
                    sd=RR(RR(1)/RR(6)*(RR(d)+old_div(RR(D),RR(d)))-RR(sqrt(D))/RR(pi)*log(d))
                    print(D, d, sd)
                    s+=sd
        if s<=0:
            print("s= {0}  D={1}".format(s,D))


def CMorbits(D):
    s=sum([2**(len(prime_divisors(old_div(D,d)))) for d in divisors(D) if is_fundamental_discriminant(-d) and Zmod(d)(old_div(D,d)).is_square()])+1
    return s

def orbittest(minD,maxD,eps=-1):
    for D in range(minD,maxD+1):
        if is_fundamental_discriminant(eps*D) and is_odd(D):
            print(D)
            print(CMorbits(eps*D)-len(Newforms(kronecker_character(eps*D),3,names='a')))

def sigma_rep(Delta,print_divisors=False):
    s=0
    for DD in divisors(Delta):
        if is_fundamental_discriminant(-DD):
            D=-DD
            for d in divisors(old_div(Delta,D)):
                s+=kronecker(D,d)
                if print_divisors:
                    print(D,d, kronecker(D,d))
    return s


def test_sigma_rep(Dmin=1,Dmax=1000,print_divisors=False):
    n=0
    sa=0
    sw=[1,0]
    sm=[1,0]
    for D in range(Dmin,Dmax):
        if is_fundamental_discriminant(-D):
            n+=1
            sr=sigma_rep(D,True)
            sa+=sr
            sm[1]=max(sr,sm[1])
            swt=old_div(sr,D)
            sw[1]=max(swt,RR(sw[1]))
            if sm[1]==sr:
                sm[0]=-D
            if sw[1]==swt:
                sw[0]=-D
            ct=RR(100)*(RR(D)**(old_div(1.0,25)))
            #if ct<=RR(sr):
            print(-D, sr, ct, sigma(D,0))
    sa=old_div(sa,n)
    print("Max: {0}".format(sm))
    print("Avg: {0}, {1}".format(sa, RR(sa)))
    print("Worst case: {0}".format(sw))

def neg_fdd_div(D):
    s=0
    for d in divisors(D):
        if is_fundamental_discriminant(-d):
            if Zmod(d)(old_div(D,d)).is_square():
                s+=1
    return s

def test_neg_fdd_div(minD,maxD):
    for D in range(minD,maxD):
        if is_fundamental_discriminant(-D):
            print("{0}, {1}, {2}".format(D, neg_fdd_div(D), old_div(sigma(D,0),2)))

def sigma_log(D):
    return sum([RR(log(d)) for d in divisors(D)])
