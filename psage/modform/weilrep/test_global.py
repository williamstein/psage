from sage.quadratic_forms.genera.genus import GenusSymbol_global_ring, Genus_Symbol_p_adic_ring, is_GlobalGenus
from psage.modules.finite_quadratic_module import FiniteQuadraticModule
from sage.matrix.matrix_space import MatrixSpace
from sage.all import ZZ, Zmod
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
        if not symbols.has_key(p):
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
    symbol._local_symbols=[Genus_Symbol_p_adic_ring(p,syms) for p,syms in symbols.iteritems()]
    symbol._signature=(r,s)
    print r,s, symbol._local_symbols
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
    8: ['2_2^+2', #'2^-4',
        '3^-1', '7^1'
        ],
    9: ['2_1^+1', '4_1^+1', '8_1^+1'
        ],
    10: ['2^+2'
        ]
}

def test_list_freitag(m=-1):
    all=True if m==-1 else False
    for n,symbols in list_freitag.iteritems():
        if all or n==m:
            for symbol in symbols:
                M=FiniteQuadraticModule(symbol)
                V=VectorValuedModularForms(M)
                k = Integer(2+n)/Integer(2)
                print "n = ", n, " ", symbol, ': Dimension('+ str(k) + ') = ', V.dimension_cusp_forms(k), ", ", is_global(M,2,n), '|M|=', M.order()

def search_for_simple_lattices(n=3,min_D=2,max_D=100):
    simple_symbols=dict()
    rank=2+n
    sign=2-n
    k = Integer(2+n)/Integer(2)
    for D in range(min_D,max_D+1):
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
            symbol._local_symbols=[Genus_Symbol_p_adic_ring(p,syms) for p,syms in sym.iteritems()]
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
                
        #global_symbols=Set(global_symbols)
        if len(global_symbols)>0:
            print "Symbols for D =", D, ": ", len(global_symbols)
        for sym in global_symbols:
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
            print symstr
            M=FiniteQuadraticModule(symstr)
            V=VectorValuedModularForms(M)
            d=V.dimension_cusp_forms(k)
            if d == 0:
                print symstr, "simple!"
                if not simple_symbols.has_key(D):
                    simple_symbols[D]=list()
                simple_symbols[D].append(symstr)
                #print sym._local_symbols[0].canonical_symbol()
    return simple_symbols
    
            
