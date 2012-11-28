from sage.quadratic_forms.genera.genus import GenusSymbol_global_ring, Genus_Symbol_p_adic_ring, is_GlobalGenus
from psage.modules.finite_quadratic_module import FiniteQuadraticModule
from sage.matrix.matrix_space import MatrixSpace
from sage.all import ZZ
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
                sym = [sym[0], sym[1], sym[2] % 8, 1, sym[3] % 8]
                if sym[1]==1:
                    if  sym[2].kronecker(2)==sym[4].kronecker(2):
                        sym[2]=sym[4]
                    else:
                        return False
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
                symbols[p].append([0,n - prank, eps,0,0])
            else:    
                symbols[p].append([0,n - prank, eps])
    symbol=GenusSymbol_global_ring(MatrixSpace(ZZ,r+s,r+s).one())
    symbol._local_symbols=[Genus_Symbol_p_adic_ring(p,syms) for p,syms in symbols.iteritems()]
    symbol._signature=(r,s)
    #print r,s, symbol._local_symbols
    isglob=is_GlobalGenus(symbol)
    if return_symbol:
        return symbol, isglob
    else:
        return isglob

list_freitag = {
    3: ['2_7^1', '2_3^-1', '2_7^+3', '2^+2.4_7^+1', '2^2.8_3^-1',
        '2_7^+1.4^2', '2^+4.4_7^+1', '2^4.8_3^-1', '2_7^1.4^2',
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
    rank=2+n
    sign=2-n
    k = Integer(2+n)/Integer(2)
    symbols=dict()
    for D in range(min_D,max_D+1):
        if True:
            D=(-1)**n*D
            fac = Integer(D).factor()
            symbols=dict()
            maxlen = 0
            symbols=list()
            for p, v in fac:
                psymbols=list()
                parts=partitions(v)
                Dp=D//(p**v)
                for j,vs in enumerate(parts):
                    #if vs.count(1)>1:
                    #    continue
                    prank=len(vs)
                    if prank <= rank:
                        l = dict()
                        if p==2:
                            l2=dict()
                            for r in range(prank):
                                l2[r]=list()
                                for t in [0,1]:
                                    if t==1 and is_odd(vs[r]):
                                        continue
                                    dets=[1,3,5,7] if t==0 else [7,3]
                                    for d in dets:
                                        l2[r].append((t,d))
                            #print 'l2=', l2
                            for r in range(prank):
                                l[r]=list()
                                for s in l2[r]:
                                    vv = vs[r] if s[0]==0 else Integer(vs[r])/Integer(2)
                                    oddity=0 if s[0]==1 else s[1]
                                    l[r].append((vv,s[0]+1,s[1],1-s[0],oddity))
                        else:
                            for r in range(prank):
                                l[r] = [(vs[r],1,eps) for eps in [-1,1]]
                        #print l
                        index=unordered_tuples(range(len(l[0])),prank)
                        #print index
                        for i in range(len(index)):
                            sym={p: [l[r][index[i][r] % len(l[r])] for r in range(prank)]}
                            #print 'sym=', sym
                            prank_s = sum([ s[1] for s in sym[p]])
                            if prank_s < rank:
                                eps= (Dp*Integer(prod([ s[2] for s in sym[p]]))).kronecker(p)
                                if p==2:
                                    sym[p].append((0,rank - prank_s, eps,0,0))
                                else:    
                                    sym[p].append((0,rank - prank_s, eps))
                            if prank_s <= rank:
                                #print vs, sym
                                psymbols.append(sym)
                if len(symbols)==0:
                    symbols=psymbols
                else:
                    symbols_new=list()
                    for sym in symbols:
                        for psym in psymbols:
                            sym.update(psym)
                            symbols_new.append(sym)
                    symbols=symbols_new
            print "Symbols for D=", D, ":"
            print len(symbols)
            for sym in symbols:
                symbol=GenusSymbol_global_ring(MatrixSpace(ZZ,rank,rank).one())
                symbol._local_symbols=[Genus_Symbol_p_adic_ring(p,syms) for p,syms in sym.iteritems()]
                symbol._signature=(2,n)
                isglob = is_GlobalGenus(symbol)
                #print sym, isglob
                if True or isglob:
                    symstr = ''
                    for p,syms in sym.iteritems():
                        for s in syms:
                            if s[0]==0:
                                continue
                            if len(symstr)!=0:
                                symstr=symstr + '.'
                            symstr = symstr + str(p**s[0])
                            if p == 2:
                                sgn = '+' if (Integer(s[2]).kronecker(2) == 1) else '-' 
                                if s[3]==1:
                                    symstr = symstr + '_' + str(s[2])
                                symstr = symstr + '^' + sgn + str(s[1])
                                #else:
                                #    symstr = symstr + '^' + str(s[1])
                            else:
                                sgn = '+' if (s[2] == 1)  else '-'
                                symstr = symstr + '^' + sgn + str(s[1])
                    #print symstr, isglob
                    #if symstr == '2^+2.4_7^+1':
                    #    return
                    M=FiniteQuadraticModule(symstr)
                    #isglob=is_global(M,2,n)
                    if isglob:
                        V=VectorValuedModularForms(M)
                        d=V.dimension_cusp_forms(k)
                        if d == 0:
                            print symstr, isglob
                    
                    
                        
                
                    
            
    
            
