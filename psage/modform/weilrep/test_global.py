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
    isglob=is_GlobalGenus(symbol)
    if return_symbol:
        return symbol, isglob
    else:
        return isglob

list_freitag = {
    3: ['2_7^1', '2_7^+3', '2^+2.4_7^+1', '2^2.8_3^-1',
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
                print "n = ", n, " ", symbol, ': Dimension('+ str(k) + ') = ', V.dimension_cusp_forms(k), ", ", is_global(M,2,n)

def search_for_simple_lattices(n=3,max_D=100):
    rank=2+n
    sign=2-n
    symbols=dict()
    for D in range(2,max_D+1):
        if D % 4 == 0 or D%4 ==1:
            fac = Integer(D).factor()
            symbols=dict()
            maxlen = 0
            for p, v in fac:
                symbols[p]=dict()
                parts=partitions(v)
                Dp=D//(p**v)
                for j,vs in enumerate(partitions(v)):
                    prank=len(vs)
                    if j > maxlen:
                        maxlen=j
                    if prank <= rank:
                        if p==2:
                            l = {r: [(vs[r],1,d,t,d) for d in [1,3,5,7] for t in [0,1]] for r in range(prank)}
                        else:
                            l = {r: [(vs[r],1,eps) for eps in [0,1]] for r in range(prank)}
                        symbols[p][j]=[[l[i][r] for i in range(len(l))] for r in range(len(l[0]))]
                        if prank < rank:
                            for sym in symbols[p][j]:
                                eps= (Dp*Integer(prod([ s[2] for s in sym ]))).kronecker(p)
                                if p==2:
                                    sym.append((0,rank - prank, eps,0,0))
                                else:    
                                    sym.append((0,rank - prank, eps))
            for j in range(maxlen):
                sym=dict()
                for p in Integer(D).prime_factors():
                    if symbols[p].has_key(j):
                        sym.append(symbols[p][j])
                
    return symbols
                        
                
                    
            
    
            
