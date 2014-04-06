def invariants_eps(FQM, TM, use_reduction = False, debug = 0):
    eps = True
    if TM != None and FQM != None:
        TMM = TM+FQM
    elif TM != None:
        TMM = TM
    else:
        TMM = FQM
        eps = False
    if debug > 1: print "FQM = {0}, TM = {1}, TMM = {2}".format(FQM, TM, TMM)
    if debug > 2:
        debug2=debug
    else:
        debug2=0
    inv = invariants(TMM, use_reduction, debug2)
    if debug > 1: print inv
    if type(inv) in [list,tuple]:
        V = inv[1]
    else:
        V = inv
    d = [0,0]
    if V.dimension() != 0:
        el = list()
        M = Matrix(V.base_ring(), V.ambient_module().dimension())
        if eps:
            for v in inv[0]:
                vv = v.c_list()
                vv[0] = -vv[0]
                vv = TMM(vv,can_coords=True)
                if inv[0].count(vv) > 0:
                    el.append(inv[0].index(vv))
                else:
                    el.append(inv[0].index(-vv))
            for a in el:
                M[el[a],a]=1
            if debug > 1: print M
            try:
                KM = (M-M.parent().one()).kernel_on(V)
                if debug > 1: print KM
                d[0] = KM.dimension()
                KM = (M+M.parent().one()).kernel_on(V)
                if debug > 1: print KM
                d[1] = KM.dimension()
            except:
                print "Error occured for ", FQM.jordan_decomposition().genus_symbol()
        else:
            d = V.dimension()
    if debug > 1: print d
    return d

def weight_one_half_dim(FQM, use_reduction, debug = 0, local=True):
    N = Integer(FQM.level())
    if not N % 4 == 0:
        return 0
    m = Integer(N/Integer(4))
    d = 0
    for l in m.divisors():
        if is_squarefree(m/l):
            if debug > 1: print "l = {0}".format(l)
            TM = FiniteQuadraticModule([2*l],[-1/Integer(4*l)])
            if local:
                dd = [0,0] # eigenvalue 1, -1 multiplicity
                for p,n in lcm(FQM.level(),4*l).factor():
                    N = None
                    TN = None
                    J = FQM.jordan_decomposition()
                    L = TM.jordan_decomposition()
                    for j in xrange(1,n+1):
                        C = J.constituent(p**j)[0]
                        D = L.constituent(p**j)[0]
                        if debug > 1: print "C = {0}, D = {1}".format(C,D)
                        if N == None and C.level() != 1:
                            N = C
                        elif C.level() != 1:
                            N = N + C
                        if TN == None and D.level() != 1:
                            TN = D
                        elif D.level() != 1:
                            TN = TN + D
                    dd1 = invariants_eps(N, TN, use_reduction, debug)
                    if debug > 1: print "dd1 = {}".format(dd1)
                    if dd1 == [0,0]:
                        dd = [0,0]
                        break
                    ddtmp = copy(dd)
                    if dd[0] == 0:
                        ddtmp[0] = dd1[0]
                    else:
                        ddtmp[0] = dd[0]*dd1[0] + dd[1]*dd1[1]
                    if dd[1] == 0:
                        ddtmp[1] = dd1[1]
                    else:
                        ddtmp[1] = dd[0]*dd1[1] + dd[1]*dd1[0]
                    dd = ddtmp
                    if debug > 1: print "dd = {0}".format(dd)
                d += dd[0]
            else:
                d += invariants_eps(FQM, TM, use_reduction, debug)[0]
    return d
