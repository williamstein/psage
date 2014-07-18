from sage.all import SageObject, Integer, RR, is_odd, next_prime, floor, RealField, ZZ, ceil, log, ComplexField, real, sqrt, exp, is_squarefree, lcm, Matrix
from psage.modules.finite_quadratic_module import FiniteQuadraticModule
from psage.modules.weil_invariants import invariants
from copy import copy

def invariants_eps(FQM, TM, use_reduction = True, proof = False, debug = 0):
    r"""
    Computes the invariants of a direct summand in the decomposition for weight
    one-half modular forms. Such a summand is of the form
    \[
      \mathbb{C}[(\mathbb{Z}/2N\mathbb{Z}, -x^2/4N)]^\epsilon \otimes \mathbb{C}[M],
    \]
    where $M$ is a given finite quadratic module and $\epsilon$ acts on the first
    factor via $\mathfrak{e}_\mu \mapsto \mathfrak{e}_{-\mu}$.

    INPUT:
        - FQM: A given finite quadratic module, referred to as $M$ above
        - TM: A cyclic module of the form as given abve (the first factor)

    NOTE:
        We do not check that TM is of the stated form. The function is an auxiliary function
        usually only called by ``weight_one_half_dim``.

    EXAMPLES:
        NONE
    """
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
    inv = invariants(TMM, use_reduction, proof, debug2)
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
            f = 1 if TMM.signature() % 4 == 0 else -1
            for v in inv[0]:
                vv = v.c_list()
                vv[0] = -vv[0]
                vv = TMM(vv,can_coords=True)
                if inv[0].count(vv) > 0:
                    el.append((inv[0].index(vv),1))
                else:
                    el.append((inv[0].index(-vv),f))
            for i in range(len(el)):
                M[el[i][0],i] = el[i][1]
            if debug > 1: print "M={0}, V={1}".format(M, V)
            try:
                KM = (M-M.parent().one()).kernel_on(V)
                if debug > 1: print "KM for ev 1 = {0}".format(KM)
                d[0] = KM.dimension()
                KM = (M+M.parent().one()).kernel_on(V)
                if debug > 1: print "KM for ev -1 = {0}".format(KM)
                d[1] = KM.dimension()
            except:
                raise RuntimeError("Error occured for {0}".format(FQM.jordan_decomposition().genus_symbol()), M, V)
        else:
            d = [V.dimension(), 0]
    if debug > 1: print d
    return d

def weight_one_half_dim(FQM, use_reduction = True, proof = False, debug = 0, local=True):
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
                    dd1 = invariants_eps(N, TN, use_reduction, proof, debug)
                    if debug > 1: print "dd1 = {}".format(dd1)
                    if dd1 == [0,0]:
                        # the result is multiplicative
                        # a single [0,0] as a local result
                        # yields [0,0] in the end
                        # and we're done here
                        dd = [0,0]
                        break
                    if dd == [0,0]:
                        # this is the first prime
                        dd = dd1
                    else:
                        # some basic arithmetic ;-)
                        # 1 = 1*1 = (-1)(-1)
                        # -1 = 1*(-1) = (-1)*1
                        ddtmp = copy(dd)
                        ddtmp[0] = dd[0]*dd1[0] + dd[1]*dd1[1]
                        ddtmp[1] = dd[0]*dd1[1] + dd[1]*dd1[0]
                        dd = ddtmp
                    if debug > 1: print "dd = {0}".format(dd)
                d += dd[0]
            else:
                d += invariants_eps(FQM, TM, use_reduction, proof, debug)[0]
    return d
