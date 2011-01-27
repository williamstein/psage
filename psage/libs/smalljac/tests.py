
def test1(B=10**4):
    from sage.all import polygen, primes, QQ, EllipticCurve
    import psage.libs.smalljac.wrapper
    x = polygen(QQ, 'x')
    J = psage.libs.smalljac.wrapper.SmallJac(x**3 + 17*x + 3)
    E = EllipticCurve([17,3])
    N = E.conductor()
    for p in primes(B):
        if N%p:
            assert E.ap(p) == J.ap(p)

def test2(B=500):
    from sage.all import polygen, primes, QQ, GF, HyperellipticCurve
    import psage.libs.smalljac.wrapper
    x = polygen(QQ, 'x')
    J = psage.libs.smalljac.wrapper.SmallJac(x**5 + 17*x + 3)
    N = 97*3749861
    for p in primes(45,B):  # 45 comes from "(2g+1)(2N-1) = 45"
        if N%p:
            x = polygen(GF(p), 'x')            
            C = HyperellipticCurve(x**5 + 17*x + 3)
            assert C.frobenius_polynomial() == J.frob(p).charpoly()

    
    
