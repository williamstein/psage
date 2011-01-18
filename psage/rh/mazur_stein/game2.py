from sage.all import prime_range, I

class J_E:
    def __init__(self, E, N, K):
        self.E = E
        self.ap = [K(ap) for ap in E.aplist(N)]
        self.N = N
        self.primes = [K(p) for p in prime_range(N)]
        self.K = K
        self.log_primes = [p.log() for p in self.primes]
        self.I = K(I)

    def F0(self, t, N):
        assert N <= self.N
        K = self.K
        s = 1 + t*self.I  # 1 is the center
        ans = 0
        i = 0
        while self.primes[i] < N:
            ap  = self.ap[i]
            p   = self.primes[i]
            X   = p**(-s)
            ans += (ap - 2*p) / (1 - ap*X + p*X*X)  * X * self.log_primes[i]
            i   += 1
        return ans


