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




class J_generic:
    def __init__(self, K, N_max=10**5):
        """
        Compute J working over the field K.
        """
        self.N_max   = N_max
        self.K = K

        from sage.all import prime_powers, factor
        PP = prime_powers(N_max+1)[1:]

        n  = len(PP)
        self.a   = [K(0)]*n
        self.s   = [K(0)]*n
        self.pv  = [0]*n
        i = 0
        for pv in PP:
            F = factor(pv)
            p, v = F[0]
            self.pv[i] = K(pv)
            logp = K(p).log()
            self.a[i] = logp/K(pv).sqrt()
            self.s[i] = v*logp
            i += 1

    def __repr__(self):
        return "The function J(t,N) over %s of the Mazur-Stein game with N_Max=%s"%(
            self.K, self.N_max)

    def F(self, t, N):
        K = self.K
        t = K(t)
        ans = K(0)
        i = 0
        if N > self.N_max:
            raise ValueError, "J not defined for N > %s"%self.N_max
        while 1:
            if self.pv[i] >= N:
                return ans
            ans += self.a[i] * (self.s[i]*t).cos()
            i += 1
            
    def e(self, t, N):
        K = self.K
        return K(N).sqrt() / (.25 + t*t) * \
               ( (K(t*K(N).log())).cos() + 2*t*(K(t*K(N).log())).sin() )

    def H(self, t, N):
        K = self.K
        ans = K(0); F = K(0)
        i = 0
        if N > self.N_max:
            raise ValueError, "J not defined for N > %s"%self.N_max
        n = 1
        while 1:
            if self.pv[i] > N:
                # time to return.  But first we have to add a few more
                # values to ans up to N (where F is constant).
                while n <= N:
                    ans += F + self.e(t,n)
                    n += 1
                # Finally return after dividing by N.
                return ans / N
            # At this point, F is equal to self.F(t, self.pv._values[i]),
            # and now we add to our running total a range of values of

            #    F(t, n) + e(t,n)
            while n <= self.pv[i]:
                ans += F + self.e(t,n)
                n += 1

            F += self.a[i] * cos(self.s[i]*t)

            # move to next prime power
            i += 1

    def J(self, t, N):
        K = self.K
        return self.H(t,N) / K(N).log()

    def __call__(self, t, N):
        return self.J(t,N)

    ######################################################
    # Cesaro of Cesaro...
    ######################################################    
    def Jc(self, t, N):
        K = self.K
        ans = K(0); F = K(0); ans2 = K(0)
        i = 0
        if N > self.N_max:
            raise ValueError, "J not defined for N > %s"%self.N_max
        n = 1
        while 1:
            if self.pv[i] > N:
                # time to return.  But first we have to add a few more
                # values to ans up to N (where F is constant).
                while n <= N:
                    ans += F + self.e(t,n)
                    if n >= 2:
                        ans2 += ans / (n*K(n).log()) 
                    #print (n, F, self.e(t,n))
                    n += 1
                # Finally return after dividing by N.
                return ans2 / N
            # At this point, F is equal to self.F(t, self.pv._values[i]),
            # and now we add to our running total a range of values of

            #    F(t, n) + e(t,n)
            while n <= self.pv[i]:
                ans += F + self.e(t,n)
                if n >= 2:
                    ans2 += ans / (n*K(n).log())
                #print (n, F, self.e(t,n))                
                n += 1

            F += self.a[i] * (self.s[i]*t).cos()

            # move to next prime power
            i += 1
