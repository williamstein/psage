June 1, 2008 smalljac version 3.0 readme.txt

The major improvements in smalljac version 3 are in genus 1 performance.

Specialized root-finding algorithms for cubics and quartics over finite
fields made the computation of 3-torsion competitive for all but very small p
(the time to factor a quartic using standard methods exceeds the time
required to compute a_p until p is quite large (outside the feasible range),
This is why we didn't use 3-torsion information previously).

Together with a number of improvements to the genus 1 point arithmetic,
performance has improved by more than a factor of 2.  3-torsion information
is also used, in a more limited way, in genus 2, via the modular polynomial
of Gaudry and Schost (this is practical only for largish p, say > 10^6).

Also included are two further optimizations suggested at ANTS 8: the use of
4-torsion information via point-halving (thanks to Noam Elkies), and a judicious
selection of the first giant step in the interval which minimizes the hamming weight
of the scalar multiple used to compute it (thanks to Dan Bernstein).

smalljac version 3 also contains several bug fixes, most of which arose only for
exceptional curves (Jacobians with extra endomorphisms) and/or small primes.

The following files are needed to compile the genus 1and genus 2 versions
of smalljac, including demonstration programs.  For genus 3, David Harvey's
zn_poly and NTL are also required and not included here.

GMP is required, as is GCC (the core genus 1 code does not need GMP).
This code is only supported on 64-bit Intel/AMD platforms (life is too short
to spend it worrying about overflowing a 32-bit integer).

demo1.c
demo2.c
demo3.c
demo4.c
lpoly.c
jgroup.c
smalljac.c
smalljac_g23.c
smalljactab.c
jacorder.c
jacstructure.c
jac.c
hecurve.c
hecurve1.c
hecurve1x.c
hecurve2.c
hecurve3.c
pointcount.c 
ffpoly.c
ffpoly_g2tor3.c
ffpolyext.c
ff.c
ffext.c
mpzutil.c
cstd.c

smalljac.h
smalljac_g23.h
smalljac_internal.h
smalljactab.h
jacorder.h
jac.h
hecurve.h
pointcount.h
polydisc.h
ffwrapper.h
ffpoly.h
g2tor3poly.h
ff.h
asm.h
mpzutil.h
cstd.h

The script builddemo compiles and links the demo programs: demo1, demo2,
demo3, demo4, lpoly, and jgroup (yes, I will eventually put together a makefile).

The main interface for applications is defined in smalljac.h, and is used by
the demo programs above.

The finite field representation is now set by default to 32 bits, meaning that
prime fields up to 2^31 in size may be used.  This can be changed to 64
bits by defining SMALLJAC_FF_64BITS to be 1 in smalljac.h.  Doing so
will decrease performance by five to ten percent but increase the maximum
field size to 2^63.

Currently, only genus 1 and 2 curves are fully supported
(but in genus 2 only curves with a single point at infinity).
Please report bugs (on 64-bit platforms) to drew@math.mit.edu.

Eventually genus 3 curves will be fully supported (they work now only in
the typical case).  I also hope to add support for genus 2 curves with zero or
two points at infinity (real hyperelliptic curves).

The following references occur in the code where specific algorithms are
implemented.  This is by no means a comprehensive bibliography, The
code does not generally cite references for well known algorithms (e.g.
Shanks' baby-steps giant-steps algorithm and many others), except
where a direct implementation of a specific version was used.

[CANT] Henri Cohen,
    "A Course in Computational Algebraic Number Theory", Springer, 1996.

[HECHECC] Henri Cohen (ed.) et al.,
    "Handbook of Elliptic Curve and Hyperelliptic Curve Cryptography",
    Chapman and Hall, 2006.
    
[KedlayaSutherland2007] Kiran S. Kedlaya and Andrew V. Sutherland,
    "Computing L-series of hyperelliptic curves", to appear in ANTS VIII,
    2008, http://arxiv.org/abs/arXiv:0801.2778v2.

[SutherlandThesis] Andrew V. Sutherland,
    "Order computations in generic groups", Ph.D. thesis, MIT, 2007,
    //groups.csail.mit.edu/cis/theses/sutherland-phd.pdf.
    
[Sutherland2007] Andrew V. Sutherland,
    "A generic approach to searching for Jacobians", to appear in Math. Comp.,
    2007, http://arxiv.org/abs/0708.3168v2.
