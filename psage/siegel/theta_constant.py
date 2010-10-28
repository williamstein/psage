"""
Theta constants
"""


def _compute_theta_char_poly(char_dict, f):
    r"""
    Return the coefficient at ``f`` of the Siegel modular form

    .. math:

        \sum_{l \in\text{char_dict}} \alpha[l] * \prod_i \theta_l[i](8Z).   

    INPUT:

    - ``char_dict`` -- a dictionary whose keys are *tuples* of theta 
      characteristics and whose values are in some ring.
    - ``f`` -- a triple `(a,b,c)` of rational numbers such that the
      quadratic form `[a,b,c]` is semi positive definite 
      (i.e. `a,c \geq 0` and `b^2-4ac \leq 0`).

    EXAMPLES::

        sage: from sage.modular.siegel.theta_constant import _multiply_theta_char
        sage: from sage.modular.siegel.theta_constant import _compute_theta_char_poly
        sage: theta_constants = {((1, 1, 0, 0), (0, 0, 1, 1), (1, 1, 0, 0), (0, 0, 1, 1)): 1}
        sage: _compute_theta_char_poly(theta_constants, [2, 0, 10])
        32
        sage: _compute_theta_char_poly(theta_constants, [2, 0, 2]) 
        8
        sage: _compute_theta_char_poly(theta_constants, [0, 0, 0])
        0
    """
    return sum(_multiply_theta_char(list(l), f)*char_dict[l] for l in char_dict)



def _multiply_theta_char(l, f):
    r"""
    Return the coefficient at ``f`` of the theta series `\prod_t \theta_t` 
    where `t` runs through the list ``l`` of theta characteristics.

    INPUT:

    - ``l`` -- a list of quadruples `t` in `{0,1}^4`
    - ``f`` -- a triple `(a,b,c)` of rational numbers such that the
      quadratic form `[a,b,c]` is semi positive definite 
      (i.e. `a,c \geq 0` and `b^2-4ac \leq 0`).

    EXAMPLES::

        sage: from sage.modular.siegel.theta_constant import _multiply_theta_char
        sage: from sage.modular.siegel.theta_constant import _compute_theta_char_poly
        sage: tc = [(1, 1, 0, 0), (0, 0, 1, 1), (1, 1, 0, 0), (0, 0, 1, 1)]
        sage: _multiply_theta_char(tc, [2, 0, 6])
        -32
        sage: _multiply_theta_char(tc, [2, 0, 10])
        32
        sage: _multiply_theta_char(tc, [0, 0, 0]) 
        0
    """
    if 0 == len(l):
        return (1 if (0, 0, 0) == f else 0)
    a, b, c = f
    m1, m2, m3, m4 = l[0]
    # if the characteristic is not even:
    if 1 == (m1*m3 + m2*m4)%2: 
        return 0
    coeff = 0
    from sage.misc.all import isqrt, xsrange
    for u in xsrange(m1, isqrt(a)+1, 2):
        for v in xsrange(m2, isqrt(c)+1, 2):
            if 0 == u and 0 == v:
                coeff += _multiply_theta_char(l[1:], (a, b, c))
                continue
            ap, bp, cp = (a-u*u, b-2*u*v, c-v*v)
            if bp*bp-4*ap*cp <= 0:
                val = (2 if  0 == (u*m3 + v*m4)%4 else -2)
                coeff += val * _multiply_theta_char(l[1:], (ap, bp, cp))
            if u != 0 and v != 0:
                ap, bp, cp = (a-u*u, b+2*u*v, c-v*v)
                if bp*bp-4*ap*cp <= 0:
                    val = (2 if  0 == (u*m3 - v*m4)%4 else -2)
                    coeff += val * _multiply_theta_char(l[1:], (ap, bp, cp))
    return coeff
