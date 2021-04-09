from __future__ import absolute_import
def SmallJac(f):
    try:
        d = f.degree()
    except TypeError:
        d = f.degree(f.variables()[0])
    g = (d - 1) // 2
    if g == 1:
        from . import wrapper1
        return wrapper1.SmallJac(f)
    elif g == 2:
        from . import wrapper2
        return wrapper2.SmallJac(f)
    elif g == 3:
        from . import wrapper3
        return wrapper3.SmallJac(f)
    raise NotImplementedError
