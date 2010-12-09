def SmallJac(f):
    x = f.variables()[0]
    g = (f.degree(x)-1) // 2
    if g == 1:
        import wrapper1
        return wrapper1.SmallJac(f)
    elif g == 2:
        import wrapper2
        return wrapper2.SmallJac(f)
    elif g == 3:
        import wrapper3
        return wrapper3.SmallJac(f)
    raise NotImplementedError
