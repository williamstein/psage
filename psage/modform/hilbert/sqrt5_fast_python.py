def ideal_characteristic(I):
    return int(I.number_field()._pari_().idealtwoelt(I._pari_())[0])

def ResidueRing(P, e):
    from sqrt5_fast import ResidueRing_nonsplit, ResidueRing_ramified_odd, ResidueRing_split
    p = ideal_characteristic(P)
    if p == 5:   # ramified
        if e % 2 == 0:
            return ResidueRing_nonsplit(P, p, e)
        else:
            return ResidueRing_ramified_odd(P, p, e)
    elif p%5 in [2,3]:  # nonsplit
        return ResidueRing_nonsplit(P, p, e)
    else: # split
        return ResidueRing_split(P, p, e)





