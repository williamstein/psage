from sage.all import *


def examples():
    r"""
    Some more examples of xxx_rankbound() on high rank curves.
    (See source for more.)

    EXAMPLES:

    We make sure that the examples all work.

    ::

        sage: from psage.ellcurve.xxx.rankbound_examples import examples
        sage: examples() # LONG TIME (lack of output indicates success)

    """
    from psage.ellcurve.xxx.rankbound import xxx_rankbound

    def check(E, Delta, answer, bad_primes = None):
        # TODO: Is there a way to quickly compute the root number
        # given the bad primes?
        if bad_primes is not None:
            # check that we actually have all of the bad primes
            D = E.discriminant()
            for p in bad_primes:
                p = Integer(p)
                assert(p.divides(D)), "p is not bad"
                while p.divides(D):
                    D = D.divide_knowing_divisible_by(p)

            assert abs(D) == 1, "Missing a bad prime."

        
        a = xxx_rankbound(E, Delta, bad_primes)

        assert abs(a - answer) < 1e-3, "FAIL"

    E20 = EllipticCurve([1, 0, 0, -431092980766333677958362095891166, 5156283555366643659035652799871176909391533088196])
    E21 = EllipticCurve([1, 1, 1, -215843772422443922015169952702159835, -19474361277787151947255961435459054151501792241320535])
    E22 = EllipticCurve([1, 0, 1, -940299517776391362903023121165864, 10707363070719743033425295515449274534651125011362])
    E23 = EllipticCurve([1, 0, 1, -19252966408674012828065964616418441723, 32685500727716376257923347071452044295907443056345614006])
    E24 = EllipticCurve([1, 0, 1, -120039822036992245303534619191166796374, 504224992484910670010801799168082726759443756222911415116])
    E28 = EllipticCurve([1, -1, 1, -20067762415575526585033208209338542750930230312178956502,
  34481611795030556467032985690390720374855944359319180361266008296291939448732243429])

    E23_bad_primes = [2, 3, 5, 11, 13, 17, 19, 23, 199, 17858193374841750901974789649, 1006218106634655545344494448610726356220703995276273]
    E24_bad_primes = [2, 3, 5, 11, 13, 17, 29, 31, 41, 458619970494582607679296750333015081, 264240973182971699094661154229360236070105974082503]
    E28_bad_primes = [2, 3, 5, 7, 11, 13, 17, 19, 48463, 20650099, 315574902691581877528345013999136728634663121, 376018840263193489397987439236873583997122096511452343225772113000611087671413]

    check(E20, 2.0, 21.6907) 

    check(E21, 2.5, 22.6727)

    check(E22, 2.0, 23.7047)

    check(E23, 2.5, 24.4834, E23_bad_primes)

    check(E24, 2.5, 25.5682, E24_bad_primes) 
    
    check(E28, 2.5, 33.4304, E28_bad_primes) # needs a larger delta to get a
                                             # good bound, but that takes a
                                             # while.

    #check(E28, 3.2, 31.2984, E28_bad_primes) # this works, but takes about
                                              # 6 minutes on my fast machine.
                                              # (Is that too long for LONG?)
