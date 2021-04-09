r"""

Algorithms for locating Maass waveforms.

"""
from __future__ import print_function
from __future__ import division

from builtins import str
from builtins import map
from builtins import range
from past.utils import old_div
import logging

from sage.all import RealField,RR,gcd,QQ,sqrt,log_b,ceil
from sage.rings.complex_number import ComplexNumber
from sage.plot.all import Graphics
from sage.functions.generalized import sign

from psage.modform.maass.maass_forms import Maasswaveform
from psage.modform.maass.maass_forms import prediction,find_Y_and_M,get_M_and_Y
from psage.modform.maass.maass_forms_alg import get_coeff_fast_cplx_dp,get_coeff_fast_real_dp
from psage.modform.maass.maass_forms_alg import get_Y_for_M_dp_old as get_Y_for_M_dp

try:
    import colorlog
    from colorlog import ColoredFormatter
    LOG_LEVEL = logging.DEBUG
    logging.root.setLevel(LOG_LEVEL)
    LOGFORMAT = "  %(log_color)s%(levelname)-10s%(filename)s:%(lineno)d%(reset)s | %(log_color)s%(message)s%(reset)s"
    formatter = ColoredFormatter(LOGFORMAT)
    stream = logging.StreamHandler()
    stream.setLevel(LOG_LEVEL)
    stream.setFormatter(formatter)
    flogger = logging.getLogger(__name__)
    logger2 = logging.getLogger('detailed log')
    #if not flogger.handlers:
    flogger.addHandler(stream)
    if not logger2.handlers:
        logger2.addHandler(stream)    
except:
    logging.basicConfig(
        format='%(asctime)s %(levelname)s: %(message)s '
        '[in %(filename)s:-%(lineno)d]',
        datefmt='%Y%m%d-%H:%M%p',
    )
    flogger = logging.getLogger(__name__)
    logger2 = logging.getLogger('detailed log')
#flogger.setLevel(logging.WARNING)
logger2.setLevel(logging.DEBUG)


def get_element_in_range(H,R1,R2,Mset=0,Yset=0,dim=1,ndigs=12,set_c=[],fnr=0,neps=10,hmax=100,method='TwoY',split=0,get_c=0,verbose=0,maxit=20,force=False,db=None,tol=1e-12,hecke_ef=True,sloppy=0,return_bad=False):
    r""" Finds element of the space Hwith R in the interval R1 and R2

    INPUT:

    - ``R1`` -- lower bound (real)
    - ``R1`` -- upper bound (real)
    - ``method`` -- int - specify which method to use (to be described later...)
    - 'force' -- if set to False we skip intervals we already investigated
    - 'db' -- MaassDB

    """
    # Dummy element
    if hasattr(R1,'prec'):
        RF = RealField(R1.prec())
    else:
        RF = RealField(53)
    R1 = RF(R1); R2=RF(R2)
    F=Maasswaveform(H,R2,dim=dim,compute=False)
    e1 = RR(10)**-RR(neps)
    if verbose>0:
        print("R2,Mset,Yset,ndigs={0}".format(R2,Mset,Yset,ndigs))
        print("set_c={0}".format(set_c))
    param=H.set_default_parameters(R2,Mset,Yset,ndigs)
    Y=param['Y']
    Q=param['Q']
    M=param['M']
    if Y==None:
        M = H.smallest_M0()
        Y = get_Y_for_M_dp(H,R2,M,e1)
        if verbose>0:
            print("Recompute Y. Got:{0}".format(Y))
        Q = M+10
    Rl=list()

    if e1<tol and tol<0.5:
        neps=abs(int(log_b(tol,10)))
    if H._verbose==-2:
        P=Graphics()
    if H._verbose>0:
        print("Y={0}".format(Y))
        print("M={0}".format(M))
        print("Q={0}".format(Q))
        print("tol={0}".format(tol))
    if split==1:
        # Split interval on zeros of the K-Bessel function
        l=H.split_interval(R1,R2)
        if verbose>0:
            print("Split into intervals:")
            print("l={0}".format(l))
            for [r1,r2,y] in l:
                print("[{0},{1}]:{2}".format(r1,r2,y))
        for [r1,r2,y] in l:
            [R,er]=find_single_ev(H,r1,r2,Yset=y,neps=neps,method=method,dim=dim,verbose=verbose,maxit=maxit,set_c=set_c,fnr=fnr,hecke_ef=hecke_ef,Mset=Mset,sloppy=sloppy,tol_in=tol)
            if R!=0 and er!=0:
                Rl.append([R,er])
    else:
        r2=0
        y=Y*0.97
        if R1>0:
            r1=R1
        else:
            r1=0.0 # 1.0E-14
            ## And we also check separately if  R=0 is an eigenvalue.
        while r1<R2:
            # Where is the next end point?
            # This is good for large r, but for small?
            h=abs(1.0/H.Weyl_law_Np(r1))
            h=min(h,min(0.5,RR.pi()/r1))
            h=min(h,0.05)
            h=h*0.5
            r2=r1+h
            flogger.info("Using stepsize: h={0}".format(h))
            if db and have_eigenvalue_in(H,r1,r2,neps,db):
                flogger.info("Continuing! We already have an eigenvalue in the database in :[{0} , {1}], y={2}, neps={3}".format(r1,r2,y,neps))
                r1=r2
                continue
            flogger.info("checking :[{0} , {1}], y={2}, neps={3}".format(r1,r2,y,neps))
            try:
                l=find_single_ev(H,r1,r2,Yset=y,neps=neps,dim=dim,hmax=hmax,method=method,verbose=verbose,maxit=maxit,set_c=set_c,fnr=fnr,hecke_ef=hecke_ef,sloppy=sloppy,db=db,tol_in=tol)
            except KeyboardInterrupt:
                return
            flogger.info("l={0}".format(l))
            if H._verbose==-2:
                R,er,pp=l
                P+=pp
            else:
                try:
                    R,er,Y,M=l
                except (ValueError, TypeError):
                    R,er,Y,M=-1,0,0,0

            if R >=0 and er>0:  # Remember that a close but non-zero is often recognized by negative value of r
                #if db:
                Rl.append(l) #[H._group._level,H.sym_type(),R,er,Y,M])
            elif return_bad:
                Rl.append(l) ## We always add also bad results in the returned values.
            if r2==r1:
                r1 = r2+1E-14
            else:
                r1=r2
    if H._verbose==-2:
        return Rl,P
    else:
        return Rl


def find_single_ev(S,R1in,R2in,Yset=None,neps=10,method='TwoY',dim=1,verbose=None,set_c={},fnr=0,get_c=0,hmax=10000,maxit=20,db=None,hecke_ef=True,Mset=None,sloppy=0,tol_in=None):
    """

    Find a single eigenvalue in an interval if it exists.

    :param S:
    :param R1in:
    :param R2in:
    :param Yset:
    :param neps:
    :param method:
    :param dim:
    :param verbose:
    :param set_c:
    :param fnr:
    :param get_c:
    :param hmax:
    :param maxit:
    :param db:
    :param hecke_ef:
    :param Mset:
    :param sloppy:
    :param tol_in:
    :return:
    """
    i=0
    if verbose==None:
        verbose=S._verbose
    else:
        verbose=verbose
    if verbose>2:
        lvl0 = flogger.getEffectiveLevel()
        flogger.setLevel(int(10))
    else:
        flogger.setLevel(int(40))

    #print "verbose=",verbose
    #Mset=None
    mi=0
    R1=R1in; R3=R2in
    # We should only need to decrease Y a couple of times
    # We do it 5 times with the same M and then change M
    # and decrease Y 5 more times.
    # if this doesn't work we give up...
    #flogger.info("find single: set_c={0}".format(set_c))
    M0=0; Y0=0; message=''
    l = None
    while i < 10:
        try:
            l=find_single_ev_1(S,R1,R3,Yset=Yset,Mset=Mset,neps=neps,
                               method=method,dim=dim,hmax=hmax,
                               verbose=verbose,get_c=get_c,maxit=maxit,set_c=set_c,fnr=fnr,hecke_ef=hecke_ef,sloppy=sloppy,tol_in=tol_in)
        #except KeyError:
        #    pass
        except ArithmeticError as aritherr:
            message=str(aritherr)
            flogger.critical("Arithmetic error! {0}".format(message))
            if 'Need smaller Y' in message:
                flogger.critical("We continue and change parameters!")
                Y = float(message.split('than')[-1])
                l = -1,R1,R3,Y,Mset
        if l and isinstance(l,(list,tuple)):
            if l[0]==-1:
                i=i+1

                R1old=l[1]; R3old=l[2]
                Yset=l[3]*0.9
                M0=Mset
                Y0=Yset
                if mi<=5:
                    Mset=l[4]; mi+=1
                else:
                    Mset=None
                    mi=0
                if R1old>R1in and R1old<R2in:
                    R1=R1old
                if R3old>R1in and R3old<R2in:
                    R3=R3old
                flogger.debug("decreasing Y:{0}".format(Yset))
                flogger.debug("keeping M:{0}".format(Mset))
                flogger.debug("checking interval:[{0},{1}]".format(R1,R3))

            else:
                return l
        elif l != None:
            return l
        else:
            return 0,0,0,0
    if M0 is None:
        M0=0
    if not hasattr(S._group,'level'):
        flogger.info("Location failed at R1,R3,M0,Y0:{0}".format(R1,R3,M0,Y0))
        return 0,0,0,0
    data={'Level':int(S._group.level()),
          'Weight':int(S._weight),
          'Char':int(S._ch),
          'Dim':int(dim),
          'Symmetry':int(S._sym_type),
          'R1':float(R1),
          'R2':float(R3),
          'M0':int(M0),
          'Y':float(Y0),
          'Cusp_evs':S._cusp_evs #mongify(S._cusp_evs)
          }
    flogger.info("Location failed and we register a problem:{0}".format(data))
    if message=='':
        data['message']='Location routine failed!'
    else:
        data['message']='Location routine failed! '+message
    print(data['message'])



def find_single_ev_1(S,R1in,R3in,Yset=None,Mset=None,neps=10,method='TwoY',verbose=0,dim=1,hmax=100,get_c=0,set_c={},
                     fnr=0,maxit=20,hecke_ef=False,sloppy=0,tol_in=None,error_test_version='normal'):
    r""" Locate a single eigenvalue on G between R1 and R2

    INPUT:(tentative)

    .- ''S''    -- space of Maass waveforms
    - ''R1in'' -- real
    - ''R1in'' -- real
    - ''Yset'' -- real (use this value of Y to compute coefficients)
    - ''neps'' -- number of desired digits
    - ''method'' -- 'TwoY' or 'Hecke'
    - ''verbose'' -- integer (default 0) level of verbosity
    - ''dim''  -- integer (default 1) look for space of this dimension
    - ''hmax'' -- integer (default 15) parameter which sets how large 'too large' derivative is
    - ''get_c'' -- integer. (default 0) set to 1 if you want to return coefficients.
    - ''maxit'' -- integer (default 20,) number of iterations before returning with an error message
    - ''hecke_ef'' -- logical (default False) if set to True we check that for Hecke eigenforms.
                             (Note: this increases the necessary precision)

    OUTPUT:

    - ''R'' --


    """
    G=S.group()
    jmax=1000  # maximal number of iterations
    if neps>=15:
        prec = ceil(neps*3.5)
    else:
        prec = 53
    RF = RealField(prec)
    R1in = RF(R1in)
    R3in = RF(R3in)
    R1v = [RF(R1in)]; R3v = [RF(R3in)]
    #R1=mpmath.mp.mpf(R1in);R3=mpmath.mp.mpf(R2in)
    if S._verbose>0:
        flogger.setLevel(10)
        #print "mpmath.mp.dps=",mpmath.mp.dps
    flogger.debug("prec={0}".format(prec))
    flogger.debug("neps={0}".format(neps))    
    flogger.debug("Yset={0}".format(Yset))
    flogger.debug("R1={0} of type {1}".format(R1in,type(R1in)))
    flogger.debug("R3={0} of type {1}".format(R3in,type(R3in)))
    flogger.debug("Finding function nr. {0}".format(fnr))
    half=RF(0.5)
    if tol_in is None:
        tol=RF(2)**RF(10-prec)
    else:
        tol = tol_in
    tol = RF(tol)
    tol_z=tol.sqrt()  ## When we check for zeros we are more relaxed
    #[Y,M]=find_Y_and_M(G,R1,neps,Yset=Yset,Mset=Mset)
    M = S.smallest_M0()
    #Y = get_Y_for_M(S,R2in,M,tol)
    flogger.debug("R3in={0}, min hieight: {1}, M={2}, tol={3}".format(R3in,S.group().minimal_height(),M,tol))
    
    M,Y=get_M_and_Y(R3in,S.group().minimal_height(),M,tol,cuspidal=1)
    # start slightly lower?
    #if Y<0: # If the required precision needs a larger M
    #    Y,M = get_Y_and_M_dp(S,R2in,tol)
    if Y<0:
        raise ArithmeticError("Can not find Y for this precision!")
    if Yset:
        if Yset < Y and Yset>0:
            Y = Yset
    Y = RF(Y)
    if Mset and Mset>0:
        M=Mset
    flogger.debug("Got M={0} and Y={1}".format(M, Y))
    signs=dict();diffs=dict()
    c=dict(); h=dict()
    is_trivial=S.multiplier().is_trivial()
    st = S.sym_type()
    if not is_trivial:
        x=S.multiplier().character()
        mod=x.modulus()
    else:
        mod=S.group().generalised_level()
    flogger.debug("Computing Maass form witt these parameters: R={0}, dim={1} and set_C={2}".format(R1in,dim,set_c))
    F=S.get_element(R1in,Yset=Y,Mset=M,dim=dim,set_c=set_c,phase2=False)
    if isinstance(F,list):
        F = F[fnr]

    c = find_appropriate_coefficients_for_checks(F,method=method,tol=tol)
                                                  
    flogger.info("+++++++++++++++++++++++ Start find_single_ev_1 +++++++++++++++++++++++++++++++++++++++++++++")
    flogger.info("R1,R3={0},{1}".format(R1in,R3in))
    flogger.debug("tol={0}\t tol_z={1}".format(tol, tol_z))
    flogger.debug("Y={0} \t M={1}".format(Y,M))
    flogger.debug("level={0}".format(S._group._level)+"dim={0}".format(dim)+' sym_type={0}'.format(S._sym_type)+" ch={0}".format(S._ch))
    flogger.debug("cusp_evs={0}".format(S._cusp_evs))
    flogger.debug("c={0}".format(c))
    flogger.debug("method={0}".format(method))

    Y1=Y*RF(0.995); Y2=RF(0.995)*Y1
    met=method
    diffs[1],h[1] = functional(S,R1in,M,Y1,Y2,signs,c,first_time=True,method=met,dim=dim,set_c=set_c,fnr=fnr,verbose=verbose,hecke_ef=hecke_ef)
    # A quick option out in case we did not give an interval but two identical points.
    if R1in == R3in:
        errest = abs(h[1])
        if errest < tol:
            return R1in,errest,Y1,M
        else:
            return 0,0,0,Y,M
    [diffs[3],h[3]] = functional(S,R3in,M,Y1,Y2,signs,c,first_time=True,method=met,dim=dim,set_c=set_c,fnr=fnr,verbose=verbose,hecke_ef=hecke_ef)

    # Debugging output
    flogger.debug("Y1={0}".format(Y1)+" Y2={0}".format(Y2)+" M={0}".format(M))
    flogger.debug("R1={0},  R3={1}".format(R1in,R3in))
    for j in [1,2,3]:
        flogger.debug("diffs[{0}]={1}".format(j,diffs.get(j)))
    flogger.debug("h={0}".format(h))
        #flogger.debug("deriv. approx:1:{0}".format(t1)
        #flogger.debug("deriv. approx:2:{0}".format(t2)
        #flogger.debug("deriv. approx:3:{0}".format(t3)
    if verbose>1:
        for n in list(c.keys()): #.sort():
            for j in list(diffs.keys()): #.sort():
                if n in diffs[j]:
                    flogger.debug("diff[{0}][{1}]={2}".format(j,c[n],diffs[j][n]))

    ### Checking results
    if abs(h[1])>hmax or abs(h[3])>hmax:
        if verbose>0:
            flogger.debug("too large h1 or h3:{0}".format(h))
        return -1,R1in,R3in,Y,M
    if is_derivative_large(R1in,R3in,diffs[1],diffs[3],met,hmax,tol,verbose):
        return -1,R1in,R3in,Y,M
    # Do we want to add a plot
    add_plot=0
    if verbose==-2:
        add_plot=1
        P = Graphics()
        #fp=open("testplot.txt{0}".format("write")
        lopts={'thickness':1}
        l=line([(R1,h[1]),(R3,h[3])],**lopts)
        P+=l
    def return_no_zero():
        if add_plot == 1:
            return 0, 0, 0, 0, P
        else:
            return 0, 0, 0, Y, M
    # Set signs and check zeros
    if met=='Hecke' and dim==1:
        if h[1]*h[3]>sqrt(tol):
            return return_no_zero()
        if h[1]>0:
            signs[1]=1
        else:
            signs[1]=-1
    # Begin by trying to find the 'direction' of the crossing
    # This is a measure of how large the differences are
    var= sum([abs(diffs[1].get(j,0))+abs(diffs[3].get(j,0)) for j in diffs[1]])

    # In the first step we want to catch more more sign changes so we relax the conditions
    # Since the precision will in general not be fine to have sign changes in all
    # coefficients we are satisfied with two out of three
    nsgn = [ diffs[1 ][j]*diffs[3][j]<tol for j in diffs[1]].count(True)
    signs = {x: 1 for x in diffs[1]}
    if nsgn<1:
        # If we did not have any sign change then we probably don't have a zero
        # at least, if the relative size is large
        if sum([abs(diffs[1][j])+abs(diffs[3][j]) for j in diffs[1]]) > 0.04*var:
            flogger.debug("No zero here!")
            flogger.debug("var={0}".format(var))
            flogger.debug("diffs[{0}]={1} - {2}".format(j,abs(diffs[1][j]),abs(diffs[3][j])))
            return return_no_zero()

    flogger.debug("Number of sign changes={0}".format(nsgn))

    # If there is a crossing and it is positive to the left we set sign to -1 or
    # If there is no crossing but the value is larger we also set the sign to -1
    signs.update({j : -1  for j in diffs[1] if (diffs[1][j]*diffs[3][j] < tol and diffs[1][j] > 0) or
                            (diffs[1][j]*diffs[3][j] > tol and abs(diffs[1][j]) > abs(diffs[3][j])) })

    # Recompute functionals using the signs
    flogger.debug("Making preliminary search")
    flogger.debug(" R1,R3={0}".format((R1in,R3in)))
    flogger.debug(" h = {0}".format(h))
    for j in [1,2,3]:
        flogger.debug("diffs[{0}]={1}".format(j,diffs.get(j)))
    flogger.debug(" signs={0}".format(signs))
    for k in [1,3]:
        h[k]=0
        for j in range(1,3+1):
            if j not in diffs[k]:
                continue
            h[k]=h[k]+signs[j]*diffs[k][j]
    Rnew=prediction(h[1 ],h[3 ],R1in,R3in)
    Rnew_old=Rnew ## Keep track of the previous estimate
    flogger.debug(" new h={0}".format(h))
    flogger.debug(" Prediction 1={0}".format(prediction(diffs[1 ][1],diffs[3 ][1],R1v[0],R3v[0])))
    if 2 in diffs[1]:
        flogger.debug(" Prediction 2={0}".format(prediction(diffs[1 ][2],diffs[3 ][2],R1v[0],R3v[0])))
    if 3 in diffs[1]:
        flogger.debug(" Prediction 3={0}".format(prediction(diffs[1 ][3],diffs[3 ][3],R1v[0],R3v[0])))
    flogger.debug(" tol={0}".format(tol))
    flogger.debug(" Rnew={0}".format(Rnew))
    flogger.debug(" R3in+2tol={0}".format(R3in+2*tol))
    flogger.debug(" R1in-2tol={0}".format(R1in-2*tol))
    if Rnew > R3in+2*tol or Rnew<R1in-2*tol:
        ## Try to take an average prediction instead
        Rnew = sum([prediction(diffs[1 ][j],diffs[3 ][j],R1in,R3in) for j in range(1,4)])
        Rnew = Rnew/RF(3)
        flogger.debug("Try a new Rnew = {0}".format(Rnew))
        if Rnew > R3in+10*tol or Rnew<R1in-10*tol:
            return 0,0,0,Y,M

    [diffs[2],h[2]]=functional(S,Rnew,M,Y1,Y2,signs,c,first_time=False,method=met,dim=dim,set_c=set_c,fnr=fnr,verbose=verbose,hecke_ef=hecke_ef)
    flogger.debug("R1,Rnew,R3={0}".format((R1in,Rnew,R3in)))
    flogger.debug("h={0}".format(h))
    ## We use two different error estimates: errest_x and errest_h
    ## estimated from values of x and of h(x) respectively
    ### Checking results
    if abs(h[2])>hmax and sloppy>1:
        if verbose>0:
            flogger.debug("got too large h 1={0}".format(h))
        return -1,R1in,R3in,Y,M
    zero_in=is_zero_in(st,h,diffs,tol_z,verbose)
    if zero_in == -1:
        R3v.insert(0,Rnew); R1v.insert(0,R1v[0]); h[3]=h[2 ]; diffs[3 ]=diffs[2 ]; errest_x=abs(Rnew-R1v[0])
        if add_plot==1:
            P+=line([(R1,h[1]),(R3,h[3])],**lopts)
    else:
        R1v.insert(0,Rnew); R3v.insert(0,R3v[0]); h[1 ]=h[2 ]; diffs[1 ]=diffs[2 ]; errest_x=abs(Rnew-R3v[0])
        if add_plot==1:
            P+=line([(R1v[0],h[1]),(R3v[0],h[3])],**lopts)
    step=0
    if is_derivative_large(R1v[0],R3v[0],diffs[1],diffs[3],met,hmax,tol_z,verbose) and sloppy>1:
        return -1,R1v[0],R3v[0],Y,M
    errest_x_old=errest_x
    errest_h = sum(map(abs,h.values()))/3.0
    errest_h_old = errest_h
    #R1v=[R1]; R3v=[R3]
    for j in range(maxit):
        flogger.info(" ----------------- step: {0} -----------------------".format(j))
        Rnew_old=Rnew
        Rnew=prediction(h[1 ],h[3 ],R1v[0],R3v[0])
        [diffs[2],h[2]]=functional(S,Rnew,M,Y1,Y2,signs,c,first_time=False,method=met,dim=dim,set_c=set_c,fnr=fnr,verbose=verbose,hecke_ef=hecke_ef)
        errest_h_old = errest_h; errest_x_old = errest_x
        errest_h = sum(map(abs,h.values()))/3.0
        if method == 'Hecke' and -1 in diffs[1]:
            errest_h = max([errest_h,abs(diffs[1][-1]),abs(diffs[3][-1]),abs(diffs[2][-1])])
            flogger.debug("R1,R3v[0],Rnew,errest_h,tol={0}".format((R1v[0],R3v[0],Rnew,errest_h,tol)))
        ## The prediction tries to find zero but now we see how close we actually were:
        flogger.debug("Rnew ={0}".format(abs(Rnew)))
        Rnew1=prediction(h[1 ],h[2 ],R1v[0],Rnew)
        Rnew2=prediction(h[2 ],h[3 ],Rnew,R3v[0])            
        flogger.debug("Rnew_1 ={0}".format(abs(Rnew1)))
        flogger.debug("Rnew_2 ={0}".format(abs(Rnew2)))
        flogger.debug("Rnew_1-Rnew ={0}".format(abs(Rnew1-Rnew)))
        flogger.debug("Rnew_2-Rnew {0}=".format(abs(Rnew2-Rnew)))        
        errest_x = abs(R1v[0]-R3v[0]) ## Too strict?
        errest_x = max( abs(R1v[0]-Rnew),abs(R3v[0]-Rnew))
        flogger.debug("[R1,R3] = [{0},{1}".format(R1v[0],R3v[0]))
        flogger.debug("\t Rnew_old = {0}".format(Rnew_old))
        flogger.debug("\t Rnew_new(pred) = {0}".format(Rnew))
        flogger.debug("Rnew - Rnew_old".format(abs(Rnew-Rnew_old)))
        flogger.debug("\t h={0}".format(h))
        flogger.debug("\t R1[0],R3[0]={0}".format((R1v[0],R3v[0])))
        for ii in range(len(R1v)):
            flogger.debug("\t R1[{0}],R3[{0}]={1},{2}".format(ii,R1v[ii],R3v[ii]))        
        flogger.debug("Current errest_x={0}".format(errest_x))
        flogger.debug("Current errest_h={0}".format(errest_h))        
        flogger.debug("Old errest_x={0}".format(errest_x_old))
        flogger.debug("Old errest_h={0}".format(errest_h_old))
        if error_test_version == 'strict':
            errest = max(errest_x,errest_h)
        else:
            errest = errest_x

        for j in [1,2,3]:
            flogger.debug("diffs[{0}]={1}".format(j,diffs.get(j)))
        flogger.debug("h={0}".format(h))
        # First check if values are out of range and need recomputation with other parameters:
        if check_if_out_of_range(R1v, R3v, Rnew, diffs, h, errest_x, errest_x_old, errest_h, errest_h_old, met, hmax,
                                 tol, tol_z, verbose, sloppy):
            return -1, R1v[1], R3v[1], Y, M
        # If not we check if we are done:
        # Check if we are done
        flogger.debug("errest={0}".format(errest))
        flogger.debug("tol={0}".format(tol))

        if errest < tol:  ## Strict
            flogger.debug("Error estimate is {0} < tol= {1}".format(errest, tol))
            if add_plot == 1:
                return Rnew, errest, P
            else:
                return Rnew, errest, Y, M
        # Check if we landed outside previous interval
        # if Rnew > R3v[0] or Rnew<R1v[0]:
        #     flogger.debug("Rnew is outside the previous interval!")
        #     ## Try to take one of the other predictions instead
        #     #return -1,R1v[1],R3v[1],Y,M
        #     try:
        #         for j in range(1,4):
        #             Rnew = prediction(diffs[1 ][j],diffs[3 ][j],R1v[0],R3v[0])
        #             flogger.debug("Try new Rnew = Prediction({0}) = {1}".format(j,Rnew))
        #             if Rnew < R3v[0] and Rnew>R1v[0]:
        #                 raise StopIteration
        #         flogger.debug("We couldn't find any prediction inside the interval!")
        #         return -1,R1v[1],R3v[1],Y,M
        #     except StopIteration:
        #         pass
        ###
        ### Checks for errors in computation or parameter choices.
        ### We return signals that new parameters need to be used.
        ### 

        # # Check if we are done
        # if errest < tol:  ## Strict
        #     flogger.debug("Error estimate is {0} < tol= {1}".format(errest,tol))
        #     if add_plot==1:
        #         return Rnew,errest,P
        #     else:
        #         return Rnew,errest,Y,M
        #
        # We now see if we can find a zero in [R1,Rnew] or [R2,Rnew]      
        zero_in=is_zero_in(st,h,diffs,tol_z,verbose)

        flogger.debug("zero_in = {0}".format(zero_in))       
        if zero_in==0 and abs(errest)<tol:
            if add_plot==1:
                return Rnew,errest,P
            return Rnew,errest,Y,M
        elif zero_in not in [1,-1]:
            if verbose>0:
                flogger.debug("No zero! Breaking!")
            return Rnew,-errest,Y,M
            #flogger.debug("break!!!!"
            #    break #raise StopIteration()
        stepz={}
        if zero_in==-1:
            flogger.debug("There is a zero in [R1,Rnew]={0}".format((R1v[0],Rnew)))
            stepz[0]=abs(Rnew-R3v[0]); stepz[1]=abs(R1v[0]-R3v[0])            
            errest_x_old = errest_x
            if R1v[0] < Rnew:
                R3v.insert(0,Rnew); R1v.insert(0,R1v[0]); h[3 ]=h[2 ]; diffs[3 ]=diffs[2 ]
            errest_x=abs(Rnew-R1v[0]) # strict error estimate
            if hecke_ef:
                errest_h = max(errest_h,diffs[2][-1])

        elif zero_in==1:
            flogger.debug("There is a zero in [Rnew,R3]={0}".format((Rnew,R3v[0])))
            stepz[0]=abs(Rnew-R1v[0]); stepz[1]=abs(R1v[0]-R3v[0])                        
            errest_x_old = errest_x
            if R3v[0] > Rnew:
                R1v.insert(0,Rnew); R3v.insert(0,R3v[0]); h[1 ]=h[2 ]; diffs[1 ]=diffs[2 ]            
            errest_x=abs(Rnew-R3v[0]) # strict error estimate
            if hecke_ef:
                errest_h = max(errest_h,diffs[2][-1])
        # Check if the error estimate in x gets worse (shouldn't happen)
        if 2*errest_x_old < errest_x and errest_x_old<1e-10:
            flogger.warning("New error estimate in x {0} is larger than previous {1}. Improve parameters!".format(errest_x,errest_x_old))
            return -1,R1v[1],R3v[1],Y,M
        if add_plot==1:
            P+=line([(R1v[0],h[1]),(R3v[0],h[3])],**lopts)
        # If we have gone in the same direction too many times we need to modify our approach
        step=step+zero_in
        flogger.debug("step={0}".format(step))
        if step>2 or step < -2:    # Have changed the left end point too many times. Need to change the right also.
            ss = sign(step)
            if step > 2: 
                flogger.debug("Have moved the left end point too many times!")
            else:
                flogger.debug("Have moved the right end point too many times!")
            ## continue to decrease R3 as long as we detect a zero.
            try:
                for zj in range(2):
                    for jj in range(1,5):
                        Rtest=Rnew +ss*half**jj*stepz[zj]  # Need to test if this modified R3 work:
                        flogger.debug("diff={0}".format(ss*stepz[zj]))
                        if abs(stepz[zj]) < tol:
                            flogger.debug("Stepz ={0} < tol={1}".format(stepz[zj],tol))
                            raise StopIteration
                        if step > 2:
                            flogger.debug("Rtest1(R){0}={1}".format(jj,Rtest))
                        else:
                            flogger.debug("Rtest1(L){0}={1}".format(jj,Rtest))
                        [diffs[2 ],h[2 ]]=functional(S,Rtest,M,Y1,Y2,signs,c,False,method=met,dim=dim,set_c=set_c,fnr=fnr,verbose=verbose,hecke_ef=hecke_ef)
                        flogger.debug("h={0}".format(h))
                        for j in [1,2,3]:
                            flogger.debug("diffs[{0}]={1}".format(j,diffs.get(j)))
                        flogger.debug("R1,Rtest,R3={0}".format((R1v[0],Rtest,R3v[0])))
                        t1 = is_derivative_large(R1v[0],Rtest,diffs[1],diffs[2],met,hmax,tol_z,verbose)
                        t2 = is_derivative_large(R3v[0], Rtest, diffs[3], diffs[2], met, hmax, tol_z, verbose)
                        if t1 and t2 and sloppy>1:
                            return -1,R1v[1],R3v[1],Y,M
                        if abs(h[2])>hmax:
                            return -1,R1v[1],R3v[1],Y,M
                        t = is_zero_in(st,h,diffs,tol_z,verbose)
                        flogger.debug("zero in:{0}".format(t))
                        # if zj==2:
                        #     raise ValueError
                        if t*step > 0:
                            continue
                        if t == -1: # have zero in [R1v[0],Rtest]
                            R3v.insert(0,Rtest); R1v.insert(0,R1v[0]); h[3]=h[2]; diffs[3]=diffs[2]; step=step-2
                            flogger.debug("Set R3: R1,R3={0}".format((R1v[0],R3v[0])))
                        elif t==1: # have zero in [Rtest,R3v[0]]
                            R1v.insert(0,Rtest); R3v.insert(0,R3v[0]); h[1]=h[2]; diffs[1]=diffs[2]; step=step+2
                            flogger.debug("Set R1: R1,R3={0},{1}".format(R1v[0],R3v[0]))
                        else:
                            continue
                        # Check if we caught the zero or if we jumed too far
                        step = 0
                        if ss*t > 0:
                            flogger.debug("Get back to main loop!")
                            raise StopIteration
                        elif ss*t < 0:
                            flogger.debug("Get back to main loop with new (good) prediction!")
                            Rnew = Rtest
                            raise StopIteration
                        
                        # Also test if by now the difference at R3 is small enough
                        #i,er = is_minimum_at(h,diffs,tol_z,verbose)
                        #flogger.debug("i,er={0},{1}".format(i,er))
                        #if er<tol and errest_x < tol:
                        #    if i==1:
                        #        return R1v[0],er,Y,M
                        #    if i==2:
                        #        return Rtest,er,Y,M
                        #    if i==3:
                        #        return R3v[0],er,Y,M
                        if t == 1:
                            flogger.debug("Continue to move the left endpoint again!")
                        else:
                            flogger.debug("Continue to move the right endpoint again!")
                        
            except StopIteration:
                pass
            if add_plot==1:
                P+=line([(R1v[0],h[1]),(R3v[0],h[3])],**lopts)
            
    errest=min(errest_x,errest_h)
    if j>=maxit-1:
        flogger.info("Warning: may be problems in the locating function!")
        return Rnew,-errest,Y,M
    elif errest>0 and errest<1.0:
        return Rnew_old,-errest,Y,M        

def find_appropriate_coefficients_for_checks(F,method='TwoY',tol=1e-10,wanted_number_of_coefficients=3):
    """
    Find a set of coefficients that can be used when constructing a functional to search for Maass waveforms.
    :param F: a Maass waveform object (not necessarily a 'True' Maass form).
    :return:
    """
    set_c = F._set_c
    setc = set_c.get(0,{})
    if not setc:
        setc = {(0,1):1} # this is the default set coefficient
    usei = 1
    used_coefficients = []
    if method == 'Hecke':
        usei = 1
        a = F.space().get_primitive_p()
        b = F.space().get_primitive_p(a)
        c = a * b
        used_coefficients = [a, b, c]
    else:
        a = 1
    for (ci, n),v in setc.items():
        if v ==0: # Can't use a coefficient that is set to 0
            continue
        missing_number_of_coefficients = wanted_number_of_coefficients - len(used_coefficients)
        if not missing_number_of_coefficients:
            break
        this_test = lambda i: (0, i) not in setc and i not in used_coefficients and (i % n != 0 or n==1) \
                              and (F.C(0, i) is None or abs(F.C(0, i)) > tol ** 0.5) \
                              and gcd(i, F.generalised_level()) == 1 and gcd(i,a) == 1
        used_coefficients += [i for i in range(2,100) if this_test(i) ][0:missing_number_of_coefficients]
        # If we did not get sufficient number of coefficients from the set coefficients (or if there were no set)

    if len(used_coefficients) < wanted_number_of_coefficients:
        raise ArithmeticError("Could not get sufficiently many coefficients!")
    else:
        return { i+1: j for (i,j) in enumerate(used_coefficients)}

def get_coefficients_for_functional(S,r,Y,M,Q,dim=1,fnr=0,do_parallel=False,ncpus=1,set_c={},do_cplx=True,Norm=None,
                                    hecke_ef=False,used_coefficients=[],diffsx={},prec=53):
    r"""
    Get coefficients for using in the functional
    :param S:
    :param r:
    :param Y:
    :param M:
    :param dim:
    :param do_parallel:
    :param ncpus:
    :param do_cplx:
    :return:
    """
    logger2.info("Complex={0}".format(do_cplx))
    if not do_cplx:
        logger2.info("GEt coeff real fast: {0},{1},{2},{3},norm={4}".format(RR(r), RR(Y), int(M), int(Q), Norm))
        C1 = get_coeff_fast_real_dp(S, RR(r), RR(Y), int(M), int(Q), Norm)
        return C1[fnr][0]
    logger2.info("Get coefficients: r,Y1,M,Q={0}".format((r, Y, M, Q)))
    logger2.debug("Norm={0}".format(Norm))
    #logger2.debug("c = {0}, signs={1}".format(c, signs))
    if S.is_congruence():
        p = S.get_primitive_p()
        logger2.info("Hecke prime={0}".format(p) + " fnr={0}".format(fnr) + "dim={0}".format(dim))
        logger2.info("set_c={0}".format(set_c))
        # Ctmp=S.Hecke_eigenfunction_from_coeffs(C1,p)
        F = S.get_element(r, Y=Y, M0=M, dim=dim, do_par=do_parallel, ncpus=ncpus, set_c=set_c)
        if isinstance(F,list):
            logger2.info("len(F)={0}".format(len(F)))
            C1 =  F[0]._coeffs[0]
        else:
            C1 = F._coeffs[0]
        # If we are using Hecke tests then we add them also here
        if hecke_ef == 1:
            diffsx[-1] = S.test_Hecke_relation(C1, signed=True)
            logger2.info("Hecke C1(a)C1(b)-C1(ab)={0}".format(diffsx[-1]))
            logger2.info("C1.keys={0}".format(list(C1[0].keys())))
        if (used_coefficients[1] * used_coefficients[2]) in C1[0]:
            logger2.info("C1({0})={1}".format(used_coefficients[1] * used_coefficients[2], C1[0][used_coefficients[1] * used_coefficients[2]]))
        for j in used_coefficients.keys():
            logger2.debug("C1_Hecke0[0][{0}]={1}".format(used_coefficients[j], C1[0][used_coefficients[j]]))
        for j in used_coefficients.keys():
            sq = S.character()(used_coefficients[j])
            if sq != 1:
                sq = sq.complex_embedding(prec).sqrt()
            C1[0][used_coefficients[j]] = C1[0][used_coefficients[j]] / sq
        if p not in used_coefficients.values():
            logger2.debug("C1_Hecke[0][{0}]={1}".format(p, C1[0][p]))
        for j in used_coefficients.keys():
            logger2.debug("C1_Hecke[0][{0}]={1}".format(used_coefficients[j], C1[0][used_coefficients[j]]))
    else:
        C1 = get_coeff_fast_cplx_dp(S, r, Y, M, Q, Norm)
        C1 = C1[fnr]
    return C1

def functional(S,r,M,Y1,Y2,signs,c,first_time=False,method='TwoY',dim=1,set_c=[],fnr=0,do_parallel=False,ncpus=1,hecke_ef=1,verbose=0):
    r"""
    Computes the functional we use as an indicator of an eigenvalue.

    INPUT:

    -''S'' -- space of Maass waveforms
    -''r'' -- real
    -''M'' -- integer
    -''Y1''-- real
    -''Y2''-- real
    -''signs'' -- dict
    -''c'' -- set which coefficients to use
    -''ST''-- normalization
    -''first_time'' --
    -''method'' -- string in ['Hecke','TwoY','Innerprod']
    -''ndigs'' -- integer (number of digits wanted)

    OUTPUT:

    -  list of real values



    """
    diffsx=dict()
    h=0
    logger2.info("r,Y1,Y2={0}".format((r,Y1,Y2)))
    logger2.info("prec={0}".format(r.prec()))
    logger2.info("dim={0}".format(dim)+'\t method={0}'.format(method))
    logger2.info("set_c={0}".format(set_c))
    prec = r.prec()
    sym_type=S.sym_type()
    cusp_evs=S.atkin_lehner_eigenvalues() ## a dict
    Norm=S.set_norm(dim,set_c=set_c)
    if not S.group().is_congruence():
        hecke_ef = 0
    Q=M+20
    do_cplx=1
    if S.group().level()>1:
        do_cplx= True
    else:
        do_cplx = False
    if do_cplx:
        chi = S.multiplier()._character
        mod   = chi.modulus()
        sqch   = {}
        o = chi.order()
        for j in range(mod):
            if chi(j) != 1:
                sqch[j]=chi(j).complex_embedding(prec).sqrt()
            else:
                sqch[j]=1
    if method == 'Innerprod':
        raise NotImplementedError
    C1 = get_coefficients_for_functional(S,r,Y1,M,Q,dim=dim,fnr=fnr,do_parallel=do_parallel,ncpus=ncpus,set_c=set_c,do_cplx=do_cplx,
                                             Norm=Norm,hecke_ef=hecke_ef,used_coefficients=c,diffsx=diffsx)
    if method=='TwoY':
        diffsx2 = {}
        C2 = get_coefficients_for_functional(S, r, Y2, M, Q, dim=dim, fnr=fnr, do_parallel=do_parallel, ncpus=ncpus,
                                             set_c=set_c, do_cplx=do_cplx,
                                             Norm=Norm, hecke_ef=hecke_ef, used_coefficients=c,diffsx=diffsx2)
        if diffsx.get(-1):
            diffsx[-2]=diffsx[-1]

        # if do_cplx:
        #     logger2.info("do_cplx: r,Y2,M,Q={0}".format(r,Y2,M,Q,Norm,cusp_evs))
        #         #logger2.info("sqch={0}".format(sqch
        #     C2=get_coeff_fast_cplx_dp(S,RR(r),float(Y2),int(M),int(Q),Norm,cusp_ev=cusp_evs)
        #     C2 = C2[fnr]
        #     if dim>1:
        #         ## If we have a congruence group we try to do a Hecke eigenfunction.
        #         #if verbose>0:
        #         #    for k in range(dim):
        #         #        for j in c.keys():
        #         #            logger2.info("C2[{0}][{1}]={2}".format(k,c[j],C2[k][0][c[j]])
        #         if S.group().is_congruence():
        #             p = S.get_primitive_p()
        #             logger2.info("Hecke prime={0}".format(p)+" fnr={0}".format(fnr))
        #             F = S.get_element(r,Yset=Y2,Mset=M,dim=dim,do_par=do_parallel,ncpus=ncpus,set_c=set_c)
        #             logger2.info("len(F)={0}".format(len(F)))
        #             if len(F)<>dim:
        #                 raise ArithmeticError,"Error in dimension!"
        #             Ctmp2={}
        #             for j in range(len(F)):
        #                 Ctmp2[j]=F[j]._coeffs[0]
        #             C2=Ctmp2[fnr]
        #             #Ctmp=S.Hecke_eigenfunction_from_coeffs(C2,p)
        #             #C2=Ctmp[fnr]
        #             if hecke_ef==1:
        #                 diffsx[-2]=S.test_Hecke_relation(C2,signed=True)
        #                 logger2.info("Hecke C2(a)C2(b)-C2(ab)={0}".format(diffsx[-1]))
        #             for j in c.keys():
        #                 C2[0][c[j]]=C2[0][c[j]]/sqch[ c[j] % mod]
        #             if p not in c.values():
        #                 logger2.info("C2_Hecke[0][{0}]={1}".format(p,C2[0][p]))
        #             for j in c.keys():
        #                 logger2.info("C2_Hecke[0][{0}]={1}".format(c[j],C2[0][c[j]]))
        #         else:
        #             C2=C2[fnr]
        #
        #         #logger2.info("C2.keys(){0}".format(C2.keys()
        #         #logger2.info("C2[0].keys(){0}".format(C2[0].keys()
        #         #return
        #     else:
        #         for j in c.keys():
        #             logger2.debug("C2_Hecke0[0][{0}]={1}".format(c[j],C2[0][c[j]]))
        #         if C2[0].has_key(c[1]*c[2]):
        #             logger2.debug("C2({0})={1}".format(c[1]*c[2],C2[0][c[1]*c[2]]))
        #         for j in c.keys():
        #             C2[0][c[j]]=C2[0][c[j]]/sqch[ c[j] % mod]
        #         for j in c.keys():
        #             logger2.info("C2[0][{0}]={1}".format(c[j],C2[0][c[j]]))
        #
        #
        #         if hecke_ef==1:
        #             diffsx[-2]=S.test_Hecke_relation(C2,signed=True)
        #             logger2.info("Hecke C2(a)C2(b)-C2(ab)={0}".format(diffsx[-2]))
        #
        # else:
        #     C2=get_coeff_fast_real_dp(S,RR(r),RR(Y2),int(M),int(Q),Norm)
        #     #if dim>1:
        #     C2=C2[fnr]
        # #if dim>1:
        # #    C1=C1[fnr]
        # #s    C2=C2[fnr]
        # print "C1=",C1
        # print "C2=",C2
        if sym_type == -1:
            numt = len(c.keys())
            for j in c.keys():  #
                den = min(1,abs(C1[0][c[j]])+abs(C2[0][c[j]]))
                if isinstance(C1[0][c[j]],ComplexNumber):
                    diffsx[j]=(C1[0][c[j]]-C2[0][c[j]]).real()/den
                    diffsx[numt+j]=(C1[0][c[j]]-C2[0][c[j]]).imag()/den
                    logger2.debug("abs(C1[{0}])+abs(C2[{0}])={1}".format(c[j],den))
                    logger2.debug("diffs0={0}".format((C1[0][c[j]]-C2[0][c[j]]).real()))
                    logger2.debug( "diffs/den={0}".format(diffsx[j]))
                elif isinstance(C1[0][c[j]],complex):
                    diffsx[j]=(C1[0][c[j]]-C2[0][c[j]]).real/den
                    diffsx[numt+j]=(C1[0][c[j]]-C2[0][c[j]]).imag/den
                    logger2.debug("abs(C1[{0}])+abs(C2[{0}])={1}".format(c[j],den))
                    logger2.debug("diffs0={0}".format((C1[0][c[j]]-C2[0][c[j]]).real))
                    logger2.debug( "diffs/den={0}".format(diffsx[j]))
                else:
                    raise ArithmeticError("Got coefficients of unknown type {0}".format(type(C1[0][c[j]])))
        else:
            for j in c.keys():
                den = abs(C1[0][c[j]])+abs(C2[0][c[j]])
                if isinstance(C1[0][c[j]],ComplexNumber):
                   diffsx[j]=(C1[0][c[j]]-C2[0][c[j]]).real()/den
                elif isinstance(C1[0][c[j]],complex):
                    diffsx[j]=(C1[0][c[j]]-C2[0][c[j]]).real
                else:
                    raise ArithmeticError("Got coefficients of unknown type {0}".format(type(C1[0][c[j]])))
        h=0.0
        for j in [1,2,3]:
            logger2.debug("C1:2[{0}]={1} : {2}".format(c[j],C1[0][c[j]],C2[0][c[j]]))
            logger2.debug("diffsx[{1}]={0}".format(diffsx[j],j))
        logger2.debug("c={0}".format(c))
        logger2.debug("signs={0}".format(signs))
        for j in c.keys():
            #print "h0={0}".format(h
            if not first_time and list(signs.keys()).count(j)>0:
                #print "h={0}".format(h,"+{0}".format(signs[j]*diffsx[j]
                h=h+signs[j]*diffsx[j]
            else:
                #print "signs={0}".format(signs
                #print "diffsx={0}".format(diffsx
                if not signs.keys():
                    # We use the sign of the first test.
                    for k in c.keys():
                        if diffsx[1] != 0:
                            signs[k]=sign(diffsx[k])/sign(diffsx[1])
                        else:
                            signs[k]=sign(diffsx[k])
                        #print "sgns[{0}".format(k,"]={0}".format(signs[k]
                h=h+signs[j]*diffsx[j]
            #else:
            #    h=h+abs(diffsx[j])
            #print "h1={0}".format(h
        logger2.debug("diffsx={0}".format(diffsx))
        return [diffsx,h]
    elif method=='Hecke':
        if dim>1:
            ### We need to make a Hecke eigenform if possible
            #raise NotImplementedError,"Use TwoY for dimension>1!"
            for j in range(min(3,dim)):
                if S._sym_type==-1:
                    diffsx[j+1]=abs(real(Ctmp[j][0][-1]))-1
                    logger2.info("C[{0}][0](-1)={1}".format(j,Ctmp[j][0][-1]))
                elif cusp_evs[j]==0:
                    diffsx[j+1]=abs(real(Ctmp[j][1][1]))-1
                    logger2.info("C[{0}][1](1)={1}".format(Ctmp[j][1][1]))
                else:
                    ## Try the Hecke relation
                    diffsx[j+1]=diffsx[-1]
            if dim==2:
                if S._character.order()==222:
                    N = S._group._level
                    diffsx[3]=abs(Ctmp[0][0][N])-1
                    logger2.info("C[0][0]({0})={1}".format(N,Ctmp[j][0][N]))
                else:
                    # Try another Hecke relation
                    a = S.get_primitive_p(p)
                    b = S.get_primitive_p(a)
                    diffsx[3]=S.test_Hecke_relation(Ctmp[0],a,b,signed=True)
                    logger2.info("C[0][0]({0})={1}".format(a,Ctmp[0][0][a]))
                    logger2.info("C[0][0]({0})={1}".format(b,Ctmp[0][0][b]))
                    logger2.info("C[0][0]({0})={1}".format(a*b,Ctmp[0][0][a*b]))
            logger2.info("diffsx({0})={1}".format(r,diffsx))
            h=0
            for j in range(1,4):
                if not first_time and list(signs.keys()).count(j)>0:
                    htmp=signs[j]*diffsx[j]
                else:
                    if not signs.keys():
                        # We use the sign of the first test.
                        for k in c.keys():
                            if diffsx[1] != 0:
                                signs[k]=sign(diffsx[k])/sign(diffsx[1])
                            else:
                                signs[k]=sign(diffsx[k])
                        #logger2.debug("sgns[{0}".format(k,"]={0}".format(signs[k]
                    htmp=signs[j]*diffsx[j]
                #logger2.debug("h({0}".format(r,j,")={0}".format(h
                h = h+htmp
            logger2.info("h({0})={1}".format(r,h))
        elif dim==1:
            difftmp=C1[0][c[1]]*C1[0][c[2]]-C1[0][c[3]]*S._multiplier._character(c[3])
            if isinstance(difftmp,(float,complex)):
                diffsx[1]=difftmp.real
            else:
                diffsx[1]=difftmp.real()
                logger2.info("C[{0}]={1}".format(c[1],C1[0][c[1]]))
                logger2.info("C[{0}]={1}".format(c[2],C1[0][c[2]]))
                logger2.info("C[{0}]={1}".format(c[3],C1[0][c[3]]))
                logger2.info("difftmp={0}".format(difftmp))
                logger2.info("diffsx={0}".format(diffsx))
            if not first_time and list(signs.keys()).count(1)>0:
                h=signs[1]*diffsx[1]
            else:
                h = diffsx[1]
        return [diffsx,h]

    elif method=='Innerprod':
        raise NotImplementedError

def check_if_out_of_range(R1v,R3v,Rnew,diffs,h,errest_x,errest_x_old,errest_h,errest_h_old,met,hmax,tol,tol_z,
                          verbose,sloppy):
    """
    Check if computed values are 'too large' and we need to recompute.

    :param R1v:
    :param R3v:
    :param Rnew:
    :param diffs:
    :param errest_x:
    :param errest_x_old:
    :param errest_h:
    :param errest_h_old:
    :param met:
    :param hmac:
    :param tol:
    :param tol_z:
    :param hmax:
    :param verbose:
    :param sloppy:
    :return:
    """
    if Rnew > R3v[0] or Rnew < R1v[0]:
        flogger.debug("Rnew is outside the previous interval!")
        return True
    if is_derivative_large(R1v[0],Rnew,diffs[1],diffs[2],met,hmax,tol_z,verbose) and sloppy >1:
        flogger.debug("The derivative in [R1,Rnew] is too large!")
        return True
    if is_derivative_large(Rnew,R3v[0],diffs[2],diffs[3],met,hmax,tol_z,verbose) and sloppy >1:
        flogger.debug("The derivative in [Rnew,R3] is too large!")
        return True
    ## See if the functional value is too large
    if abs(h[2])>hmax:
        return True
    # Check if we got larger error estimate in x or h (shouldn't happen)
    if errest_x > errest_x_old: # or errest_h > errest_h_old:
        if verbose>0:
            flogger.debug("Got larger error estimate: {0} > {1}".format(errest_x,errest_x_old))
#            flogger.debug("Got larger error estimate: {0} > {1}".format(errest_h,errest_h_old))
        return True
    # Check if the R doesn't move but we have too large functional values.
    # Normally the function is much smaller.
    if errest_x < tol and errest_h > 10*tol:
        return True
    return False

def find_single_ev_noncong(S,R1in,R3in,Yset=None,Mset=None,dim=1,tol=1e-7,neps=10,method='TwoY',verbose=0,hmax=100,get_c=0,set_c=[],fnr=0,maxit=20,hecke_ef=False,sloppy=0):
    r""" Locate a single eigenvalue on G between R1 and R2

    INPUT:(tentative)

    .- ''S''    -- space of Maass waveforms
    - ''R1in'' -- real
    - ''R1in'' -- real
    - ''Yset'' -- real (use this value of Y to compute coefficients)
    - ''neps'' -- number of desired digits
    - ''method'' -- 'TwoY' or 'Hecke'
    - ''verbose'' -- integer (default 0) level of verbosity
    - ''dim''  -- integer (default 1) look for space of this dimension
    - ''hmax'' -- integer (default 15) parameter which sets how large 'too large' derivative is
    - ''get_c'' -- integer. (default 0) set to 1 if you want to return coefficients.
    - ''maxit'' -- integer (default 20,) number of iterations before returning with an error message
    - ''hecke_ef'' -- logical (default False) if set to True we check that for Hecke eigenforms.
                             (Note: this increases the necessary precision)

    OUTPUT:

    - ''R'' --


    """
    G=S.group()
    jmax=1000  # maximal number of iterations
    if hasattr(R1in,'prec'):
        prec = R1in.prec()
    else:
        prec = 53
    RF = RealField(prec)
    R1in = RF(R1in)
    R3in = RF(R3in)
    R1v = [RF(R1in)]; R3v = [RF(R3in)]
    #R1=mpmath.mp.mpf(R1in);R3=mpmath.mp.mpf(R2in)
    if verbose>0:
        flogger.setLevel(10)
        #print "mpmath.mp.dps=",mpmath.mp.dps
    flogger.debug("prec={0}".format(prec))
    flogger.debug("Yset={0}".format(Yset))
    flogger.debug("R1={0} of type {1}".format(R1in,type(R1in)))
    flogger.debug("R3={0} of type {1}".format(R3in,type(R3in)))
    flogger.debug("Finding function nr. {0}".format(fnr))
    half=RF(0.5)
    tol = RF(tol)
    tol_z=tol.sqrt()  ## When we check for zeros we are more relaxed
    M = S.smallest_M0()
    flogger.debug("R3in={0}, min hieight: {1}, M={2}, tol={3}".format(R3in,S.group().minimal_height(),M,tol))
    M,Y=get_M_and_Y(R3in,S.group().minimal_height(),M,tol)
    flogger.debug("Got M={0} and Y={1}".format(M,Y))
    if Y<0:
        raise ArithmeticError("Can not find Y for this precision!")
    if Yset!=None:
        if Yset < Y and Yset>0:
            Y = Yset
    Y = RF(Y)
    if Mset and Mset>0:
        M=Mset
    signs=dict();diffs=dict()
    c=dict(); h=dict()
    mod=S.group().generalised_level()
    flogger.debug("Computing Maass form with these parameters: R={0}, dim={1} and set_C={2}".format(R1in,dim,set_c))
    F=S.get_element(R1in,Yset=Y,Mset=M,dim=dim,set_c=set_c,phase2=False)
    if isinstance(F,list):
        F = F[fnr]
    if method=='TwoY':
        #if dim==1:
        #    c[1]=2 ; c[2]=3 ; c[3]=4
        #else:
        usei=1
        if len(set_c)>fnr:
            setc = set_c[fnr]
            for (ci,n) in setc:
                if setc[(ci,n)] != 0:
                    for i in range(2,100):
                        if (0,i) in setc:
                            continue
                        if i % n == 0:  # We might want to use this                    
                            if abs(F.C(0,i))>tol**0.5:
                                if verbose>0:
                                    flogger.debug("Using these coeffs: F.C(0,{0})={1}".format(i,F.C(0,i)))
                                c[usei] = i
                                usei+=1
                        if usei>3:
                            break
                if usei>3:
                    break
                
        if usei < 3:
            for i in range(dim+1,100):
                if gcd(i,mod)>1:
                    continue
                cont = 0
                if len(set_c)>fnr:
                    setc = set_c[fnr]
                    if (0,i) in setc:
                        cont = 1                        
                if i>M-2:
                    cont = 1
                if cont:
                    continue
                c[usei]=i
                usei+=1
                if usei>3:
                    break
                                                  
    flogger.info("+++++++++++++++++++++++ Start find_single_ev_1 +++++++++++++++++++++++++++++++++++++++++++++")
    flogger.info("R1,R3={0},{1}".format(R1in,R3in))
    flogger.debug("tol={0}\t tol_z={1}".format(tol, tol_z))
    flogger.debug("Y={0} \t M={1}".format(Y,M))
    flogger.debug("c={0}".format(c))

    Y1=Y*RF(0.995); Y2=RF(0.995)*Y1
    flogger.debug("for diffs1: R1in={0},M={1}, Y1={2}, Y2={3}".format(R1in,M,Y1,Y2))
    met=method
    [diffs[1 ],h[1 ]]=functional_noncong(S,R1in,M,Y1,Y2,signs,c,first_time=True,method=met,dim=dim,set_c=set_c,fnr=fnr,verbose=verbose,hecke_ef=hecke_ef)
    if R1in == R3in:
        errest = abs(h[1])
        if errest < tol:
            return R1in,errest,Y1,M
        else:
            return 0,0,0,Y,M
    [diffs[3 ],h[3 ]]=functional_noncong(S,R3in,M,Y1,Y2,signs,c,first_time=True,method=met,dim=dim,set_c=set_c,fnr=fnr,verbose=verbose,hecke_ef=hecke_ef)
    for j in [1,2,3]:
        flogger.debug("diffs[{0}]={1}".format(j,diffs.get(j)))
    flogger.debug("Y1={0}".format(Y1)+" Y2={0}".format(Y2)+" M={0}".format(M))
    flogger.debug("h={0}".format(h))
    flogger.debug("R1={0},  R3={1}".format(R1in,R3in))
    if verbose>1:
        for n in list(c.keys()): #.sort():
            for j in list(diffs.keys()): #.sort():
                if n in diffs[j]:
                    flogger.debug("diff[{0}][{1}]={2}".format(j,c[n],diffs[j][n]))
    ### Checking results
    if abs(h[1])>hmax or abs(h[3])>hmax:
        if verbose>0:
            flogger.debug("too large h1 or h3:{0}".format(h))
        return -1,R1in,R3in,Y,M
    if is_derivative_large(R1in,R3in,diffs[1],diffs[3],met,hmax,tol,verbose):
        return -1,R1in,R3in,Y,M
    add_plot=0
    if verbose==-2:
        add_plot=1
        P = Graphics()
        #fp=open("testplot.txt{0}".format("write")
        lopts={'thickness':1}
        l=line([(R1,h[1]),(R3,h[3])],**lopts)
        P+=l
    # Begin by trying to find the 'direction' of the crossing
    var=0.0
    for j in range(1 ,3 +1 ):
        var+=abs(diffs[1 ].get(j,0))+abs(diffs[3 ].get(j,0))
    # In the first step we want to catch more more sign changes so we relax the conditions
    # Since the precision will in general not be fine to have sign changes in all
    # coefficients we are satisfied with two out of three
    nsgn=0
    for j in range(1 ,3 +1 ):
        signs[j]=1
        if j not in diffs[1]:
            continue
        if diffs[1 ][j]*diffs[3][j]<tol:
            nsgn+=1
    if nsgn<1:
        # If we did not have any sign change then we probably don't have a zero
        # at least, if the relative size is large
        if sum([abs(diffs[1][j])+abs(diffs[3][j]) for j in range(1,4)]) > 0.04*var:
            flogger.debug("No zero here!")
            flogger.debug("var={0}".format(var))
            flogger.debug("diffs[{0}]={1} - {2}".format(j,abs(diffs[1][j]),abs(diffs[3][j])))
            if add_plot==1:
                return 0,0,0 ,0,P
            else:
                return 0,0 ,0,Y,M
    flogger.debug("Number of sign changes={0}".format(nsgn))

    for j in range(1,4):
        if j not in diffs[1]:
            continue
        if diffs[1 ][j]*diffs[3 ][j]<tol:
            if diffs[1 ][j]>0:
                signs[j]=-1
            else:
                signs[j]=1
        else:
            if abs(diffs[1][j])>abs(diffs[3][j]):
                signs[j]=-1
            else:
                signs[j]=1
    # Recompute functionals using the signs
    flogger.debug("Making preliminary search")
    flogger.debug(" R1,R3={0}".format((R1in,R3in)))
    flogger.debug(" h = {0}".format(h))
    for j in [1,2,3]:
        flogger.debug("diffs[{0}]={1}".format(j,diffs.get(j)))
    flogger.debug(" signs={0}".format(signs))
    for k in [1,3]:
        h[k]=0
        for j in range(1,3+1):
            if j not in diffs[k]:
                continue
            h[k]=h[k]+signs[j]*diffs[k][j]
    Rnew=prediction(h[1 ],h[3 ],R1in,R3in)
    Rnew_old=Rnew ## Keep track of the previous estimate
    flogger.debug(" new h={0}".format(h))
    flogger.debug(" Prediction 1={0}".format(prediction(diffs[1 ][1],diffs[3 ][1],R1v[0],R3v[0])))
    if 2 in diffs[1]:
        flogger.debug(" Prediction 2={0}".format(prediction(diffs[1 ][2],diffs[3 ][2],R1v[0],R3v[0])))
    if 3 in diffs[1]:
        flogger.debug(" Prediction 3={0}".format(prediction(diffs[1 ][3],diffs[3 ][3],R1v[0],R3v[0])))
    flogger.debug(" tol={0}".format(tol))
    flogger.debug(" Rnew={0}".format(Rnew))
    flogger.debug(" R3in+2tol={0}".format(R3in+2*tol))
    flogger.debug(" R1in-2tol={0}".format(R1in-2*tol))
    if Rnew > R3in+2*tol or Rnew<R1in-2*tol:
        ## Try to take an average prediction instead
        Rnew = sum([prediction(diffs[1][j],diffs[3][j],R1in,R3in) for j in range(1,4) if diffs[1][j] != diffs[3][j]])
        Rnew = Rnew/RF(3)
        flogger.debug("Try a new Rnew = {0}".format(Rnew))
        if Rnew > R3in+10*tol or Rnew<R1in-10*tol:
            return 0,0,0,Y,M

    [diffs[2],h[2]]=functional_noncong(S,Rnew,M,Y1,Y2,signs,c,first_time=False,method=met,dim=dim,set_c=set_c,fnr=fnr,verbose=verbose,hecke_ef=hecke_ef)
    flogger.debug("R1,Rnew,R3={0}".format((R1in,Rnew,R3in)))
    flogger.debug("h={0}".format(h))
    ## We use two different error estimates: errest_x and errest_h
    ## estimated from values of x and of h(x) respectively
    ### Checking results
    st = -1
    if abs(h[2])>hmax and sloppy>1:
        if verbose>0:
            flogger.debug("got too large h 1={0}".format(h))
        return -1,R1in,R3in,Y,M
    zero_in=is_zero_in(st,h,diffs,tol_z,verbose)
    if zero_in == -1:
        R3v.insert(0,Rnew); R1v.insert(0,R1v[0]); h[3]=h[2 ]; diffs[3 ]=diffs[2 ]; errest_x=abs(Rnew-R1v[0])
        if add_plot==1:
            P+=line([(R1,h[1]),(R3,h[3])],**lopts)
    else:
        R1v.insert(0,Rnew); R3v.insert(0,R3v[0]); h[1 ]=h[2 ]; diffs[1 ]=diffs[2 ]; errest_x=abs(Rnew-R3v[0])
        if add_plot==1:
            P+=line([(R1v[0],h[1]),(R3v[0],h[3])],**lopts)
    step=0
    if is_derivative_large(R1v[0],R3v[0],diffs[1],diffs[3],met,hmax,tol_z,verbose) and sloppy>1:
        return -1,R1v[0],R3v[0],Y,M
    errest_x_old=errest_x
    errest_h = sum(map(abs,h.values()))/3.0
    errest_h_old = errest_h
    #R1v=[R1]; R3v=[R3]
    for j in range(maxit):
        flogger.info(" ----------------- step: {0} -----------------------".format(j))
        Rnew_old=Rnew
        Rnew=prediction(h[1 ],h[3 ],R1v[0],R3v[0])
        [diffs[2 ],h[2 ]]=functional_noncong(S,Rnew,M,Y1,Y2,signs,c,first_time=False,method=met,dim=dim,set_c=set_c,fnr=fnr,verbose=verbose,hecke_ef=hecke_ef)
        errest_h_old = errest_h; errest_x_old = errest_x
        errest_h = sum(map(abs,h.values()))/3.0
        if method=='Hecke' and -1 in diffs[1]:
            errest_h = max([errest_h,abs(diffs[1][-1]),abs(diffs[3][-1]),abs(diffs[2][-1])])
            flogger.debug("R1,R3v[0],Rnew,errest_h,tol={0}".format((R1v[0],R3v[0],Rnew,errest_h,tol)))

        ## The prediction tries to find zero but now we see how close we actually were:
        Rnew1=prediction(h[1 ],h[2 ],R1v[0],Rnew)
        Rnew2=prediction(h[2 ],h[3 ],Rnew,R3v[0])            
        flogger.debug("Rnew_1 ={0}".format(abs(Rnew1)))
        flogger.debug("Rnew_2 ={0}".format(abs(Rnew2)))
        flogger.debug("Rnew_1-Rnew ={0}".format(abs(Rnew1-Rnew)))
        flogger.debug("Rnew_2-Rnew {0}=".format(abs(Rnew2-Rnew)))        
        errest_x = abs(R1v[0]-R3v[0]) ## Too strict?
        errest_x = max( abs(Rnew1-Rnew),abs(Rnew2-Rnew))
        flogger.debug("[R1,R3] = [{0},{1}".format(R1v[0],R3v[0]))
        flogger.debug("\t Rnew_old = {0}".format(Rnew_old))
        flogger.debug("\t Rnew_new(pred) = {0}".format(Rnew))
        flogger.debug("Rnew - Rnew_old".format(abs(Rnew-Rnew_old)))
        flogger.debug("\t h={0}".format(h))
        flogger.debug("\t R1[0],R3[0]={0}".format((R1v[0],R3v[0])))
        for ii in range(len(R1v)):
            flogger.debug("\t R1[{0}],R3[{0}]={1},{2}".format(ii,R1v[ii],R3v[ii]))        
        flogger.debug("Current errest_x={0}".format(errest_x))
        flogger.debug("Current errest_h={0}".format(errest_h))        
        flogger.debug("Old errest_x={0}".format(errest_x_old))
        flogger.debug("Old errest_h={0}".format(errest_h_old))
        errest = max(errest_x,errest_h)
        
        for j in [1,2,3]:
            flogger.debug("diffs[{0}]={1}".format(j,diffs.get(j)))
        flogger.debug("h={0}".format(h))
          
        # Check if we landed outside previous interval
        if Rnew > R3v[0] or Rnew<R1v[0]:
            flogger.debug("Rnew is outside the previous interval!")
            ## Try to take one of the other predictions instead
            #return -1,R1v[1],R3v[1],Y,M
            try: 
                for j in range(1,4):
                    Rnew = prediction(diffs[1 ][j],diffs[3 ][j],R1v[0],R3v[0])
                    flogger.debug("Try new Rnew = Prediction({0}) = {1}".format(j,Rnew))
                    if Rnew < R3v[0] and Rnew>R1v[0]:
                        raise StopIteration
                flogger.debug("We couldn't find any prediction inside the interval!")
                return -1,R1v[1],R3v[1],Y,M
            except StopIteration:
                pass
        ###
        ### Checks for errors in computation or parameter choices.
        ### We return signals that new parameters need to be used.
        ### 
        ## See if the derivative has been too large:
        if is_derivative_large(R1v[0],Rnew,diffs[1],diffs[2],met,hmax,tol_z,verbose) and sloppy >1:
            flogger.debug("The derivative in [R1,Rnew] is too large!")
            return -1,R1v[1],R3v[1],Y,M
        if is_derivative_large(Rnew,R3v[0],diffs[2],diffs[3],met,hmax,tol_z,verbose) and sloppy >1:
            flogger.debug("The derivative in [Rnew,R3] is too large!")
            return -1,R1v[1],R3v[1],Y,M
        ## See if the functional value is too large
        if abs(h[2])>hmax:
            return -1,R1v[1],R3v[1],Y,M
        # Check if we got larger error estimate in x or h (shouldn't happen)
        if errest_x > errest_x_old or errest_h > errest_h_old:
            if verbose>0:
                flogger.debug("Got larger error estimate: {0} > {1}".format(errest_x,errest_x_old))
                flogger.debug("Got larger error estimate: {0} > {1}".format(errest_h,errest_h_old))
            return -1,R1v[1],R3v[1],Y,M
        # Check if the R doesn't move but we have too large functional values. 
        # Normally the function is much smaller.
        if errest_x < tol and errest_h > 10*tol:
            return -1,R1v[1],R3v[1],Y,M

        # Check if we are done
        if errest < tol:  ## Strict
            flogger.debug("Error estimate is {0} < tol= {1}".format(errest,tol))
            if add_plot==1:
                return Rnew,errest,P
            else:
                return Rnew,errest,Y,M
      
        # We now see if we can find a zero in [R1,Rnew] or [R2,Rnew]      
        zero_in=is_zero_in(st,h,diffs,tol_z,verbose)

        flogger.debug("zero_in = {0}".format(zero_in))       
        if zero_in==0 and abs(errest)<tol:
            if add_plot==1:
                return Rnew,errest,P
            return Rnew,errest,Y,M
        elif zero_in not in [1,-1]:
            if verbose>0:
                flogger.debug("No zero! Breaking!")
            return Rnew,-errest,Y,M
            #flogger.debug("break!!!!"
            #    break #raise StopIteration()
        stepz={}
        if zero_in==-1:
            flogger.debug("There is a zero in [R1,Rnew]={0}".format((R1v[0],Rnew)))
            stepz[0]=abs(Rnew-R3v[0]); stepz[1]=abs(R1v[0]-R3v[0])            
            errest_x_old = errest_x
            if R1v[0] < Rnew:
                R3v.insert(0,Rnew); R1v.insert(0,R1v[0]); h[3 ]=h[2 ]; diffs[3 ]=diffs[2 ]
            errest_x=abs(Rnew-R1v[0]) # strict error estimate
            if hecke_ef:
                errest_h = max(errest_h,diffs[2][-1])

        elif zero_in==1:
            flogger.debug("There is a zero in [Rnew,R3]={0}".format((Rnew,R3v[0])))
            stepz[0]=abs(Rnew-R1v[0]); stepz[1]=abs(R1v[0]-R3v[0])                        
            errest_x_old = errest_x
            if R3v[0] > Rnew:
                R1v.insert(0,Rnew); R3v.insert(0,R3v[0]); h[1 ]=h[2 ]; diffs[1 ]=diffs[2 ]            
            errest_x=abs(Rnew-R3v[0]) # strict error estimate
            if hecke_ef:
                errest_h = max(errest_h,diffs[2][-1])
        # Check if the error estimate in x gets worse (shouldn't happen)
        if 2*errest_x_old < errest_x and errest_x_old<1e-10:
            flogger.warning("New error estimate in x {0} is larger than previous {1}. Improve parameters!".format(errest_x,errest_x_old))
            return -1,R1v[1],R3v[1],Y,M
        if add_plot==1:
            P+=line([(R1v[0],h[1]),(R3v[0],h[3])],**lopts)
        # If we have gone in the same direction too many times we need to modify our approach
        step=step+zero_in
        flogger.debug("step={0}".format(step))
        if step>2 or step < -2:    # Have changed the left end point too many times. Need to change the right also.
            ss = sign(step)
            if step > 2: 
                flogger.debug("Have moved the left end point too many times!")
            else:
                flogger.debug("Have moved the right end point too many times!")
            ## continue to decrease R3 as long as we detect a zero.
            try:
                for zj in range(2):
                    for jj in range(1,100):
                        Rtest=Rnew +ss*half**jj*stepz[zj]  # Need to test if this modified R3 work:
                        if step > 2:
                            flogger.debug("Rtest1(R){0}={1}".format(jj,Rtest))
                        else:
                            flogger.debug("Rtest1(L){0}={1}".format(jj,Rtest))
                        [diffs[2 ],h[2 ]]=functional_noncong(S,Rtest,M,Y1,Y2,signs,c,False,method=met,dim=dim,set_c=set_c,fnr=fnr,verbose=verbose,hecke_ef=hecke_ef)
                        flogger.debug("h={0}".format(h))
                        for j in [1,2,3]:
                            flogger.debug("diffs[{0}]={1}".format(j,diffs.get(j)))
                        flogger.debug("R1,Rtest,R3={0}".format((R1v[0],Rtest,R3v[0])))
                        if is_derivative_large(R1v[0],Rtest,diffs[1],diffs[2],met,hmax,tol_z,verbose) and sloppy>1:
                            return -1,R1v[1],R3v[1],Y,M
                        if abs(h[2])>hmax:
                            return -1,R1v[1],R3v[1],Y,M
                        t = is_zero_in(st,h,diffs,tol_z,verbose)
                        flogger.debug("zero in:{0}".format(t))
                        if t == -1: # have zero in [R1v[0],Rtest]
                            R3v.insert(0,Rtest); R1v.insert(0,R1v[0]); h[3]=h[2]; diffs[3]=diffs[2]; step=step-2
                            flogger.debug("Set R3: R1,R3={0}".format((R1v[0],R3v[0])))
                        elif t==1: # have zero in [Rtest,R3v[0]]
                            R1v.insert(0,Rtest); R3v.insert(0,R3v[0]); h[1]=h[2]; diffs[1]=diffs[2]; step=step+2
                            flogger.debug("Set R1: R1,R3={0},{1}".format(R1v[0],R3v[0]))
                        else:
                            continue
                        # Check if we caught the zero or if we jumed too far
                        step = 0
                        if ss*t > 0:
                            flogger.debug("Get back to main loop!")
                            raise StopIteration
                        elif ss*t < 0:
                            flogger.debug("Get back to main loop with new (good) prediction!")
                            Rnew = Rtest
                            raise StopIteration
                        
                        # Also test if by now the difference at R3 is small enough
                        #i,er = is_minimum_at(h,diffs,tol_z,verbose)
                        #flogger.debug("i,er={0},{1}".format(i,er))
                        #if er<tol and errest_x < tol:
                        #    if i==1:
                        #        return R1v[0],er,Y,M
                        #    if i==2:
                        #        return Rtest,er,Y,M
                        #    if i==3:
                        #        return R3v[0],er,Y,M
                        if t == 1:
                            flogger.debug("Continue to move the left endpoint again!")
                        else:
                            flogger.debug("Continue to move the right endpoint again!")
                        
            except StopIteration:
                pass
            if add_plot==1:
                P+=line([(R1v[0],h[1]),(R3v[0],h[3])],**lopts)
            
    errest=min(errest_x,errest_h)
    if j>=maxit-1:
        flogger.info("Warning: may be problems in the locating function!")
        return Rnew,-errest,Y,M
    elif errest>0 and errest<1.0:
        return Rnew_old,-errest,Y,M     

def functional_noncong(S,r,M,Y1,Y2,signs,c,first_time=False,method='TwoY',dim=1,set_c=[],fnr=0,hecke_ef=1,verbose=0):
    r"""
    Computes the functional we use as an indicator of an eigenvalue.

    INPUT:

    -''S'' -- space of Maass waveforms
    -''r'' -- real
    -''M'' -- integer
    -''Y1''-- real
    -''Y2''-- real
    -''signs'' -- dict
    -''c'' -- set which coefficients to use
    -''ST''-- normalization
    -''first_time'' --
    -''method'' -- string in ['TwoY','Innerprod']
    -''ndigs'' -- integer (number of digits wanted)

    OUTPUT:

    -  list of real values



    """
    diffsx=dict()
    h=0
    logger2.info("r,Y1,Y2={0}".format((r,Y1,Y2)))
    logger2.info("prec={0}".format(r.prec()))
    logger2.info("dim={0}".format(dim)+'\t method={0}'.format(method))
    logger2.info("set_c={0}".format(set_c))
    prec = r.prec()
    sym_type=S.sym_type()
    cusp_evs=S.atkin_lehner_eigenvalues() ## a dict
    Norm=S.set_norm(dim,set_c=set_c)
    #if S.group().is_congruence():
    #    raise ValueError,"Method for noncongruence groups"
    Q=M+20
    logger2.info("r,Y1,M,Q={0}".format((r,Y1,M,Q)))
    logger2.debug("Norm={0}".format(Norm))
    logger2.debug("c = {0}, signs={1}".format(c,signs))
    C1=get_coeff_fast_cplx_dp(S,r,Y1,M,Q,Norm)
    C1=C1[fnr]
    logger2.info("do_cplx: r,Y2,M,Q={0}".format(r,Y2,M,Q,Norm,cusp_evs))
    C2=get_coeff_fast_cplx_dp(S,RR(r),float(Y2),M,Q,Norm)
    C2 = C2[fnr]
    numt = len(c.keys())
    for j in c.keys():  #
        ## To avoid 'blowups' we sometimes Look at relative difference(?)
        den = max(1,abs(C1[0][c[j]])+abs(C2[0][c[j]]))
        if isinstance(C1[0][c[j]],ComplexNumber):
                diffsx[j]=(C1[0][c[j]]-C2[0][c[j]]).real()/den
                diffsx[numt+j]=(C1[0][c[j]]-C2[0][c[j]]).imag()/den
        elif isinstance(C1[0][c[j]],complex):
                diffsx[j]=(C1[0][c[j]]-C2[0][c[j]]).real/den
                diffsx[numt+j]=(C1[0][c[j]]-C2[0][c[j]]).imag/den
        else:
            raise TypeError("Got coefficients of unknown type {0}".format(type(C1[0][c[j]])))
        logger2.debug("C1[{0}]={1}".format(c[j],C1[0][c[j]]))
        logger2.debug("C2[{0}]={1}".format(c[j],C2[0][c[j]]))
        logger2.debug("abs(C1[{0}])+abs(C2[{0}])={1}".format(c[j],den))
        logger2.debug("diffs(real)={0}".format((C1[0][c[j]]-C2[0][c[j]]).real))
        logger2.debug( "diffs/den={0}".format(diffsx[j]))
    h=0.0
    for j in [1,2,3]:
        logger2.debug("C1:2[{0}]={1} : {2}".format(c[j],C1[0][c[j]],C2[0][c[j]]))
        logger2.debug("diffsx[1]={0}".format(diffsx[1]))
    logger2.debug("c={0}".format(c))
    logger2.debug("signs={0}".format(signs))
    ## Now we add up all contributions using the same alignment of signs as previous point..
    abh=0
    for j in diffsx.keys():
        if not first_time and list(signs.keys()).count(j)>0:
            h=h+signs[j]*diffsx[j]
        else:
            if not signs.keys():
                # We use the sign of the first test.
                for k in diffsx.keys():
                    if diffsx[1] != 0:
                        signs[k]=sign(diffsx[k])/sign(diffsx[1])
                    else:
                        signs[k]=sign(diffsx[k])
                    #print "sgns[{0}".format(k,"]={0}".format(signs[k]
            h=h+signs[j]*diffsx[j]
        abh = abh+abs(diffsx[j])
        logger2.debug("abs(h)={0}".format(abh))
    return [diffsx,h]

    
def is_derivative_large(R1,R2,d1,d2,method='TwoY',hmax=100,tol=1E-10,verbose=0):
    t=dict()
    if abs(R2-R1)<tol:
        return 0
    for j in d1.keys():
        if j<0:
            continue
        t[j]=(d1[j]-d2[j])/(R2-R1)
        #if verbose>0:
        # print "diff(",j,":",R1,")=",(d1[j]-d2[j]),"/",(R2-R1),"=",t[j]
    #print "t=",t
    #print "t.values()=",t.values()
    #print "abs(t)=",map(abs,t.values())
    #print "max(abs(t)=",max(map(abs,t.values()))
    mmax = max(map(abs,t.values()))
    if mmax>hmax:
        if verbose>0:
            flogger.debug("too large derivative! max={0}".format(mmax))
        return 1
    return 0


def is_minimum_at(h,diffs,tol=1E-12,verbose=0):
    a1=sum(map(abs,diffs[1].values()))
    if a1<tol:
        return 1,a1
    a2=sum(map(abs,diffs[1].values()))
    if a2<tol:
        return 2,a2
    a3=sum(map(abs,diffs[1].values()))
    if a3<tol:
        return 3,a3
    return 0,0


def is_zero_in(st,h,diffs,tol=1E-12,verbose=0):
    r"""
    Tells which interval contains changes of sign.
    and make sure that the individual tests are also having sign-changes in the same interval.

    INPUT:

    - ''h'' -- dictionary h[j]=[h1,h2,h3]
    - `diffs` -- dictionary with differences at points
    OUTPUT:

    - integer  (-1 0 1)

    EXAMPLES::


    """
    if verbose>2:
        lvl0 = logger2.getEffectiveLevel()
        logger2.setLevel(int(10))
    else:
        logger2.setLevel(int(40))
    #for x in diffs.keys():
    #    logger2.debug("diffs[{0}]={1}".format(x,diffs[x]))
    zi=dict(); i=0
    ok=1
    a1=0; a2=0; a3=0
    for j in diffs[1].keys():
        if j <=0:
            continue
        a1+=abs(diffs[1][j])
        a2+=abs(diffs[2][j])
        a3+=abs(diffs[3][j])
        
    logger2.debug("a1={0} \t a2={1}\t a3={2}".format(a1,a2,a3))

    ## Try to decide which interval: [R1,R2] or [R2,R3] has the best chance of containing a zero
    
    ## First look for simultaneous sign changes:
    ## First actual sign change and then approximate...
    logger2.debug("tol :{0}".format(tol))
    zero_count={}
    ##  If we don't use symmetry we are using tests for both real and imaginary parts of coefficients.
    if st == -1:
        l = 0
        for j in diffs[1].keys():
            if j>0:
                l+=1
        l = QQ(l)/QQ(2)
#    print "keys=",diffs[1].keys()
#    print "l=",l
    for toll in [0,tol]:
        zero_count[toll]=[{},{}]
        for j in diffs[1].keys():
            zero_count[toll][0][j]=0
            zero_count[toll][1][j]=0                       
        ok=1; i=0
        zero_in_1=1; zero_in_2=1;
        for j in diffs[1].keys():
            if j<=0: continue
            if diffs[1][j]*diffs[2][j]<toll:
                zero_count[toll][0][j]=1
        for j in diffs[2].keys():
            if j<=0: continue
            if diffs[2][j]*diffs[3][j]<toll:
                zero_count[toll][1][j]=1
        if st==-1:
            for j in diffs[1].keys():
                if j>l: continue
                if zero_count[toll][0][j]==0 and zero_count[toll][0][j+l]==0:
                    zero_in_1=0
                    break
            for j in diffs[1].keys():
                if j>l: continue
                if zero_count[toll][1][j]==0 and zero_count[toll][1][j+l]==0:
                    zero_in_2=0
                    break
        else:
            for j in diffs[1].keys():
                if j<=0: continue
                if zero_count[toll][0][j]==0:
                    zero_in_1=0
                if zero_count[toll][1][j]==0:
                    zero_in_2=0
        logger2.debug("zero_in_1,2= {0} - {1}".format(zero_in_1,zero_in_2))
        logger2.debug("zero_count={0}".format(zero_count))
#        if zero_in_1==1 and zero_in_2 == 1:
#            if zero_count[toll][1].count(1)>zero_count[toll][0].count(1):
#                zero_in_1 = 0
#       if verbose>0:
#            print "zero_in_1,2=",zero_in_1,zero_in_2
        if zero_in_1==1 and zero_in_2==1:
            # possible zero in both intervals
            ## See if there is a true zero in one of the intervals
            if list(zero_count[0][0].values()).count(1)>list(zero_count[0][1].values()).count(1):
                i=-1
            elif list(zero_count[0][0].values()).count(1)<list(zero_count[0][1].values()).count(1):
                i=1
            elif a1<a2:
                i=-1
            else:
                i=1
        elif zero_in_1==1:
            i=-1
        elif zero_in_2==1:
            i=1
        else:
            i=0
        if toll==0 and i != 0:
            break
        if toll>0 and i != 0:
            ## make extra check
            logger2.debug("Check extra:i={0}".format(i))
            if i==-1:
                a1 = max(map(abs,diffs[1].values()))
                a2 = max(map(abs,diffs[2].values()))
                logger2.debug("a1={0},a2={1}".format(a1,a2))

                if max(a1,a2)>sqrt(tol):
                    i=0
            else:
                a2 = max(list(map(abs,list(diffs[2].values()))))
                a3 = max(list(map(abs,list(diffs[3].values()))))
                logger2.debug("a2={0},a3={1}".format(a2,a3))
                if max(a2,a3)>sqrt(tol):
                    i=0
    if verbose>0:
        logger2.debug("zero in :".format(i))
    if i==0: ## We make a more loose test for zero before giving up
        a1=sum(map(abs,diffs[1].values()))
        a2=sum(map(abs,diffs[2].values()))
        a3=sum(map(abs,diffs[3].values()))
        if verbose>0:
            logger2.debug("a1={0}\t a2={1}\t a3={2}".format(a1,a2,a3))
    ## If the differences are really small at both endpoints of the interval
    ## we "probably" have a zero...
        aa1=max(a1,a2)
        aa2=max(a2,a3)
        if aa1<aa2 and (aa1)<tol:
            if verbose>0:
                logger2.debug("aa1<aa2 and aa1={0} < tol={1}".format(aa1,tol))
                logger2.setLevel(lvl0)
            return -1
        if aa2<aa1 and (aa2)<tol:
            if verbose>0:
                logger2.debug("aa2<aa1 and aa2={0} < tol={1}".format(aa2,tol))
                logger2.setLevel(lvl0)

            return 1
    #if h[1]*h[2] < 0:
    #    zi[-1]=1; i=-1
    #if h[3]*h[2] < 0:
    #    zi[1]=1; i=1
    #if zi.values().count(1) >1: # need to split
    #    return -2
    #print "i1=",i
    #s="Need to split! Not implemented!"
    #raise ValueError,s
    if verbose>2:
        logger2.setLevel(lvl0)

    return i
