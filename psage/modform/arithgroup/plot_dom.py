
""" Excerpt from my MySubgroup class. Contains routines to draw fundamental domains.

r"""
from __future__ import print_function
from __future__ import division
from builtins import str
from builtins import range
from builtins import object
import matplotlib
import matplotlib.patches as patches
import matplotlib.path as path
from matplotlib.backends.backend_agg import FigureCanvasAgg, FigureCanvas
from sage.all import deepcopy,hyperbolic_arc
from sage.all import I,Gamma0,Gamma1,Gamma,SL2Z,ZZ,RR,ceil,sqrt,CC,line,text,latex,exp,pi,infinity
from sage.all import deepcopy
from sage.plot.all import Graphics
from sage.plot.plot import (Graphics,line)
from sage.functions.trig import (cos,sin)
from sage.plot.plot import line
from sage.functions.trig import arcsin


from matplotlib import pyplot

from sage.plot.all import Graphics,arc
from sage.plot.arc import Arc
from sage.plot.line import Line
from sage.plot.point import Point
from sage.plot.primitive import GraphicPrimitive
from sage.plot.colors import to_mpl_color
from sage.plot.plot import options
from sage.misc.decorators import rename_keyword
from sage.plot.all import text
from sage.plot.hyperbolic_polygon import HyperbolicPolygon, hyperbolic_triangle
from sage.plot.hyperbolic_arc import hyperbolic_arc
from sage.plot.bezier_path import BezierPath

def draw_fundamental_domain(N,group='Gamma0',model="H",axes=None,filename=None,**kwds):
    r""" Draw fundamental domain
    INPUT:
        - ''model'' -- (default ''H'')
        - ''H'' -- Upper halfplane
        - ''D'' -- Disk model
        - ''filename''-- filename to print to
        - ''**kwds''-- additional arguments to matplotlib 
        - ''axes''  -- set geometry of output
        =[x0,x1,y0,y1] -- restrict figure to [x0,x1]x[y0,y1]

    EXAMPLES::

        sage: G=MySubgroup(Gamma0(3))
        sage: G.draw_fundamental_domain()

    """
    G=eval(group+'('+str(N)+')')
    #print G
    name ="$"+latex(G)+"$"
    ## need a "nice" set of coset representatives to draw a connected fundamental domain. Only implemented for Gamma_0(N)
    coset_reps = nice_coset_reps(G)
    if(model=="D"):
        g=draw_funddom_d(coset_reps,format,I)
    else:
        g=draw_funddom(coset_reps,format)
    if axes!=None:
        [x0,x1,y0,y1]=axes
    elif model=="D":
        x0=-1 ; x1=1 ; y0=-1.1 ; y1=1 
    else:
        # find the width of the fundamental domain
        w=0  #self.cusp_width(Cusp(Infinity))
        wmin=0 ; wmax=1
        max_x = RR(0.55)
        rho = CC( exp(2*pi*I/3))
        for V in coset_reps:
            ## we also compare the real parts of where rho and infinity are mapped
            r1 = (V.acton(rho)).real()
            if V[1,0] != 0:
                inf1 = RR(V[0,0] / V[1,0])
            else:
                inf1 = 0
            if(V[1 ,0 ]==0  and V[0 ,0 ]==1 ):
                if(V[0 ,1 ]>wmax):
                    wmax=V[0 ,1 ]
                if(V[0 ,1 ]<wmin):
                    wmin=V[0 ,1 ]
            if( max(r1,inf1) > max_x):
                max_x = max(r1,inf1)
            #print "wmin,wmax=",wmin,wmax
            #x0=-1; x1=1; y0=-0.2; y1=1.5
        x0=RR(-max_x) ; x1=RR(max_x) ; y0=RR(-0.15) ; y1=RR(1.5) 
    ## Draw the axes manually (since  can't figure out how to do it automatically)
    ax = line([[x0,0.0],[x1,0.0]],color='black')
        #ax = ax + line([[0.0,y0],[0.0,y1]],color='black')
        ## ticks
    ax = ax + line([[-0.5,-0.01],[-0.5,0.01]],color='black')
    ax = ax + line([[0.5,-0.01],[0.5,0.01]],color='black')
    g = g + ax
    if model=="H":
        t = text(name, (0, -0.1), fontsize=16, color='black')
    else:
        t = text(name, (0, -1.1), fontsize=16, color='black')
        g = g + t
        g.set_aspect_ratio(1)
        g.set_axes_range(x0,x1,y0,y1)
        g.axes(False)
        if(filename!=None):
            fig = g.matplotlib()
            fig.set_canvas(FigureCanvasAgg(fig))
            axes = fig.get_axes()[0]
            axes.minorticks_off()
            axes.set_yticks([])
            fig.savefig(filename,**kwds)
        else:
            return g
#        g.show(figsize=[5,5])

def draw_funddom(coset_reps,format="S"):
    r""" Draw a fundamental domain for G.
    
    INPUT:
    
    - ``format``  -- (default 'Disp') How to present the f.d.
    -   ``S`` -- Display directly on the screen
    
    EXAMPLES::        


    sage: G=MySubgroup(Gamma0(3))
    sage: G._draw_funddom()
        
    """
    pi=RR.pi()
    pi_3 = pi / RR(3.0)
    g=Graphics()
    x1=RR(-0.5) ; y1=RR(sqrt(3 )/2 )
    x2=RR(0.5) ; y2=RR(sqrt(3 )/2 )
    xmax=RR(20.0) 
    l1 = line([[x1,y1],[x1,xmax]])
    l2 = line([[x2,y2],[x2,xmax]])
    l3 = line([[x2,xmax],[x1,xmax]]) # This is added to make a closed contour
    c0=_circ_arc(RR(pi/3.0) ,RR(2.0*pi)/RR(3.0) ,0 ,1 ,100 )
    tri=c0+l1+l3+l2
    g=g+tri
    for A in coset_reps:
        if list(A)==[1,0,0,1]:
            continue

        tri=draw_transformed_triangle_H(A,xmax=xmax)
        g=g+tri
    return g

def draw_transformed_triangle_H(A,xmax=20):
    r"""
    Draw the modular triangle translated by A=[a,b,c,d]
    """
    #print "A=",A,type(A)
    pi=RR.pi()
    pi_3 = pi / RR(3.0)
    x1=RR(-0.5) ; y1=RR(sqrt(3 )/2 )
    x2=RR(0.5) ; y2=RR(sqrt(3 )/2 )
    a,b,c,d = A #[0,0]; b=A[0,1]; c=A[1,0]; d=A[1,1]
    if a<0:
        a=RR(-a); b=RR(-b); c=RR(-c); d=RR(-d) 
    else:
        a=RR(a); b=RR(b); c=RR(c); d=RR(d) 
    if c==0: # then this is easier
        if a*d != 0:
            a=a/d; b=b/d; 
        L0 = [[a*cos(pi_3*RR(i/100.0))+b,a*sin(pi_3*RR(i/100.0))] for i in range(100 ,201 )]
        L1 = [[a*x1+b,a*y1],[a*x1+b,xmax]]
        L2 = [[a*x2+b,a*y2],[a*x2+b,xmax]]
        L3 = [[a*x2+b,xmax],[a*x1+b,xmax]]
        c0=line(L0); l1=line(L1); l2=line(L2); l3=line(L3)
        tri=c0+l1+l3+l2
    else:
        den=(c*x1+d)**2 +c**2 *y1**2 
        x1_t=(a*c*(x1**2 +y1**2 )+(a*d+b*c)*x1+b*d)/den
        y1_t=y1/den
        den=(c*x2+d)**2 +c**2 *y2**2 
        x2_t=(a*c*(x2**2 +y2**2 )+(a*d+b*c)*x2+b*d)/den
        y2_t=y2/den
        inf_t=a/c
        #print "A=",A
        #print "arg1=",x1_t,y1_t,x2_t,y2_t
        c0=_geodesic_between_two_points(x1_t,y1_t,x2_t,y2_t)
        #print "arg1=",x1_t,y1_t,inf_t
        c1=_geodesic_between_two_points(x1_t,y1_t,inf_t,0. )
        #print "arg1=",x2_t,y2_t,inf_t
        c2=_geodesic_between_two_points(x2_t,y2_t,inf_t,0.0)
        tri=c0+c1+c2
    return tri
    

class HyperbolicTriangle(HyperbolicPolygon):
    """
    Primitive class for hyberbolic triangle type. See ``hyperbolic_triangle?``
    for information about plotting a hyperbolic triangle in the complex plane.

      INPUT:

    - ``a,b,c`` - coordinates of the hyperbolic triangle in the upper
      complex plane

    - ``options`` - dict of valid plot options to pass to constructor

    EXAMPLES:

    Note that constructions should use ``hyperbolic_``::

         sage: from sage.plot.hyperbolic_triangle import HyperbolicTriangle
         sage: print(HyperbolicTriangle(0, 1/2, I, {}))
         Hyperbolic triangle (0.000000000000000, 0.500000000000000, 1.00000000000000*I)

    Note: The main reason for this extension is that we sometimes want to draw only certain sides of the triangle.

    """
    def __init__(self, A, B, C, options):
        """
        Initialize HyperbolicTriangle:

        Examples::

            sage: from sage.plot.hyperbolic_triangle import HyperbolicTriangle
            sage: print(HyperbolicTriangle(0, 1/2, I, {}))
            Hyperbolic triangle (0.000000000000000, 0.500000000000000, 1.00000000000000*I)

        """
        A, B, C = (CC(A), CC(B), CC(C))
        self.path = []
        sides = options.pop('sides',[1,2,3])
        verbose=options.pop('verbose',0)
        sides.sort()
        if sides == [1]:
            if verbose>0:
                print("Drawing A - B!")
            self._hyperbolic_arc(A, B, True)
        elif sides == [2]:
            if verbose>0:
                print("Drawing B - C!")                
            self._hyperbolic_arc(B, C, True)
        elif sides == [3]:
            if verbose>0:
                print("Drawing C - A!")                
            self._hyperbolic_arc(C, A, True)
        elif sides == [1,2]:
            if verbose>0:
                print("Drawing A - B! & B - C!")
            self._hyperbolic_arc(A, B, True);
            self._hyperbolic_arc(B, C, False)
        elif sides == [1,3]:
            if verbose>0:
                print("Drawing C - A! & A - B")
            self._hyperbolic_arc(C, A,True)
            self._hyperbolic_arc(A, B, False)
        elif sides == [2,3]:
            if verbose>0:
                print("Drawing B - C! & C - A")
            self._hyperbolic_arc(B, C,True)
            self._hyperbolic_arc(C, A, False)            
        else:
            self._hyperbolic_arc(A, B,True)            
            self._hyperbolic_arc(B, C,False)
            self._hyperbolic_arc(C, A, False)            
            
        BezierPath.__init__(self, self.path, options)
        self.A, self.B, self.C = (A, B, C)
        self._pts = [A,B,C]

class HyperbolicTriangleDisc(object): #]GraphicPrimitive):
    r"""
    Hyperbolic triangles in the disc model of the hyperbolic plane.
    Note that we are given coordinates in the upper half-plane and map them to the disc.
    
    """
    @options(center=CC(0,1))
    def __init__(self, A, B, C,**options):
        """
        Initialize HyperbolicTriangle under the map (z-z0)/(z-\bar(z0)):
        
        Examples::
        
            sage: from sage.plot.hyperbolic_triangle import HyperbolicTriangle
            sage: print(HyperbolicTriangle(0, 1/2, I, {}))
            Hyperbolic triangle (0.000000000000000, 0.500000000000000, 1.00000000000000*I)
        """
        A, B, C = (CC(A), CC(B), CC(C))
        self.path = []
        self._options = {}
        self._graphics = Graphics()
        self._verbose = options.pop('verbose',None)
        Z0 = options['center'] #.get('center',C(0,1))
        options.pop('center',None)
        self._npts = options.pop('npts',10)
        sides = options.pop('sides',[1,2,3])
        #options.pop('fill',None)
        self._options.update(options)
        self._z0=CC(Z0); self._z0bar=CC(Z0).conjugate()        
        verbose=self._verbose
        sides.sort()
        if sides == [1]:
            if verbose>0:
                print("Drawing A - B!")
            self._hyperbolic_arc_d(A, B, True)
        elif sides == [2]:
            if verbose>0:
                print("Drawing B - C!")                
            self._hyperbolic_arc_d(B, C, True)
        elif sides == [3]:
            if verbose>0:
                print("Drawing C - A!")                
            self._hyperbolic_arc_d(C, A, True)
        elif sides == [1,2]:
            if verbose>0:
                print("Drawing A - B! & B - C!")
            self._hyperbolic_arc_d(A, B, True)
            self._hyperbolic_arc_d(B, C,False)
        elif sides == [1,3]:
            if verbose>0:
                print("Drawing C - A! & A - B")
            self._hyperbolic_arc_d(C, A,True)
            self._hyperbolic_arc_d(A, B, False)
        elif sides == [2,3]:
            if verbose>0:
                print("Drawing B - C! & C - A")
            self._hyperbolic_arc_d(B, C,True)
            self._hyperbolic_arc_d(C, A, False)            
        else:
            self._hyperbolic_arc_d(A, B,True)            
            self._hyperbolic_arc_d(B, C,False)
            self._hyperbolic_arc_d(C, A, False)
        #self._hyperbolic_arc_d(A, B, True);
        #self._hyperbolic_arc_d(B, C);
        #self._hyperbolic_arc_d(C, A);
        #BezierPath.__init__(self, self.path, options)


        #super(HyperbolicTriangleDisc,self).__init__(options)
        self.A, self.B, self.C = (A, B, C)

    def _cayley_transform(self,z):
        #print "z=",z,z==infyinity,type(z),type(infinity)
        if z==infinity or z==CC(infinity):
            return CC(1,0)
        return (CC(z)-self._z0)/(CC(z)-self._z0bar)

    def _hyperbolic_arc_d(self, z0, z3, first=False):
        """
        Function to construct Bezier path as an approximation to
        the hyperbolic arc between the complex numbers z0 and z3 in the
        hyperbolic plane.
        """

        w0 = self._cayley_transform(z0)
        w3 = self._cayley_transform(z3)
        if self._verbose>0:
            print("in plane z0,z3={0},{1}".format(z0,z3))
            print("in disc: {0},{1}".format(w0,w3))
        npts = self._npts
        if z0 == infinity or z0==CC(infinity):
            zm = [z3 + CC(0,j+0.5) for j in range(npts-2)]
            wm = [self._cayley_transform(x) for x in zm]
            pts = [w3]
            pts.extend(wm)
            pts.append(w0)
            opt = self._options
            opt['fill']=False
            self._graphics.add_primitive(BezierPath([[(x.real(),x.imag()) for x in pts ]],opt))
            return 
        if z3 == infinity  or z3==CC(infinity):
            zm = [z0 + CC(0,j+0.5) for j in range(npts-2)]
            wm = [self._cayley_transform(x) for x in zm]
            pts = [w0]
            pts.extend(wm)
            pts.append(w3)
            opt = self._options
            opt['fill']=False
            self._graphics.add_primitive(BezierPath([[(x.real(),x.imag()) for x in pts ]],opt))
            #self._graphics.add_primitive(Line([w1.real(),w2.real()],[w1.imag(),w2.imag()],self._options))
            #self.path.append([(0,w0.imag()),CC(0,y), (0,w3.imag())])
            return
        x0=z0.real(); y0 = z0.imag()
        x3=z3.real(); y3 = z3.imag()
        if y0 == 0 and y3 == 0:
            p = (z0.real()+z3.real())/2
            r = abs(z0-p)
            zm = CC(p, r)
            self._hyperbolic_arc_d(z0, zm, first)
            self._hyperbolic_arc_d(zm, z3)
            return                   
        else:
            if abs(x0-x3)<1e-10:  ## on the same vertical line
                xmid = (x0+x3)*0.5; h = y3-y0
                zm = [ CC(xmid,y0+t*h/(npts-1)) for t in range(npts) ]
            else:
                p = RR((x0+x3)*(x3-x0)+(y0+y3)*(y3-y0))/(2*(x3-x0)) 
                r = RR((p - x0)**2 + y0**2).sqrt()  # radius of the circle in H
                zm = ((z0+z3)/2-p)/abs((z0+z3)/2-p)*r+p  # midpoint (at least approximately) of geodesic between z0 and z3
                t0 = CC(z0 - p).argument()
                t3 = CC(z3 - p).argument()
                if self._verbose>1:
                    print("x0,x3={0},{1}".format(x0,x3))
                    print("t0,t3={0},{1}".format(t0,t3))
                    print("r={0}".format(r))
                    print("opt={0}".format(self._options))
                if x0 <= x3:
                    zm = [p + r*CC(0,(t0+t*(t3-t0)/(npts-1))).exp() for t in range(npts)]
                else:
                    zm = [p + r*CC(0,(t0+t*(t3-t0)/(npts-1))).exp() for t in range(npts)]            
                #print "zm=",zm
                #zm.insert(0,z0)
                #zm.append(z3)
            pts = [self._cayley_transform(x) for x in zm]
            opt = self._options
            opt['fill']=False
            #w0 = self._cayley_transform(z0)
            #w3 = self._cayley_transform(z3)

            if self._verbose>2:
                print("pts={0}".format(pts))
            self._graphics.add_primitive(BezierPath([[(x.real(),x.imag()) for x in pts ]],opt))
            return 
            #print "z0_test=",(p+r*exp(t0*I))
                #print "z3_test=",(p+r*exp(t3*I))
#t = (8*zm-4*(z0+z3)).imag()/3/(z3-z0).real()
            # I have no idea what these points should represent....
            #z1 = z0 + t*CC(z0.imag(), (p-z0.real()))
            #z2 = z3 - t*CC(z3.imag(), (p-z3.real()))
            wm = [self._cayley_transform(x) for x in zm]
            pp = self._cayley_transform(CC(p,0))
            w1 = self._cayley_transform(z0)
            w2 = self._cayley_transform(z3)
            c = self._cayley_transform(CC(p,0)) # center of circle on the unit disk.
            if self._verbose>2:
                print("p,r={0},{1}".format(p,r))
                print("zm={0}".format(zm))
                #print "t=",t
                #print "tt=",(8*zm-4*(z0+z3)).imag()/3/(z3-z0).real()

                print("C(c)={0}".format(pp))
                print("C(zm)={0}".format(wm))
                print("C(z0)={0}".format(w1))
                print("C(z3)={0}".format(w2))
                #print "z2=",z2
            
            r = abs(w1-c) # radius
            t0 = CC(w0 - pp).argument()
            t3 = CC(w3 - pp).argument()
            t = abs(t0-t3)
        if self._verbose>0:
            print("adding a:rc {0}".format((zm.real(),zm.imag(),r,r,t,t0,t3)))
        self._graphics.add_primitive(Line([w1.real(),w2.real(),wm.real()],[w1.imag(),w2.imag(),wm.imag()],{'thickness':2,'alpha':1,
                                                                                       'rgbcolor':'blue',
                                                                                                         'legend_label':""}))
        self._graphics.add_primitive(Point([w1.real(),w2.real(),wm.real()],[w1.imag(),w2.imag(),wm.imag()],{'size':10,'alpha':1,
                                                                                                            'faceted':True,
                                                                                                            'rgbcolor':'red',
                                                                                                            'legend_label':""}))
        
        #self._graphics.add_primitive(Arc(pp.real(),pp.imag(),r,r,t,t0,t3,self._options))
        #self._graphics. Arc(zm.real(),zm.imag(),r,r,abs(t),-abs(t),abs(t),self._options)

    def __call__(self):
        return self._graphics
        # if first:
        #     self.path.append([(w0.real(), w0.imag()),
        #                       (w1.real(), w1.imag()),
        #                       (w2.real(), w2.imag()),
        #                       (w3.real(), w3.imag())]);
        #     first = False
        # else:
        #     self.path.append([(w1.real(), w1.imag()),
        #                       (w2.real(), w2.imag()),
        #                       (w3.real(), w3.imag())]);


@rename_keyword(color='rgbcolor')
#@options(alpha=1,fill=False,  thickness=1, rgbcolor="blue", zorder=2, linestyle='solid', model='H')
@options(alpha=1,  thickness=1, rgbcolor="black", zorder=2, linestyle='solid', model='H',fill=False)
def my_hyperbolic_triangle(a, b, c, **options):
    r"""
    Return a hyperbolic triangle in the complex hyperbolic plane with points
    (a, b, c). Type ``?hyperbolic_triangle`` to see all options.

    INPUT:

    - ``a, b, c`` - complex numbers in the upper half complex plane

    OPTIONS:
    
    - ``alpha`` - default: 1
       
    - ``fill`` - default: False

    - ``thickness`` - default: 1

    - ``rgbcolor`` - default: 'blue'

    - ``linestyle`` - default: 'solid'

    EXAMPLES:

    Show a hyperbolic triangle with coordinates 0, `1/2+i\sqrt{3}/2` and
    `-1/2+i\sqrt{3}/2`::
    
         sage: hyperbolic_triangle(0, -1/2+I*sqrt(3)/2, 1/2+I*sqrt(3)/2)

    A hyperbolic triangle with coordinates 0, 1 and 2+i::
    
         sage: hyperbolic_triangle(0, 1, 2+i, fill=true, rgbcolor='red')
    """
    g = Graphics()
    g._set_extra_kwds(g._extract_kwds_for_show(options))
    model = options['model']
    npts = options.get('npts',10)
    sides = options.pop('sides',[1,2,3]) ## Which of the sides we will draw.
    verbose = options.get('verbose',0)
    options.pop('model',None); options.pop('method',None)
    if model=="D":
        #g.add_primitive(HyperbolicTriangleDisc(a, b, c, **options))
        #print "options=",options
        options['sides']=sides
        H = HyperbolicTriangleDisc(a, b, c, **options)
        g += H()
    else:
        options.pop('npts',None)
        if sides == [1,2,3]:
            if verbose>0:
                print("adding HyperbolicTriangle({0}, {1}, {2},options={3})".format(a,b,c,options))
            options.pop('verbose',0)
            ## check if We need my class or the original class
            g.add_primitive(HyperbolicTriangle(a, b, c, options))
        else:
            options['sides']=sides
            if verbose>0:
                print("adding HyperbolicTriangle({0}, {1}, {2},options={3})".format(a,b,c,options))
            g.add_primitive(HyperbolicTriangle(a, b, c, options))                           
    g.set_aspect_ratio(1)
    return g


def draw_funddom_d(coset_reps,format="MP",z0=I,verbose=0):
    r""" Draw a fundamental domain for self in the circle model
    INPUT:
    - ''format''  -- (default 'Disp') How to present the f.d.
    =  'S'  -- Display directly on the screen
    - z0          -- (default I) the upper-half plane is mapped to the disk by z-->(z-z0)/(z-z0.conjugate())
    EXAMPLES::
        

    sage: G=MySubgroup(Gamma0(3))
    sage: G._draw_funddom_d()
        
    """
    # The fundamental domain consists of copies of the standard fundamental domain
    pi=RR.pi()
    from sage.plot.plot import (Graphics,line)
    g=Graphics()
    bdcirc=_circ_arc(0 ,2 *pi,0 ,1 ,1000 )
    g=g+bdcirc
    # Corners
    x1=-RR(0.5); y1=RR(sqrt(3 )/2)
    x2=RR(0.5); y2=RR(sqrt(3 )/2)
    z_inf=1 
    l1 = _geodesic_between_two_points_d(x1,y1,x1,infinity)
    l2 = _geodesic_between_two_points_d(x2,y2,x2,infinity)
    c0 = _geodesic_between_two_points_d(x1,y1,x2,y2)
    tri=c0+l1+l2
    g=g+tri
    for A in coset_reps:
        a,b,c,d=A
        if a==1  and b==0  and c==0  and d==1:
            continue
        if a<0:
            a=-a; b=-b; c=-c; d=-d 
        if verbose>0:
                print("a,b,c,d={0},{1},{2},{3}".format(a,b,c,d))
        if c==0: # then this is easier
            l1 = _geodesic_between_two_points_d(x1+b,y1,x1+b,infinity)
            l2 = _geodesic_between_two_points_d(x2+b,y2,x2+b,infinity)
            c0 = _geodesic_between_two_points_d(x1+b,y1,x2+b,y2)
            # c0=line(L0); l1=line(L1); l2=line(L2); l3=line(L3)
            tri=c0+l1+l2
            g=g+tri
        else:
            den=(c*x1+d)**2 +c**2 *y1**2 
            x1_t=(a*c*(x1**2 +y1**2 )+(a*d+b*c)*x1+b*d)/den
            y1_t=y1/den
            den=(c*x2+d)**2 +c**2 *y2**2 
            x2_t=(a*c*(x2**2 +y2**2 )+(a*d+b*c)*x2+b*d)/den
            y2_t=y2/den            
            inf_t=a/c
            if verbose>0:
                    print("x1_t=",x1_t)
                    print("y1_t=",y1_t)
                    print("x2_t=",x2_t)
                    print("y2_t=",y2_t)
                    print("inf_t=",inf_t)
            c0=_geodesic_between_two_points_d(x1_t,y1_t,x2_t,y2_t)
            c1=_geodesic_between_two_points_d(x1_t,y1_t,inf_t,1.0 )
            c2=_geodesic_between_two_points_d(x2_t,y2_t,inf_t,1.0 )
            tri=c0+c1+c2
            g=g+tri
    g.xmax(1 )
    g.ymax(1 )
    g.xmin(-1 )
    g.ymin(-1 )
    g.set_aspect_ratio(1 )
    return g


#### Methods not dependent explicitly on the group 
def _geodesic_between_two_points(x1,y1,x2,y2):
    r""" Geodesic path between two points hyperbolic upper half-plane

    INPUTS:
    
    - ''(x1,y1)'' -- starting point (0<y1<=infinity)
    - ''(x2,y2)'' -- ending point   (0<y2<=infinity)
    - ''z0''  -- (default I) the point in the upper corresponding
                 to the point 0 in the disc. I.e. the transform is
                 w -> (z-I)/(z+I)
    OUTPUT:

    - ''ca'' -- a polygonal approximation of a circular arc centered
    at c and radius r, starting at t0 and ending at t1

    
    EXAMPLES::


        sage: l=_geodesic_between_two_points(0.1,0.2,0.0,0.5)
    
    """
    pi=RR.pi()
    #print "z1=",x1,y1
    #print "z2=",x2,y2
    if( abs(x1-x2)< 1E-10):
        # The line segment [x=x1, y0<= y <= y1]
        return line([[x1,y1],[x2,y2]])  #[0,0,x0,infinity]
    c=RR(y1**2 -y2**2 +x1**2 -x2**2 )/RR(2 *(x1-x2))
    r=RR(sqrt(y1**2 +(x1-c)**2 ))
    r1=RR(y1/r); r2=RR(y2/r)
    if(abs(r1-1 )< 1E-12 ):
        r1=RR(1.0)
    elif(abs(r2+1 )< 1E-12 ):
        r2=-RR(1.0)
    if(abs(r2-1 )< 1E-12 ):
        r2=RR(1.0)
    elif(abs(r2+1 )<1E-12 ):
        r2=-RR(1.0)
    if(x1>=c):
        t1 = RR(arcsin(r1))
    else:
        t1 = RR(pi)-RR(arcsin(r1))
    if(x2>=c):
        t2 = RR(arcsin(r2))
    else:
        t2 = RR(pi)-arcsin(r2)
    tmid = (t1+t2)*RR(0.5)
    a0=min(t1,t2)
    a1=max(t1,t2)
    #print "c,r=",c,r
    #print "t1,t2=",t1,t2
    return _circ_arc(t1,t2,c,r)

def _geodesic_between_two_points_d(x1,y1,x2,y2,z0=I):
    r""" Geodesic path between two points represented in the unit disc
         by the map w = (z-I)/(z+I)
    INPUTS:
    - ''(x1,y1)'' -- starting point (0<y1<=infinity)
    - ''(x2,y2)'' -- ending point   (0<y2<=infinity)
    - ''z0''  -- (default I) the point in the upper corresponding
                 to the point 0 in the disc. I.e. the transform is
                 w -> (z-I)/(z+I)
    OUTPUT:
    - ''ca'' -- a polygonal approximation of a circular arc centered
    at c and radius r, starting at t0 and ending at t1

    
    EXAMPLES::

        sage: l=_geodesic_between_two_points_d(0.1,0.2,0.0,0.5)
    
    """
    pi=RR.pi()
    from sage.plot.plot import line
    from sage.functions.trig import (cos,sin)    
    # First compute the points
    if(y1<0  or y2<0 ):
        raise ValueError("Need points in the upper half-plane! Got y1=%s, y2=%s" %(y1,y2))
    if(y1==infinity):
        P1=CC(1 )
    else:
        P1=CC((x1+I*y1-z0)/(x1+I*y1-z0.conjugate()))
    if(y2==infinity):
        P2=CC(1 )
    else:
        P2=CC((x2+I*y2-z0)/(x2+I*y2-z0.conjugate()))
        # First find the endpoints of the completed geodesic in D
    if(x1==x2):
        a=CC((x1-z0)/(x1-z0.conjugate()))
        b=CC(1 )
    else:
        c=RR(y1**2 -y2**2 +x1**2 -x2**2 )/RR(2 *(x1-x2))
        r=RR(sqrt(y1**2 +(x1-c)**2 ))
        a=c-r
        b=c+r
        a=CC((a-z0)/(a-z0.conjugate()))
        b=CC((b-z0)/(b-z0.conjugate()))
    if( abs(a+b) < 1E-10 ): # On a diagonal
        return line([[P1.real(),P1.imag()],[P2.real(),P2.imag()]])
    th_a=a.argument()
    th_b=b.argument()
    # Compute the center of the circle in the disc model
    if( min(abs(b-1 ),abs(b+1 ))< 1E-10  and  min(abs(a-1 ),abs(a+1 ))>1E-10 ):
        c=b+I*(1 -b*cos(th_a))/sin(th_a)
    elif( min(abs(b-1 ),abs(b+1 ))> 1E-10  and  min(abs(a-1 ),abs(a+1 ))<1E-10 ):
        c=a+I*(1 -a*cos(th_b))/RR(sin(th_b))
    else:
        cx=(sin(th_b)-sin(th_a))/sin(th_b-th_a)
        c=cx+I*(1 -cx*cos(th_b))/RR(sin(th_b))
    # First find the endpoints of the completed geodesic
    r=abs(c-a)
    t1=CC(P1-c).argument()
    t2=CC(P2-c).argument()
    #print "t1,t2=",t1,t2
    return _circ_arc(t1,t2,c,r)


def _circ_arc(t0,t1,c,r,num_pts=5000 ):
    r""" Circular arc
    INPUTS:
    - ''t0'' -- starting parameter
    - ''t1'' -- ending parameter
    - ''c''  -- center point of the circle
    - ''r''  -- radius of circle
    - ''num_pts''  -- (default 100) number of points on polygon
    OUTPUT:
    - ''ca'' -- a polygonal approximation of a circular arc centered
    at c and radius r, starting at t0 and ending at t1

    
    EXAMPLES::

        sage: ca=_circ_arc(0.1,0.2,0.0,1.0,100)
    
    """
    from sage.plot.plot import line,parametric_plot
    from sage.functions.trig import (cos,sin)
    from sage.all import var
    t00=t0; t11=t1
    ## To make sure the line is correct we reduce all arguments to the same branch,
    ## e.g. [0,2pi]
    pi=RR.pi()
    while(t00<0.0):
        t00=t00+RR(2.0*pi)
    while(t11<0):
        t11=t11+RR(2.0*pi)
    while(t00>2*pi):
        t00=t00-RR(2.0*pi)
    while(t11>2*pi):
        t11=t11-RR(2.0*pi)

    xc=CC(c).real()
    yc=CC(c).imag()
    num_pts=3
    t = var('t')
    if t11>t00:
        ca = parametric_plot((r*cos(t)+xc, r*sin(t)+yc), (t, t00,t11))
    else:
        ca = parametric_plot((r*cos(t)+xc, r*sin(t)+yc), (t, t11,t00))
    #L0 = [[RR(r*cos(t00+i*(t11-t00)/num_pts))+xc,RR(r*sin(t00+i*(t11-t00)/num_pts))+yc] for i in range(0 ,num_pts)]
    #ca=line(L0)
    return ca


def nice_coset_reps(G):
        r"""
        Compute a better/nicer list of right coset representatives [V_j]
        i.e. SL2Z = \cup G V_j
        Use this routine for known congruence subgroups.

        EXAMPLES::


            sage: G=MySubgroup(Gamma0(5))
            sage: G._get_coset_reps_from_G(Gamma0(5))
            [[1 0]
            [0 1], [ 0 -1]
            [ 1  0], [ 0 -1]
            [ 1  1], [ 0 -1]
            [ 1 -1], [ 0 -1]
            [ 1  2], [ 0 -1]
            [ 1 -2]]
    
        """
        cl=list()
        S,T=SL2Z.gens()
        lvl=G.generalised_level()
        # Start with identity rep.
        cl.append(SL2Z([1 ,0 ,0 ,1 ]))
        if(not S in G):
            cl.append(S)
        # If the original group is given as a Gamma0 then
        # the reps are not the one we want
        # I.e. we like to have a fundamental domain in
        # -1/2 <=x <= 1/2 for Gamma0, Gamma1, Gamma
        for j in range(1 , ZZ( ceil(RR(lvl/2.0))+2)):
            for ep in [1 ,-1 ]:
                if(len(cl)>=G.index()):
                    break
                # The ones about 0 are all of this form
                A=SL2Z([0 ,-1 ,1 ,ep*j])
                # just make sure they are inequivalent
                try:
                    for V in cl:
                        if((A!=V and A*V**-1  in G) or cl.count(A)>0 ):
                            raise StopIteration()
                    cl.append(A)
                except StopIteration:
                    pass
        # We now add the rest of the "flips" of these reps.
        # So that we end up with a connected domain
        i=1 
        while(True):
            lold=len(cl)
            for V in cl:
                for A in [S,T,T**-1 ]:
                    B=V*A
                    try:
                        for W in cl:
                            if( (B*W**-1  in G) or cl.count(B)>0 ):
                                raise StopIteration()
                        cl.append(B)
                    except StopIteration:
                        pass
            if(len(cl)>=G.index() or lold>=len(cl)):
                # If we either did not add anything or if we added enough
                # we exit
                break
        # If we missed something (which is unlikely)        
        if len(cl) != G.index():
            print("cl=",cl)
            raise ValueError("Problem getting coset reps! Need %s and got %s" %(G.index(),len(cl)))
        return cl

def get_contour(G,version=1,model='D',standalone=False,as_patch=True,translates=1,**kwds):
    r"""
    Compute the contour of a fundamental domain of G.
    :param G:
    :param version:
    :param model:
    :param standalone:
    :param as_patch:
    :param as_translated: if positive gives a number of translates of the given fundamental domain to include
    :param kwds:
    :return:
    """
    from .mysubgroups_alg import SL2Z_elt
    ymax = kwds.get('ymax',None)
    if translates > 1:
        original_reps = G.coset_reps()
        T = SL2Z_elt(1, 1, 0, 1)
        G._coset_reps_v0 = sum(
            [[(T ** j) * x for x in original_reps] for j in range(int(-translates / 2.), int(translates / 2. + 1))], [])
    if ymax is None:
        P=G.draw_fundamental_domain(version=version,method='a',model=model,fill=False,show_tesselation=False,contour=True,draw_circle=False,rgbcolor=kwds.get('color','red'),as_arcs=True)
    else:
        P=G.draw_fundamental_domain(version=version,method='a',model=model,fill=False,show_tesselation=False,contour=True,draw_circle=False,rgbcolor=kwds.get('color','red'),as_arcs=True,ymax=ymax)
    if translates > 1:
        # reset the representatives.
        G._coset_reps_v0 = original_reps
        #    #return P
    l=build_connected_path(P,**kwds) #,model=model)
    if standalone:
        fig=pyplot.Figure()
        canvas = FigureCanvas(fig)
        ax=fig.add_subplot(111)
        ax.add_patch(patches.PathPatch(l,facecolor='none'))
        #ax.plot(*zip(*l.vertices))
        return fig
    if as_patch:
        return patches.PathPatch(l,facecolor='none')
    return l


        
def build_connected_path(P,**kwds):
    paths = []
    ymax = P._axes_range.get('ymax',kwds.get('ymax',10000))
    xmax = P._axes_range.get('xmax',kwds.get('xmax',0))
    xmin = P._axes_range.get('xmin',kwds.get('xmin',0))
    new_paths = []
    if len(P)==1: ### SL2Z
        A = P[0].A; B=P[0].B; C=P[0].C
        new_paths =  [hyperbolic_arc(CC(A),CC(B))[0],
                hyperbolic_arc(CC(B),CC(C))[0],
                hyperbolic_arc(CC(C),CC(A))[0]]
    else:
        for x in P:
            if x.A.imag() >= 10000:
                ### these have to be treated specially..
                ### We truncate to the maximum y height
                ### and set the x-coordinate to the same as the other endpoint.
                A = CC(x.B.real(),ymax)
                B = CC(x.B.real(),x.B.imag())
                paths.append(hyperbolic_arc(A,B)[0])
                xmin = x.B.real()
            elif x.B.imag() >= 10000:
                ### these have to be treated specially..
                ### We truncate to the maximum y height
                ### and set the x-coordinate to the same as the other endpoint.
                A = CC(x.A.real(),x.A.imag())
                B = CC(x.A.real(),ymax)
                paths.append(hyperbolic_arc(A,B)[0])
                xmax = x.A.real()
            else:
                paths.append(x)
        ## Add a 'closing' path between the two vertical sides.
        paths.append(hyperbolic_arc(CC(xmin,ymax),CC(xmax,ymax))[0])
        ## first find the left most:
        ## an arc has a .A and .B
        As = [x.A for x in P]
        minA = min(As)
        mini = As.index(minA)
        new_paths = [P[mini]]
        paths.remove(P[mini])
        eps = 1e-15
        current = P[mini].B
        while paths != []:
            current = new_paths[-1].B
            #print "current=",current
            try:
                for p in paths:
                    #print "p.A=",p.A
                    #print "p.B=",p.B
                    if abs(p.A-current)<eps:
                        #print "appending p"
                        new_paths.append(p)
                        paths.remove(p)
                        raise StopIteration                    
                    elif abs(p.B-current)<eps: ## we reverse it 
                        #print "appending p reversed"
                        pnew = hyperbolic_arc(p.B,p.A)
                        new_paths.append(pnew[0])
                        paths.remove(p)
                        raise StopIteration
                    else:
                        continue

                print(paths)
                raise ArithmeticError("Could not connect from {0}".format(current))
            except StopIteration:
                pass
    ### new paths is now a list of hyperbolic arcs
    res = []
    import matplotlib.patches as patches
    from matplotlib.path import Path
    import numpy
    i = 0
    vertices = []
    codes = []
    for p in new_paths:
        new_codes = p.codes
        if i > 0:
            new_codes[0]=2
        i+=1
        if i==len(new_paths):
            new_codes[-1] = 79
        for v in p.vertices:
            vertices.append(v)
        codes = codes + new_codes
        #pt = patches.Path(p.vertices,codes)
        #res.append(pt)
    vertices = numpy.array(vertices, float)
    res = path = Path(vertices, codes)

#    res = patches.Path.make_compound_path(*res)
    return res
        
def build_connected_path2(P,verbose=0,model='H'):
    r"""
    Takes a path P consisting of a list of segments (not necessarily in correct order) and constructs a connected path.
    
    """
    ## This is the only place where we allow for coddes to be = 1 (i.e. starting a path)
    if len(P)==0:
        return P
    print("P=",P)
    print("P0=",P[0])
    ## Begin by locating left most path
    if model == 'H':
        xmin = P[0].vertices.min()
    else:
        xmin = None
    imin = 0
    for i in range(0,len(P)):
        try:
            if model=='H':
                testmin = P[i].vertices.min()
            else:
                v = P[i].vertices[0]
                w = (CC(v[1],-v[0]) - CC(0,1))/(CC(v[0],v[1]) - CC(1,0))
                testmin = w.real()
                for v in P[i].vertices:
                    if v[1]==0 and v[0]==1.0: # correspond to infinity
                        continue 
                    w = (CC(v[1],-v[0]) - CC(0,1))/(CC(v[0],v[1]) - CC(1,0))
                    if w.real()< testmin:
                        testmin = w.real()
                if verbose>0:
                    print("testmin[{0}]={1}".format(P[i],testmin))
            if xmin==None:
                xmin = testmin
            if testmin < xmin:
                imin = i
                xmin = testmin
            if verbose>0:
                print("xmin=",xmin)
        except AttributeError:
            pass
    codes = deepcopy(P[imin].codes)
    for j in range(len(codes)):
        if codes[j]==10:
            codes[j]=4
#        if codes[j]==79:
#            codes[j]=0
    pt = patches.Path(P[imin].vertices,codes)
    res = [pt]
    used = [imin]
    current = imin
    if verbose>0:
        print("left most path is nr. {0}".format(imin))
    ii = 0
    x,y = P[imin].A
    eps = 1e-12
    fc='orange'
    lw=2
    while len(res)<len(P) and ii<len(P):
        ii+=1
        if verbose>0:
            print("ii=",ii)
            print("current=",current)
        if ii == len(P):
            x1,y1=P[current].A
            if verbose > 0:
                print("x,y=",x,y)
                print("x1,y1=",x1,y1)
            if x==x1 and y==y1:
                x,y=P[current].B
                ii = 0
        if verbose>0:
            print("x,y=",x,y)
            print("res=",res)
        for j in range(len(P)):
            if verbose>0:
                print("++++++++++++++++++++++++++++++++++++++j=",j)
            if j in used:
                continue
            if not hasattr(P[j],"vertices"):
                used.append(j)
                continue
            xx,yy = P[j].A #vertices[0]
            #if verbose>0:p
            #    print "xx,yy(0)=",xx,yy
            #if abs(xx-x) > eps or abs(yy-y) > eps:
            #    xx,yy = P[j].B
            if abs(xx-x)<eps and abs(yy-y)<eps:
                if verbose>0:
                    print("connect prevsious segment with {0}".format(j))
                    print("segment =",P[j])
                #print "xx=x,yy=y"
                vertices = P[j].vertices; codes = deepcopy(P[j].codes)
                if not isinstance(vertices,list):
                    vertices = vertices.tolist()
                if not isinstance(codes,list):
                    codes = codes.tolist()
                for jjj in range(len(codes)):
                    if codes[jjj]==10:
                        if j<3:
                            codes[jjj]=2
                        else:
                            codes[jjj]=4                        
                    #print "vertex=",vertices[jjj]
                    if model=='D' and abs(vertices[jjj][0]+0.34)<1e-2 and abs(vertices[jjj][1]+0.52)<1e-2:
                        print("==========Changing to 2")
                        codes[jjj]=2

                if ii < len(P) and model=='H':
                    vertices.pop(); codes.pop()
                if verbose>0 and j==3:
                    print("vertices=",vertices)
                    print("codes=",codes)
                pt = patches.Path(vertices,codes)
                res.append(pt)
                used.append(j)
                current = j
                x,y = P[j].A #vertices[-1]
                break
            
            xx,yy = P[j].B #vertices[-1]
            #if verbose>0:
            #    print "xx,yy(-1)=",xx,yy
            if abs(xx-x)<eps and abs(yy-y)<eps:
                if verbose>0:
                    print("connect previous segment with {0} reversed".format(j))                
                    print("segment =",P[j])
                vertices = P[j].vertices; codes = deepcopy(P[j].codes)
                if not isinstance(vertices,list):
                    vertices = vertices.tolist()
                vertices.reverse()
                if not isinstance(codes,list):
                    codes = codes.tolist()
                if codes[0]!=1 and codes[-1]==1:
                    codes.reverse()
                for jjj in range(len(codes)):
                    if codes[jjj]==10:
                        if j==2 and (jjj==7 or jjj==8):
                            #print "THIS IS HERE!"
                            codes[jjj]=3
                        else:
                            codes[jjj]=3
                if verbose>0 and j==2:
                    print("vertices=",vertices)
                    print("codes=",codes)
                    #for jjj in range(len(vertices)-1):
                    #    print "dy/dx({0},{1})={2}".format(jjj,jjj+1,(vertices[jjj][1]-vertices[jjj+1][1])/(vertices[jjj][0]-vertices[jjj+1][0]))
#                    if model=='D' and abs(vertices[jjj][0]+0.098)<1e-3 and abs(vertices[jjj][1]+0.328)<1e-2:
#                        print "Changing to 2"
#                        codes[jjj]=2
                if codes[0]!=1:
                    raise ValueError("Startpoint of path needs a code 1. got {0}".format(codes[0]))
                #print "codes[0]=",codes[0]
                #print "codes[-1]=",codes[-1]
                if ii < len(P) and model=='H':
                    vertices.pop(); codes.pop()
                pt = patches.Path(vertices,codes)
                res.append(pt)
                used.append(j)
                current = j
                x,y = P[j].vertices[0]
                break
    if len(used) < len(P):
        print("Could not connect all paths!")
    res = patches.Path.make_compound_path(*res)
    # Now we have to make this into a closed path.
    #codes = [max(res.codes[j],13) for j in range(len(res.codes))]
    for i in range(len(res.codes)):
        if res.codes[i]==1 and i>0:
            res.codes[i]=2
    if '1.3' in matplotlib.__version__:
        res.codes[-1]=13
    else:
        res.codes[-1] = patches.Path.CLOSEPOLY
    if verbose>0:
        print("vertices=",res.vertices)
        print("codes=",res.codes)
        
    #return patches.PathPatch(res)
    #res.codes = codes
    return res
        
