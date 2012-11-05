r"""
Cython wrapper of K-Bessel routine written in C.

Authors:
Andrew Booker (original code)
Fredrik Stroemberg (wrapper)


"""

#include <stdio.h>
#include <math.h>
cdef extern from "math.h":
    double fabs(double)
    double fmax(double,double)
    int ceil(double) 
    double sqrt(double)
    double exp(double)
    #double sin(double)
    #double cos(double)
    double power(double,double)

cdef extern from "complex.h":
    ### "Hack" suggested by Robert Bradshaw in 2009.
    ### TODO: Is there now a better way for pure c double complex?
    ctypedef double cdouble "double complex"
    cdef double creal(cdouble)
    cdef double cimag(cdouble)
    cdef cdouble _Complex_I
    cdef double carg(cdouble)
    cdef double cabs(cdouble)
    cdef cdouble cexp(cdouble)
    cdef cdouble cpow(cdouble,cdouble)
    
cdef float ystep =  0.125
cdef int nsamples= 320
cdef int ntaylor = 16
cdef int nasymp  = 14
cdef int wdegree = 19
cdef int cdegree = 27
cdef float ymin = 0.5
cdef float ymax
ymax = (ymin+nsamples*ystep)

#typedef struct {
#       double re,im;
#} complex;

cdef complex C
cdef double R,R2,asymp[nasymp],taylor[nsamples][ntaylor]
#if 0
#extern struct { long double re,im; } Cdata[cdegree+1];
#extern long double Wdata[nsamples][2][wdegree+1];
#else
#include "wdata"
cdef extern from "wdata":
#static struct {long double re,im} Cdata
    ctypedef struct ccomplex:  #double cdouble "double complex"
        pass
    cdef ccomplex Cdata[cdegree+1]
    long double Wdata[nsamples][2][wdegree+1];
    #cdef struct {long double re,im} 
#endif

cdef double W(double y):
    cdef double temp,s
    cdef complex z
    cdef int i,j
    if y <= ymin:
        if (R < 1.0e-16):
            z.re = 1.0
            z.im = log(y)
        else:
            temp = log(y)*R
            #/* z.im will always be the actual imaginary part divided by R */
            #/* note we don\'t lose accuracy here, thanks to floating point */
            z.re = cos(temp)
            z.im = sin(temp)/R
            
        temp = sqrt(y);
        z.re *= temp; z.im *= temp;
        temp = z.re*C.re-z.im*C.im*R2;
        z.im = z.re*C.im+z.im*C.re;
        z.re = temp;
        y *= y*0.25;
        s = z.im;
        j = 1;
        while (fabs(z.re)+fabs(z.im) > 1.0e-17):
            temp = y/(j*(j*j+R2));
            z.re *= temp; z.im *= temp;
            temp = z.re*j+z.im*R2;
            z.im = z.im*j-z.re;
            z.re = temp;
            s += z.im;
            j+=1;
        return s;
    if (y < ymax):
        i = (int)((y-ymin)/ystep);
        temp = y-((i+0.5)*ystep+ymin);
        j = ntaylor-1;
        s = taylor[i][j];
        while (--j >= 0):
            s = s*temp+taylor[i][j];
        return s;

    temp = 1.0/y
    s = asymp[nasymp-1];
    #for (j=nasymp-2;j>=0;j--):
    for j from nasymp-2>=j >=0:
        s = s*temp+asymp[j];
        j-=1
    return s*exp(-y);

cdef long double halfpi    = 1.5707963267948966192313217
cdef long double sqrthalfpi= 1.2533141373155002512078826
cdef double set_r = 0
cdef void Winit(double r):
    #cdef long double lambda,s,t1,t2,temp[ntaylor],w[ntaylor]
    cdef long double xlambda,s,t1,t2,temp[ntaylor],w[ntaylor]
    cdef int i,j,l
    if set_r == r:
        return # already inited
    set_r = r
    R=fabs(r); t1=<long double>(R*R);
    R2=<double>t1; xlambda=t1+<long double>0.25;
    t1 = Cdata[cdegree].re; t2 = Cdata[cdegree].im;
    #for (j=cdegree-1;j>=0;j--):
    for j from cdegree-1>=j>=0: #(j=cdegree-1;j>=0;j--):
        t1 = t1*R+Cdata[j].re;
        t2 = t2*R+Cdata[j].im;
        j-=1
    C.re = <double>t1, C.im = <double>t2;
    for i from 0<=i<nsamples: #(i=0;i<nsamples;i++):
        for l from 0<=l<2: #(l=0;l<2;l++):
            s = Wdata[i][l][wdegree]
            for j from wdegre>=j>=0: #(j=wdegree-1;j>=0;j--):
                s = s*<long double>R+Wdata[i][l][j]
                j-=1
            w[l] = s;
            l+=1
        i+=1
        t1 = -<long double>1.0/((i+<long double>0.5)*ystep+ymin) #; /* -1/y */
        #for (j=1,t2=xlambda*t1;j<ntaylor;j++,t2*=t1):
        t2=xlambda*t1
        for j from 1<=j<ntaylor:
            temp[j]=(j-1)*t2;
            j+=1
            t2*=t1
        for l from 2<=l<ntaylor by -1:
            #(;l<ntaylor;l++):
            #for (j=l-2,s=w[j];j>=0;j--):
            s=w[l-2]
            for j from l-1>=j>=0 by -1: # (j=l-2,s=w[j];j>=0;j--):
                s -= temp[l-j]*w[j]
                #j-=1
            w[l] = s/(l*l-l);
            
        #for (l=0;l<ntaylor;l++):
        for l from 0<=l<ntaylor: #(l=0;l<ntaylor;l++):
            taylor[i][l] = <double>w[l];
        t1 = sqrthalfpi*exp(halfpi*<long double>R);
        for j from 0<=j<nasymp: #(j=0;j<nasymp;j++):
            asymp[j] = <double>t1;
            t1 *= -(xlambda+j*(j+1))/(2*j+2);

#cdef int main(void):
cdef double kbessel2(r,y):
    if r<>set_r:
        Winit(r)
    return W(y)
    ## cdef char buf[256]
    ## cdef unsigned int last=0,temp1,temp2
    ## cdef int i,imax=0
    ## cdef double r,y,w,w2
    ## while(fgets(buf,sizeof(buf),stdin)):
    ##     sscanf(buf,"%u, %u, %lf",&temp1,&temp2,&w)
    ##     y = <double>(temp2+1)*50.0/(1LL<<31)
    ##     if (temp1 != last):
    ##         r = <double>temp1/(1LL<<31)
    ##         Winit(r)
    ##         last = temp1
    ##     w2 = W(y)
    ##     #i = *(int *)&w-*(int *)&w2
    ##     i = (<int *>&w)[0]-(<int *>&w2)[0]
    ##     if (i < 0): i = -i
    ##     if (i > imax):
    ##         print "r=%.15f, y=%.15f, w=%g, error=%d\n".format(r,y,w,i)
    ##         imax = i
    ## return 0
