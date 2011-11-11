#include <complex>
#include <cmath>

#include <iostream>

using namespace std;

const double log_2Pi = 1.83787706640935;


inline double my_norm(complex<double> z)
{
    return(real(z*conj(z)));
}

const double bernoulli [] = {
    1.0000000000000000000000000000,
    -0.50000000000000000000000000000,
    0.16666666666666666666666666667,
    0.00000000000000000000000000000,
    -0.033333333333333333333333333333,
    0.00000000000000000000000000000,
    0.023809523809523809523809523810,
    0.00000000000000000000000000000,
    -0.033333333333333333333333333333,
    0.00000000000000000000000000000,
    0.075757575757575757575757575758,
    0.00000000000000000000000000000,
    -0.25311355311355311355311355311,
    0.00000000000000000000000000000,
    1.1666666666666666666666666667,
    0.00000000000000000000000000000,
    -7.0921568627450980392156862745,
    0.00000000000000000000000000000,
    54.971177944862155388471177945,
    0.00000000000000000000000000000,
    -529.12424242424242424242424242,
    0.00000000000000000000000000000,
    6192.1231884057971014492753623,
    0.00000000000000000000000000000,
    -86580.253113553113553113553114,
    0.00000000000000000000000000000,
    1.4255171666666666666666666667e6,
    0.00000000000000000000000000000,
    -2.7298231067816091954022988506e7,
    0.00000000000000000000000000000,
    6.0158087390064236838430386817e8,
    0.00000000000000000000000000000,
    -1.5116315767092156862745098039e10,
    0.00000000000000000000000000000,
    4.2961464306116666666666666667e11,
    0.00000000000000000000000000000,
    -1.3711655205088332772159087949e13,
    0.00000000000000000000000000000,
    4.8833231897359316666666666667e14,
    0.00000000000000000000000000000,
    -1.9296579341940068148632668145e16,
    0.00000000000000000000000000000,
    8.4169304757368261500055370986e17,
    0.00000000000000000000000000000,
    -4.0338071854059455413076811594e19,
    0.00000000000000000000000000000,
    2.1150748638081991605601453901e21,
    0.00000000000000000000000000000,
    -1.2086626522296525934602731194e23,
    0.00000000000000000000000000000
};

complex<double> log_GAMMA (complex<double> z,int n)
{
    const long DIGITS = 15;
    double M;
    complex<double> log_G,r,r2,y;

    double xx=z.real(),yy=z.imag();
    if(xx<0) xx=-xx;

    double x;
    int i,m;


    //assume the remainder stopping at the mth term is roughly bernoulli(2m)/(z+M)^(2m).
    //Now bernoulli(2m) is roughly (2m)!/(2Pi)^(2m). So remainder is more or less
    //(2m/(2ePi(z+M))^(2m). Setting 2m = Digits, we see that we need m/(ePi(z+M))
    //to be roughly 1/10 in size, i.e. |z+M| should be roughly 10*m/(ePi)=10*Digits/(2ePi).
    //squaring, gives us how large M should be.

    //n==0 leads to |z+M| >10*2m/(2*Pi*e) with 2m=Digits. However, when we differentiate
    //we end up with extra powers of (z+M) in the denominators, but extra (2m-1)...(2m+n-2)
    //in the numerators. So assume further that |z+M| > 2m+n = Digits+n

    if(n==0){
        //.343 = 100./(4*Pi*Pi*exp(2.))
        if((xx*xx+yy*yy)> .343*DIGITS*DIGITS) M=0;
        else M=ceil(sqrt((DIGITS*DIGITS*.343-yy*yy))-xx+1);
    }
    else{
        if((xx*xx+yy*yy)> .343*(DIGITS+n)*(DIGITS+n)) M=0;
        else M=ceil(sqrt(((DIGITS+n)*(DIGITS+n)*.343-yy*yy))-xx+1);
    }

    if(n==0)
       log_G=log_2Pi/2+(z+M-.5)*log(z+M)-(z+M);
    else if(n==1)
       log_G=log(z+M)-.5/(z+M);
    else{
       r=1;
       for(m=1;m<=n-1;m++){
          r=-r*(double)m/(z+M);
       }
       log_G=log_G-r/((double)n-1)-.5*r/(z+M);
    }

    r=1;
    for(m=1;m<=n;m++){
       r=-r*(double)m/(z+M);
    }
    r=r/(z+M);

    r2=1.0/((z+M)*(z+M));
    m=2;
    x=my_norm(r);
    do{
        y=bernoulli[m]*r;
        log_G=log_G+y/((double)m*(m-1));

        r=r*((double)m+n-1)*(m+(double)n)*r2/((m-1)*(double)m);
        m=m+2;
    //}while(m<=DIGITS&&my_norm(r)*tolerance_sqrd<x);
    }while(m<=DIGITS);

    for(m=0;m<=M-1;m++){
       if(n==0){
           log_G=log_G-log(z+(double)m); //XXX might be faster to multiply and take one log,
                                 //but careful about overflow errors
       }
       else{
           r=1;
           for(i=1;i<=n;i++){
              r=-r*(double)i/(z+(double)m);
           }
           log_G=log_G+r/(double)n;
       }
    }

    return log_G;
}

double sinc(double u, double Delta) {
    // return (sin(pi * u * Delta)/(pi * u * Delta))
    //
    // If the input is very small, we use the power series expansion.
    // Otherwise we will compute it directly.

    double x = u * Delta * M_PI;

    if(abs(x) < .00001) {
        return -pow(x, 6)/5040 + pow(x, 4)/120 - pow(x, 2)/6 + 1;
    }
    else {
        return sin(x)/x;
    }
}

double sinc_square(double u, double Delta) {
    // return (sin(pi * u * Delta)/(pi * u * Delta))^2
    //
    // If the input is very small, we use the power series expansion.
    // Otherwise we will compute it directly.

    double x = u * Delta * M_PI;

    if(abs(x) < .00001) {
        return 1 - 1.0/3.0 * pow(x, 2) + 2.0/45.0 * pow(x, 4) - 1.0/315.0 * pow(x, 6);
    }
    else {
        double a = sin(x)/x;
        return a*a;
    }
}

double triangle(double u, double Delta) {
    return (1 - abs(u/Delta))/Delta;
}
