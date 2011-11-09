/*
 * rankbound.cc
 *
 * Code to analytically bound the rank of an elliptic curve (assuming BSD
 * + RH). More specifically, this program computes
 *
 * sum f(\gamma)
 *
 * where 1/2 + i\gamma runs over the nontrivial zeros of the L function
 * of an elliptic curve (using the "analytic" normalization here). We use
 *
 * f = (sin(delta pi x)/(delta pi x))^2
 *
 * so that f(0) = 1 and f(x) > 0 if x is real. So if all of the
 * gammas are real, then the resulting sum is an upper bound for
 * the number of zeros of L(s, E) at s = 1/2.
 *
 */


/*
 * Basic workings of the program:
 *
 * We compute the sum over zeros using the explicit formula, which does not
 * require us to actually find the zeros. The main() function initializes
 * some variables, computes some quantities that are independent of the
 * curve or that only rely on the conductor, and then creates the curve
 * as a smalljac object.
 *
 * We register a callback with smalljac to compute successive terms in the
 * "sum over primes" part of the explicit formula as smalljac computes
 * each ap for us. Since smalljac doesn't compute ap when E has bad reduction
 * at p, these p and ap need to be specified in advance.
 *
 */


#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include "precomputation.h"
#include "mathlib.h"

#include <complex>
using namespace std;

extern "C" {
#include "smalljac.h"
}

// We are only going to be doing one curve at a time,
// so we don't mind making the following global variables.

long * bad_primes;
int * bad_prime_aps;
int number_of_bad_primes;

double delta = 3.6; 
double endpoint;
double sum;
long count;
const complex<double> I = complex<double>(0.0, 1.0);

int verbose = 0;

inline complex<double> digamma(complex<double> z) {
    return log_GAMMA(z, 1);
}

int smalljac_callback(smalljac_Qcurve_t curve, unsigned long p, int good, long a[], int n, void *arg) {
    count++;

    double p_endpoint = log(endpoint)/log(p);
    if(!good) {
        long ap = 2;
        for(int n = 0; n < number_of_bad_primes; n++) {
            if(p == bad_primes[n])
                ap = bad_prime_aps[n];
            continue;
        }
        if(ap == 2) {
            cerr << "warning: missing or bad ap for bad prime " << p << endl;
        }
        for(int k = 1; k <= p_endpoint; k++) {
            sum -= log(p) * ap/pow((double)p, k) * triangle(k * log(p)/(2.0 * M_PI), delta)/M_PI;
        }
    }
    else {
        double ap = -a[0]/sqrt(p);
        complex<double> alpha = ap/2.0 + I * sqrt(1 - ap/2.0);
        complex<double> beta = ap/2.0 - I * sqrt(1 - ap/2.0);

        
        for(int k = 1; k <= p_endpoint; k++) {
            double cn = (pow(alpha, k) + pow(beta, k)).real();
            double z = log(p) * ((pow(alpha, k) + pow(beta, k)).real()/pow(p, k/2.0) *
                                                triangle(k * log(p)/(2.0 * M_PI), delta)/M_PI);
            sum -= z;
        }
        
        if( verbose && ((count % 100000 == 0) || (count < 10000 && count % 10 == 0)))
            cout << p << " " << (double)p/(double)endpoint << " " << sum << " " << endl;
    }
    return 1;
}

extern "C"
double rankbound(char * curve_string, double logN, long * bad_primes_, int * ap_, int len_bad_primes_, double delta_, int verbose_) {
    verbose = verbose_;
    delta = delta_;
    bad_primes = bad_primes_;
    bad_prime_aps = ap_;
    number_of_bad_primes = len_bad_primes_;

    endpoint = exp(2 * M_PI * delta);
    double conductor_term = triangle(0, delta)/(2 * M_PI) * logN;
    sum = gamma_terms(delta) - triangle(0, delta) * log(M_PI)/M_PI + conductor_term;

    smalljac_Qcurve_t curve;

    long result;
    int err;
    long maxp = endpoint;

    count = 0;

    curve = smalljac_Qcurve_init(curve_string, &err);

    result = smalljac_Lpolys(curve, 1, maxp, 0, smalljac_callback, (void * )NULL);
    smalljac_Qcurve_clear(curve);
    
    return sum;
}

int main(int argc, char *argv[]) {
    if(argc < 3) {
        cout << "Usage: rankbound input_file delta" << endl;
        cout << "Example: rankbound rank20_example 2.5" << endl;
        return 1;
    }

    double delta;

    delta = strtod(argv[2], NULL);

    cout << setprecision(17);
    
    ifstream infile(argv[1]);

    //
    // we expect infile to contain one line, which looks like
    //
    // [a1,a2,a3,a4,a6] logN p ap p ap ...
    //

    string curve_string;
    infile >> curve_string;

    double logN;
    infile >> logN;

    vector<long> bad_primes;
    vector<int> bad_prime_aps;

    while(!infile.eof()) {
        unsigned long p;
        int ap;
        infile >> p;
        bad_primes.push_back(p);
        infile >> ap;
        bad_prime_aps.push_back(ap);

    }

    verbose = 1;

    double bound = rankbound(const_cast<char *>(curve_string.c_str()),
                             logN,
                             &bad_primes[0],
                             &bad_prime_aps[0],
                             bad_primes.size(),
                             delta,
                             verbose);

    cout << curve_string << " " << bound << endl;


}
