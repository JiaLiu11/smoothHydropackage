// Ver 1.1
// Note that all calculations are done at a given particle rapidity y; and all
// "y_minus_y_minus_eta_s" appearences in the code are y-y_minus_eta_s.

#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "arsenal.h"
#include "GammaDistribution.h"
#include "Stopwatch.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

int main()
{
    double alpha = 1.0;
    double beta = 0.5;
    GammaDistribution test(alpha, beta);
    Stopwatch sw;
    long sum;
    
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    ofstream check("checkGammaDistribution.dat");
    for (long i=1; i<=1000000; i++)
    {
        check << gsl_ran_gamma(r, alpha, 2.0) << endl;
    }
    sw.toc();

    gsl_rng_free (r);

    // sw.tic();
    // sum = 0;
    // for (long i=1; i<=10000; i++) sum += binomial_coefficient(i, i/2);
    // cout << "sum=" << sum << endl;
    // sw.toc();
    // cout << "1 takes time: " << sw.takeTime() << endl;

    // sw.tic();
    // sum = 0;
    // for (long i=1; i<=10000; i++) sum += 1.0/(i+1)/beta_function(i-i/2+1, i/2+1);
    // cout << "sum=" << sum << endl;
    // sw.toc();
    // cout << "2 takes time: " << sw.takeTime() << endl;


    // for (long i=0; i<10; i++) cout << i << "  " << nbd.pdf(i) << endl;
    // for (double i=0; i<10; i+=0.25) cout << setw(10) << i << "  " << nbd.envelopPdfTab->map(i) << endl;
    //for (double i=0; i<0.1; i+=0.003641) cout << setw(10) << i << "  " << nbd.envelopInvCDFTab->map(i) << endl;
    // nbd.envelopPdfTab->printFunction();
    
    //ofstream check("checkGammaDistribution.dat");
    //for (long i=1; i<=1000000; i++)
    //{
    //    check << test.rand() << endl;
    //}
    //sw.toc();
    cout << "takes time: " << sw.takeTime() << endl;

}


 
