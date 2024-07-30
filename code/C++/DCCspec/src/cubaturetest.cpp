#include <iostream>
#include <cmath>
#include <complex>
#include <functional>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <map>
#include <boost/math/interpolators/pchip.hpp>
#include <boost/program_options.hpp>
#include <filesystem>
#include <ctime>   // std::time, std::local_time
#include <iomanip> // std::put_time
#include <chrono>

#include <algorithm>
#include <iterator>

#include "savedata.h"       // write to file
#include "processdata.h"    // initial conditions and freezeout functions
#include "constants.h"      // pion mass, GeV to inverse fm
#include "cubature.h"
#include "cuba.h"

double func(double x, double y)
{
    double u = (1 + x*y)/sqrt(2+x*x);

    if(u <= 1)  return 0;
    else        return 1/(sqrt(u*u-1) * sqrt(1-y*y));
}

const double XMIN(0), XMAX(1),YMIN(-1),YMAX(1);
static int Integrand(const int *ndim, const double xx[],
  const int *ncomp, double ff[], void *userdata) {
    double x = XMIN + (XMAX - XMIN) * xx[0];
    double y = YMIN + (YMAX - YMIN) * xx[1];

    ff[0] = (XMAX - XMIN) * (YMAX - YMIN) * func(x,y);

    return 0;
}

int main(int ac, char* av[])
{
    auto start = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    // CUBA LIBRARY

    #define NDIM 2
    #define NCOMP 1
    #define USERDATA NULL
    #define NVEC 1
    #define EPSREL 1e-3
    #define EPSABS 1e-7
    #define VERBOSE 0 //2
    #define LAST 4
    #define SEED 0
    #define MINEVAL 0
    #define MAXEVAL 1e6

    #define NSTART 1000
    #define NINCREASE 500
    #define NBATCH 1000
    #define GRIDNO 0
    #define STATEFILE NULL
    #define SPIN NULL

    #define NNEW 1000
    #define NMIN 2
    #define FLATNESS 25.

    #define KEY1 47
    #define KEY2 1
    #define KEY3 1
    #define MAXPASS 5
    #define BORDER 0.
    #define MAXCHISQ 10.
    #define MINDEVIATION .25
    #define NGIVEN 0
    #define LDXGIVEN NDIM
    #define NEXTRA 0

    #define KEY 0

    int comp, nregions, neval, fail;
    double integral[NCOMP], error[NCOMP], prob[NCOMP];

    printf("-------------------- Vegas test --------------------\n");

    start = std::chrono::steady_clock::now();
    Vegas(NDIM, NCOMP, Integrand, USERDATA, NVEC,
        EPSREL, EPSABS, VERBOSE, SEED,
        MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
        GRIDNO, STATEFILE, SPIN,
        &neval, &fail, integral, error, prob);
    end = std::chrono::steady_clock::now();
    elapsed_seconds = end - start;

    printf("VEGAS RESULT:\tneval %d\tfail %d\n",
        neval, fail);
    for( comp = 0; comp < NCOMP; ++comp )
        printf("VEGAS RESULT:\t%.8f +- %.8f\tp = %.3f\n%.5f seconds\n",
        (double)integral[comp], (double)error[comp], (double)prob[comp],elapsed_seconds.count());
        
    printf("\n-------------------- Suave test --------------------\n");

    start = std::chrono::steady_clock::now();
    Suave(NDIM, NCOMP, Integrand, USERDATA, NVEC,
        EPSREL, EPSABS, VERBOSE | LAST, SEED,
        MINEVAL, MAXEVAL, NNEW, NMIN, FLATNESS,
        STATEFILE, SPIN,
        &nregions, &neval, &fail, integral, error, prob);
    end = std::chrono::steady_clock::now();
    elapsed_seconds = end - start;

    printf("SUAVE RESULT:\tnregions %d\tneval %d\tfail %d\n",
        nregions, neval, fail);
    for( comp = 0; comp < NCOMP; ++comp )
        printf("SUAVE RESULT:\t%.8f +- %.8f\tp = %.3f\n%.5f seconds\n",
        (double)integral[comp], (double)error[comp], (double)prob[comp],elapsed_seconds.count());

    printf("\n------------------- Divonne test -------------------\n");

    start = std::chrono::steady_clock::now();
    Divonne(NDIM, NCOMP, Integrand, USERDATA, NVEC,
        EPSREL, EPSABS, VERBOSE, SEED,
        MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
        BORDER, MAXCHISQ, MINDEVIATION,
        NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
        STATEFILE, SPIN,
        &nregions, &neval, &fail, integral, error, prob);
    end = std::chrono::steady_clock::now();
    elapsed_seconds = end - start;

    printf("DIVONNE RESULT:\tnregions %d\tneval %d\tfail %d\n",
        nregions, neval, fail);
    for( comp = 0; comp < NCOMP; ++comp )
        printf("DIVONNE RESULT:\t%.8f +- %.8f\tp = %.3f\n%.5f seconds\n",
        (double)integral[comp], (double)error[comp], (double)prob[comp],elapsed_seconds.count());

    printf("\n-------------------- Cuhre test --------------------\n");

    start = std::chrono::steady_clock::now();
    Cuhre(NDIM, NCOMP, Integrand, USERDATA, NVEC,
        EPSREL, EPSABS, VERBOSE | LAST,
        MINEVAL, MAXEVAL, KEY,
        STATEFILE, SPIN,
        &nregions, &neval, &fail, integral, error, prob);
    end = std::chrono::steady_clock::now();
    elapsed_seconds = end - start;

    printf("CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n",
        nregions, neval, fail);
    for( comp = 0; comp < NCOMP; ++comp )
        printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n%.5f seconds\n",
        (double)integral[comp], (double)error[comp], (double)prob[comp],elapsed_seconds.count());

    // CUBATURE LIBRARY

    printf("\n-------------------- Cubature test --------------------\n");

    auto f = [](
        unsigned ndim,
        const double *xx,
        void *fdata,
        unsigned fdim,
        double *fval)
    {
        double x = xx[0];
        double y = xx[1];
        fval[0] = func(x,y);
        return 0;
    };

    double xmin[2] = {XMIN,YMIN}, xmax[2] = {XMAX,YMAX};
    int XDIM(2), FDIM(1);
    double val, err;

    start = std::chrono::steady_clock::now();
    hcubature(FDIM, f, USERDATA, NDIM, xmin, xmax, MAXEVAL, EPSABS, EPSREL, ERROR_INDIVIDUAL, &val, &err);
    end = std::chrono::steady_clock::now();
    elapsed_seconds = end - start;

    printf("(H)CUBATURE RESULT:\t%.8f +- %.8f\n%.5f seconds\n",
        val, err,elapsed_seconds.count());
}