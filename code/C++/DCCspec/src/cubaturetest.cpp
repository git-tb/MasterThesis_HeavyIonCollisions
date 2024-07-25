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

#include <algorithm>
#include <iterator>

#include "savedata.h"       // write to file
#include "processdata.h"    // initial conditions and freezeout functions
#include "constants.h"      // pion mass, GeV to inverse fm
#include "cubature.h"

int main(int ac, char* av[])
{
    // auto f = [](
    //     unsigned ndim,
    //     const double *x,
    //     void *fdata,
    //     unsigned fdim,
    //     double *fval)
    // {
    //     fval[0] = exp(-x[0]*x[0]-x[1]*x[1]);
    //     return 0;
    // };

    // double xmin[2] = {-10, -10}, xmax[2] = {10,10};
    // int XDIM(2), FDIM(1), MAXEVAL(10000);
    // double ERRABS(1e-7), ERRREL(1e-7);
    // double val, err;
    // void* fdata = NULL;
    // hcubature(FDIM, f, fdata, XDIM, xmin, xmax, MAXEVAL, ERRABS, ERRREL, ERROR_INDIVIDUAL, &val, &err);
    // printf("Computed integral = %0.10g +/- %g\n", val, err);

    auto f = [](
        unsigned ndim,
        const double *x,
        void *fdata,
        unsigned fdim,
        double *fval)
    {
        fval[0] = 2/sqrt(1-x[0]*x[0]);
        return 0;
    };

    double xmin[1] = {-1}, xmax[1] = {1};
    int XDIM(1), FDIM(1), MAXEVAL(10000);
    double ERRABS(1e-7), ERRREL(1e-7);
    double val, err;
    void* fdata = NULL;
    hcubature(FDIM, f, fdata, XDIM, xmin, xmax, MAXEVAL, ERRABS, ERRREL, ERROR_INDIVIDUAL, &val, &err);
    printf("Computed integral = %0.10g +/- %g\n", val, err);
}

/* EXAMPLE CALL TO CUBATURE LIBRARY

int f(unsigned ndim, const double *x, void *fdata,
      unsigned fdim, double *fval);

int hcubature(unsigned fdim, integrand f, void *fdata,
              unsigned dim, const double *xmin, const double *xmax, 
              size_t maxEval, double reqAbsError, double reqRelError,
              error_norm norm,
              double *val, double *err);

int f(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {
    double sigma = *((double *) fdata); // we can pass σ via fdata argument
    double sum = 0;
    unsigned i;
    for (i = 0; i < ndim; ++i) sum += x[i] * x[i];
    // compute the output value: note that fdim should == 1 from below
    fval[0] = exp(-sigma * sum);
    return 0; // success*
}

double xmin[3] = {-2,-2,-2}, xmax[3] = {2,2,2}, sigma = 0.5, val, err;
hcubature(1, f, &sigma, 3, xmin, xmax, 0, 0, 1e-4, ERROR_INDIVIDUAL, &val, &err);
printf("Computed integral = %0.10g +/- %g\n", val, err);

Computed integral = 13.69609043 +/- 0.00136919

*/