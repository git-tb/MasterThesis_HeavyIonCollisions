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

double Q(1.0);
double ma(1), mb(1), mc(1);
double p_abc = 1/(2*ma) * sqrt(
    (pow(ma+mb,2)-pow(mc,2))*(pow(ma-mb,2)-pow(mc,2))
);
double E_abc(sqrt(mb*mb + p_abc*p_abc));

std::function<double(double,double)> w = [](double p, double m) { return sqrt(pow(p,2)+pow(m,2));};
std::function<double(double, double,double,double)> absgradg = [](double t, double eta, double phi, double p) {
    return sqrt(
        pow(-w(p,mb)*Q*Q*t*cosh(eta)/w(Q*t,ma) + Q*p*cos(phi),2)+
        pow(w(Q*t,ma)*w(p,mb)*sinh(eta),2)+
        pow(t*Q*p*sin(phi),2)
    );
};
std::function<double(double,double,double)> u_star = [](double t, double phi, double p){
    return (ma * E_abc + t*Q*p*cos(phi))/(w(Q*t,ma) * w(p,mb));
};
std::function<double(double, double,double)> areafactor = [](double t, double phi, double p) {
    double u = u_star(t,phi,p);
    return sqrt(
        1/(u*u-1) * pow(Q*p/(w(t*Q,ma)*w(p,mb)),2) (
            pow(cos(phi)*(1 - pow(t*Q/w(t*Q,ma),2)),2) + pow(t*sin(phi),2)
        )
    );
};

std::function<double(double)> primespec;
std::function<double(double,double,void*)> restrictfunc = [](double t, double phi,double p){
    double u = u_star(t,v,p);
    if(u < 1) return 0;
    return t * primespec(t*Q) * areafactor(t,v,p)/absgradg(t,u_,v,p);
};

struct args
{
    double p;
};
std::function myintegrand = [](
        unsigned ndim,
        const double *x,
        void *fdata,
        unsigned fdim,
        double *fval)
{
    double t(x[0]), v(x[1]);
    args myargs = *(struct args *)fdata;
    fval[0] = restrictfunc(t,v,myargs.p);
    return 0;
};

std::vector<double> decayspec(
    std::vector<double> ps,
    std::function<double(double)> primespec)
{
    std::vector<double> finalspec(ps.size());

    double ERRABS(0), ERRREL(1e-2);
    int MAXEVAL(1000), FDIM(1), XDIM(2);
    args myargs;

    for (int i = 0; i < ps.size(); i++)
    {
        std::cout << i + 1 << " / " << ps.size() << std::endl;

        myargs.p = ps[i];
        hcubature(FDIM, f, &myargs, XDIM, xmin, xmax, MAXEVAL, ERRABS, ERRREL, ERROR_INDIVIDUAL, &val, &err);
    }
    return finalspec;
}
/* #endregion */

int main(int ac, char* av[])
{
    auto f = [](
        unsigned ndim,
        const double *x,
        void *fdata,
        unsigned fdim,
        double *fval)
    {
        fval[0] = exp(-x[0]*x[0]-x[1]*x[1]);
        return 0;
    };

    double xmin[2] = {-10, -10}, xmax[2] = {10,10};
    int XDIM(2), FDIM(1), MAXEVAL(10000);
    double ERRABS(1e-7), ERRREL(1e-7);
    double val, err;
    void* fdata = NULL;
    hcubature(FDIM, f, fdata, XDIM, xmin, xmax, MAXEVAL, ERRABS, ERRREL, ERROR_INDIVIDUAL, &val, &err);
    printf("Computed integral = %0.10g +/- %g\n", val, err);


    // /* #region COMMAND LINE OPTIONS */
    // // DEFINE VARIABLES TO BE SET
    // double pmin(0), pmax(1.0);
    // int Nps(100);
    // std::string primespec_data;

    // // DECLARE SUPPORTED OPTIONS
    // namespace po = boost::program_options;
    // po::options_description desc("Allowed options");
    // desc.add_options()
    //     ("help", "produce help message")
    //     ("pTmax", po::value<double>()->default_value(1.0), "spectrum is computed on [0,pTmax]")
    //     ("NpT", po::value<int>()->default_value(100), "number of sample points within [0,pTmax]")
    //     ("primespecpath",po::value<std::string>(),"csv file containing primary spectrum");

    // po::variables_map vm;
    // po::store(po::parse_command_line(ac, av, desc), vm);
    // po::notify(vm);

    // // PROCESS CMDLINE OPTIONS
    // if (vm.count("help"))    {
    //     std::cout << desc << "\n";
    //     return 1;
    // }
    // pmax = vm["pTmax"].as<double>();
    // Nps = vm["NpT"].as<int>();
    // primespecpath = vm["primespecpath"].as<std::string>();
    // /* #endregion */

    // // CREATE DIRECTORY TO SAVE FILES
    // std::time_t t = std::time(nullptr);
    // std::tm tm = *std::localtime(&t);
    // std::stringstream timestamp_sstr;
    // timestamp_sstr << std::put_time(&tm, "%Y%m%d_%H%M%S");
    // std::string timestamp = timestamp_sstr.str();
    // std::string pathname = "data/spec_"+timestamp;
    // std::filesystem::create_directories(pathname);
   
    
    // // DEFINE PION FIELD AND DERIVATIVE ON FREEZOUT SURFACE
    // std::function<double(double)> primespecc;

    // std::vector<double> myspectr_abs2 = decayspec(ps, primespec);

    // writeSamplesToFile(pathname+"/spectr.txt", ps, myspectr_abs2,{"pT","abs2Re","abs2Im"},{timestamp});
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