#include <iostream>
#include <cmath>
#include <complex>
#include <functional>
#include <gsl/gsl_integration.h>
#include <map>

#include <libInterpolate/Interpolate.hpp>

#include "filereader.h"

using namespace std::complex_literals;

double MPI = 0.14;
double GeVtoIfm = 5.0677;
double IfmtoGeV = 1 / GeVtoIfm;
double fmtoIGeV = GeVtoIfm;
double IGeVtofm = 1 / fmtoIGeV;

_1D::CubicSplineInterpolator<double> tau,r,Dtau,Dr,ur,utau;

double w(double p) { return sqrt(p * p + MPI * MPI); }
double J0rp(double alpha, double p) { return std::cyl_bessel_j(0, r(alpha) * p * GeVtoIfm); }
double J1rp(double alpha, double p) { return std::cyl_bessel_j(1, r(alpha) * p * GeVtoIfm); }
double Y0tw(double alpha, double p) { return std::cyl_neumann(0, tau(alpha) * w(p) * GeVtoIfm); }
double Y1tw(double alpha, double p) { return std::cyl_neumann(1, tau(alpha) * w(p) * GeVtoIfm); }
double J0tw(double alpha, double p) { return std::cyl_bessel_j(0, tau(alpha) * w(p) * GeVtoIfm); }
double J1tw(double alpha, double p) { return std::cyl_bessel_j(1, tau(alpha) * w(p) * GeVtoIfm); }

std::complex<double> H1(double alpha, double p) { return J0rp(alpha, p) * (-Y0tw(alpha, p) + 1i * J0tw(alpha, p)); }
std::complex<double> H2(double alpha, double p)
{
    return J1rp(alpha, p) * (-Y0tw(alpha, p) + 1i * J0tw(alpha, p)) * Dtau(alpha) * fmtoIGeV * p +
           J0rp(alpha, p) * (-Y1tw(alpha, p) + 1i * J1tw(alpha, p)) * Dr(alpha) * fmtoIGeV * w(p);
}

struct args
{
    double p;
    std::function<double(double)> func;
    std::function<double(double)> Dfunc;
};
std::complex<double> spectr(double p, std::function<double(double)> func, std::function<double(double)> Dfunc)
{
    auto integrand_re = [](double alpha, void *params)
    {
        args myargs = *(struct args *)params;
        return tau(alpha) * r(alpha) * fmtoIGeV * fmtoIGeV * (myargs.Dfunc(alpha) * H1(alpha, myargs.p).real() + myargs.func(alpha) * H2(alpha, myargs.p).real());
    };
    auto integrand_im = [](double alpha, void *params)
    {
        args myargs = *(struct args *)params;
        return tau(alpha) * r(alpha) * fmtoIGeV * fmtoIGeV * (myargs.Dfunc(alpha) * H1(alpha, myargs.p).imag() + myargs.func(alpha) * H2(alpha, myargs.p).imag());
    };

    args myargs = {p, func, Dfunc};

    gsl_function F_re, F_im;
    F_re.function = integrand_re;
    F_im.function = integrand_im;
    F_re.params = &myargs;
    F_im.params = &myargs;

    double result_re, result_im, error_re, error_im;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(10000);

    gsl_integration_qags(&F_re, 0, M_PI_2, 0, 1e-3, 10000, w, &result_re, &error_re);
    gsl_integration_workspace_free(w);

    w = gsl_integration_workspace_alloc(1000);
    gsl_integration_qags(&F_im, 0, M_PI_2, 0, 1e-3, 10000, w, &result_im, &error_im);

    return result_re + 1i * result_im;
}

int main()
{
    // READ IN FREEZOUT DATA
    csvdata mydata = readcsv("./../../Mathematica/data/ExampleFreezeOut.csv", ",", true, true);
    std::ofstream debug("debug.txt");
    for (int i = 0; i < mydata.header.size(); i++)
    {
        debug << mydata.header[i] << (i == mydata.header.size() - 1 ? "" : ", ") << std::endl;
        std::cout << mydata.header[i] << (i == mydata.header.size() - 1 ? "" : ", ") << std::endl;
        for (int j = 0; j < mydata.data[i].size(); j++)
            debug << mydata.data[i][j] << std::endl;
        debug << std::endl;
    }
    debug.close();

    // DEFINE INTERPOLATING FUNCTIONS
    tau.setData(mydata.data[0],mydata.data[1]);
    r.setData(mydata.data[0],mydata.data[2]);
    Dtau.setData(mydata.data[0],mydata.data[3]);
    Dr.setData(mydata.data[0],mydata.data[4]);
    ur.setData(mydata.data[0],mydata.data[5]);
    utau.setData(mydata.data[0],mydata.data[6]);

    // DEFINE PION FIELD AND DERIVATIVE ON FREEZOUT SURFACE
    std::function<double(double)> func = [](double alpha)
    { return 1; };
    std::function<double(double)> Dfunc = [](double alpha)
    { return 0; };

    // COMPUTE SPECTRUM AT P = 1
    spectr(1, func, Dfunc);

}