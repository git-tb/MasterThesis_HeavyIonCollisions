#include <iostream>
#include <cmath>
#include <complex>
#include <functional>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <map>
#include <boost/math/interpolators/pchip.hpp>

#include "filereader.h"

using namespace std::complex_literals;

double MPI = 0.14;
double GeVtoIfm = 5.0677;
double IfmtoGeV = 1 / GeVtoIfm;
double fmtoIGeV = GeVtoIfm;
double IGeVtofm = 1 / fmtoIGeV;

std::function<double(double)> tau, r, Dtau, Dr, ur, utau;

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
    int ITERATIONS(1000);
    double EPSABS(0), EPSREL(1e-3);

    gsl_set_error_handler_off();

    gsl_integration_workspace *w_re = gsl_integration_workspace_alloc(ITERATIONS);
    int status = gsl_integration_qags(&F_re, 0, M_PI_2, EPSABS, EPSREL, ITERATIONS, w_re, &result_re, &error_re);
    if(status) 
    {
        std::cout << gsl_strerror (status) << " at p = " << myargs.p << " | estimated error (re): " << error_re << std::endl;
    }
    gsl_integration_workspace_free(w_re);

    gsl_integration_workspace *w_im = gsl_integration_workspace_alloc(ITERATIONS);
    gsl_integration_qags(&F_im, 0, M_PI_2, EPSABS, EPSREL, ITERATIONS, w_im, &result_im, &error_im);
    if(status) {
        std::cout << gsl_strerror (status) << " at p = " << myargs.p << " | estimated error (imag): " << error_im << std::endl;
    }
    gsl_integration_workspace_free(w_im);

    return 2 * M_PI * M_PI * (result_re + 1i * result_im);
}

void writeSamplesToFile(std::string path, std::vector<double> x, std::vector<double> y)
{
    std::ofstream output(path);

    if (!output.is_open())
    {
        std::cerr << "Error opening the file!" << std::endl;
        return;
    }

    for(int i = 0; i < x.size(); i++)
    {
        output << x[i] << ";" << y[i] << std::endl;
    }

    output.close();
}

void writeFuncToFile(std::string path, std::function<double(double)> func, double a, double b, int Nsamples)
{
    std::ofstream output(path);

    if (!output.is_open())
    {
        std::cerr << "Error opening the file!" << std::endl;
        return;
    }

    double dx = (b-a)/(Nsamples-1);
    for(int i = 0; i < Nsamples; i++)
    {
        output << i*dx << ";" << func(i*dx) << std::endl;
    }

    output.close();
}

void linearExtrapolate(std::vector<double> &x, std::vector<double> &y, double x_extr)
{
    std::vector<double> newx, newy;
    newx.push_back(x_extr);
    if(x_extr < x[0])
    {
        double slope = (y[1]-y[0])/(x[1]-x[0]);
        newy.push_back(y[0] + slope * (x_extr - x[0]));
        x.insert(x.begin(), newx.begin(), newx.end());
        y.insert(y.begin(), newy.begin(), newy.end());
    }else if (x_extr > x[0])
    {
        int n = x.size();
        double slope = (y[n-1]-y[n-2])/(x[n-1]-x[n-2]);
        newy.push_back(y[n-1] + slope * (x_extr - x[n-1]));
        x.insert(x.end(), newx.begin(), newx.end());
        y.insert(y.end(), newy.begin(), newy.end());        
    }else
    {
        std::cout << "ERROR! Extrapolation point must lie outside data range" << std::endl;
    }
    
}

int main()
{
    // READ IN FREEZOUT DATA
    csvdata mydata = readcsv("./../../Mathematica/data/ExampleFreezeOut.csv", ",", true, true);
    std::ofstream debug("debug.txt");
    for (int i = 0; i < mydata.header.size(); i++)
    {
        debug << mydata.header[i] << (i == mydata.header.size() - 1 ? "" : ", ") << std::endl;
        // std::cout << mydata.header[i] << (i == mydata.header.size() - 1 ? "" : ", ") << std::endl;
        for (int j = 0; j < mydata.data[i].size(); j++)
            debug << mydata.data[i][j] << std::endl;
        debug << std::endl;
    }
    debug.close();

    // DEFINE INTERPOLATING FUNCTIONS
    std::vector<double> x = mydata.data[0];
    std::vector<double> x1(x), x2(x), x3(x), x4(x), x5(x), x6(x);
    std::vector<double> y1(mydata.data[1]),
                        y2(mydata.data[2]),
                        y3(mydata.data[3]),
                        y4(mydata.data[4]),
                        y5(mydata.data[5]),
                        y6(mydata.data[6]);

    linearExtrapolate(x1,y1,0-0.01);
    linearExtrapolate(x1,y1,M_PI/2.0+0.01);
    linearExtrapolate(x2,y2,0-0.01);
    linearExtrapolate(x2,y2,M_PI/2.0+0.01);
    linearExtrapolate(x3,y3,0-0.01);
    linearExtrapolate(x3,y3,M_PI/2.0+0.01);
    linearExtrapolate(x4,y4,0-0.01);
    linearExtrapolate(x4,y4,M_PI/2.0+0.01);
    linearExtrapolate(x5,y5,0-0.01);
    linearExtrapolate(x5,y5,M_PI/2.0+0.01);
    linearExtrapolate(x6,y6,0-0.01);
    linearExtrapolate(x6,y6,M_PI/2.0+0.01);

    using boost::math::interpolators::pchip;
    auto tau_spline = pchip(
        std::move(x1),
        std::move(y1));
    auto r_spline = pchip(
        std::move(x2),
        std::move(y2));
    auto Dtau_spline = pchip(
        std::move(x3),
        std::move(y3));
    auto Dr_spline = pchip(
        std::move(x4),
        std::move(y4));
    auto ur_spline = pchip(
        std::move(x5),
        std::move(y5));
    auto utau_spline = pchip(
        std::move(x6),
        std::move(y6));

    tau = tau_spline;
    r = r_spline;
    Dtau = Dtau_spline;
    Dr = Dr_spline;
    ur = ur_spline;
    utau = utau_spline;

    // SAVE FOR INSPECTION
    writeSamplesToFile("data/tau_samp.txt",mydata.data[0],mydata.data[1]);
    writeSamplesToFile("data/r_samp.txt",mydata.data[0],mydata.data[2]);
    writeSamplesToFile("data/Dtau_samp.txt",mydata.data[0],mydata.data[3]);
    writeSamplesToFile("data/Dr_samp.txt",mydata.data[0],mydata.data[4]);
    writeSamplesToFile("data/ur_samp.txt",mydata.data[0],mydata.data[5]);
    writeSamplesToFile("data/utau_samp.txt",mydata.data[0],mydata.data[6]);
    
    int NSAMPLE = 1000;
    writeFuncToFile("data/tau_interp.txt",tau,0,M_PI/2.0,NSAMPLE);
    writeFuncToFile("data/r_interp.txt",r,0,M_PI/2.0,NSAMPLE);
    writeFuncToFile("data/Dtau_interp.txt",Dtau,0,M_PI/2.0,NSAMPLE);
    writeFuncToFile("data/Dr_interp.txt",Dr,0,M_PI/2.0,NSAMPLE);
    writeFuncToFile("data/ur_interp.txt",ur,0,M_PI/2.0,NSAMPLE);
    writeFuncToFile("data/utau_interp.txt",utau,0,M_PI/2.0,NSAMPLE);

    // DEFINE PION FIELD AND DERIVATIVE ON FREEZOUT SURFACE
    std::function<double(double)> func = [](double alpha)
    { return sin(50*alpha); };
    std::function<double(double)> Dfunc = [](double alpha)
    { return 0; };

    // COMPUTE SPECTRUM
    std::function<double(double)> spectrfun = [&func, &Dfunc](double p)
    { return std::norm(spectr(p, func, Dfunc)); }; // this is |...|^2
    writeFuncToFile("data/spectr.txt",spectrfun, 0, 2, 100);

}