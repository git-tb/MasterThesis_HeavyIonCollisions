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

using namespace std::complex_literals;


/* #region HELPER FUNCTIONS FOR SPECTRUM COMPUTATION */
double w(double p) { return sqrt(p * p + m_pion * m_pion); }
double J0rp(double alpha, double p) { return std::cyl_bessel_j(0, r(alpha) * p * GeVtoIfm); }
double J1rp(double alpha, double p) { return std::cyl_bessel_j(1, r(alpha) * p * GeVtoIfm); }
double Y0tw(double alpha, double p) { return std::cyl_neumann(0, tau(alpha) * w(p) * GeVtoIfm); }
double Y1tw(double alpha, double p) { return std::cyl_neumann(1, tau(alpha) * w(p) * GeVtoIfm); }
double J0tw(double alpha, double p) { return std::cyl_bessel_j(0, tau(alpha) * w(p) * GeVtoIfm); }
double J1tw(double alpha, double p) { return std::cyl_bessel_j(1, tau(alpha) * w(p) * GeVtoIfm); }

std::complex<double> G0(double alpha, double p) { return J0rp(alpha, p) * (-Y0tw(alpha, p) + 1i * J0tw(alpha, p)); }
std::complex<double> G0_ANTI(double alpha, double p) { return J0rp(alpha, p) * (-Y0tw(alpha, p) - 1i * J0tw(alpha, p)); }
std::complex<double> G1(double alpha, double p) {
    return J1rp(alpha, p) * (-Y0tw(alpha, p) + 1i * J0tw(alpha, p)) * Dtau(alpha) * fmtoIGeV * p +
           J0rp(alpha, p) * (-Y1tw(alpha, p) + 1i * J1tw(alpha, p)) * Dr(alpha) * fmtoIGeV * w(p);
}
std::complex<double> G1_ANTI(double alpha, double p)
{
    return J1rp(alpha, p) * (-Y0tw(alpha, p) - 1i * J0tw(alpha, p)) * Dtau(alpha) * fmtoIGeV * p +
           J0rp(alpha, p) * (-Y1tw(alpha, p) - 1i * J1tw(alpha, p)) * Dr(alpha) * fmtoIGeV * w(p);
}

// struct args
// {
//     double p;
//     std::function<double(double)> func;
//     std::function<double(double)> Dfunc;
// };

// auto integrand_re = [](double alpha, void *params)
// {
//     args myargs = *(struct args *)params;
//     return tau(alpha) * r(alpha) * fmtoIGeV * fmtoIGeV * (myargs.Dfunc(alpha) * H1(alpha, myargs.p).real() + myargs.func(alpha) * H2(alpha, myargs.p).real());
// };
// auto integrand_im = [](double alpha, void *params)
// {
//     args myargs = *(struct args *)params;
//     return tau(alpha) * r(alpha) * fmtoIGeV * fmtoIGeV * (myargs.Dfunc(alpha) * H1(alpha, myargs.p).imag() + myargs.func(alpha) * H2(alpha, myargs.p).imag());
// };

struct args
{
    double p;
    std::function<std::complex<double>(double)> func, Dfunc;
};

std::complex<double> (*integrand)(double, void*) = [](double alpha, void* params)
{
    args myargs = *(struct args *)params;
    return tau(alpha) * r(alpha) * fmtoIGeV * fmtoIGeV * (
        myargs.Dfunc(alpha) * G0(alpha, myargs.p) + myargs.func(alpha) * G1(alpha, myargs.p)
    );
};

std::complex<double> (*integrand_anti)(double, void*) = [](double alpha, void* params)
{
    args myargs = *(struct args *)params;
    return tau(alpha) * r(alpha) * fmtoIGeV * fmtoIGeV * (
        myargs.Dfunc(alpha) * G0_ANTI(alpha, myargs.p) + myargs.func(alpha) * G1_ANTI(alpha, myargs.p)
    );
};
/* #endregion */

/* #region SPECTRUM COMPUTATION AT SINGLE P-VALUE OR LIST OF P-VALUES */
std::vector<std::complex<double>> spectr(
    std::vector<double> ps,
    std::function<std::complex<double>(double)> func,
    std::function<std::complex<double>(double)> Dfunc,
    bool anti = false)
{
    // FOR EVALUATION ON ARRAY OF PS, WE DON'T NEED TO ALLOCATE AND FREE
    //  THE MEMORY FOR THE INTEGRATION WORKSPACE OVER AND OVER AGAIN
    std::vector<std::complex<double>> result(ps.size());

    double EPSABS(0), EPSREL(1e-2);
    int ITERATIONS(1000);
    int KEY(6);
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(ITERATIONS);
    gsl_set_error_handler_off();

    for (int i = 0; i < ps.size(); i++)
    {
        std::cout << i + 1 << " / " << ps.size() << std::endl;

        double p = ps[i];

        args myargs = {p, func, Dfunc};

        gsl_function F_re, F_im;
        if(!anti)
        {
            F_re.function = [](double alpha, void* params){ return std::real(integrand(alpha, params)); };
            F_im.function = [](double alpha, void* params){ return std::imag(integrand(alpha, params)); };
        }
        else
        {
            F_re.function = [](double alpha, void* params){ return std::real(integrand_anti(alpha, params)); };
            F_im.function = [](double alpha, void* params){ return std::imag(integrand_anti(alpha, params)); };
        }
        F_re.params = &myargs;
        F_im.params = &myargs;

        double result_re, result_im, error_re, error_im;

        // int status = gsl_integration_qags(&F_re, 0, M_PI_2, EPSABS, EPSREL, ITERATIONS, workspace, &result_re, &error_re);
        int status = gsl_integration_qag(&F_re, 0, M_PI_2, EPSABS, EPSREL, ITERATIONS, KEY, workspace, &result_re, &error_re);
        if (status)
            std::cout << gsl_strerror(status) << " at p = " << myargs.p << " | estimated error (re): " << error_re << std::endl;

        // status = gsl_integration_qags(&F_im, 0, M_PI_2, EPSABS, EPSREL, ITERATIONS, workspace, &result_im, &error_im);
        status = gsl_integration_qag(&F_im, 0, M_PI_2, EPSABS, EPSREL, ITERATIONS, KEY, workspace, &result_im, &error_im);
        if (status)
            std::cout << gsl_strerror(status) << " at p = " << myargs.p << " | estimated error (imag): " << error_im << std::endl;

        result[i] = 2 * M_PI * M_PI * (result_re + 1i * result_im);
    }

    gsl_integration_workspace_free(workspace);
    return result;
}
/* #endregion */

int main(int ac, char* av[])
{
    /* #region COMMAND LINE OPTIONS */
    // DEFINE VARIABLES TO BE SET
    double pmin(0), pmax(1.0);
    int Nps(100);
    std::string initdata;

    // DECLARE SUPPORTED OPTIONS
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("pTmax", po::value<double>()->default_value(1.0), "spectrum is computed on [0,pTmax]")
        ("NpT", po::value<int>()->default_value(100), "number of sample points within [0,pTmax]")
        ("initpath",po::value<std::string>(),"csv file containing initial field data");

    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);

    // PROCESS CMDLINE OPTIONS
    if (vm.count("help"))    {
        std::cout << desc << "\n";
        return 1;
    }
    pmax = vm["pTmax"].as<double>();
    Nps = vm["NpT"].as<int>();
    initdata = vm["initpath"].as<std::string>();
    /* #endregion */

    // CREATE DIRECTORY TO SAVE FILES
    std::time_t t = std::time(nullptr);
    std::tm tm = *std::localtime(&t);
    std::stringstream timestamp_sstr;
    timestamp_sstr << std::put_time(&tm, "%Y%m%d_%H%M%S");
    std::string timestamp = timestamp_sstr.str();
    std::string pathname = "data/spec_"+timestamp;
    std::filesystem::create_directories(pathname);

    csvdata mydata = ProcessFreezeoutData();
    csvdata initialdata = ProcessInitialData(initdata);

    // SAVE FOR INSPECTION
    //  THE RAW SAMPLES
    writeSamplesToFile(pathname+"/tau_samp.txt", mydata.data[0], mydata.data[1],{"alpha","tauRe","tauIm"},{timestamp});
    writeSamplesToFile(pathname+"/r_samp.txt", mydata.data[0], mydata.data[2],{"alpha","rRe","rIm"},{timestamp});
    writeSamplesToFile(pathname+"/Dtau_samp.txt", mydata.data[0], mydata.data[3],{"alpha","DtauRe","DtauIm"},{timestamp});
    writeSamplesToFile(pathname+"/Dr_samp.txt", mydata.data[0], mydata.data[4],{"alpha","DrRe","DrIm"},{timestamp});
    writeSamplesToFile(pathname+"/ur_samp.txt", mydata.data[0], mydata.data[5],{"alpha","urRe","urIm"},{timestamp});
    writeSamplesToFile(pathname+"/utau_samp.txt", mydata.data[0], mydata.data[6],{"alpha","utauRe","utauIm"},{timestamp});
    {
        std::vector<std::complex<double>> f0dat, Df0dat;    // REAL AND IMAGINARY PART ARE READ IN AS 2 SEPARATE REAL-VALUED ARRAYS
                                                            //  BUT SHOULD BE PRINTED AS 1 COMBINED COMPLEX-VALUED ARRAY
        for(int i  = 0; i < initialdata.data[0].size(); i++)
        {
            f0dat.push_back(initialdata.data[1][i] + 1i * initialdata.data[2][i]);
            Df0dat.push_back(initialdata.data[3][i] + 1i * initialdata.data[4][i]);
        }
        writeSamplesToFile(pathname+"/f0_samp.txt", initialdata.data[0], f0dat,{"alpha","f0Re","f0Im"},{timestamp});
        writeSamplesToFile(pathname+"/Df0_samp.txt", initialdata.data[0], Df0dat,{"alpha","Df0Re","Df0Im"},{timestamp});
    }

    // THE INTERPOLATED FUNCTIONS, SAMPLED AT SOME PRESCRIBED RESOLUTION
    int NSAMPLE = 1000;
    writeFuncToFile(pathname+"/tau_interp.txt", tau, 0, M_PI / 2.0, NSAMPLE,{"alpha","tauRe","tauIm"},{timestamp});
    writeFuncToFile(pathname+"/r_interp.txt", r, 0, M_PI / 2.0, NSAMPLE,{"alpha","rRe","rIm"},{timestamp});
    writeFuncToFile(pathname+"/Dtau_interp.txt", Dtau, 0, M_PI / 2.0, NSAMPLE,{"alpha","DtauRe","DtauIm"},{timestamp});
    writeFuncToFile(pathname+"/Dr_interp.txt", Dr, 0, M_PI / 2.0, NSAMPLE,{"alpha","DrRe","DrIm"},{timestamp});
    writeFuncToFile(pathname+"/ur_interp.txt", ur, 0, M_PI / 2.0, NSAMPLE,{"alpha","urRe","urIm"},{timestamp});
    writeFuncToFile(pathname+"/utau_interp.txt", utau, 0, M_PI / 2.0, NSAMPLE,{"alpha","utauRe","utauIm"},{timestamp});

    writeFuncToFile(pathname+"/f0_interp.txt", f0, 0, M_PI / 2.0, NSAMPLE,{"alpha","f0Re","f0Im"},{timestamp});
    writeFuncToFile(pathname+"/Df0_interp.txt", Df0, 0, M_PI / 2.0, NSAMPLE,{"alpha","Df0Re","Df0Im"},{timestamp});    
    
    // DEFINE PION FIELD AND DERIVATIVE ON FREEZOUT SURFACE
    std::function<std::complex<double>(double)> func = [](double alpha)
    { return f0(alpha); };
    std::function<std::complex<double>(double)> Dfunc = [](double alpha)
    { return Df0(alpha); };

    writeFuncToFile(pathname+"/field0.txt", func, 0, M_PI / 2.0, 1000,{"alpha","field0Re","field0Im"},{timestamp});
    writeFuncToFile(pathname+"/field0_deriv.txt", Dfunc, 0, M_PI / 2.0, 1000,{"alpha","Dfield0Re","Dfield0Im"},{timestamp});

    // COMPUTE SPECTRUM
    // USING SLIGHTLY OPTIMIZED VERSION FOR ARRAY-LIKE ARGUMENTS
    std::vector<double> ps(Nps);
    for (int i = 0; i < ps.size(); i++)
        ps[i] = pmin + i * (pmax - pmin) / (Nps - 1);

    std::vector<std::complex<double>> myspectr = spectr(ps, func, Dfunc,false);
    std::vector<std::complex<double>> myspectr_anti = spectr(ps, func, Dfunc,true);

    std::vector<double> myspectr_abs2(myspectr.size());
    for (int i = 0; i < myspectr.size(); i++)
        myspectr_abs2[i] = (1 / std::pow(2 * M_PI, 3)) * std::norm(myspectr[i]);

    std::vector<double> myspectr_anti_abs2(myspectr_anti.size());
    for (int i = 0; i < myspectr_anti.size(); i++)
        myspectr_anti_abs2[i] = (1 / std::pow(2 * M_PI, 3)) * std::norm(myspectr_anti[i]);

    writeSamplesToFile(pathname+"/spectr.txt", ps, myspectr_abs2,{"pT","abs2Re","abs2Im"},{timestamp});
    writeSamplesToFile(pathname+"/spectr_anti.txt", ps, myspectr_anti_abs2,{"pT","abs2Re","abs2Im"},{timestamp});
}