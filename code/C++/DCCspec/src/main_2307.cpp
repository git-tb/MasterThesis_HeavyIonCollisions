#include <iostream>
#include <cmath>
#include <complex>
#include <functional>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <map>
#include <boost/math/interpolators/pchip.hpp>
#include <boost/program_options.hpp>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include <algorithm>
#include <iterator>

#include "filereader.h"
#include "spectr.h"

using namespace std::complex_literals;

double m_pion = 0.14;
double GeVtoIfm = 5.0677;
double IfmtoGeV = 1 / GeVtoIfm;
double fmtoIGeV = GeVtoIfm;
double IGeVtofm = 1 / fmtoIGeV;

std::function<double(double)> tau, r, Dtau, Dr, ur, utau;

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
// std::complex<double> spectr(double p, std::function<double(double)> func, std::function<double(double)> Dfunc)
// {
//     args myargs = {p, func, Dfunc};

//     gsl_function F_re, F_im;
//     F_re.function = integrand_re;
//     F_im.function = integrand_im;
//     F_re.params = &myargs;
//     F_im.params = &myargs;

//     double EPSABS(0), EPSREL(1e-2);
//     int ITERATIONS(1000);
//     int KEY(6);

//     gsl_set_error_handler_off();

//     double result_re, result_im, error_re, error_im;
//     gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(ITERATIONS);
//     int status = gsl_integration_qag(&F_re, 0, M_PI_2, EPSABS, EPSREL, ITERATIONS, KEY, workspace, &result_re, &error_re);
//     if (status)
//         std::cout << gsl_strerror(status) << " at p = " << myargs.p << " | estimated error (re): " << error_re << std::endl;

//     status = gsl_integration_qag(&F_im, 0, M_PI_2, EPSABS, EPSREL, ITERATIONS, KEY, workspace, &result_im, &error_im);
//     if (status)
//         std::cout << gsl_strerror(status) << " at p = " << myargs.p << " | estimated error (imag): " << error_im << std::endl;

//     gsl_integration_workspace_free(workspace);

//     return 2 * M_PI * M_PI * (result_re + 1i * result_im);
// }

// std::vector<std::complex<double>> spectr(
//     std::vector<double> ps,
//     std::function<double(double)> func,
//     std::function<double(double)> Dfunc)
// {
//     // FOR EVALUATION ON ARRAY OF PS, WE DON'T NEED TO ALLOCATE AND FREE
//     //  THE MEMORY FOR THE INTEGRATION WORKSPACE OVER AND OVER AGAIN
//     std::vector<std::complex<double>> result(ps.size());

//     double EPSABS(0), EPSREL(1e-2);
//     int ITERATIONS(1000);
//     int KEY(6);
//     gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(ITERATIONS);
//     gsl_set_error_handler_off();

//     for (int i = 0; i < ps.size(); i++)
//     {
//         std::cout << i + 1 << " / " << ps.size() << std::endl;

//         double p = ps[i];

//         args myargs = {p, func, Dfunc};

//         gsl_function F_re, F_im;
//         F_re.function = integrand_re;
//         F_im.function = integrand_im;
//         F_re.params = &myargs;
//         F_im.params = &myargs;

//         double result_re, result_im, error_re, error_im;

//         // int status = gsl_integration_qags(&F_re, 0, M_PI_2, EPSABS, EPSREL, ITERATIONS, workspace, &result_re, &error_re);
//         int status = gsl_integration_qag(&F_re, 0, M_PI_2, EPSABS, EPSREL, ITERATIONS, KEY, workspace, &result_re, &error_re);
//         if (status)
//             std::cout << gsl_strerror(status) << " at p = " << myargs.p << " | estimated error (re): " << error_re << std::endl;

//         // status = gsl_integration_qags(&F_im, 0, M_PI_2, EPSABS, EPSREL, ITERATIONS, workspace, &result_im, &error_im);
//         status = gsl_integration_qag(&F_im, 0, M_PI_2, EPSABS, EPSREL, ITERATIONS, KEY, workspace, &result_im, &error_im);
//         if (status)
//             std::cout << gsl_strerror(status) << " at p = " << myargs.p << " | estimated error (imag): " << error_im << std::endl;

//         result[i] = 2 * M_PI * M_PI * (result_re + 1i * result_im);
//     }

//     gsl_integration_workspace_free(workspace);
//     return result;
// }

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

/* #region HELPER FUNCTIONS TO STORE ARRAYS AND FUNCTION VALUES */
template <typename T>
void writeSamplesToFile(std::string path, std::vector<double> x, std::vector<T> y)
{
    std::ofstream output(path);

    if (!output.is_open())
    {
        std::cerr << "Error opening the file!" << std::endl;
        return;
    }

    for (int i = 0; i < x.size(); i++)
    {
        output << x[i] << ";" << std::real(y[i]) << ";" << std::imag(y[i]) << std::endl;
    }

    output.close();
}

template <typename T>
void writeFuncToFile(std::string path, std::function<T(double)> func, double a, double b, int Nsamples)
{
    std::ofstream output(path);

    if (!output.is_open())
    {
        std::cerr << "Error opening the file!" << std::endl;
        return;
    }

    double dx = (b - a) / (Nsamples - 1);
    for (int i = 0; i < Nsamples; i++)
    {
        T y = func(i * dx);
        output << i * dx << ";" << std::real(y) << ";" << std::imag(y) << std::endl;
    }

    output.close();
}
/* #endregion */

/* #region SIMPLE EXTRAPOLATION TO PREVENT DOMAIN ERRORS IN INTERPOLATION FUNCTIONS */
void linearExtrapolate(std::vector<double> &x, std::vector<double> &y, double x_extr)
{
    std::vector<double> newx, newy;
    newx.push_back(x_extr);
    if (x_extr < x[0])
    {
        double slope = (y[1] - y[0]) / (x[1] - x[0]);
        newy.push_back(y[0] + slope * (x_extr - x[0]));
        x.insert(x.begin(), newx.begin(), newx.end());
        y.insert(y.begin(), newy.begin(), newy.end());
    }
    else if (x_extr > x[0])
    {
        int n = x.size();
        double slope = (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]);
        newy.push_back(y[n - 1] + slope * (x_extr - x[n - 1]));
        x.insert(x.end(), newx.begin(), newx.end());
        y.insert(y.end(), newy.begin(), newy.end());
    }
    else
    {
        std::cout << "ERROR! Extrapolation point must lie outside data range" << std::endl;
    }
}
/* #endregion */

/* #region ODE SYSTEM FOR INITIAL CONDITIONS */
int func_pi0(double alpha, const double y[], double f[], void *params)
{
    double deriv = (-Dtau(alpha) * utau(alpha) + Dr(alpha) * ur(alpha));
    f[0] = y[1] * deriv;
    f[1] = -m_pion * m_pion * y[0] * deriv;
    return GSL_SUCCESS;
};

int jac_pi0(double alpha, const double y[], double *dfdy, double dfdt[], void *params)
{
    double deriv = (-Dtau(alpha) * utau(alpha) + Dtau(alpha) * ur(alpha));

    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 2, 2);
    gsl_matrix *m = &dfdy_mat.matrix;
    gsl_matrix_set(m, 0, 0, 0.0);
    gsl_matrix_set(m, 0, 1, deriv);
    gsl_matrix_set(m, 1, 0, -m_pion * m_pion * deriv);
    gsl_matrix_set(m, 1, 1, 0);
    dfdt[0] = 0.0; // not yet implemented
    dfdt[1] = 0.0; // not yet implemented
    return GSL_SUCCESS;
};

int func_piplus(double alpha, const double y[], double f[], void *params)
{
    double deriv = (-Dtau(alpha) * utau(alpha) + Dr(alpha) * ur(alpha));
    f[0] = m_pion * deriv; // chi = m_pion
    return GSL_SUCCESS;
};

struct InitialConditions_Pi0
{
    std::function<double(double)> pi0;
    std::function<double(double)> chi; // =Sqrt[-(\partial_\mu\pi^0)(\partial^\mu\pi^0)]
};
InitialConditions_Pi0 getInitialCondition_Pi0(double pi0_start, double chi_start)
{
    double EPSREL(1e-3), EPSABS(0);
    double alphastart(0.0), alphaend(M_PI / 2.0);
    double stepsize(1e-6);                // initial value
    double y[2] = {pi0_start, chi_start}; // initial value

    const gsl_odeiv2_step_type *TYPE = gsl_odeiv2_step_rk8pd;
    gsl_odeiv2_step *STEP = gsl_odeiv2_step_alloc(TYPE, 2); // dimension = 2
    gsl_odeiv2_control *CONTROL = gsl_odeiv2_control_y_new(EPSABS, EPSREL);
    gsl_odeiv2_evolve *EVOLVE = gsl_odeiv2_evolve_alloc(2); // dimension = 2
    gsl_odeiv2_system SYS({func_pi0, NULL, 2, NULL});          // this is {func, jac, dimension, parameters}

    double alpha(alphastart);
    std::vector<double> alphas({alpha}), pi0s({y[0]}), chis({y[1]});
    while (alpha < alphaend)
    {
        int status = gsl_odeiv2_evolve_apply(EVOLVE, CONTROL, STEP, &SYS, &alpha, alphaend, &stepsize, y);
        if (status != GSL_SUCCESS)
            break;

        alphas.push_back(alpha);
        pi0s.push_back(y[0]);
        chis.push_back(y[1]);
    }

    gsl_odeiv2_control_free(CONTROL);
    gsl_odeiv2_evolve_free(EVOLVE);
    gsl_odeiv2_step_free(STEP);

    std::vector<double> alphas1(alphas), alphas2(alphas);
    auto pi0_spline = boost::math::interpolators::pchip(
        std::move(alphas1),
        std::move(pi0s));
    auto chi_spline = boost::math::interpolators::pchip(
        std::move(alphas2),
        std::move(chis));

    return InitialConditions_Pi0({pi0_spline,chi_spline});
};

struct InitialConditions_PiPlus
{
    std::function<std::complex<double>(double)> piplus;
    std::function<double(double)> theta; // =Sqrt[-(\partial_\mu\pi^0)(\partial^\mu\pi^0)]
};
InitialConditions_PiPlus getInitialCondition_PiPlus(double thetastart)
{
    double EPSREL(1e-3), EPSABS(0);
    double alphastart(0.0), alphaend(M_PI / 2.0);
    double stepsize(1e-6);                // initial value
    double y_piplus[1] = {thetastart}; // initial value

    const gsl_odeiv2_step_type *TYPE = gsl_odeiv2_step_rk8pd;
    gsl_odeiv2_step *STEP = gsl_odeiv2_step_alloc(TYPE, 1); // dimension = 1
    gsl_odeiv2_control *CONTROL = gsl_odeiv2_control_y_new(EPSABS, EPSREL);
    gsl_odeiv2_evolve *EVOLVE = gsl_odeiv2_evolve_alloc(1); // dimension = 1
    gsl_odeiv2_system SYS({func_piplus, NULL, 1, NULL});          // this is {func, jac, dimension, parameters}

    double alpha = alphastart;
    std::vector<double> alphas({alpha}), thetas({thetastart});
    while (alpha < alphaend)
    {
        int status = gsl_odeiv2_evolve_apply(EVOLVE, CONTROL, STEP, &SYS, &alpha, alphaend, &stepsize, y_piplus);
        if (status != GSL_SUCCESS)
            break;
                
        thetas.push_back(y_piplus[0]);
        alphas.push_back(alpha);
    }

    gsl_odeiv2_control_free(CONTROL);
    gsl_odeiv2_evolve_free(EVOLVE);
    gsl_odeiv2_step_free(STEP);

    auto theta_spline = boost::math::interpolators::pchip(
        std::move(alphas),
        std::move(thetas));
    std::function<std::complex<double>(double)> piplus = [](double alpha)
    {
        std::cout << "PIPLUS CANNOT BE SET BY INITIALIZATION FUNCTION, YOU NEED TO CONSTRUCT IT YOURSELF FROM THETA" << std::endl;
        return 0.0 + 1i * 0.0;
    };

    return InitialConditions_PiPlus({piplus, theta_spline});
};
/* #endregion */

int main(int ac, char* av[])
{
    /* #region COMMAND LINE OPTIONS */
    // Declare the supported options.
    // namespace po = boost::program_options;
    // po::options_description desc("Allowed options");
    // desc.add_options()
    //     ("help", "produce help message")
    //     ("pTmax", po::value<double>()->default_value(1.0), "set pTmax");

    // po::variables_map vm;
    // po::store(po::parse_command_line(ac, av, desc), vm);
    // po::notify(vm);

    // if (vm.count("help"))
    // {
    //     std::cout << desc << "\n";
    //     return 1;
    // }
    // double pTmax(vm["pTmax"].as<double>());
    // std::cout << pTmax << std::endl;
    // return 1;
    /* #endregion */

    /* #region INTERPOLATING FUNCTIONS FROM THE FREEZOUT DATA */
    // READ IN FREEZOUT DATA
    csvdata mydata = readcsv("./../../Mathematica/data/ExampleFreezeOut.csv", ",", true, true);
    csvdata initialdata = readcsv("./data/initialfield.csv", ",", true, true);
    std::ofstream debug("debug.txt");
    for (int i = 0; i < mydata.header.size(); i++)
    {
        debug << mydata.header[i] << (i == mydata.header.size() - 1 ? "" : ", ") << std::endl;
        for (int j = 0; j < mydata.data[i].size(); j++)
            debug << mydata.data[i][j] << std::endl;
        debug << std::endl;
    }
    for (int i = 0; i < initialdata.header.size(); i++)
    {
        debug << initialdata.header[i] << (i == initialdata.header.size() - 1 ? "" : ", ") << std::endl;
        for (int j = 0; j < initialdata.data[i].size(); j++)
            debug << initialdata.data[i][j] << std::endl;
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

    linearExtrapolate(x1, y1, 0 - 0.01);
    linearExtrapolate(x1, y1, M_PI / 2.0 + 0.01);
    linearExtrapolate(x2, y2, 0 - 0.01);
    linearExtrapolate(x2, y2, M_PI / 2.0 + 0.01);
    linearExtrapolate(x3, y3, 0 - 0.01);
    linearExtrapolate(x3, y3, M_PI / 2.0 + 0.01);
    linearExtrapolate(x4, y4, 0 - 0.01);
    linearExtrapolate(x4, y4, M_PI / 2.0 + 0.01);
    linearExtrapolate(x5, y5, 0 - 0.01);
    linearExtrapolate(x5, y5, M_PI / 2.0 + 0.01);
    linearExtrapolate(x6, y6, 0 - 0.01);
    linearExtrapolate(x6, y6, M_PI / 2.0 + 0.01);

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
    writeSamplesToFile("data/tau_samp.txt", mydata.data[0], mydata.data[1]);
    writeSamplesToFile("data/r_samp.txt", mydata.data[0], mydata.data[2]);
    writeSamplesToFile("data/Dtau_samp.txt", mydata.data[0], mydata.data[3]);
    writeSamplesToFile("data/Dr_samp.txt", mydata.data[0], mydata.data[4]);
    writeSamplesToFile("data/ur_samp.txt", mydata.data[0], mydata.data[5]);
    writeSamplesToFile("data/utau_samp.txt", mydata.data[0], mydata.data[6]);

    int NSAMPLE = 1000;
    writeFuncToFile("data/tau_interp.txt", tau, 0, M_PI / 2.0, NSAMPLE);
    writeFuncToFile("data/r_interp.txt", r, 0, M_PI / 2.0, NSAMPLE);
    writeFuncToFile("data/Dtau_interp.txt", Dtau, 0, M_PI / 2.0, NSAMPLE);
    writeFuncToFile("data/Dr_interp.txt", Dr, 0, M_PI / 2.0, NSAMPLE);
    writeFuncToFile("data/ur_interp.txt", ur, 0, M_PI / 2.0, NSAMPLE);
    writeFuncToFile("data/utau_interp.txt", utau, 0, M_PI / 2.0, NSAMPLE);
    /* #endregion */

    // INITIAL CONDITIONS
    double epsilon = 0.001 * 0.160054; // energy density of condensate
                                       // epsilon = Epot + Ekin
                                       //  Epot = (1/2) m^2 phi^2
                                       //  Ekin = (1/2) chi^2
    double ratio(0.5);                 // Epot/epsilon
    int sign(-1);                      // +-1, relative sign between pi0(alpha=0) and \partial pi0(alpha=0)
    double Epot(ratio * epsilon), Ekin((1 - ratio) * epsilon);

    std::cout << "generatin pi0 initial conditions..." << std::flush;
    InitialConditions_Pi0 ic_pi0 = getInitialCondition_Pi0(sqrt(2*Epot)/m_pion, sqrt(2*Ekin));
    std::cout << " done" << std::endl;
    std::cout << "generatin pi+ initial conditions..." << std::flush;
    InitialConditions_PiPlus ic_piplus = getInitialCondition_PiPlus(0.6);
    std::cout << " done" << std::endl;
    double n(1.0);
    ic_piplus.piplus = [&](double alpha) { return sqrt(n) * (cos(ic_piplus.theta(alpha)) + 1i * sin(ic_piplus.theta(alpha))); };

    std::cout << "saving initial conditions..." << std::flush;
    writeFuncToFile("data/pi0_initial.txt", ic_pi0.pi0, 0, M_PI / 2.0, 1000);
    writeFuncToFile("data/chi_initial.txt", ic_pi0.chi, 0, M_PI / 2.0, 1000);
    writeFuncToFile<std::complex<double>>("data/piplus_initial.txt", ic_piplus.piplus, 0, M_PI / 2.0, 1000);
    writeFuncToFile("data/theta_initial.txt", ic_piplus.theta, 0, M_PI / 2.0, 1000);
    std::cout << " done" << std::endl;
    
    
    // DEFINE PION FIELD AND DERIVATIVE ON FREEZOUT SURFACE
    std::function<std::complex<double>(double)> func = [&](double alpha)
    { return ic_piplus.piplus(alpha); };
    std::function<std::complex<double>(double)> Dfunc = [&](double alpha)
    { return 1i * ic_piplus.piplus(alpha) * m_pion * (-Dr(alpha) * utau(alpha) + Dtau(alpha) * ur(alpha)); };
    // std::function<std::complex<double>(double)> func = [&](double alpha)
    // { return pi0_spline(alpha); };
    // std::function<std::complex<double>(double)> Dfunc = [&](double alpha)
    // { return chi_spline(alpha) * (-Dr(alpha) * utau(alpha) + Dtau(alpha) * ur(alpha)); };

    writeFuncToFile("data/field0.txt", func, 0, M_PI / 2.0, 1000);
    writeFuncToFile("data/field0_deriv.txt", Dfunc, 0, M_PI / 2.0, 1000);

    // COMPUTE SPECTRUM
    // USING SLIGHTLY OPTIMIZED VERSION FOR ARRAY-LIKE ARGUMENTS
    double pmin(0), pmax(1.0);
    int Nps(100);

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

    writeSamplesToFile("data/spectr.txt", ps, myspectr_abs2);
    writeSamplesToFile("data/spectr_anti.txt", ps, myspectr_anti_abs2);
}