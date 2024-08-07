#include <iostream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <filesystem>
#include <ctime>   // std::time, std::local_time
#include <iomanip> // std::put_time
#include <boost/program_options.hpp>

#include "savedata.h"
#include "constants.h"
#include "processdata.h"

double M_particle(0.14);

/* #region GSL ODE STUFF */

/*
    For real fields $\phi(x)$ with mass $m_\phi$ given by M_particle, solve the differential equation

        $\epsilon = 1/2 \times ( m_\phi^2 \phi^2 + \chi^2 )$

    where $\chi^2 = -(\partial_\mu\phi)(\partial^\mu\phi)$ and $\partial_\mu\phi\sim u_\mu$ with the 
    4-velociy $u_\mu$. $\epsilon$, $\phi$ and $\chi$ are functions on the freezeout surface and depend
    on 1 parameter $\alpha$ Assume $\epsilon=\text{const.}$, then

        $0 = d\epsilon = m_\phi^2 \phi d\phi + \chi d\chi$

    and further use $d\phi = (\partial_\mu\phi) d^\mu s = u_\mu \chi d^\mu s$ where $d^\mu s$ is the 
    displacement along the freezeout surface, to find

        $0 =m_\phi^2 \phi \chi u_\mu d^\mu s + \chi d\chi$

    Therefore we need to solve 

        $d\phi = \chi u_\mu d^\mu s$
        $d\chi = -m_\phi^2 \phi u_\mu d^\mu s$

    In the following, $\phi=y[0]$ and $\chi=y[1]$.
*/
int func_realfield(double alpha, const double y[], double f[], void *params)
{
    double deriv = (-Dtau(alpha) * utau(alpha) + Dr(alpha) * ur(alpha)); // = u_\mu d^\mu s
    f[0] = y[1] * deriv;
    f[1] = -M_particle * M_particle * y[0] * deriv;
    return GSL_SUCCESS;
};

int jac_realfield(double alpha, const double y[], double *dfdy, double dfdt[], void *params)
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

int func_complexfield(double alpha, const double y[], double f[], void *params)
{
    // ???
    return GSL_SUCCESS;
};

int jac_realfield(double alpha, const double y[], double *dfdy, double dfdt[], void *params)
{
    // ???

    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 2, 2);
    gsl_matrix *m = &dfdy_mat.matrix;
    gsl_matrix_set(m, 0, 0, 0.0);   // not yet implemented
    gsl_matrix_set(m, 0, 1, 0.0);   // not yet implemented
    gsl_matrix_set(m, 1, 0, 0.0);   // not yet implemented
    gsl_matrix_set(m, 1, 1, 0.0);   // not yet implemented
    dfdt[0] = 0.0;                  // not yet implemented
    dfdt[1] = 0.0;                  // not yet implemented
    return GSL_SUCCESS;
};

struct InitialConditions_RealField
{
    std::function<double(double)> phi;
    std::function<double(double)> chi; // =Sqrt[-(\partial_\mu\phi)(\partial^\mu\phi)]
};
InitialConditions_RealField getInitialCondition_RealField(double phi_start, double chi_start)
{
    double EPSREL(1e-3), EPSABS(0);
    double alphastart(0.0), alphaend(M_PI / 2.0);
    double stepsize(1e-6);                // initial value
    double y[2] = {phi_start, chi_start}; // initial value

    const gsl_odeiv2_step_type *TYPE = gsl_odeiv2_step_rk8pd;
    gsl_odeiv2_step *STEP = gsl_odeiv2_step_alloc(TYPE, 2); // dimension = 2
    gsl_odeiv2_control *CONTROL = gsl_odeiv2_control_y_new(EPSABS, EPSREL);
    gsl_odeiv2_evolve *EVOLVE = gsl_odeiv2_evolve_alloc(2); // dimension = 2
    gsl_odeiv2_system SYS({func_realfield, NULL, 2, NULL});       // this is {func, jac, dimension, parameters}

    double alpha(alphastart);
    std::vector<double> alphas({alpha}), phis({y[0]}), chis({y[1]});
    while (alpha < alphaend)
    {
        int status = gsl_odeiv2_evolve_apply(EVOLVE, CONTROL, STEP, &SYS, &alpha, alphaend, &stepsize, y);
        if (status != GSL_SUCCESS)
            break;

        alphas.push_back(alpha);
        phis.push_back(y[0]);
        chis.push_back(y[1]);
    }

    gsl_odeiv2_control_free(CONTROL);
    gsl_odeiv2_evolve_free(EVOLVE);
    gsl_odeiv2_step_free(STEP);

    std::vector<double> alphas1(alphas), alphas2(alphas);
    auto phi_spline = boost::math::interpolators::pchip(
        std::move(alphas1),
        std::move(phis));
    auto chi_spline = boost::math::interpolators::pchip(
        std::move(alphas2),
        std::move(chis));

    return InitialConditions_RealField({phi_spline, chi_spline});
};

struct InitialConditions_ComplexField
{
    std::function<std::complex<double>(double)> phi;
    std::function<double(double)> theta;
};
InitialConditions_PiPlus getInitialCondition_ComplexField(double thetastart)
{
    double EPSREL(1e-3), EPSABS(0);
    double alphastart(0.0), alphaend(M_PI / 2.0);
    double stepsize(1e-6);             // initial value
    double y[1] = {thetastart}; // initial value

    const gsl_odeiv2_step_type *TYPE = gsl_odeiv2_step_rk8pd;
    gsl_odeiv2_step *STEP = gsl_odeiv2_step_alloc(TYPE, 1); // dimension = 1
    gsl_odeiv2_control *CONTROL = gsl_odeiv2_control_y_new(EPSABS, EPSREL);
    gsl_odeiv2_evolve *EVOLVE = gsl_odeiv2_evolve_alloc(1); // dimension = 1
    gsl_odeiv2_system SYS({func_complexfield, NULL, 1, NULL});    // this is {func, jac, dimension, parameters}

    double alpha = alphastart;
    std::vector<double> alphas({alpha}), thetas({thetastart});
    while (alpha < alphaend)
    {
        int status = gsl_odeiv2_evolve_apply(EVOLVE, CONTROL, STEP, &SYS, &alpha, alphaend, &stepsize, y);
        if (status != GSL_SUCCESS)
            break;

        thetas.push_back(y[0]);
        alphas.push_back(alpha);
    }

    gsl_odeiv2_control_free(CONTROL);
    gsl_odeiv2_evolve_free(EVOLVE);
    gsl_odeiv2_step_free(STEP);

    auto theta_spline = boost::math::interpolators::pchip(
        std::move(alphas),
        std::move(thetas));
    std::function<std::complex<double>(double)> phi_spline = [](double alpha)
    {
        std::cout << "COMPLEXFIELD CANNOT BE SET BY INITIALIZATION FUNCTION, YOU NEED TO CONSTRUCT IT YOURSELF FROM THETA" << std::endl;
        return 0.0 + 1i * 0.0;
    };

    return InitialConditions_ComplexField({phi_spline, theta_spline});
};
/* #endregion */

int main(int ac, char* av[])
{
    csvdata freezeoutdata = ProcessFreezeoutData();

    // CREATE DIRECTORY TO SAVE FILES
    std::time_t t = std::time(nullptr);
    std::tm tm = *std::localtime(&t);
    std::stringstream timestamp_sstr;
    timestamp_sstr << std::put_time(&tm, "%Y%m%d_%H%M%S");
    std::string timestamp = timestamp_sstr.str();
    std::string pathname = "data/init_"+timestamp;
    std::filesystem::create_directories(pathname);

    /* #region COMMAND LINE OPTIONS */
    // DEFINE VARIABLES TO BE SET
    double epsilon, ratio;
    int sign;

    // DECLARE SUPPORTED OPTIONS
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("M", po::value<double>()->default_value(0.14), "particle mass")
        ("eps", po::value<double>()->default_value(0.001 * 0.160054), "initial (constant) energy density")
        ("ratio", po::value<double>()->default_value(1.0), "Epot/eps at alpha=0")
        ("sign",po::value<int>()->default_value(1),"sign of phi/dphi at alpha=0");

    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);

    // PROCESS CMDLINE OPTIONS
    if (vm.count("help"))    {
        std::cout << desc << "\n";
        return 1;
    }
    M_particle = vm["M"].as<double>();
    epsilon = vm["eps"].as<double>();
    ratio = vm["ratio"].as<double>();
    sign = vm["sign"].as<int>();
    /* #endregion */

    // SPECIFIY ODE INITIAL CONDITIONS
    double Epot(ratio * epsilon), Ekin((1 - ratio) * epsilon);

    // INTEGRATING THE ODE SYSTEMS
    InitialConditions_RealField ic_pi0 = getInitialCondition_RealField(sqrt(2 * Epot) / M_particle, sqrt(2 * Ekin));
    InitialConditions_ComplexField ic_piplus = getInitialCondition_ComplexField(0.6);
    double n(1.0);
    ic_piplus.piplus = [&](double alpha)
    { return sqrt(n) * (cos(ic_piplus.theta(alpha)) + 1i * sin(ic_piplus.theta(alpha))); };

    // SAVE TO FILE
    std::stringstream mass_ss;
    mass_ss << "m: " << M_particle;
    std::vector<std::string> comments({
        timestamp, mass_ss.str()
    });

    writeFuncToFile(pathname+"/pi0_initial.txt", ic_pi0.pi0, 0, M_PI / 2.0, 1000, {"alpha", "pi0Re", "pi0Im"},comments);
    writeFuncToFile(pathname+"/chi_initial.txt", ic_pi0.chi, 0, M_PI / 2.0, 1000, {"alpha", "chiRe", "chiIm"},comments);

    std::function<std::complex<double>(double)> Dpi0 = [&](double alpha)
    { return ic_pi0.chi(alpha) * (-Dr(alpha) * utau(alpha) + Dtau(alpha) * ur(alpha)); };
    writeFuncToFile(pathname+"/Dpi0_initial.txt", Dpi0, 0, M_PI / 2.0, 1000, {"alpha", "Dpi0Re", "Dpi0Im"},comments);
    writeFunctionsToFile<std::complex<double>>(
        pathname+"/initialfields_pi0.csv",
        {ic_pi0.pi0, Dpi0},
        0, M_PI / 2.0, 1000,
        {"alpha", "pi0Re", "pi0Im", "Dpi0Re", "Dpi0Im"},comments);

    writeFuncToFile<std::complex<double>>(pathname+"/piplus_initial.txt", ic_piplus.piplus, 0, M_PI / 2.0, 1000, {"alpha", "piplusRe", "piplusIm"},comments);
    writeFuncToFile(pathname+"/theta_initial.txt", ic_piplus.theta, 0, M_PI / 2.0, 1000, {"alpha", "thetaRe", "thetaIm"},comments);
    
    std::function<std::complex<double>(double)> Dpiplus = [&](double alpha)
    { return 1i * ic_piplus.piplus(alpha) * m_pion * (-Dr(alpha) * utau(alpha) + Dtau(alpha) * ur(alpha)); }; // THIS ASSUMES CHI = MPION
    writeFuncToFile<std::complex<double>>(pathname+"/Dpiplus_initial.txt", Dpiplus, 0, M_PI / 2.0, 1000, {"alpha", "DpiplusRe", "DpiplusIm"},comments);
    writeFunctionsToFile<std::complex<double>>(
        pathname+"/initialfields_piplus.csv",
        {ic_piplus.piplus, Dpiplus},
        0, M_PI / 2.0, 1000,
        {"alpha", "piplusRe", "piplusIm", "DpiplusRe", "DpiplusIm"},
        comments);

    return 0;
}