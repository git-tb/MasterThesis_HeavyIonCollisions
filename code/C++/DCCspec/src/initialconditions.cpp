#include <iostream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include "savedata.h"
#include "constants.h"
#include "processdata.h"

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

int main()
{
    csvdata freezeoutdata = ProcessFreezeoutData();

    double epsilon = 0.001 * 0.160054; // energy density of condensate
                                       // epsilon = Epot + Ekin
                                       //  Epot = (1/2) m^2 phi^2
                                       //  Ekin = (1/2) chi^2
    double ratio(0.5);                 // Epot/epsilon
    int sign(-1);                      // +-1, relative sign between pi0(alpha=0) and \partial pi0(alpha=0)
    double Epot(ratio * epsilon), Ekin((1 - ratio) * epsilon);

    // INTEGRATING THE ODE SYSTEMS
    InitialConditions_Pi0 ic_pi0 = getInitialCondition_Pi0(sqrt(2*Epot)/m_pion, sqrt(2*Ekin));
    InitialConditions_PiPlus ic_piplus = getInitialCondition_PiPlus(0.6);
    double n(1.0);
    ic_piplus.piplus = [&](double alpha) { return sqrt(n) * (cos(ic_piplus.theta(alpha)) + 1i * sin(ic_piplus.theta(alpha))); };

    // SAVE TO FILE
    writeFuncToFile("data/pi0_initial.txt", ic_pi0.pi0, 0, M_PI / 2.0, 1000,{"alpha","pi0Re","pi0Im"});
    writeFuncToFile("data/chi_initial.txt", ic_pi0.chi, 0, M_PI / 2.0, 1000,{"alpha","chiRe","chiIm"});
    std::function<std::complex<double>(double)> Dpi0 = [&](double alpha)
    { return ic_pi0.chi(alpha) * (-Dr(alpha) * utau(alpha) + Dtau(alpha) * ur(alpha)); };
    writeFuncToFile("data/Dpi0_initial.txt", Dpi0, 0, M_PI / 2.0, 1000,{"alpha","Dpi0Re","Dpi0Im"});
    writeFunctionsToFile<std::complex<double>>(
        "data/initialfields_pi0.csv",
        {ic_pi0.pi0, Dpi0},
        0, M_PI / 2.0, 1000,
        {"alpha","pi0Re","pi0Im","Dpi0Re","Dpi0Im"});

    writeFuncToFile<std::complex<double>>("data/piplus_initial.txt", ic_piplus.piplus, 0, M_PI / 2.0, 1000,{"alpha","piplusRe","piplusIm"});
    writeFuncToFile("data/theta_initial.txt", ic_piplus.theta, 0, M_PI / 2.0, 1000,{"alpha","thetaRe","thetaIm"});
    std::function<std::complex<double>(double)> Dpiplus = [&](double alpha)
    { return 1i * ic_piplus.piplus(alpha) * m_pion * (-Dr(alpha) * utau(alpha) + Dtau(alpha) * ur(alpha)); }; // THIS ASSUMES CHI = MPION
    writeFuncToFile<std::complex<double>>("data/Dpiplus_initial.txt", Dpiplus, 0, M_PI / 2.0, 1000,{"alpha","DpiplusRe","DpiplusIm"});
    writeFunctionsToFile<std::complex<double>>(
        "data/initialfields_piplus.csv",
        {ic_piplus.piplus, Dpiplus},
        0, M_PI / 2.0, 1000,
        {"alpha","piplusRe","piplusIm","DpiplusRe","DpiplusIm"});

    return 0;
}