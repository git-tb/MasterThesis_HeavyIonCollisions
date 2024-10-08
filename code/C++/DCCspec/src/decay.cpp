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
#include <limits>   // signaling NaN

#include <algorithm>
#include <iterator>

#include "savedata.h"       // write to file
#include "processdata.h"    // initial conditions and freezeout functions
#include "constants.h"      // pion mass, GeV to inverse fm
#include "cubature.h"
#include "cuba.h"

double B(1.0);
double Q(1.0);
double ma(0.6), mb(0.14), mc(0.14);
double p_abc(std::numeric_limits<double>::signaling_NaN());
double E_abc(std::numeric_limits<double>::signaling_NaN());

std::function<double(double,double)> w = [](double p, double m) { return sqrt(p*p + m*m);};
std::function<double(double, double,double,double)> absgradg = [](double t, double u, double v, double p) {
    return sqrt(
        pow(-t*u*Q*Q*w(p,mb)/w(Q*t,ma) + Q*p*v,2)
        + pow(-w(Q*t,ma)*w(p,mb),2)
        + pow(t*Q*p,2)
    );
};
std::function<double(double,double,double)> u_star = [](double t, double v, double p){
    return (ma * E_abc + t*Q*p*v)/(w(Q*t,ma) * w(p,mb));
};
std::function<double(double, double,double)> areafactor = [](double t, double v, double p) {
    return sqrt(
        pow(
            ma*Q*(-t*Q*E_abc + v*ma*p) / (w(p,mb) * pow(w(t*Q,ma),3))
            ,2)
        +pow(-1,2)
        +pow(
            t*Q*p/(w(Q*t,ma)*w(p,mb))
            ,2)
    );
};

std::function<double(double)> primespec = [](double x){ std::cout << "PRIMARY SPECTRUM FUNCTION UNINITIALIZED!"; return 0.0;};
std::function<double(double,double,double)> restrictfunc = [](double t, double v, double p){
    // double u = (1 + t*v)/sqrt(2+t*t);
    // if(u <= 1)  return 0.0;
    // else        return 1/(sqrt(u*u-1) * sqrt(1-v*v));

    double u = u_star(t,v,p);
    if(u <= 1) return 0.0;
    return (1/sqrt(u*u-1)) * (1/sqrt(1-v*v)) * areafactor(t,v,p)/absgradg(t,u,v,p) * t * primespec(t*Q);
    // return (1/sqrt(u*u-1)) * (1/sqrt(1-v*v)) * areafactor(t,v,p)/absgradg(t,u,v,p);
};

struct argsCUBATURE
{
    double p;
};

struct argsCUBA
{
    double p, tmin, tmax, vmin, vmax;
};

static int integrandCUBA(
    const int *ndim, const double xx[],
    const int *ncomp, double ff[], void *userdata)
{    
    argsCUBA myargs = *(struct argsCUBA*)userdata;

    double t = myargs.tmin + (myargs.tmax - myargs.tmin) * xx[0];
    double v = myargs.vmin + (myargs.vmax - myargs.vmin) * xx[1];

    ff[0] = (myargs.tmax - myargs.tmin) * (myargs.vmax - myargs.vmin) * restrictfunc(t,v,myargs.p);

    return 0;
}

int integrandCUBATURE(
        unsigned ndim, const double *x,
        void *fdata,
        unsigned fdim, double *fval)
{
    double t(x[0]), v(x[1]);
    argsCUBATURE myargs = *(struct argsCUBATURE *)fdata;
    fval[0] = restrictfunc(t,v,myargs.p);
    return 0;
}

std::vector<double> decayspec(
    std::vector<double> ps,
    std::function<double(double)> primespec,
    double qmax,
    std::function<void(double, double)> callback= [](double p, double val){ return; })
{
    std::vector<double> finalspec(ps.size());

    double ERRABS(1e-7), ERRREL(1e-3);
    int MAXEVAL(1e6), FDIM(1), XDIM(2);
    argsCUBATURE myargsCUBATURE;
    argsCUBA myargsCUBA;

    for (int i = 0; i < ps.size(); i++)
    {
        std::cout << i + 1 << " / " << ps.size() << std::endl;

        double vmin(-1);
        double tmin(0);
        double tmax(qmax/Q);

        // double  A(ma*E_abc),
        //         B(ps[i]*Q),
        //         C(Q),
        //         D(ma),
        //         F(ps[i]),
        //         G(mb);

        // if(-A*A*D*D + D*D*D*D*(F*F + G*G) >= 0)
        //     vmin = A*C*(-A + sqrt(F*F + G*G) * sqrt(D*D*D*D*(F*F + G*G)/(A*A)))/
        //                 (B*sqrt(-A*A*D*D + D*D*D*D*(F*F + G*G)));
            
        if(-E_abc*E_abc + ps[i]*ps[i] + mb*mb >= 0)
            vmin = sqrt(-E_abc*E_abc + ps[i]*ps[i] + mb*mb)/ps[i];
        // if(E_abc >= mb) // this is always true
        // tmin = ma * (E_abc*ps[i] - w(ps[i],mb) * sqrt(E_abc*E_abc - mb*mb))/(mb*mb);
        // tmax = ma * (E_abc*ps[i] + w(ps[i],mb) * sqrt(E_abc*E_abc - mb*mb))/(mb*mb);
        tmin = ma * (E_abc*ps[i] - w(ps[i],mb) * p_abc)/(mb*mb);
        tmax = ma * (E_abc*ps[i] + w(ps[i],mb) * p_abc)/(mb*mb);

        tmin = tmin > 0 ? tmin : 0;
        tmax = tmax < qmax/Q ? tmax : qmax/Q;

        // ------------------ WITH CUBATURE LIBRARY -------------------------
        // const double    XMIN[2] = {0,vmin},
        //                 XMAX[2] = {qmax/Q,1};

        // myargsCUBATURE.p = ps[i];
        // double val, err;
        // hcubature(FDIM, integrandCUBATURE, &myargsCUBATURE, XDIM, XMIN, XMAX, MAXEVAL, ERRABS, ERRREL, ERROR_INDIVIDUAL, &val, &err);
        // finalspec[i] = val;


        // ------------------ WITH CUBA LIBRARY -------------------------
        const int   NDIM(2), 
                    NCOMP(1), 
                    NVEC(1), 
                    FLAGS(0), 
                    MINEVAL(0), 
                    MAXEVAL(1e6), 
                    KEY(0);
        void* USERDATA(&myargsCUBA);
        void* SPIN(NULL);
        const double EPSREL(1e-7), EPSABS(1e-3);
        const char* STATEFILE(NULL);

        int nregions, neval, fail;
        double integral[NCOMP], error[NCOMP], prob[NCOMP];

        myargsCUBA.tmin = tmin;
        myargsCUBA.tmax = tmax;
        myargsCUBA.vmin = vmin;
        myargsCUBA.vmax = 1;
        myargsCUBA.p = ps[i];

        // This performs the numerical integraion
        Cuhre(NDIM, NCOMP, integrandCUBA, USERDATA, NVEC,
            EPSREL, EPSABS, FLAGS,
            MINEVAL, MAXEVAL, KEY,
            STATEFILE, SPIN,
            &nregions, &neval, &fail, integral, error, prob);

        finalspec[i] = B * integral[0] * Q*Q*ma / (p_abc * M_PI);
        callback(ps[i],finalspec[i]);
    }
    return finalspec;
}

int main(int ac, char* av[])
{
    /* #region COMMAND LINE OPTIONS */
    // DEFINE VARIABLES TO BE SET
    double pmin(0), pmax(1.0);
    int Nps(100);
    std::string primespecpath;

    // DECLARE SUPPORTED OPTIONS
    namespace po = boost::program_options;
    po::options_description desc(   "// ============================================================= \\\\\n"
                                    "|| This program computes the induced pT-spectrum of a particle b ||\n"
                                    "|| from the decay of a primary particle a via the decay process  ||\n"
                                    "||                       a -> b +c                               ||\n"
                                    "\\\\ ============================================================= //\n\n"
                                    "Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("ma", po::value<double>()->default_value(0.8), "mass of particle a")
        ("mb", po::value<double>()->default_value(0.14), "mass of particle b")
        ("mc", po::value<double>()->default_value(0.14), "mass of particle c")
        ("pTmax", po::value<double>()->default_value(2.0), "spectrum is computed on [0,pTmax]")
        ("NpT", po::value<int>()->default_value(200), "number of sample points within [0,pTmax]")
        ("primespecpath",po::value<std::string>(),"csv file containing primary spectrum")
        ("B", po::value<double>()->default_value(1.0), "branching ratio of the decay")
        ("Q", po::value<double>()->default_value(1.0), "dimensionless scale (should not influence the result)")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);

    // PROCESS CMDLINE OPTIONS
    if (vm.count("help"))    {
        std::cout << desc << "\n";
        return 1;
    }
    ma = vm["ma"].as<double>();
    mb = vm["mb"].as<double>();
    mc = vm["mc"].as<double>();
    pmax = vm["pTmax"].as<double>();
    Nps = vm["NpT"].as<int>();
    primespecpath = vm["primespecpath"].as<std::string>();
    B = vm["B"].as<double>();
    Q = vm["Q"].as<double>();

    p_abc = 1/(2*ma) * sqrt( (pow(ma+mb,2)-pow(mc,2)) * (pow(ma-mb,2)-pow(mc,2)) );
    E_abc = sqrt(mb*mb + p_abc*p_abc);
    /* #endregion */

    // CREATE DIRECTORY TO SAVE FILES
    std::time_t t = std::time(nullptr);
    std::tm tm = *std::localtime(&t);
    std::stringstream timestamp_sstr;
    timestamp_sstr << std::put_time(&tm, "%Y%m%d_%H%M%S");
    std::string timestamp = timestamp_sstr.str();
    std::string pathname = "data/decayspec_"+timestamp;
    // std::filesystem::create_directories(pathname); // create later, such that in case of read-in error & abort no empty directory is created

    // INTERPOLATE PRIMESPEC FROM DATA
    csvdata primespecdata = readcsv(primespecpath);
    std::vector<double> x(primespecdata.data[0]);
    std::vector<double> y(primespecdata.data[1]);
    double qmin(0), qmax(x[x.size()-1]);

    std::vector<double> log10y(y.size()); // logarithmics data is probably interpolated better
    for(int i = 0; i < y.size(); i++) log10y[i] = log10(y[i]);
    linearExtrapolate(x, log10y, 0 - 0.01);

    using boost::math::interpolators::pchip;
    auto logspecspline = pchip(
        std::move(x),
        std::move(log10y));
    primespec = [&](double p){ return pow(10, logspecspline(p)); };

    int NSAMPLE(1000);
    std::filesystem::create_directories(pathname); // CREATE DIRECTORY HERE
    writeFuncToFile(pathname+"/primespec_interp.txt", primespec, qmin, qmax, NSAMPLE, {"q","primespecRe","primespecIm"},{timestamp});

    // AGAIN, WRITE TO FILE DURING COMPUTATION, SO PREPARE THE FILE BEFOREHAND
    std::stringstream qmax_ss, primespec_ss, userinput_ss;
    qmax_ss << "qmax:\t" << qmax;
    primespec_ss << "primespec:\t" << primespecpath;
    userinput_ss    << "ma:\t" << ma << "\n# "  // the other comments are just 1 line per comment, therefore add the comment symbol "#" manually....
                    << "mb:\t" << mb << "\n# "
                    << "mc:\t" << mc << "\n# "
                    << "B:\t" << B << "\n# "
                    << "Q:\t" << Q << "\n# "
                    << "[DEBUG] pabc:\t" << p_abc << "\n# "
                    << "[DEBUG] Eabc:\t" << E_abc;
    std::vector<std::string> comments({
        timestamp,
        primespec_ss.str(),
        qmax_ss.str(),
    });
    std::vector<std::string> headers({"p","finalspecRe","finalspecIm"});

    std::string path = pathname+"/decayspec.txt";
    std::ofstream decayspec_output(path);
    if (!decayspec_output.is_open()) { std::cerr << "Error opening the file: " << path << " to save to" << std::endl; return -1; }
    for(int i = 0; i < comments.size(); i++) decayspec_output << "# " << comments[i] << std::endl;
    for(int i = 0; i < headers.size(); i++) { decayspec_output << headers[i]; if(i != headers.size()-1) decayspec_output << ","; }
    decayspec_output << std::endl;

    // COMPUTE DECAY SPEC
    std::vector<double> ps(Nps);
    for (int i = 0; i < ps.size(); i++) ps[i] = pmin + i * (pmax - pmin) / (Nps - 1);

    std::function<void(double,double)> callback = [&decayspec_output](double p, double value)
    {
        decayspec_output << p << "," << std::real(value) << "," << std::imag(value) << std::endl;
    }; 
    std::vector<double> finalspec = decayspec(ps, primespec,qmax,callback);
    
    // std::stringstream qmax_ss, primespec_ss;
    // qmax_ss << "gmax:\t" << qmax;
    // primespec_ss << "primespec:\t" << primespecpath;
    // writeSamplesToFile(pathname+"/decayspec.txt", ps, finalspec, {"p","finalspecRe","finalspecIm"},
    // {
    //     timestamp,
    //     primespec_ss.str(),
    //     qmax_ss.str(),
    // });
}