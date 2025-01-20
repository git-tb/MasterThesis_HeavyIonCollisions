#include <iostream>
#include <cmath>
#include <complex>
#include <functional>
#include <boost/math/interpolators/pchip.hpp>
#include <boost/program_options.hpp>
#include <filesystem>
#include <ctime>   // std::time, std::local_time
#include <iomanip> // std::put_time

#include <algorithm>

#include "savedata.h"       // write to file
#include "filereader.h"     // readcsv
#include "constants.h"      // GeV to inverse fm
#include "cuba.h"

#include "debugmsg.h"

// DEFINE HELPER FUNCTIONS
double w(double p, double m) { return sqrt(p*p + m*m);};
double absgradg(double t, double u, double v, double p, double ma, double mb, double Q) {
    return sqrt(
        pow(-t*u*Q*Q*w(p,mb)/w(Q*t,ma) + Q*p*v,2)
        + pow(-w(Q*t,ma)*w(p,mb),2)
        + pow(t*Q*p,2)
    );
};
double u_star(double t, double v, double p, double ma, double mb, double E_abc, double Q){
    return (ma * E_abc + t*Q*p*v)/(w(Q*t,ma) * w(p,mb));
};
double areafactor(double t, double v, double p, double ma, double mb, double E_abc, double Q) {
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

// THIS IS THE FULL INTEGRAND UP TO A FACTOR of B*Q*Q*ma/(p_abc * M_PI);
double myintegrand(std::function<double(double)> primespec, double t, double v, double p, double ma, double mb, double E_abc, double Q){
    double u = u_star(t,v,p, ma, mb, E_abc, Q);
    if(u <= 1) return 0.0;
    return (1/sqrt(u*u-1)) * (1/sqrt(1-v*v)) * areafactor(t,v,p,ma,mb,E_abc,Q)/absgradg(t,u,v,p,ma,mb,Q) * t * primespec(t*Q);
};

// TO CAST THE INTEGRAND INTO THE FORM REQUIREY THE INTEGRATION LIBRARY, INTRODUCE A STRUCT THAT HOLDS ALL PARAMETERS EXCEPT 
//  THE INTEGRATION VARIABLES t AND v

struct argsCUBA
{
    // MAYBE THINK ABOUT WHETHER TO DELETE TO DEFAULT CONSTRUCTOR FOR SAFETY...THOUGH IT MAKES THINGS A LITTLE MORE MESSY
    // argsCUBA() = delete; // NOT ACCIDENTAL INITIALIZATION
    // argsCUBA(double p_, double tmin_, double tmax_, double vmin_, double vmax_, std::function<double(double)> primespec_, double ma_, double mb_, double E_abc_, double Q_) :
    //    p(p_), tmin(tmin_), tmax(tmax_), vmin(vmin_), vmax(vmax_), primespec(primespec_), ma(ma_), mb(mb_), E_abc(E_abc_), Q(Q_)
    // {};

    // MAYBE IT'S EASIER TO UPDATE ALL RELEVANT VARIABLES IN A FUNCTION CALL IF NECESSARY
    void update(double p_, double tmin_, double tmax_, double vmin_, double vmax_) {
        p = p_;
        tmin = tmin_;
        tmax = tmax_;
        vmin = vmin_;
        vmax = vmax_;
    }

    double p, tmin, tmax, vmin, vmax;
    std::function<double(double)> primespec;
    double ma, mb, E_abc, Q;
};
// struct argsCUBATURE
// {
//     double p;
// };

static int integrandCUBA(
    const int *ndim, const double xx[],
    const int *ncomp, double ff[], void *userdata)
{    
    argsCUBA* myargs = (struct argsCUBA*)userdata;

    // CUBA INTEGRATES OVER [0,1]^D, THEREFORE WE NEED TO SCALE THE ARGUMENTS
    // BUT ALSO THE INTEGRAL MEASURE (OR EQUIVALENTLY, THE INTEGRAND)
    double t = myargs->tmin + (myargs->tmax - myargs->tmin) * xx[0];
    double v = myargs->vmin + (myargs->vmax - myargs->vmin) * xx[1];

    ff[0] = (myargs->tmax - myargs->tmin) * (myargs->vmax - myargs->vmin) * myintegrand(myargs->primespec,
                                                                                        t,
                                                                                        v,
                                                                                        myargs->p,
                                                                                        myargs->ma,
                                                                                        myargs->mb,
                                                                                        myargs->E_abc,
                                                                                        myargs->Q);

    return 0;
}
// int integrandCUBATURE(
//         unsigned ndim, const double *x,
//         void *fdata,
//         unsigned fdim, double *fval)
// {
//     double t(x[0]), v(x[1]);
//     argsCUBATURE myargs = *(struct argsCUBATURE *)fdata;
//     fval[0] = restrictfunc(t,v,myargs.p);
//     return 0;
// }

std::vector<double> decayspec(
    std::vector<double> pTs,
    double ma,
    double mb,
    double mc,
    double B,
    std::function<double(double)> primespec,
    double qmax,
    double Q,
    double epsabs,
    double epsrel,
    int iterations,
    std::function<void(double, double)> callback= [](double p, double val){ return; })
{
    double  p_abc = 1/(2*ma) * sqrt( (pow(ma+mb,2)-pow(mc,2)) * (pow(ma-mb,2)-pow(mc,2)) ),
            E_abc = sqrt(mb*mb + p_abc*p_abc);
    std::vector<double> finalspec(pTs.size());

    // argsCUBATURE myargsCUBATURE;
    argsCUBA myargsCUBA;
    myargsCUBA.E_abc = E_abc;
    myargsCUBA.ma = ma;
    myargsCUBA.mb = mb;
    myargsCUBA.primespec = primespec;
    myargsCUBA.Q = Q;
    myargsCUBA.tmax = qmax/Q;
    myargsCUBA.tmin = 0;
    myargsCUBA.vmax = 1;
    myargsCUBA.vmin = -1;
    myargsCUBA.p = 0.0;

    const int   NDIM(2), 
                NCOMP(1), 
                NVEC(1), 
                FLAGS(0), 
                MINEVAL(0), 
                MAXEVAL(iterations), 
                KEY(0);
    void* USERDATA(&myargsCUBA);
    void* SPIN(NULL);
    const double EPSREL(epsrel), EPSABS(epsabs);
    const char* STATEFILE(NULL);

    int nregions, neval, fail;
    double integral[NCOMP], error[NCOMP], prob[NCOMP];
    for (int i = 0; i < pTs.size(); i++)
    {
        std::cout << i + 1 << " / " << pTs.size() << std::endl;

        double vmin(-1);
        double tmin(0);
        double vmax(1);
        double tmax(qmax/Q);
            
        if(-E_abc*E_abc + pTs[i]*pTs[i] + mb*mb >= 0) // if w_p >= E_abc
            vmin = sqrt(-E_abc*E_abc + pTs[i]*pTs[i] + mb*mb)/pTs[i];
        // if(E_abc >= mb) // this is always true
        // tmin = ma * (E_abc*pTs[i] - w(pTs[i],mb) * sqrt(E_abc*E_abc - mb*mb))/(mb*mb);
        // tmax = ma * (E_abc*pTs[i] + w(pTs[i],mb) * sqrt(E_abc*E_abc - mb*mb))/(mb*mb);
        tmin = ma * (E_abc*pTs[i] - w(pTs[i],mb) * p_abc)/(Q*mb*mb);
        tmax = ma * (E_abc*pTs[i] + w(pTs[i],mb) * p_abc)/(Q*mb*mb);

        tmin = tmin > 0 ? tmin : 0;
        tmax = tmax < qmax/Q ? tmax : qmax/Q;

        myargsCUBA.update(pTs[i], tmin, tmax, vmin, vmax);

        // ------------------ WITH CUBATURE LIBRARY -------------------------
        // const double    XMIN[2] = {tmin,vmin},
        //                 XMAX[2] = {tmax,vmax};

        // myargsCUBATURE.p = pTs[i];
        // double val, err;
        // hcubature(FDIM, integrandCUBATURE, &myargsCUBATURE, XDIM, XMIN, XMAX, MAXEVAL, ERRABS, ERRREL, ERROR_INDIVIDUAL, &val, &err);
        // finalspec[i] = val;


        // ------------------ WITH CUBA LIBRARY -------------------------
        // This performs the numerical integraion
        Cuhre(NDIM, NCOMP, integrandCUBA, USERDATA, NVEC,
            EPSREL, EPSABS, FLAGS,
            MINEVAL, MAXEVAL, KEY,
            STATEFILE, SPIN,
            &nregions, &neval, &fail, integral, error, prob);

        finalspec[i] = B * integral[0] * Q*Q*ma / (p_abc * M_PI);

        callback(pTs[i],finalspec[i]);
    }
    return finalspec;
}

int main(int ac, char* av[])
{
    /* #region COMMAND LINE OPTIONS */
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
        ("epsabs",          po::value<double>()->default_value(0),      "absolute integration error goal")
        ("epsrel",          po::value<double>()->default_value(1e-3),   "relative integration error goal")
        ("iter",            po::value<int>()->default_value(1e4),       "maximum integration iterations")
        ("primespecpath",po::value<std::string>(),"csv file containing primary spectrum")
        ("parentdir",       po::value<std::string>()->default_value("Data"),
            "data folder, in which a subfolder for the results is created")
        ("B", po::value<double>()->default_value(1.0), "branching ratio of the decay")
        ("Q", po::value<double>()->default_value(1.0), "dimensionful pT-scale (should not influence the result)")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);

    // PROCESS CMDLINE OPTIONS
    if (vm.count("help"))    {
        std::cout << desc << "\n";
        return 1;
    }

    DEBUGMSG("start command line processing");
    double  ma                  = vm["ma"].as<double>(),
            mb                  = vm["mb"].as<double>(),
            mc                  = vm["mc"].as<double>(),
            pTmax               = vm["pTmax"].as<double>(),
            B                   = vm["B"].as<double>(),
            Q                   = vm["Q"].as<double>(),
            epsabs              = vm["epsabs"].as<double>(),
            epsrel              = vm["epsrel"].as<double>();
    int NpTs                    = vm["NpT"].as<int>(),
        iter                    = vm["iter"].as<int>();
    std::string primespecpath   = vm["primespecpath"].as<std::string>(),
                parentdir       = vm["parentdir"].as<std::string>();
    double  p_abc               = 1/(2*ma) * sqrt( (pow(ma+mb,2)-pow(mc,2)) * (pow(ma-mb,2)-pow(mc,2)) ),
            E_abc               = sqrt(mb*mb + p_abc*p_abc);
    DEBUGMSG("command line processing completed");
    /* #endregion */

    // CREATE DIRECTORY TO SAVE FILES
    std::time_t t = std::time(nullptr);
    std::tm tm = *std::localtime(&t);
    std::stringstream timestamp_sstr;
    timestamp_sstr << std::put_time(&tm, "%Y%m%d_%H%M%S");
    std::string timestamp = timestamp_sstr.str();
    std::string pathname = parentdir+"/decayspec_"+timestamp;
    // std::filesystem::create_directories(pathname); // create later, such that in case of read-in error & abort no empty directory is created

    // INTERPOLATE PRIMESPEC FROM DATA
    csvdata primespecdata = readcsv(primespecpath);
    std::vector<double> x(primespecdata.data[0]);
    std::vector<double> y(primespecdata.data[1]);
    double qmin(0), qmax(x[x.size()-1]);

    std::vector<double> log10y(y.size()); // logarithmics data is probably interpolated better
    for(int i = 0; i < y.size(); i++) log10y[i] = log10(y[i]);

    using boost::math::interpolators::pchip;
    auto logspecspline = pchip(
        std::move(x),
        std::move(log10y));
    std::function<double(double)> primespec = [&](double p){ return pow(10, logspecspline(p)); };

    int NSAMPLE(1000);
    std::filesystem::create_directories(pathname); // CREATE DIRECTORY HERE
    writeFuncToFile(pathname+"/primespec.txt", primespec, qmin, qmax, NSAMPLE, {"q","primespecRe","primespecIm"},{timestamp});

    // AGAIN, WRITE TO FILE DURING COMPUTATION, SO PREPARE THE FILE BEFOREHAND
    std::stringstream qmax_ss, primespec_ss, userinput_ss;
    qmax_ss         << "qmax:\t" << qmax;
    primespec_ss    << "primespec:\t" << primespecpath;
    userinput_ss    << "ma:\t" << ma << "\n# "  // the other comments are just 1 line per comment, therefore add the comment symbol "#" manually....
                    << "mb:\t" << mb << "\n# "
                    << "mc:\t" << mc << "\n# "
                    << "[=>] pabc:\t" << p_abc << "\n# "
                    << "[=>] Eabc:\t" << E_abc << "\n# "
                    << "B:\t" << B << "\n# "
                    << "Q:\t" << Q;
    std::vector<std::string> comments({
        timestamp,
        primespec_ss.str(),
        qmax_ss.str(),
        userinput_ss.str(),
    });
    std::vector<std::string> headers({"p","finalspecRe","finalspecIm"});

    std::string path = pathname+"/decayspec.txt";
    std::ofstream decayspec_output(path);
    if (!decayspec_output.is_open()) { std::cerr << "Error opening the file: " << path << " to save to" << std::endl; return -1; }
    for(int i = 0; i < comments.size(); i++) decayspec_output << "# " << comments[i] << std::endl;
    for(int i = 0; i < headers.size(); i++) { decayspec_output << headers[i]; if(i != headers.size()-1) decayspec_output << ","; }
    decayspec_output << std::endl;

    // COMPUTE DECAY SPEC
    std::vector<double> pTs(NpTs);
    for (int i = 0; i < pTs.size(); i++) pTs[i] = i * pTmax / (NpTs - 1);

    std::function<void(double,double)> callback = [&decayspec_output](double p, double value)
    {
        decayspec_output << p << "," << std::real(value) << "," << std::imag(value) << std::endl;
    }; 

    std::vector<double> finalspec = decayspec(pTs, ma, mb, mc, B, primespec, qmax, Q, epsabs, epsrel, iter, callback);
    // for(int i = 0; i < pTs.size(); i++)
    // {
    //     decayspec_output << pTs[i] << "," << std::real(finalspec[i]) << "," << std::imag(finalspec[i]) << std::endl;
    // }
    decayspec_output.close();
}