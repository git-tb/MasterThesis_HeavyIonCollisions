#ifndef PROCESSDATA_H
#define PRCOESSDATA_H

#include <iostream>
#include <functional>       // std::function
#include <string>           // std::string
#include <cmath>            // M_PI
#include <complex>          // std::complex<double>
#include <boost/math/interpolators/pchip.hpp>   // pchip interpolation

#include "filereader.h"     // csvdata, readcsv()
#include "extrapolate.h"    // linearExtrapolate()

using namespace std::complex_literals; // imaginary unit is 1i

std::function<double(double)> UNINITIALIZED_FUNCTION = [](double x){
    std::cout << "THIS FUNCTION HAS NOT BEEN INITIALIZED" << std::endl;
    return 0.0;
};
std::function<double(double)> 
    tau(UNINITIALIZED_FUNCTION),
    r(UNINITIALIZED_FUNCTION),
    Dtau(UNINITIALIZED_FUNCTION),
    Dr(UNINITIALIZED_FUNCTION),
    ur(UNINITIALIZED_FUNCTION),
    utau(UNINITIALIZED_FUNCTION);
csvdata ProcessFreezeoutData(std::string datafile = "./../../Mathematica/data/ExampleFreezeOut.csv")
{
    std::ofstream debug("debug_freezout.txt");

    // READ IN DATA
    csvdata mydata = readcsv(datafile, ",", true, true); // FLUID INFO
    for (int i = 0; i < mydata.header.size(); i++)    {
        debug << mydata.header[i] << (i == mydata.header.size() - 1 ? "" : ", ") << std::endl;
        for (int j = 0; j < mydata.data[i].size(); j++)
            debug << mydata.data[i][j] << std::endl;
        debug << std::endl;
    }

    // INTERPOLATE
    std::vector<double> x = mydata.data[0];
    std::vector<double> x1(x), x2(x), x3(x), x4(x), x5(x), x6(x);
    std::vector<double> 
        y1(mydata.data[1]),
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

    tau = std::move(tau_spline);    // std::move doesn't seem to be necessary, but it in principle it makes sense because we don't
                                        //  want the splines to be destroyed at the end of the scope
    r = std::move(r_spline);
    Dtau = std::move(Dtau_spline);
    Dr = std::move(Dr_spline);
    ur = std::move(ur_spline);
    utau = std::move(utau_spline);

    return mydata;
}

std::function<double(double)>
    f0Re(UNINITIALIZED_FUNCTION),
    f0Im(UNINITIALIZED_FUNCTION),
    Df0Re(UNINITIALIZED_FUNCTION),
    Df0Im(UNINITIALIZED_FUNCTION);
std::function<std::complex<double>(double)> f0 = [](double alpha) { return f0Re(alpha) + 1i * f0Im(alpha);};
std::function<std::complex<double>(double)> Df0 = [](double alpha) { return Df0Re(alpha) + 1i * Df0Im(alpha);};
csvdata ProcessInitialData(std::string datafile = "./data/initialfields.csv")
{ 
    std::ofstream debug("debug_initial.txt");

    // READ DATA
    csvdata initialdata = readcsv(datafile, ",", true, true);              // INITIAL FIELDS
    for (int i = 0; i < initialdata.header.size(); i++)    {
        debug << initialdata.header[i] << (i == initialdata.header.size() - 1 ? "" : ", ") << std::endl;
        for (int j = 0; j < initialdata.data[i].size(); j++)
            debug << initialdata.data[i][j] << std::endl;
        debug << std::endl;
    }
    debug.close();

    // INTERPOLATE
    std::vector<double> x = initialdata.data[0];
    std::vector<double> x1(x), x2(x), x3(x), x4(x);
    std::vector<double> 
        y1(initialdata.data[1]),
        y2(initialdata.data[2]),
        y3(initialdata.data[3]),
        y4(initialdata.data[4]);

    linearExtrapolate(x1, y1, 0 - 0.01);
    linearExtrapolate(x1, y1, M_PI / 2.0 + 0.01);
    linearExtrapolate(x2, y2, 0 - 0.01);
    linearExtrapolate(x2, y2, M_PI / 2.0 + 0.01);
    linearExtrapolate(x3, y3, 0 - 0.01);
    linearExtrapolate(x3, y3, M_PI / 2.0 + 0.01);
    linearExtrapolate(x4, y4, 0 - 0.01);
    linearExtrapolate(x4, y4, M_PI / 2.0 + 0.01);

    using boost::math::interpolators::pchip;
    auto f0Re_spline = pchip(
        std::move(x1),
        std::move(y1));
    auto f0Im_spline = pchip(
        std::move(x2),
        std::move(y2));
    auto Df0Re_spline = pchip(
        std::move(x3),
        std::move(y3));
    auto Df0Im_spline = pchip(
        std::move(x4),
        std::move(y4));
            
    f0Re = std::move(f0Re_spline);
    f0Im = std::move(f0Im_spline);
    Df0Re = std::move(Df0Re_spline);
    Df0Im = std::move(Df0Im_spline);

    return initialdata;
}

#endif