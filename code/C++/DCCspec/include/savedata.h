#ifndef SAVEDATA_H
#define SAVEDATA_H

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <string>
#include <functional>

template <typename T>
void writeSamplesToFile(std::string path, std::vector<double> x, std::vector<T> y, std::vector<std::string> headers)
{
    std::ofstream output(path);

    if (!output.is_open())
    {
        std::cerr << "Error opening the file: " << path << std::endl;
        return;
    }

    for(int i = 0; i < headers.size(); i++)
    {
        output << headers[i];
        if(i != headers.size()-1) output << ",";
    }
    output << std::endl;

    for (int i = 0; i < x.size(); i++)
    {
        output << x[i] << "," << std::real(y[i]) << "," << std::imag(y[i]) << std::endl;
    }

    output.close();
}

template <typename T>
void writeFuncToFile(std::string path, std::function<T(double)> func, double a, double b, int Nsamples, std::vector<std::string> headers)
{
    std::ofstream output(path);

    if (!output.is_open())
    {
        std::cerr << "Error opening the file: " << path << std::endl;
        return;
    }

    for(int i = 0; i < headers.size(); i++)
    {
        output << headers[i];
        if(i != headers.size()-1) output << ",";
    }
    output << std::endl;

    double dx = (b - a) / (Nsamples - 1);
    for (int i = 0; i < Nsamples; i++)
    {
        T y = func(i * dx);
        output << i * dx << "," << std::real(y) << "," << std::imag(y) << std::endl;
    }

    output.close();
}

template <typename T>
void writeFunctionsToFile(std::string path, std::vector<std::function<T(double)>> funcs, double a, double b, int Nsamples, std::vector<std::string> headers)
{
    std::ofstream output(path);

    if (!output.is_open())
    {
        std::cerr << "Error opening the file: " << path << std::endl;
        return;
    }

    double dx = (b - a) / (Nsamples - 1);
    for(int i = 0; i < headers.size(); i++)
    {
        output << headers[i];
        if(i != headers.size()-1) output << ",";
    }
    output << std::endl;
    for (int i = 0; i < Nsamples; i++)
    {
        output << i * dx ;
        for(int j = 0; j < funcs.size(); j++)
        {
            std::function<T(double)> func = funcs[j];
            T y = func(i * dx);
            output << "," << std::real(y) << "," << std::imag(y);
        }
        output << std::endl;
    }

    output.close();
}

#endif