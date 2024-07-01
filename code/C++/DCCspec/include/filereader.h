#ifndef FILEREADER_H
#define FILEREADER_H

#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <regex>

struct csvdata
{
    std::vector<std::string> header;
    std::vector<std::vector<double>> data;
};
csvdata readcsv(std::string source, std::string delimiter = ",", bool has_header = true, bool transpose = true)
{
    std::ifstream input(source);

    if (!input.is_open())
    {
        std::cerr << "Error opening the file!" << std::endl;
        return {{"ERROR"}, {{-1}}};
    }

    csvdata mydata;

    std::string line;

    if (has_header)
    {
        std::getline(input, line);
        size_t pos = 0;
        std::string token;
        while ((pos = line.find(delimiter)) != std::string::npos)
        {
            mydata.header.push_back(line.substr(0, pos));
            line.erase(0, pos + delimiter.length());
        }
        mydata.header.push_back(line.substr(0, pos));
    }

    while (std::getline(input, line))
    {
        mydata.data.push_back(std::vector<double>());

        size_t pos = 0;
        std::string token;
        while ((pos = line.find(delimiter)) != std::string::npos)
        {
            mydata.data.back().push_back(std::stof(line.substr(0, pos)));
            line.erase(0, pos + delimiter.length());
        }
        mydata.data.back().push_back(std::stof(line.substr(0, pos)));
    }
    input.close();

    if (transpose)
    {
        // transpose the data
        std::vector<std::vector<double>> transposeddata(mydata.data[0].size());
        for (int i = 0; i < mydata.data[0].size(); i++)
        {
            for (int j = 0; j < mydata.data.size(); j++)
            {
                transposeddata[i].push_back(mydata.data[j][i]);
            }
        }

        mydata.data = transposeddata;
    }

    return mydata;
}

#endif