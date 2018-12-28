#include <reference_reader.hpp>
#include <string>
#include <fstream>
#include <iostream>

std::string ReferenceReader::read_reference(std::string path)
{
    std::ifstream file(path);
    std::string reference;

    if (file.is_open())
    {
        std::string line;
        getline(file, line);
        while (getline(file, line))
        {
            reference += line.c_str();
        }
        file.close();
    }

    return reference;
}