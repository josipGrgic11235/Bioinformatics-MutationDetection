// Created by Leon

#include <reverse_complement.hpp>
#include <string>
#include <stdexcept>

char ReverseComplement::get_complement_base(char c)
{
    switch (c)
    {
    case 'T':
        return 'A';
    case 'A':
        return 'T';
    case 'G':
        return 'C';
    case 'C':
        return 'G';
    default:
        throw std::invalid_argument("Unknown base");
    }
}

std::string ReverseComplement::get_reverse_complement(std::string input)
{
    std::string reverse_complement;
    for (int i = input.size() - 1; i >= 0; i--)
    {
        reverse_complement += get_complement_base(input[i]);
    }
    return reverse_complement;
}