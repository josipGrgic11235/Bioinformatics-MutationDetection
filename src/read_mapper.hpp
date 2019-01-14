// Created by Leon

#include <kmer_index.hpp>
#include <utility>

#ifndef READ_MAPPER
#define READ_MAPPER

/**
 * Structure that represents region with start and end index.
 **/
struct Region
{
    int start_index;
    int end_index;

    Region(int start_index, int end_index) : start_index(start_index), end_index(end_index) {}
};

class ReadMapper
{
    KmerIndex &kmerIndex;
    const int region_divider;
    const int k;

  public:
    /**
     * Constructor
     * param: region_divider = reference genome region divider
     * param: k = k-mer length
     **/
    ReadMapper(KmerIndex &kmerIndex, int region_divider, int k) : kmerIndex(kmerIndex), region_divider(region_divider), k(k){};

    /**
     * Maps given string to the reference and returns start and end index of it's position.
     * param: input = string that we want to map
     * returns: Region = Region structure with start and end index where input is positionited in reference
     **/
    Region map(std::string &input);

    int get_region_divider()
    {
        return region_divider;
    }
};

#endif