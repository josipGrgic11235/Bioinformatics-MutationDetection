#include <kmer_index.hpp>
#include <utility>

#ifndef READ_MAPPER
#define READ_MAPPER

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
    ReadMapper(KmerIndex &kmerIndex, int region_divider, int k) : kmerIndex(kmerIndex), region_divider(region_divider), k(k){};

    Region map(std::string &input);

    int get_region_divider()
    {
        return region_divider;
    }
};

#endif