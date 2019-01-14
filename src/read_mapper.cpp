#include <read_mapper.hpp>
#include <numeric>
#include <string>

Region ReadMapper::map(std::string &input)
{
    std::vector<int> region(kmerIndex.get_reference_size() / region_divider + 1, 0);

    for (int i = 0; i <= input.size() - k; i += 2)
    {
        std::string kmer = input.substr(i, k);
        auto index_list = kmerIndex.find(kmer);

        for (auto i = index_list.begin(); i != index_list.end(); i++)
        {
            region[*i / region_divider]++;
        }
    }

    int sequence_region_length = (input.size() / region_divider) + 1;

    int max = 0;
    int region_id = 0;
    for (int i = 0; i <= region.size() - sequence_region_length; i++)
    {
        int region_count = std::accumulate(region.begin() + i, region.begin() + i + sequence_region_length, 0);
        if (region_count > max)
        {
            max = region_count;
            region_id = i;
        }
    }

    return Region(region_id, region_id + sequence_region_length);
}