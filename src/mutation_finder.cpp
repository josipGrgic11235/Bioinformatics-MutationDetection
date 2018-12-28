#include <mutation_finder.hpp>
#include <kmer_index.hpp>
#include <read_mapper.hpp>
#include <local_alignment.hpp>
#include <reverse_complement.hpp>
#include <execution_timer.hpp>
#include <iostream>
#include <algorithm>
#include <sstream>

template <typename KeyType, typename ValueType>
std::pair<KeyType, ValueType> get_max(const std::map<KeyType, ValueType> &x)
{
    using pairtype = std::pair<KeyType, ValueType>;
    return *std::max_element(x.begin(), x.end(), [](const pairtype &p1, const pairtype &p2) {
        return p1.second < p2.second;
    });
}

std::string MutationFinder::find_mutations(int k, int region_divider, int local_align_k, int match_score, int change_score, int gap_score, int confirmation_count)
{
    ExecutionTimer timer;
    timer.start();

    KmerIndex kmerIndex(reference, k);
    std::map<int, std::map<std::string, int>> result_map;
    ReadMapper readMapper(kmerIndex, region_divider, k);
    LocalAlignment localAlignment(readMapper, reference, result_map, local_align_k, region_divider, match_score, change_score, gap_score);

    for (int i = 0; i < read_data.size(); i++)
    {
        if (!localAlignment.align(read_data[i]))
        {
            std::string reverse_complement = ReverseComplement::get_reverse_complement(read_data[i]);
            localAlignment.align(reverse_complement);
        }
        if ((i + 1) % 100 == 0)
        {
            timer.printExecutionTime();
            std::cout << ", processed " << i + 1 << "/" << read_data.size() << std::endl;
        }
    }

    std::string csv_result = parse_results(result_map, confirmation_count);

    timer.printExecutionTime();

    return csv_result;
}

std::string MutationFinder::parse_results(std::map<int, std::map<std::string, int>> result_map, int confirmation_count)
{
    std::ostringstream result_stream;
    for (int i = 0; i < reference.size(); i++)
    {
        auto it = result_map.find(i);
        if (it != result_map.end())
        {
            auto result = get_max(it->second);
            if (result.first[0] != 'M' && result.second > confirmation_count)
            {
                result_stream << result.first[0] << "," << i << "," << result.first[1] << std::endl;
            }
        }
    }

    return result_stream.str();
}