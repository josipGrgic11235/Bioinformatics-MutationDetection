#include <mutation_finder.hpp>
#include <kmer_index.hpp>
#include <read_mapper.hpp>
#include <local_alignment.hpp>
#include <reverse_complement.hpp>
#include <get_rss.hpp>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <fstream>

template <typename KeyType, typename ValueType>
std::pair<KeyType, ValueType> get_max(const std::map<KeyType, ValueType> &x)
{
    using pairtype = std::pair<KeyType, ValueType>;
    return *std::max_element(x.begin(), x.end(), [](const pairtype &p1, const pairtype &p2) {
        return p1.second < p2.second;
    });
}

void MutationFinder::find_mutations(int k, int region_divider, int local_align_k, int match_score, int change_score, int gap_score, int confirmation_count, std::string output_file, ExecutionTimer &timer)
{
    KmerIndex kmerIndex(reference, k);
    std::cout << "Done with generating reference genome k-mers." << std::endl;

    std::map<int, std::map<std::string, int>> result_map;
    ReadMapper readMapper(kmerIndex, region_divider, k);
    LocalAlignment localAlignment(readMapper, reference, result_map, local_align_k, region_divider, match_score, change_score, gap_score);

    std::cout << "Starting read alignment." << std::endl;
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
            std::cout << ", processed "
                      << i + 1 << "/"
                      << read_data.size()
                      << ", peak ram usage "
                      << getPeakRSS() / (1024 * 1024)
                      << " MB"
                      << std::endl;
        }
    }
    std::cout << "Done with read alignment." << std::endl;

    parse_results(result_map, confirmation_count, output_file);

    std::cout << "Execution completed." << std::endl;
    timer.printExecutionTime();
    std::cout << ", peak ram usage "
              << getPeakRSS() / (1024 * 1024)
              << " MB"
              << std::endl;
}

void MutationFinder::parse_results(std::map<int, std::map<std::string, int>> result_map, int confirmation_count, std::string output_path)
{
    std::ofstream output_file;
    output_file.open(output_path);

    for (int i = 0; i < reference.size(); i++)
    {
        auto it = result_map.find(i);
        if (it != result_map.end())
        {
            auto result = get_max(it->second);
            if (result.first[0] != 'M' && result.second > confirmation_count)
            {
                output_file << result.first[0] << "," << i << "," << result.first[1] << std::endl;
            }
        }
    }
    output_file.close();
}