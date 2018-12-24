#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <utility>
#include <numeric>
#include <ctype.h>
#include <chrono>

#define REGION_DIVIDER 100

void generate_kmers(std::map<std::string, std::vector<int>> &kmer_map, std::string &input, int k)
{
    for (int i = 0; i <= input.length() - k; i++)
    {
        std::string kmer = input.substr(i, k);
        auto it = kmer_map.find(kmer);
        if (it != kmer_map.end())
        {
            it->second.push_back(i);
        }
        else
        {
            std::vector<int> index_list = {i};
            kmer_map.insert(std::pair<std::string, std::vector<int>>(kmer, index_list));
        }
    }
}

std::string read_reference(std::string path)
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

void read_sequence_data(std::vector<std::string> &sequence_list, std::string path)
{
    std::ifstream file(path);

    if (file.is_open())
    {
        std::string line;
        while (getline(file, line))
        {
            getline(file, line);
            sequence_list.push_back(line);
        }
        file.close();
    }
}

std::pair<int, int> map_to_reference(std::map<std::string, std::vector<int>> &reference_kmer_map, std::string &reference, std::string &input, int k)
{
    std::vector<int> region(reference.size() / REGION_DIVIDER + 1, 0);

    for (int i = 0; i <= input.length() - k; i += k)
    {
        std::string kmer = input.substr(i, k);
        auto it = reference_kmer_map.find(kmer);
        if (it != reference_kmer_map.end())
        {
            auto index_list = it->second;
            for (auto i = index_list.begin(); i != index_list.end(); i++)
            {
                region[*i / REGION_DIVIDER]++;
            }
        }
    }

    int sequence_region_length = (input.size() / REGION_DIVIDER);

    int max = 0;
    int region_id = -1;
    for (int i = 0; i <= region.size() - sequence_region_length; i++)
    {
        int region_count = std::accumulate(region.begin() + i, region.begin() + i + sequence_region_length, 0);
        if (region_count > max)
        {
            max = region_count;
            region_id = i;
        }
    }

    return std::pair<int, int>(region_id, region_id + sequence_region_length);
}

char get_complement_base(char c)
{
    switch (std::toupper(c))
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
        throw std::invalid_argument("Unknown base " + c);
    }
}

void append_reverse_complements(std::vector<std::string> &sequence_list)
{
    int sequence_length = sequence_list.size();
    for (int i = 0; i < sequence_length; i++)
    {
        std::string reverse_complement;
        for (int j = sequence_list[i].size() - 1; j >= 0; j--)
        {
            reverse_complement += get_complement_base(sequence_list[i][j]);
        }

        sequence_list.push_back(reverse_complement);
    }
}

int main(int argc, char *argv[])
{
    auto start = std::chrono::high_resolution_clock::now();

    std::string args_reference_genome = argv[1];
    std::string args_sequenced_reads = argv[2];
    int k = atoi(argv[3]);

    std::ios::sync_with_stdio(false);
    std::string reference = read_reference(args_reference_genome);
    std::map<std::string, std::vector<int>> kmer_map;
    generate_kmers(kmer_map, reference, k);

    std::vector<std::string> sequence_list;
    read_sequence_data(sequence_list, args_sequenced_reads);
    append_reverse_complements(sequence_list);

    for (int i = 0; i < sequence_list.size(); i++)
    {
        auto result = map_to_reference(kmer_map, reference, sequence_list[i], k);
        std::cout << i << " " << result.first << "-" << result.second << std::endl;
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::cout << std::endl
              << "Execution time: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(finish - start).count() / 1e6
              << "ms\n";
    return 0;
}