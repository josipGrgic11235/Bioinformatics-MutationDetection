#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <numeric>

#define KMER_SIZE 8
#define REGION_DIVIDER 100

void generate_kmers(std::map<std::string, std::vector<int>> &kmer_map, std::string input, int k)
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

std::string read_reference()
{
    std::ifstream file("C:\\Users\\leon\\Documents\\Bioinformatika\\Bioinformatics-MutationDetection\\train_data\\lambda.fasta");
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

void read_sequence_data(std::vector<std::string> &sequence_list)
{
    std::ifstream file("C:\\Users\\leon\\Documents\\Bioinformatika\\Bioinformatics-MutationDetection\\train_data\\lambda_simulated_reads.fasta");

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

int main()
{
    std::ios::sync_with_stdio(false);
    std::string reference = read_reference();
    std::map<std::string, std::vector<int>> kmer_map;
    generate_kmers(kmer_map, reference, KMER_SIZE);

    std::vector<std::string> sequence_list;
    read_sequence_data(sequence_list);

    std::map<std::string, std::vector<int>> sequence_kmer_map;

    std::vector<int> region(reference.size() / REGION_DIVIDER + 1, 0);

    auto interest_list = sequence_list[22];
    for (int i = 0; i <= interest_list.length() - KMER_SIZE; i++)
    {
        std::string kmer = interest_list.substr(i, KMER_SIZE);
        auto it = kmer_map.find(kmer);
        if (it != kmer_map.end())
        {
            auto index_list = it->second;
            for (auto i = index_list.begin(); i != index_list.end(); i++)
            {
                region[*i / REGION_DIVIDER]++;
            }
        }
    }

    int sequence_region_length = (interest_list.size() / REGION_DIVIDER + 1);

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
    std::cout << region_id << "-" << region_id + sequence_region_length << " " << max << std::endl;

    return 0;
}