// Created by Leon

#include <kmer_index.hpp>
#include <iostream>

KmerIndex::KmerIndex(std::string &input, int k)
{
    reference_size = input.size();
    generate_kmers(input, k);
}

void KmerIndex::generate_kmers(std::string &input, int k)
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

std::vector<int> KmerIndex::find(std::string &kmer)
{
    auto it = kmer_map.find(kmer);
    if (it != kmer_map.end())
    {
        return it->second;
    }
    return {};
}

int KmerIndex::get_reference_size()
{
    return reference_size;
}