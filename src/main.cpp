#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <utility>
#include <numeric>
#include <ctype.h>
#include <chrono>
#include <algorithm>
#include <tuple>
#include <cmath>
#include <kmer_index.hpp>
#include <reference_reader.hpp>
#include <sequence_reader.hpp>
#include <read_mapper.hpp>
#include <reverse_complement.hpp>
#include <local_alignment.hpp>

#define REGION_DIVIDER 100
#define INPUT_DIVIDER 4000
#define LOCAL_ALLIGN_K 50
#define CONFIRMATION_COUNT 3
#define M 5
#define X -4
#define G -7

template <typename KeyType, typename ValueType>
std::pair<KeyType, ValueType> get_max(const std::map<KeyType, ValueType> &x)
{
    using pairtype = std::pair<KeyType, ValueType>;
    return *std::max_element(x.begin(), x.end(), [](const pairtype &p1, const pairtype &p2) {
        return p1.second < p2.second;
    });
}

int main(int argc, char *argv[])
{
    auto start = std::chrono::high_resolution_clock::now();

    /*std::string args_reference_genome = argv[1];
    std::string args_sequenced_reads = argv[2];
    std::string args_output_file_name = argv[3];
    int k = atoi(argv[4]);*/

    std::string args_reference_genome = "C:\\Users\\leon\\Documents\\Bioinformatika\\Bioinformatics-MutationDetection\\train_data\\lambda.fasta";
    std::string args_sequenced_reads = "C:\\Users\\leon\\Documents\\Bioinformatika\\Bioinformatics-MutationDetection\\train_data\\lambda_simulated_reads.fasta";
    std::string args_output_file_name = "train_data\\lambda_result.csv";
    int k = 8;

    std::ios::sync_with_stdio(false);

    std::string reference = ReferenceReader::read_reference(args_reference_genome);
    KmerIndex kmerIndex(reference, k);
    std::vector<std::string> sequence_list = SequenceReader::read_sequence_data(args_sequenced_reads);

    std::map<int, std::map<std::string, int>> result_map;
    ReadMapper readMapper(kmerIndex, REGION_DIVIDER, k);
    LocalAlignment localAlignment(readMapper, reference, result_map, LOCAL_ALLIGN_K, REGION_DIVIDER, M, X, G);

    for (int i = 0; i < sequence_list.size(); i++)
    {
        if (!localAlignment.align(sequence_list[i]))
        {
            std::string reverse_complement = ReverseComplement::get_reverse_complement(sequence_list[i]);
            localAlignment.align(reverse_complement);
        }
        if ((i + 1) % 100 == 0)
        {
            auto finish = std::chrono::high_resolution_clock::now();
            std::cout << std::endl
                      << "Execution time: "
                      << std::chrono::duration_cast<std::chrono::nanoseconds>(finish - start).count() / 1e6
                      << "ms, ";
            std::cout << "processed " << i + 1 << "/" << sequence_list.size() << std::endl;
        }
    }

    std::ofstream output_file;
    output_file.open(args_output_file_name);
    for (int i = 0; i < reference.size(); i++)
    {
        auto it = result_map.find(i);
        if (it != result_map.end())
        {
            auto result = get_max(it->second);
            if (result.first[0] != 'M' && result.second > CONFIRMATION_COUNT)
            {
                output_file << result.first[0] << "," << i << "," << result.first[1] << std::endl;
            }
        }
    }
    output_file.close();

    auto finish = std::chrono::high_resolution_clock::now();
    std::cout << std::endl
              << "Execution time: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(finish - start).count() / 1e6
              << "ms\n";

    /*std::map<int, std::map<std::string, int>> result_map;
    std::string reference = "TATATGCGGCGTTT";
    std::string input = "GGTATGCTGGCGCTA";

    //apply_local_allign(reference, input);
    apply_local_allign_optimized(reference, input, 3, result_map, 0);*/

    return 0;
}

/*

./main.exe C:\\Users\\leon\\Documents\\Bioinformatika\\Bioinformatics-MutationDetection\\train_data\\lambda.fasta C:\\Users\\leon\\Documents\\Bioinformatika\\Bioinformatics-MutationDetection\\train_data\\lambda_simulated_reads.fasta train_data\\lambda_result.csv 8
./a.exe C:\\Users\\leon\\Documents\\Bioinformatika\\Bioinformatics-MutationDetection\\train_data\\ecoli.fasta C:\\Users\\leon\\Documents\\Bioinformatika\\Bioinformatics-MutationDetection\\train_data\\ecoli_simulated_reads.fasta train_data\\ecoli_result.csv 8

python train_data/jaccard.py -b train_data/lambda_mutated.csv -a train_data/lambda_result.csv
python train_data/jaccard.py -b train_data/ecoli_mutated.csv -a train_data/ecoli_result.csv

 g++ -I src/ src/*.cpp -o main.exe

*/