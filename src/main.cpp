#include <reference_reader.hpp>
#include <sequence_reader.hpp>
#include <mutation_finder.hpp>
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

#define REGION_DIVIDER 100
#define INPUT_DIVIDER 4000
#define LOCAL_ALLIGN_K 50
#define CONFIRMATION_COUNT 3
#define M 5
#define X -4
#define G -7

int main(int argc, char *argv[])
{
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
    std::vector<std::string> sequence_list = SequenceReader::read_sequence_data(args_sequenced_reads);

    MutationFinder mutationFinder(reference, sequence_list);

    std::ofstream output_file;
    output_file.open(args_output_file_name);
    output_file << mutationFinder.find_mutations(k, REGION_DIVIDER, LOCAL_ALLIGN_K, M, X, G, CONFIRMATION_COUNT);
    output_file.close();

    return 0;
}

/*

./main.exe C:\\Users\\leon\\Documents\\Bioinformatika\\Bioinformatics-MutationDetection\\train_data\\lambda.fasta C:\\Users\\leon\\Documents\\Bioinformatika\\Bioinformatics-MutationDetection\\train_data\\lambda_simulated_reads.fasta train_data\\lambda_result.csv 8
./a.exe C:\\Users\\leon\\Documents\\Bioinformatika\\Bioinformatics-MutationDetection\\train_data\\ecoli.fasta C:\\Users\\leon\\Documents\\Bioinformatika\\Bioinformatics-MutationDetection\\train_data\\ecoli_simulated_reads.fasta train_data\\ecoli_result.csv 8

python train_data/jaccard.py -b train_data/lambda_mutated.csv -a train_data/lambda_result.csv
python train_data/jaccard.py -b train_data/ecoli_mutated.csv -a train_data/ecoli_result.csv

 g++ -I src/ src/*.cpp -o main.exe

*/