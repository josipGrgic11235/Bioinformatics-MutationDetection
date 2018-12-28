#include <arg_parser.hpp>
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

/* 
g++ -I src/ src/*.cpp -o main.exe
python train_data/jaccard.py -b train_data/lambda_mutated.csv -a train_data/lambda_result.csv
python train_data/jaccard.py -b train_data/ecoli_mutated.csv -a train_data/ecoli_result.csv

*/
int main(int argc, char *argv[])
{
    auto args = ArgParser::parse(argc, argv);

    std::ios::sync_with_stdio(false);

    std::string reference = ReferenceReader::read_reference(args.reference_genome_path);
    std::vector<std::string> sequence_list = SequenceReader::read_sequence_data(args.sequenced_reads_path);

    MutationFinder mutationFinder(reference, sequence_list);

    std::ofstream output_file;
    output_file.open(args.output_file_path);
    output_file << mutationFinder.find_mutations(args.kmer_k, args.region_divider, args.local_align_k, args.match_score, args.change_score, args.gap_score, args.confirmation_count);
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