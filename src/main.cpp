#include <arg_parser.hpp>
#include <reference_reader.hpp>
#include <sequence_reader.hpp>
#include <mutation_finder.hpp>
#include <iostream>
#include <fstream>
#include <vector>

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