#include <arg_parser.hpp>
#include <reference_reader.hpp>
#include <sequence_reader.hpp>
#include <mutation_finder.hpp>
#include <execution_timer.hpp>
#include <iostream>
#include <vector>

int main(int argc, char *argv[])
{
    ExecutionTimer timer;
    timer.start();
    auto args = ArgParser::parse(argc, argv);

    std::ios::sync_with_stdio(false);

    std::cout << "Reading input files." << std::endl;
    std::string reference = ReferenceReader::read_reference(args.reference_genome_path);
    std::vector<std::string> sequence_list = SequenceReader::read_sequence_data(args.sequenced_reads_path);
    std::cout << "Done with reading input files." << std::endl;
    timer.printExecutionTime();
    std::cout << std::endl;

    MutationFinder mutationFinder(reference, sequence_list);
    mutationFinder.find_mutations(args.kmer_k, args.region_divider, args.local_align_k, args.match_score, args.change_score, args.gap_score, args.confirmation_count, args.output_file_path, timer);

    return 0;
}