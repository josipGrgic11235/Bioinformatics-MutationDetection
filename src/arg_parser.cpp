#include <arg_parser.hpp>
#include <algorithm>
#include <iostream>

char *get_arg(char **begin, char **end, const std::string &option)
{
    char **itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

bool arg_present(char **begin, char **end, const std::string &option)
{
    return std::find(begin, end, option) != end;
}

void print_help(char *program_name)
{
    std::cerr << "Usage:"
              << "\n\t-r   [reference genome file path] (file_path)"
              << "\n\t-s   [sequenced reads file path] (file_path) "
              << "\n\t-o   [output file path] (file_path)"
              << "\n\t-k   [kmer length, default=" << K << "] (int)"
              << "\n\t-rd  [region divider, default=" << REGION_DIVIDER << "] (int)"
              << "\n\t-lak [local alignment k value, default=" << LOCAL_ALLIGN_K << "] (int)"
              << "\n\t-cc  [confirmation count, default=" << CONFIRMATION_COUNT << "] (int)"
              << "\n\t-m   [match score, default=" << M << "] (int)"
              << "\n\t-x   [mismatch score, default=" << X << "] (int)"
              << "\n\t-g   [gap score, default=" << G << "] (int)"
              << std::endl;
    exit(-1);
}

Args ArgParser::parse(int argc, char *argv[])
{
    Args args;

    if (arg_present(argv, argv + argc, "-h"))
    {
        print_help(argv[0]);
    }
    if (arg_present(argv, argv + argc, "-k"))
    {
        args.kmer_k = std::atoi(get_arg(argv, argv + argc, "-k"));
    }
    if (arg_present(argv, argv + argc, "-rd"))
    {
        args.region_divider = std::atoi(get_arg(argv, argv + argc, "-rd"));
    }
    if (arg_present(argv, argv + argc, "-lak"))
    {
        args.local_align_k = std::atoi(get_arg(argv, argv + argc, "-lak"));
    }
    if (arg_present(argv, argv + argc, "-cc"))
    {
        args.confirmation_count = std::atoi(get_arg(argv, argv + argc, "-cc"));
    }
    if (arg_present(argv, argv + argc, "-m"))
    {
        args.match_score = std::atoi(get_arg(argv, argv + argc, "-m"));
    }
    if (arg_present(argv, argv + argc, "-x"))
    {
        args.change_score = std::atoi(get_arg(argv, argv + argc, "-x"));
    }
    if (arg_present(argv, argv + argc, "-g"))
    {
        args.gap_score = std::atoi(get_arg(argv, argv + argc, "-g"));
    }

    bool allMandatoryArgsPresent = true;
    if (arg_present(argv, argv + argc, "-r"))
    {
        args.reference_genome_path = get_arg(argv, argv + argc, "-r");
    }
    else
    {
        std::cerr << "Missing mandatory param [-r] (reference genome path)" << std::endl;
        allMandatoryArgsPresent &= false;
    }

    if (arg_present(argv, argv + argc, "-s"))
    {
        args.sequenced_reads_path = get_arg(argv, argv + argc, "-s");
    }
    else
    {
        std::cerr << "Missing mandatory param [-s] (sequenced reads path)" << std::endl;
        allMandatoryArgsPresent &= false;
    }

    if (arg_present(argv, argv + argc, "-o"))
    {
        args.output_file_path = get_arg(argv, argv + argc, "-o");
    }
    else
    {
        std::cerr << "Missing mandatory param [-o] (output file path)" << std::endl;
        allMandatoryArgsPresent &= false;
    }

    if (!allMandatoryArgsPresent)
    {
        exit(-1);
    }
    return args;
}