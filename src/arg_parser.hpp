// Created by Leon

#include <string>

#define REGION_DIVIDER 100
#define LOCAL_ALLIGN_K 50
#define CONFIRMATION_COUNT 3
#define M 5
#define X -4
#define G -7
#define K 8

/**
 * Structure which holds arguments that are passed to the program.
 * Initial values are assigned as following:
 * region_divider = 100
 * local_align_k = 50
 * confirmation_count = 3
 * match_score = 5
 * change_score = -4
 * gap_score = -7
 * kmer_k = 8
**/
struct Args
{
    int kmer_k = K;
    int region_divider = REGION_DIVIDER;
    int local_align_k = LOCAL_ALLIGN_K;
    int confirmation_count = CONFIRMATION_COUNT;
    int match_score = M;
    int change_score = X;
    int gap_score = G;

    std::string reference_genome_path;
    std::string sequenced_reads_path;
    std::string output_file_path;
};

class ArgParser
{
  public:
    /**
     * Method that parses given arguments and returns Args structure.
     * params: argc = argument count
     * params: argv = arguments
     * returns: Args structure filled with values from arguments
    **/
    static Args parse(int argc, char *argv[]);
};