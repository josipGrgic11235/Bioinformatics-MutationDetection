// Created by Leon

#include <execution_timer.hpp>
#include <string>
#include <vector>
#include <map>

/**
 * Maps sequenced reads to the reference genome, aligns them and finds mutations. 
 **/
class MutationFinder
{
  std::string &reference;
  std::vector<std::string> &read_data;
  void parse_results(std::map<int, std::map<std::string, int>> result_map, int confirmation_count, std::string output_path);

public:
  /**
   * Constructor
   * param: reference = reference genome
   * param: read_data = sequenced reads
   **/
  MutationFinder(std::string &reference, std::vector<std::string> &read_data) : reference(reference), read_data(read_data) {}

  /**
   * Finds mutations and writes them in csv format.
   * param: k = kmer length
   * param: region_divider = reference genome divider used in mapping phase
   * param: local_align_k = boundary length during local alignment
   * param: match_score = score for matching bases
   * param: change_score = score for supstitution mutation
   * param: gap_score = score for indel (insertion and deletion)
   * param: confirmation_conut = number of required confirmations
   * param: output_path = result file path
   **/
  void find_mutations(int k, int region_divider, int local_align_k, int match_score, int change_score, int gap_score, int confirmation_count, std::string output_path, ExecutionTimer &executionTimer);
};