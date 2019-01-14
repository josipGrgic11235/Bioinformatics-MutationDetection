// Created by Leon

#include <map>
#include <vector>

#ifndef KMER_INDEX
#define KMER_INDEX

/**
 * Generates k-mer index for given input.
 **/
class KmerIndex
{
  std::map<std::string, std::vector<int>> kmer_map;
  int reference_size;

  void generate_kmers(std::string &input, int k);

public:
  /**
   * Constructor
   * params: input = input string
   * params: k = k-mer length
   **/
  KmerIndex(std::string &input, int k);

  /**
   * Find index of given k-mer.
   * params: kmer = k-mer that we want to find in input string
   * returns: list of indexes where given k-mer is present in the input string
   **/
  std::vector<int> find(std::string &kmer);

  /**
   * returns: number of characters in input string
   **/
  int get_reference_size();
};
#endif
