#include <map>
#include <vector>

#ifndef KMER_INDEX
#define KMER_INDEX

class KmerIndex
{
  std::map<std::string, std::vector<int>> kmer_map;
  int reference_size;

  void generate_kmers(std::string &input, int k);

public:
  KmerIndex(std::string &input, int k);
  std::vector<int> find(std::string &kmer);
  int get_reference_size();
};
#endif
