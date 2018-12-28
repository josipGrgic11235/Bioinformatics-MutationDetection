#include <string>
#include <vector>
#include <map>

class MutationFinder
{
  std::string &reference;
  std::vector<std::string> &read_data;
  std::string parse_results(std::map<int, std::map<std::string, int>> result_map, int confirmation_count);

public:
  MutationFinder(std::string &reference, std::vector<std::string> &read_data) : reference(reference), read_data(read_data) {}
  std::string find_mutations(int k, int region_divider, int local_align_k, int match_score, int change_score, int gap_score, int confirmation_count);
};