// Created by Leon

#include <string>
#include <vector>

class SequenceReader
{
public:
  /**
   * Reads sequence from given path.
   * returns: vector<string> = list of sequenced reads
   **/
  static std::vector<std::string> read_sequence_data(std::string path);
};