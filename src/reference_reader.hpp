// Created by Leon

#include <string>

class ReferenceReader
{
public:
  /**
   * Reads the reference genome from the file from given path and stores it into a string.
   **/
  static std::string read_reference(std::string path);
};