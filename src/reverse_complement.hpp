// Created by Leon

#include <string>

class ReverseComplement
{
public:
  /**
   * Returns complement of the given base.
   * param: c = base for which we want to get it's complement.
   * A-T, G-C
   * returns: string = complement base
   **/
  static char get_complement_base(char c);

  /**
   * Generates reverse complement for input string.
   * returns: string = reverse complement
   **/
  static std::string get_reverse_complement(std::string input);
};