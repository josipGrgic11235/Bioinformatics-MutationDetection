#include <read_mapper.hpp>
#include <string>
#include <map>

/*
* The possible actions for the local align algorithm
*/
enum Action
{
    Match,
    Insertion,
    Deletion,
    None
};

/*
* The cell types for the generated matrix as part of the local align algorithm
*/
struct Score
{
    int score;
    Action action;
    Score() : score(0), action(None) {}
    Score(int score, Action action) : score(score), action(action) {}
};

/*
* Relevant data about the cell with the highest value inside the generated matrix as part of the local align algorithm
*/
struct Max_Cell
{
    int i;
    int j;
    int value;
};

class LocalAlignment
{
    unsigned const char A = 'A';
    unsigned const char C = 'C';
    unsigned const char G = 'G';
    unsigned const char T = 'T';
    unsigned const char X = 'X';
    unsigned const char INDEL = '-';

    ReadMapper &readMapper;
    std::string &reference;
    std::map<int, std::map<std::string, int>> &result_map;
    int region_divider;
    int local_align_k;
    int match_score;
    int change_score;
    int gap_score;
    int scoring_matrix['T' + 1]['T' + 1];
    Score *none = new Score(0, None);

    /**
    * The local align algorithm (Smith-Waterman).
    * Optimized using the bounded dynamic programming; using the -lak parameter passed from the command line (local_align_k class member) 
    * params: reference_region = the suspected region within the reference genome the input is matched against
    * params: input = the seqeunce which is to be matched with the reference genome 
    * params: reference_offset = the start index of the reference_region parameter inside the reference genome
    * returns: True if the similarity between the most similar matched substrings returned from the local align function is greater than 80%, false otherwise
    **/
    bool apply_local_allign(std::string &reference_region, std::string &input, int reference_offset);

    /*
    * params: input = the original string
    * params: start = the start index of the substring 
    * params: length = the desired length of the substring
    * returns: the string with the specified length
    * remarks: The desired substring length can be unfeasible. Please check the LocalAlignment::get_substr_length specification for more info
    */
    std::string get_substr(std::string &input, int start, int length);

    /*
    * Calculates the length of the substring
    * params: max_lenght = the maximum length of the substirng (= the length of the original string)
    * params: start = the start index of the substring 
    * params: length = the desired lenght of the substring
    * returns: the length of the substring, if the desired length is out of bounds in respect to the original string and the start index, only the length from the start index to the end of the string is returned 
    */
    int get_substr_length(int max_length, int start, int length);

    /*
    * The backtrack procedure - used for indentifying and constructing the most simmilar substrings in the reference genome region and the sequence
    * params: matrix = the constructed matrix from the apply_local_allign method
    * params: max_cell = the i (horizontal) and j (vertical) coordinates, and the value of the cell with the greatest value inisde the local align matrix
    * returns: True if the similarity between the most similar matched substrings returned from the local align function is greater than 80%, false otherwise
    */
    bool backtrack(Score **matrix, Max_Cell max_cell, int columns, std::string &input, std::string &reference_region, int reference_offset);

    /*
    * Calculates the cell count at the specified row based on the local_align_k value
    * params: row = the row index
    * params: max_columns = the column count of the local align matrix
    * returns: the number of needed cells at the specified row
    */
    int get_array_size_at_row(int row, int max_columns);

    /*
    * Calculates the Hamming distance between two strings of the same length
    * params: s1 = the first string
    * params: s2 = the second string
    * returns: the number of different characters on the same index inside the passed strings
    */
    static int get_distance(std::string &s1, std::string &s2);

  public:
    // constructor method
    LocalAlignment(ReadMapper &readMapper, std::string &reference, std::map<int, std::map<std::string, int>> &result_map, int local_align_k, int region_divider, int match_score, int change_score, int gap_score);

    /**
    * Main function for the align procedure
    * params: read_input = the sequence we wish to allign with the reference genome passed via the ctor
    * returns: True if the similarity between the most similar matched substrings returned from the local align function is greater than 80%, false otherwise
    */
    bool align(std::string &read_input);
};
