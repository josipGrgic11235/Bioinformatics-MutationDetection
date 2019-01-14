#include <read_mapper.hpp>
#include <string>
#include <map>

enum Action
{
    Match,
    Insertion,
    Deletion,
    None
};

struct Score
{
    int score;
    Action action;
    Score() : score(0), action(None) {}
    Score(int score, Action action) : score(score), action(action) {}
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

    bool apply_local_allign(std::string &reference_region, std::string &input, int reference_offset);
    std::string get_substr(std::string &input, int start, int length);
    int get_substr_length(int max_length, int start, int length);
    bool backtrack(Score **matrix, int max_i, int max_j, int columns, std::string &input, std::string &reference_region, int reference_offset);
    int get_array_size_at_row(int row, int max_columns);
    void print_local_allign_matrix(Score **matrix, int rows, int columns, std::string &reference_region, std::string &input);

    static int get_distance(std::string &s1, std::string &s2);

  public:
    LocalAlignment(ReadMapper &readMapper, std::string &reference, std::map<int, std::map<std::string, int>> &result_map, int local_align_k, int region_divider, int match_score, int change_score, int gap_score);

    bool align(std::string &read_input);
};
