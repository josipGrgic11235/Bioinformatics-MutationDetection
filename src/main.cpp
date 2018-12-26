#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <utility>
#include <numeric>
#include <ctype.h>
#include <chrono>
#include <algorithm>
#include <tuple>
#include <cmath>

#define REGION_DIVIDER 100
#define INPUT_DIVIDER 100
#define M 5
#define X -4
#define G -5

void generate_kmers(std::map<std::string, std::vector<int>> &kmer_map, std::string &input, int k)
{
    for (int i = 0; i <= input.length() - k; i++)
    {
        std::string kmer = input.substr(i, k);
        auto it = kmer_map.find(kmer);
        if (it != kmer_map.end())
        {
            it->second.push_back(i);
        }
        else
        {
            std::vector<int> index_list = {i};
            kmer_map.insert(std::pair<std::string, std::vector<int>>(kmer, index_list));
        }
    }
}

std::string read_reference(std::string path)
{
    std::ifstream file(path);
    std::string reference;

    if (file.is_open())
    {
        std::string line;
        getline(file, line);
        while (getline(file, line))
        {
            reference += line.c_str();
        }
        file.close();
    }

    return reference;
}

void read_sequence_data(std::vector<std::string> &sequence_list, std::string path)
{
    std::ifstream file(path);

    if (file.is_open())
    {
        std::string line;
        while (getline(file, line))
        {
            getline(file, line);
            sequence_list.push_back(line);
        }
        file.close();
    }
}

std::pair<int, int> map_to_reference(std::map<std::string, std::vector<int>> &reference_kmer_map, std::string &reference, std::string &input, int k)
{
    std::vector<int> region(reference.size() / REGION_DIVIDER + 1, 0);

    for (int i = 0; i <= input.size() - k; i += 2)
    {
        std::string kmer = input.substr(i, k);
        auto it = reference_kmer_map.find(kmer);
        if (it != reference_kmer_map.end())
        {
            auto index_list = it->second;
            for (auto i = index_list.begin(); i != index_list.end(); i++)
            {
                region[*i / REGION_DIVIDER]++;
            }
        }
    }

    int sequence_region_length = (input.size() / REGION_DIVIDER) + 1;

    int max = 0;
    int region_id = 0;
    for (int i = 0; i <= region.size() - sequence_region_length; i++)
    {
        int region_count = std::accumulate(region.begin() + i, region.begin() + i + sequence_region_length, 0);
        if (region_count > max)
        {
            max = region_count;
            region_id = i;
        }
    }

    return std::pair<int, int>(region_id, region_id + sequence_region_length);
}

char get_complement_base(char c)
{
    switch (c)
    {
    case 'T':
        return 'A';
    case 'A':
        return 'T';
    case 'G':
        return 'C';
    case 'C':
        return 'G';
    default:
        throw std::invalid_argument("Unknown base " + c);
    }
}

std::string get_reverse_complement(std::string input)
{
    std::string reverse_complement;
    for (int i = input.size() - 1; i >= 0; i--)
    {
        reverse_complement += get_complement_base(input[i]);
    }
    return reverse_complement;
}

int get_base_index(char c)
{
    switch (c)
    {
    case 'A':
        return 0;
    case 'C':
        return 1;
    case 'G':
        return 2;
    case 'T':
        return 3;
    case '-':
        return 4;
    default:
        throw std::invalid_argument("Unknown base " + c);
    }
}

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

bool compare_scores(const Score &left, const Score &right)
{
    return left.score < right.score;
}

Score max(Score &a, Score &b, Score &c, Score &d)
{
    return std::max(std::max(a, b, compare_scores), std::max(c, d, compare_scores), compare_scores);
}

int get_distance(std::string s1, std::string s2)
{
    int distance = 0;
    for (int i = 0; i < s1.size(); i++)
    {
        if (s1[i] != s2[i])
        {
            distance++;
        }
    }
    return distance;
}

void insert_or_increment(std::map<int, std::map<std::string, int>> &result_map, std::string change_id, int index)
{
    result_map[index][change_id]++;
}

bool backtrack(Score **matrix, int max_i, int max_j, int columns, std::string input, std::string reference, std::map<int, std::map<std::string, int>> &result_map, int reference_offset)
{
    int i = max_i;
    int j = max_j;

    std::string input_result;
    std::string reference_result;
    while (i > 0 && j > 0)
    {
        Score score = matrix[i][j];
        if (score.action == Insertion)
        {
            j--;
            input_result += '-';
            reference_result += reference[j];
        }
        else if (score.action == Deletion)
        {
            i--;
            input_result += input[i];
            reference_result += '-';
        }
        else if (score.action == Match)
        {
            i--;
            j--;
            input_result += input[i];
            reference_result += reference[j];
        }
        else
        {
            break;
        }
    }
    //std::cout << i << " " << j << std::endl;
    std::reverse(input_result.begin(), input_result.end());
    std::reverse(reference_result.begin(), reference_result.end());

    //std::cout << reference_result << std::endl;
    //std::cout << input_result << std::endl;

    double similarity = (1 - (float)get_distance(reference_result, input_result) / reference_result.size());
    if (similarity < 0.8)
    {
        return false;
    }

    int deletion_count = 0;
    int insertion_count = 0;
    // TODO deletion_count?????
    int offset = j;
    for (int i = 0; i < reference_result.size(); i++)
    {
        int corrected_index = offset + i - insertion_count + reference_offset;
        if (reference_result[i] == '-')
        {
            //std::cout << "Insertion " << input_result[i] << " at index: " << corrected_index << std::endl;
            result_map[corrected_index][std::string("I") + input_result[i]]++;
            insertion_count++;
        }
        else if (input_result[i] == '-')
        {
            //std::cout << "Deletion " << reference_result[i] << " at index: " << corrected_index << std::endl;
            result_map[corrected_index][std::string("D") + input_result[i]]++;
            deletion_count++;
        }
        else if (input_result[i] != reference_result[i])
        {
            //std::cout << "Change " << input_result[i] << " at index: " << corrected_index << std::endl;
            result_map[corrected_index][std::string("X") + input_result[i]]++;
        }
        else
        {
            //std::cout << "Match " << input_result[i] << " at index: " << corrected_index << std::endl;
            result_map[corrected_index][std::string("M") + input_result[i]]++;
        }
    }

    return true;
}

bool apply_local_allign(std::string const &reference, std::string const &input, std::map<int, std::map<std::string, int>> &result_map, int reference_offset)
{
    int scoring_matrix[5][5] = {
        M, X, X, X, G,
        X, M, X, X, G,
        X, X, M, X, G,
        X, X, X, M, G,
        G, G, G, G, 0};

    int rows = input.length() + 1;
    int columns = reference.length() + 1;
    int cell_count = rows * columns;

    Score **matrix = new Score *[rows];

    for (int i = 0; i < rows; i++)
    {
        matrix[i] = new Score[columns];
    }

    int max_value = 0;
    int max_index_i = 0;
    int max_index_j = 0;

    Score *none = new Score(0, None);

    int minus_index = get_base_index('-');

    Score *horizontal_distance = new Score(0, Insertion);
    Score *vertical_distance = new Score(0, Deletion);
    Score *diagonal_distance = new Score(0, Match);

    for (int i = 1; i < rows; i++)
    {
        char current_input_character = input[i - 1];
        int current_input_index = get_base_index(current_input_character);

        for (int j = 1; j < columns; j++)
        {
            char current_reference_character = reference[j - 1];
            int current_reference_index = get_base_index(current_reference_character);

            horizontal_distance->score = matrix[i][j - 1].score + scoring_matrix[current_input_index][minus_index];
            vertical_distance->score = matrix[i - 1][j].score + scoring_matrix[minus_index][current_reference_index];
            diagonal_distance->score = matrix[i - 1][j - 1].score + scoring_matrix[current_input_index][current_reference_index];

            matrix[i][j] = max(*horizontal_distance, *vertical_distance, *diagonal_distance, *none);

            if (matrix[i][j].score > max_value)
            {
                max_index_i = i;
                max_index_j = j;
                max_value = matrix[i][j].score;
            }
        }
    }

    bool backtrack_result = backtrack(matrix, max_index_i, max_index_j, columns, input, reference, result_map, reference_offset);

    for (int i = 0; i < rows; i++)
    {
        delete matrix[i];
    }
    delete[] matrix;

    return backtrack_result;
}

template <typename KeyType, typename ValueType>
std::pair<KeyType, ValueType> get_max(const std::map<KeyType, ValueType> &x)
{
    using pairtype = std::pair<KeyType, ValueType>;
    return *std::max_element(x.begin(), x.end(), [](const pairtype &p1, const pairtype &p2) {
        return p1.second < p2.second;
    });
}

int get_substr_length(int max_length, int start, int length)
{
    return start + length > max_length ? max_length - start : length;
}

std::string get_substr(std::string input, int start, int length)
{
    int input_start = start;
    int input_length = get_substr_length(input.size(), input_start, length);
    return input.substr(input_start, input_length);
}

bool align(std::map<std::string, std::vector<int>> &kmer_map, std::map<int, std::map<std::string, int>> &result_map, std::string &reference, std::string &read_input, int k)
{
    bool align_result = false;
    for (int j = 0; j < std::ceil(read_input.size() / (float)INPUT_DIVIDER); j++)
    {
        std::string input = get_substr(read_input, j * INPUT_DIVIDER, INPUT_DIVIDER);

        // TODO prevent index out of bounds
        if (input.size() < k)
            continue;

        auto result = map_to_reference(kmer_map, reference, input, k);

        // Ain't nobody got time for that
        /*int region_start = std::max(0, (result.first - 2) * REGION_DIVIDER);
        std::string reference_substring = get_substr(reference, region_start, (result.second - result.first + 2) * REGION_DIVIDER);*/
        int region_start = result.first * REGION_DIVIDER;
        std::string reference_substring = get_substr(reference, region_start, (result.second - result.first) * REGION_DIVIDER);

        align_result |= apply_local_allign(reference_substring, input, result_map, region_start);
    }

    return align_result;
}
int main(int argc, char *argv[])
{
    auto start = std::chrono::high_resolution_clock::now();

    std::string args_reference_genome = argv[1];
    std::string args_sequenced_reads = argv[2];
    std::string args_output_file_name = argv[3];
    int k = atoi(argv[4]);

    std::ios::sync_with_stdio(false);
    std::string reference = read_reference(args_reference_genome);

    std::map<std::string, std::vector<int>> kmer_map;
    generate_kmers(kmer_map, reference, k);

    std::vector<std::string> sequence_list;
    read_sequence_data(sequence_list, args_sequenced_reads);

    std::map<int, std::map<std::string, int>> result_map;

    for (int i = 0; i < sequence_list.size(); i++)
    {
        if (!align(kmer_map, result_map, reference, sequence_list[i], k))
        {
            std::string reverse_complement = get_reverse_complement(sequence_list[i]);
            align(kmer_map, result_map, reference, reverse_complement, k);
        }
        std::cout << "Processed " << i + 1 << "/" << sequence_list.size() << std::endl;
    }

    std::ofstream output_file;
    output_file.open(args_output_file_name);
    for (int i = 0; i < reference.size(); i++)
    {
        auto it = result_map.find(i);
        if (it != result_map.end())
        {
            auto result = get_max(it->second);
            if (result.first[0] != 'M' && result.second > 2)
            {
                output_file << result.first[0] << "," << i << "," << result.first[1] << std::endl;
            }
        }
    }
    output_file.close();

    auto finish = std::chrono::high_resolution_clock::now();
    std::cout << std::endl
              << "Execution time: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(finish - start).count() / 1e6
              << "ms\n";
    return 0;
}

/*

./a.exe C:\\Users\\leon\\Documents\\Bioinformatika\\Bioinformatics-MutationDetection\\train_data\\lambda.fasta C:\\Users\\leon\\Documents\\Bioinformatika\\Bioinformatics-MutationDetection\\train_data\\lambda_simulated_reads.fasta train_data\\lambda_result.csv 8

python train_data/jaccard.py -b train_data/lambda_mutated.csv -a train_data/lambda_result.csv
*/