#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <utility>
#include <numeric>
#include <ctype.h>
#include <chrono>
#include <algorithm>

#define REGION_DIVIDER 100

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

    for (int i = 0; i <= input.length() - k; i += k)
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

    int sequence_region_length = (input.size() / REGION_DIVIDER);

    int max = 0;
    int region_id = -1;
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
    switch (std::toupper(c))
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

void append_reverse_complements(std::vector<std::string> &sequence_list)
{
    int sequence_length = sequence_list.size();
    for (int i = 0; i < sequence_length; i++)
    {
        std::string reverse_complement;
        for (int j = sequence_list[i].size() - 1; j >= 0; j--)
        {
            reverse_complement += get_complement_base(sequence_list[i][j]);
        }

        sequence_list.push_back(reverse_complement);
    }
}

int max(int a, int b, int c, int d)
{
    return std::max(std::max(a, b), std::max(c, d));
}

int get_index(int i, int j, int columns)
{
    return i * columns + j;
}

int get_base_index(char c)
{
    switch (std::toupper(c))
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

void apply_local_allign(std::string const &reference, std::string const &input)
{
    int scoring_matrix[5][5] = {
        2, -4, -4, -4, -6,
        -4, 2, -4, -4, -6,
        -4, -4, 2, -4, -6,
        -4, -4, -4, 2, -6,
        -6, -6, -6, -6, 0};

    int rows = input.length() + 1;
    int columns = reference.length() + 1;
    int cell_count = rows * columns;

    int *matrix = new int[cell_count];

    for (int i = 0; i < cell_count; i++)
    {
        matrix[i] = 0;
    }

    int max_value = 0;
    int max_index_i = 0;
    int max_index_j = 0;

    for (int i = 1; i < rows; i++)
    {
        char current_input_character = input[i - 1];
        int current_input_index = get_base_index(current_input_character);

        for (int j = 1; j < columns; j++)
        {
            char current_reference_character = reference[j - 1];
            int current_reference_index = get_base_index(current_reference_character);

            int horizontal_distance = matrix[get_index(i, j - 1, columns)] + scoring_matrix[current_input_index][get_base_index('-')];
            int vertical_distance = matrix[get_index(i - 1, j, columns)] + scoring_matrix[get_base_index('-')][current_reference_index];
            int diagonal_distance = matrix[get_index(i - 1, j - 1, columns)] + scoring_matrix[current_input_index][current_reference_index];

            matrix[get_index(i, j, columns)] = max(horizontal_distance, vertical_distance, diagonal_distance, 0);

            if (matrix[get_index(i, j, columns)] > max_value)
            {
                max_index_i = i;
                max_index_j = j;
            }
        }
    }
}

int get_array_size_at_row(int row, int max_columns, int k)
{
    return std::min(max_columns - 1, row + k) - std::max(0, row - k) + 1;
}

int get_first_index(int row, int k)
{
    return std::max(0, row - k);
}

int get_last_index(int row, int max_columns, int k)
{
    return std::min(max_columns - 1, row + k);
}

void apply_local_allign_optimized(std::string const &reference, std::string const &input, int k)
{
    int scoring_matrix[5][5] = {
        {2, -4, -4, -4, -6},
        {-4, 2, -4, -4, -6},
        {-4, -4, 2, -4, -6},
        {-4, -4, -4, 2, -6},
        {-6, -6, -6, -6, 0}};

    // Matrix init
    int rows = input.length() + 1;
    int columns = reference.length() + 1;

    int **matrix = new int *[rows];
    for (int i = 0; i < rows; i++)
    {
        matrix[i] = new int[get_array_size_at_row(i, columns, k)];
        std::cout << "Array size at i = " << i << ": " << get_array_size_at_row(i, columns, k) << std::endl;
    }
    std::cout << std::endl;

    for (int i = 0; i <= k; i++)
    {
        matrix[0][i] = 0;
        matrix[i][0] = 0;
    }

    int max_value = 0;
    int max_index_i = 0;
    int max_index_j = 0;

    for (int i = 1; i < rows; i++)
    {
        char current_input_character = input[i - 1];
        int current_input_index = get_base_index(current_input_character);

        int start_reference_index = std::max(0, i - k - 1);
        int offset = (i <= k) ? 0 : 1;
        int array_size = get_array_size_at_row(i, columns, k);
        int start_array_index = offset == 0 ? 1 : 0;

        for (int j = start_array_index; j < array_size; j++)
        {
            char current_reference_character = reference[start_reference_index + j - start_array_index];
            int current_reference_index = get_base_index(current_reference_character);

            int horizontal_distance = j == 0 ? 0 : (matrix[i][j - 1] + scoring_matrix[current_input_index][get_base_index('-')]);
            int vertical_distance = matrix[i - 1][j + offset] + scoring_matrix[get_base_index('-')][current_reference_index];
            int diagonal_distance = matrix[i - 1][j - 1 + offset] + scoring_matrix[current_input_index][current_reference_index];

            int value = max(horizontal_distance, vertical_distance, diagonal_distance, 0);
            std::cout << "(" << i << ", " << j << ") = " << value << std::endl;

            matrix[i][j] = value;

            if (matrix[i][j] > max_value)
            {
                max_value = matrix[i][j];
                max_index_i = i;
                max_index_j = j;
            }
        }
    }

    printf("        ");
    for (int i = 0; i < reference.length(); i++)
    {
        printf("%4c", reference[i]);
    }
    std::cout << std::endl;

    for (int i = 0; i < rows; i++)
    {
        if (i == 0)
        {
            printf("    ");
        }
        else
        {
            printf("%4c", input[i - 1]);
        }

        int offset = std::max(0, i - k);
        for (int j = 0; j < offset; j++)
        {
            printf("    ");
        }
        for (int j = 0; j < get_array_size_at_row(i, columns, k); j++)
        {
            int element = matrix[i][j];
            printf("%4d", element);
        }
        std::cout << std::endl;
    }
}

int main(int argc, char *argv[])
{
    auto start = std::chrono::high_resolution_clock::now();

    std::string args_reference_genome = argv[1];
    std::string args_sequenced_reads = argv[2];
    int k = atoi(argv[3]);

    std::ios::sync_with_stdio(false);
    std::string reference = read_reference(args_reference_genome);
    std::map<std::string, std::vector<int>> kmer_map;
    generate_kmers(kmer_map, reference, k);

    std::vector<std::string> sequence_list;
    read_sequence_data(sequence_list, args_sequenced_reads);
    append_reverse_complements(sequence_list);

    for (int i = 0; i < sequence_list.size(); i++)
    {
        auto result = map_to_reference(kmer_map, reference, sequence_list[i], k);
        std::cout << i << " " << result.first << "-" << result.second << std::endl;
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::cout << std::endl
              << "Execution time: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(finish - start).count() / 1e6
              << "ms\n";

    /*std::string reference = "TATATGCGGCGTTT";
    std::string input = "GGTATGCTGGCGCTA";

    apply_local_allign(reference, input);
    apply_local_allign_optimized(reference, input, 3);
    */

    return 0;
}