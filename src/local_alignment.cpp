// Create by Josip

#include <local_alignment.hpp>
#include <string>
#include <algorithm>
#include <iostream>

LocalAlignment::LocalAlignment(ReadMapper &readMapper, std::string &reference, std::map<int, std::map<std::string, int>> &result_map, int local_align_k, int region_divider, int match_score, int change_score, int gap_score)
    : readMapper(readMapper),
      reference(reference),
      result_map(result_map),
      region_divider(region_divider),
      local_align_k(local_align_k),
      match_score(match_score),
      change_score(change_score),
      gap_score(gap_score)
{
    scoring_matrix[A][A] = match_score;
    scoring_matrix[A][C] = change_score;
    scoring_matrix[A][G] = change_score;
    scoring_matrix[A][T] = change_score;

    scoring_matrix[C][A] = change_score;
    scoring_matrix[C][C] = match_score;
    scoring_matrix[C][G] = change_score;
    scoring_matrix[C][T] = change_score;

    scoring_matrix[G][A] = change_score;
    scoring_matrix[G][C] = change_score;
    scoring_matrix[G][G] = match_score;
    scoring_matrix[G][T] = change_score;

    scoring_matrix[T][A] = change_score;
    scoring_matrix[T][C] = change_score;
    scoring_matrix[T][G] = change_score;
    scoring_matrix[T][T] = match_score;

    scoring_matrix[A][INDEL] = gap_score;
    scoring_matrix[C][INDEL] = gap_score;
    scoring_matrix[G][INDEL] = gap_score;
    scoring_matrix[T][INDEL] = gap_score;

    scoring_matrix[INDEL][A] = gap_score;
    scoring_matrix[INDEL][C] = gap_score;
    scoring_matrix[INDEL][G] = gap_score;
    scoring_matrix[INDEL][T] = gap_score;
}

bool LocalAlignment::align(std::string &read_input)
{
    auto result = readMapper.map(read_input);

    int region_start = result.start_index * region_divider;
    std::string reference_region = get_substr(reference, region_start, (result.end_index - result.start_index) * region_divider);

    return apply_local_allign(reference_region, read_input, region_start);
}

bool LocalAlignment::apply_local_allign(std::string &reference_region, std::string &input, int reference_offset)
{
    // Matrix initialization
    int rows = input.length() + 1;
    int columns = reference_region.length() + 1;

    Score **matrix = new Score *[rows];
    for (int i = 0; i < rows; i++)
    {
        matrix[i] = new Score[get_array_size_at_row(i, columns)];
    }

    // Used for finding the cell element with the maximum score
    Max_Cell max_cell = {0, 0, 0};

    Score *horizontal_distance = new Score(0, Insertion);
    Score *vertical_distance = new Score(0, Deletion);
    Score *diagonal_distance = new Score(0, Match);

    int prev_index_array_size = get_array_size_at_row(0, columns);

    // Starting from 1 because the first row contains zeros
    for (int i = 1; i < rows; i++)
    {
        unsigned char current_input_character = input[i - 1];

        int start_reference_index = std::max(0, i - local_align_k - 1);
        // offset: indicates if the first row element is horizontally moved in relation to the first one in the previous row
        int offset = (i <= local_align_k) ? 0 : 1;
        int array_size = get_array_size_at_row(i, columns);
        int start_array_index = offset == 0 ? 1 : 0;

        for (int j = start_array_index; j < array_size; j++)
        {
            // current character from the reference string => start index + offset
            unsigned char current_reference_character = reference_region[start_reference_index + j - start_array_index];

            horizontal_distance->score = j == 0 ? none->score : (matrix[i][j - 1].score + scoring_matrix[current_input_character][INDEL]);

            // no vertical value if it's the last element in the current row and the previous row length is smaller or equal than the current one
            bool no_vertical_value = j == (array_size - 1) && prev_index_array_size <= array_size;
            vertical_distance->score = no_vertical_value ? 0 : matrix[i - 1][j + offset].score + scoring_matrix[INDEL][current_reference_character];

            diagonal_distance->score = matrix[i - 1][j - 1 + offset].score + scoring_matrix[current_input_character][current_reference_character];

            if (diagonal_distance->score > vertical_distance->score && diagonal_distance->score > horizontal_distance->score)
            {
                matrix[i][j] = *diagonal_distance;
            }
            else if (vertical_distance->score > horizontal_distance->score)
            {
                matrix[i][j] = *vertical_distance;
            }
            else
            {
                matrix[i][j] = *horizontal_distance;
            }

            // update the info about the max cell if the current cell is the biggest so far
            if (matrix[i][j].score > max_cell.value)
            {
                max_cell.value = matrix[i][j].score;
                max_cell.i = i;
                max_cell.j = j;
            }
        }
        prev_index_array_size = array_size;
    }

    bool backtrack_result = backtrack(matrix, max_cell, columns, input, reference_region, reference_offset);

    // free the memory
    for (int i = 0; i < rows; i++)
    {
        delete matrix[i];
    }
    delete[] matrix;

    return backtrack_result;
}

bool LocalAlignment::backtrack(Score **matrix, Max_Cell max_cell, int columns, std::string &input, std::string &reference_region, int reference_offset)
{
    // Initial values: the (i,j) coordinates from the cell with the greatest value
    int i = max_cell.i;
    int j = max_cell.j;

    std::string input_result;
    input_result.reserve(i + j);
    std::string reference_result;
    reference_result.reserve(i + j);

    // constructing the most similar strings - reversed
    int row_offset = (i <= local_align_k) ? 0 : 1;
    while (i > 0 && j > 0)
    {
        Score score = matrix[i][j];
        if (score.action == Insertion)
        {
            j--;
            input_result += '-';
            reference_result += reference_region[std::max(0, i - local_align_k) + j];
        }
        else if (score.action == Deletion)
        {
            i--;
            row_offset = (i <= local_align_k) ? 0 : 1;
            j += row_offset;
            input_result += input[i];
            reference_result += '-';
        }
        else if (score.action == Match)
        {
            i--;
            j = j - 1 + row_offset;
            row_offset = (i <= local_align_k) ? 0 : 1;
            input_result += input[i];
            reference_result += reference_region[std::max(0, i - local_align_k) + j];
        }
        else
        {
            break;
        }
    }

    // reversing the reversed stirngs
    std::reverse(input_result.begin(), input_result.end());
    std::reverse(reference_result.begin(), reference_result.end());

    double similarity = (1 - (float)get_distance(reference_result, input_result) / reference_result.size());
    // Only accept the matched stirng if their similarity is greater than 80%
    if (similarity < 0.8)
    {
        return false;
    }

    int insertion_count = 0;
    int offset = j;
    // process the differences between the matched strings and update the result map
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

int LocalAlignment::get_array_size_at_row(int row, int max_columns)
{
    return std::max(0, std::min(max_columns - 1, row + local_align_k) - std::max(0, row - local_align_k) + 1);
}

int LocalAlignment::get_distance(std::string &s1, std::string &s2)
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

int LocalAlignment::get_substr_length(int max_length, int start, int length)
{
    return start + length > max_length ? max_length - start : length;
}

std::string LocalAlignment::get_substr(std::string &input, int start, int length)
{
    int input_start = start;
    int input_length = get_substr_length(input.size(), input_start, length);
    return input.substr(input_start, input_length);
}
