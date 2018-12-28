#include <sequence_reader.hpp>
#include <fstream>

std::vector<std::string> SequenceReader::read_sequence_data(std::string path)
{
    std::vector<std::string> sequence_list;
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

    return sequence_list;
}