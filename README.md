# Bioinformatics-MutationDetection
 Project assignment as part of the course https://www.fer.unizg.hr/en/course/bio

### Usage:

Compile Windows: `g++ -I src/ src/*.c src/*.cpp -O3 -lpsapi -m64 -std=c++11 -o main.exe`

Compile Linux/OSX: `g++ -I src/ src/*.c src/*.cpp -O3 -m64 -std=c++11 -o main.out`

Run: `./main.exe -r {PATH_TO_REFERENCE_GENOME} -s {PATH_TO_SEQUENCED_READS} -o {PATH_TO_OUTPUT_FILE}`

Jaccard: `python jaccard.py -b {PATH_TO_OUTPUT_FILE} -a {PATH_TO_REFERENCE_FILE}`

Help: `./main.exe -h`

### Available flags:
```
        -r   [reference genome file path] (file_path)
        -s   [sequenced reads file path] (file_path)
        -o   [output file path] (file_path)
        -k   [kmer length, default=8] (int)
        -rd  [region divider, default=100] (int)
        -lak [local alignment k value, default=100] (int)
        -cc  [confirmation count, default=3] (int)
        -m   [match score, default=5] (int)
        -x   [mismatch score, default=-4] (int)
        -g   [gap score, default=-7] (int)
```