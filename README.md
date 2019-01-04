# Bioinformatics-MutationDetection
 Project assignment as part of the course https://www.fer.unizg.hr/en/course/bio

### Usage:

Compile Windows: `g++ -I src/ src/*.c src/*.cpp -O3 -lpsapi -m64 -o main.exe`

Compile Linux/OSX: `g++ -I src/ src/*.c src/*.cpp -O3 -m64 -o main.exe`

Run: `./main.exe -r train_data/lambda.fasta -s train_data/lambda_simulated_reads.fasta -o train_data/lambda_result.csv`

Jaccard: `python train_data/jaccard.py -b train_data/lambda_mutated.csv -a train_data/lambda_result.csv`

Help: `./main.exe -h`

### Available flags:
```
        -r   [reference genome file path] (file_path)
        -s   [sequenced reads file path] (file_path)
        -o   [output file path] (file_path)
        -k   [kmer length, default=8] (int)
        -rd  [region divider, default=100] (int)
        -lak [local alignment k value, default=50] (int)
        -cc  [confirmation count, default=3] (int)
        -m   [match score, default=5] (int)
        -x   [mismatch score, default=-4] (int)
        -g   [gap score, default=-7] (int)
```
