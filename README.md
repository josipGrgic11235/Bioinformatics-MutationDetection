# Bioinformatics-MutationDetection
 Project assignment as part of the course https://www.fer.unizg.hr/en/course/bio

### Usage:

Compile: `g++ -I src/ src/*.cpp -o main.exe`

Run: `./main.exe -r train_data/lambda.fasta -s train_data/lambda_simulated_reads.fasta -o train_data/lambda_result.csv`

Jaccard: `python train_data/jaccard.py -b train_data/lambda_mutated.csv -a train_data/lambda_result.csv`

Help: `./main.exe -h`