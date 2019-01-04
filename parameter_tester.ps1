For ($lak = 45; $lak -lt 56; $lak++ ) {
    For ($k = 4; $k -lt 12; $k++ ) {
        For ($rd = 80; $rd -lt 130; $rd += 10) {
            Write-Output "Local alignment k value: " $lak
            Write-Output "K-mer length: " $k
            Write-Output "Region divider: " $rd

            ./main.exe -r train_data/lambda.fasta -s train_data/lambda_simulated_reads.fasta -o train_data/lambda_result.csv -lak $lak -k $k -rd $rd
            python train_data/jaccard.py -b train_data/lambda_mutated.csv -a train_data/lambda_result.csv

            Write-Output "--------------------------------------------------------"
        }
    }
}