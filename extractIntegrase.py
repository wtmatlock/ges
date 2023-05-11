import os

# open the .tsv file and read in the data
with open('sys.argv[1], 'r') as f:
    data = [line.strip().split('\t') for line in f.readlines()]

# loop over each row in the .tsv file
for row in data:
    # extract the data
    name = row[0]
    start = int(row[1])
    end = int(row[2])

    # open the corresponding .fasta file and read in the sequence
    with open(f'./{name}.fasta', 'r') as f:
        fasta_data = f.readlines()
    seq = ''.join(fasta_data[1:]).replace('\n', '')

    # extract the desired subsequence
    subseq = seq[start-1:end]

    # Write the subsequence to a new .fasta file
    with open(f'./contigs/{name}-integrase.fasta', 'w') as f:
        f.write(f'>{name}_integrase\n{subseq}\n')
