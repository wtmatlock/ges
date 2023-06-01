import os
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein

# open the TSV file and read in the data
with open('/Users/willmatlock/Desktop/integrase-positions.tsv', 'r') as f:
    data = [line.strip().split('\t') for line in f.readlines()]

# loop over each row in the TSV file
for row in data:
    # extract the relevant data
    name = row[0]
    start = int(row[1])
    end = int(row[2])
    strand = int(row[3])

    # open the corresponding FASTA file and read in the sequence
    with open(f'/Users/willmatlock/Desktop/ges/rescomp-backup/contigs/{name}.fasta', 'r') as f:
        fasta_data = f.readlines()
    seq = ''.join(fasta_data[1:]).replace('\n', '')

    # extract the desired subsequence
    subseq = seq[start-1:end]

    # reverse complement the subsequence if the strand is -1
    if strand == -1:
        subseq = str(Seq(subseq).reverse_complement())

    # write the subsequence to a new FASTA file
    with open(f'/Users/willmatlock/Desktop/integrase-sequences/{name}-integrase.fasta', 'w') as f:
        f.write(f'>{name}-integrase\n{subseq}\n')

    # translate the subsequence to protein and write it to a new FASTA file
    prot_seq = Seq(subseq, generic_dna).translate()
    with open(f'/Users/willmatlock/Desktop/integrase-sequences/{name}-integrase-prot.fasta', 'w') as f:
        f.write(f'>{name}-integrase-prot\n{prot_seq}\n')
