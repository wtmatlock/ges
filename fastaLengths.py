import sys
from Bio.SeqIO.FastaIO import SimpleFastaParser

fastaFile = open(sys.argv[1], 'r')
lengthsFile= open('./lengths.tsv', 'w')

for name, seq in SimpleFastaParser(fastaFile):
    seqLength = len(seq)
    lengthsFile.write(name + '\t' + str(seqLength) + '\n')

fastaFile.close()
lengthsFile.close()