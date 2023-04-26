# GES-5 flank analysis workflow

## Data curation

- We retrived 1,375 contigs from NCBI's Pathogen Detection [MicroBIGG-E](https://www.ncbi.nlm.nih.gov/pathogens/microbigge) using the query `element_symbol:blaGES*` (as of 27/02/2023). The full metadata table was also downloaded.
- Contigs were annotated for GES-5 using the NCBI GES-5 protein reference sequence [WP_012658785.1](https://www.ncbi.nlm.nih.gov/protein/WP_012658785.1) and blastx. We filtered by 100% identity/coverage and a single hit, giving 431 contigs.
- Contigs were then deduplicated as follows: first, we took the longest contig from each BioSample, discarding the rest. In the case of ties, we chose a random representative. Then, within each BioProject, we discarded any contigs that were perfectly contained in another. If this resulted in all contigs from a BioProject being discarded, we chose a random, longest representative.
- The script `deduplication.R` was used to wrangle the various inputs (NCBI metadata, contig lengths, BioProject accessions, mash screen) then deduplicate the contigs.

### Contig lengths

Contig sequence lengths were calculated using `fastaLengths.py`. Usage is `python3 fastaLengths.py contigs.fa` which outputs `lengths.tsv`. Requires `SimpleFastaParser` from [Biopython](https://github.com/biopython/biopython).

### BioProject accessions

BioProject accessions were retreived using NCBI's [E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK179288/). The following code writes a tab-delimited table of BioProjects for a given list of BioSamples:
```
cat BioSamples.txt | { 
while read line
do
echo "$line" | tr '\n' '\t'
esearch -db biosample -query $line < /dev/null | elink -target bioproject | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc | tr '\n' '\t' && echo ""
sleep 1
done
} > BioProjects.tsv
```
BioProject accessions were manually curated to remove any meta-analyses which reused contigs from earlier studies. 

### Sequence containment
Pairwise contig containment used [Mash](https://github.com/marbl/Mash) (v. 2.2):

```
mash sketch -s 1000000 ./contigs/*.fasta -o contigs
echo contigs/*.fasta | xargs -n 1 mash screen contigs.msh > mash-output.tsv
```
With this command, the Mash output will only list the 'screened' contigs, not the contig that we are scoring containment in. For example,
```
x
a.fasta
b.fasta
c.fasta
a.fasta
b.fasta
c.fasta
a.fasta
b.fasta
c.fasta
```
really means scoring the containment of x in y:
``` 
x y
a.fasta a.fasta
b.fasta a.fasta
c.fasta a.fasta
a.fasta b.fasta
b.fasta b.fasta
c.fasta b.fasta
a.fasta c.fasta
b.fasta c.fasta
c.fasta c.fasta
```
A sense check here is that given *m* contigs, every *(n-1)m+n*-th row should score 1 since we are screening a contig against itself.

## Extracting flanking sequences

- Flanking sequences were extracted using the python tool [Flanker](https://github.com/wtmatlock/flanker) with a custom [Abricate](https://github.com/tseemann/abricate) database containing only the NCBI GES-5 nucleotide reference sequence [NG_049137.1](https://www.ncbi.nlm.nih.gov/nuccore/NG_049137.1). Flanker was run by the Slurm scheduler script `runFlanker.sh` within the Conda environment described [here](https://flanker.readthedocs.io/en/latest/#installation). For a contig called `foo.fasta` and `--window bar`, the output file will be called `foo_GES-5_bar_both_flank.fasta`.

## Sequence annotations

Contig annotations were pre-generated by NCBI's [Prokaryotic Genome Annotation Pipeline](https://www.ncbi.nlm.nih.gov/genome/annotation_prok/), and retrived as follows:

```
mkdir gffs
cat contigAccessions.txt | { 
while read line
do
wget "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=${line}" -O ./gffs/"$line".gff
sleep 1
done
}
```
Additional, integron-specific annoations were generated using [Integron Finder](https://github.com/gem-pasteur/Integron_Finder) (v. 2.0.2) with the command
```
integron_finder ./contigs/*.fasta --gbk --local-max --promoter-attI --lin
```
For each contig, this wrote a `.integron` file, which we combined with the NCBI annotations for plotting.

## Clustering flanking sequences
Pangraph... 

## Bringing it all together
- Plotting the data...
- Contig metadata countries were mapped to WHO region using this [table](https://github.com/lukes/ISO-3166-Countries-with-Regional-Codes/blob/master/all/all.csv).


