# GES-5 analysis workflow

## Data curation

- We retrived 1,375 contigs from NCBI's Pathogen Detection [MicroBIGG-E](https://www.ncbi.nlm.nih.gov/pathogens/microbigge) using the query `element_symbol:blaGES*` (as of 27/02/2023).
- Contigs were annotated for GES-5 using the NCBI protein reference sequence [WP_012658785.1](https://www.ncbi.nlm.nih.gov/protein/WP_012658785.1) with blastx. We filtered by 100% identity/coverage and a single hit, giving 431 contigs.
- Contigs were then deduplicated using (i) BioProject and (ii) sequence containment. For contigs from the same BioProject, contained in one another, the longest, container contig was always chosen, to maximse information per study. This approach is less conservative than keeping all unique contigs per BioProject, but makes the assumption that the contained contig likely continued homologously to the container contig.

BioProject accessions were retreived using [E-utilites](https://www.ncbi.nlm.nih.gov/books/NBK179288/):
`cat BioSamples.txt | { 
while read line\m
do
echo "$line" | tr '\n' '\t'
esearch -db biosample -query $line < /dev/null | elink -target bioproject | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc | tr '\n' '\t' && echo ""
sleep 1
done
} < BioProjects.tsv`
