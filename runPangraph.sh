#!/bin/bash
#SBATCH -J pangraph
#SBATCH -A bag.prj
#SBATCH -p short

echo "****************************************************"
echo "Slurm Job ID: $SLURM_JOB_ID"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "****************************************************"

cd /well/bag/fnd111/ges
mkdir -p pangraph-output
samples=./flanker-output-5000/*.fasta

pangraph=/well/bag/fnd111/pangraph/bin/pangraph
blast=/well/bag/fnd111/miniconda3/bin/makeblastdb

$pangraph build $samples > ./pangraph-output/pangraph-build.json
$pangraph export --edge-minimum-length 0 ./pangraph-output/pangraph-build.json -p pangraph-export -o ./pangraph-output/

$blast -in ./pangraph-output/pangraph-export.fa -dbtype 'nucl'
geneBlock=$(blastn -query ./NG_049137.1.fasta -db ./pangraph-output/pangraph-export.fa -outfmt 7)
echo $geneBlock

python ./scripts/pangraphGFA.py ./pangraph-output/pangraph-export.gfa

echo "Finished at: "`date`
echo "****************************************************"

echo "*****************************************************"
echo "Finished!"
echo "*****************************************************"
