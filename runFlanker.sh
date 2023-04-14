#!/bin/bash
#SBATCH -J flanker
#SBATCH -A bag.prj
#SBATCH -p short
#SBATCH --array 1-115:1

echo "****************************************************"
echo "Slurm Job ID: $SLURM_JOB_ID"
echo "Slurm Task ID: $SLURM_ARRAY_TASK_ID"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "****************************************************"

cd /well/bag/fnd111/ges/
samples=./contigs-blastx-dedup.txt
f=$(sed -n "$SLURM_ARRAY_TASK_ID"p $samples)
echo "Working on ${f}"

# the '--window' parameter controls the flank length extracted from the contig
flanker --gene GES-5 --database GES-5 --fasta_file ./contigs-blastx/"$f".fasta --flank both --window 5000 --include_gene

echo "Finished at: "`date`
echo "****************************************************"

echo "*****************************************************"
echo "Finished!"
echo "*****************************************************"
