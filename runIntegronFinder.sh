#!/bin/bash
#SBATCH -J integronFinder
#SBATCH -A bag.prj
#SBATCH -p short
#SBATCH --array 1-109:1

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

integron_finder ./contigs/"$f".fasta --keep-tmp --local-max --promoter-attI --lin

echo "Finished at: "`date`
echo "****************************************************"

echo "*****************************************************"
echo "Finished!"
echo "*****************************************************"
