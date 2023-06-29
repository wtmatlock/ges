#!/bin/bash
#SBATCH -J mlst
#SBATCH -A bag.prj
#SBATCH -p short
#SBATCH --array 1-115:1

# conda activate mlst_env

echo "****************************************************"
echo "Slurm Job ID: $SLURM_JOB_ID"
echo "Slurm Task ID: $SLURM_ARRAY_TASK_ID"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "****************************************************"

echo "****************************************************"
echo "Starting at: "`date`

cd /well/bag/fnd111/ges
mlst=/well/bag/fnd111/miniconda3/envs/mlst_env/bin/mlst

samples=./contigs-blastx-dedup.txt
f=$(sed -n "$SLURM_ARRAY_TASK_ID"p $samples)
echo "Working on ${f}"

mkdir -p ./mlst-output/"$f"

$mlst ./contigs-blastx/"$f".fasta > ./mlst-output/"$f"/"$f"_mlst.tab

echo "Finished at: "`date`
echo "****************************************************"

echo "*****************************************************"
echo "Finished!"
