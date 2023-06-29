#!/bin/bash
#SBATCH -J mobtyper
#SBATCH -A bag.prj
#SBATCH -p short
#SBATCH --array 1-115:1

# conda activate mobtyper_env

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
mobtyper=/well/bag/fnd111/miniconda3/envs/mobtyper_env/bin/mob_typer

samples=./contigs-blastx-dedup.txt
f=$(sed -n "$SLURM_ARRAY_TASK_ID"p $samples)
echo "Working on ${f}"

mkdir -p ./mobtyper-output/"$f"

$mobtyper --infile ./contigs-blastx/"$f".fasta --out_file ./mobtyper-output/"$f"/"$f"_mobtyper.txt

echo "Finished at: "`date`
echo "****************************************************"

echo "*****************************************************"
echo "Finished!"
echo "*****************************************************"
