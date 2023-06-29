#!/bin/bash
#SBATCH -J plasmidfinder
#SBATCH -A bag.prj
#SBATCH -p short
#SBATCH --array 1-115:1

# conda activate flanker_env

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
abricate=/well/bag/fnd111/miniconda3/envs/flanker_env/bin/abricate

samples=./contigs-blastx-dedup.txt
f=$(sed -n "$SLURM_ARRAY_TASK_ID"p $samples)
echo "Working on ${f}"

mkdir -p ./plasmidfinder-output/"$f"

$abricate ./contigs-blastx/"$f".fasta --db plasmidfinder > ./plasmidfinder-output/"$f"/"$f"_plasmidfinder.tab
$abricate --summary ./plasmidfinder-output/*/*_plasmidfinder.tab > ./plasmidfinder-output/plasmidfinder_summary.tab

echo "Finished at: "`date`
echo "****************************************************"

echo "*****************************************************"
echo "Finished!"
echo "*****************************************************"
