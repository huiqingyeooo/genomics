#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --mem=250G

echo "starting run at: `date`"

####### Set environment variables ###############
source ~/software/init-conda-sl
conda activate megahit

list=/work/soghigian_lab/huiqing.yeo/genome_assembly/ahe_sample.txt
trimmed=/work/soghigian_lab/huiqing.yeo/trimmed/ahe/

cd /scratch/$SLURM_JOB_ID
echo results in $SLURM_JOB_ID

for sample in $(cat ${list})
do
megahit -1 ${trimmed}${sample}_outfile1.fq.gz -2 ${trimmed}${sample}_outfile2.fq.gz \
-r ${trimmed}${sample}_s.fq.gz -o ${sample}
done

echo "Job finised with exit code $? at: `date`"
