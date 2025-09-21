#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --mem=250G

echo "starting run at: `date`"

####### Set environment variables ###############
PATH=/work/soghigian_lab/apps/SPAdes-3.15.5-Linux/bin/:$PATH

list=/work/soghigian_lab/huiqing.yeo/genome_assembly/ahe_sample.txt
trimmed=/work/soghigian_lab/huiqing.yeo/trimmed/ahe/

cd /scratch/$SLURM_JOB_ID
echo results in $SLURM_JOB_ID

for sample in $(cat ${list})
do
spades.py -1 ${trimmed}${sample}_outfile1.fq.gz -2 ${trimmed}${sample}_outfile2.fq.gz \
-s ${trimmed}${sample}_s.fq.gz -o ${sample} -t 16 -m 250
done

echo "Job finised with exit code $? at: `date`"
