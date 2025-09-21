#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=60:00:00
#SBATCH --mem=250G

echo "starting run at: `date`"

####### Set environment variables ###############
module load gcc/10.2.0
module load lib/zlib/1.2.11

PATH=/work/soghigian_lab/apps/samtools-1.15/:$PATH
PATH=/work/soghigian_lab/apps/picard.jar/:$PATH
source ~/software/init-conda-sl
conda activate discovardenovo

file_list=/work/soghigian_lab/huiqing.yeo/genome_assembly/ahe_sample.txt
trimmed=/work/soghigian_lab/huiqing.yeo/trimmed/ahe/

cd /scratch/$SLURM_JOB_ID
echo results in $SLURM_JOB_ID

#Discovar assembly
for sample in $(cat $file_list)
do
DiscovarDeNovo READS=${trimmed}${sample}_outfile1.fq.gz,${trimmed}${sample}_outfile2.fq.gz OUT_DIR=${sample} NUM_THREADS=16
done

echo "Job finised with exit code $? at: `date`"
