#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=72:00:00
#SBATCH --mem=250G

echo "starting run at: `date`"

####### Set environment variables ###############
module load gcc/10.2.0
module load lib/zlib/1.2.11
module load perl/5.32.0

PATH=/work/soghigian_lab/apps/bedtools2/bin:$PATH
PATH=/work/soghigian_lab/apps/samtools-1.15/bin:$PATH
PATH=/work/soghigian_lab/apps/gatb-minia-pipeline:$PATH
PATH=/work/soghigian_lab/apps/bwa:$PATH

source ~/software/init-conda-sl
conda activate /work/soghigian_lab/apps/conda/envs/gatb_minia

file_list=/work/soghigian_lab/huiqing.yeo/genome_assembly/ahe_sample.txt
trimmed=/work/soghigian_lab/huiqing.yeo/trimmed/ahe/

cd /scratch/$SLURM_JOB_ID
echo results in $SLURM_JOB_ID

#gatb minia
for sample in $(cat $file_list)
do
gatb -1 ${trimmed}${sample}_outfile1.fq.gz -2 ${trimmed}${sample}_outfile2.fq.gz \
-s ${trimmed}${sample}_s.fq.gz -o ${sample} \
--max-memory 250000 --nb-cores 16 --no-scaffolding
done

echo "Job finised with exit code $? at: `date`"
