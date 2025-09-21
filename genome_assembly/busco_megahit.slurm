#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --mem=100G

echo "starting run at: `date`"

####### Set environment variables ###############
#Run this first:
source ~/software/init-conda-sl
conda activate /work/soghigian_lab/apps/miniconda3/envs/busco_env

assembly=/work/soghigian_lab/huiqing.yeo/genome_assembly/megahit/
list=/work/soghigian_lab/huiqing.yeo/genome_assembly/all_samples.txt

#try in batch mode; do it in scratch
cd /scratch/${SLURM_JOB_ID}
echo check /scratch/${SLURM_JOB_ID}

####### Busco ###############
#run busco in batch mode (autodetects)
for file in $(cat $list)
do
busco -m genome -c 16 --offline --augustus -f \
-i ${assembly}${file}_final.contigs.fa \
-o ${file}_busco \
-l /work/soghigian_lab/apps/busco_lineages/diptera_odb10
done

echo "Job finised with exit code $? at: `date`"
