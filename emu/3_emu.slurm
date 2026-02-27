#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=1-00:00:00
#SBATCH --mem=100G

####### Set environment variables ###############
source ~/software/init-conda-sl
conda activate /work/soghigian_lab/apps/conda/envs/emu

sample=sample5_vsearch
db_name=ref_nem_cerv
database=/work/soghigian_lab/huiqing.yeo/filiaroidDiv/emu/ref_database

echo "starting run at: `date`"
####### script ##################################

# fasta file
#emu abundance ${sample}.fasta \
#--db ${database}/${db_name} --threads 24 --output-dir ${sample}_${db_name}_t24 \
#--type map-ont --keep-counts --keep-read-assignments

# fastq file
emu abundance --type sr ${sample}.fastq \
--db ${database}/${db_name} --threads 24 --output-dir ${sample}_${db_name}_t24 \
--type map-ont --keep-counts --keep-read-assignments

echo "Job finised with exit code $? at: `date`"
