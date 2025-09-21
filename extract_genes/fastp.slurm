#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:00
#SBATCH --mem=50G

####### Set environment variables ###############
source ~/software/init-conda-sl
cd /work/soghigian_lab/apps/aTRAM
conda activate aTRAM
PATH=/work/soghigian_lab/apps/fastp:$PATH

#env PERL5LIB= PERL_LOCAL_LIB_ROOT= parallel sbatch fastp_se.slurm :::: list_se.txt


####### Script ######################################
raw_reads=/work/soghigian_lab/data/ahe/reads/single_end
trimmed=/work/soghigian_lab/huiqing.yeo/geneminer/reads

taxname=`echo ${1}_R1.fq.gz | awk -F"_" '{print $1}'`

fastp -i ${raw_reads}/${1}*.fastq.gz -o ${trimmed}/${1}_outfile.fq.gz --detect_adapter_for_pe --cut_front -g -x -w 12 -h ${1}.html

echo "Job finised with exit code $? at: `date`"
