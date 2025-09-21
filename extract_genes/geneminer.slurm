#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:00
#SBATCH --mem=50G

####### Set environment variables ###############
PATH=/home/huiqing.yeo/software/ncbi-blast-2.14.1+/bin:$PATH
trimmed=/work/soghigian_lab/huiqing.yeo/geneminer/reads

#env PERL5LIB= PERL_LOCAL_LIB_ROOT= parallel sbatch geneminer.slurm :::: list_pe.txt ::: ak coi

echo ${1}

/work/soghigian_lab/apps/GeneMiner/geneminer.py \
-1 ${trimmed}/${1}_outfile1.fq.gz \
-2 ${trimmed}/${1}_outfile2.fq.gz \
-rtfa /work/soghigian_lab/huiqing.yeo/atram/${2}.fasta \
-o ${2}_${1}

echo "Job finised with exit code $? at: `date`"
