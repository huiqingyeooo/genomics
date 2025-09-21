#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:00
#SBATCH --mem=50G

####### Set environment variables ###############
source ~/software/init-conda-sl
conda activate /work/soghigian_lab/apps/conda/envs/atram
atramPath=/work/soghigian_lab/apps/conda/envs/atram/aTRAM
raw_reads=/work/soghigian_lab/data/ahe/reads/paired_end
trimmed=/work/soghigian_lab/huiqing.yeo/atram/reads

#env PERL5LIB= PERL_LOCAL_LIB_ROOT= parallel sbatch atram.slurm :::: list.txt :::: ak
cd /work/soghigian_lab/huiqing.yeo/atram/

${atramPath}/atram_preprocessor.py \
--blast-db=/work/soghigian_lab/huiqing.yeo/atram/blastdb/${1} \
--end-1=${trimmed}/${1}_outfile1.fq.gz \
--end-2=${trimmed}/${1}_outfile2.fq.gz \
--gzip

${atramPath}/atram.py \
--blast-db=/work/soghigian_lab/huiqing.yeo/atram/blastdb/${1} \
--query=/work/soghigian_lab/huiqing.yeo/atram/${2}.fasta \
--assembler=spades \
--path=/work/soghigian_lab/apps/SPAdes-4.0.0-Linux/bin \
-t tmp \
--keep-temp-dir \
--output-prefix=/work/soghigian_lab/huiqing.yeo/atram/${1}_${2}

echo "Job finised with exit code $? at: `date`"
