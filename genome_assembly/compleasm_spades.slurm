#!/bin/bash
####### Reserve computing resources #############
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=12:00:00
#SBATCH --mem=10G

####### Set environment variables ###############
module load bioconda/conda3
source activate /work/soghigian_lab/apps/conda/envs/compleasm

list=/work/soghigian_lab/huiqing.yeo/genome_assembly/all_samples.txt
assembly=/work/soghigian_lab/huiqing.yeo/genome_assembly/spades/

cd ${assembly}

####### Script ###############
for file in $(cat $list)
do
compleasm run -a ${assembly}${file}_scaffolds.fasta -o ${file}_compleasm -t 12 \
-l diptera -L /work/soghigian_lab/apps/busco_lineages/diptera_odb10 
done
